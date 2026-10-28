#pragma once

#include "0_basic/Error.h"
#include "0_basic/MPI_WRAPPER.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace HALO_EXEC
{
struct Message
{
    int peer = -1;
    int tag = -1;
    int length = 0;
    std::size_t buffer_offset = 0;
    bool enabled = true;
};

struct Plan
{
    std::string name;
    std::uint64_t generation = 0;
    std::vector<Message> sends;
    std::vector<Message> recvs;
    std::size_t send_doubles = 0;
    std::size_t recv_doubles = 0;
};

struct Statistics
{
    std::uint64_t executions = 0;
    std::uint64_t workspace_growth_events = 0;
    std::uint64_t sent_messages = 0;
    std::uint64_t received_messages = 0;
    std::uint64_t sent_doubles = 0;
    std::uint64_t received_doubles = 0;
};

inline void finalize_plan(Plan &plan)
{
    plan.send_doubles = 0;
    for (std::size_t i = 0; i < plan.sends.size(); ++i)
    {
        Message &message = plan.sends[i];
        if (message.length < 0)
            ERROR::Abort("[Halo] negative prepared send length: " + plan.name);
        message.buffer_offset = plan.send_doubles;
        plan.send_doubles += static_cast<std::size_t>(message.length);
    }
    plan.recv_doubles = 0;
    for (std::size_t i = 0; i < plan.recvs.size(); ++i)
    {
        Message &message = plan.recvs[i];
        if (message.length < 0)
            ERROR::Abort("[Halo] negative prepared receive length: " + plan.name);
        message.buffer_offset = plan.recv_doubles;
        plan.recv_doubles += static_cast<std::size_t>(message.length);
    }
}

// A Halo deliberately shares this workspace among every prepared exchange,
// including Coupling. Halo calls are serialized; nested/concurrent use is a
// contract violation rather than a reason to allocate a second workspace.
class PreparedExchangeExecutor
{
public:
    void reset(std::uint64_t generation)
    {
        if (in_use_)
            ERROR::Abort("[Halo] cannot reset an exchange workspace in use");
        generation_ = generation;
        ready_ = false;
        max_send_doubles_ = max_recv_doubles_ = 0;
        max_send_messages_ = max_recv_messages_ = 0;
    }

    void invalidate() noexcept
    {
        ready_ = false;
        generation_ = 0;
    }

    void reserve_plan(const Plan &plan)
    {
        if (plan.generation != generation_)
            ERROR::Abort("[Halo] exchange generation mismatch: " + plan.name);
        max_send_doubles_ = std::max(max_send_doubles_, plan.send_doubles);
        max_recv_doubles_ = std::max(max_recv_doubles_, plan.recv_doubles);
        max_send_messages_ = std::max(max_send_messages_, plan.sends.size());
        max_recv_messages_ = std::max(max_recv_messages_, plan.recvs.size());
        grow_workspace_();
    }

    void activate()
    {
        grow_workspace_();
        ready_ = true;
    }

    bool ready(std::uint64_t generation) const noexcept
    {
        return ready_ && generation_ == generation;
    }

    const Statistics &statistics() const noexcept { return statistics_; }

    template <class Pack, class Unpack>
    void execute(const Plan &plan, Pack pack, Unpack unpack)
    {
        if (!ready_ || plan.generation != generation_)
            ERROR::Abort("[Halo] stale/inactive exchange plan: " + plan.name);
        if (in_use_)
            ERROR::Abort("[Halo] nested/concurrent exchange: " + plan.name);
        if (plan.send_doubles > send_buffer_.size() ||
            plan.recv_doubles > recv_buffer_.size() ||
            plan.sends.size() > send_requests_.size() ||
            plan.recvs.size() > recv_requests_.size())
            ERROR::Abort("[Halo] unreserved exchange plan: " + plan.name);

        struct Guard
        {
            explicit Guard(bool &flag) : flag_(flag) { flag_ = true; }
            ~Guard() { flag_ = false; }
            bool &flag_;
        } guard(in_use_);

        // Clear persistent sends so an unallocated block cannot leak values
        // left by an earlier exchange.
        std::fill(send_buffer_.begin(),
                  send_buffer_.begin() + plan.send_doubles, 0.0);
        for (std::size_t i = 0; i < plan.sends.size(); ++i)
        {
            const Message &m = plan.sends[i];
            send_requests_[i] = MPI_REQUEST_NULL;
            if (m.enabled)
                pack(i, send_pointer_(m.buffer_offset), m.length);
        }

        // One deterministic schedule: receives, sends, receive wait, send
        // wait, unpack. Zero-length messages remain valid MPI messages.
        for (std::size_t i = 0; i < plan.recvs.size(); ++i)
        {
            const Message &m = plan.recvs[i];
            recv_requests_[i] = MPI_REQUEST_NULL;
            if (m.enabled)
                PARALLEL::mpi_data_recv(m.peer, m.tag,
                    recv_pointer_(m.buffer_offset), m.length,
                    &recv_requests_[i]);
        }
        for (std::size_t i = 0; i < plan.sends.size(); ++i)
        {
            const Message &m = plan.sends[i];
            if (m.enabled)
                PARALLEL::mpi_data_send(m.peer, m.tag,
                    send_pointer_(m.buffer_offset), m.length,
                    &send_requests_[i]);
        }

        int nrecv = static_cast<int>(plan.recvs.size());
        int nsend = static_cast<int>(plan.sends.size());
        if (nrecv > 0)
            PARALLEL::mpi_wait(nrecv, recv_requests_.data(),
                               recv_statuses_.data());
        if (nsend > 0)
            PARALLEL::mpi_wait(nsend, send_requests_.data(),
                               send_statuses_.data());

        for (std::size_t i = 0; i < plan.recvs.size(); ++i)
        {
            const Message &m = plan.recvs[i];
            if (m.enabled)
                unpack(i, recv_pointer_(m.buffer_offset), m.length);
        }

        ++statistics_.executions;
        for (std::size_t i = 0; i < plan.sends.size(); ++i)
            if (plan.sends[i].enabled)
            {
                ++statistics_.sent_messages;
                statistics_.sent_doubles += plan.sends[i].length;
            }
        for (std::size_t i = 0; i < plan.recvs.size(); ++i)
            if (plan.recvs[i].enabled)
            {
                ++statistics_.received_messages;
                statistics_.received_doubles += plan.recvs[i].length;
            }
    }

private:
    double *send_pointer_(std::size_t offset) noexcept
    {
        return send_buffer_.empty() ? nullptr
                                    : send_buffer_.data() + offset;
    }

    double *recv_pointer_(std::size_t offset) noexcept
    {
        return recv_buffer_.empty() ? nullptr
                                    : recv_buffer_.data() + offset;
    }

    void grow_workspace_()
    {
        bool grew = false;
        if (send_buffer_.size() < max_send_doubles_)
        {
            send_buffer_.resize(max_send_doubles_);
            grew = true;
        }
        if (recv_buffer_.size() < max_recv_doubles_)
        {
            recv_buffer_.resize(max_recv_doubles_);
            grew = true;
        }
        if (send_requests_.size() < max_send_messages_)
        {
            send_requests_.resize(max_send_messages_, MPI_REQUEST_NULL);
            send_statuses_.resize(max_send_messages_);
            grew = true;
        }
        if (recv_requests_.size() < max_recv_messages_)
        {
            recv_requests_.resize(max_recv_messages_, MPI_REQUEST_NULL);
            recv_statuses_.resize(max_recv_messages_);
            grew = true;
        }
        if (grew)
            ++statistics_.workspace_growth_events;
    }

    std::uint64_t generation_ = 0;
    bool ready_ = false;
    bool in_use_ = false;
    std::size_t max_send_doubles_ = 0, max_recv_doubles_ = 0;
    std::size_t max_send_messages_ = 0, max_recv_messages_ = 0;
    std::vector<double> send_buffer_, recv_buffer_;
    std::vector<MPI_Request> send_requests_, recv_requests_;
    std::vector<MPI_Status> send_statuses_, recv_statuses_;
    Statistics statistics_;
};
}
