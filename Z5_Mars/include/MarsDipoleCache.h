#pragma once
#include <cstdint>
#include <cstddef>
#include <string>
#include <vector>

namespace MARS
{
// Deterministic fingerprint for accidental cache staleness/corruption detection.
struct DipoleFingerprint
{
    std::uint64_t value=14695981039346656037ULL;
    void Bytes(const void *data,std::size_t size)
    {
        const auto *p=static_cast<const unsigned char*>(data);
        for(std::size_t i=0;i<size;++i){value^=p[i];value*=1099511628211ULL;}
    }
    template<class T> void Add(const T &v){Bytes(&v,sizeof(v));}
};
// Native-endian binary cache, versioned and checked against the complete setup
// fingerprint and expected payload length. A miss never modifies values.
bool ReadDipoleCache(const std::string &path,std::uint64_t signature,
                     std::vector<double> &values);
bool WriteDipoleCache(const std::string &path,std::uint64_t signature,
                      const std::vector<double> &values);
}
