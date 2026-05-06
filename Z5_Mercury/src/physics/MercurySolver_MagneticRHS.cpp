#include "MercurySolver.h"
#include <vector>

void MercurySolver::AssembleRHS_Induction_CT_()
{
    const int nb = fld_->num_blocks();

    // 0) E_edge 清零
    for (int ib = 0; ib < nb; ++ib)
    {
        auto &Exi = fld_->field(fid_.fid_E.xi, ib);
        auto &Eeta = fld_->field(fid_.fid_E.eta, ib);
        auto &Eze = fld_->field(fid_.fid_E.zeta, ib);
        if (!Exi.is_allocated())
            continue;

        auto zero_electric = [&](FieldBlock &E)
        {
            Int3 lo = E.get_lo(), hi = E.get_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        E(i, j, k, 0) = 0.0;
        };
        zero_electric(Exi);
        zero_electric(Eeta);
        zero_electric(Eze);
    }

    // Explicit Advance
    {
        // 1) Calculate the E fields
        Build_E_explicit_edge_();

        // 2) E_edge 也要进边界/halo（否则 curl 更新 B_face 时边界附近会乱）
        mercury_bound_.Sync("Eedge"); // 你需要加一个 group：fields={E_xi,E_eta,E_zeta}
        FilterPoleNearAxisEedge_();
        mercury_bound_.Sync("Eedge");

        // 3) curl(E_edge) -> RHS_Bface，然后 Bface += dt*RHS
        for (int ib = 0; ib < nb; ++ib)
        {
            auto &Exi = fld_->field(fid_.fid_E.xi, ib);
            auto &Eeta = fld_->field(fid_.fid_E.eta, ib);
            auto &Eze = fld_->field(fid_.fid_E.zeta, ib);

            auto &RHSBxi = fld_->field(fid_.fid_RHS_b.xi, ib);
            auto &RHSBeta = fld_->field(fid_.fid_RHS_b.eta, ib);
            auto &RHSBze = fld_->field(fid_.fid_RHS_b.zeta, ib);
            if (!Exi.is_allocated())
                continue;

            // multiper 的符号你要和你的取向一致：
            // 一般 CT 是 B^{n+1} = B^n - dt * curl(E)
            CTOperators::CurlEdgeToFace(ib, Exi, Eeta, Eze, RHSBxi, RHSBeta, RHSBze, /*multiper=*/-1.0);
        }
    }
}

void MercurySolver::FilterPoleNearAxisEedge_()
{
    constexpr int shift_layers = 4;
    constexpr double filter_floor = 0.0;

    auto loc_delta = [](int edge_axis) -> Int3
    {
        if (edge_axis == 0)
            return {1, 0, 0};
        if (edge_axis == 1)
            return {0, 1, 0};
        return {0, 0, 1};
    };

    auto get_comp = [](const Int3 &a, int ax) -> int
    {
        if (ax == 0)
            return a.i;
        if (ax == 1)
            return a.j;
        return a.k;
    };

    auto set_comp = [](Int3 &a, int ax, int v)
    {
        if (ax == 0)
            a.i = v;
        else if (ax == 1)
            a.j = v;
        else
            a.k = v;
    };

    auto filter_one_edge = [&](FieldBlock &E,
                               int edge_axis,
                               const TOPO::PhysicalPatch &p)
    {
        if (!E.is_allocated())
            return;

        const int pole_dir = std::abs(p.direction);
        if (pole_dir != 1 && pole_dir != 2)
            return;

        const int norm_axis = pole_dir - 1;
        const bool high_side = (p.direction > 0);
        const Int3 delta = loc_delta(edge_axis);

        const Int3 elo = E.inner_lo();
        const Int3 ehi = E.inner_hi();

        Int3 base_lo = elo;
        Int3 base_hi = ehi;

        for (int ax = 0; ax < 3; ++ax)
        {
            if (ax == norm_axis)
                continue;

            const int d = get_comp(delta, ax);
            const int patch_lo = get_comp(p.this_box_node.lo, ax);
            const int patch_hi = get_comp(p.this_box_node.hi, ax) - d;

            set_comp(base_lo, ax, std::max(get_comp(elo, ax), patch_lo));
            set_comp(base_hi, ax, std::min(get_comp(ehi, ax), patch_hi));
        }

        const int nlo = get_comp(elo, norm_axis);
        const int nhi = get_comp(ehi, norm_axis);
        const int max_layers = std::min(shift_layers, nhi - nlo);

        if (max_layers <= 0)
            return;

        for (int layer = 0; layer < max_layers; ++layer)
        {
            Int3 lo = base_lo;
            Int3 hi = base_hi;

            const int n = high_side ? (nhi - 1 - layer) : (nlo + layer);
            set_comp(lo, norm_axis, n);
            set_comp(hi, norm_axis, n + 1);

            if (!(lo.i < hi.i && lo.j < hi.j && lo.k < hi.k))
                continue;

            double weight = 1.0 - static_cast<double>(layer) / static_cast<double>(max_layers);
            weight = std::max(filter_floor, std::min(1.0, weight));

            if (weight <= 0.0)
                continue;

            const int unique_k_hi = (edge_axis == 2) ? hi.k : (hi.k - 1);
            const int nk = unique_k_hi - lo.k;

            if (nk <= 1)
                continue;

            std::vector<double> old(static_cast<std::size_t>(nk));
            std::vector<double> filtered(static_cast<std::size_t>(nk));

            for (int i = lo.i; i < hi.i; ++i)
            {
                for (int j = lo.j; j < hi.j; ++j)
                {
                    for (int q = 0; q < nk; ++q)
                        old[static_cast<std::size_t>(q)] = E(i, j, lo.k + q, 0);

                    for (int q = 0; q < nk; ++q)
                    {
                        const int qm = (q + nk - 1) % nk;
                        const int qp = (q + 1) % nk;

                        const double smooth =
                            0.25 * old[static_cast<std::size_t>(qm)] +
                            0.50 * old[static_cast<std::size_t>(q)] +
                            0.25 * old[static_cast<std::size_t>(qp)];

                        filtered[static_cast<std::size_t>(q)] =
                            (1.0 - weight) * old[static_cast<std::size_t>(q)] +
                            weight * smooth;
                    }

                    for (int q = 0; q < nk; ++q)
                        E(i, j, lo.k + q, 0) = filtered[static_cast<std::size_t>(q)];

                    if (edge_axis != 2 && hi.k > unique_k_hi)
                        E(i, j, unique_k_hi, 0) = filtered[0];
                }
            }
        }
    };

    const int nb = fld_->num_blocks();
    for (int ib = 0; ib < nb; ++ib)
    {
        FieldBlock &Exi = fld_->field(fid_.fid_E.xi, ib);
        FieldBlock &Eeta = fld_->field(fid_.fid_E.eta, ib);
        FieldBlock &Eze = fld_->field(fid_.fid_E.zeta, ib);

        for (const auto &p : topo_->physical_patches)
        {
            if (p.this_block != ib)
                continue;
            if (p.bc_name != "Pole")
                continue;

            filter_one_edge(Exi, 0, p);
            filter_one_edge(Eeta, 1, p);
            filter_one_edge(Eze, 2, p);
        }
    }
}
