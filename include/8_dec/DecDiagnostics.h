#pragma once

#include <string>

#include "3_field/Field.h"
#include "3_field/FieldBlock.h"
#include "8_dec/DecTypes.h"

namespace DEC
{
    double max_abs_curl_grad(FieldBlock &node, int node_comp = 0);

    double max_abs_curl_grad(Field &fields,
                             const std::string &node_name,
                             int node_comp = 0);

    double max_abs_div_curl(FieldBlock &edge_xi,
                            FieldBlock &edge_eta,
                            FieldBlock &edge_zeta,
                            int edge_comp = 0);

    double max_abs_div_curl(Field &fields,
                            const EdgeFormNames &edge_names,
                            int edge_comp = 0);
}
