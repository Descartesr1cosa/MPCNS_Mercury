#pragma once

class Field;

namespace mpcns::compat
{

// Legacy Field-backed reconstruction kernels. They return false when their
// required staggered samples or cached reconstruction data are unavailable.
bool reconstruct_legacy_face_2form_to_cell(Field &fields, int fid_xi, int fid_eta, int fid_zeta, int iblock, int i, int j, int k, double out_xyz[3]);
bool reconstruct_legacy_edge_1form_to_cell(Field &fields, int fid_xi, int fid_eta, int fid_zeta, int iblock, int i, int j, int k, double out_xyz[3]);
bool reconstruct_legacy_face_2form_to_node(Field &fields, int fid_xi, int fid_eta, int fid_zeta, int iblock, int i, int j, int k, double out_xyz[3]);
bool reconstruct_legacy_edge_1form_to_node(Field &fields, int fid_xi, int fid_eta, int fid_zeta, int iblock, int i, int j, int k, double out_xyz[3]);

} // namespace mpcns::compat
