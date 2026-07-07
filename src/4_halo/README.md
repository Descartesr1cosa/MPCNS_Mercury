# Halo Module Layout

`Halo` is organized by synchronization family first, then by implementation
detail:

- `core/`: registration, sync dispatch, debug helpers, and legacy facade calls.
- `build/`: halo region construction for 1D face, 2D edge, and 3D vertex layers.
- `component_copy/`: same-field copy exchange used by Cell, Node, and ordinary
  vector fields.
- `forms/`: orientation-aware triplet exchange for edge 1-forms and face
  2-forms.
- `coupling/`: coupling-buffer transfer and form-aware coupling helpers.
- `owner/`: owner-to-alias synchronization helpers used by equivalent
  Node/Edge/Face entities.

The runtime registry classifies each `FieldHaloRequest` into one of three sync
families:

- component copy: Cell/Node scalar, Cartesian vector, conserved vector, and any
  non-orientation-aware field.
- edge 1-form triplet: `E_xi/E_eta/E_zeta` style covariant edge forms.
- face 2-form triplet: `B_xi/B_eta/B_zeta` style contravariant face forms.

Layer names map to geometric halo depth:

- `FaceOnly`: 1D face halo.
- `Edge`: 2D edge/corner halo, including `FaceOnly`.
- `Vertex`: 3D vertex/corner halo, including `FaceOnly` and `Edge`.
