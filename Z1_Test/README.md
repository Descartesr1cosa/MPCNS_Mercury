# Z1_Test

`Z1_Test` collects small framework tests.

Current tests:

- `topo_test`: reads a case and grid, builds topology, then checks:
  1. local `boundary(boundary(face)) == 0`
  2. local `boundary(boundary(cell)) == 0`
  3. global `D1 * D0 == 0`
  4. global `D2 * D1 == 0`
  5. owner/alias uniqueness and sign consistency
