# Z1_Test

`Z1_Test` collects small framework tests. Each subdirectory owns one focused
test executable and its own `CMakeLists.txt`.

Current tests:

- `topo_test`: reads a case and grid, builds topology, then checks local and
  global topological boundary closure.
