# Z1_Test

`Z1_Test` collects small framework tests. Each subdirectory owns one focused
test executable and its own `CMakeLists.txt`.

Current tests:

- `topo_test`: reads a case and grid, builds topology, then checks local and
  global topological boundary closure.
- `halo_form_test`: owner/alias sign checks for edge 1-forms and face 2-forms,
  with interface consistency placeholders.
- `dec_operator_test`: DEC exterior derivative smoke checks, with
  codifferential adjointness placeholder.
- `metric_hodge_test`: metric/Hodge finite and positivity diagnostics.
