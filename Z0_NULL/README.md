# Z0_NULL

`Z0_NULL` is the minimal topology construction example. It reads a case,
preprocesses the grid, builds `TOPO::Topology` including entity equivalence
and owner information, registers one cell-centered scalar field, constructs
the halo/boundary/solver skeleton, runs one no-op step, and exits normally.

It deliberately does not implement any PDE. Use it as the empty application
template before adding real field catalogs, boundary handlers, and solver
operators.
