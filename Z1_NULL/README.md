# Z1_NULL

`Z1_NULL` is the minimal topology construction example. It reads a case,
preprocesses the grid, builds `TOPO::Topology` including entity equivalence
and owner information, and then releases its local objects and exits.

It deliberately does not build field or halo objects and does not solve any
PDE. For framework tests, use `Z0_CoreDebug`; for Mercury use `Z4_Mercury`.
