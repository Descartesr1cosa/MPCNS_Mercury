# Z1_NULL

`Z1_NULL` is a full empty physics template. It demonstrates the standard
structure for adding a new solver module, including field registration,
coupling/halo setup, boundary hooks, initialization, diagnostics, output hooks,
and a no-op time loop.

It does not solve any PDE. For low-level framework tests, use `Z0_CoreDebug`.
For real Mercury MHD/Hall-MHD, use `Z4_Mercury`.
