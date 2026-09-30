# SGS comparison

This setup compares the SGS closures that are active in the present CryosLES
source:

- `static` — `MODEL=1`, static Smagorinsky (`C_s=0.16`);
- `dynamic` — `MODEL=2`, plane-averaged dynamic Smagorinsky;
- `lagrangian_scale_dependent` — `MODEL=5`, Lagrangian scale-dependent dynamic Smagorinsky.

The cases share the same baseline (`../20_PARAMETERS.py_0`) and differ only in
the SGS selector. `prepare_sgs_case.sh` stages the selected runtime inputs:

```bash
cd model
./sgs_comparison/prepare_sgs_case.sh static
```

It writes the runtime-only `PARAMETERS.py_0`, `patches.in`, and
`hetero_srf_props.in`. The staged setup uses a `4*pi x 4*pi x 4*pi` m domain,
the compiled `32 x 32 x 32` grid, and a homogeneous surface with
`z0 = 0.01 m`. Compile and launch with 16 MPI ranks after changing this branch.


At 20 Hz, `inst_slice` writes only the vertical centreline to
`output_0/ts1_centerline/`: `u_Z.c<rank>`, `v_Z.c<rank>`, `w_Z.c<rank>`, and
`t_Z.c<rank>`. Each direct-access record contains the local `NZ-1` levels; for
this 32³ setup, concatenate rank files `c0` through `c15` in rank order to
obtain the full 32-level profile at each sample.
`MODEL=3` and `MODEL=4` remain documented in an old input-file comment, but
they have no active closure call in this revision, so they are intentionally
excluded.
