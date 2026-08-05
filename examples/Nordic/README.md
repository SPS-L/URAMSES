# Nordic test system (regression case)

The Nordic voltage-collapse case, used by `tools/nordic_gate.sh` to prove that
a `dynsim` built here still reproduces RAMSES trajectory for trajectory.

These files are a copy of the same case in `stepss-ramses`, which is where the
baseline in `tests/baselines/nordic_baseline.npz` was generated. The case does
not exercise the example models in `custom_models/`, so a correct `dynsim` must
match that baseline exactly; any difference means the module kit, the model
routers or the build changed the simulator's behaviour.

| File | Role |
|---|---|
| `dyn_A.dat` | Network, machines and controls (77 buses, 23 generators) |
| `volt_rat_A.dat` | Voltage ratings |
| `settings1.dat` | Solver settings |
| `obs.dat` | Observable definitions (1417 columns) |
| `short_trip_branch.dst` | Disturbance: branch trip leading to voltage collapse |
| `cmd_ci.txt` | Command file wiring the above together |

Run it by hand with:

```sh
make -f build/Makefile.linux all
bash tools/nordic_gate.sh Release_l/dynsim . tests/baselines/nordic_baseline.npz
```

The run is expected to exit 255 (127 under MSYS2): the case ends in a
by-design `sim_minmaxvolt` trip, not a clean stop. The gate checks for that
exit code, for the `**Simulation finished**` marker, and then compares the
trajectory.

Keep in step with `stepss-ramses/examples/Nordic/`; the two are meant to be
the same case. Licensed under the Apache License 2.0, see `LICENSE`.
