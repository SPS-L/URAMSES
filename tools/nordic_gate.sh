#!/usr/bin/env bash
# Nordic voltage-collapse regression gate for CI.
#
#   nordic_gate.sh <dynsim-binary> <repo-root> <baseline-npz>
#
# Ported from stepss-ramses, where the same script gates the `ramses` binary
# before a release is published. Here it gates `dynsim`, which is the same
# simulator linked against a pre-compiled kit plus this repository's model
# routers, so it proves a refreshed kit still reproduces RAMSES trajectory
# for trajectory rather than merely linking. The baseline is RAMSES' own:
# the Nordic case does not exercise the example models in custom_models/, so
# a correct dynsim must match it exactly.
#
# Keep in step with stepss-ramses/tools/nordic_gate.sh; the two are meant to
# apply the same assertions to the same case.
#
# Runs the Nordic case (dyn_A + volt_rat_A + short_trip_branch.dst) in a
# temp dir and asserts:
#   1. exit code is exactly 255 (Linux) or 127 (Windows/MSYS2) for the
#      by-design sim_minmaxvolt trip; never 0. The trip path exits via
#      exit(-1); glibc masks it to 255 while MSYS2 maps the 32-bit Windows
#      code to 127. A genuine launch failure can also read 127, but is
#      caught by the log-marker and comparator assertions that follow.
#   2. the log contains '**Simulation finished**' and a nonzero time-step
#      count,
#   3. tools/compare_trj.py accepts obs.trj against the baseline.
set -u
if [ $# -ne 3 ]; then
    echo "usage: $0 <ramses-binary> <repo-root> <baseline-npz>" >&2
    exit 2
fi
BIN="$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
ROOT="$(cd "$2" && pwd)"
BASE="$(cd "$(dirname "$3")" && pwd)/$(basename "$3")"
PY="$(command -v python3 || command -v python)"
[ -x "$BIN" ] || { echo "FAIL: ramses binary not found: $BIN"; exit 1; }
[ -f "$BASE" ] || { echo "FAIL: baseline not found: $BASE"; exit 1; }
[ -n "$PY" ] || { echo "FAIL: no python interpreter on PATH"; exit 1; }

RUN="$(mktemp -d)"
trap 'rm -rf "$RUN"' EXIT
cp "$ROOT"/examples/Nordic/{dyn_A.dat,volt_rat_A.dat,settings1.dat,obs.dat,short_trip_branch.dst,cmd_ci.txt} "$RUN/" \
    || { echo "FAIL: could not copy Nordic inputs from $ROOT/examples/Nordic"; exit 1; }
cd "$RUN" || exit 1

"$BIN" -t cmd_ci.txt > run.log 2>&1
rc=$?
# The by-design sim_minmaxvolt trip exits via exit(-1) (src/main.f90 ->
# ramses() returning diagno=-1): glibc masks that to 255, while MSYS2 bash
# maps the Windows exit code 0xFFFFFFFF to 127. A genuine launch failure
# can also yield 127, but then run.log lacks the markers checked next.
if [ "$rc" -ne 255 ] && [ "$rc" -ne 127 ]; then
    echo "FAIL: expected exit code 255 (Linux) or 127 (MSYS2) for the by-design collapse trip, got $rc"
    tail -40 run.log
    exit 1
fi
if ! grep -qF '**Simulation finished**' run.log; then
    echo "FAIL: '**Simulation finished**' marker missing from log"
    tail -40 run.log
    exit 1
fi
steps="$(awk '/^Time steps/{sub(/\r$/,""); print $3}' run.log)"
if [ -z "$steps" ] || [ "$steps" -le 0 ]; then
    echo "FAIL: zero or missing time-step count"
    tail -40 run.log
    exit 1
fi
echo "run OK: exit $rc, $steps time steps"
"$PY" "$ROOT/tools/compare_trj.py" compare obs.trj "$BASE"
