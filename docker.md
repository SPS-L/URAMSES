## Dockerfile

This route always builds the **Linux** artifact (`ramses.so` via
`Makefile.gfortran` and `modules_lin/`), whatever the host. It is the
recommended option on Intel Macs, for which no pre-compiled RAMSES kit is
published. For a native build, use `Makefile.macos` (Apple Silicon) or
`Makefile.mingw` (Windows/MSYS2) instead.

Place this at the root of your `stepss-uramses` clone (alongside `Makefile.gfortran`):

```dockerfile
# docker/Dockerfile
FROM ubuntu:24.04

LABEL maintainer="SPS-L <info@sps-lab.org>"
LABEL description="stepss-uramses build environment – compiles ramses.so for stepss-pyramses"

# ── System deps ────────────────────────────────────────────────────────────────
RUN apt-get update && apt-get install -y --no-install-recommends \
        gfortran \
        libopenblas-dev \
        make \
        git \
    && rm -rf /var/lib/apt/lists/*

# ── Working directory (mapped from host at run time) ──────────────────────────
WORKDIR /uramses

# ── Default: build ramses.so and copy it to /output ──────────────────────────
CMD make -f Makefile.gfortran dll \
    && mkdir -p /output \
    && cp Release_gnu_l/ramses.so /output/ramses.so \
    && echo "✅  ramses.so written to /output/ramses.so"
```


***

## docker-compose.yml

For repeatable one-command usage:

```yaml
services:
  uramses-build:
    build:
      context: .
      dockerfile: docker/Dockerfile
    volumes:
      # repo root → container working dir (edits to my_models/ visible immediately)
      - .:/uramses
      # output dir on host receives the compiled ramses.so
      - ./output:/output
```


***

## Workflow for the User

### One-time setup

```bash
git clone https://github.com/SPS-L/stepss-uramses.git
cd stepss-uramses
docker compose build          # ~2 min, only once
```


### Add / edit a model → rebuild

```bash
# Edit or drop new .f90 files into my_models/
nano my_models/my_new_model.f90

# Rebuild ramses.so
docker compose run --rm uramses-build
# → output/ramses.so is updated
```


### Use ramses.so with stepss-pyramses

```bash
# Point pyramses to the freshly built .so
cp output/ramses.so /path/to/stepss-pyramses/lib/ramses.so
```

Or set the path via the `RAMSES_SO` env var if stepss-pyramses supports it.

***

## Key Design Decisions

| Choice | Rationale |
| :-- | :-- |
| `ubuntu:24.04` base | Matches your target; `gfortran-13` ships in noble |
| `libopenblas-dev` | Satisfies `-lopenblas` in `Makefile.gfortran` |
| Bind-mount `.:/uramses` | User edits `my_models/` on host, no rebuild of the image needed |
| `/output` volume | Clean separation; `ramses.so` lands on the host automatically |
| `make dll` only | Only builds `ramses.so`, skipping `dynsim` – that's all pyramses needs |
| No `COPY` in `Dockerfile` | Image stays generic; works for any clone of the repo |


***

## Optional: `build.sh` Helper Script

For users who don't want to think about Docker commands:

```bash
#!/usr/bin/env bash
# build.sh – rebuild ramses.so inside Docker
set -e
docker compose run --rm uramses-build
echo "ramses.so is at: $(pwd)/output/ramses.so"
```

```bash
chmod +x build.sh
./build.sh
```


***

## Adding a New Model (end-to-end)

1. Create `my_models/inj_MyDevice.f90` implementing the RAMSES injector interface.
2. Register it in `src/usr_inj_models.f90` (add `USE inj_MyDevice` and the dispatch case).
3. Run `./build.sh` → `output/ramses.so` is rebuilt with the new model linked in.
4. Copy `ramses.so` to your `stepss-pyramses` environment and run your Python simulations normally.

No compiler installation, no Intel oneAPI setup, no Visual Studio – just Docker.
<span style="display:none">[^1][^2]</span>

<div align="center">⁂</div>

[^1]: STEPSS-Repositories.md

[^2]: propose-the-best-framework-pla-k8Y5ppWKS5KWp1sU.lYQRA.md

