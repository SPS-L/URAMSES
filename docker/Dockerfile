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
CMD make -f Makefile.linux dll \
    && mkdir -p /output \
    && cp Release_gnu_l/ramses.so /output/ramses.so \
    && echo "✅  ramses.so written to /output/ramses.so"
