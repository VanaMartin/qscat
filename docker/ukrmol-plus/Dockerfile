# syntax=docker/dockerfile:1
# Pinned external-engine toolchain/Psi4 layer; pilot retains reference binaries.
FROM zdenekmasin/ukrmol_plus@sha256:1cef2ee6aeaca2a4d5b813ff31e7d28bef4436a4282398c8933a1c047a8ec76c AS pilot

USER root
RUN apt-get update \
    && apt-get install -y --no-install-recommends perl numdiff ca-certificates \
    && rm -rf /var/lib/apt/lists/* \
    && /opt/ukrmolp/bin/python3 -m pip install --no-cache-dir --no-deps \
        numpy==2.3.4 scipy==1.16.2
ENV SSL_CERT_FILE=/etc/ssl/certs/ca-certificates.crt
RUN test -s "$SSL_CERT_FILE" \
    && perl -MStorable -MTime::HiRes -e 'print "Perl runtime ready\n"' \
    && /opt/ukrmolp/entrypoint.sh bash -c \
        'eval "$(psi4 --psiapi)"; python3 -c "import psi4, numpy, scipy; print(psi4.__version__)"'

USER ukrmol
WORKDIR /work
ENTRYPOINT ["/opt/ukrmolp/entrypoint.sh"]
CMD ["bash"]

# UKRmol-in includes its compatible GBTOlib sources; no independent GBTO tag.
FROM pilot AS build
USER root
RUN apt-get update \
    && apt-get install -y --no-install-recommends curl ca-certificates make binutils libc6-dev \
    && rm -rf /var/lib/apt/lists/*
USER ukrmol
ENV OMPI_FC=/opt/compiler/bin/gfortran OMPI_CC=/opt/compiler/bin/gcc
COPY docker/ukrmol-plus/source-build.sh /home/ukrmol/source-build.sh
ARG BUILD_JOBS=4
RUN /opt/ukrmolp/entrypoint.sh bash /home/ukrmol/source-build.sh "${BUILD_JOBS}"

FROM build AS test
# These upstream decks compare against distributed reference outputs via numdiff.
# A single worker bounds memory use; each parallel deck itself uses two ranks.
RUN --mount=type=tmpfs,target=/dev/shm,size=1073741824 \
    /opt/ukrmolp/entrypoint.sh bash -c \
    'export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1; \
     ctest --test-dir /home/ukrmol/source/build --output-on-failure -j1 \
       -R "C2v_target_HF.*_(serial|parallel)$" --no-tests=error'

FROM pilot AS runtime
USER root
# Remove the reference double engine before copying the source-built one.
RUN rm -rf /opt/ukrmolp/bin.double /opt/ukrmolp/lib.double
COPY --from=test /home/ukrmol/source/install/bin/ /opt/ukrmolp/bin.double/
COPY --from=test /home/ukrmol/source/install/lib/ /opt/ukrmolp/lib.double/
COPY --from=test /home/ukrmol/source/provenance.json /opt/ukrmolp/source-provenance.json
COPY --from=test /home/ukrmol/source/build/Testing/ /opt/ukrmolp/source-tests/
USER ukrmol

# Experiment layers follow the engine gate, so deck edits do not rebuild Fortran.
FROM pilot AS co-pilot
COPY projects/__init__.py /opt/qmodeling/projects/__init__.py
COPY projects/ukrmol_co /opt/qmodeling/projects/ukrmol_co
WORKDIR /opt/qmodeling
ENTRYPOINT ["/opt/ukrmolp/entrypoint.sh", "python3", "-m", "projects.ukrmol_co.run"]
CMD ["--workdir", "/work/runs/co-eq-sep20"]

FROM runtime AS co-source
COPY projects/__init__.py /opt/qmodeling/projects/__init__.py
COPY projects/ukrmol_co /opt/qmodeling/projects/ukrmol_co
WORKDIR /opt/qmodeling
ENTRYPOINT ["/opt/ukrmolp/entrypoint.sh", "python3", "-m", "projects.ukrmol_co.run"]
CMD ["--workdir", "/work/runs/co-eq-sep20"]
