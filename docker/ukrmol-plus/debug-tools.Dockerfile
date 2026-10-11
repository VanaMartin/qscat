# Diagnostic tools over the exact successful CAS12 runtime/library image.
FROM sha256:6b3b0ededa85494076229b111efe969985407023e404bff888f17cafb4630666
USER root
RUN apt-get update \
    && apt-get install -y --no-install-recommends gdb binutils make libc6-dev \
    && rm -rf /var/lib/apt/lists/*
USER ukrmol
ENTRYPOINT ["/opt/ukrmolp/entrypoint.sh"]
