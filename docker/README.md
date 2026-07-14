# docker/

CPU-only images. AWS deployment (later) builds on the `runtime` stage.

Setup inside the image uses the canonical `uv sync --all-packages`, which installs
`qscat` and builds the Rust `qscat_kernels` in one step.

```bash
# Run the test stage (fails the build if tests fail):
docker build --target test -f docker/Dockerfile .

# Build and run the runtime image:
docker build --target runtime -t qmodeling:latest -f docker/Dockerfile .
docker run --rm qmodeling:latest
```
