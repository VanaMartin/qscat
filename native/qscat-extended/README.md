# Experimental extended residual kernel

Rust/PyO3 checks matrix buffer dimensions and a small C shim evaluates long-double
dot products. The Python array interface and retained NumPy reference live in
`projects/ukrmol_co/numerov_native_residual.py`; the extension's byte interface is
an internal transport ABI. No cast to binary64 occurs.

Qualified native format: little-endian x86-64,16-byte storage,64-significand-bit
longdouble, maximum exponent16384. Other hosts can build/import the capability
probe and use the NumPy reference; native arithmetic rejects unsupported formats.
The C shim is compiled without reassociation or multiply/add contraction.

```sh
uv run --no-sync maturin develop --release --manifest-path native/qscat-extended/Cargo.toml
uv run --no-sync pytest -q projects/ukrmol_co/test_numerov_native_residual.py
cargo test --locked --manifest-path native/qscat-extended/Cargo.toml
cargo clippy --locked --manifest-path native/qscat-extended/Cargo.toml --all-targets -- -D warnings
cargo bench --locked --manifest-path native/qscat-extended/Cargo.toml --bench residual
```

If both Conda and a virtualenv are active, unset `CONDA_PREFIX` for maturin.
The published qualification wheel uses abi3-py311 to run in the pinned external
tools container's Python3.11 as well as the layered build image's Python3.12.
The repository's ordinary Python requirement remains>=3.12. Wheel extraction in
the external tools image is a declared numerical experiment boundary.

[Design](../../docs/physics/co-numerov-native-residual.md) ·
[Qualification](../../projects/ukrmol_co/NUMEROV_NATIVE_RESIDUAL.md).
