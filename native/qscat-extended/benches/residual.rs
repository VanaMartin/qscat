use criterion::{Criterion, criterion_group, criterion_main};
use qscat_extended_kernels::{benchmark_fixture, capability, checked_residual};
use std::hint::black_box;
use std::time::Duration;

fn residual(c: &mut Criterion) {
    if !capability().3 {
        eprintln!("Native long-double benchmark unsupported on this host");
        return;
    }
    for n in [64, 298] {
        let a = benchmark_fixture(n * n);
        let b = benchmark_fixture(n * n);
        let x = benchmark_fixture(n * n);
        c.bench_function(&format!("extended_residual_{n}"), |bench| {
            bench.iter(|| black_box(checked_residual(&a, &b, &x, n, n, n).unwrap()))
        });
    }
}
criterion_group! { name = benches; config = Criterion::default().sample_size(10)
.warm_up_time(Duration::from_secs(1)).measurement_time(Duration::from_secs(2)); targets = residual }
criterion_main!(benches);
