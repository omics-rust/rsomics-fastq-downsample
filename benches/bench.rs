use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;
use std::path::PathBuf;
use std::process::Command;

fn bench_fastq_downsample(c: &mut Criterion) {
    let bin = env!("CARGO_BIN_EXE_rsomics-fastq-downsample");
    let manifest = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let fq = manifest.join("tests/golden/small.fq");
    c.bench_function("rsomics-fastq-downsample golden", |b| {
        b.iter(|| {
            let out = Command::new(black_box(bin))
                .args([fq.to_str().unwrap(), "-f", "0.5", "--seed", "42"])
                .output()
                .unwrap();
            assert!(out.status.success());
        });
    });
}

criterion_group!(benches, bench_fastq_downsample);
criterion_main!(benches);
