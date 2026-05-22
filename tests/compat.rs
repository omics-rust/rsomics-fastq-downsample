use std::path::PathBuf;
use std::process::{Command, Stdio};

fn ours() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-fastq-downsample"))
}

fn fixture() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/golden/small.fq")
}

fn seqtk_available() -> bool {
    Command::new("seqtk")
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .is_ok()
}

fn seqkit_available() -> bool {
    Command::new("seqkit")
        .arg("version")
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .map(|s| s.success())
        .unwrap_or(false)
}

fn count_reads(fq: &str) -> usize {
    fq.lines().filter(|l| l.starts_with('@')).count()
}

#[test]
fn fraction_1_keeps_all() {
    let out = Command::new(ours())
        .arg(fixture())
        .args(["-f", "1.0", "--seed", "42"])
        .output()
        .unwrap();
    assert!(out.status.success());
    let s = String::from_utf8_lossy(&out.stdout);
    assert_eq!(count_reads(&s), 4, "fraction=1.0 should keep all 4 reads");
}

#[test]
fn fraction_0_keeps_none() {
    let out = Command::new(ours())
        .arg(fixture())
        .args(["-f", "0.0", "--seed", "42"])
        .output()
        .unwrap();
    assert!(out.status.success());
    let s = String::from_utf8_lossy(&out.stdout);
    assert_eq!(count_reads(&s), 0, "fraction=0.0 should keep no reads");
}

#[test]
fn deterministic() {
    let run = |seed: &str| -> String {
        let out = Command::new(ours())
            .arg(fixture())
            .args(["-f", "0.5", "--seed", seed])
            .output()
            .unwrap();
        assert!(out.status.success());
        String::from_utf8(out.stdout).unwrap()
    };
    let r1 = run("99");
    let r2 = run("99");
    assert_eq!(r1, r2, "same seed must produce identical output");
}

#[test]
fn output_is_valid_fastq() {
    let out = Command::new(ours())
        .arg(fixture())
        .args(["-f", "1.0", "--seed", "1"])
        .output()
        .unwrap();
    assert!(out.status.success());
    let s = String::from_utf8_lossy(&out.stdout);
    let lines: Vec<&str> = s.lines().collect();
    assert_eq!(lines.len() % 4, 0, "FASTQ must have lines divisible by 4");
    for chunk in lines.chunks(4) {
        assert!(chunk[0].starts_with('@'), "header line must start with @");
        assert_eq!(chunk[2], "+", "separator line must be +");
        assert_eq!(
            chunk[1].len(),
            chunk[3].len(),
            "seq and qual must have equal length"
        );
    }
}

#[test]
fn count_comparable_to_seqtk() {
    if !seqtk_available() {
        eprintln!("skipping: seqtk not found");
        return;
    }
    let ours_out = Command::new(ours())
        .arg(fixture())
        .args(["-f", "0.5", "--seed", "42"])
        .output()
        .unwrap();
    assert!(ours_out.status.success());
    let our_count = count_reads(&String::from_utf8_lossy(&ours_out.stdout));

    let seqtk_out = Command::new("seqtk")
        .args(["sample", "-s", "42", fixture().to_str().unwrap(), "2"])
        .output()
        .unwrap();
    let seqtk_count = count_reads(&String::from_utf8_lossy(&seqtk_out.stdout));

    assert!(our_count <= 4 && our_count > 0);
    assert!(seqtk_count <= 4 && seqtk_count > 0);
}

// Rate-based compat vs `seqkit sample -p` (the perfgate upstream): random
// samplers can't match exact reads across RNGs, so verify both keep a
// comparable fraction. Uses a larger generated fixture for a stable rate.
#[test]
fn rate_comparable_to_seqkit() {
    if !seqkit_available() {
        eprintln!("skipping: seqkit not found");
        return;
    }
    let dir = std::env::temp_dir().join("rsomics-downsample-seqkit");
    std::fs::create_dir_all(&dir).unwrap();
    let big = dir.join("reads.fq");
    let mut s = String::new();
    for i in 0..2000 {
        s.push_str(&format!("@r{i}\nACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIII\n"));
    }
    std::fs::write(&big, &s).unwrap();

    let our_count = count_reads(&String::from_utf8_lossy(
        &Command::new(ours())
            .arg(&big)
            .args(["-f", "0.3", "--seed", "42"])
            .output()
            .unwrap()
            .stdout,
    ));
    let sk_count = count_reads(&String::from_utf8_lossy(
        &Command::new("seqkit")
            .args(["sample", "-p", "0.3", "-s", "42"])
            .arg(&big)
            .output()
            .unwrap()
            .stdout,
    ));
    let our_frac = our_count as f64 / 2000.0;
    let sk_frac = sk_count as f64 / 2000.0;
    assert!(
        (our_frac - 0.3).abs() < 0.05,
        "ours kept {our_frac}, want ~0.3"
    );
    assert!(
        (our_frac - sk_frac).abs() < 0.05,
        "ours {our_frac} vs seqkit {sk_frac} rate mismatch"
    );
}
