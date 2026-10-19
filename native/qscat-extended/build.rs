fn main() {
    println!("cargo:rerun-if-changed=src/extended.c");
    cc::Build::new()
        .file("src/extended.c")
        .opt_level(3)
        .flag_if_supported("-fno-fast-math")
        .flag_if_supported("-ffp-contract=off")
        .warnings(true)
        .compile("qscat_extended");
}
