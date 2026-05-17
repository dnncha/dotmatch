use std::fs;

use dotmatch_workbench::{
    build_dotmatch_command, canonical_workspace, resolve_workspace_path, DotmatchAction,
};

#[test]
fn resolve_workspace_path_rejects_absolute_paths_and_parent_traversal() {
    let tmp = tempfile::tempdir().expect("tempdir");
    let workspace = canonical_workspace(tmp.path()).expect("workspace");

    assert!(resolve_workspace_path(&workspace, "/tmp/reads.fastq").is_err());
    assert!(resolve_workspace_path(&workspace, "../reads.fastq").is_err());
}

#[test]
fn resolve_workspace_path_rejects_symlink_escape() {
    let tmp = tempfile::tempdir().expect("tempdir");
    let outside = tempfile::tempdir().expect("outside");
    let workspace = canonical_workspace(tmp.path()).expect("workspace");
    let link = tmp.path().join("outside-link");

    #[cfg(unix)]
    std::os::unix::fs::symlink(outside.path(), &link).expect("symlink");

    #[cfg(unix)]
    assert!(resolve_workspace_path(&workspace, "outside-link/file.fastq").is_err());
}

#[test]
fn resolve_workspace_path_allows_missing_files_when_parent_is_inside_workspace() {
    let tmp = tempfile::tempdir().expect("tempdir");
    fs::create_dir(tmp.path().join("inputs")).expect("inputs");
    let workspace = canonical_workspace(tmp.path()).expect("workspace");

    let path = resolve_workspace_path(&workspace, "inputs/new_assay.toml").expect("path");

    assert!(path.ends_with("inputs/new_assay.toml"));
}

#[test]
fn build_dotmatch_command_uses_explicit_argv_for_allowed_actions() {
    let argv = build_dotmatch_command(
        "/usr/local/bin/dotmatch",
        DotmatchAction::AssayPlan,
        vec!["assay.toml".to_string()],
    )
    .expect("argv");

    assert_eq!(argv, vec!["/usr/local/bin/dotmatch", "assay", "plan", "assay.toml"]);
}

#[test]
fn build_dotmatch_command_rejects_shell_metacharacters_in_extra_args() {
    let err = build_dotmatch_command(
        "/usr/local/bin/dotmatch",
        DotmatchAction::AssayCheck,
        vec!["assay.toml; rm -rf /".to_string()],
    )
    .expect_err("rejected");

    assert!(err.contains("unsafe argument"));
}
