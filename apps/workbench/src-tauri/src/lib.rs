use serde::{Deserialize, Serialize};
use std::env;
use std::ffi::OsString;
use std::fs;
use std::path::{Component, Path, PathBuf};
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};

#[derive(Clone, Debug, Deserialize, Serialize)]
pub enum DotmatchAction {
    Version,
    Dist,
    AssayInfer,
    AssayCheck,
    AssayPlan,
    AssayRun,
    AssayAutopsy,
}

#[derive(Debug, Deserialize)]
pub struct WorkbenchCommandRequest {
    pub workspace: String,
    pub dotmatch_path: Option<String>,
    pub action: DotmatchAction,
    pub args: Vec<String>,
    pub log_name: Option<String>,
}

#[derive(Debug, Serialize)]
pub struct CommandResult {
    pub argv: Vec<String>,
    pub exit_code: i32,
    pub stdout: String,
    pub stderr: String,
    pub log_path: String,
}

#[derive(Debug, Serialize)]
pub struct DoctorReport {
    pub dotmatch_path: String,
    pub checks: Vec<CommandResult>,
}

pub fn canonical_workspace<P: AsRef<Path>>(workspace: P) -> Result<PathBuf, String> {
    let path = workspace
        .as_ref()
        .canonicalize()
        .map_err(|err| format!("workspace does not exist: {err}"))?;
    if !path.is_dir() {
        return Err("workspace must be a directory".to_string());
    }
    Ok(path)
}

pub fn resolve_workspace_path(workspace: &Path, relative: &str) -> Result<PathBuf, String> {
    if relative.trim().is_empty() {
        return Err("path must not be empty".to_string());
    }
    let rel_path = Path::new(relative);
    if rel_path.is_absolute() {
        return Err("absolute paths are not allowed; use a path relative to the workspace".to_string());
    }
    if rel_path
        .components()
        .any(|component| matches!(component, Component::ParentDir | Component::RootDir | Component::Prefix(_)))
    {
        return Err("paths must stay inside the workspace".to_string());
    }

    let candidate = workspace.join(rel_path);
    let existing_anchor = existing_anchor(&candidate)?;
    let canonical_anchor = existing_anchor
        .canonicalize()
        .map_err(|err| format!("could not canonicalize path: {err}"))?;

    if !canonical_anchor.starts_with(workspace) {
        return Err("path resolves outside the workspace".to_string());
    }

    if candidate.exists() {
        return candidate
            .canonicalize()
            .map_err(|err| format!("could not canonicalize path: {err}"));
    }

    let suffix = candidate
        .strip_prefix(&existing_anchor)
        .map_err(|_| "could not resolve workspace-relative path".to_string())?;
    Ok(canonical_anchor.join(suffix))
}

pub fn build_dotmatch_command(
    dotmatch_path: &str,
    action: DotmatchAction,
    args: Vec<String>,
) -> Result<Vec<String>, String> {
    for arg in &args {
        validate_safe_arg(arg)?;
    }

    let mut argv = vec![dotmatch_path.to_string()];
    match action {
        DotmatchAction::Version => {
            if !args.is_empty() {
                return Err("version check does not accept extra arguments".to_string());
            }
            argv.push("--version".to_string());
        }
        DotmatchAction::Dist => {
            if args.is_empty() {
                argv.extend(["dist".to_string(), "ACGT".to_string(), "AGGT".to_string()]);
            } else {
                argv.push("dist".to_string());
                argv.extend(args);
            }
        }
        DotmatchAction::AssayInfer => {
            argv.extend(["assay".to_string(), "infer".to_string()]);
            argv.extend(args);
        }
        DotmatchAction::AssayCheck => {
            argv.extend(["assay".to_string(), "check".to_string()]);
            argv.extend(args);
        }
        DotmatchAction::AssayPlan => {
            argv.extend(["assay".to_string(), "plan".to_string()]);
            argv.extend(args);
        }
        DotmatchAction::AssayRun => {
            argv.extend(["assay".to_string(), "run".to_string()]);
            argv.extend(args);
        }
        DotmatchAction::AssayAutopsy => {
            argv.extend(["assay".to_string(), "autopsy".to_string()]);
            argv.extend(args);
        }
    }
    Ok(argv)
}

#[tauri::command]
pub fn doctor(dotmatch_path: Option<String>) -> Result<DoctorReport, String> {
    let dotmatch = find_dotmatch(dotmatch_path)?;
    let checks = vec![
        execute_argv(
            build_dotmatch_command(&dotmatch, DotmatchAction::Version, vec![])?,
            None,
            None,
        )?,
        execute_argv(
            build_dotmatch_command(&dotmatch, DotmatchAction::Dist, vec![])?,
            None,
            None,
        )?,
    ];
    Ok(DoctorReport {
        dotmatch_path: dotmatch,
        checks,
    })
}

#[tauri::command]
pub fn run_workbench_command(request: WorkbenchCommandRequest) -> Result<CommandResult, String> {
    let workspace = canonical_workspace(&request.workspace)?;
    validate_workspace_args(&workspace, &request.action, &request.args)?;
    let dotmatch = find_dotmatch(request.dotmatch_path)?;
    let argv = build_dotmatch_command(&dotmatch, request.action, request.args)?;
    let log_path = log_path(&workspace, request.log_name.as_deref())?;
    execute_argv(argv, Some(&workspace), Some(&log_path))
}

#[tauri::command]
pub fn read_text_artifact(workspace: String, relative_path: String) -> Result<String, String> {
    let workspace = canonical_workspace(workspace)?;
    let path = resolve_workspace_path(&workspace, &relative_path)?;
    fs::read_to_string(path).map_err(|err| format!("could not read artifact: {err}"))
}

#[tauri::command]
pub fn write_text_artifact(workspace: String, relative_path: String, contents: String) -> Result<String, String> {
    let workspace = canonical_workspace(workspace)?;
    let path = resolve_workspace_path(&workspace, &relative_path)?;
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent).map_err(|err| format!("could not create output directory: {err}"))?;
    }
    fs::write(&path, contents).map_err(|err| format!("could not write artifact: {err}"))?;
    display_workspace_path(&workspace, &path)
}

pub fn run() {
    tauri::Builder::default()
        .plugin(tauri_plugin_dialog::init())
        .plugin(tauri_plugin_opener::init())
        .invoke_handler(tauri::generate_handler![
            doctor,
            run_workbench_command,
            read_text_artifact,
            write_text_artifact
        ])
        .run(tauri::generate_context!())
        .expect("failed to run DotMatch Workbench");
}

fn existing_anchor(candidate: &Path) -> Result<PathBuf, String> {
    let mut current = candidate.to_path_buf();
    while !current.exists() {
        if !current.pop() {
            return Err("workspace path has no existing parent".to_string());
        }
    }
    Ok(current)
}

fn validate_safe_arg(arg: &str) -> Result<(), String> {
    if arg.chars().any(|ch| ch.is_control() || matches!(ch, ';' | '&' | '|' | '<' | '>' | '`' | '$')) {
        return Err(format!("unsafe argument rejected: {arg}"));
    }
    Ok(())
}

fn validate_workspace_args(workspace: &Path, action: &DotmatchAction, args: &[String]) -> Result<(), String> {
    let path_flags = [
        "--targets",
        "--barcodes",
        "--left-targets",
        "--right-targets",
        "--reads",
        "--out",
        "--report",
        "--candidates",
        "--out-dir",
    ];
    let mut i = 0;
    while i < args.len() {
        if path_flags.contains(&args[i].as_str()) {
            let value = args.get(i + 1).ok_or_else(|| format!("{} requires a value", args[i]))?;
            resolve_workspace_path(workspace, value)?;
            i += 2;
            continue;
        }
        i += 1;
    }

    match action {
        DotmatchAction::AssayCheck | DotmatchAction::AssayPlan | DotmatchAction::AssayRun => {
            let spec = args.first().ok_or_else(|| "assay command requires a spec path".to_string())?;
            resolve_workspace_path(workspace, spec)?;
        }
        DotmatchAction::AssayAutopsy => {
            let spec = args.first().ok_or_else(|| "autopsy requires a spec path".to_string())?;
            resolve_workspace_path(workspace, spec)?;
        }
        DotmatchAction::Version | DotmatchAction::Dist | DotmatchAction::AssayInfer => {}
    }
    Ok(())
}

fn find_dotmatch(configured: Option<String>) -> Result<String, String> {
    if let Some(path) = configured.filter(|value| !value.trim().is_empty()) {
        return executable_path(path);
    }
    if let Ok(path) = env::var("DOTMATCH_WORKBENCH_DOTMATCH") {
        if !path.trim().is_empty() {
            return executable_path(path);
        }
    }
    find_in_path("dotmatch").ok_or_else(|| {
        "could not find dotmatch; install DotMatch, add it to PATH, set DOTMATCH_WORKBENCH_DOTMATCH, or configure an executable path".to_string()
    })
}

fn executable_path(path: String) -> Result<String, String> {
    let candidate = PathBuf::from(path);
    if !candidate.exists() {
        return Err("configured dotmatch executable does not exist".to_string());
    }
    Ok(candidate.to_string_lossy().to_string())
}

fn find_in_path(binary: &str) -> Option<String> {
    let path_var = env::var_os("PATH")?;
    env::split_paths(&path_var).find_map(|dir| {
        let candidate = dir.join(binary);
        if candidate.exists() {
            Some(candidate.to_string_lossy().to_string())
        } else {
            None
        }
    })
}

fn execute_argv(argv: Vec<String>, cwd: Option<&Path>, log_path: Option<&Path>) -> Result<CommandResult, String> {
    let program = argv.first().ok_or_else(|| "empty command".to_string())?;
    let args: Vec<OsString> = argv.iter().skip(1).map(OsString::from).collect();
    let mut command = Command::new(program);
    command.args(args);
    if let Some(cwd) = cwd {
        command.current_dir(cwd);
    }
    let output = command.output().map_err(|err| format!("could not run dotmatch: {err}"))?;
    let result = CommandResult {
        argv,
        exit_code: output.status.code().unwrap_or(-1),
        stdout: String::from_utf8_lossy(&output.stdout).to_string(),
        stderr: String::from_utf8_lossy(&output.stderr).to_string(),
        log_path: log_path
            .map(|path| path.to_string_lossy().to_string())
            .unwrap_or_default(),
    };
    if let Some(path) = log_path {
        if let Some(parent) = path.parent() {
            fs::create_dir_all(parent).map_err(|err| format!("could not create log directory: {err}"))?;
        }
        let json = serde_json::to_string_pretty(&result).map_err(|err| format!("could not serialize log: {err}"))?;
        fs::write(path, json).map_err(|err| format!("could not write log: {err}"))?;
    }
    Ok(result)
}

fn log_path(workspace: &Path, label: Option<&str>) -> Result<PathBuf, String> {
    let timestamp = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_err(|err| format!("clock error: {err}"))?
        .as_secs();
    let clean_label = label
        .unwrap_or("dotmatch-command")
        .chars()
        .map(|ch| if ch.is_ascii_alphanumeric() || ch == '-' || ch == '_' { ch } else { '-' })
        .collect::<String>();
    Ok(workspace
        .join(".dotmatch")
        .join("workbench")
        .join("logs")
        .join(format!("{timestamp}-{clean_label}.json")))
}

fn display_workspace_path(workspace: &Path, path: &Path) -> Result<String, String> {
    path.strip_prefix(workspace)
        .map(|path| path.to_string_lossy().to_string())
        .map_err(|_| "path is outside workspace".to_string())
}
