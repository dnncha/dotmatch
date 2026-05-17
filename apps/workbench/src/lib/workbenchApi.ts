import { invoke } from "@tauri-apps/api/core";

export type DotmatchAction =
  | "Version"
  | "Dist"
  | "AssayInfer"
  | "AssayCheck"
  | "AssayPlan"
  | "AssayRun"
  | "AssayAutopsy";

export type CommandResult = {
  argv: string[];
  exit_code: number;
  stdout: string;
  stderr: string;
  log_path: string;
};

export type DoctorReport = {
  dotmatch_path: string;
  checks: CommandResult[];
};

export function doctor(dotmatchPath?: string): Promise<DoctorReport> {
  return invoke("doctor", { dotmatchPath });
}

export function runWorkbenchCommand(input: {
  workspace: string;
  dotmatchPath?: string;
  action: DotmatchAction;
  args: string[];
  logName?: string;
}): Promise<CommandResult> {
  return invoke("run_workbench_command", {
    request: {
      workspace: input.workspace,
      dotmatch_path: input.dotmatchPath ?? null,
      action: input.action,
      args: input.args,
      log_name: input.logName ?? null
    }
  });
}

export function readTextArtifact(workspace: string, relativePath: string): Promise<string> {
  return invoke("read_text_artifact", { workspace, relativePath });
}

export function writeTextArtifact(workspace: string, relativePath: string, contents: string): Promise<string> {
  return invoke("write_text_artifact", { workspace, relativePath, contents });
}
