"use client";
import { useState, useRef, useEffect } from "react";
import { publishedVersion } from "./site-metadata";
import styles from "./home.module.css";
export function InstallCommand() {
  const [copied, setCopied] = useState(false);
  const [failed, setFailed] = useState(false);
  const timer = useRef<ReturnType<typeof setTimeout> | null>(null);
  const command = `python3 -m pip install dotmatch==${publishedVersion}`;
  useEffect(
    () => () => {
      if (timer.current) clearTimeout(timer.current);
    },
    [],
  );
  async function copy() {
    try {
      await navigator.clipboard.writeText(command);
      setCopied(true);
      setFailed(false);
      if (timer.current) clearTimeout(timer.current);
      timer.current = setTimeout(() => setCopied(false), 2200);
    } catch {
      setFailed(true);
    }
  }
  return (
    <div>
      <div className={styles.installCommand}>
        <code>{command}</code>
        <button
          type="button"
          onClick={copy}
          aria-label="Copy installation command"
        >
          {copied ? "Copied" : "Copy"}
        </button>
      </div>
      <span className={styles.copyStatus} role="status">
        {failed
          ? "Select the command above to copy it manually."
          : copied
            ? "Installation command copied."
            : ""}
      </span>
    </div>
  );
}
