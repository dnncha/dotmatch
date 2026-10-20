import { defineConfig, globalIgnores } from "eslint/config";
import js from "@eslint/js";
import tseslint from "typescript-eslint";
import reactHooks from "eslint-plugin-react-hooks";

export default defineConfig([
  js.configs.recommended,
  ...tseslint.configs.recommended,
  reactHooks.configs.flat.recommended,
  globalIgnores([
    ".next/**", "out/**", "build/**", "next-env.d.ts",
    "repos/**", "packages/**", "public/**", "docs/**",
  ]),
]);
