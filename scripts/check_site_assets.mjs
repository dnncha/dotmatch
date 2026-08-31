import { existsSync, readFileSync, readdirSync, statSync } from "node:fs";
import { dirname, join, resolve } from "node:path";
import { fileURLToPath } from "node:url";

const root = resolve(dirname(fileURLToPath(import.meta.url)), "..");

const requiredFiles = [
  "README.md",
  "app/page.tsx",
  "app/layout.tsx",
  "app/robots.ts",
  "app/sitemap.ts",
  "app/globals.css",
  "docs/index.md",
  "docs/getting-started.md",
  "docs/command-reference.md",
  "docs/agent-guide.md",
  "docs/agent-crispr.md",
  "docs/agent-perturb-seq.md",
  "public/llms.txt",
  "public/llms-full.txt",
  "public/agent-capabilities.json",
  "public/agent-capabilities.schema.json",
  "public/agent-tools.json",
  "public/agent-tools.schema.json",
  "public/agent-reference-crispr.json",
  "public/dotmatch-read-assignment-v2.webp",
  "public/dotmatch-read-assignment-mobile-v2.webp",
  "public/dotmatch-og.png",
  "public/dotmatch-twitter.png"
];

for (const path of requiredFiles) {
  if (!existsSync(join(root, path))) {
    console.error(`Missing public file: ${path}`);
    process.exit(1);
  }
}

const read = (path) => readFileSync(join(root, path), "utf8");
const page = read("app/page.tsx");
const layout = read("app/layout.tsx");
const readme = read("README.md");
const docsIndex = read("docs/index.md");
const packageMetadata = JSON.parse(read("package.json"));
const normalizedPage = page.replace(/\s+/g, " ");

for (const anchor of [
  'id="top"',
  'id="agent-workflow"',
  'id="failure-modes"',
  'id="workflow"',
  'id="use-cases"',
  'id="evidence"',
  'id="install"'
]) {
  if (!page.includes(anchor)) {
    console.error(`Missing homepage section: ${anchor}`);
    process.exit(1);
  }
}

if ((page.match(/<h1\b/g) ?? []).length !== 1) {
  console.error("Homepage must contain exactly one H1.");
  process.exit(1);
}

for (const phrase of [
  "DotMatch",
  "Known-target sequencing read assignment",
  "Match reads without hiding uncertainty.",
  "unique",
  "ambiguous",
  "none",
  "invalid",
  "CRISPR guides",
  "inline barcodes",
  "It is not a genome aligner or basecaller.",
  "Run with a local agent",
  "dotmatch agent tools --json",
  "structured verdict remains <code>failed</code>",
  "Published ${publishedVersion}; candidate ${releaseVersion}",
  "python3 -m pip install dotmatch",
  "getting-started.html",
  "https://dotmatch.readthedocs.io/en/latest/"
]) {
  if (!normalizedPage.includes(phrase)) {
    console.error(`Homepage is missing required user-facing content: ${phrase}`);
    process.exit(1);
  }
}

if (page.includes("github.com/dnncha/dotmatch/blob/main/docs/")) {
  console.error("Homepage documentation links must point to the rendered documentation site.");
  process.exit(1);
}

for (const image of [
  "dotmatch-read-assignment-v2.webp",
  "dotmatch-read-assignment-mobile-v2.webp"
]) {
  if (!page.includes(image)) {
    console.error(`Homepage does not reference its responsive explainer image: ${image}`);
    process.exit(1);
  }
}

if (!page.includes("<picture>") || !page.includes('media="(max-width: 520px)"')) {
  console.error("Homepage explainer must provide a dedicated mobile image source.");
  process.exit(1);
}

if (!layout.includes('applicationName: "DotMatch"') ||
    !layout.includes('siteName: "DotMatch"') ||
    !layout.includes("known short DNA targets")) {
  console.error("Site metadata must describe DotMatch consistently.");
  process.exit(1);
}

if (!page.includes('type="application/ld+json"') || !page.includes("SoftwareApplication")) {
  console.error("Homepage must include SoftwareApplication structured data.");
  process.exit(1);
}

if (!layout.includes('rel="describedby"') || !layout.includes("llms.txt")) {
  console.error("Homepage metadata must point agents to llms.txt.");
  process.exit(1);
}

if (!page.includes("featureList") || !page.includes("agent-capabilities.json")) {
  console.error("Homepage structured data must expose agent task vocabulary and capability JSON.");
  process.exit(1);
}

if (!page.includes('import packageMetadata from "../package.json"') ||
    !page.includes("softwareVersion: publishedVersion") ||
    !page.includes('const publishedVersion = "0.3.1"') ||
    !page.includes("Published ${publishedVersion}; candidate ${releaseVersion}") ||
    !page.includes("Version ${releaseVersion} source candidate")) {
  console.error("Homepage must separate candidate source metadata from the current published package.");
  process.exit(1);
}

if (typeof packageMetadata.version !== "string" || !/^\d+\.\d+\.\d+/.test(packageMetadata.version)) {
  console.error("package.json must declare a semantic release version.");
  process.exit(1);
}

const css = read("app/globals.css");
for (const accessibilityRule of [
  ".skip-link",
  ":focus-visible",
  "prefers-reduced-motion",
  ".install-copy a"
]) {
  if (!css.includes(accessibilityRule)) {
    console.error(`Homepage accessibility rule is missing: ${accessibilityRule}`);
    process.exit(1);
  }
}

const markdownLinkPattern = /!?\[[^\]]*\]\(([^)\s]+)(?:\s+[^)]*)?\)/g;
for (const [label, text] of [["README.md", readme], ["docs/index.md", docsIndex]]) {
  for (const match of text.matchAll(markdownLinkPattern)) {
    const destination = match[1].replace(/^<|>$/g, "");
    if (destination.startsWith("#") ||
        destination.startsWith("https://") ||
        destination.startsWith("http://") ||
        destination.startsWith("mailto:")) {
      continue;
    }
    const sourcePath = join(root, label);
    const pathPart = destination.split("#", 1)[0];
    const resolved = resolve(dirname(sourcePath), pathPart);
    if (!existsSync(resolved)) {
      console.error(`${label} contains a missing local link: ${destination}`);
      process.exit(1);
    }
  }
}

for (const match of readme.matchAll(markdownLinkPattern)) {
  const destination = match[1].replace(/^<|>$/g, "");
  if (!destination.startsWith("https://") && !destination.startsWith("http://") && !destination.startsWith("#")) {
    console.error(`README link is not safe when rendered on PyPI: ${destination}`);
    process.exit(1);
  }
}

const publicLanguageFiles = [
  "README.md",
  "app/page.tsx",
  "app/layout.tsx",
  "pyproject.toml",
  "codemeta.json",
  ".zenodo.json",
  "CITATION.cff"
];

function markdownFiles(directory) {
  const results = [];
  for (const name of readdirSync(join(root, directory))) {
    const relative = join(directory, name);
    const full = join(root, relative);
    if (name === "_build") continue;
    if (statSync(full).isDirectory()) results.push(...markdownFiles(relative));
    else if (name.endsWith(".md")) results.push(relative);
  }
  return results;
}

publicLanguageFiles.push(...markdownFiles("docs"));

const forbiddenPhrases = [
  "adoption evidence",
  "adoption trust",
  "AI slop",
  "big wins",
  "evidence-bounded",
  "industry exposure",
  "massive industry impact",
  "next wins",
  "pilot conversations",
  "private feedback",
  "quote-approved",
  "turning private evaluation into public adoption evidence",
  "without turning private feedback into public evidence",
  "game-changing",
  "revolutionary",
  "world-class",
  "best-in-class",
  "enterprise-grade",
  "just works"
];

for (const path of publicLanguageFiles) {
  const text = read(path).toLowerCase();
  for (const phrase of forbiddenPhrases) {
    if (text.includes(phrase.toLowerCase())) {
      console.error(`${path} contains internal or inflated public language: ${phrase}`);
      process.exit(1);
    }
  }
}

console.log("Public site, documentation links, PyPI README links, and language checks passed");
