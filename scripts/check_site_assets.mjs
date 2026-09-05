import { existsSync, readFileSync, readdirSync, statSync } from "node:fs";
import { dirname, join, resolve } from "node:path";
import { fileURLToPath } from "node:url";
import assert from "node:assert/strict";
const root = resolve(dirname(fileURLToPath(import.meta.url)), "..");
const read = path => readFileSync(join(root, path), "utf8");
const required = ["README.md", "app/page.tsx", "app/layout.tsx", "app/robots.ts", "app/sitemap.ts", "app/globals.css", "app/site-metadata.ts", "app/research-shell.tsx", "app/research.module.css", "app/crispr-guide-counting/page.tsx", "app/tools/library-safety/page.tsx", "app/tools/library-safety/explorer.tsx", "lib/library-safety.ts", "tests/site/library-safety.test.mjs", "tests/site/metadata.test.mjs", "docs/index.md", "docs/getting-started.md", "docs/command-reference.md", "docs/agent-guide.md", "docs/agent-crispr.md", "docs/agent-perturb-seq.md", "public/llms.txt", "public/llms-full.txt", "public/agent-capabilities.json", "public/agent-capabilities.schema.json", "public/agent-tools.json", "public/agent-tools.schema.json", "public/agent-reference-crispr.json", "public/dotmatch-read-assignment-v2.webp", "public/dotmatch-read-assignment-mobile-v2.webp", "public/dotmatch-og.png", "public/dotmatch-twitter.png"];
for (const path of required) assert(existsSync(join(root, path)), `Missing public file: ${path}`);
const home = read("app/page.tsx");
for (const anchor of ["top", "workflow", "failure-modes", "use-cases", "evidence", "install", "agent-workflow"]) assert(home.includes(`id="${anchor}"`), `Missing homepage section: ${anchor}`);
for (const path of ["app/page.tsx", "app/crispr-guide-counting/page.tsx", "app/tools/library-safety/page.tsx"]) {
  const page = read(path);
  assert.equal((page.match(/<h1\b/g) ?? []).length, 1, `${path}: exactly one H1 required`);
  assert(!page.includes("github.com/dnncha/dotmatch/blob/main/docs/"), `${path}: use rendered documentation links`);
}
for (const phrase of ["CRISPR", "ambiguous", "unique", "none", "invalid", "dotmatch agent tools --json", "python3 -m pip install dotmatch", "getting-started.html", "SoftwareApplication", "softwareVersion: publishedVersion", "featureList", 'type="application/ld+json"']) assert(home.includes(phrase), `Missing scientific/task content: ${phrase}`);
for (const image of ["dotmatch-read-assignment-v2.webp", "dotmatch-read-assignment-mobile-v2.webp"]) assert(home.includes(image), `Missing responsive image: ${image}`);
assert(home.includes("<picture>") && home.includes('media="(max-width: 520px)"'), "Responsive explainer source missing");
const layout = read("app/layout.tsx"), metadata = read("app/site-metadata.ts"), sitemap = read("app/sitemap.ts");
assert(layout.includes('applicationName: "DotMatch"') && metadata.includes('siteName: "DotMatch"'), "Metadata identity missing");
assert(layout.includes('rel="describedby"') && layout.includes("llms.txt"), "Agent discovery link missing");
assert(metadata.includes('publishedVersion = "0.4.1"'), "Review the published release and its install evidence before changing the public version");
assert(home.includes("packageMetadata.version") && home.includes("Website source version"), "Keep published and source versions distinct");
assert(/^\d+\.\d+\.\d+/.test(JSON.parse(read("package.json")).version), "Source version must be semantic");
for (const route of ["crispr-guide-counting", "tools/library-safety"]) {
  assert(sitemap.includes(route), `Sitemap missing ${route}`);
  assert(read(`app/${route}/page.tsx`).includes("pageMetadata("), `Page-specific canonical missing: ${route}`);
}
assert(!sitemap.includes("new Date("), "Do not manufacture lastModified from build time");
for (const rule of [".skip-link", ":focus-visible", "prefers-reduced-motion", ".install-copy a"]) assert(read("app/globals.css").includes(rule), `Missing accessibility rule: ${rule}`);
const linkPattern = /!?\[[^\]]*\]\(([^)\s]+)(?:\s+[^)]*)?\)/g;
for (const file of ["README.md", "docs/index.md"]) {
  for (const match of read(file).matchAll(linkPattern)) {
    const destination = match[1].replace(/^<|>$/g, "");
    if (/^(#|https?:\/\/)/.test(destination)) continue;
    assert(file !== "README.md", `README link not safe on PyPI: ${destination}`);
    if (destination.startsWith("mailto:")) continue;
    assert(existsSync(resolve(dirname(join(root, file)), destination.split("#", 1)[0])), `${file}: broken local link ${destination}`);
  }
}
function markdownFiles(dir) {
  return readdirSync(join(root, dir)).flatMap(name => {
    const path = join(dir, name);
    if (name === "_build") return [];
    return statSync(join(root, path)).isDirectory() ? markdownFiles(path) : name.endsWith(".md") ? [path] : [];
  });
}
const languageFiles = ["README.md", "app/page.tsx", "app/layout.tsx", "app/crispr-guide-counting/page.tsx", "app/tools/library-safety/page.tsx", "pyproject.toml", "codemeta.json", ".zenodo.json", "CITATION.cff", ...markdownFiles("docs")];
const forbidden = ["adoption evidence", "adoption trust", "AI slop", "big wins", "evidence-bounded", "industry exposure", "massive industry impact", "next wins", "pilot conversations", "private feedback", "quote-approved", "turning private evaluation into public adoption evidence", "without turning private feedback into public evidence", "game-changing", "revolutionary", "world-class", "best-in-class", "enterprise-grade", "just works"];
for (const path of languageFiles) for (const phrase of forbidden) assert(!read(path).toLowerCase().includes(phrase.toLowerCase()), `${path}: internal or inflated language: ${phrase}`);
console.log("Public routes, accessibility hooks, release identity, documentation links, and scientific language checks passed");
