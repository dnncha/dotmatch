import { readFileSync } from "node:fs";

const page = readFileSync(new URL("../app/page.tsx", import.meta.url), "utf8");
const css = readFileSync(new URL("../app/globals.css", import.meta.url), "utf8");
const layout = readFileSync(new URL("../app/layout.tsx", import.meta.url), "utf8");
const nextConfig = readFileSync(new URL("../next.config.ts", import.meta.url), "utf8");
const renderSocialImages = readFileSync(
  new URL("../scripts/render_social_images.py", import.meta.url),
  "utf8"
);
const pageNormalized = page.replace(/\s+/g, " ");
const layoutNormalized = layout.replace(/\s+/g, " ");

const requiredSnippets = [
  ["GitHub source link", "https://github.com/dnncha/dotmatch"],
  ["install section anchor", "id=\"install\""],
  ["citation section anchor", "id=\"cite\""],
  ["packaging docs link", "docs/packaging.md"],
  ["methods docs link", "docs/methods-and-citation.md"],
  ["benchmark evidence link", "docs/benchmarks/public_crispr/README.md"],
  ["citation file link", "CITATION.cff"],
  ["plain maintainer voice", "DotMatch is a small C/Python tool"],
  ["auditable benchmark framing", "These rows are not a leaderboard."],
  ["honest install framing", "The documented install path is a source build or local Python package install from the repository."],
  ["human benchmark framing", "Benchmarks you can inspect, not just quote."],
  ["biology-first hero", "Known-target DNA assignment for guide counts and barcode reads."],
  ["validated use-case lead", "Current public benchmark: CRISPR guide counting"],
  ["decision box inclusion", "Use DotMatch for"],
  ["decision box exclusion", "Use other tools for"],
  ["hamming translation", "allow one mismatch, no indels"],
  ["levenshtein translation", "allow one substitution, insertion, or deletion"],
  ["plain benchmark sentence", "DotMatch Hamming k=1 processed about 331k reads/s"],
  ["validation sample size", "2,000 checked reads"],
  ["distribution status", "package-channel availability is documented without presenting unpublished channels as install paths"],
  ["Bioconda PR link", "https://github.com/bioconda/bioconda-recipes/pull/65367"],
  ["ambiguity example setup", "Some tools may pick or double-count."],
  ["ambiguity example result", "DotMatch reports: ambiguous"],
  ["workflow status table", "Public benchmark"],
  ["generated hero visual", "dotmatch-hero-workflow.png"],
  ["human hero visual caption", "Reads move into known targets; ambiguous and unmatched lanes stay visible."]
];

const missing = requiredSnippets.filter(([, snippet]) => !pageNormalized.includes(snippet));
const bannedPhrases = [
  "launch path",
  "evidence trail",
  "current scope",
  "boundary",
  "workflow-ready",
  "strongest public claim",
  "checked artifacts",
  "current support",
  "claim gates",
  "coming next",
  "ci green",
  "validated now",
  "supported, bounded",
  "best-supported workflow"
];
const checkedCopyLower = `${page}\n${layout}`.toLowerCase();
const banned = bannedPhrases.filter((phrase) => checkedCopyLower.includes(phrase));

if (missing.length > 0) {
  console.error("Missing launch-surface affordances:");
  for (const [label, snippet] of missing) {
    console.error(`- ${label}: ${snippet}`);
  }
  process.exit(1);
}

if (banned.length > 0) {
  console.error("Copy still contains release-note or machine-like phrasing:");
  for (const phrase of banned) {
    console.error(`- ${phrase}`);
  }
  process.exit(1);
}

if (!css.includes(".sequence-rail {\n  position: relative;")) {
  console.error("The hero sequence rail must stay in normal flow to avoid overlapping metric text.");
  process.exit(1);
}

if (css.includes("min-height: 360px;")) {
  console.error("Launch cards should size to their content; fixed tall cards create empty mobile space.");
  process.exit(1);
}

if (!css.includes(".launch-card {\n  min-width: 0;")) {
  console.error("Launch cards need min-width: 0 so long commands cannot force mobile overflow.");
  process.exit(1);
}

if (!css.includes("overflow-wrap: anywhere;")) {
  console.error("Launch command text should wrap on mobile instead of hiding the repository URL.");
  process.exit(1);
}

if (!nextConfig.includes("devIndicators: false")) {
  console.error("Disable the local Next.js dev indicator so it does not look like part of the site.");
  process.exit(1);
}

const requiredCss = [
  [".hero-visual", "The generated workflow image should stay visible in the hero."],
  [".decision-grid", "The near-hero decision box should stay styled and visible."],
  [".translation-grid", "Jargon translations need a distinct scannable layout."],
  [".example-layout", "The biological example needs a stable two-column desktop layout."],
  [".status-table", "Workflow maturity should remain separated from benchmark charts."],
  [".ambiguity-example", "The ambiguity story should remain concrete."]
];

const missingCss = requiredCss.filter(([selector]) => !css.includes(selector));
if (missingCss.length > 0) {
  console.error("Missing adoption-focused CSS hooks:");
  for (const [selector, message] of missingCss) {
    console.error(`- ${selector}: ${message}`);
  }
  process.exit(1);
}

if (!layoutNormalized.includes("CRISPR guide counts, barcode splits, and QC tables")) {
  console.error("Metadata should describe practical user outcomes, not just implementation mechanics.");
  process.exit(1);
}

const requiredMetadataSnippets = [
  "DotMatch - Known-Target DNA Assignment",
  "openGraph: {",
  "siteName: \"DotMatch\"",
  "locale: \"en_US\"",
  "secureUrl: socialImageUrl",
  "type: \"image/png\"",
  "DotMatch social preview showing CRISPR guide-count assignment into known DNA target rows"
];

const missingMetadata = requiredMetadataSnippets.filter((snippet) => !layout.includes(snippet));
if (missingMetadata.length > 0) {
  console.error("Missing rich Open Graph/Twitter metadata:");
  for (const snippet of missingMetadata) {
    console.error(`- ${snippet}`);
  }
  process.exit(1);
}

if (!renderSocialImages.includes("public\" / \"dotmatch-social-art.png\"")) {
  console.error("Social image renderer should use the project-local generated background source.");
  process.exit(1);
}

function readPngDimensions(imagePath) {
  const png = readFileSync(new URL(imagePath, import.meta.url));
  if (png.length < 24) {
    console.error(`${imagePath} is too small to be a valid PNG.`);
    process.exit(1);
  }
  return {
    width: png.readUInt32BE(16),
    height: png.readUInt32BE(20)
  };
}

for (const imagePath of ["../public/dotmatch-og.png", "../public/dotmatch-twitter.png"]) {
  const { width, height } = readPngDimensions(imagePath);
  if (width !== 1200 || height !== 630) {
    console.error(`${imagePath} must be 1200x630 for Open Graph/Twitter previews; saw ${width}x${height}.`);
    process.exit(1);
  }
}

{
  const imagePath = "../public/dotmatch-hero-workflow.png";
  const { width, height } = readPngDimensions(imagePath);
  if (width < 1200 || height < 650) {
    console.error(`${imagePath} should be a substantial generated hero visual; saw ${width}x${height}.`);
    process.exit(1);
  }
}
