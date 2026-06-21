import { existsSync, readFileSync } from "node:fs";

const requiredFiles = [
  "../app/page.tsx",
  "../app/layout.tsx",
  "../app/globals.css",
  "../next.config.ts",
  "../public/dotmatch-read-assignment.svg",
  "../public/dotmatch-og.png",
  "../public/dotmatch-twitter.png",
  "../scripts/render_social_images.py"
];

for (const path of requiredFiles) {
  if (!existsSync(new URL(path, import.meta.url))) {
    console.error(`Missing site file: ${path}`);
    process.exit(1);
  }
}

const page = readFileSync(new URL("../app/page.tsx", import.meta.url), "utf8");
const normalizedPage = page.replace(/\s+/g, " ");
const css = readFileSync(new URL("../app/globals.css", import.meta.url), "utf8");
const layout = readFileSync(new URL("../app/layout.tsx", import.meta.url), "utf8");
const socialRenderer = readFileSync(new URL("../scripts/render_social_images.py", import.meta.url), "utf8");
const nextConfig = readFileSync(new URL("../next.config.ts", import.meta.url), "utf8");

for (const anchor of ['id="top"', 'id="failure-modes"', 'id="workflow"', 'id="evidence"', 'id="install"']) {
  if (!page.includes(anchor)) {
    console.error(`Missing site section anchor: ${anchor}`);
    process.exit(1);
  }
}

for (const selector of [
  ".hero",
  ".positioning",
  ".assignment-figure",
  ".outcome-grid",
  ".failure-grid",
  ".workflow-grid",
  ".context-rail",
  ".evidence-layout",
  ".terminal",
  ".site-footer"
]) {
  if (!css.includes(selector)) {
    console.error(`Missing site CSS selector: ${selector}`);
    process.exit(1);
  }
}

const h1Matches = page.match(/<h1\b/g) ?? [];
if (h1Matches.length !== 1 || !page.includes("Know which read assignments you can trust.")) {
  console.error("Homepage must have exactly one required H1.");
  process.exit(1);
}

for (const phrase of [
  "Assignment reliability for known-target sequencing assays.",
  "unique, ambiguous, none, or invalid",
  "CRISPR guides",
  "inline barcodes",
  "feature tags",
  "primers / panels",
  "whitelists",
  "pip install dotmatch"
]) {
  if (!normalizedPage.includes(phrase)) {
    console.error(`Missing repositioning copy: ${phrase}`);
    process.exit(1);
  }
}

for (const phrase of [
  "docs/scientific-claims.md",
  "docs/evidence-gallery/README.md",
  "docs/methods-and-citation.md",
  "docs/packaging.md"
]) {
  if (!page.includes(phrase)) {
    console.error(`Missing evidence boundary link: ${phrase}`);
    process.exit(1);
  }
}

const aboveEvidence = page.split('id="evidence"')[0];
for (const forbidden of [
  "reads" + "/s",
  "Ham" + "ming",
  "Leven" + "shtein",
  "SIMD",
  "G" + "PU",
  "throughput",
  "benchmark"
]) {
  if (aboveEvidence.includes(forbidden)) {
    console.error(`Performance or algorithm claim appears before evidence section: ${forbidden}`);
    process.exit(1);
  }
}

for (const stale of [
  "Fast FASTQ " + "matching",
  "Design panels. " + "Count guides",
  "guide counts and barcode QC",
  "331k reads" + "/s"
]) {
  if (page.includes(stale) || layout.includes(stale) || socialRenderer.includes(stale)) {
    console.error(`Stale homepage positioning remains: ${stale}`);
    process.exit(1);
  }
}

if (!layout.includes("export const metadata") || !layout.includes("openGraph") || !layout.includes("twitter")) {
  console.error("Site metadata must include Open Graph and Twitter metadata objects.");
  process.exit(1);
}

if (!layout.includes("Assignment Reliability") || !layout.includes("Know which read assignments you can trust")) {
  console.error("Site metadata must match the new assignment reliability positioning.");
  process.exit(1);
}

if (!layout.includes("export const viewport")) {
  console.error("Site layout must export viewport metadata for mobile rendering.");
  process.exit(1);
}

if (!nextConfig.includes("devIndicators: false")) {
  console.error("Next.js dev indicator should be disabled for local screenshots.");
  process.exit(1);
}

function readPngDimensions(imagePath) {
  const png = readFileSync(new URL(imagePath, import.meta.url));
  if (png.length < 24 || png.toString("ascii", 1, 4) !== "PNG") {
    console.error(`${imagePath} is not a valid PNG.`);
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
    console.error(`${imagePath} must be 1200x630; saw ${width}x${height}.`);
    process.exit(1);
  }
}

const svg = readFileSync(new URL("../public/dotmatch-read-assignment.svg", import.meta.url), "utf8");
if (!svg.startsWith("<svg ") || !svg.includes('role="img"') || !svg.includes("<title")) {
  console.error("Workflow SVG should be a valid image asset with title metadata.");
  process.exit(1);
}

for (const outcome of ["unique", "ambiguous", "none", "invalid"]) {
  if (!svg.includes(outcome)) {
    console.error(`Workflow SVG must show the ${outcome} outcome.`);
    process.exit(1);
  }
}
