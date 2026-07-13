import { existsSync, readFileSync } from "node:fs";

const requiredFiles = [
  "../app/page.tsx",
  "../app/layout.tsx",
  "../app/robots.ts",
  "../app/sitemap.ts",
  "../app/globals.css",
  "../next.config.ts",
  "../docs/industry-exposure.md",
  "../docs/industry-exposure-plan.json",
  "../docs/industry-next-wins.md",
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
const robots = readFileSync(new URL("../app/robots.ts", import.meta.url), "utf8");
const sitemap = readFileSync(new URL("../app/sitemap.ts", import.meta.url), "utf8");
const exposureKit = readFileSync(new URL("../docs/industry-exposure.md", import.meta.url), "utf8");
const nextWinsDoc = readFileSync(new URL("../docs/industry-next-wins.md", import.meta.url), "utf8");
const exposurePlan = JSON.parse(readFileSync(new URL("../docs/industry-exposure-plan.json", import.meta.url), "utf8"));
const socialRenderer = readFileSync(new URL("../scripts/render_social_images.py", import.meta.url), "utf8");
const nextConfig = readFileSync(new URL("../next.config.ts", import.meta.url), "utf8");

for (const anchor of [
  'id="top"',
  'id="failure-modes"',
  'id="workflow"',
  'id="industry-routes"',
  'id="evidence"',
  'id="exposure"',
  'id="next-wins"',
  'id="install"'
]) {
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
  ".audience-grid",
  ".evidence-layout",
  ".exposure-list",
  ".next-wins-grid",
  ".next-wins-note",
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
  "Core facilities",
  "CRISPR screen teams",
  "Workflow maintainers",
  "Assay developers",
  "Outreach and integration kit",
  "Workflow submission pack",
  "Next 10 exposure wins",
  "Decision tree",
  "Release calendar",
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
  "docs/packaging.md",
  "docs/industry-exposure.md",
  "docs/industry-next-wins.md",
  "docs/workflow-submissions.md",
  "docs/adopters/README.md"
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

if (!page.includes('type="application/ld+json"') || !page.includes("SoftwareApplication")) {
  console.error("Homepage must include structured SoftwareApplication JSON-LD.");
  process.exit(1);
}

if (!robots.includes("MetadataRoute.Robots") || !robots.includes("sitemap.xml")) {
  console.error("Robots route must expose the sitemap.");
  process.exit(1);
}

for (const [label, route] of [["robots", robots], ["sitemap", sitemap]]) {
  if (!route.includes('export const dynamic = "force-static"')) {
    console.error(`${label} route must be static-export safe.`);
    process.exit(1);
  }
}

if (!sitemap.includes("MetadataRoute.Sitemap") || !sitemap.includes("changeFrequency")) {
  console.error("Sitemap route must expose the public homepage.");
  process.exit(1);
}

for (const phrase of [
  "The Big 5 Distribution Moves",
  "Next 10 Wins",
  "Workflow distribution handoff",
  "Citation and methods flywheel",
  "Evidence-first launch packet",
  "Public adopter record",
  "Copy-Paste Outreach"
]) {
  if (!exposureKit.includes(phrase)) {
    console.error(`Exposure kit is missing required section: ${phrase}`);
    process.exit(1);
  }
}

if (!Array.isArray(exposurePlan.items) || exposurePlan.items.length !== 10) {
  console.error("Industry exposure plan must contain exactly 10 next-win items.");
  process.exit(1);
}

for (const item of exposurePlan.items) {
  for (const field of ["id", "title", "primary_audience", "asset", "done_when"]) {
    if (!item[field]) {
      console.error(`Industry exposure plan item missing ${field}.`);
      process.exit(1);
    }
  }
  if (!nextWinsDoc.includes(item.title)) {
    console.error(`Next-win doc missing plan title: ${item.title}`);
    process.exit(1);
  }
}

for (const phrase of [
  "Evaluator Decision Tree",
  "Persona One-Pagers",
  "Integration Target Tracker",
  "Reviewer Evidence Packet",
  "Conference Abstracts",
  "Social And Forum Pack",
  "Maintainer Issue Templates",
  "Pilot Scorecard",
  "Adoption KPI Dashboard Spec",
  "Release Communications Calendar"
]) {
  if (!nextWinsDoc.includes(phrase)) {
    console.error(`Next-win doc is missing required section: ${phrase}`);
    process.exit(1);
  }
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
