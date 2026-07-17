import { existsSync, readFileSync } from "node:fs";

const requiredFiles = [
  "../app/page.tsx",
  "../app/layout.tsx",
  "../app/robots.ts",
  "../app/sitemap.ts",
  "../app/globals.css",
  "../next.config.ts",
  "../docs/bioinformatics-evaluation.md",
  "../docs/external-review-packet.md",
  "../docs/integration-targets.json",
  "../docs/pilot-program.md",
  "../docs/reviewer-readiness.json",
  "../docs/workflow-integration-kit.md",
  "../docs/workflow-integration-plan.json",
  "../docs/workflow-integration-roadmap.md",
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
const evaluationPacket = readFileSync(new URL("../docs/bioinformatics-evaluation.md", import.meta.url), "utf8");
const reviewPacket = readFileSync(new URL("../docs/external-review-packet.md", import.meta.url), "utf8");
const integrationTargets = JSON.parse(readFileSync(new URL("../docs/integration-targets.json", import.meta.url), "utf8"));
const reviewerReadiness = JSON.parse(readFileSync(new URL("../docs/reviewer-readiness.json", import.meta.url), "utf8"));
const integrationKit = readFileSync(new URL("../docs/workflow-integration-kit.md", import.meta.url), "utf8");
const integrationRoadmap = readFileSync(new URL("../docs/workflow-integration-roadmap.md", import.meta.url), "utf8");
const integrationPlan = JSON.parse(readFileSync(new URL("../docs/workflow-integration-plan.json", import.meta.url), "utf8"));
const socialRenderer = readFileSync(new URL("../scripts/render_social_images.py", import.meta.url), "utf8");
const nextConfig = readFileSync(new URL("../next.config.ts", import.meta.url), "utf8");
const readme = readFileSync(new URL("../README.md", import.meta.url), "utf8");

for (const anchor of [
  'id="top"',
  'id="failure-modes"',
  'id="workflow"',
  'id="industry-routes"',
  'id="evidence"',
  'id="evaluation"',
  'id="ecosystem"',
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
  ".evaluation-list",
  ".evaluation-links",
  ".ecosystem-grid",
  ".ecosystem-note",
  ".terminal",
  ".site-footer"
]) {
  if (!css.includes(selector)) {
    console.error(`Missing site CSS selector: ${selector}`);
    process.exit(1);
  }
}

const h1Matches = page.match(/<h1\b/g) ?? [];
if (h1Matches.length !== 1 || !page.includes("Design the assay. Trust the assignment.")) {
  console.error("Homepage must have exactly one required H1.");
  process.exit(1);
}

for (const phrase of [
  "Assay compilation and reliability for known-target sequencing.",
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
  "Bioinformatics evaluation packet",
  "External review packet",
  "Package channels",
  "Validated scope",
  "Output contracts",
  "Workflow status",
  "Distribution status record",
  "Workflow adoption status",
  "Integration target tracker",
  "DotMatch evaluation protocol",
  "Reviewer readiness record",
  "Workflow submission pack",
  "nf-core modules",
  "MultiQC module",
  "Galaxy / IUC",
  "bio.tools record",
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
  "docs/bioinformatics-evaluation.md",
  "docs/external-review-packet.md",
  "docs/integration-targets.json",
  "docs/pilot-program.md",
  "docs/reviewer-readiness.json",
  "docs/workflow-integration-kit.md",
  "docs/workflow-submissions.md",
  "docs/adopters/README.md",
  "docs/workflow-adoption.json",
  "docs/distribution-release.json"
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
  "331k reads" + "/s",
  "Next 10 exposure wins",
  "big wins",
  "massive industry impact",
  "Adoption trust plan",
  "External integration kit",
  "industry-" + "exposure.md",
  "industry-" + "next-wins.md",
  "adoption-" + "trust-plan.json"
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

if (!layout.includes("AssayCode") || !layout.includes("Design the assay and trust the assignment")) {
  console.error("Site metadata must match the AssayCode positioning.");
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
  "Priority Integration Work",
  "Workflow Integration Roadmap",
  "Workflow distribution handoff",
  "Methods and citation artifacts",
  "Reviewer packet",
  "Public use records",
  "Copy-Paste Outreach"
]) {
  if (!integrationKit.includes(phrase)) {
    console.error(`Workflow integration kit is missing required section: ${phrase}`);
    process.exit(1);
  }
}

if (!Array.isArray(integrationPlan.items) || integrationPlan.items.length !== 10) {
  console.error("Workflow integration plan must contain exactly 10 reviewer/integration items.");
  process.exit(1);
}

for (const item of integrationPlan.items) {
  for (const field of ["id", "title", "primary_audience", "asset", "done_when"]) {
    if (!item[field]) {
      console.error(`Workflow integration plan item missing ${field}.`);
      process.exit(1);
    }
  }
  if (!integrationRoadmap.includes(item.title)) {
    console.error(`Workflow integration roadmap missing plan title: ${item.title}`);
    process.exit(1);
  }
}

for (const phrase of [
  "Bioinformatics Evaluation Packet",
  "Current Package Surface",
  "Validated Assay Scope",
  "Minimum Local Evaluation",
  "Outputs To Inspect",
  "Workflow Integration Status",
  "Language Rules For Public Surfaces"
]) {
  if (!evaluationPacket.includes(phrase)) {
    console.error(`Evaluation packet is missing required section: ${phrase}`);
    process.exit(1);
  }
}

const credibilitySurfaces = [
  ["homepage", page],
  ["layout metadata", layout],
  ["README", readme],
  ["evaluation packet", evaluationPacket],
  ["external review packet", reviewPacket],
  ["integration kit", integrationKit],
  ["integration roadmap", integrationRoadmap]
];

for (const [label, content] of credibilitySurfaces) {
  for (const phrase of [
    "massive industry impact",
    "AI slop",
    "game-changing",
    "revolutionary",
    "world-class",
    "best-in-class",
    "enterprise-grade",
    "production-ready",
    "just works",
    "magic",
    "adoption evidence",
    "adoption trust",
    "evidence-bounded",
    "industry exposure",
    "next wins",
    "pilot conversations",
    "private feedback",
    "private pilot",
    "quote-approved",
    "turning private",
    "without turning"
  ]) {
    if (content.toLowerCase().includes(phrase.toLowerCase())) {
      console.error(`${label} contains hype or unsupported credibility language: ${phrase}`);
      process.exit(1);
    }
  }
}

if (!Array.isArray(reviewerReadiness.items) || reviewerReadiness.items.length !== 10) {
  console.error("Reviewer readiness record must contain exactly 10 concrete assets.");
  process.exit(1);
}

if (!Array.isArray(integrationTargets.targets) || integrationTargets.targets.length !== 5) {
  console.error("Integration target tracker must contain exactly five ecosystem targets.");
  process.exit(1);
}

for (const phrase of [
  "Reviewer Decision Tree",
  "Persona One-Pagers",
  "Integration Target Tracker",
  "Reviewer Packet",
  "Conference Abstracts",
  "Technical Communication Pack",
  "Maintainer Issue Templates",
  "Evaluation Scorecard",
  "Integration Tracking Metrics",
  "Release Publication Checklist"
]) {
  if (!integrationRoadmap.includes(phrase)) {
    console.error(`Workflow integration roadmap is missing required section: ${phrase}`);
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
