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
  ["autopsy section anchor", "id=\"autopsy\""],
  ["install section anchor", "id=\"install\""],
  ["citation section anchor", "id=\"cite\""],
  ["packaging docs link", "docs/packaging.md"],
  ["methods docs link", "docs/methods-and-citation.md"],
  ["barcode science link", "docs/barcode-science-readiness.md"],
  ["barcode benchmark link", "docs/benchmarks/barcode_demux/README.md"],
  ["benchmark evidence link", "docs/benchmarks/public_crispr/README.md"],
  ["citation file link", "CITATION.cff"],
  ["plain maintainer voice", "DotMatch is a small C/Python tool"],
  ["barcode hero", "Barcode autopsy for fixed-window FASTQ assays."],
  ["barcode flagship command", "dotmatch barcode autopsy"],
  ["autopsy report artifact", "report.html"],
  ["autopsy findings artifact", "findings.tsv"],
  ["autopsy provenance artifact", "provenance.json"],
  ["bounded replacement claim", "DotMatch does not replace BCL Convert, genome aligners, or general adapter trimming."],
  ["auditable benchmark framing", "Speed is shown only after the comparator semantics are documented."],
  ["honest install framing", "Use the source install until the public package channels finish publication."],
  ["human caveat framing", "We are keeping the claims narrow"],
  ["fixed-window story", "post-FASTQ fixed-window assignment"],
  ["science gate command", "make repository-ready"],
  ["decision box inclusion", "Use DotMatch when you have"],
  ["decision box exclusion", "Do not use DotMatch for"],
  ["hamming translation", "allow one mismatch, no indels"],
  ["levenshtein translation", "allow one substitution, insertion, or deletion"],
  ["plain benchmark sentence", "DotMatch Hamming k=1 processed about 331k reads/s"],
  ["validation sample size", "2,000 reads"],
  ["distribution maturity now", "Current distribution: source install and repository release artifacts, with Bioconda review tracked in PR #65367."],
  ["distribution maturity claim gate", "Channel availability is claimed only after the package appears on that channel."],
  ["Bioconda PR link", "https://github.com/bioconda/bioconda-recipes/pull/65367"],
  ["ambiguity example setup", "Some tools may pick or double-count."],
  ["ambiguity example result", "DotMatch reports: ambiguous"],
  ["workflow status table", "Comparator-backed, bounded"],
  ["generated barcode visual", "dotmatch-barcode-autopsy.png"],
  ["human hero visual caption", "FASTQ reads become unique, ambiguous, none, and invalid outcomes"]
];

const missing = requiredSnippets.filter(([, snippet]) => !pageNormalized.includes(snippet));
const bannedPhrases = [
  "utterly dominate",
  "replace cutadapt",
  "100% scientifically accurate",
  "launch path",
  "evidence trail",
  "current scope",
  "boundary",
  "workflow-ready",
  "strongest public claim",
  "checked artifacts",
  "current support"
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
  [".autopsy-layout", "The barcode autopsy demo should have a prominent two-column layout."],
  [".artifact-grid", "Barcode autopsy outputs should be scannable."],
  [".finding-list", "Autopsy findings should stay visible as scientific diagnosis examples."],
  [".report-table", "The scientist-readable report preview should remain visible."],
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

if (!layoutNormalized.toLowerCase().includes("barcode splits, crispr guide counts, and qc reports")) {
  console.error("Metadata should describe practical user outcomes, not just implementation mechanics.");
  process.exit(1);
}

const requiredMetadataSnippets = [
  "DotMatch - Barcode Autopsy for Fixed-Window FASTQs",
  "openGraph: {",
  "siteName: \"DotMatch\"",
  "locale: \"en_US\"",
  "secureUrl: socialImageUrl",
  "type: \"image/png\"",
  "DotMatch social preview showing fixed-window barcode and guide assignment with auditable outcomes"
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
  const imagePath = "../public/dotmatch-barcode-autopsy.png";
  const { width, height } = readPngDimensions(imagePath);
  if (width < 1200 || height < 900) {
    console.error(`${imagePath} should be a substantial generated barcode autopsy visual; saw ${width}x${height}.`);
    process.exit(1);
  }
}
