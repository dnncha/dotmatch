import { existsSync, readFileSync } from "node:fs";

const requiredFiles = [
  "../app/page.tsx",
  "../app/layout.tsx",
  "../app/globals.css",
  "../next.config.ts",
  "../public/dotmatch-read-assignment.svg",
  "../public/dotmatch-panel-certificate.png",
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
const css = readFileSync(new URL("../app/globals.css", import.meta.url), "utf8");
const layout = readFileSync(new URL("../app/layout.tsx", import.meta.url), "utf8");
const nextConfig = readFileSync(new URL("../next.config.ts", import.meta.url), "utf8");

for (const anchor of ['id="top"', 'id="panel-design"', 'id="barcode-qc"', 'id="benchmarks"', 'id="install"', 'id="cite"', 'id="use-cases"']) {
  if (!page.includes(anchor)) {
    console.error(`Missing site section anchor: ${anchor}`);
    process.exit(1);
  }
}

for (const selector of [
  ".hero",
  ".hero-visual",
  ".panel-design-layout",
  ".panel-output-grid",
  ".metric-grid",
  ".autopsy-layout",
  ".artifact-grid",
  ".report-table",
  ".decision-grid",
  ".example-layout",
  ".status-table",
  ".terminal"
]) {
  if (!css.includes(selector)) {
    console.error(`Missing site CSS selector: ${selector}`);
    process.exit(1);
  }
}

if (!layout.includes("export const metadata") || !layout.includes("openGraph") || !layout.includes("twitter")) {
  console.error("Site metadata must include Open Graph and Twitter metadata objects.");
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
