import { readFileSync, existsSync } from 'node:fs';
import { join } from 'node:path';
import assert from 'node:assert/strict';
const site = (process.env.NEXT_PUBLIC_SITE_URL || 'https://dnncha.github.io/dotmatch').replace(/\/+$/, '');
const base = (process.env.NEXT_PUBLIC_BASE_PATH || '').replace(/\/+$/, '');
for (const path of ['', 'crispr-guide-counting', 'tools/library-safety', 'assignment-sensitivity']) {
  const html = readFileSync(join('out', path, 'index.html'), 'utf8');
  const expected = `${site}/${path ? `${path}/` : ''}`;
  assert.equal((html.match(/<h1(?:\s|>)/g) || []).length, 1, `${path}: one rendered H1`);
  const canonicals = html.match(/<link\b[^>]*rel="canonical"[^>]*>/g) || [];
  assert.equal(canonicals.length, 1, `${path}: exactly one canonical`);
  assert(canonicals[0].includes(`href="${expected}"`), `${path}: wrong canonical ${canonicals[0]}`);
  assert(html.includes(`href="${base}/tools/library-safety/"`), `${path}: checker not internally linked`);
  assert(!/<meta\b[^>]*name="robots"[^>]*noindex/.test(html), `${path}: unexpected noindex`);
}
const sitemap = readFileSync('out/sitemap.xml', 'utf8');
for (const path of ['', 'crispr-guide-counting/', 'tools/library-safety/', 'assignment-sensitivity/']) assert(sitemap.includes(`<loc>${site}/${path}</loc>`), `Sitemap missing ${path}`);
console.log('Four exported directory routes, rendered headings, internal links, canonicals, and sitemap verified');

const demoReport = 'examples/assignment-review/report.html';
const demoArchive = 'examples/assignment-review/dotmatch-review-example.zip';
for (const path of [demoReport, demoArchive, 'examples/assignment-review/demo-manifest.json']) {
  assert(existsSync(join('out', path)), `Missing native-generated review asset: ${path}`);
}
const sensitivityPage = readFileSync('out/assignment-sensitivity/index.html', 'utf8');
assert(sensitivityPage.includes(`href="${base}/${demoReport}"`), 'Working viewer entry is not linked');
assert(sensitivityPage.includes(`href="${base}/${demoArchive}"`), 'Portable example download is not linked');
const demoHtml = readFileSync(join('out', demoReport), 'utf8');
assert(demoHtml.includes('Synthetic example — not biological data.'), 'Demo must be visibly labelled');
assert(demoHtml.includes('id="review-data"'), 'Demo must use the real report renderer');
console.log('Native viewer report, portable download and base-path links verified');
