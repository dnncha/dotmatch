import test from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
const read = path => readFileSync(new URL('../../' + path, import.meta.url), 'utf8');

test('report sections use scientific labels rather than promotional headings', () => {
  const source = read('python/dotmatch/sensitivity_review_assets.py');
  for (const [id, title] of Object.entries({
    'outcomes-title': 'Read outcomes',
    'guides-title': 'Guide counts',
    'transitions-title': 'Read-state transitions',
    'provenance-title': 'Provenance',
  })) assert(source.includes(`<h2 id="${id}">${title}</h2>`), title);
  assert(source.includes('$("headline").textContent="Assignment sensitivity";'));
});

test('public viewer entry points do not restore removed slogans', () => {
  const paths = ['app/page.tsx', 'app/assignment-sensitivity/page.tsx',
    'app/assignment-demo.tsx', 'python/dotmatch/sensitivity_review_assets.py', 'README.md'];
  const removed = [/every read[, ]+(?:has an|accounted)/i,
    /account for every read/i, /find the difference\. follow the evidence/i,
    /same assigned total\. different guide counts/i, /a small step in your workflow/i,
    /count your guides/i, /no forced call/i];
  for (const path of paths) for (const phrase of removed)
    assert(!phrase.test(read(path)), `${path}: ${phrase}`);
});

test('the report can be opened from the homepage without scrolling to a promotion', () => {
  const home = read('app/page.tsx');
  const top = home.slice(home.indexOf('id="top"'), home.indexOf('id="install"'));
  assert(top.includes('siteAsset("examples/assignment-review/report.html")'));
  assert(top.includes('Open example report'));
  assert(top.includes('sitePath("crispr-guide-counting")'));
  assert(top.includes('sitePath("tools/library-safety")'));
});

test('shorter copy retains scientific and privacy qualifications', () => {
  const source = read('python/dotmatch/sensitivity_review_assets.py');
  assert(source.includes('These checks do not authenticate the report or validate the assay.'));
  assert(source.includes('This is not the selected pair’s changed-read count.'));
  assert(source.includes('These comparisons do not measure biological accuracy.'));
  assert(source.includes('Optional attached read IDs remain in memory until cleared'));
  assert(source.includes('A base-level alignment or a specific mismatch explanation cannot be inferred'));
});
