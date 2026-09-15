import test from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { execFileSync } from 'node:child_process';
const read = path => readFileSync(new URL('../../' + path, import.meta.url), 'utf8');
const moduleUrl = new URL('../../app/site-metadata.ts', import.meta.url).href;
for (const base of ['', '/dotmatch']) {
  test(`viewer file links preserve extensions without a directory slash (${base || 'root'})`, () => {
    const value = execFileSync(process.execPath, ['--experimental-strip-types', '--input-type=module', '-e', `import {siteAsset} from ${JSON.stringify(moduleUrl)}; console.log(siteAsset('examples/assignment-review/report.html'));`], {env: {...process.env, NEXT_PUBLIC_BASE_PATH: base}, encoding:'utf8'}).trim();
    assert.equal(value, `${base}/examples/assignment-review/report.html`);
  });
}
test('entry page distinguishes the released engine from the candidate viewer', () => {
  const page = read('app/assignment-sensitivity/page.tsx');
  assert(page.includes('siteAsset("examples/assignment-review/report.html")'));
  assert(page.includes('siteAsset("examples/assignment-review/dotmatch-review-example.zip")'));
  assert(page.includes('unreleased upgrade'));
  assert(!page.includes('not in published {publishedVersion}'));
  assert(page.includes('bundle/read_changes.tsv'));
});
test('public install surfaces use the confirmed 0.5.0 release, not stale 0.4.1', () => {
  for (const path of ['app/site-metadata.ts', 'public/llms.txt', 'public/llms-full.txt']) {
    const text = read(path); assert(text.includes('0.5.0'), path); assert(!text.includes('0.4.1'), path);
  }
});
test('site build generates the real demo before export rather than shipping a mock', () => {
  const pkg = JSON.parse(read('package.json'));
  assert.equal(pkg.scripts.prebuild, 'python3 scripts/build_sensitivity_demo.py --site');
  assert(read('scripts/build_sensitivity_demo.py').includes('run_sensitivity('));
});
