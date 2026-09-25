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
test('entry page distinguishes the public website review from the installed CLI report', () => {
  const page = read('app/assignment-sensitivity/page.tsx');
  assert(page.includes('siteAsset("examples/assignment-review/report.html")'));
  assert(page.includes('siteAsset("examples/assignment-review/dotmatch-review-example.zip")'));
  assert(page.includes('published {publishedVersion}'));
  assert(page.includes('still writes its static report'));
  assert(!page.includes('unreleased upgrade'));
  assert(page.includes('bundle/read_changes.tsv'));
});
test('public install surfaces match the authorized release candidate', () => {
  const version = JSON.parse(read('package.json')).version;
  const published = read('app/site-metadata.ts').match(/publishedVersion = "([^"]+)"/)?.[1];
  assert.equal(published, version);
  const request = JSON.parse(read('.github/release-request.json'));
  const record = JSON.parse(read('docs/distribution-release.json'));
  assert.equal(request.authorized, true);
  assert.equal(request.version, version);
  assert.equal(request.tag, `v${version}`);
  assert.equal(record.publication_authorized, true);
  assert.equal(record.release_version, version);
  assert.equal(record.release_tag, `v${version}`);

  const capabilities = JSON.parse(read('public/agent-capabilities.json'));
  assert.equal(capabilities.generated_for_version, version);
  assert.equal(capabilities.install.recommended, `python3 -m pip install dotmatch==${version}`);
  assert.equal(JSON.parse(read('public/agent-tools.json')).generated_for_version, version);
  assert.equal(JSON.parse(read('public/agent-reference-crispr.json')).dotmatch_version, version);

  for (const path of ['public/llms.txt', 'public/llms-full.txt']) assert(read(path).includes(version), path);
});
test('site build generates the real demo before export rather than shipping a mock', () => {
  const pkg = JSON.parse(read('package.json'));
  assert.equal(pkg.scripts.prebuild, 'python3 scripts/build_sensitivity_demo.py --site');
  assert(read('scripts/build_sensitivity_demo.py').includes('run_sensitivity('));
});
