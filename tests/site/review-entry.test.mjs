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
  assert(page.includes('GitHub release {releaseVersion}'));
  assert(page.includes('dotmatch sensitivity-review'));
  assert(page.includes('It produces the count tables and an interactive offline report'));
  assert(page.includes('fall back to a static summary'));
  assert(!page.includes('unreleased upgrade'));
  assert(page.includes('bundle/read_changes.tsv'));
});
test('public install surfaces match the authorized release candidate', () => {
  const version = JSON.parse(read('package.json')).version;
  const metadata = read('app/site-metadata.ts');
  const release = metadata.match(/releaseVersion = "([^"]+)"/)?.[1];
  const published = metadata.match(/publishedVersion = "([^"]+)"/)?.[1];
  assert.equal(release, version);
  const request = JSON.parse(read('.github/release-request.json'));
  const record = JSON.parse(read('docs/distribution-release.json'));
  assert.equal(request.authorized, true);
  assert.equal(request.version, version);
  assert.equal(request.tag, `v${version}`);
  assert.deepEqual(request.publication_channels, ['GitHub Releases']);
  assert.equal(record.release_version, published);
  assert.equal(record.release_tag, `v${published}`);
  assert.equal(record.candidate_version, version);
  assert(['prepared_not_published', 'github_released'].includes(record.candidate_status));
  assert.equal(record.github_release.version, version);
  assert.equal(record.github_release.status, record.candidate_status === 'github_released' ? 'verified' : 'prepared');
  assert.equal(record.github_release.public_url, `https://github.com/dnncha/dotmatch/releases/tag/v${version}`);
  assert.equal(record.publication_authorized, false);
  assert.equal(record.channels.find(channel => channel.id === 'pypi').status, 'verified');
  assert.equal(record.channels.find(channel => channel.id === 'pypi').public_url, `https://pypi.org/project/dotmatch/${published}/`);

  const capabilities = JSON.parse(read('public/agent-capabilities.json'));
  assert.equal(capabilities.generated_for_version, version);
  assert.equal(capabilities.install.recommended, `python3 -m pip install https://github.com/dnncha/dotmatch/releases/download/v${version}/dotmatch-${version}.tar.gz`);
  assert(metadata.includes('export const releaseInstallCommand ='));
  for (const path of ['app/install-command.tsx', 'app/crispr-guide-counting/page.tsx', 'app/assignment-sensitivity/page.tsx']) assert(read(path).includes('releaseInstallCommand'), path);
  assert.equal(JSON.parse(read('public/agent-tools.json')).generated_for_version, version);
  assert.equal(JSON.parse(read('public/agent-reference-crispr.json')).dotmatch_version, version);

  for (const path of ['public/llms.txt', 'public/llms-full.txt']) assert(read(path).includes(version), path);
});
test('site build generates the real demo before export rather than shipping a mock', () => {
  const pkg = JSON.parse(read('package.json'));
  assert.equal(pkg.scripts.prebuild, 'python3 scripts/build_sensitivity_demo.py --site');
  assert(read('scripts/build_sensitivity_demo.py').includes('run_sensitivity('));
});
test('every homepage example id exists in the checked demo fixture', () => {
  const source = read('app/assignment-demo.tsx');
  const block = source.match(/const examples = \[([\s\S]*?)\] as const;/);
  assert(block, 'app/assignment-demo.tsx: homepage example list not found');
  const ids = [...block[1].matchAll(/\n\s*\["([^"]+)",/g)].map(match => match[1]);
  assert.deepEqual(ids, ['exact_isolated', 'one_mismatch', 'two_candidates', 'unmatched', 'short']);
  const records = new Set(JSON.parse(read('public/assignment-demo.json')).records.map(row => row.id));
  for (const id of ids) assert(records.has(id), `public/assignment-demo.json has no ${id} record for the homepage to look up`);
});
