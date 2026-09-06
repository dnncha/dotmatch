import test from 'node:test';
import assert from 'node:assert/strict';
import { execFileSync } from 'node:child_process';
const moduleUrl = new URL('../../app/site-metadata.ts', import.meta.url).href;
function evaluate(env) {
  return JSON.parse(execFileSync(process.execPath, ['--experimental-strip-types', '--input-type=module', '-e', `import { sitePath, canonicalUrl, pageMetadata } from ${JSON.stringify(moduleUrl)};console.log(JSON.stringify({home:sitePath(),path:sitePath('/tools/library-safety/'),canonical:canonicalUrl('/tools/library-safety/'),metadata:pageMetadata('Library audit','Description','tools/library-safety')}))`], { env: { ...process.env, ...env }, encoding: 'utf8', stdio: ['ignore','pipe','ignore'] }));
}
test('GitHub Pages links preserve subpath exactly once and use page-specific canonical', () => {
  const result = evaluate({ NEXT_PUBLIC_BASE_PATH: '/dotmatch/', NEXT_PUBLIC_SITE_URL: 'https://dnncha.github.io/dotmatch/' });
  assert.equal(result.home, '/dotmatch/'); assert.equal(result.path, '/dotmatch/tools/library-safety/');
  assert.equal(result.canonical, 'https://dnncha.github.io/dotmatch/tools/library-safety/');
  assert.equal(result.metadata.alternates.canonical, result.canonical); assert.equal(result.metadata.openGraph.url, result.canonical);
});
test('root-hosted build does not inherit a project base path', () => {
  const result = evaluate({ NEXT_PUBLIC_BASE_PATH: '', NEXT_PUBLIC_SITE_URL: 'https://example.org/' });
  assert.equal(result.home, '/'); assert.equal(result.path, '/tools/library-safety/');
  assert.equal(result.canonical, 'https://example.org/tools/library-safety/'); assert.equal(result.metadata.twitter.title, 'Library audit');
});
