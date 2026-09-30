const fs = require('fs');
const os = require('os');
const path = require('path');
const http = require('http');
const assert = require('node:assert/strict');
const { describe, it, before, after, mock } = require('node:test');
const { EventEmitter } = require('node:events');
const childProcess = require('node:child_process');

const tmp = fs.mkdtempSync(path.join(os.tmpdir(), 'hu-rnmp-browser-'));
const resultsDir = path.join(tmp, 'results');
const runA = path.join(resultsDir, 'runA');
const runB = path.join(resultsDir, 'runB');

fs.mkdirSync(path.join(runA, 'sub'), { recursive: true });
fs.mkdirSync(runB, { recursive: true });
fs.writeFileSync(path.join(runA, 'report.txt'), 'hello report');
fs.writeFileSync(path.join(runA, 'sub', 'nested.txt'), 'nested ok');
fs.writeFileSync(path.join(runB, 'private.txt'), 'runB private');

fs.mkdirSync(path.join(tmp, 'results-old'));
fs.writeFileSync(path.join(tmp, 'results-old', 'secret.txt'), 'TOP SECRET');
fs.writeFileSync(path.join(tmp, 'outside.txt'), 'outside data');

fs.symlinkSync(path.join(tmp, 'outside.txt'), path.join(runA, 'leak.txt'));
fs.symlinkSync(path.join(tmp, 'results-old'), path.join(resultsDir, 'linked-escape'));

process.env.RESULTS_DIR = resultsDir;

const app = require('../server');

const spawnCalls = [];
let spawnPlan = { code: 0, stdout: 'launch ok\n' };

function makeChild(plan) {
  const child = new EventEmitter();
  child.stdout = new EventEmitter();
  child.stderr = new EventEmitter();
  child.stdin = new EventEmitter();
  child.kill = () => {};
  process.nextTick(() => {
    if (plan.stdout) child.stdout.emit('data', Buffer.from(plan.stdout));
    if (plan.stderr) child.stderr.emit('data', Buffer.from(plan.stderr));
    if (plan.error) child.emit('error', plan.error);
    child.emit('close', plan.error ? null : plan.code);
  });
  return child;
}

function request(method, urlPath, { body, raw, headers } = {}) {
  return new Promise((resolve, reject) => {
    const payload = raw !== undefined ? Buffer.from(raw) : body === undefined ? null : Buffer.from(JSON.stringify(body));
    const req = http.request(
      {
        hostname: '127.0.0.1',
        port,
        path: urlPath,
        method,
        headers: {
          ...(payload
            ? {
                'Content-Type': 'application/json',
                'Content-Length': payload.length,
              }
            : {}),
          ...headers,
        },
      },
      (res) => {
        const chunks = [];
        res.on('data', (chunk) => chunks.push(chunk));
        res.on('end', () => {
          const text = Buffer.concat(chunks).toString('utf8');
          let parsed = null;
          try {
            parsed = text ? JSON.parse(text) : null;
          } catch {
            parsed = null;
          }
          resolve({
            status: res.statusCode,
            body: parsed,
            raw: text,
            contentType: String(res.headers['content-type'] || ''),
          });
        });
      },
    );
    req.setTimeout(4000, () => {
      req.destroy(new Error(`timeout ${method} ${urlPath}`));
    });
    req.on('error', reject);
    if (payload) req.end(payload);
    else req.end();
  });
}

let server;
let port;
const crashes = [];

function onCrash(err) {
  crashes.push(err);
}

describe('results browser security', () => {
  before(async () => {
    mock.method(childProcess, 'spawn', (cmd, args, opts) => {
      spawnCalls.push({ cmd, args, opts });
      return makeChild(spawnPlan);
    });
    process.on('uncaughtException', onCrash);
    server = app.listen(0, '127.0.0.1');
    await new Promise((resolve) => server.on('listening', resolve));
    port = server.address().port;
  });

  after(async () => {
    process.off('uncaughtException', onCrash);
    mock.restoreAll();
    if (server) {
      await new Promise((resolve, reject) => {
        server.close((err) => (err ? reject(err) : resolve()));
      });
    }
    fs.rmSync(tmp, { recursive: true, force: true });
  });

  it('serves health and a file that stays inside the run', async () => {
    const health = await request('GET', '/healthz');
    assert.equal(health.status, 200);
    assert.equal(health.body.status, 'ok');

    const file = await request('GET', '/api/runs/runA/file/report.txt');
    assert.equal(file.status, 200);
    assert.equal(file.raw, 'hello report');

    const nested = await request('GET', '/api/runs/runA/file/sub/nested.txt');
    assert.equal(nested.status, 200);
    assert.equal(nested.raw, 'nested ok');

    const listing = await request('GET', '/api/runs/runA/ls');
    assert.equal(listing.status, 200);
    const names = listing.body.files.map((entry) => entry.path).sort();
    assert.deepEqual(names, ['report.txt', 'sub/nested.txt']);
  });

  it('rejects dot-dot run names and does not list the parent directory', async () => {
    const mark = spawnCalls.length;
    const encoded = await request('GET', '/api/runs/%2E%2E/ls');
    assert.equal(encoded.status, 403);
    assert.equal(encoded.body.error, 'Forbidden');
    assert.equal(encoded.raw.includes('secret.txt'), false);
    assert.equal(encoded.raw.includes('outside.txt'), false);
    assert.equal(encoded.raw.includes('TOP SECRET'), false);

    const rawDots = await request('GET', '/api/runs/../ls');
    assert.notEqual(rawDots.status, 200);
    assert.equal(rawDots.raw.includes('TOP SECRET'), false);
    assert.equal(rawDots.raw.includes('secret.txt'), false);
    assert.equal(spawnCalls.length, mark);
  });

  it('rejects a sibling directory whose path shares the results prefix', async () => {
    const run = encodeURIComponent('../results-old');
    const listing = await request('GET', `/api/runs/${run}/ls`);
    assert.equal(listing.status, 403);
    assert.equal(listing.body.error, 'Forbidden');
    assert.equal(listing.raw.includes('TOP SECRET'), false);
    assert.equal(listing.raw.includes('secret.txt'), false);

    const file = await request('GET', `/api/runs/${run}/file/secret.txt`);
    assert.equal(file.status, 403);
    assert.equal(file.body.error, 'Forbidden');
    assert.equal(file.raw.includes('TOP SECRET'), false);
  });

  it('rejects encoded traversal in the file path, including another run', async () => {
    const escaped = encodeURIComponent('../../results-old/secret.txt');
    const sibling = await request('GET', `/api/runs/runA/file/${escaped}`);
    assert.equal(sibling.status, 403);
    assert.equal(sibling.raw.includes('TOP SECRET'), false);

    const slashes = await request('GET', '/api/runs/runA/file/../../results-old/secret.txt');
    assert.equal(slashes.status, 403);
    assert.equal(slashes.raw.includes('TOP SECRET'), false);

    const otherRun = await request('GET', '/api/runs/runA/file/../runB/private.txt');
    assert.equal(otherRun.status, 403);
    assert.equal(otherRun.raw.includes('runB private'), false);
  });

  it('rejects symlinks that point outside the results directory', async () => {
    const leaked = await request('GET', '/api/runs/runA/file/leak.txt');
    assert.equal(leaked.status, 403);
    assert.equal(leaked.raw.includes('outside data'), false);

    const linked = await request('GET', '/api/runs/linked-escape/ls');
    assert.equal(linked.status, 403);
    assert.equal(linked.raw.includes('TOP SECRET'), false);
    assert.equal(linked.raw.includes('secret.txt'), false);

    const runs = await request('GET', '/api/runs');
    assert.equal(runs.status, 200);
    assert.deepEqual(runs.body.runs.map((run) => run.name).sort(), ['runA', 'runB']);
  });

  it('returns JSON 404 for a missing file without crashing', async () => {
    const before = crashes.length;
    const missing = await request('GET', '/api/runs/runA/file/missing.txt');
    assert.equal(missing.status, 404);
    assert.equal(missing.body.error, 'File not found');
    assert.equal(missing.contentType.includes('application/json'), true);

    const health = await request('GET', '/healthz');
    assert.equal(health.status, 200);
    assert.equal(crashes.length, before);
  });

  it('rejects shell injection in accession and genome before spawn', async () => {
    const payloads = [
      { accession: 'SRR123456; id' },
      { accession: 'SRR123456 && id' },
      { accession: 'SRR123456 | id' },
      { accession: '$(id)' },
      { accession: '`id`' },
      { accession: 'SRR123456\nid' },
      { accession: 'SRR123456; cat /etc/passwd' },
      { accession: 'GSE123; curl http://evil' },
      { accession: '../SRR123456' },
      { accession: 'SRR123456', genome: 'GRCh38; id' },
      { accession: 'SRR123456', genome: 'sacCer3 && id' },
      { accession: 'SRR123456', genome: 'GRCh38`id`' },
      { accession: 'SRR123456', genome: '/refs/evil' },
      { accession: 'SRR123456', genome: 'GRCh38\n--profile' },
      { accession: { toString: 'SRR123456' } },
    ];

    for (const body of payloads) {
      const mark = spawnCalls.length;
      const before = crashes.length;
      const res = await request('POST', '/api/run', { body });
      assert.equal(res.status, 400, `expected 400 for ${JSON.stringify(body)}`);
      assert.equal(spawnCalls.length, mark, `spawn called for ${JSON.stringify(body)}`);
      assert.equal(crashes.length, before);
      assert.equal(res.raw.includes('/etc/passwd'), false);
    }
  });

  it('spawns nextflow with an argument array and shell disabled', async () => {
    spawnPlan = { code: 0, stdout: 'launch ok\n' };
    const mark = spawnCalls.length;
    const res = await request('POST', '/api/run', {
      body: { accession: 'SRP123456', genome: 'sacCer3' },
    });
    assert.equal(res.status, 200);
    assert.equal(res.body.status, 'queued');
    assert.equal(res.body.stdout.includes('launch ok'), true);

    const call = spawnCalls[mark];
    assert.equal(call.cmd, 'nextflow');
    assert.equal(call.opts.shell, false);
    assert.deepEqual(call.args, [
      'run',
      '/app/nextflow/main.nf',
      '-profile',
      'docker',
      '--accession',
      'SRP123456',
      '--genome',
      'sacCer3',
    ]);
    assert.equal(call.args.includes('SRP123456'), true);
    assert.equal(call.args.some((arg) => /[;&|`$]/.test(arg)), false);

    const defaultGenome = await request('POST', '/api/run', {
      body: { accession: 'PRJNA123456' },
    });
    assert.equal(defaultGenome.status, 200);
    const genomeArg = spawnCalls[mark + 1].args;
    assert.equal(genomeArg[genomeArg.indexOf('--genome') + 1], 'GRCh38');
  });

  it('returns 500 JSON when the pipeline fails to start or exits non-zero', async () => {
    const before = crashes.length;
    spawnPlan = {
      error: Object.assign(new Error('spawn nextflow ENOENT'), { code: 'ENOENT' }),
    };
    const started = await request('POST', '/api/run', {
      body: { accession: 'ERR12345', genome: 'GRCh38' },
    });
    assert.equal(started.status, 500);
    assert.equal(started.body.error, 'failed to start pipeline');
    assert.equal(typeof started.body.error, 'string');

    spawnPlan = { code: 1, stderr: 'nextflow failed\n' };
    const failed = await request('POST', '/api/run', {
      body: { accession: 'GSE12345' },
    });
    assert.equal(failed.status, 500);
    assert.equal(failed.body.error, 'pipeline failed');
    assert.equal(failed.body.stderr.includes('nextflow failed'), true);
    assert.equal(crashes.length, before);

    const health = await request('GET', '/healthz');
    assert.equal(health.status, 200);
    assert.equal(health.body.status, 'ok');
  });
});
