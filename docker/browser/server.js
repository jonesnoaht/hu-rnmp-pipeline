const express = require('express');
const fs = require('fs');
const path = require('path');
const morgan = require('morgan');
const childProcess = require('child_process');
const RateLimit = require('express-rate-limit');

const app = express();
const PORT = process.env.PORT || 8050;
const RESULTS_DIR = process.env.RESULTS_DIR || '/results';

// One proxy hop (ingress / Authentik). Needed so the rate limiter keys on the
// client address instead of rejecting X-Forwarded-For.
app.set('trust proxy', 1);

// Public SRA / ENA / GEO accessions. Each alternative is a fixed prefix plus
// one bounded run of digits, so the match stays linear in the input length.
const ACCESSION_RE = /^(?:SRR|ERR|DRR|SRP|ERP|DRP|GSE)\d{1,12}$|^PRJ[EDN][A-Z]\d{1,12}$/;
// Values nextflow/main.nf actually branches on (params.genome).
const ALLOWED_GENOMES = new Set(['GRCh38', 'sacCer3']);
const NEXTFLOW_BIN = 'nextflow';
const NEXTFLOW_SCRIPT = '/app/nextflow/main.nf';
const OUTPUT_LIMIT = 1000;

app.use(morgan('combined'));
app.use(express.static(path.join(__dirname, 'public')));
app.use(express.json());

// File reads and pipeline launches are expensive. Cap API traffic per client.
const apiLimiter = RateLimit({
  windowMs: 60 * 1000,
  max: 120,
  standardHeaders: true,
  legacyHeaders: false,
});

function resultsRoot() {
  return path.resolve(RESULTS_DIR);
}

// True when candidate is root itself or a path strictly under it.
// path.relative rejects the "/results" vs "/results-old" prefix trick that
// String.prototype.startsWith allows.
function isInside(root, candidate) {
  const resolvedRoot = path.resolve(root);
  const resolved = path.resolve(candidate);
  const rel = path.relative(resolvedRoot, resolved);
  if (rel === '') return true;
  if (rel === '..' || rel.startsWith(`..${path.sep}`) || path.isAbsolute(rel)) return false;
  return true;
}

function isSafeRunName(run) {
  if (typeof run !== 'string' || run.length === 0 || run.includes('\0')) return false;
  if (run.includes('/') || run.includes('\\') || run === '.' || run === '..') return false;
  return run === path.basename(run);
}

function appendCapped(current, chunk, max) {
  if (current.length >= max) return current;
  const next = current + chunk.toString();
  return next.length > max ? next.slice(0, max) : next;
}

app.get('/healthz', (req, res) => {
  res.json({ status: 'ok', timestamp: new Date().toISOString() });
});

app.use(apiLimiter);

app.get('/api/runs', (req, res) => {
  try {
    const root = resultsRoot();
    const entries = fs.readdirSync(root, { withFileTypes: true });
    const runs = entries
      .filter((entry) => entry.isDirectory())
      .map((dir) => ({
        name: dir.name,
        path: path.join(root, dir.name),
        mtime: fs.statSync(path.join(root, dir.name)).mtime,
      }));
    res.json({ runs });
  } catch {
    res.status(404).json({ error: 'Results directory not found or empty' });
  }
});

function resolveRunDir(run) {
  if (!isSafeRunName(run)) return { error: 403 };
  const lexical = path.resolve(resultsRoot(), run);
  if (!isInside(resultsRoot(), lexical)) return { error: 403 };
  let real;
  try {
    real = fs.realpathSync(lexical);
  } catch (err) {
    if (err.code === 'ENOENT') return { error: 404 };
    return { error: 403 };
  }
  if (!isInside(resultsRoot(), real)) return { error: 403 };
  return { path: real };
}

app.get('/api/runs/:run/ls', (req, res) => {
  const located = resolveRunDir(req.params.run);
  if (located.error) {
    const message = located.error === 404 ? 'Run not found' : 'Forbidden';
    return res.status(located.error).json({ error: message });
  }
  try {
    const files = walkDir(located.path);
    res.json({ run: req.params.run, files });
  } catch (err) {
    if (err.code === 'ENOENT' || err.code === 'ENOTDIR') {
      return res.status(404).json({ error: 'Run not found' });
    }
    res.status(500).json({ error: 'Failed to list run' });
  }
});

function resolveResultFile(run, filepath) {
  if (!isSafeRunName(run)) return { error: 403 };
  if (typeof filepath !== 'string' || filepath.length === 0 || filepath.includes('\0')) {
    return { error: 403 };
  }
  if (path.isAbsolute(filepath)) return { error: 403 };

  const root = resultsRoot();
  const runLexical = path.resolve(root, run);
  const relRun = path.relative(root, runLexical);
  if (path.isAbsolute(relRun) || relRun === '..' || relRun.startsWith('..' + path.sep)) {
    return { error: 403 };
  }

  let realRun;
  try {
    realRun = fs.realpathSync(runLexical);
  } catch (err) {
    if (err.code === 'ENOENT') return { error: 404 };
    return { error: 403 };
  }
  const relRealRun = path.relative(root, realRun);
  if (path.isAbsolute(relRealRun) || relRealRun === '..' || relRealRun.startsWith('..' + path.sep)) {
    return { error: 403 };
  }

  const lexical = path.resolve(realRun, filepath);
  // path.relative + a ".." rejection is the containment check the file sink
  // has to see directly. A helper call is not enough for that data flow.
  const relLexical = path.relative(realRun, lexical);
  if (
    lexical === realRun ||
    path.isAbsolute(relLexical) ||
    relLexical === '..' ||
    relLexical.startsWith('..' + path.sep)
  ) {
    return { error: 403 };
  }

  try {
    const real = fs.realpathSync(lexical);
    const relReal = path.relative(realRun, real);
    if (
      path.isAbsolute(relReal) ||
      relReal === '..' ||
      relReal.startsWith('..' + path.sep)
    ) {
      return { error: 403 };
    }
    let stat;
    try {
      stat = fs.statSync(real);
    } catch (err) {
      if (err.code === 'ENOENT') return { error: 404 };
      return { error: 403 };
    }
    if (!stat.isFile()) return { error: 404 };
    return { path: real };
  } catch (err) {
    if (err.code !== 'ENOENT') return { error: 403 };
    // Path is confined but missing. sendFile reports that asynchronously.
    return { path: lexical };
  }
}

app.get('/api/runs/:run/file/*', (req, res) => {
  const raw = req.params[0];
  const filepath = Array.isArray(raw) ? raw.join('/') : raw;
  const located = resolveResultFile(req.params.run, filepath);
  if (located.error) {
    const message = located.error === 404 ? 'File not found' : 'Forbidden';
    return res.status(located.error).json({ error: message });
  }

  const root = resultsRoot();
  const rel = path.relative(root, located.path);
  if (path.isAbsolute(rel) || rel === '..' || rel.startsWith('..' + path.sep)) {
    return res.status(403).json({ error: 'Forbidden' });
  }

  // root confines the send, and the callback catches async failures (ENOENT,
  // EISDIR) that a synchronous try/catch around sendFile never sees.
  res.sendFile(rel, { root }, (err) => {
    if (!err || res.headersSent) return;
    const status = Number(err.status || err.statusCode) || 404;
    if (status === 403) {
      res.status(403).json({ error: 'Forbidden' });
      return;
    }
    res.status(404).json({ error: 'File not found' });
  });
});

app.post('/api/run', (req, res) => {
  const body = req.body;
  if (!body || typeof body !== 'object' || Array.isArray(body)) {
    return res.status(400).json({ error: 'JSON object body required' });
  }

  const { accession, genome } = body;
  if (typeof accession !== 'string' || accession.length > 16 || !ACCESSION_RE.test(accession)) {
    return res.status(400).json({ error: 'invalid accession' });
  }

  const genomeArg = genome == null || genome === '' ? 'GRCh38' : genome;
  if (typeof genomeArg !== 'string' || !ALLOWED_GENOMES.has(genomeArg)) {
    return res.status(400).json({ error: 'invalid genome' });
  }

  // Argument vector, never a shell string. shell:false is explicit so a
  // future edit cannot turn this back into /bin/sh -c interpolation.
  const args = [
    'run',
    NEXTFLOW_SCRIPT,
    '-profile',
    'docker',
    '--accession',
    accession,
    '--genome',
    genomeArg,
  ];

  let child;
  try {
    child = childProcess.spawn(NEXTFLOW_BIN, args, { shell: false });
  } catch {
    return res.status(500).json({ error: 'failed to start pipeline' });
  }

  let stdout = '';
  let stderr = '';
  let settled = false;

  function finish(status, payload) {
    if (settled || res.headersSent) return;
    settled = true;
    res.status(status).json(payload);
  }

  if (child.stdout) {
    child.stdout.on('data', (chunk) => {
      stdout = appendCapped(stdout, chunk, OUTPUT_LIMIT);
    });
  }
  if (child.stderr) {
    child.stderr.on('data', (chunk) => {
      stderr = appendCapped(stderr, chunk, OUTPUT_LIMIT);
    });
  }

  // 'error' (spawn ENOENT) and 'close' (non-zero exit) can both fire.
  // Respond once; an uncaught throw here would take down the process.
  child.on('error', () => {
    finish(500, { error: 'failed to start pipeline' });
  });

  child.on('close', (code) => {
    if (code !== 0) {
      finish(500, { error: 'pipeline failed', stderr: stderr.slice(0, OUTPUT_LIMIT) });
      return;
    }
    finish(200, { status: 'queued', stdout: stdout.slice(0, OUTPUT_LIMIT) });
  });
});

function walkDir(dir, base = '', runRoot = dir) {
  let results = [];
  const entries = fs.readdirSync(dir, { withFileTypes: true });
  for (const entry of entries) {
    const rel = base ? `${base}/${entry.name}` : entry.name;
    const full = path.join(dir, entry.name);
    if (entry.isSymbolicLink()) {
      try {
        const real = fs.realpathSync(full);
        if (!isInside(runRoot, real)) continue;
        const stat = fs.statSync(real);
        if (!stat.isFile()) continue;
        results.push({ path: rel, size: stat.size, mtime: stat.mtime });
      } catch {
        continue;
      }
      continue;
    }
    if (entry.isDirectory()) {
      results = results.concat(walkDir(full, rel, runRoot));
      continue;
    }
    if (!entry.isFile()) continue;
    try {
      const stat = fs.statSync(full);
      results.push({ path: rel, size: stat.size, mtime: stat.mtime });
    } catch {
      continue;
    }
  }
  return results;
}

app.use((err, req, res, next) => {
  if (res.headersSent) return next(err);
  if (err.type === 'entity.parse.failed' || (err instanceof SyntaxError && err.status === 400)) {
    return res.status(400).json({ error: 'Invalid JSON' });
  }
  const status = err.status || err.statusCode || 500;
  res.status(status).json({ error: status === 404 ? 'Not found' : 'Internal error' });
});

if (require.main === module) {
  app.listen(PORT, () => {
    console.log(`HU-rNMP browser on :${PORT}`);
  });
}

module.exports = app;
