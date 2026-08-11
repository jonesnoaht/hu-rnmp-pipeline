const express = require('express');
const fs = require('fs');
const path = require('path');
const morgan = require('morgan');
const app = express();
const PORT = process.env.PORT || 8050;
const RESULTS_DIR = process.env.RESULTS_DIR || '/results';

app.use(morgan('combined')); // logs via k8s
app.use(express.static(path.join(__dirname, 'public')));

// JSON body parsing for job submission
app.use(express.json());

// Healthz (k8s probes)
app.get('/healthz', (req, res) => {
  res.json({ status: 'ok', timestamp: new Date().iso() });
});

// List all results runs
app.get('/api/runs', (req, res) => {
  try {
    const entries = fs.readdirSync(RESULTS_DIR, { withFileTypes: true });
    const runs = entries
      .filter(e => e.isDirectory())
      .map(dir => ({
        name: dir.name,
        path: path.join(RESULTS_DIR, dir.name),
        mtime: fs.statSync(path.join(RESULTS_DIR, dir.name)).mtime,
      }));
    res.json({ runs });
  } catch {
    res.status(404).json({ error: 'Results directory not found or empty' });
  }
});

// Get a specific run (list files)
app.get('/api/runs/:run/ls', (req, res) => {
  const runDir = path.join(RESULTS_DIR, req.params.run);
  try {
    const files = walkDir(runDir);
    res.json({ run: req.params.run, files });
  } catch {
    res.status(404).json({ error: 'Run not found' });
  }
});

// Serve a specific file (data download)
app.get('/api/runs/:run/file/:filepath', (req, res) => {
  const filePath = path.join(RESULTS_DIR, req.params.run, req.params.filepath);
  if (!filePath.startsWith(RESULTS_DIR)) {
    return res.status(403).json({ error: 'Forbidden' });
  }
  try {
    res.sendFile(filePath);
  } catch {
    res.status(404).json({ error: 'File not found' });
  }
});

// Run pipeline on accession (submit to nextflow via exec)
const { exec } = require('child_process');
app.post('/api/run', (req, res) => {
  const { accession, genome } = req.body;
  if (!accession) return res.status(400).json({ error: 'accession required' });
  const genomeArg = genome || 'GRCh38';
  const cmd = `nextflow run /app/nextflow/main.nf -profile docker \
    --accession ${accession} --genome ${genomeArg}`;
  exec(cmd, (err, stdout, stderr) => {
    if (err) return res.status(500).json({ error, stderr });
    res.json({ status: 'queued', stdout: stdout.slice(0, 1000) });
  });
});

function walkDir(dir, base = '') {
  let results = [];
  const entries = fs.readdirSync(dir, { withFileTypes: true });
  for (const entry of entries) {
    const rel = base ? `${base}/${entry.name}` : entry.name;
    if (entry.isDirectory()) {
      results = results.concat(walkDir(path.join(dir, entry.name), rel));
    } else {
      const stat = fs.statSync(path.join(dir, entry.name));
      results.push({ path: rel, size: stat.size, mtime: stat.mtime });
    }
  }
  return results;
}

app.listen(PORT, () => {
  console.log(`HU-rNMP browser on :${PORT}`);
});
