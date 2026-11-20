import { createHash } from 'node:crypto';
import { readFileSync, writeFileSync, mkdirSync } from 'node:fs';

const features = [
  'Minimum Degree', 'Molecular Weight', 'Number of H-Bond Donors',
  'Number of Rings', 'Number of Rotatable Bonds', 'Polar Surface Area',
];
const target = 'measured log solubility in mols per litre';
const url = 'https://raw.githubusercontent.com/deepchem/deepchem/master/datasets/delaney-processed.csv';
const mean = xs => xs.reduce((sum, x) => sum + x, 0) / xs.length;
const dot = (a, b) => a.reduce((sum, x, i) => sum + x * b[i], 0);
const hash = value => createHash('sha256').update(value).digest('hex');

function csv(line) {
  return Array.from(line.matchAll(/(?:^|,)(?:"((?:[^"]|"")*)"|([^",]*))/g),
    match => match[1] === undefined ? match[2] : match[1].replaceAll('""', '"'));
}

function ridge(train) {
  const yMean = mean(train.map(row => row.y));
  const centers = features.map((_, j) => mean(train.map(row => row.x[j])));
  const scales = centers.map((center, j) =>
    Math.sqrt(mean(train.map(row => (row.x[j] - center) ** 2))) || 1);
  const transform = row => row.x.map((x, j) => (x - centers[j]) / scales[j]);
  const x = train.map(transform);
  const size = features.length;
  // Solve the regularized normal equations with pivoted elimination.
  const matrix = centers.map((_, i) => [
    ...centers.map((__, j) => mean(x.map(row => row[i] * row[j])) + (i === j ? 0.01 : 0)),
    mean(x.map((row, k) => row[i] * (train[k].y - yMean))),
  ]);
  for (let i = 0; i < size; i++) {
    let pivot = i;
    for (let j = i + 1; j < size; j++) {
      if (Math.abs(matrix[j][i]) > Math.abs(matrix[pivot][i])) pivot = j;
    }
    [matrix[i], matrix[pivot]] = [matrix[pivot], matrix[i]];
    matrix[i] = matrix[i].map(value => value / matrix[i][i]);
    for (let j = 0; j < size; j++) {
      if (j === i) continue;
      const factor = matrix[j][i];
      matrix[j] = matrix[j].map((value, k) => value - factor * matrix[i][k]);
    }
  }
  const weights = matrix.map(row => row[size]);
  const scaledTrain = x.map((values, i) => ({ values, y: train[i].y }));
  return row => {
    const values = transform(row);
    const neighbors = scaledTrain.map(item => ({
      y: item.y,
      distance: item.values.reduce((sum, xj, j) => sum + (xj - values[j]) ** 2, 0),
    })).sort((a, b) => a.distance - b.distance).slice(0, 8);
    return { mean: yMean, ridge: yMean + dot(values, weights), knn: mean(neighbors.map(item => item.y)) };
  };
}

function metrics(rows, name) {
  const errors = rows.map(row => row.predictions[name] - row.y);
  return {
    n: rows.length,
    mae: mean(errors.map(Math.abs)),
    rmse: Math.sqrt(mean(errors.map(error => error ** 2))),
  };
}

const [input, output = 'results'] = process.argv.slice(2);
const response = input ? null : await fetch(url);
if (response && !response.ok) throw new Error(`Dataset download failed: ${response.status}`);
const raw = input ? readFileSync(input, 'utf8') : await response.text();
const lines = raw.trimEnd().split(/\r?\n/);
const headers = csv(lines.shift());
const groupField = headers.includes('scaffold') ? 'scaffold' : 'smiles';
for (const field of ['Compound ID', 'smiles', target, groupField, ...features]) {
  if (!headers.includes(field)) throw new Error(`Missing column: ${field}`);
}
const rows = lines.map((line, index) => {
  const fields = csv(line);
  if (fields.length !== headers.length) throw new Error(`Malformed row ${index + 2}`);
  const get = name => fields[headers.indexOf(name)].trim();
  const y = Number(get(target));
  const x = features.map(name => Number(get(name)));
  if ([target, ...features].some(name => !get(name)) || ![y, ...x].every(Number.isFinite)) {
    throw new Error(`Invalid numeric value in row ${index + 2}`);
  }
  if (!get('smiles') || !get(groupField)) throw new Error(`Missing molecule or group in row ${index + 2}`);
  return { id: get('Compound ID'), smiles: get('smiles'), group: get(groupField), y, x };
});
const grouped = new Map();
for (const row of rows) {
  if (!grouped.has(row.group)) grouped.set(row.group, []);
  grouped.get(row.group).push(row);
}
const folds = Array.from({ length: 5 }, () => []);
const groups = [...grouped.entries()].map(([group, members]) => ({ members, key: hash('esol-v1:' + group) }));
groups.sort((a, b) => b.members.length - a.members.length || a.key.localeCompare(b.key));
for (const { members } of groups) {
  const smallest = folds.reduce((best, fold) => fold.length < best.length ? fold : best);
  smallest.push(...members);
}
if (folds.some(fold => !fold.length)) throw new Error('Five nonempty folds are required');
for (const [index, fold] of folds.entries()) {
  const predict = ridge(folds.filter((_, j) => j !== index).flat());
  for (const row of fold) {
    row.fold = index + 1;
    row.predictions = predict(row);
  }
}
const models = ['mean', 'ridge', 'knn'];
const scored = rows.map(row => ({
  ...row, absError: Math.abs(row.predictions.ridge - row.y),
}));
const summary = {
  rows: rows.length, groupField, groups: groups.length, foldSizes: folds.map(fold => fold.length),
  dataSha256: hash(raw.replaceAll('\r\n', '\n').trimEnd()), features,
  models: Object.fromEntries(models.map(name => [name, metrics(rows, name)])),
  donorSlices: [
    { label: '0 donors', values: rows.filter(row => row.x[2] === 0) },
    { label: '1 donor', values: rows.filter(row => row.x[2] === 1) },
    { label: '2+ donors', values: rows.filter(row => row.x[2] >= 2) },
  ].filter(slice => slice.values.length).map(slice => ({ label: slice.label, ridge: metrics(slice.values, 'ridge'), knn: metrics(slice.values, 'knn') })),
  topRidgeErrors: [...scored].sort((a, b) => b.absError - a.absError).slice(0, 5),
};
mkdirSync(output, { recursive: true });
writeFileSync(`${output}/summary.json`, JSON.stringify(summary, null, 2) + '\n');
const quote = value => '"' + String(value).replaceAll('"', '""') + '"';
const predictions = ['compound,smiles,group,fold,measured,mean,ridge,knn'];
for (const row of rows) {
  predictions.push([quote(row.id), quote(row.smiles), quote(row.group), row.fold, row.y,
    ...models.map(name => row.predictions[name])].join(','));
}
writeFileSync(`${output}/predictions.csv`, predictions.join('\n') + '\n');
console.log(JSON.stringify(summary, null, 2));
