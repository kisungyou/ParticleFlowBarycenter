#!/usr/bin/env python3
"""Check archived evidence without rerunning costly experiments or importing dependencies."""
import csv
import hashlib
import itertools
import json
import math
from pathlib import Path
import statistics
import sys

ROOT = Path(__file__).resolve().parent

def rows(path):
    with (ROOT / path).open(newline='') as handle:
        return list(csv.DictReader(handle))

def close(left, right, tolerance=1e-10):
    if not math.isclose(float(left), float(right), rel_tol=tolerance, abs_tol=tolerance):
        raise ValueError(f'Numerical mismatch: {left} != {right}')

def grid_check(records, fields, expected):
    actual = [tuple(float(row[field]) for field in fields) for row in records]
    if len(set(actual)) != len(actual) or set(actual) != set(itertools.product(*expected)):
        raise ValueError(f'Incomplete or duplicated grid for {fields}')

def stopped(records):
    for row in records:
        if row['status'] != 'converged':
            raise ValueError('A checked fit has an unexpected termination status')
        if not (0 <= float(row['delta']) <= 1e-8 and 0 <= float(row['residual']) <= 1e-4):
            raise ValueError('A checked fit fails the recorded objective/residual tolerance')

def main():
    manifest = json.loads((ROOT / 'archive-sha256.json').read_text())
    for name, expected in manifest.items():
        path = ROOT / name
        if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != expected:
            raise ValueError(f'Archived file is missing or changed: {name}')
    print(f'Archive integrity: {len(manifest)} files match')

    gauss = rows('revision-work/refresh/gaussian_resolution.csv')
    grid_check(gauss, ['rep', 'support', 'alpha'], [range(1, 4), [10, 25, 50, 100, 200, 500], [1, .5]])
    stopped(gauss)
    wasp = rows('revision-work/refresh/wasp_metrics.csv')
    grid_check(wasp, ['rep', 'nsplit', 'support', 'alpha'], [range(1, 4), [5, 10, 20], [50, 100, 200], [1, .5]])
    stopped(wasp)
    sensitivity = rows('revision-work/refresh/draw_sensitivity.csv')
    grid_check(sensitivity, ['rep', 'draws'], [range(1, 4), [25, 50, 100, 200]])
    stopped(sensitivity)
    print('Checked fits: 36 Gaussian, 54 WASP, 12 posterior-sensitivity configurations')

    paired = rows('revision-work/gaussian/paired_results.csv')
    methods = {'r_exact', 'pot_exact', 'sinkhorn_0.1', 'sinkhorn_0.3', 'sinkhorn_1'}
    expected_ids = {f'rep{rep:02d}_m{support:03d}' for rep in range(1, 11) for support in [50, 100, 200]}
    by_id = {}
    for row in paired:
        block = by_id.setdefault(row['id'], {})
        if row['method'] in block:
            raise ValueError('Duplicate paired solver output')
        block[row['method']] = row
        if row['status'] not in {'converged', 'iteration_limit'}:
            raise ValueError('Unexpected paired solver status')
    if set(by_id) != expected_ids:
        raise ValueError('Missing paired Gaussian problems')
    for block in by_id.values():
        if set(block) != methods:
            raise ValueError('Methods are not paired on every Gaussian problem')
        for field in ['objective', 'semi_w2', 'mean_error', 'cov_error']:
            close(block['r_exact'][field], block['pot_exact'][field], tolerance=1e-8)
        if float(block['pot_exact']['pot_native_max_difference']) > 1e-8:
            raise ValueError('POT native implementation check exceeds tolerance')
    capped = sum(row['status'] == 'iteration_limit' for row in paired)
    if capped != 7:
        raise ValueError('Unexpected number of archived iteration-limit exits')
    print('Paired Gaussian: 30 shared problems, 150 outputs, 7 capped fits explicitly retained')

    medoid = rows('revision-work/mnist/split_metrics.csv')
    summary = json.loads((ROOT / 'revision-work/mnist/summary.json').read_text())
    if len(medoid) != 3 or summary['completed_splits'] != 3:
        raise ValueError('MNIST summary must include three balanced splits')
    for field in ['Accuracy', 'MacroPrecision', 'MacroRecall', 'MacroF1', 'classification_sec']:
        values = [float(row[field]) for row in medoid]
        close(statistics.mean(values), summary[field + '_mean'])
        close(statistics.stdev(values), summary[field + '_sd'])
    close(summary['Accuracy_mean'], .5968)
    print('MNIST restricted medoid: three-split metrics reproduce the reported summary')

    dvq = rows('revision-work/dvq/paired_metrics.csv')
    pair_groups = {}
    for row in dvq:
        key = row['data'], int(row['S']), int(row['repeat_id'])
        block = pair_groups.setdefault(key, {})
        if row['method'] in block:
            raise ValueError('Duplicate DVQ pair member')
        block[row['method']] = row
    if set(pair_groups) != set(itertools.product(['pbmc', 'fashion'], [5, 10], range(1, 4))):
        raise ValueError('Incomplete DVQ paired grid')
    differences = {}
    for (dataset, subsets, repeat), block in pair_groups.items():
        if set(block) != {'DVQ', 'Summary k-means'}:
            raise ValueError('DVQ comparator is missing')
        left, right = block['DVQ'], block['Summary k-means']
        for field in ['seed', 'n', 'd', 'k', 'transmitted_coordinates']:
            if left[field] != right[field]:
                raise ValueError(f'DVQ pair differs in {field}')
        if int(left['transmitted_coordinates']) != subsets * 10 * 20:
            raise ValueError('DVQ summary budget differs from the manuscript')
        if left['converged'] != 'TRUE' or not 3 <= int(left['iterations']) <= 5:
            raise ValueError('DVQ final stopping status differs from the manuscript')
        record = {
            'distortion_relative_difference_pct': 100 * (float(left['distortion']) / float(right['distortion']) - 1),
            'ARI_difference': float(left['ARI']) - float(right['ARI']),
            'NMI_difference': float(left['NMI']) - float(right['NMI']),
        }
        differences.setdefault((dataset, subsets), []).append(record)
    for summary_row in rows('revision-work/dvq/paired_difference_summary.csv'):
        values = differences[(summary_row['data'], int(summary_row['S']))]
        for field in values[0]:
            close(statistics.mean(row[field] for row in values), summary_row[field + '_mean'])
            close(statistics.stdev(row[field] for row in values), summary_row[field + '_sd'])
    print('DVQ: 12 matched-budget pairs reproduce the reported differences')

    stems = ['fig-sim-gauss-2', 'fig-sim-gauss-3', 'fig-sim-gauss-4',
             'fig-sim-normal-1', 'fig-sim-normal-2', 'fig-sim-normal-3',
             'fig-real-digits-1', 'fig-real-digits-2',
             'fig-real-clustering-1', 'fig-real-clustering-2']
    for stem in stems:
        for suffix in ['pdf', 'pgf']:
            if not (ROOT / 'manuscript-figures/published' / f'{stem}.{suffix}').is_file():
                raise ValueError(f'Missing manuscript figure {stem}.{suffix}')
    print('Figures: all ten preserved PDF/PGF pairs are present')
    print('PASS. These checks verify archived records, not unavailable historical experiments.')

if __name__ == '__main__':
    try:
        main()
    except (OSError, ValueError, KeyError) as error:
        print(f'FAIL: {error}', file=sys.stderr)
        sys.exit(1)
