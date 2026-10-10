"""Side experiment: transfer_rates_baomi without the x100 on kbh_a2.

Checks whether the empirical x100 in transfer_rates_baomi is still needed after
the Penman-Monteith unit fix. Uses the same virtual experiments and benchmarks as
run_validation.py / analyze.py; writes validation/results/experiments/.

Run from the repository root after validation/build_inputs.py:
    python validation/exp_baomi_no_x100.py
"""
import os
import sys
import warnings

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
warnings.filterwarnings('ignore')

import analyze  # noqa: E402
import run_validation  # noqa: E402
from dmitte.transfer_equotions import transfer_rates_baomi  # noqa: E402

OUT = os.path.join(HERE, 'results', 'experiments', 'baomi_no_x100')
K_BH_A2 = run_validation.K_NAMES.index('kbh_a2')


def baomi_no_x100(*args):
    k = transfer_rates_baomi(*args).copy()
    k[:, K_BH_A2] /= 100
    return k


def main():
    run_validation.VARIANTS = {'baomi': baomi_no_x100}
    run_validation.OUT = OUT
    run_validation.main()
    acute = pd.read_csv(os.path.join(OUT, 'acute_exposures.csv'), dtype={'exposure': str})
    chronic = pd.read_csv(os.path.join(OUT, 'chronic_exposure.csv'))
    bench = pd.read_csv(os.path.join(HERE, 'benchmarks', 'literature_benchmarks.csv')).set_index('id')
    bench['id'] = bench.index
    analyze.VARIANT_ORDER = ['baomi']
    comp = analyze.comparison(acute, chronic, bench)
    comp.to_csv(os.path.join(OUT, 'benchmark_comparison.csv'), index=False)
    print(comp[['benchmark', 'crop', 'metric', 'bench_low', 'bench_high', 'model_median', 'verdict']].to_string(index=False))
    print(comp.verdict.value_counts())


if __name__ == '__main__':
    main()
