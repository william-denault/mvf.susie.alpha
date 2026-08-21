import json
import math
from collections import defaultdict

import numpy as np


RNG = np.random.default_rng(20260821)
T = 256
P = 200
REPS = 40


def col_scale(x):
    x = np.asarray(x, dtype=float)
    center = x.mean(axis=0)
    scale = x.std(axis=0, ddof=1)
    scale[scale == 0] = 1.0
    return (x - center) / scale


def haar_rows(x):
    current = np.asarray(x, dtype=float)
    details = []
    inv_sqrt_2 = 1.0 / math.sqrt(2.0)
    while current.shape[1] > 1:
        left = current[:, 0::2]
        right = current[:, 1::2]
        details.append((left - right) * inv_sqrt_2)
        current = (left + right) * inv_sqrt_2
    return np.concatenate(details + [current], axis=1)


def pipeline(raw):
    return col_scale(haar_rows(col_scale(raw)))


def genotype_matrix(n):
    maf = RNG.uniform(0.1, 0.4, size=P)
    return col_scale(RNG.binomial(2, maf, size=(n, P)))


def fixed_mixture_log_bf(z):
    # A deliberately conservative fixed prior: 99% spike at zero and three
    # equal non-null components with prior-variance/SE-variance ratios 1,4,16.
    weights = np.array([0.99, 0.01 / 3, 0.01 / 3, 0.01 / 3])
    ratios = np.array([0.0, 1.0, 4.0, 16.0])
    log_terms = []
    for weight, ratio in zip(weights, ratios):
        log_bf = -0.5 * np.log1p(ratio) + 0.5 * z * z * ratio / (1.0 + ratio)
        log_terms.append(math.log(weight) + log_bf)
    stacked = np.stack(log_terms, axis=0)
    maxima = stacked.max(axis=0)
    log_mix = maxima + np.log(np.exp(stacked - maxima).sum(axis=0))
    return log_mix.sum(axis=1)


def pairwise_square_dependence(w):
    sample = min(96, w.shape[1])
    columns = RNG.choice(w.shape[1], size=sample, replace=False)
    square = col_scale(w[:, columns] ** 2)
    corr = square.T @ square / (square.shape[0] - 1)
    off_diag = corr[np.triu_indices(sample, k=1)]
    return float(np.nanmean(np.abs(off_diag)))


def excess_kurtosis(w):
    centered = w - w.mean(axis=0)
    variance = np.mean(centered**2, axis=0)
    keep = variance > 0
    kurtosis = np.mean(centered[:, keep] ** 4, axis=0) / variance[keep] ** 2 - 3
    return float(np.median(kurtosis))


def simulate_scenario(kind, n, modalities, rate=None):
    metrics = defaultdict(list)
    for _ in range(REPS):
        x = genotype_matrix(n)
        d = n - 1
        joint_lbf = np.zeros(P)
        first_w = None
        first_raw = None
        for _m in range(modalities):
            if kind == "gaussian":
                raw = RNG.normal(size=(n, T))
            else:
                raw = RNG.poisson(rate, size=(n, T))
            w = pipeline(raw)
            z = x.T @ w / math.sqrt(d)
            joint_lbf += fixed_mixture_log_bf(z)
            if first_w is None:
                first_w = w
                first_raw = raw

        ordered = np.sort(joint_lbf)
        metrics["max_joint_logbf"].append(float(ordered[-1]))
        metrics["max_minus_second"].append(float(ordered[-1] - ordered[-2]))
        metrics["median_wc_excess_kurtosis"].append(excess_kurtosis(first_w))
        metrics["mean_abs_corr_squared_wc"].append(pairwise_square_dependence(first_w))
        metrics["fraction_zero_rows"].append(float(np.mean(first_raw.sum(axis=1) == 0)))
        metrics["fraction_zero_columns"].append(float(np.mean(first_raw.sum(axis=0) == 0)))

    summary = {}
    for name, values in metrics.items():
        array = np.asarray(values)
        summary[name] = {
            "mean": float(array.mean()),
            "median": float(np.median(array)),
            "q95": float(np.quantile(array, 0.95)),
        }
    maximum = np.asarray(metrics["max_joint_logbf"])
    summary["prob_max_logbf_gt_log10"] = float(np.mean(maximum > math.log(10)))
    summary["prob_max_logbf_gt_log100"] = float(np.mean(maximum > math.log(100)))
    return summary


scenarios = {
    "gaussian_n84_m1": ("gaussian", 84, 1, None),
    "gaussian_n84_m6": ("gaussian", 84, 6, None),
    "poisson_5_n84_m6": ("poisson", 84, 6, 5.0),
    "poisson_0.135_n84_m1": ("poisson", 84, 1, math.exp(-2)),
    "poisson_0.135_n84_m6": ("poisson", 84, 6, math.exp(-2)),
    "poisson_0.0025_n84_m1": ("poisson", 84, 1, math.exp(-6)),
    "poisson_0.0025_n84_m6": ("poisson", 84, 6, math.exp(-6)),
    "gaussian_n300_m1": ("gaussian", 300, 1, None),
    "poisson_0.135_n300_m1": ("poisson", 300, 1, math.exp(-2)),
}

output = {}
for label, arguments in scenarios.items():
    output[label] = simulate_scenario(*arguments)
    print(label, json.dumps(output[label], sort_keys=True), flush=True)
