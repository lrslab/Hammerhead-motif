#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 14/3/2025 9:21 pm
# @Author  : Runsheng
# @File    : stat.py


import numpy as np
from scipy import stats


def assess_motif_significance(motifs, pvalue_threshold=0.05):
    """
    Assess statistical significance of found motifs.

    Parameters:
    -----------
    motifs : list
        List of motif dictionaries from MEME
    pvalue_threshold : float
        P-value threshold for significance

    Returns:
    --------
    list
        List of significant motifs
    """
    significant_motifs = []

    for motif in motifs:
        # Calculate information content per position
        ic_per_pos = np.sum(motif['pwm'] * np.log2(motif['pwm'] + 1e-10), axis=1)

        # Perform statistical tests
        # 1. Chi-square test for position-specific nucleotide frequencies
        chi_square_pvals = []
        expected = np.array([0.25, 0.25, 0.25, 0.25])  # Assuming uniform background

        for pos in motif['pwm']:
            observed = pos * motif['nsites']
            expected_counts = expected * motif['nsites']
            chi2, pval = stats.chisquare(observed, expected_counts)
            chi_square_pvals.append(pval)

        # 2. Calculate overall motif significance
        combined_pval = stats.combine_pvalues(chi_square_pvals)[1]

        if combined_pval < pvalue_threshold:
            motif['combined_pval'] = combined_pval