#!/usr/bin/python3

import pytest
import numpy as np

def test_nans_in_same_place(testCOB):
    norm_expr = testCOB.expr(raw=False)
    raw_expr = testCOB.expr(raw=True).ix[norm_expr.index,norm_expr.columns]
    assert all(np.isnan(norm_expr) == np.isnan(raw_expr))
    assert all(np.isnan(raw_expr) == np.isnan(norm_expr))

def test_inplace_nansort(testCOB):
    x = np.random.rand(50000)
    for i in  np.random.randint(0,50000,500):
        x[i] = np.nan
    sorted_x = testCOB.inplace_nansort(x)
    assert all(np.isnan(x) == np.isnan(sorted_x))
