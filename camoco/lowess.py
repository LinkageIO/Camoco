"""
This module implements the Lowess function for nonparametric regression.

Functions:
lowess        Fit a smooth nonparametric regression curve to a scatterplot.

For more information, see

William S. Cleveland: "Robust locally weighted regression and smoothing
scatterplots", Journal of the American Statistical Association, December 1979,
volume 74, number 368, pp. 829-836.

William S. Cleveland and Susan J. Devlin: "Locally weighted regression: An
approach to regression analysis by local fitting", Journal of the American
Statistical Association, September 1988, volume 83, number 403, pp. 596-610.
"""

from math import ceil
import numpy as np
from scipy import linalg


def confband(xd,yd,a,b,conf=0.95,x=None):
    """
    Calculates the confidence band of the linear regression model at the desired confidence
    level, using analytical methods. The 2sigma confidence interval is 95% sure to contain 
    the best-fit regression line. This is not the same as saying it will contain 95% of 
    the data points.
    
    Arguments:
    - conf: desired confidence level, by default 0.95 (2 sigma)
    - xd,yd: data arrays
    - a,b: linear fit parameters as in y=ax+b
    - x: (optional) array with x values to calculate the confidence band. If none is provided, will
    by default generate 100 points in the original x-range of the data.
    
    Returns:
    Sequence (lcb,ucb,x) with the arrays holding the lower and upper confidence bands 
    corresponding to the [input] x array.
    
    Usage:
    >>> lcb,ucb,x=nemmen.confband(all.kp,all.lg,a,b,conf=0.95)
    calculates the confidence bands for the given input arrays
    
    >>> pylab.fill_between(x, lcb, ucb, alpha=0.3, facecolor='gray')
    plots a shaded area containing the confidence band
    
    References:
    1. http://en.wikipedia.org/wiki/Simple_linear_regression, see Section Confidence intervals
    2. http://www.weibull.com/DOEWeb/confidence_intervals_in_simple_linear_regression.htm
    
    Author: Rodrigo Nemmen
    v1 Dec. 2011
    v2 Jun. 2012: corrected bug in computing dy
    """
    alpha=1.-conf   # significance
    n=xd.size   # data sample size

    if x==None: x=numpy.linspace(xd.min(),xd.max(),100)

    # Predicted values (best-fit model)
    y=a*x+b

    # Auxiliary definitions
    sd=scatterfit(xd,yd,a,b)    # Scatter of data about the model
    sxd=numpy.sum((xd-xd.mean())**2)
    sx=(x-xd.mean())**2 # array

    # Quantile of Student's t distribution for p=1-alpha/2
    q=scipy.stats.t.ppf(1.-alpha/2.,n-2)

    # Confidence band
    dy=q*sd*numpy.sqrt( 1./n + sx/sxd )
    ucb=y+dy    # Upper confidence band
    lcb=y-dy    # Lower confidence band

    return lcb,ucb,x

def lowess(x, y, f=2./3., iter=3):
    """
    lowess(x, y, f=2./3., iter=3) -> yest

    Lowess smoother: Robust locally weighted regression.
    The lowess function fits a nonparametric regression curve to a scatterplot.
    The arrays x and y contain an equal number of elements; each pair
    (x[i], y[i]) defines a data point in the scatterplot. The function returns
    the estimated (smooth) values of y.

    The smoothing span is given by f. A larger value for f will result in a
    smoother curve. The number of robustifying iterations is given by iter. The
    function will run faster with a smaller number of iterations.
    """
    n = len(x)
    r = int(ceil(f*n))
    h = [np.sort(np.abs(x - x[i]))[r] for i in range(n)]
    w = np.clip(np.abs((x[:,None] - x[None,:]) / h), 0.0, 1.0)
    w = (1 - w**3)**3
    yest = np.zeros(n)
    delta = np.ones(n)
    for iteration in range(iter):
        for i in range(n):
            weights = delta * w[:,i]
            b = np.array([np.sum(weights*y), np.sum(weights*y*x)])
            A = np.array([[np.sum(weights), np.sum(weights*x)],
                   [np.sum(weights*x), np.sum(weights*x*x)]])
            beta = linalg.solve(A, b)
            yest[i] = beta[0] + beta[1]*x[i]

        residuals = y - yest
        s = np.median(np.abs(residuals))
        delta = np.clip(residuals / (6.0 * s), -1, 1)
        delta = (1 - delta**2)**2

    return yest

if __name__ == '__main__':
    import math
    n = 100
    x = np.linspace(0, 2 * math.pi, n)
    y = np.sin(x) + 0.3*np.random.randn(n)

    f = 0.25
    yest = lowess(x, y, f=f, iter=3)

    import pylab as pl
    pl.clf()
    pl.plot(x, y, label='y noisy')
    pl.plot(x, yest, label='y pred')
    pl.legend()
    pl.show()
