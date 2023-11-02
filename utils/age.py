import numpy as np
from scipy.stats import chi2

# t squared test from here https://en.wikipedia.org/wiki/Hotelling%27s_T-squared_distribution
def hotteling_t2_test(xnew, mu, S):
    "Compute Mahalanobis distance and corresponding p-value for a given point."
    p = mu.shape[0]
    d = (xnew - mu) @ np.linalg.inv(S) @ (xnew - mu) #squared mahalanobis distance
    pval = 1 - chi2.cdf(d, p)
    return d, pval

#tmp 
def compute_age_univariate(x, means, cov):
    mu = means[1] + cov[0,1] / cov[0,0] * (x - means[0]) 
    s = np.sqrt(cov[1,1] - cov[0,1]**2 / cov[0,0]) 
    return mu, s

#compute conditional mean and std for age, i.e. mu, sigma from P(y|X)
#or inference most probable age for the corresponding x and its prediction uncertainty (epistemic)
def inference_age_multivariate(x, means, cov):
    mu_x = means[:-1]
    mu_y = means[-1]
    sx = cov[:-1, :-1]
    sy = cov[-1, -1] 
    sxy = cov[:-1, -1]
    y_mean_new = mu_y + sxy @ np.linalg.pinv(sx) @ (x - mu_x)
    s_y_new = np.sqrt(sy - sxy @ np.linalg.pinv(sx) @ sxy.T) 
    return y_mean_new, s_y_new