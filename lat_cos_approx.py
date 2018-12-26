# Function to sample from a cosine distribution approximated by a normal distribution for speed
def lat_cos_approx(n_samples):
	import numpy as np
	from scipy import stats
	lower, upper = -np.pi/2.0, np.pi/2.0
	mu, sigma = 0, 0.75*np.sqrt( (np.pi**2/3.0) - 2.0)
	X = stats.truncnorm(
	    (lower - mu) / sigma, (upper - mu) / sigma, loc=mu, scale=sigma)
	return X.rvs(n_samples)
