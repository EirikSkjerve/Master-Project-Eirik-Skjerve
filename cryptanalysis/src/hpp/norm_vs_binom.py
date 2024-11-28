import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import binom, norm


def compare_distributions(n=16, p=0.5, mu=0, sigma=2, size=1000000):
    """
    Compares a binomial distribution and a rounded normal distribution.

    Parameters:
    - n (int): Number of trials for the binomial distribution.
    - p (float): Probability of success for the binomial distribution.
    - mu (float): Mean of the normal distribution.
    - sigma (float): Standard deviation of the normal distribution.
    - size (int): Number of samples for the normal distribution.
    """

    n = 4*sigma**2

    # Generate Binomial Distribution
    x_binom = np.arange(0, n+1)
    y_binom = binom.pmf(x_binom, n, p)

    # Generate Rounded Normal Distribution
    normal_samples = np.random.normal(mu, sigma, size)
    rounded_samples = np.round(normal_samples).astype(int)
    unique, counts = np.unique(rounded_samples, return_counts=True)
    y_normal = counts / size  # Convert counts to probabilities

    # Align x values for visualization
    x_normal = unique

    # Plot
    plt.figure(figsize=(10, 6))
    plt.bar(x_binom - n*p, y_binom, width=0.4, label='Binomial Distribution (n={}, p={})'.format(n, p), alpha=0.7)
    plt.bar(x_normal, y_normal, width=0.4, label='Rounded Normal (µ={}, sigma={})'.format(mu, sigma), alpha=0.7)
    plt.xlabel('Value')
    plt.ylabel('Probability')
    plt.title('Comparison of Binomial and Rounded Normal Distributions')
    plt.legend()
    plt.show()


# Example use of the function
compare_distributions()
