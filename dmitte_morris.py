# %%
import numpy as np
from SALib.sample.morris import sample
from SALib.analyze.morris import analyze

# Define the model function
def model_function(X):
    return X[:, 0]**2 + np.sin(X[:, 1]) + X[:, 2]

# Problem definition
problem = {
    'num_vars': 3,  # Number of parameters
    'names': ['x1', 'x2', 'x3'],  # Names of the parameters
    'bounds': [[-10, 10], [-np.pi, np.pi], [-10, 10]]  # Ranges of the parameters
}

# Generate samples
param_values = sample(problem, N=1000, num_levels=4, optimal_trajectories=None)

# Run model (evaluate the output of the model at each set of parameters)
Y = np.array([model_function(param.reshape(1, -1)) for param in param_values])

# Perform Morris analysis
Si = analyze(problem, param_values, Y, print_to_console=False, num_levels=4, num_resamples=100)

# Print results
print("Elementary Effects (mu):", Si['mu'])
print("Standard deviations (sigma):", Si['sigma'])
print("Means of absolute deviations (mu_star):", Si['mu_star'])

# Optionally, visualize the results
# plot_horizontal_bar(Si)

# %%
