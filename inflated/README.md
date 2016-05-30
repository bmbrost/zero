# inflated

MCMC algorithms for implementing zero-inflated mixture models for count data

#### Contents

- *zi.poisson*: Standard zero-inflated Poisson model
- *zi.poisson.varying.coef.3*: Three-level hierarchical zero-inflated Poisson model with coefficients that vary by 'group'
- *zi.poisson.varying.coef.4*: Four-level hierarchical zero-inflated Poisson model with coefficients that vary by 'group' and 'subgroup'

###### Naming conventions
- **.pdf.R*: Model statement and full conditional distributions
- **.sim.R*: Simulate data and fit model
- **.mcmc.R*: MCMC algorithm for parameter estimation
