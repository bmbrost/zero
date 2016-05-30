# truncated

MCMC algorithms for implementing zero-truncated Poisson and negative binomial models for count data

#### Contents

- *zt.poisson*: Standard zero-truncated Poisson model
- *zt.nb*: Standard zero-truncated negative binomial model
- *zt.poisson.varying.coef.3*: Three-level hierarchical zero-truncated Poisson model with coefficients that vary by 'group'
- *zt.nb.varying.coef.3*: Three-level hierarchical zero-truncated negative binomial model with coefficients that vary by 'group'

###### Naming conventions
- **.pdf.R*: Model statement and full conditional distributions
- **.sim.R*: Simulate data and fit model
- **.mcmc.R*: MCMC algorithm for parameter estimation
