# truncated

MCMC algorithms for implementing zero-truncated Poisson or negative binomial models for count data

##### Standard zero-truncated Poisson model
- *zt.poisson.pdf.R*: Model statement and full conditional distributions
- *zt.poisson.sim.R*: Simulate data and fit model
- *zt.poisson.mcmc.R*: MCMC algorithm for parameter estimation

##### Standard zero-truncated negative binomial model
- *zt.nb.pdf.R*: Model statement and full conditional distributions
- *zt.nb.sim.R*: Simulate data and fit model
- *zt.nb.mcmc.R*: MCMC algorithm for parameter estimation

##### Three-level hierarchical zero-truncated Poisson model with coefficients that vary by 'group'
- *zt.poisson.varying.coef.3.pdf.R*: Model statement and full conditional distributions
- *zt.poisson.varying.coef.3.sim.R*: Simulate data and fit model
- *zt.poisson.varying.coef.3.mcmc.R*: MCMC algorithm for parameter estimation

##### Three-level hierarchical zero-truncated negative binomial model with coefficients that vary by 'group'
- *zt.nb.varying.coef.3.pdf.R*: Model statement and full conditional distributions
- *zt.nb.varying.coef.3.sim.R*: Simulate data and fit model
- *zt.nb.varying.coef.3.mcmc.R*: MCMC algorithm for parameter estimation
