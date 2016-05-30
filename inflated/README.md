# inflated

MCMC algorithms for implementing zero-inflated mixture models for count data

##### Standard zero-inflated Poisson model
- *zi.poisson.pdf.R*: Model statement and full conditional distributions
- *zi.poisson.sim.R*: Simulate data and fit model
- *zi.poisson.mcmc.R*: MCMC algorithm for parameter estimation

##### Three-level hierarchical zero-inflated Poisson model with coefficients that vary by 'group'
- *zi.poisson.varying.coef.3.pdf.R*: Model statement and full conditional distributions
- *zi.poisson.varying.coef.3.sim.R*: Simulate data and fit model
- *zi.poisson.varying.coef.3.mcmc.R*: MCMC algorithm for parameter estimation

##### Four-level hierarchical zero-inflated Poisson model with coefficients that vary by 'group' and 'subgroup'
- *zi.poisson.varying.coef.4.pdf.R*: Model statement and full conditional distributions
- *zi.poisson.varying.coef.4.sim.R*: Simulate data and fit model
- *zi.poisson.varying.coef.4.mcmc.R*: MCMC algorithm for parameter estimation