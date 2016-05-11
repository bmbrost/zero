# inflated

MCMC algorithms for implementing zero-inflated mixture models for count data

##### Standard zero-inflated Poisson model
- *zi.poisson.pdf.R*: Model statement and full conditional distributions
- *zi.poisson.sim.R*: Simulate data and fit model
- *zi.poisson.mcmc.R*: MCMC algorithm for parameter estimation

##### Three-level hierarchical zero-inflated Poisson model with coefficients that vary by 'group'
- *zi.poisson.varying.coef.3.pdf.R*: Model statement and full conditional distributions
- *binomial.logit.glm.sim.R*: Simulate data and fit model
- *binomial.logit.glm.mcmc.R*: MCMC algorithm for parameter estimation