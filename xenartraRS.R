source('./b_random_skewers.R')

load("./xenartraMCMC.samples.Rdata")

MCMC.R.proj = list()
MCMC.R.proj[['25']] = R.proj(Ps[['25']])
MCMC.R.proj[['28']] = R.proj(Ps[['28']])
MCMC.R.proj[['32']] = R.proj(Ps[['32']])
MCMC.R.proj[['35']] = R.proj(Ps[['35']])
save(MCMC.R.proj, file = 'rs.proj.Rdata')
load('rs.proj.Rdata')

rs.plots = list()
rs.plots[['25']] = PlotBayesianRS(MCMC.R.proj[['25']], Ps[['25']])
rs.plots[['29']] = PlotBayesianRS(MCMC.R.proj[['28']], Ps[['28']], 4)
rs.plots[['32']] = PlotBayesianRS(MCMC.R.proj[['32']], Ps[['32']], 4)
rs.plots[['35']] = PlotBayesianRS(MCMC.R.proj[['35']], Ps[['35']])
save(rs.plots, file = 'rs.Rdata')


