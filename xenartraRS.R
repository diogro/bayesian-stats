library(MCMCglmm)
source('./b_random_skewers.R')

load("./xenartraMCMCsamples.Rdata")

MCMC.R.proj = c('25', '28', '32', '35')
MCMC.R.proj = alply(MCMC.R.proj, 1, function(traits) R.proj(Ps[[traits]]), .parallel = TRUE)
names(MCMC.R.proj) <- c('25', '28', '32', '35')
save(MCMC.R.proj, file = 'rs.proj.Rdata')
#load('rs.proj.Rdata')

rs.plots = list()
rs.plots[['25']] = PlotBayesianRS(MCMC.R.proj[['25']], Ps[['25']])
rs.plots[['28']] = PlotBayesianRS(MCMC.R.proj[['28']], Ps[['28']], 4)
rs.plots[['32']] = PlotBayesianRS(MCMC.R.proj[['32']], Ps[['32']], 4)
rs.plots[['35']] = PlotBayesianRS(MCMC.R.proj[['35']], Ps[['35']])
save(rs.plots, file = 'rs.Rdata')
#load('rs.Rdata')

for(i in names(rs.plots)) ggsave(paste0("~/Desktop/RS.plot.", i,"D.tiff"), rs.plots[[i]], width = 40, height = 20, units = "cm")
