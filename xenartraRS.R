source('./b_random_skewers.R')



# save( MCMC.R.proj.25, 
#       MCMC.R.proj.28,
#       MCMC.R.proj.32, 
#       MCMC.R.proj.35, 
#       MCMC.R.proj.32.noMyr,
#       MCMC.R.proj.35.noMyr,
#       file = 'rs.proj.Rdata')

load('rs.proj.Rdata')

load('./xenartraMCMC.25.samples.Rdata')

#MCMC.R.proj.25 = R.proj(Ps)
rs.plots.25 = PlotBayesianRS(MCMC.R.proj.25, Ps)

load('./xenartraMCMC.28.samples.Rdata')

#MCMC.R.proj.28 = R.proj(Ps)
rs.plots.28 = PlotBayesianRS(MCMC.R.proj.28, Ps, 4)
graphics.off()

load('./xenartraMCMC.32.samples.Rdata')

#MCMC.R.proj.32 = R.proj(Ps)
rs.plots.32 = PlotBayesianRS(MCMC.R.proj.32, Ps, 4)
graphics.off()

load('./xenartraMCMC.35.samples.Rdata')

#MCMC.R.proj.35 = R.proj(Ps)
rs.plots.35 = PlotBayesianRS(MCMC.R.proj.35, Ps)
graphics.off()

source('xenartraSamples.R')
x$Myrmecophaga <- NULL

#Ps.35.noMyr = generateMCMCArray(35)
#MCMC.R.proj.35.noMyr = R.proj(Ps.35.noMyr)
rs.plots.35.noMyr = PlotBayesianRS(MCMC.R.proj.35.noMyr, Ps.35.noMyr)
graphics.off()

#Ps.32.noMyr = generateMCMCArray(32)
#MCMC.R.proj.32.noMyr = R.proj(Ps.32.noMyr)
rs.plots.32.noMyr = PlotBayesianRS(MCMC.R.proj.32.noMyr, Ps.32.noMyr, 4)
graphics.off()



rs.plots = list('25' = rs.plots.25,
                '28' = rs.plots.28,
                '32' = rs.plots.32,
                '35' = rs.plots.35)#,
                #'32.noMyr' = rs.plots.32.noMyr,
                #'35.noMyr' = rs.plots.35.noMyr)

save(rs.plots, file = 'rs.Rdata')


