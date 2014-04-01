source("xenartraSamples.R")

Ps = generateMCMCArray(25)
Ps = Ps[,,,901:1000]
save(Ps, file = 'xenartraMCMC.25.samples.Rdata')

Ps = generateMCMCArray(28)
Ps = Ps[,,,901:1000]
save(Ps, file = 'xenartraMCMC.28.samples.Rdata')

Ps = generateMCMCArray(32)
Ps = Ps[,,,901:1000]
save(Ps, file = 'xenartraMCMC.32.samples.Rdata')

Ps = generateMCMCArray(35)
Ps = Ps[,,,901:1000]
save(Ps, file = 'xenartraMCMC.35.samples.Rdata')
