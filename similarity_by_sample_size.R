library(Morphometrics)
library(gdata)
load('./Rdatas/maindata.RData')

cov.matrices = llply(x, function(x) x$vcv)
mats = llply(names(cov.matrices[[1]]), function(x) llply(cov.matrices, function(y) y[[x]]))
names(mats) = names(cov.matrices[[1]])
mats = llply(mats, function(x) x[!is.na(x)])
comparisons_RS = llply(mats, function(x) RandomSkewers(x)[[1]])
comparisons_Krz = llply(mats, function(x) KrzCor(x))
mats = llply(mats, function(x) llply(x,function(y) cov2cor(y)))
comparisons_Mantel = llply(mats, function(x) MantelCor(x, iterations = 1)[[1]])
comparisons_krz_corr = llply(mats, function(x) KrzCor(x))
comparisons_RS_corr = llply(mats, function(x) RandomSkewers(x)[[1]])
sample.size = ldply(x, function(x) as.data.frame(x$df))
CrossSamples = function(name){
    samples = sample.size[[name]]
    ids = sample.size[['.id']][!is.na(samples)]
    samples = samples[!is.na(samples)]
    m = length(samples)
    cross = matrix(NA, m, m)
    for(i in 1:m)
        for(j in 1:m)
            cross[i,j] = 1/mean(1/c(samples[i], samples[j]))
    dimnames(cross) = list(ids, ids)
    return(cross)
}
cross = llply(names(mats), CrossSamples)
names(cross) = names(mats)
DfComps = function(name, comps, cross){
    comps = lowerTriangle(comps)
    cross = lowerTriangle(cross)
    df = data.frame(comparisons = comps, samples = cross)
    df$trait = name
    return(df)
}
dat_RS       = ldply(names(mats), function(x) DfComps(x, comparisons_RS[[x]], cross[[x]]))
dat_Krz      = ldply(names(mats), function(x) DfComps(x, comparisons_Krz[[x]], cross[[x]]))
dat_Mantel   = ldply(names(mats), function(x) DfComps(x, comparisons_Mantel[[x]], cross[[x]]))
dat_Krz_corr = ldply(names(mats), function(x) DfComps(x, comparisons_krz_corr[[x]], cross[[x]]))
dat_RS_corr  = ldply(names(mats), function(x) DfComps(x, comparisons_RS_corr[[x]], cross[[x]]))

library(ggplot2)
ggplot(dat_RS, aes(samples, comparisons)) + geom_point() + theme_bw() + facet_wrap(~trait, ncol = 2) + labs(x = 'Sample size', y = 'Structural similarity')
ggsave("~/Dropbox/labbio/articles/VCV/Aguirre Stuff/RS_SampleSize.tiff")
ggplot(dat_RS_corr, aes(samples, comparisons)) + geom_point() + theme_bw() + facet_wrap(~trait, ncol = 2) + labs(x = 'Sample size', y = 'Structural similarity')
ggsave("~/Dropbox/labbio/articles/VCV/Aguirre Stuff/RS_corr_SampleSize.tiff")
ggplot(dat_Krz, aes(samples, comparisons)) + geom_point() + theme_bw() + facet_wrap(~trait, ncol = 2) + labs(x = 'Sample size', y = 'Structural similarity')
ggsave("~/Dropbox/labbio/articles/VCV/Aguirre Stuff/Krz_SampleSize.tiff")
ggplot(dat_Mantel, aes(samples, comparisons)) + geom_point() + theme_bw() + facet_wrap(~trait, ncol = 2) + labs(x = 'Sample size', y = 'Structural similarity')
ggsave("~/Dropbox/labbio/articles/VCV/Aguirre Stuff/Mantel_SampleSize.tiff")
ggplot(dat_Krz_corr, aes(samples, comparisons)) + geom_point() + theme_bw() + facet_wrap(~trait, ncol = 2) + labs(x = 'Sample size', y = 'Structural similarity')
ggsave("~/Dropbox/labbio/articles/VCV/Aguirre Stuff/KrzCor_SampleSize.tiff")
