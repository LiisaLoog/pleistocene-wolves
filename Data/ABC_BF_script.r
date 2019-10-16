ABC.results <- read.csv("ABC_Marginal_Densities.csv",header = T, sep = ",")

number.retrained = ABC.results[,3] -1
americas.filter.correction <- array(1, dim = length(ABC.results[,1]))
americas.filter.correction[c(8,9,15,16)] = 5
retrained.corrected = number.retrained * americas.filter.correction
likelihood = ABC.results[,2] * retrained.corrected
bayes.factor = likelihood / max(likelihood)


cols <- c("white","white","dodgerblue4","thistle1", "goldenrod2","violetred1","deeppink4","black","mistyrose4", "dodgerblue4","thistle1", "goldenrod2","violetred1","deeppink4","black","mistyrose4")

x <- bayes.factor
pdf("WOLVES_BAYESFACTOR_BAR_PLOT.pdf", width = 16, height = 7)
barplot(x, cex.names = 0.2, col = cols, space = 0.5)
dev.off()

ABC.BF.results =  cbind(ABC.results[,c(1,2)], retrained.corrected, likelihood,bayes.factor)

write.csv(ABC.BF.results, "ABC.BF.results.csv")
