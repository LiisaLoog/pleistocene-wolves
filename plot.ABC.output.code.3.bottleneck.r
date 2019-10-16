############################################ PLOT RESULTS ########################################################

setwd("3_Bottleneck/ABC") # Provide a full path

prior <- read.table("ABC_GLM_TruncatedPrior_Obs0.txt", sep = "	", head = T)
posterior <- read.table("ABC_GLM_PosteriorEstimates_Obs0.txt", sep = "	", head = T)

pdf("3_Bottleneck_ABC_log10.mig.rate.1.0.pdf")
plot(posterior[,2],posterior[,3], type = "l", col = "red", main = "Model = 3_Bottleneck_2 ", xlab = "log10.mig.rate", ylab="Posterior density")
# lines(prior[,2], prior[,3], col = "blue")
# legend(0.5,1,legend = c("prior","posterior"), col = c("blue", "red"), lty = 1)
dev.off()

pdf("3_Bottleneck_ABC_log10.Ka.1.0.pdf")
plot(posterior[,4],posterior[,5], type = "l", col = "red", main = "Model = 3_Bottleneck_2 ", xlab = "log10.Ka", ylab="Posterior density")
# lines(prior[,4], prior[,5], col = "blue")
# legend(1,0.8,legend = c("prior","posterior"), col = c("blue", "red"), lty = 1)
dev.off()

pdf("3_Bottleneck_ABC_log10.Kb.1.0.pdf")
plot(posterior[,6],posterior[,7], type = "l", col = "red",  main = "Model = 3_Bottleneck_2 ", xlab = "log10.Kb", ylab="Posterior density")
# lines(prior[,6], prior[,7], col = "blue")
# legend(1,0.4,legend = c("prior","posterior"), col = c("blue", "red"), lty = 1)
dev.off()

pdf("3_Bottleneck_ABC_log10.Kc.1.0.pdf")
plot(posterior[,8],posterior[,9], type = "l", col = "red",  main = "Model = 3_Bottleneck_2 ", xlab = "log10.Kc", ylab="Posterior density")
# lines(prior[,8], prior[,9], col = "blue")
# legend(0.5,0.6,legend = c("prior","posterior"), col = c("blue", "red"), lty = 1)
dev.off()

setwd("../..")
