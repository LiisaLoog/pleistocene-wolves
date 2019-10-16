############################################ PLOT RESULTS #######################################################
setwd("/1_Static/ABC") # Provide a full path

prior <- read.table("ABC_GLM_TruncatedPrior_Obs0.txt", sep = "	", head = T)
posterior <- read.table("ABC_GLM_PosteriorEstimates_Obs0.txt", sep = "	", head = T)

pdf("ABC_log10.mig.rate.1.0.pdf")
plot(posterior[,2],posterior[,3], type = "l", col = "red", main = "Model = 1_Static ", xlab = "log10.mig.rate", ylab="Posterior density")
# lines(prior[,2], prior[,3], col = "blue")
# legend(0.5,0.8,legend = c("prior","posterior"), col = c("blue", "red"), lty = 1)
dev.off()


pdf("ABC_log10.K.1.0.pdf")
plot(posterior[,4],posterior[,5], type = "l", col = "red",  main = "Model = 1_Static ", xlab = "log10.K", ylab="Posterior density")
# lines(prior[,4], prior[,5], col = "blue")
# legend(-1.0,0.8,legend = c("prior","posterior"), col = c("blue", "red"), lty = 1)
dev.off()

setwd("../..")
