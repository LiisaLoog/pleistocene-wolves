
############################################ PLOT RESULTS ########################################################

for(n in 1:7){
	print(n)
	setwd(paste("2_Expansion/ABC", n,  sep = ""))

	prior <- read.table("ABC_GLM_TruncatedPrior_Obs0.txt", sep = "	", head = T)
	posterior <- read.table("ABC_GLM_PosteriorEstimates_Obs0.txt", sep = "	", head = T)

	pdf("2_Expansion_ABC_log10.mig.rate.1.0.pdf")
	plot(posterior[,2],posterior[,3], type = "l", col = "red", main = paste("Model = 2_Expansion", n, sep = "_"), xlab = "log10.mig.rate", ylab="Posterior density")
	# lines(prior[,2], prior[,3], col = "blue")
	# legend(0.5,1,legend = c("prior","posterior"), col = c("blue", "red"), lty = 1)
	dev.off()


	pdf("2_Expansion_ABC_log10.K.1.0.pdf")
	plot(posterior[,4],posterior[,5], type = "l", col = "red", main = paste("Model = 2_Expansion", n, sep = "_"), xlab = "log10.K", ylab="Posterior density")
	# lines(prior[,4], prior[,5], col = "blue")
	# legend(-1,0.8,legend = c("prior","posterior"), col = c("blue", "red"), lty = 1)
	dev.off()

	pdf("2_Expansion_ABC_x.bottel.1.0.pdf")
	plot(posterior[,6],posterior[,7], type = "l", col = "red",  main = paste("Model = 2_Expansion", n, sep = "_"), xlab = "x.bottle", ylab="Posterior density")
	# lines(prior[,6], prior[,7], col = "blue")
	# legend(0.4,1.4,legend = c("prior","posterior"), col = c("blue", "red"), lty = 1)
	dev.off()

	pdf("2_Expansion_ABC_exp.start.1.0.pdf")
	plot(posterior[,8],posterior[,9], type = "l", col = "red",  main = paste("Model = 2_Expansion", n, sep = "_"), xlab = "exp.start", ylim = c(0,0.11), ylab="Posterior density")
	# lines(prior[,8], prior[,9], col = "blue")
	# legend(10,0.08,legend = c("prior","posterior"), col = c("blue", "red"), lty = 1)
	dev.off()

	pdf("2_Expansion_ABC_exp.interval.1.0.pdf")
	plot(posterior[,10],posterior[,11], type = "l", col = "red",  main = paste("Model = 2_Expansion", n, sep = "_"), xlab = "exp.interval", ylab="Posterior density")
	# lines(prior[,10], prior[,11], col = "blue")
	# legend(40,0.01,legend = c("prior","posterior"), col = c("blue", "red"), lty = 1)
	dev.off()

	setwd("../..")
}
