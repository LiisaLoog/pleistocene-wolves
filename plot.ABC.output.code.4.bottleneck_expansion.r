
############################################ PLOT RESULTS #######################################

for(n in 1:7){
	print(n)
	setwd(paste("4_Bottleneck_Expansion_2/ABC", n,sep = "")) # Provide a full path

	prior <- read.table("ABC_GLM_TruncatedPrior_Obs0.txt", sep = "	", head = T)
	posterior <- read.table("ABC_GLM_PosteriorEstimates_Obs0.txt", sep = "	", head = T)

	pdf(paste(n,"4_Bottleneck_Expansion_ABC_log10.mig.rate.1.0.pdf",sep = "_"))
	plot(posterior[,2],posterior[,3], type = "l", col = "red", main = paste("Model = ", n," 4_Bottleneck_Expansion " ,sep = "-"), xlab = "log10.mig.rate", ylab="Posterior density",  xlim = c(min(c(prior[,2], posterior[,2])), max(c(prior[,2], posterior[,2]))), ylim = c(min(c(prior[,3], posterior[,3])), max(c(prior[,3], posterior[,3]))))
	# lines(prior[,2], prior[,3], col = "blue")
	# legend(min(c(prior[,2], posterior[,2])),max(c(prior[,3], posterior[,3])),legend = c("prior","posterior"), col = c("blue", "red"), lty = 1)
	dev.off()

	pdf(paste(n,"4_Bottleneck_Expansion_ABC_log10.Ka.1.0.pdf",sep = "_"))
	plot(posterior[,4],posterior[,5], type = "l", col = "red", main = paste("Model = ", n," 4_Bottleneck_Expansion " ,sep = "-"), xlab = "log10.Ka", ylab="Posterior density", xlim = c(min(c(prior[,4], posterior[,4])), max(c(prior[,4], posterior[,4]))), ylim = c(min(c(prior[,5], posterior[,5])), max(c(prior[,5], posterior[,5]))))
	# lines(prior[,4], prior[,5], col = "blue")
	# legend(min(c(prior[,4], posterior[,4])),max(c(prior[,5], posterior[,5])),legend = c("prior","posterior"), col = c("blue", "red"), lty = 1)
	dev.off()


	pdf(paste(n,"4_Bottleneck_Expansion_ABC_log10.Kb.1.0.pdf", sep = "_"))
	plot(posterior[,6],posterior[,7], type = "l", col = "red", main = paste("Model = ", n," 4_Bottleneck_Expansion " ,sep = "-"), xlab = "log10.Kb", ylab="Posterior density", xlim = c(min(c(prior[,6], posterior[,6])), max(c(prior[,6], posterior[,6]))), ylim = c(min(c(prior[,7], posterior[,7])), max(c(prior[,7], posterior[,7]))))
	# lines(prior[,6], prior[,7], col = "blue")
	# legend(min(c(prior[,6], posterior[,6])),max(c(prior[,7], posterior[,7])),legend = c("prior","posterior"), col = c("blue", "red"), lty = 1)
	dev.off()

	pdf(paste(n,"4_Bottleneck_Expansion_ABC_log10.Kc.1.0.pdf",sep = "_"))
	plot(posterior[,8],posterior[,9], type = "l", col = "red", main = paste("Model = ", n," 4_Bottleneck_Expansion " ,sep = "-"), xlab = "log10.Kc", ylab="Posterior density", xlim = c(min(c(prior[,8], posterior[,8])), max(c(prior[,8], posterior[,8]))), ylim = c(min(c(prior[,9], posterior[,9])), max(c(prior[,9], posterior[,9]))))
	# lines(prior[,8], prior[,9], col = "blue")
	# legend(min(c(prior[,8], posterior[,8])),max(c(prior[,9], posterior[,9])),legend = c("prior","posterior"), col = c("blue", "red"), lty = 1)
	dev.off()

	pdf(paste(n,"4_Bottleneck_Expansion_ABC_x.bottel.1.0.pdf",sep = "_"))
	plot(posterior[,10],posterior[,11], type = "l", col = "red",  main = paste("Model = ", n," 4_Bottleneck_Expansion " ,sep = "-"), xlab = "x.bottle", ylab="Posterior density", xlim = c(min(c(prior[,10], posterior[,10])), max(c(prior[,10], posterior[,10]))), ylim = c(min(c(prior[,11], posterior[,11])), max(c(prior[,11], posterior[,11]))))
	# lines(prior[,10], prior[,11], col = "blue")
	# legend(min(c(prior[,10], posterior[,10])),max(c(prior[,11], posterior[,11])),legend = c("prior","posterior"), col = c("blue", "red"), lty = 1)
	dev.off()

	pdf(paste(n,"4_Bottleneck_Expansion_ABC_exp.start.1.0.pdf",sep = "_"))
	plot(posterior[,12],posterior[,13], type = "l", col = "red",  main = paste("Model = ", n," 4_Bottleneck_Expansion " ,sep = "-"), xlab = "exp.start", ylab="Posterior density", xlim = c(min(c(prior[,12], posterior[,12])), max(c(prior[,12], posterior[,12]))), ylim = c(min(c(prior[,13], posterior[,13])), max(c(prior[,13], posterior[,13]))))
	# lines(prior[,12], prior[,13], col = "blue")
	# legend(min(c(prior[,12], posterior[,12])),max(c(prior[,13], posterior[,13])),legend = c("prior","posterior"), col = c("blue", "red"), lty = 1)
	dev.off()

	pdf(paste(n,"4_Bottleneck_Expansion_ABC_exp.interval.1.0.pdf",sep = "_"))
	plot(posterior[,14],posterior[,15], type = "l", col = "red",  main = paste("Model = ", n," 4_Bottleneck_Expansion " ,sep = "-"), xlab = "exp.interval", ylab="Posterior density", xlim = c(min(c(prior[,14], posterior[,14])), max(c(prior[,14], posterior[,14]))), ylim = c(min(c(prior[,15], posterior[,15])), max(c(prior[,15], posterior[,15]))))
	# lines(prior[,14], prior[,15], col = "blue")
	# legend(min(c(prior[,14], posterior[,14])),max(c(prior[,15], posterior[,15])),legend = c("prior","posterior"), col = c("blue", "red"), lty = 1)
	dev.off()

	setwd("../..")
}
