setwd("Data")

wolf.data <- read.table("Wolf_Info_For_Simulations.3.3.csv", head = T, sep = "\t")

EUR.pos <- which(wolf.data[,6] == "Europe" & wolf.data[,7] == 1)
CNE.pos <- which(wolf.data[,6] == "Central_North_Eurasia" & wolf.data[,7] == 1)
BER.pos <- which(wolf.data[,6] == "Beringia" & wolf.data[,7] == 1)
ME.pos <- which(wolf.data[,6] == "Middle_East" & wolf.data[,7] == 1)
SEA.pos <- which(wolf.data[,6] == "South_East_Asia"& wolf.data[,7] == 1)
NA.pos <- which(wolf.data[,6] == "North_America"& wolf.data[,7] == 1)
ANA.pos <- which(wolf.data[,6] == "Arctic_North_America"& wolf.data[,7] == 1)

DEMES5 <- list(as.numeric(EUR.pos),as.numeric(ME.pos),c(as.numeric(CNE.pos),as.numeric(SEA.pos),as.numeric(BER.pos)),c(as.numeric(NA.pos),as.numeric(ANA.pos)))

DEMES <- DEMES5
ndemes = length(DEMES)

############################################ Observed TMRC by Deme ####################################################

observed.tmrc.mat <- read.csv("TMRCA.mat.WOLVES.for_simulation.3.0.csv", head = F)
observed.tmrc.mat <- as.matrix(observed.tmrc.mat)

observed.tmrc.by.deme <- matrix(NA, nrow = ndemes, ncol = ndemes)
for(i in 1: ndemes){
	paste(i)
	for(j in 1: ndemes){
		paste(j)
		observed.tmrc.by.deme[i,j] <- mean(observed.tmrc.mat[as.array(unlist(DEMES[i])),as.array(unlist(DEMES[j]))])/1000
	}
}

abc.out.observed <- c(observed.tmrc.by.deme[1,1],observed.tmrc.by.deme[1,2], observed.tmrc.by.deme[1,3], observed.tmrc.by.deme[1,4], observed.tmrc.by.deme[2,3], observed.tmrc.by.deme[2,4], observed.tmrc.by.deme[3,3], observed.tmrc.by.deme[3,4], observed.tmrc.by.deme[4,4])

abc.out.observed.table <- matrix(abc.out.observed, ncol = 9, nrow = 1)
colnames(abc.out.observed.table)  <- c("tEUR-EUR", "tEUR-ME", "tEUR-NEA", "tEUR-AM", "tME-NEA","tME-AM", "tNEA-NEA", "tNEA-AM","tAM-AM")

setwd("../2_Expansion") # Provide a full path
write.table(abc.out.observed.table, file = "abc.out.observed.1.0.txt", sep = " ", row.names = F, quote = F)

################################################################################################
######################################## 2_EXPANSION ##############################################
################################################################################################

for(n in 1:7){
	print(n)
	name <- paste("a", n, "b.txt.gz", sep = "")

	output <- as.matrix(read.table(gzfile(name), sep = "\t", skip = 19))
    #    output <- as.matrix(read.table(name, sep = "\t", skip = 19))

##################################Ice-Sheet Filter Americas #####################################

	if(n == 6 | n == 7){
        print(paste("Number original", length(output[,1])))
		output = output[output[,7] >9,]
		output = output[output[,7] <16,]
        print(paste("Number retrained", length(output[,1])))
	}

################################################################################################

	sim.nr <- output[,1]
	r.sq <- output[,2]
	mig.rate <- output[,4]
	K <- output[,5]
	x.bottel <- output[,6]
	exp.start <- output[,7]
	exp.interval <- output[,8]
	nsims <- length(sim.nr)

	abc.out.bottleneck <- cbind(log10(mig.rate), log10(K), x.bottel, exp.start, exp.interval, output[,9:length(output[1,])])

	dir.create(paste("ABC",n,sep = ""))
	setwd(paste("ABC",n,sep = ""))
	write.table(abc.out.bottleneck, file = paste("abc.out_",n,"expansion.","1.0.txt", sep = "_"), sep = " ", col.names = c("log10(mig.rate)", "log10(K)", "x.bottel","exp.start","exp.interval","tEUR-EUR", "tEUR-ME", "tEUR-NEA","tME-NEA", "tNEA-NEA","tEUR-AM","tME-AM","tNEA-AM","tAM-AM"), row.names = F, quote = F)
	setwd("..")
}

setwd("..")
