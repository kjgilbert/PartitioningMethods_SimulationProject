library(scales)


setwd("~/Documents/UNIL/PartitioningMethods_SimulationProject/FunctioningPipeline/AnalysisResults")


make.merge.id <- function(dat){
	id.col <- rep(NA, dim(dat)[1])
	for(i in 1:dim(dat)[1]){
		if(length(grep("lopho", dat$file[i])) > 0) spp <- "lopho_"
		if(length(grep("myria", dat$file[i])) > 0) spp <- "myria_"
		
		if(length(grep("10", dat$file[i])) > 0) samp <- "10_"
		if(length(grep("20", dat$file[i])) > 0) samp <- "20_"
		if(length(grep("40", dat$file[i])) > 0) samp <- "40_"
		if(length(grep("80", dat$file[i])) > 0) samp <- "80_"
		
		if(length(grep("rcluster", dat$file[i])) > 0) algo <- "rcluster_"
		if(length(grep("greedy", dat$file[i])) > 0) algo <- "greedy_"
		if(length(grep("noPart", dat$file[i])) > 0) algo <- "noPart_"
	
		if(length(grep("AIC", dat$file[i])) > 0) crit <- "AIC"
		if(length(grep("AICc", dat$file[i])) > 0) crit <- "AICc"
		if(length(grep("BIC", dat$file[i])) > 0) crit <- "BIC"
		if(length(grep("noPart", dat$file[i])) > 0) crit <- "noPart"
		
		id.col[i] <- paste(c(spp, samp, algo, crit), collapse="")
	}
	return(id.col)
}


make.merge.id2 <- function(dat){
	id.col <- rep(NA, dim(dat)[1])
	for(i in 1:dim(dat)[1]){
		if(length(grep("lopho", dat$file[i])) > 0) spp <- "lopho_"
		if(length(grep("myria", dat$file[i])) > 0) spp <- "myria_"
		
		if(length(grep("10", dat$file[i])) > 0) samp <- "10_"
		if(length(grep("20", dat$file[i])) > 0) samp <- "20_"
		if(length(grep("40", dat$file[i])) > 0) samp <- "40_"
		if(length(grep("80", dat$file[i])) > 0) samp <- "80_"
		
		if(length(grep("r-", dat$file[i])) > 0) algo <- "rcluster_"
		if(length(grep("greedy", dat$file[i])) > 0) algo <- "greedy_"
		if(length(grep("noPart", dat$file[i])) > 0) algo <- "noPart_"
	
		if(length(grep("-a", dat$file[i])) > 0) crit <- "AIC"
		if(length(grep("-c", dat$file[i])) > 0) crit <- "AICc"
		if(length(grep("-b", dat$file[i])) > 0) crit <- "BIC"
		if(length(grep("noPart", dat$file[i])) > 0) crit <- "noPart"
		
		id.col[i] <- paste(c(spp, samp, algo, crit), collapse="")
	}
	return(id.col)
}




####################################

get.rep.difference <- function(dat, rep, samp.size, spp, partitioning.criterion, distance){
	temp.dat <- dat[dat$replicate == rep & dat$samp.size == samp.size & dat$species == spp , ]
	no.part.result <- temp.dat[temp.dat$search.algo == "noPart" ,]
	rclust.AIC.result <- temp.dat[temp.dat$search.algo == "rcluster" & temp.dat$part.IC == "AIC" ,]
	rclust.BIC.result <- temp.dat[temp.dat$search.algo == "rcluster" & temp.dat$part.IC == "BIC" ,]
	rclust.AICc.result <- temp.dat[temp.dat$search.algo == "rcluster" & temp.dat$part.IC == "AICc",]
	
	if(partitioning.criterion == "AIC"){
		if(dim(rclust.AIC.result[distance])[1] == 0 | dim(no.part.result[distance])[1] == 0){ difference <- NA}
		else{ difference <- no.part.result[distance] - rclust.AIC.result[distance] }
	}
	if(partitioning.criterion == "BIC"){
		if(dim(rclust.BIC.result[distance])[1] == 0 | dim(no.part.result[distance])[1] == 0){ difference <- NA}
		else{ difference <- no.part.result[distance] - rclust.BIC.result[distance] }
	}
	if(partitioning.criterion == "AICc"){
		if(dim(rclust.AICc.result[distance])[1] == 0 | dim(no.part.result[distance])[1] == 0){ difference <- NA}
		else{ difference <- no.part.result[distance] - rclust.AICc.result[distance] }
	}
	return(as.numeric(difference))			
}

####################################

get.boxplot.dat <- function(dat, samp.size, spp, partitioning.criterion, distance){	
	samp.size.dists <- NULL
	for(j in unique(dat$replicate)){
		dist <- get.rep.difference(dat, rep=j, samp.size=samp.size, spp=spp, partitioning.criterion=partitioning.criterion, distance=distance)
		samp.size.dists <- c(samp.size.dists, dist)
	}
	return(samp.size.dists)
}

####################################

get.partition.dat <- function(dat, samp.size, spp, partitioning.criterion){	
	partition.numbers <- NULL
	for(j in unique(dat$replicate)){
		
		temp.dat <- dat[dat$replicate == j & dat$samp.size == samp.size & dat$species == spp , ]
		rclust.AIC.result <- temp.dat[temp.dat$search.algo == "rcluster" & temp.dat$part.IC == "AIC" ,]
		rclust.BIC.result <- temp.dat[temp.dat$search.algo == "rcluster" & temp.dat$part.IC == "BIC" ,]
		rclust.AICc.result <- temp.dat[temp.dat$search.algo == "rcluster" & temp.dat$part.IC == "AICc",]
	
		if(partitioning.criterion == "AIC"){
			if(length(rclust.AIC.result$num.partitions)[1] == 0){ num.parts <- NA}
			else{ num.parts <- rclust.AIC.result$num.partitions }
		}
		if(partitioning.criterion == "BIC"){
			if(length(rclust.BIC.result$num.partitions)[1] == 0){ num.parts <- NA}
			else{ num.parts <- rclust.BIC.result$num.partitions }
		}
		if(partitioning.criterion == "AICc"){
			if(length(rclust.AICc.result$num.partitions)[1] == 0){ num.parts <- NA}
			else{ num.parts <- rclust.AICc.result$num.partitions }
		}
		partition.numbers <- c(partition.numbers, num.parts)
	}
	return(partition.numbers)
}

####################################



# mut rate 1000

infiles.a <- c("AllStatsOutput_New_Rep1_100-1000.csv",
				"AllStatsOutput_New_Rep2_100-1000.csv",
				"AllStatsOutput_New_Rep3_100-1000.csv",
				"AllStatsOutput_New_Rep4_100-1000.csv",
				"AllStatsOutput_New_Rep5_100-1000.csv",
				"AllStatsOutput_New_Rep6_100-1000.csv",
				"AllStatsOutput_Newest_Rep7_100-1000.csv",
				"AllStatsOutput_Newest_Rep8_100-1000.csv",
				"AllStatsOutput_Newest_Rep9_100-1000.csv",
				"AllStatsOutput_Newest_Rep10_100-1000.csv",
				"AllStatsOutput_Newest_Rep11_100-1000.csv",
				"AllStatsOutput_Newest_Rep12_100-1000.csv",
				"AllStatsOutput_Newest_Rep13_100-1000.csv",
				"AllStatsOutput_Newest_Rep14_100-1000.csv",
				"AllStatsOutput_Newest_Rep15_100-1000.csv",
				"AllStatsOutput_Newest_Rep16_100-1000.csv",
				"AllStatsOutput_Newest_Rep17_100-1000.csv",
				"AllStatsOutput_Newest_Rep18_100-1000.csv",
				"AllStatsOutput_Newest_Rep19_100-1000.csv",
				"AllStatsOutput_Newest_Rep20_100-1000.csv",
				"AllStatsOutput_Newest_Rep21_100-1000.csv",
				"AllStatsOutput_Newest_Rep22_100-1000.csv",
				"AllStatsOutput_Newest_Rep23_100-1000.csv",
				"AllStatsOutput_Newest_Rep24_100-1000.csv",
				"AllStatsOutput_Newest_Rep25_100-1000.csv",
				"AllStatsOutput_Newest_Rep26_100-1000.csv",
				"AllStatsOutput_Newest_Rep27_100-1000.csv",
				"AllStatsOutput_Newest_Rep28_100-1000.csv",
				"AllStatsOutput_Newest_Rep29_100-1000.csv",
				"AllStatsOutput_Newest_Rep30_100-1000.csv",
				"AllStatsOutput_Newest_Rep31_100-1000.csv",
				"AllStatsOutput_Newest_Rep32_100-1000.csv",
				"AllStatsOutput_Newest_Rep33_100-1000.csv",
				"AllStatsOutput_Newest_Rep34_100-1000.csv",
				"AllStatsOutput_Newest_Rep35_100-1000.csv"
				)
infiles.a2 <- c(	"AllStatsOutput_Rep1.csv",
				"AllStatsOutput_Rep2.csv",
				"AllStatsOutput_Rep3.csv",
				"AllStatsOutput_Rep4.csv",
				"AllStatsOutput_Rep5.csv",
				"AllStatsOutput_Rep6.csv",
				"AllStatsOutput_Rep7.csv",
				"AllStatsOutput_Rep8.csv",
				"AllStatsOutput_Rep9.csv",
				"AllStatsOutput_Rep10.csv"
				)
				
# IQtree results
infiles.b <- c("RF_Eucl_results_New_Rep1_100-1000_wConsensus.txt",
				"RF_Eucl_results_New_Rep2_100-1000_wConsensus.txt",
				"RF_Eucl_results_New_Rep3_100-1000_wConsensus.txt",
				"RF_Eucl_results_New_Rep4_100-1000_wConsensus.txt",
				"RF_Eucl_results_New_Rep5_100-1000_wConsensus.txt",
				"RF_Eucl_results_New_Rep6_100-1000_wConsensus.txt",
				"RF_Eucl_results_New_Rep7_100-1000_wConsensus.txt",
				"RF_Eucl_results_New_Rep8_100-1000_wConsensus.txt",
				"RF_Eucl_results_New_Rep9_100-1000_wConsensus.txt",
				"RF_Eucl_results_New_Rep10_100-1000_wConsensus.txt",
				"RF_Eucl_results_New_Rep11_100-1000_wConsensus.txt",
				"RF_Eucl_results_New_Rep12_100-1000_wConsensus.txt",
				"RF_Eucl_results_New_Rep13_100-1000_wConsensus.txt",
				"RF_Eucl_results_New_Rep14_100-1000_wConsensus.txt",
				"RF_Eucl_results_New_Rep15_100-1000_wConsensus.txt",
				"RF_Eucl_results_New_Rep16_100-1000_wConsensus.txt",
				"RF_Eucl_results_New_Rep17_100-1000_wConsensus.txt",
				"RF_Eucl_results_New_Rep18_100-1000_wConsensus.txt",
				"RF_Eucl_results_New_Rep19_100-1000_wConsensus.txt",
				"RF_Eucl_results_New_Rep20_100-1000_wConsensus.txt",
				"RF_Eucl_results_New_Rep21_100-1000_wConsensus.txt",
				"RF_Eucl_results_New_Rep22_100-1000_wConsensus.txt",
				"RF_Eucl_results_New_Rep23_100-1000_wConsensus.txt",
				"RF_Eucl_results_New_Rep24_100-1000_wConsensus.txt",
				"RF_Eucl_results_New_Rep25_100-1000_wConsensus.txt",
				"RF_Eucl_results_New_Rep26_100-1000_wConsensus.txt",
				"RF_Eucl_results_New_Rep27_100-1000_wConsensus.txt",
				"RF_Eucl_results_New_Rep28_100-1000_wConsensus.txt",
				"RF_Eucl_results_New_Rep29_100-1000_wConsensus.txt",
				"RF_Eucl_results_New_Rep30_100-1000_wConsensus.txt",
				"RF_Eucl_results_New_Rep31_100-1000_wConsensus.txt",
				"RF_Eucl_results_New_Rep32_100-1000_wConsensus.txt",
				"RF_Eucl_results_New_Rep33_100-1000_wConsensus.txt",
				"RF_Eucl_results_New_Rep34_100-1000_wConsensus.txt",
				"RF_Eucl_results_New_Rep35_100-1000_wConsensus.txt"
				)
infiles.b2 <- c(	"RF_Eucl_results_N10_Rep1_wConsensus.txt",
				"RF_Eucl_results_N10_Rep2_wConsensus.txt",
				"RF_Eucl_results_N10_Rep3_wConsensus.txt",
				"RF_Eucl_results_N10_Rep4_wConsensus.txt",
				"RF_Eucl_results_N10_Rep5_wConsensus.txt",
				"RF_Eucl_results_N10_Rep6_wConsensus.txt",
				"RF_Eucl_results_N10_Rep7_wConsensus.txt",
				"RF_Eucl_results_N10_Rep8_wConsensus.txt",
				"RF_Eucl_results_N10_Rep9_wConsensus.txt",
				"RF_Eucl_results_N10_Rep10_wConsensus.txt"
				)
				
## RAxML results
infiles.c <- c("RF_Eucl_results_New_Rep1_100-1000_RAxML.txt",
				"RF_Eucl_results_New_Rep2_100-1000_RAxML.txt",
				"RF_Eucl_results_New_Rep3_100-1000_RAxML.txt",
				"RF_Eucl_results_New_Rep4_100-1000_RAxML.txt",
				"RF_Eucl_results_New_Rep5_100-1000_RAxML.txt",
				"RF_Eucl_results_New_Rep6_100-1000_RAxML.txt",
				"RF_Eucl_results_New_Rep7_100-1000_RAxML.txt",
				"RF_Eucl_results_New_Rep8_100-1000_RAxML.txt",
				"RF_Eucl_results_New_Rep9_100-1000_RAxML.txt",
				"RF_Eucl_results_New_Rep10_100-1000_RAxML.txt",
				"RF_Eucl_results_New_Rep11_100-1000_RAxML.txt",
				"RF_Eucl_results_New_Rep12_100-1000_RAxML.txt",
				"RF_Eucl_results_New_Rep13_100-1000_RAxML.txt",
				"RF_Eucl_results_New_Rep14_100-1000_RAxML.txt",
				"RF_Eucl_results_New_Rep15_100-1000_RAxML.txt",
				"RF_Eucl_results_New_Rep16_100-1000_RAxML.txt",
				"RF_Eucl_results_New_Rep17_100-1000_RAxML.txt",
				"RF_Eucl_results_New_Rep18_100-1000_RAxML.txt",
				"RF_Eucl_results_New_Rep19_100-1000_RAxML.txt",
				"RF_Eucl_results_New_Rep20_100-1000_RAxML.txt",
				"RF_Eucl_results_New_Rep21_100-1000_RAxML.txt",
				"RF_Eucl_results_New_Rep22_100-1000_RAxML.txt",
				"RF_Eucl_results_New_Rep23_100-1000_RAxML.txt",
				"RF_Eucl_results_New_Rep24_100-1000_RAxML.txt",
				"RF_Eucl_results_New_Rep25_100-1000_RAxML.txt",
				"RF_Eucl_results_New_Rep26_100-1000_RAxML.txt",
				"RF_Eucl_results_New_Rep27_100-1000_RAxML.txt",
				"RF_Eucl_results_New_Rep28_100-1000_RAxML.txt",
				"RF_Eucl_results_New_Rep29_100-1000_RAxML.txt",
				"RF_Eucl_results_New_Rep30_100-1000_RAxML.txt",
				"RF_Eucl_results_New_Rep31_100-1000_RAxML.txt",
				"RF_Eucl_results_New_Rep32_100-1000_RAxML.txt",
				"RF_Eucl_results_New_Rep33_100-1000_RAxML.txt",
				"RF_Eucl_results_New_Rep34_100-1000_RAxML.txt",
				"RF_Eucl_results_New_Rep35_100-1000_RAxML.txt"
				)
infiles.c2 <- c("RF_Eucl_results_N10_Rep1_RAxML.txt",
				"RF_Eucl_results_N10_Rep2_RAxML.txt",
				"RF_Eucl_results_N10_Rep3_RAxML.txt",
				"RF_Eucl_results_N10_Rep4_RAxML.txt",
				"RF_Eucl_results_N10_Rep5_RAxML.txt",
				"RF_Eucl_results_N10_Rep6_RAxML.txt",
				"RF_Eucl_results_N10_Rep7_RAxML.txt",
				"RF_Eucl_results_N10_Rep8_RAxML.txt",
				"RF_Eucl_results_N10_Rep9_RAxML.txt",
				"RF_Eucl_results_N10_Rep10_RAxML.txt"
				)
				
dat <- NULL
for(i in c(1:length(infiles.a))){ #**************UNCOMMENT ME FOR ALL REPS
	temp.dat1 <- read.csv(infiles.a[i], header=TRUE, stringsAsFactors=FALSE)
	temp.dat1$uniqueID <- paste(temp.dat1$species, temp.dat1$samp.size, temp.dat1$search.algo, temp.dat1$part.IC, sep="_")
	
	temp.dat2 <- read.table(infiles.b[i], header=TRUE, stringsAsFactors=FALSE)
	# because also have consenseus trees now:
	temp.dat2 <- cbind(temp.dat2[1:32,], temp.dat2[33:64, 2:3])
	names(temp.dat2) <- c("filename", "rf_dist", "eucl_dist", "cons_rf_dist", "cons_eucl_dist")
	temp.dat2$uniqueID <- make.merge.id(temp.dat2)
	
	temp.dat3 <- read.table(infiles.c[i], header=TRUE, stringsAsFactors=FALSE)
	names(temp.dat3) <- c("filename_raxml", "rf_dist_raxml", "eucl_dist_raxml")
	temp.dat3$uniqueID <- make.merge.id2(temp.dat3)
	
	temp.temp.dat <- merge(temp.dat1, temp.dat2, by="uniqueID")
	temp.dat <- merge(temp.temp.dat, temp.dat3, by="uniqueID")
	temp.dat$gene.size <- 100
	temp.dat$mut.rate <- 1000
	temp.dat$replicate <- i

	dat <- rbind(dat, temp.dat)
}

# additional 10 reps at only 10 species and all else same, plus 10 species with rate variation and no indel model
dat2 <- NULL
for(i in c(1:length(infiles.a2))){ #**************UNCOMMENT ME FOR ALL REPS
	temp.dat1 <- read.csv(infiles.a2[i], header=TRUE, stringsAsFactors=FALSE)
	temp.dat1$uniqueID <- paste(temp.dat1$species, temp.dat1$samp.size, temp.dat1$search.algo, temp.dat1$part.IC, sep="_")
	
	temp.dat2 <- read.table(infiles.b2[i], header=TRUE, stringsAsFactors=FALSE)
	# because also have consenseus trees now:
	temp.dat2 <- cbind(temp.dat2[1:32,], temp.dat2[33:64, 2:3])
	names(temp.dat2) <- c("filename", "rf_dist", "eucl_dist", "cons_rf_dist", "cons_eucl_dist")
	temp.dat2$uniqueID <- make.merge.id(temp.dat2)
	
	temp.dat3 <- read.table(infiles.c2[i], header=TRUE, stringsAsFactors=FALSE)
	names(temp.dat3) <- c("filename_raxml", "rf_dist_raxml", "eucl_dist_raxml")
	temp.dat3$uniqueID <- make.merge.id2(temp.dat3)
	
	temp.temp.dat <- merge(temp.dat1, temp.dat2, by="uniqueID")
	temp.dat <- merge(temp.temp.dat, temp.dat3, by="uniqueID")
	temp.dat$gene.size <- 100
	temp.dat$mut.rate <- 1000
	temp.dat$replicate <- i
	
	dat2 <- rbind(dat2, temp.dat)
}
dat2$uniqueID <- gsub('lopho','TenSpp', dat2$uniqueID)
dat2$uniqueID <- gsub('myria','RateVar', dat2$uniqueID)
dat2$species <- gsub('lopho','TenSpp', dat2$species)
dat2$species <- gsub('myria','RateVar', dat2$species)

dat <- rbind(dat, dat2)

no.part.dat <- dat[dat$part.IC == "noPart" ,]


## IS there a statistically significant difference between the RF/eucl distance to true tree of no partitioning vs with partitioning??
# compare at each sample size, for each species, for both RF and euclidean distance

spp <- c("myria", "lopho", "TenSpp", "RateVar")
samps <- c(10,20,40,80)
stat.data.set.BIC <- data.frame(matrix(NA, ncol=20, nrow=8))
names(stat.data.set.BIC) <- c("species", "samp.size", "rf_twoside_p", "rf_less_p", "rf_greater_p", "eucl_twoside_p", "eucl_less_p", "eucl_greater_p", "cons_rf_twoside_p", "cons_rf_less_p", "cons_rf_greater_p", "cons_eucl_twoside_p", "cons_eucl_less_p", "cons_eucl_greater_p","raxml_rf_twoside_p", "raxml_rf_less_p", "raxml_rf_greater_p", "raxml_eucl_twoside_p", "raxml_eucl_less_p", "raxml_eucl_greater_p")

iterate <- 1
for(j in spp){
	for(k in samps){
		temp.dat.no.part <- dat[dat$samp.size == k & dat$species == j & dat$part.IC == "noPart" , ]
		temp.dat.part <- dat[dat$samp.size == k & dat$species == j & dat$part.IC == "BIC" , ]
		
		rf_ts_p <- wilcox.test(temp.dat.no.part$rf_dist, temp.dat.part$rf_dist, paired=TRUE, alternative="two.sided")$p.value # is there a significant diff?
		rf_less_p <- wilcox.test(temp.dat.no.part$rf_dist, temp.dat.part$rf_dist, paired=TRUE, alternative="less")$p.value # is the rf distance smaller in no part vs partitioned?
		rf_greater_p <- wilcox.test(temp.dat.no.part$rf_dist, temp.dat.part$rf_dist, paired=TRUE, alternative="greater")$p.value # is the rf distance smaller in partitioned trees?
		eucl_ts_p <- wilcox.test(temp.dat.no.part$eucl_dist, temp.dat.part$eucl_dist, paired=TRUE, alternative="two.sided")$p.value
		eucl_less_p <- wilcox.test(temp.dat.no.part$eucl_dist, temp.dat.part$eucl_dist, paired=TRUE, alternative="less")$p.value # is the euclidean distance smaller in no part vs partitioned?
		eucl_greater_p <- wilcox.test(temp.dat.no.part$eucl_dist, temp.dat.part$eucl_dist, paired=TRUE, alternative="greater")$p.value # is the euclidean distance smaller in partitioned trees?

		# also for the consensus tree
		cons_rf_ts_p <- wilcox.test(temp.dat.no.part$cons_rf_dist, temp.dat.part$cons_rf_dist, paired=TRUE, alternative="two.sided")$p.value # is there a significant diff?
		cons_rf_less_p <- wilcox.test(temp.dat.no.part$cons_rf_dist, temp.dat.part$cons_rf_dist, paired=TRUE, alternative="less")$p.value # is the rf distance smaller in no part vs partitioned?
		cons_rf_greater_p <- wilcox.test(temp.dat.no.part$cons_rf_dist, temp.dat.part$cons_rf_dist, paired=TRUE, alternative="greater")$p.value # is the rf distance smaller in partitioned trees?
		cons_eucl_ts_p <- wilcox.test(temp.dat.no.part$cons_eucl_dist, temp.dat.part$cons_eucl_dist, paired=TRUE, alternative="two.sided")$p.value
		cons_eucl_less_p <- wilcox.test(temp.dat.no.part$cons_eucl_dist, temp.dat.part$cons_eucl_dist, paired=TRUE, alternative="less")$p.value # is the euclidean distance smaller in no part vs partitioned?
		cons_eucl_greater_p <- wilcox.test(temp.dat.no.part$cons_eucl_dist, temp.dat.part$cons_eucl_dist, paired=TRUE, alternative="greater")$p.value # is the euclidean distance smaller in partitioned trees?
		
		# and for raxml
		rax_rf_ts_p <- wilcox.test(temp.dat.no.part$rf_dist_raxml, temp.dat.part$rf_dist_raxml, paired=TRUE, alternative="two.sided")$p.value
		rax_rf_less_p <- wilcox.test(temp.dat.no.part$rf_dist_raxml, temp.dat.part$rf_dist_raxml, paired=TRUE, alternative="less")$p.value
		rax_rf_greater_p <- wilcox.test(temp.dat.no.part$rf_dist_raxml, temp.dat.part$rf_dist_raxml, paired=TRUE, alternative="greater")$p.value
		rax_eucl_ts_p <- wilcox.test(temp.dat.no.part$eucl_dist_raxml, temp.dat.part$eucl_dist_raxml, paired=TRUE, alternative="two.sided")$p.value
		rax_eucl_less_p <- wilcox.test(temp.dat.no.part$eucl_dist_raxml, temp.dat.part$eucl_dist_raxml, paired=TRUE, alternative="less")$p.value
		rax_eucl_greater_p <- wilcox.test(temp.dat.no.part$eucl_dist_raxml, temp.dat.part$eucl_dist_raxml, paired=TRUE, alternative="greater")$p.value 


		stat.data.set.BIC[iterate, ] <- c(j, k, rf_ts_p, rf_less_p, rf_greater_p, eucl_ts_p, eucl_less_p, eucl_greater_p, cons_rf_ts_p, cons_rf_less_p, cons_rf_greater_p, cons_eucl_ts_p, cons_eucl_less_p, cons_eucl_greater_p, rax_rf_ts_p, rax_rf_less_p, rax_rf_greater_p, rax_eucl_ts_p, rax_eucl_less_p, rax_eucl_greater_p)
		iterate <- iterate + 1
	}
}





## AIC STATS

stat.data.set.AIC <- data.frame(matrix(NA, ncol=20, nrow=8))
names(stat.data.set.AIC) <- c("species", "samp.size", "rf_twoside_p", "rf_less_p", "rf_greater_p", "eucl_twoside_p", "eucl_less_p", "eucl_greater_p", "cons_rf_twoside_p", "cons_rf_less_p", "cons_rf_greater_p", "cons_eucl_twoside_p", "cons_eucl_less_p", "cons_eucl_greater_p","raxml_rf_twoside_p", "raxml_rf_less_p", "raxml_rf_greater_p", "raxml_eucl_twoside_p", "raxml_eucl_less_p", "raxml_eucl_greater_p")

iterate <- 1
for(j in spp){
	for(k in samps){
		temp.dat.no.part <- dat[dat$samp.size == k & dat$species == j & dat$part.IC == "noPart" , ]
		temp.dat.part <- dat[dat$samp.size == k & dat$species == j & dat$part.IC == "AIC" , ]
		
		rf_ts_p <- wilcox.test(temp.dat.no.part$rf_dist, temp.dat.part$rf_dist, paired=TRUE, alternative="two.sided")$p.value # is there a significant diff?
		rf_less_p <- wilcox.test(temp.dat.no.part$rf_dist, temp.dat.part$rf_dist, paired=TRUE, alternative="less")$p.value # is the rf distance smaller in no part vs partitioned?
		rf_greater_p <- wilcox.test(temp.dat.no.part$rf_dist, temp.dat.part$rf_dist, paired=TRUE, alternative="greater")$p.value # is the rf distance smaller in partitioned trees?
		eucl_ts_p <- wilcox.test(temp.dat.no.part$eucl_dist, temp.dat.part$eucl_dist, paired=TRUE, alternative="two.sided")$p.value
		eucl_less_p <- wilcox.test(temp.dat.no.part$eucl_dist, temp.dat.part$eucl_dist, paired=TRUE, alternative="less")$p.value # is the euclidean distance smaller in no part vs partitioned?
		eucl_greater_p <- wilcox.test(temp.dat.no.part$eucl_dist, temp.dat.part$eucl_dist, paired=TRUE, alternative="greater")$p.value # is the euclidean distance smaller in partitioned trees?

		# also for the consensus tree
		cons_rf_ts_p <- wilcox.test(temp.dat.no.part$cons_rf_dist, temp.dat.part$cons_rf_dist, paired=TRUE, alternative="two.sided")$p.value
		cons_rf_less_p <- wilcox.test(temp.dat.no.part$cons_rf_dist, temp.dat.part$cons_rf_dist, paired=TRUE, alternative="less")$p.value
		cons_rf_greater_p <- wilcox.test(temp.dat.no.part$cons_rf_dist, temp.dat.part$cons_rf_dist, paired=TRUE, alternative="greater")$p.value
		cons_eucl_ts_p <- wilcox.test(temp.dat.no.part$cons_eucl_dist, temp.dat.part$cons_eucl_dist, paired=TRUE, alternative="two.sided")$p.value
		cons_eucl_less_p <- wilcox.test(temp.dat.no.part$cons_eucl_dist, temp.dat.part$cons_eucl_dist, paired=TRUE, alternative="less")$p.value
		cons_eucl_greater_p <- wilcox.test(temp.dat.no.part$cons_eucl_dist, temp.dat.part$cons_eucl_dist, paired=TRUE, alternative="greater")$p.value
		
		# and for raxml
		rax_rf_ts_p <- wilcox.test(temp.dat.no.part$rf_dist_raxml, temp.dat.part$rf_dist_raxml, paired=TRUE, alternative="two.sided")$p.value
		rax_rf_less_p <- wilcox.test(temp.dat.no.part$rf_dist_raxml, temp.dat.part$rf_dist_raxml, paired=TRUE, alternative="less")$p.value
		rax_rf_greater_p <- wilcox.test(temp.dat.no.part$rf_dist_raxml, temp.dat.part$rf_dist_raxml, paired=TRUE, alternative="greater")$p.value
		rax_eucl_ts_p <- wilcox.test(temp.dat.no.part$eucl_dist_raxml, temp.dat.part$eucl_dist_raxml, paired=TRUE, alternative="two.sided")$p.value
		rax_eucl_less_p <- wilcox.test(temp.dat.no.part$eucl_dist_raxml, temp.dat.part$eucl_dist_raxml, paired=TRUE, alternative="less")$p.value
		rax_eucl_greater_p <- wilcox.test(temp.dat.no.part$eucl_dist_raxml, temp.dat.part$eucl_dist_raxml, paired=TRUE, alternative="greater")$p.value 

		stat.data.set.AIC[iterate, ] <- c(j, k, rf_ts_p, rf_less_p, rf_greater_p, eucl_ts_p, eucl_less_p, eucl_greater_p, cons_rf_ts_p, cons_rf_less_p, cons_rf_greater_p, cons_eucl_ts_p, cons_eucl_less_p, cons_eucl_greater_p, rax_rf_ts_p, rax_rf_less_p, rax_rf_greater_p, rax_eucl_ts_p, rax_eucl_less_p, rax_eucl_greater_p)
		iterate <- iterate + 1
	}
}




### AICc STATS

stat.data.set.AICc <- data.frame(matrix(NA, ncol=20, nrow=8))
names(stat.data.set.AICc) <- c("species", "samp.size", "rf_twoside_p", "rf_less_p", "rf_greater_p", "eucl_twoside_p", "eucl_less_p", "eucl_greater_p", "cons_rf_twoside_p", "cons_rf_less_p", "cons_rf_greater_p", "cons_eucl_twoside_p", "cons_eucl_less_p", "cons_eucl_greater_p", ,"raxml_rf_twoside_p", "raxml_rf_less_p", "raxml_rf_greater_p", "raxml_eucl_twoside_p", "raxml_eucl_less_p", "raxml_eucl_greater_p")

iterate <- 1
for(j in spp){
	for(k in samps){
		temp.dat.no.part <- dat[dat$samp.size == k & dat$species == j & dat$part.IC == "noPart" , ]
		temp.dat.part <- dat[dat$samp.size == k & dat$species == j & dat$part.IC == "AICc" , ]
		
		rf_ts_p <- wilcox.test(temp.dat.no.part$rf_dist, temp.dat.part$rf_dist, paired=TRUE, alternative="two.sided")$p.value # is there a significant diff?
		rf_less_p <- wilcox.test(temp.dat.no.part$rf_dist, temp.dat.part$rf_dist, paired=TRUE, alternative="less")$p.value # is the rf distance smaller in no part vs partitioned?
		rf_greater_p <- wilcox.test(temp.dat.no.part$rf_dist, temp.dat.part$rf_dist, paired=TRUE, alternative="greater")$p.value # is the rf distance smaller in partitioned trees?
		eucl_ts_p <- wilcox.test(temp.dat.no.part$eucl_dist, temp.dat.part$eucl_dist, paired=TRUE, alternative="two.sided")$p.value
		eucl_less_p <- wilcox.test(temp.dat.no.part$eucl_dist, temp.dat.part$eucl_dist, paired=TRUE, alternative="less")$p.value # is the euclidean distance smaller in no part vs partitioned?
		eucl_greater_p <- wilcox.test(temp.dat.no.part$eucl_dist, temp.dat.part$eucl_dist, paired=TRUE, alternative="greater")$p.value # is the euclidean distance smaller in partitioned trees?

		# also for the consensus tree
		cons_rf_ts_p <- wilcox.test(temp.dat.no.part$cons_rf_dist, temp.dat.part$cons_rf_dist, paired=TRUE, alternative="two.sided")$p.value
		cons_rf_less_p <- wilcox.test(temp.dat.no.part$cons_rf_dist, temp.dat.part$cons_rf_dist, paired=TRUE, alternative="less")$p.value
		cons_rf_greater_p <- wilcox.test(temp.dat.no.part$cons_rf_dist, temp.dat.part$cons_rf_dist, paired=TRUE, alternative="greater")$p.value
		cons_eucl_ts_p <- wilcox.test(temp.dat.no.part$cons_eucl_dist, temp.dat.part$cons_eucl_dist, paired=TRUE, alternative="two.sided")$p.value
		cons_eucl_less_p <- wilcox.test(temp.dat.no.part$cons_eucl_dist, temp.dat.part$cons_eucl_dist, paired=TRUE, alternative="less")$p.value
		cons_eucl_greater_p <- wilcox.test(temp.dat.no.part$cons_eucl_dist, temp.dat.part$cons_eucl_dist, paired=TRUE, alternative="greater")$p.value
		
		# and for raxml
		rax_rf_ts_p <- wilcox.test(temp.dat.no.part$rf_dist_raxml, temp.dat.part$rf_dist_raxml, paired=TRUE, alternative="two.sided")$p.value
		rax_rf_less_p <- wilcox.test(temp.dat.no.part$rf_dist_raxml, temp.dat.part$rf_dist_raxml, paired=TRUE, alternative="less")$p.value
		rax_rf_greater_p <- wilcox.test(temp.dat.no.part$rf_dist_raxml, temp.dat.part$rf_dist_raxml, paired=TRUE, alternative="greater")$p.value
		rax_eucl_ts_p <- wilcox.test(temp.dat.no.part$eucl_dist_raxml, temp.dat.part$eucl_dist_raxml, paired=TRUE, alternative="two.sided")$p.value
		rax_eucl_less_p <- wilcox.test(temp.dat.no.part$eucl_dist_raxml, temp.dat.part$eucl_dist_raxml, paired=TRUE, alternative="less")$p.value
		rax_eucl_greater_p <- wilcox.test(temp.dat.no.part$eucl_dist_raxml, temp.dat.part$eucl_dist_raxml, paired=TRUE, alternative="greater")$p.value 

		stat.data.set.AICc[iterate, ] <- c(j, k, rf_ts_p, rf_less_p, rf_greater_p, eucl_ts_p, eucl_less_p, eucl_greater_p, cons_rf_ts_p, cons_rf_less_p, cons_rf_greater_p, cons_eucl_ts_p, cons_eucl_less_p, cons_eucl_greater_p, rax_rf_ts_p, rax_rf_less_p, rax_rf_greater_p, rax_eucl_ts_p, rax_eucl_less_p, rax_eucl_greater_p)
		iterate <- iterate + 1
	}
}


