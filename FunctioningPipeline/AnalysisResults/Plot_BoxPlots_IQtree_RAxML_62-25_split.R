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




get.rep.difference <- function(dat, rep, samp.size, spp, partitioning.criterion, distance){
	temp.dat <- dat[dat$replicate == rep & dat$samp.size == samp.size & dat$species == spp , ]
	no.part.result <- temp.dat[temp.dat$search.algo == "noPart" ,]
	rclust.AIC.result <- temp.dat[temp.dat$search.algo == "rcluster" & temp.dat$part.IC == "AIC",]
	rclust.BIC.result <- temp.dat[temp.dat$search.algo == "rcluster" & temp.dat$part.IC == "BIC",]
	rclust.AICc.result <- temp.dat[temp.dat$search.algo == "rcluster" & temp.dat$part.IC == "AICc",]
	
	if(partitioning.criterion == "AIC"){
		num.partitions.inferred <- rclust.AIC.result$num.partitions
		if(dim(rclust.AIC.result[distance])[1] == 0 | dim(no.part.result[distance])[1] == 0){ difference <- NA}
		else{ difference <- no.part.result[distance] - rclust.AIC.result[distance] }
		part.dist <- rclust.AIC.result[distance]
	}
	if(partitioning.criterion == "BIC"){
		num.partitions.inferred <- rclust.BIC.result$num.partitions
		if(dim(rclust.BIC.result[distance])[1] == 0 | dim(no.part.result[distance])[1] == 0){ difference <- NA}
		else{ difference <- no.part.result[distance] - rclust.BIC.result[distance] }
		part.dist <- rclust.BIC.result[distance]
	}
	if(partitioning.criterion == "AICc"){
		num.partitions.inferred <- rclust.AICc.result$num.partitions
		if(dim(rclust.AICc.result[distance])[1] == 0 | dim(no.part.result[distance])[1] == 0){ difference <- NA}
		else{ difference <- no.part.result[distance] - rclust.AICc.result[distance] }
		part.dist <- rclust.AICc.result[distance]
	}
	return(c(no.part.result[distance], part.dist, as.numeric(difference), num.partitions.inferred))		
}

####################################

get.boxplot.dat <- function(dat, samp.size, spp, partitioning.criterion, distance){	
	samp.size.dists <- data.frame(matrix(NA, ncol=4))
	names(samp.size.dists) <- c("no.part.raw", "part.raw", "distance", "num.parts")
	for(j in unique(dat$replicate)){
		result <- get.rep.difference(dat, rep=j, samp.size=samp.size, spp=spp, partitioning.criterion=partitioning.criterion, distance=distance)
		#dist <- result[3]
		#parts <- result[4]
		samp.size.dists[j,] <- c(result[1], result[2], result[3], result[4])
	}
	return(samp.size.dists)
}




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















####################################
samp.size.col <- c("steelblue4", "darkorange2", "green4", "red3")
outline.col <- "gray30"
samp.sizes <- c(10,20,40,80)
####################################
ylim1 <- c(-0.2775, 0.3)
#ylim1 <- c(-0.1, 0.1)
#ylim2 <- c(-1.35, 1.35)
ylim2 <- c(-1.35, 5)
ylim3 <- c(0.5,8)
box.width <- 0.1
line.type <- 1
line.width <- 1.75
ytext1 <- 0.25
ytext2 <- 1.4
p.height1 <- 0.3
#p.height1 <- 0.1
#p.height2 <- 1.375
p.height2 <- 5.125
text.size <- 0.825


#criteria <- c("AIC", "BIC", "AICc")
#for(criterion in criteria){
criterion <- "BIC"


if(criterion == "AIC") stat.data.set <- stat.data.set.AIC
if(criterion == "BIC") stat.data.set <- stat.data.set.BIC
if(criterion == "AICc") stat.data.set <- stat.data.set.AICc
part.crit <- criterion

pdf(paste(c("Results_Plots/Boxplots_", part.crit, "_IQtree-RAxML_split.pdf"), collapse=""), width=8, height=10)
par(mar=c(2.5,4.25,1.75,0.25))
layout(matrix(c(1,1,3,3,1,1,3,3,1,1,3,3,2,2,4,4,2,2,4,4,2,2,4,4,5,5,5,5), ncol=4, byrow=TRUE))


#_________________________________________________________________________________________________________#

plot(0, type="n", xlim=c(0.5,9), ylim=ylim1, xlab="", xaxt="n", ylab="difference in distance from true tree", main="Robinson Fould distance - IQtree", mgp=c(3,0.9,0))
axis(side=1, at=c(2.5, 7), c("Lophotrochozoa (62 spp)", "Myriapoda (25 spp)"))
mtext(side=2, "(no partitioning - partitioning)", cex=0.75, line=2)
abline(h=0, col="gray75", lty=3, lwd=1.5)
text(9,ytext1, expression(italic("partitioning better")), cex=1, col="green4")
text(9,-ytext1, expression(italic("no partitioning better")), cex=1, col="red3")

# lopho (62 spp)
for(i in 1:4){
	x.spot1 <- i-0.2
	x.spot2 <- i+0.2
	plot.dat <- get.boxplot.dat(dat=dat, samp.size=samp.sizes[i], spp="lopho", partitioning.criterion=part.crit, distance="rf_dist")
	one.partition.dat <- plot.dat[plot.dat$num.parts == 1 ,]
	more.partitions.dat <- plot.dat[plot.dat$num.parts > 1 ,]
	if(dim(one.partition.dat)[1] > 0) boxplot(one.partition.dat$distance, add=TRUE, at=x.spot1, pch=23, col=samp.size.col[i], border=outline.col, width=box.width, lty=line.type, lwd=line.width, axes=FALSE)
	if(dim(more.partitions.dat)[1] > 0) boxplot(more.partitions.dat$distance, add=TRUE, at=x.spot2, pch=23, col=samp.size.col[i], border=outline.col, width=box.width, lty=line.type, lwd=line.width, axes=FALSE)
	
	if(dim(one.partition.dat)[1] > 0) text(x=x.spot1, y=p.height1, format(as.numeric(wilcox.test(one.partition.dat$no.part.raw, one.partition.dat$part.raw, paired=TRUE, alternative="two.sided")$p.value), digits=2), cex=text.size, srt=45)
	if(dim(more.partitions.dat)[1] > 0) text(x=x.spot2, y=p.height1, format(as.numeric(wilcox.test(more.partitions.dat$no.part.raw, more.partitions.dat$part.raw, paired=TRUE, alternative="two.sided")$p.value), digits=2), cex=text.size, srt=45)
}

# myria (25 spp)
for(i in 1:4){
	x.spot1 <- i + 4.5 - 0.2
	x.spot2 <- i + 4.5 + 0.2
	plot.dat <- get.boxplot.dat(dat=dat, samp.size=samp.sizes[i], spp="myria", partitioning.criterion=part.crit, distance="rf_dist")
	one.partition.dat <- plot.dat[plot.dat$num.parts == 1 ,]
	more.partitions.dat <- plot.dat[plot.dat$num.parts > 1 ,]
	if(dim(one.partition.dat)[1] > 0) boxplot(one.partition.dat$distance, add=TRUE, at=x.spot1, pch=23, col=samp.size.col[i], border=outline.col, width=box.width, lty=line.type, lwd=line.width, axes=FALSE)
	if(dim(more.partitions.dat)[1] > 0) boxplot(more.partitions.dat$distance, add=TRUE, at=x.spot2, pch=23, col=samp.size.col[i], border=outline.col, width=box.width, lty=line.type, lwd=line.width, axes=FALSE)
	
	if(dim(one.partition.dat)[1] > 0) text(x=x.spot1, y=p.height1, format(as.numeric(wilcox.test(one.partition.dat$no.part.raw, one.partition.dat$part.raw, paired=TRUE, alternative="two.sided")$p.value), digits=2), cex=text.size, srt=45)
	if(dim(more.partitions.dat)[1] > 0) text(x=x.spot2, y=p.height1, format(as.numeric(wilcox.test(more.partitions.dat$no.part.raw, more.partitions.dat$part.raw, paired=TRUE, alternative="two.sided")$p.value), digits=2), cex=text.size, srt=45)
}

legend("bottomright", c("10", "20", "40", "80"), pch=15, col=samp.size.col, title = "num. genes", cex=0.875, pt.cex=1.5) # , bg="white"


#_________________________________________________________________________________________________________#


# eucl distance
plot(0, type="n", xlim=c(0.5,9), ylim=ylim2, xlab="", xaxt="n", ylab="difference in distance from true tree", main="Euclidean distance - IQtree", mgp=c(3,0.9,0))
axis(side=1, at=c(2.5, 7), c("Lophotrochozoa (62 spp)", "Myriapoda (25 spp)"))
mtext(side=2, "(no partitioning - partitioning)", cex=0.75, line=2)
abline(h=0, col="gray75", lty=3, lwd=1.5)
text(9,ytext1, expression(italic("partitioning better")), cex=1, col="green4")
text(9,-ytext1, expression(italic("no partitioning better")), cex=1, col="red3")

# lopho (62 spp)
for(i in 1:4){
	x.spot1 <- i-0.2
	x.spot2 <- i+0.2
	plot.dat <- get.boxplot.dat(dat=dat, samp.size=samp.sizes[i], spp="lopho", partitioning.criterion=part.crit, distance="eucl_dist")
	one.partition.dat <- plot.dat[plot.dat$num.parts == 1 ,]
	more.partitions.dat <- plot.dat[plot.dat$num.parts > 1 ,]
	if(dim(one.partition.dat)[1] > 0) boxplot(one.partition.dat$distance, add=TRUE, at=x.spot1, pch=23, col=samp.size.col[i], border=outline.col, width=box.width, lty=line.type, lwd=line.width, axes=FALSE)
	if(dim(more.partitions.dat)[1] > 0) boxplot(more.partitions.dat$distance, add=TRUE, at=x.spot2, pch=23, col=samp.size.col[i], border=outline.col, width=box.width, lty=line.type, lwd=line.width, axes=FALSE)
	
	if(dim(one.partition.dat)[1] > 0) text(x=x.spot1, y=p.height2, format(as.numeric(wilcox.test(one.partition.dat$no.part.raw, one.partition.dat$part.raw, paired=TRUE, alternative="two.sided")$p.value), digits=2), cex=text.size, srt=45)
	if(dim(more.partitions.dat)[1] > 0) text(x=x.spot2, y=p.height2, format(as.numeric(wilcox.test(more.partitions.dat$no.part.raw, more.partitions.dat$part.raw, paired=TRUE, alternative="two.sided")$p.value), digits=2), cex=text.size, srt=45)
}

# myria (25 spp)
for(i in 1:4){
	x.spot1 <- i + 4.5 - 0.2
	x.spot2 <- i + 4.5 + 0.2
	plot.dat <- get.boxplot.dat(dat=dat, samp.size=samp.sizes[i], spp="myria", partitioning.criterion=part.crit, distance="eucl_dist")
	one.partition.dat <- plot.dat[plot.dat$num.parts == 1 ,]
	more.partitions.dat <- plot.dat[plot.dat$num.parts > 1 ,]
	if(dim(one.partition.dat)[1] > 0) boxplot(one.partition.dat$distance, add=TRUE, at=x.spot1, pch=23, col=samp.size.col[i], border=outline.col, width=box.width, lty=line.type, lwd=line.width, axes=FALSE)
	if(dim(more.partitions.dat)[1] > 0) boxplot(more.partitions.dat$distance, add=TRUE, at=x.spot2, pch=23, col=samp.size.col[i], border=outline.col, width=box.width, lty=line.type, lwd=line.width, axes=FALSE)
	
	if(dim(one.partition.dat)[1] > 0) text(x=x.spot1, y=p.height2, format(as.numeric(wilcox.test(one.partition.dat$no.part.raw, one.partition.dat$part.raw, paired=TRUE, alternative="two.sided")$p.value), digits=2), cex=text.size, srt=45)
	if(dim(more.partitions.dat)[1] > 0) text(x=x.spot2, y=p.height2, format(as.numeric(wilcox.test(more.partitions.dat$no.part.raw, more.partitions.dat$part.raw, paired=TRUE, alternative="two.sided")$p.value), digits=2), cex=text.size, srt=45)
}

legend("bottomright", c("10", "20", "40", "80"), pch=15, col=samp.size.col, title = "num. genes", cex=0.875, pt.cex=1.5) # , bg="white"

dev.off()
#_________________________________________________________________________________________________________#



# RAxML RF distance
plot(0, type="n", xlim=c(0.5,18), ylim=ylim1, xlab="", xaxt="n", ylab="difference in distance from true tree", main="Robinson Fould distance - RAxML", mgp=c(3,0.9,0))
mtext(side=2, "(no partitioning - partitioning)", cex=0.75, line=2)
axis(side=1, at=c(2.5, 7, 11.5, 16), c("Lophotrochozoa (62 spp)", "Myriapoda (25 spp)", "10 spp", "10 spp w/ rate variation"))
abline(h=0, col="gray75", lty=3, lwd=1.5)
text(9,ytext1, expression(italic("partitioning better")), cex=1, col="green4")
text(9,-ytext1, expression(italic("no partitioning better")), cex=1, col="red3")
for(i in 1:4){
	plot.dat <- get.boxplot.dat(dat=dat, samp.size=samp.sizes[i], spp="lopho", partitioning.criterion=part.crit, distance="rf_dist_raxml")
	boxplot(plot.dat, add=TRUE, at=i, pch=23, col=samp.size.col[i], border=outline.col, width=box.width, lty=line.type, lwd=line.width, axes=FALSE)
}
text(x=c(0.35,1:4), y=rep(p.height1, 5), cex=text.size, c(expression(italic("p")*" = "), 
	format(as.numeric(stat.data.set$raxml_rf_twoside_p[stat.data.set$species == "lopho" & stat.data.set$samp.size == 10][1]), digits=2), 
	format(as.numeric(stat.data.set$raxml_rf_twoside_p[stat.data.set$species == "lopho" & stat.data.set$samp.size == 20][1]), digits=2),
	format(as.numeric(stat.data.set$raxml_rf_twoside_p[stat.data.set$species == "lopho" & stat.data.set$samp.size == 40][1]), digits=2),
	format(as.numeric(stat.data.set$raxml_rf_twoside_p[stat.data.set$species == "lopho" & stat.data.set$samp.size == 80][1]), digits=2)))

for(i in 1:4){
	plot.dat <- get.boxplot.dat(dat=dat, samp.size=samp.sizes[i], spp="myria", partitioning.criterion=part.crit, distance="rf_dist_raxml")
	boxplot(plot.dat, add=TRUE, at=i+4.5, pch=23, col=samp.size.col[i], border=outline.col, width=box.width, lty=line.type, lwd=line.width, axes=FALSE)
}
text(x=1:4+4.5, y=rep(p.height1, 4), cex=text.size, c(
	format(as.numeric(stat.data.set$raxml_rf_twoside_p[stat.data.set$species == "myria" & stat.data.set$samp.size == 10][1]), digits=2), 
	format(as.numeric(stat.data.set$raxml_rf_twoside_p[stat.data.set$species == "myria" & stat.data.set$samp.size == 20][1]), digits=2),
	format(as.numeric(stat.data.set$raxml_rf_twoside_p[stat.data.set$species == "myria" & stat.data.set$samp.size == 40][1]), digits=2),
	format(as.numeric(stat.data.set$raxml_rf_twoside_p[stat.data.set$species == "myria" & stat.data.set$samp.size == 80][1]), digits=2)))

# 10 spp
for(i in 1:4){
	plot.dat <- get.boxplot.dat(dat=dat, samp.size=samp.sizes[i], spp="TenSpp", partitioning.criterion=part.crit, distance="rf_dist_raxml")
	boxplot(plot.dat, add=TRUE, at=i+9, pch=23, col=samp.size.col[i], border=outline.col, width=box.width, lty=line.type, lwd=line.width, axes=FALSE)
}
text(x=1:4+9, y=rep(p.height1, 4), cex=text.size, c(
	format(as.numeric(stat.data.set$raxml_rf_twoside_p[stat.data.set$species == "TenSpp" & stat.data.set$samp.size == 10][1]), digits=2), 
	format(as.numeric(stat.data.set$raxml_rf_twoside_p[stat.data.set$species == "TenSpp" & stat.data.set$samp.size == 20][1]), digits=2),
	format(as.numeric(stat.data.set$raxml_rf_twoside_p[stat.data.set$species == "TenSpp" & stat.data.set$samp.size == 40][1]), digits=2),
	format(as.numeric(stat.data.set$raxml_rf_twoside_p[stat.data.set$species == "TenSpp" & stat.data.set$samp.size == 80][1]), digits=2)))

# rate variation
for(i in 1:4){
	plot.dat <- get.boxplot.dat(dat=dat, samp.size=samp.sizes[i], spp="RateVar", partitioning.criterion=part.crit, distance="rf_dist_raxml")
	boxplot(plot.dat, add=TRUE, at=i++13.5, pch=23, col=samp.size.col[i], border=outline.col, width=box.width, lty=line.type, lwd=line.width, axes=FALSE)
}
text(x=1:4+13.5, y=rep(p.height1, 4), cex=text.size, c(
	format(as.numeric(stat.data.set$raxml_rf_twoside_p[stat.data.set$species == "RateVar" & stat.data.set$samp.size == 10][1]), digits=2), 
	format(as.numeric(stat.data.set$raxml_rf_twoside_p[stat.data.set$species == "RateVar" & stat.data.set$samp.size == 20][1]), digits=2),
	format(as.numeric(stat.data.set$raxml_rf_twoside_p[stat.data.set$species == "RateVar" & stat.data.set$samp.size == 40][1]), digits=2),
	format(as.numeric(stat.data.set$raxml_rf_twoside_p[stat.data.set$species == "RateVar" & stat.data.set$samp.size == 80][1]), digits=2)))


#_________________________________________________________________________________________________________#



# RAxML eucl distance
plot(0, type="n", xlim=c(0.5,18), ylim=ylim2, xlab="", xaxt="n", ylab="difference in distance from true tree", main="Euclidean distance - RAxML", mgp=c(3,0.9,0))
mtext(side=2, "(no partitioning - partitioning)", cex=0.75, line=2)
axis(side=1, at=c(2.5, 7, 11.5, 16), c("Lophotrochozoa (62 spp)", "Myriapoda (25 spp)", "10 spp", "10 spp w/ rate variation"))
abline(h=0, col="gray75", lty=3, lwd=1.5)
text(9,ytext2, expression(italic("partitioning better")), cex=1, col="green4")
text(9,-ytext2, expression(italic("no partitioning better")), cex=1, col="red3")
for(i in 1:4){
	plot.dat <- get.boxplot.dat(dat=dat, samp.size=samp.sizes[i], spp="lopho", partitioning.criterion=part.crit, distance="eucl_dist_raxml")
	boxplot(plot.dat, add=TRUE, at=i, pch=23, col=samp.size.col[i], border=outline.col, width=box.width, lty=line.type, lwd=line.width, axes=FALSE)
}
text(x=c(0.35,1:4), y=rep(p.height2, 5), cex=text.size, c(expression(italic("p")*" = "), 
	format(as.numeric(stat.data.set$raxml_eucl_twoside_p[stat.data.set$species == "lopho" & stat.data.set$samp.size == 10][1]), digits=2), 
	format(as.numeric(stat.data.set$raxml_eucl_twoside_p[stat.data.set$species == "lopho" & stat.data.set$samp.size == 20][1]), digits=2),
	format(as.numeric(stat.data.set$raxml_eucl_twoside_p[stat.data.set$species == "lopho" & stat.data.set$samp.size == 40][1]), digits=2),
	format(as.numeric(stat.data.set$raxml_eucl_twoside_p[stat.data.set$species == "lopho" & stat.data.set$samp.size == 80][1]), digits=2)))

for(i in 1:4){
	plot.dat <- get.boxplot.dat(dat=dat, samp.size=samp.sizes[i], spp="myria", partitioning.criterion=part.crit, distance="eucl_dist_raxml")
	boxplot(plot.dat, add=TRUE, at=i+4.5, pch=23, col=samp.size.col[i], border=outline.col, width=box.width, lty=line.type, lwd=line.width, axes=FALSE)
}
text(x=1:4+4.5, y=rep(p.height2, 4), cex=text.size, c(
	format(as.numeric(stat.data.set$raxml_eucl_twoside_p[stat.data.set$species == "myria" & stat.data.set$samp.size == 10][1]), digits=2), 
	format(as.numeric(stat.data.set$raxml_eucl_twoside_p[stat.data.set$species == "myria" & stat.data.set$samp.size == 20][1]), digits=2),
	format(as.numeric(stat.data.set$raxml_eucl_twoside_p[stat.data.set$species == "myria" & stat.data.set$samp.size == 40][1]), digits=2),
	format(as.numeric(stat.data.set$raxml_eucl_twoside_p[stat.data.set$species == "myria" & stat.data.set$samp.size == 80][1]), digits=2)))

# 10 spp
for(i in 1:4){
	plot.dat <- get.boxplot.dat(dat=dat, samp.size=samp.sizes[i], spp="TenSpp", partitioning.criterion=part.crit, distance="eucl_dist_raxml")
	boxplot(plot.dat, add=TRUE, at=i+9, pch=23, col=samp.size.col[i], border=outline.col, width=box.width, lty=line.type, lwd=line.width, axes=FALSE)
}
text(x=1:4+9, y=rep(p.height2, 4), cex=text.size, c(
	format(as.numeric(stat.data.set$raxml_eucl_twoside_p[stat.data.set$species == "TenSpp" & stat.data.set$samp.size == 10][1]), digits=2), 
	format(as.numeric(stat.data.set$raxml_eucl_twoside_p[stat.data.set$species == "TenSpp" & stat.data.set$samp.size == 20][1]), digits=2),
	format(as.numeric(stat.data.set$raxml_eucl_twoside_p[stat.data.set$species == "TenSpp" & stat.data.set$samp.size == 40][1]), digits=2),
	format(as.numeric(stat.data.set$raxml_eucl_twoside_p[stat.data.set$species == "TenSpp" & stat.data.set$samp.size == 80][1]), digits=2)))

# rate variation
for(i in 1:4){
	plot.dat <- get.boxplot.dat(dat=dat, samp.size=samp.sizes[i], spp="RateVar", partitioning.criterion=part.crit, distance="eucl_dist_raxml")
	boxplot(plot.dat, add=TRUE, at=i++13.5, pch=23, col=samp.size.col[i], border=outline.col, width=box.width, lty=line.type, lwd=line.width, axes=FALSE)
}
text(x=1:4+13.5, y=rep(p.height2, 4), cex=text.size, c(
	format(as.numeric(stat.data.set$raxml_eucl_twoside_p[stat.data.set$species == "RateVar" & stat.data.set$samp.size == 10][1]), digits=2), 
	format(as.numeric(stat.data.set$raxml_eucl_twoside_p[stat.data.set$species == "RateVar" & stat.data.set$samp.size == 20][1]), digits=2),
	format(as.numeric(stat.data.set$raxml_eucl_twoside_p[stat.data.set$species == "RateVar" & stat.data.set$samp.size == 40][1]), digits=2),
	format(as.numeric(stat.data.set$raxml_eucl_twoside_p[stat.data.set$species == "RateVar" & stat.data.set$samp.size == 80][1]), digits=2)))





#_________________________________________________________________________________________________________#


# number partitions
par(mar=c(2.5,4.25,0.25,0.25))
plot(0, type="n", xlim=c(0.5,18), ylim=ylim3, xlab="", xaxt="n", ylab="number partitions inferred", main="", mgp=c(3,0.9,0))
axis(side=1, at=c(2.5, 7, 11.5, 16), c("Lophotrochozoa (62 spp)", "Myriapoda (25 spp)", "10 spp", "10 spp w/ rate variation"))
abline(h=1, col="gray50", lty=2)
for(i in 1:4){
	part.dat <- get.partition.dat(dat=dat, samp.size=samp.sizes[i], spp="lopho", partitioning.criterion=part.crit)
	boxplot(part.dat, add=TRUE, at=i, pch=23, col=samp.size.col[i], border=outline.col, width=box.width, lty=line.type, lwd=line.width, axes=FALSE)
}
for(i in 1:4){
	part.dat <- get.partition.dat(dat=dat, samp.size=samp.sizes[i], spp="myria", partitioning.criterion=part.crit)
	boxplot(part.dat, add=TRUE, at=i+4.5, pch=23, col=samp.size.col[i], border=outline.col, width=box.width, lty=line.type, lwd=line.width, axes=FALSE)
}
for(i in 1:4){
	part.dat <- get.partition.dat(dat=dat, samp.size=samp.sizes[i], spp="TenSpp", partitioning.criterion=part.crit)
	boxplot(part.dat, add=TRUE, at=i+9, pch=23, col=samp.size.col[i], border=outline.col, width=box.width, lty=line.type, lwd=line.width, axes=FALSE)
}
for(i in 1:4){
	part.dat <- get.partition.dat(dat=dat, samp.size=samp.sizes[i], spp="RateVar", partitioning.criterion=part.crit)
	boxplot(part.dat, add=TRUE, at=i+13.5, pch=23, col=samp.size.col[i], border=outline.col, width=box.width, lty=line.type, lwd=line.width, axes=FALSE)
}



dev.off()



#}
