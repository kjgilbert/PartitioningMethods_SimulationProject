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


####################################
species <- c("lopho", "myria")
spp.pch <- c(1,2)
samp.sizes <- c(10,20,40,80)
samp.size.col <- c("steelblue4", "orange", "green4", "red")
pt.wid <- 2
opac <- 0.95

analyze.rf.datset <- function(dat, criterion){
	spp.it <- 1
	for(i in species){
		col.it <- 1
		for(j in samp.sizes){
			no.part.result <- dat[dat$species == i & dat$samp.size == j & dat$search.algo == "noPart" ,]
			rclust.AIC.result <- dat[dat$species == i & dat$samp.size == j & dat$search.algo == "rcluster" & dat$part.IC == "AIC" ,]
			rclust.BIC.result <- dat[dat$species == i & dat$samp.size == j & dat$search.algo == "rcluster" & dat$part.IC == "BIC" ,]
			rclust.AICc.result <- dat[dat$species == i & dat$samp.size == j & dat$search.algo == "rcluster" & dat$part.IC == "AICc",]
			
			if(criterion == "AIC") points(no.part.result$rf_dist - rclust.AIC.result$rf_dist, no.part.result$max.log.lik - rclust.AIC.result$max.log.lik, col=alpha(samp.size.col[col.it], opac), pch=spp.pch[spp.it], lwd=pt.wid)
			if(criterion == "BIC") points(no.part.result$rf_dist - rclust.BIC.result$rf_dist, no.part.result$max.log.lik - rclust.BIC.result$max.log.lik, col=alpha(samp.size.col[col.it], opac), pch=spp.pch[spp.it], lwd=pt.wid)
			if(criterion == "AICc") points(no.part.result$rf_dist - rclust.AICc.result$rf_dist, no.part.result$max.log.lik - rclust.AICc.result$max.log.lik, col=alpha(samp.size.col[col.it], opac), pch=spp.pch[spp.it], lwd=pt.wid)			
			
			col.it <- col.it + 1	
		}
		spp.it <- spp.it + 1	
	}
}


analyze.euc.datset <- function(dat, criterion){
	spp.it <- 1
	for(i in species){
		col.it <- 1
		for(j in samp.sizes){
			no.part.result <- dat[dat$species == i & dat$samp.size == j & dat$search.algo == "noPart" ,]
			rclust.AIC.result <- dat[dat$species == i & dat$samp.size == j & dat$search.algo == "rcluster" & dat$part.IC == "AIC" ,]
			rclust.BIC.result <- dat[dat$species == i & dat$samp.size == j & dat$search.algo == "rcluster" & dat$part.IC == "BIC" ,]
			rclust.AICc.result <- dat[dat$species == i & dat$samp.size == j & dat$search.algo == "rcluster" & dat$part.IC == "AICc",]
			
			if(criterion == "AIC") points(no.part.result$eucl_dist - rclust.AIC.result$eucl_dist, no.part.result$max.log.lik - rclust.AIC.result$max.log.lik, col=alpha(samp.size.col[col.it], opac), pch=spp.pch[spp.it], lwd=pt.wid)
			if(criterion == "BIC") points(no.part.result$eucl_dist - rclust.BIC.result$eucl_dist, no.part.result$max.log.lik - rclust.BIC.result$max.log.lik, col=alpha(samp.size.col[col.it], opac), pch=spp.pch[spp.it], lwd=pt.wid)
			if(criterion == "AICc") points(no.part.result$eucl_dist - rclust.AICc.result$eucl_dist, no.part.result$max.log.lik - rclust.AICc.result$max.log.lik, col=alpha(samp.size.col[col.it], opac), pch=spp.pch[spp.it], lwd=pt.wid)			
			
			col.it <- col.it + 1	
		}
		spp.it <- spp.it + 1	
	}
}

# can use any of these scores from the treebuilding in IQtree:
#   max.log.lik
#   uncon.log.lik
#   AIC
#   BIC
#   AICc
#   cons.log.lik (consensus tree LL)
#   sum.branch.lens
#   sum.int.branch.lens
euc.size.by.partitions <- function(dat, criterion, stat="max.log.lik"){
	spp.it <- 1
	for(i in species){
		col.it <- 1
		for(j in samp.sizes){
			no.part.result <- dat[dat$species == i & dat$samp.size == j & dat$search.algo == "noPart" ,]
			rclust.AIC.result <- dat[dat$species == i & dat$samp.size == j & dat$search.algo == "rcluster" & dat$part.IC == "AIC" ,]
			rclust.BIC.result <- dat[dat$species == i & dat$samp.size == j & dat$search.algo == "rcluster" & dat$part.IC == "BIC" ,]
			rclust.AICc.result <- dat[dat$species == i & dat$samp.size == j & dat$search.algo == "rcluster" & dat$part.IC == "AICc",]
			
			if(stat == "max.log.lik"){
				if(criterion == "AIC") points(no.part.result$eucl_dist - rclust.AIC.result$eucl_dist, no.part.result$max.log.lik - rclust.AIC.result$max.log.lik, col=alpha(samp.size.col[col.it], opac), pch=spp.pch[spp.it], lwd=pt.wid, cex=rclust.AIC.result$num.partitions/2)
				if(criterion == "BIC") points(no.part.result$eucl_dist - rclust.BIC.result$eucl_dist, no.part.result$max.log.lik - rclust.BIC.result$max.log.lik, col=alpha(samp.size.col[col.it], opac), pch=spp.pch[spp.it], lwd=pt.wid, cex=rclust.BIC.result$num.partitions/2)
				if(criterion == "AICc") points(no.part.result$eucl_dist - rclust.AICc.result$eucl_dist, no.part.result$max.log.lik - rclust.AICc.result$max.log.lik, col=alpha(samp.size.col[col.it], opac), pch=spp.pch[spp.it], lwd=pt.wid, cex=rclust.AICc.result$num.partitions/2)			
			}
			if(stat == "uncon.log.lik"){
				if(criterion == "AIC") points(no.part.result$eucl_dist - rclust.AIC.result$eucl_dist, no.part.result$uncon.log.lik - rclust.AIC.result$uncon.log.lik, col=alpha(samp.size.col[col.it], opac), pch=spp.pch[spp.it], lwd=pt.wid, cex=rclust.AIC.result$num.partitions/2)
				if(criterion == "BIC") points(no.part.result$eucl_dist - rclust.BIC.result$eucl_dist, no.part.result$uncon.log.lik - rclust.BIC.result$uncon.log.lik, col=alpha(samp.size.col[col.it], opac), pch=spp.pch[spp.it], lwd=pt.wid, cex=rclust.BIC.result$num.partitions/2)
				if(criterion == "AICc") points(no.part.result$eucl_dist - rclust.AICc.result$eucl_dist, no.part.result$uncon.log.lik - rclust.AICc.result$uncon.log.lik, col=alpha(samp.size.col[col.it], opac), pch=spp.pch[spp.it], lwd=pt.wid, cex=rclust.AICc.result$num.partitions/2)			
			}
			if(stat == "AIC"){
				if(criterion == "AIC") points(no.part.result$eucl_dist - rclust.AIC.result$eucl_dist, no.part.result$AIC - rclust.AIC.result$AIC, col=alpha(samp.size.col[col.it], opac), pch=spp.pch[spp.it], lwd=pt.wid, cex=rclust.AIC.result$num.partitions/2)
				if(criterion == "BIC") points(no.part.result$eucl_dist - rclust.BIC.result$eucl_dist, no.part.result$AIC - rclust.BIC.result$AIC, col=alpha(samp.size.col[col.it], opac), pch=spp.pch[spp.it], lwd=pt.wid, cex=rclust.BIC.result$num.partitions/2)
				if(criterion == "AICc") points(no.part.result$eucl_dist - rclust.AICc.result$eucl_dist, no.part.result$AIC - rclust.AICc.result$AIC, col=alpha(samp.size.col[col.it], opac), pch=spp.pch[spp.it], lwd=pt.wid, cex=rclust.AICc.result$num.partitions/2)			
			}
			if(stat == "BIC"){
				if(criterion == "AIC") points(no.part.result$eucl_dist - rclust.AIC.result$eucl_dist, no.part.result$BIC - rclust.AIC.result$BIC, col=alpha(samp.size.col[col.it], opac), pch=spp.pch[spp.it], lwd=pt.wid, cex=rclust.AIC.result$num.partitions/2)
				if(criterion == "BIC") points(no.part.result$eucl_dist - rclust.BIC.result$eucl_dist, no.part.result$BIC - rclust.BIC.result$BIC, col=alpha(samp.size.col[col.it], opac), pch=spp.pch[spp.it], lwd=pt.wid, cex=rclust.BIC.result$num.partitions/2)
				if(criterion == "AICc") points(no.part.result$eucl_dist - rclust.AICc.result$eucl_dist, no.part.result$BIC - rclust.AICc.result$BIC, col=alpha(samp.size.col[col.it], opac), pch=spp.pch[spp.it], lwd=pt.wid, cex=rclust.AICc.result$num.partitions/2)			
			}
			if(stat == "AICc"){
				if(criterion == "AIC") points(no.part.result$eucl_dist - rclust.AIC.result$eucl_dist, no.part.result$AICc - rclust.AIC.result$AICc, col=alpha(samp.size.col[col.it], opac), pch=spp.pch[spp.it], lwd=pt.wid, cex=rclust.AIC.result$num.partitions/2)
				if(criterion == "BIC") points(no.part.result$eucl_dist - rclust.BIC.result$eucl_dist, no.part.result$AICc - rclust.BIC.result$AICc, col=alpha(samp.size.col[col.it], opac), pch=spp.pch[spp.it], lwd=pt.wid, cex=rclust.BIC.result$num.partitions/2)
				if(criterion == "AICc") points(no.part.result$eucl_dist - rclust.AICc.result$eucl_dist, no.part.result$AICc - rclust.AICc.result$AICc, col=alpha(samp.size.col[col.it], opac), pch=spp.pch[spp.it], lwd=pt.wid, cex=rclust.AICc.result$num.partitions/2)			
			}			
			
			
			col.it <- col.it + 1	
		}
		spp.it <- spp.it + 1	
	}
}
####################################




# mut rate 500
dat1 <- read.csv("AllStatsOutput_New_Rep1_100-500.csv", header=TRUE, stringsAsFactors=FALSE)
dat1$uniqueID <- paste(dat1$species, dat1$samp.size, dat1$search.algo, dat1$part.IC, sep="_")
dat1b <- read.table("RF_Eucl_results_New_Rep1_100-500.txt", header=TRUE, stringsAsFactors=FALSE)
dat1b$uniqueID <- make.merge.id(dat1b)
dat <- merge(dat1, dat1b, by="uniqueID")
dat$gene.size <- 100
dat$mut.rate <- 500
dat1 <- dat

dat2 <- read.csv("AllStatsOutput_New_Rep2_100-500.csv", header=TRUE, stringsAsFactors=FALSE)
dat2$uniqueID <- paste(dat2$species, dat2$samp.size, dat2$search.algo, dat2$part.IC, sep="_")
dat2b <- read.table("RF_Eucl_results_New_Rep2_100-500.txt", header=TRUE, stringsAsFactors=FALSE)
dat2b$uniqueID <- make.merge.id(dat2b)
dat <- merge(dat2, dat2b, by="uniqueID")
dat$gene.size <- 100
dat$mut.rate <- 500
dat2 <- dat

dat3 <- read.csv("AllStatsOutput_New_Rep3_100-500.csv", header=TRUE, stringsAsFactors=FALSE)
dat3$uniqueID <- paste(dat3$species, dat3$samp.size, dat3$search.algo, dat3$part.IC, sep="_")
dat3b <- read.table("RF_Eucl_results_New_Rep3_100-500.txt", header=TRUE, stringsAsFactors=FALSE)
dat3b$uniqueID <- make.merge.id(dat3b)
dat <- merge(dat3, dat3b, by="uniqueID")
dat$gene.size <- 100
dat$mut.rate <- 500
dat3 <- dat

dat4 <- read.csv("AllStatsOutput_New_Rep4_100-500.csv", header=TRUE, stringsAsFactors=FALSE)
dat4$uniqueID <- paste(dat4$species, dat4$samp.size, dat4$search.algo, dat4$part.IC, sep="_")
dat4b <- read.table("RF_Eucl_results_New_Rep4_100-500.txt", header=TRUE, stringsAsFactors=FALSE)
dat4b$uniqueID <- make.merge.id(dat4b)
dat <- merge(dat4, dat4b, by="uniqueID")
dat$gene.size <- 100
dat$mut.rate <- 500
dat4 <- dat

dat5 <- read.csv("AllStatsOutput_New_Rep5_100-500.csv", header=TRUE, stringsAsFactors=FALSE)
dat5$uniqueID <- paste(dat5$species, dat5$samp.size, dat5$search.algo, dat5$part.IC, sep="_")
dat5b <- read.table("RF_Eucl_results_New_Rep5_100-500.txt", header=TRUE, stringsAsFactors=FALSE)
dat5b$uniqueID <- make.merge.id(dat5b)
dat <- merge(dat5, dat5b, by="uniqueID")
dat$gene.size <- 100
dat$mut.rate <- 500
dat5 <- dat


# mut rate 1000
dat6 <- read.csv("AllStatsOutput_New_Rep1_100-1000.csv", header=TRUE, stringsAsFactors=FALSE)
dat6$uniqueID <- paste(dat6$species, dat6$samp.size, dat6$search.algo, dat6$part.IC, sep="_")
dat6b <- read.table("RF_Eucl_results_New_Rep1_100-1000.txt", header=TRUE, stringsAsFactors=FALSE)
dat6b$uniqueID <- make.merge.id(dat6b)
dat <- merge(dat6, dat6b, by="uniqueID")
dat$gene.size <- 100
dat$mut.rate <- 1000
dat6 <- dat

dat7 <- read.csv("AllStatsOutput_New_Rep2_100-1000.csv", header=TRUE, stringsAsFactors=FALSE)
dat7$uniqueID <- paste(dat7$species, dat7$samp.size, dat7$search.algo, dat7$part.IC, sep="_")
dat7b <- read.table("RF_Eucl_results_New_Rep2_100-1000.txt", header=TRUE, stringsAsFactors=FALSE)
dat7b$uniqueID <- make.merge.id(dat7b)
dat <- merge(dat7, dat7b, by="uniqueID")
dat$gene.size <- 100
dat$mut.rate <- 1000
dat7 <- dat

dat8 <- read.csv("AllStatsOutput_New_Rep3_100-1000.csv", header=TRUE, stringsAsFactors=FALSE)
dat8$uniqueID <- paste(dat8$species, dat8$samp.size, dat8$search.algo, dat8$part.IC, sep="_")
dat8b <- read.table("RF_Eucl_results_New_Rep3_100-1000.txt", header=TRUE, stringsAsFactors=FALSE)
dat8b$uniqueID <- make.merge.id(dat8b)
dat <- merge(dat8, dat8b, by="uniqueID")
dat$gene.size <- 100
dat$mut.rate <- 1000
dat8 <- dat

dat9 <- read.csv("AllStatsOutput_New_Rep4_100-1000.csv", header=TRUE, stringsAsFactors=FALSE)
dat9$uniqueID <- paste(dat9$species, dat9$samp.size, dat9$search.algo, dat9$part.IC, sep="_")
dat9b <- read.table("RF_Eucl_results_New_Rep4_100-1000.txt", header=TRUE, stringsAsFactors=FALSE)
dat9b$uniqueID <- make.merge.id(dat9b)
dat <- merge(dat9, dat9b, by="uniqueID")
dat$gene.size <- 100
dat$mut.rate <- 1000
dat9 <- dat

dat10 <- read.csv("AllStatsOutput_New_Rep5_100-1000.csv", header=TRUE, stringsAsFactors=FALSE)
dat10$uniqueID <- paste(dat10$species, dat10$samp.size, dat10$search.algo, dat10$part.IC, sep="_")
dat10b <- read.table("RF_Eucl_results_New_Rep5_100-1000.txt", header=TRUE, stringsAsFactors=FALSE)
dat10b$uniqueID <- make.merge.id(dat10b)
dat <- merge(dat10, dat10b, by="uniqueID")
dat$gene.size <- 100
dat$mut.rate <- 1000
dat10 <- dat


########


##plot(0, type="n", xlim=c(-0.5,0.5), ylim=c(-50,50))
##analyze.rf.datset(dat1, "AIC")

########


pdf("Results_Plots/CompareLogLikelihoods.pdf", width=8, height=4)
par(mfrow=c(1,3))

opac <- 0.75
pt.lwd <- 2
xlimits <- c(-0.9,0.9)
ylimits <- c(-52.5,3)
plot.stat <- "max.log.lik"
plot(0, type="n", xlim=xlimits, ylim=ylimits, xlab="no part - part eucl dist from true tree", ylab="no part - part max log likelihood", main="partitioning w/ AIC")
abline(h=0, col="gray75", lty=2); abline(v=0, col="gray75", lty=2)
#analyze.euc.datset(dat1, "AIC")
text(-0.7,3.5, "no partitioning", cex=0.6)
text(-0.7,2, "higher max LL", cex=0.6)
text(-0.7,-2, "partitioning", cex=0.6)
text(-0.7,-3.5, "higher max LL", cex=0.6)
euc.size.by.partitions(dat1, "AIC", stat=plot.stat)
euc.size.by.partitions(dat2, "AIC", stat=plot.stat)
euc.size.by.partitions(dat3, "AIC", stat=plot.stat)
euc.size.by.partitions(dat4, "AIC", stat=plot.stat)
euc.size.by.partitions(dat5, "AIC", stat=plot.stat)
euc.size.by.partitions(dat6, "AIC", stat=plot.stat)
euc.size.by.partitions(dat7, "AIC", stat=plot.stat)
euc.size.by.partitions(dat8, "AIC", stat=plot.stat)
euc.size.by.partitions(dat9, "AIC", stat=plot.stat)
euc.size.by.partitions(dat10, "AIC", stat=plot.stat)
legend("bottomleft", c("lopho", "myria"), species, pch=spp.pch, col="black", bg="white", bty="n")

plot(0, type="n", xlim=xlimits, ylim=ylimits, xlab="no part - part eucl dist from true tree", ylab="no part - part max log likelihood", main="partitioning w/ BIC")
abline(h=0, col="gray75", lty=2); abline(v=0, col="gray75", lty=2)
text(-0.5,-49, "no partitioning", cex=0.7)
text(-0.5,-51, "closer to true tree", cex=0.7)
text(0.5,-49, "partitioning", cex=0.7)
text(0.5,-51, "closer to true tree", cex=0.7)
euc.size.by.partitions(dat1, "BIC", stat=plot.stat)
euc.size.by.partitions(dat2, "BIC", stat=plot.stat)
euc.size.by.partitions(dat3, "BIC", stat=plot.stat)
euc.size.by.partitions(dat4, "BIC", stat=plot.stat)
euc.size.by.partitions(dat5, "BIC", stat=plot.stat)
euc.size.by.partitions(dat6, "BIC", stat=plot.stat)
euc.size.by.partitions(dat7, "BIC", stat=plot.stat)
euc.size.by.partitions(dat8, "BIC", stat=plot.stat)
euc.size.by.partitions(dat9, "BIC", stat=plot.stat)
euc.size.by.partitions(dat10, "BIC", stat=plot.stat)

plot(0, type="n", xlim=xlimits, ylim=ylimits, xlab="no part - part eucl dist from true tree", ylab="no part - part max log likelihood", main="partitioning w/ AICc")
abline(h=0, col="gray75", lty=2); abline(v=0, col="gray75", lty=2)
euc.size.by.partitions(dat1, "AICc", stat=plot.stat)
euc.size.by.partitions(dat2, "AICc", stat=plot.stat)
euc.size.by.partitions(dat3, "AICc", stat=plot.stat)
euc.size.by.partitions(dat4, "AICc", stat=plot.stat)
euc.size.by.partitions(dat5, "AICc", stat=plot.stat)
euc.size.by.partitions(dat6, "AICc", stat=plot.stat)
euc.size.by.partitions(dat7, "AICc", stat=plot.stat)
euc.size.by.partitions(dat8, "AICc", stat=plot.stat)
euc.size.by.partitions(dat9, "AICc", stat=plot.stat)
euc.size.by.partitions(dat10, "AICc", stat=plot.stat)

legend("bottomleft", c("10", "20", "40", "80"), pch=15, col=samp.size.col, title = "num. genes", bg="white", cex=0.9, pt.cex=1.5)
dev.off()
