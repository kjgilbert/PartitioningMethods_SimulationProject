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
samp.size.col <- c("steelblue4", "darkorange2", "green4", "red3")
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
				if(criterion == "AIC") points(no.part.result$max.log.lik - rclust.AIC.result$max.log.lik, no.part.result$eucl_dist - rclust.AIC.result$eucl_dist, col=alpha(samp.size.col[col.it], opac), pch=spp.pch[spp.it], lwd=pt.wid, cex=rclust.AIC.result$num.partitions/2)
				if(criterion == "BIC") points(no.part.result$max.log.lik - rclust.BIC.result$max.log.lik, no.part.result$eucl_dist - rclust.BIC.result$eucl_dist, col=alpha(samp.size.col[col.it], opac), pch=spp.pch[spp.it], lwd=pt.wid, cex=rclust.BIC.result$num.partitions/2)
				if(criterion == "AICc") points(no.part.result$max.log.lik - rclust.AICc.result$max.log.lik, no.part.result$eucl_dist - rclust.AICc.result$eucl_dist, col=alpha(samp.size.col[col.it], opac), pch=spp.pch[spp.it], lwd=pt.wid, cex=rclust.AICc.result$num.partitions/2)			
			}
			if(stat == "uncon.log.lik"){
				if(criterion == "AIC") points(no.part.result$uncon.log.lik - rclust.AIC.result$uncon.log.lik, no.part.result$eucl_dist - rclust.AIC.result$eucl_dist, col=alpha(samp.size.col[col.it], opac), pch=spp.pch[spp.it], lwd=pt.wid, cex=rclust.AIC.result$num.partitions/2)
				if(criterion == "BIC") points(no.part.result$uncon.log.lik - rclust.BIC.result$uncon.log.lik, no.part.result$eucl_dist - rclust.BIC.result$eucl_dist, col=alpha(samp.size.col[col.it], opac), pch=spp.pch[spp.it], lwd=pt.wid, cex=rclust.BIC.result$num.partitions/2)
				if(criterion == "AICc") points(no.part.result$uncon.log.lik - rclust.AICc.result$uncon.log.lik, no.part.result$eucl_dist - rclust.AICc.result$eucl_dist, col=alpha(samp.size.col[col.it], opac), pch=spp.pch[spp.it], lwd=pt.wid, cex=rclust.AICc.result$num.partitions/2)			
			}		
			
			col.it <- col.it + 1	
		}
		spp.it <- spp.it + 1	
	}
}


rf_dist.size.by.partitions <- function(dat, criterion, stat="max.log.lik"){
	spp.it <- 1
	for(i in species){
		col.it <- 1
		for(j in samp.sizes){
			no.part.result <- dat[dat$species == i & dat$samp.size == j & dat$search.algo == "noPart" ,]
			rclust.AIC.result <- dat[dat$species == i & dat$samp.size == j & dat$search.algo == "rcluster" & dat$part.IC == "AIC" ,]
			rclust.BIC.result <- dat[dat$species == i & dat$samp.size == j & dat$search.algo == "rcluster" & dat$part.IC == "BIC" ,]
			rclust.AICc.result <- dat[dat$species == i & dat$samp.size == j & dat$search.algo == "rcluster" & dat$part.IC == "AICc",]
			
			if(stat == "max.log.lik"){
				if(criterion == "AIC") points(no.part.result$max.log.lik - rclust.AIC.result$max.log.lik, no.part.result$rf_dist - rclust.AIC.result$rf_dist, col=alpha(samp.size.col[col.it], opac), pch=spp.pch[spp.it], lwd=pt.wid, cex=rclust.AIC.result$num.partitions/2)
				if(criterion == "BIC") points(no.part.result$max.log.lik - rclust.BIC.result$max.log.lik, no.part.result$rf_dist - rclust.BIC.result$rf_dist, col=alpha(samp.size.col[col.it], opac), pch=spp.pch[spp.it], lwd=pt.wid, cex=rclust.BIC.result$num.partitions/2)
				if(criterion == "AICc") points(no.part.result$max.log.lik - rclust.AICc.result$max.log.lik, no.part.result$rf_dist - rclust.AICc.result$rf_dist, col=alpha(samp.size.col[col.it], opac), pch=spp.pch[spp.it], lwd=pt.wid, cex=rclust.AICc.result$num.partitions/2)			
			}
			if(stat == "uncon.log.lik"){
				if(criterion == "AIC") points(no.part.result$uncon.log.lik - rclust.AIC.result$uncon.log.lik, no.part.result$rf_dist - rclust.AIC.result$rf_dist, col=alpha(samp.size.col[col.it], opac), pch=spp.pch[spp.it], lwd=pt.wid, cex=rclust.AIC.result$num.partitions/2)
				if(criterion == "BIC") points(no.part.result$uncon.log.lik - rclust.BIC.result$uncon.log.lik, no.part.result$rf_dist - rclust.BIC.result$rf_dist, col=alpha(samp.size.col[col.it], opac), pch=spp.pch[spp.it], lwd=pt.wid, cex=rclust.BIC.result$num.partitions/2)
				if(criterion == "AICc") points(no.part.result$uncon.log.lik - rclust.AICc.result$uncon.log.lik, no.part.result$rf_dist - rclust.AICc.result$rf_dist, col=alpha(samp.size.col[col.it], opac), pch=spp.pch[spp.it], lwd=pt.wid, cex=rclust.AICc.result$num.partitions/2)			
			}		
			
			col.it <- col.it + 1	
		}
		spp.it <- spp.it + 1	
	}
}
####################################





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

dat11 <- read.csv("AllStatsOutput_New_Rep6_100-1000.csv", header=TRUE, stringsAsFactors=FALSE)
dat11$uniqueID <- paste(dat11$species, dat11$samp.size, dat11$search.algo, dat11$part.IC, sep="_")
dat11b <- read.table("RF_Eucl_results_New_Rep6_100-1000.txt", header=TRUE, stringsAsFactors=FALSE)
dat11b$uniqueID <- make.merge.id(dat11b)
dat <- merge(dat11, dat11b, by="uniqueID")
dat$gene.size <- 100
dat$mut.rate <- 1000
dat11 <- dat

dat12 <- read.csv("AllStatsOutput_Newest_Rep7_100-1000.csv", header=TRUE, stringsAsFactors=FALSE)
dat12$uniqueID <- paste(dat12$species, dat12$samp.size, dat12$search.algo, dat12$part.IC, sep="_")
dat12b <- read.table("RF_Eucl_results_New_Rep7_100-1000.txt", header=TRUE, stringsAsFactors=FALSE)
dat12b$uniqueID <- make.merge.id(dat12b)
dat <- merge(dat12, dat12b, by="uniqueID")
dat$gene.size <- 100
dat$mut.rate <- 1000
dat12 <- dat

dat13 <- read.csv("AllStatsOutput_Newest_Rep8_100-1000.csv", header=TRUE, stringsAsFactors=FALSE)
dat13$uniqueID <- paste(dat13$species, dat13$samp.size, dat13$search.algo, dat13$part.IC, sep="_")
dat13b <- read.table("RF_Eucl_results_New_Rep8_100-1000.txt", header=TRUE, stringsAsFactors=FALSE)
dat13b$uniqueID <- make.merge.id(dat13b)
dat <- merge(dat13, dat13b, by="uniqueID")
dat$gene.size <- 100
dat$mut.rate <- 1000
dat13 <- dat

dat14 <- read.csv("AllStatsOutput_Newest_Rep9_100-1000.csv", header=TRUE, stringsAsFactors=FALSE)
dat14$uniqueID <- paste(dat14$species, dat14$samp.size, dat14$search.algo, dat14$part.IC, sep="_")
dat14b <- read.table("RF_Eucl_results_New_Rep9_100-1000.txt", header=TRUE, stringsAsFactors=FALSE)
dat14b$uniqueID <- make.merge.id(dat14b)
dat <- merge(dat14, dat14b, by="uniqueID")
dat$gene.size <- 100
dat$mut.rate <- 1000
dat14 <- dat

dat15 <- read.csv("AllStatsOutput_Newest_Rep10_100-1000.csv", header=TRUE, stringsAsFactors=FALSE)
dat15$uniqueID <- paste(dat15$species, dat15$samp.size, dat15$search.algo, dat15$part.IC, sep="_")
dat15b <- read.table("RF_Eucl_results_New_Rep10_100-1000.txt", header=TRUE, stringsAsFactors=FALSE)
dat15b$uniqueID <- make.merge.id(dat15b)
dat <- merge(dat15, dat15b, by="uniqueID")
dat$gene.size <- 100
dat$mut.rate <- 1000
dat15 <- dat

dat16 <- read.csv("AllStatsOutput_Newest_Rep11_100-1000.csv", header=TRUE, stringsAsFactors=FALSE)
dat16$uniqueID <- paste(dat16$species, dat16$samp.size, dat16$search.algo, dat16$part.IC, sep="_")
dat16b <- read.table("RF_Eucl_results_New_Rep11_100-1000.txt", header=TRUE, stringsAsFactors=FALSE)
dat16b$uniqueID <- make.merge.id(dat16b)
dat <- merge(dat16, dat16b, by="uniqueID")
dat$gene.size <- 100
dat$mut.rate <- 1000
dat16 <- dat

dat17 <- read.csv("AllStatsOutput_Newest_Rep12_100-1000.csv", header=TRUE, stringsAsFactors=FALSE)
dat17$uniqueID <- paste(dat17$species, dat17$samp.size, dat17$search.algo, dat17$part.IC, sep="_")
dat17b <- read.table("RF_Eucl_results_New_Rep12_100-1000.txt", header=TRUE, stringsAsFactors=FALSE)
dat17b$uniqueID <- make.merge.id(dat17b)
dat <- merge(dat17, dat17b, by="uniqueID")
dat$gene.size <- 100
dat$mut.rate <- 1000
dat17 <- dat

########
########




plot.stat <- "max.log.lik"
part.text1 <- -35
part.text2 <- 0.65
part.text3 <- -0.65


pdf("Results_Plots/Flipped_CompareLogLikelihoods_BICpartitioning.pdf", width=8, height=4)
par(mfrow=c(1,2), mar=c(4,4,1,0.5)) # if just looking at BIC for partitionfinder
opac <- 0.75
pt.lwd <- 2
xlimits <- c(-60,3)




# robinson foulds distance difference

ylimits <- c(-0.25,0.25)

plot(0, type="n", xlim=xlimits, ylim=ylimits, ylab="", xlab="diff in IQtree max log lik of 'no part' - 'part'", main="")
mtext(side=2, at=0, line=2.5, "diff in RF dist from true tree: 'no part' - 'part'")
abline(h=0, col="gray75", lty=2); abline(v=0, col="gray75", lty=2)
text(part.text1, part.text2, expression(italic("partitioning better")), cex=1, col="green4")
text(part.text1, part.text3, expression(italic("no partitioning better")), cex=1, col="red3")
rf_dist.size.by.partitions(dat6, "BIC", stat=plot.stat)
rf_dist.size.by.partitions(dat7, "BIC", stat=plot.stat)
rf_dist.size.by.partitions(dat8, "BIC", stat=plot.stat)
rf_dist.size.by.partitions(dat9, "BIC", stat=plot.stat)
rf_dist.size.by.partitions(dat10, "BIC", stat=plot.stat)
rf_dist.size.by.partitions(dat11, "BIC", stat=plot.stat)
rf_dist.size.by.partitions(dat12, "BIC", stat=plot.stat)
rf_dist.size.by.partitions(dat13, "BIC", stat=plot.stat)
rf_dist.size.by.partitions(dat14, "BIC", stat=plot.stat)
rf_dist.size.by.partitions(dat15, "BIC", stat=plot.stat)
rf_dist.size.by.partitions(dat16, "BIC", stat=plot.stat)
rf_dist.size.by.partitions(dat17, "BIC", stat=plot.stat)

legend("bottomleft", c("10", "20", "40", "80"), pch=15, col=samp.size.col, title = "num. genes", bg="white", cex=0.9, pt.cex=1.5)
legend("topleft", c("lopho", "myria"), species, pch=spp.pch, col="black", bg="white", bty="n")


#__________________________________________________________________________________________#


# euclidean distance difference

ylimits <- c(-0.8,0.8)

plot(0, type="n", xlim=xlimits, ylim=ylimits, ylab="", xlab="diff in IQtree max log lik of 'no part' - 'part'", main="")
mtext(side=2, at=0, line=2.5, "diff in eucl dist from true tree: 'no part' - 'part'")
abline(h=0, col="gray75", lty=2); abline(v=0, col="gray75", lty=2)
text(part.text1, part.text2, expression(italic("partitioning better")), cex=1, col="green4")
text(part.text1, part.text3, expression(italic("no partitioning better")), cex=1, col="red3")
euc.size.by.partitions(dat6, "BIC", stat=plot.stat)
euc.size.by.partitions(dat7, "BIC", stat=plot.stat)
euc.size.by.partitions(dat8, "BIC", stat=plot.stat)
euc.size.by.partitions(dat9, "BIC", stat=plot.stat)
euc.size.by.partitions(dat10, "BIC", stat=plot.stat)
euc.size.by.partitions(dat11, "BIC", stat=plot.stat)
euc.size.by.partitions(dat12, "BIC", stat=plot.stat)
euc.size.by.partitions(dat13, "BIC", stat=plot.stat)
euc.size.by.partitions(dat14, "BIC", stat=plot.stat)
euc.size.by.partitions(dat15, "BIC", stat=plot.stat)
euc.size.by.partitions(dat16, "BIC", stat=plot.stat)
euc.size.by.partitions(dat17, "BIC", stat=plot.stat)


dev.off()


