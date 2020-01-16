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





dat1 <- read.csv("AllStatsOutput_Rep1_genes100.csv", header=TRUE, stringsAsFactors=FALSE)
dat1$uniqueID <- paste(dat1$species, dat1$samp.size, dat1$search.algo, dat1$part.IC, sep="_")
dat1b <- read.table("RF_Eucl_results_Rep1_genes100.txt", header=TRUE, stringsAsFactors=FALSE)
dat1b$uniqueID <- make.merge.id(dat1b)
dat <- merge(dat1, dat1b, by="uniqueID")
dat$gene.size <- 100
dat$mut.rate <- 250
dat1 <- dat

dat2 <- read.csv("AllStatsOutput_Rep1_genes100_mutRate500.csv", header=TRUE, stringsAsFactors=FALSE)
dat2$uniqueID <- paste(dat2$species, dat2$samp.size, dat2$search.algo, dat2$part.IC, sep="_")
dat2b <- read.table("RF_Eucl_results_Rep1_genes100_mutRate500.txt", header=TRUE, stringsAsFactors=FALSE)
dat2b$uniqueID <- make.merge.id(dat2b)
dat <- merge(dat2, dat2b, by="uniqueID")
dat$gene.size <- 100
dat$mut.rate <- 500
dat2 <- dat

dat3 <- read.csv("AllStatsOutput_Rep2_genes100.csv", header=TRUE, stringsAsFactors=FALSE)
dat3$uniqueID <- paste(dat3$species, dat3$samp.size, dat3$search.algo, dat3$part.IC, sep="_")
dat3b <- read.table("RF_Eucl_results_Rep2_genes100.txt", header=TRUE, stringsAsFactors=FALSE)
dat3b$uniqueID <- make.merge.id(dat3b)
dat <- merge(dat3, dat3b, by="uniqueID")
dat$gene.size <- 100
dat$mut.rate <- 250
dat3 <- dat

dat4 <- read.csv("AllStatsOutput_Rep2_genes100_mutRate500.csv", header=TRUE, stringsAsFactors=FALSE)
dat4$uniqueID <- paste(dat4$species, dat4$samp.size, dat4$search.algo, dat4$part.IC, sep="_")
dat4b <- read.table("RF_Eucl_results_Rep2_genes100_mutRate500.txt", header=TRUE, stringsAsFactors=FALSE)
dat4b$uniqueID <- make.merge.id(dat4b)
dat <- merge(dat4, dat4b, by="uniqueID")
dat$gene.size <- 100
dat$mut.rate <- 500
dat4 <- dat

dat5 <- read.csv("AllStatsOutput_Rep3_genes100.csv", header=TRUE, stringsAsFactors=FALSE)
dat5$uniqueID <- paste(dat5$species, dat5$samp.size, dat5$search.algo, dat5$part.IC, sep="_")
dat5b <- read.table("RF_Eucl_results_Rep3_genes100.txt", header=TRUE, stringsAsFactors=FALSE)
dat5b$uniqueID <- make.merge.id(dat5b)
dat <- merge(dat5, dat5b, by="uniqueID")
dat$gene.size <- 100
dat$mut.rate <- 250
dat5 <- dat

dat6 <- read.csv("AllStatsOutput_Rep3_genes100_mutRate500.csv", header=TRUE, stringsAsFactors=FALSE)
dat6$uniqueID <- paste(dat6$species, dat6$samp.size, dat6$search.algo, dat6$part.IC, sep="_")
dat6b <- read.table("RF_Eucl_results_Rep3_genes100_mutRate500.txt", header=TRUE, stringsAsFactors=FALSE)
dat6b$uniqueID <- make.merge.id(dat6b)
dat <- merge(dat6, dat6b, by="uniqueID")
dat$gene.size <- 100
dat$mut.rate <- 500
dat6 <- dat

dat7 <- read.csv("AllStatsOutput_Rep4_genes100.csv", header=TRUE, stringsAsFactors=FALSE)
dat7$uniqueID <- paste(dat7$species, dat7$samp.size, dat7$search.algo, dat7$part.IC, sep="_")
dat7b <- read.table("RF_Eucl_results_Rep4_genes100.txt", header=TRUE, stringsAsFactors=FALSE)
dat7b$uniqueID <- make.merge.id(dat7b)
dat <- merge(dat7, dat7b, by="uniqueID")
dat$gene.size <- 100
dat$mut.rate <- 250
dat7 <- dat

dat8 <- read.csv("AllStatsOutput_Rep5_genes100.csv", header=TRUE, stringsAsFactors=FALSE)
dat8$uniqueID <- paste(dat8$species, dat8$samp.size, dat8$search.algo, dat8$part.IC, sep="_")
dat8b <- read.table("RF_Eucl_results_Rep5_genes100.txt", header=TRUE, stringsAsFactors=FALSE)
dat8b$uniqueID <- make.merge.id(dat8b)
dat <- merge(dat8, dat8b, by="uniqueID")
dat$gene.size <- 100
dat$mut.rate <- 250
dat8 <- dat


stats.dat <- rbind(dat1, dat2, dat3, dat4, dat5, dat6, dat7, dat8)



### FOR NOW ONLY TAKE N=10 SAMPLE SIZE DAT
stats.dat <- stats.dat[stats.dat$samp.size == 10 , ]

# stats.dat <- stats.dat[-131,]

uniq.names <- names(table(stats.dat$uniqueID))

####################################


## make some point types for plotting
stats.dat$plot.pch <- NA
stats.dat$plot.pch[which(stats.dat$samp.size == 10 & stats.dat$search.algo == "rcluster")] <- 0
stats.dat$plot.pch[which(stats.dat$samp.size == 20 & stats.dat$search.algo == "rcluster")] <- 1
stats.dat$plot.pch[which(stats.dat$samp.size == 40 & stats.dat$search.algo == "rcluster")] <- 2
stats.dat$plot.pch[which(stats.dat$samp.size == 80 & stats.dat$search.algo == "rcluster")] <- 6
stats.dat$plot.pch[which(stats.dat$samp.size == 10 & stats.dat$search.algo == "greedy")] <- 15
stats.dat$plot.pch[which(stats.dat$samp.size == 20 & stats.dat$search.algo == "greedy")] <- 19
stats.dat$plot.pch[which(stats.dat$samp.size == 40 & stats.dat$search.algo == "greedy")] <- 17
stats.dat$plot.pch[which(stats.dat$samp.size == 80 & stats.dat$search.algo == "greedy")] <- 25
stats.dat$plot.pch[which(stats.dat$samp.size == 10 & stats.dat$search.algo == "noPart")] <- 0
stats.dat$plot.pch[which(stats.dat$samp.size == 20 & stats.dat$search.algo == "noPart")] <- 1
stats.dat$plot.pch[which(stats.dat$samp.size == 40 & stats.dat$search.algo == "noPart")] <- 2
stats.dat$plot.pch[which(stats.dat$samp.size == 80 & stats.dat$search.algo == "noPart")] <- 6



# for when plotting only n=10
stats.dat$plot.pch[which(stats.dat$samp.size == 10 & stats.dat$search.algo == "rcluster" & stats.dat$species == "lopho")] <- 0
stats.dat$plot.pch[which(stats.dat$samp.size == 10 & stats.dat$search.algo == "rcluster" & stats.dat$species == "myria")] <- 15
stats.dat$plot.pch[which(stats.dat$samp.size == 10 & stats.dat$search.algo == "greedy" & stats.dat$species == "lopho")] <- 1
stats.dat$plot.pch[which(stats.dat$samp.size == 10 & stats.dat$search.algo == "greedy" & stats.dat$species == "myria")] <- 19


AIC.dat <- stats.dat[stats.dat$part.IC == "AIC", ]
AICc.dat <- stats.dat[stats.dat$part.IC == "AICc", ]
BIC.dat <- stats.dat[stats.dat$part.IC == "BIC", ]
noPart.dat <- stats.dat[stats.dat$part.IC == "noPart", ]





stats.dat$part.IC.score[stats.dat$part.IC.score == "noPart"] <- NA
stats.dat$part.IC.score <- as.numeric(stats.dat$part.IC.score)



pdf("AICc_treebuilding.pdf", width=5, height=4)
par(mar=c(4,4,0.85,0.5)); size <- 1.75
plot(stats.dat$num.partitions, stats.dat$AICc, xlab="'best' number of partitions", ylab="treebuilding AICc", col="white", log="y")
points(AIC.dat$num.partitions+0.1, AIC.dat$AICc, col="red", pch=AIC.dat$plot.pch, lwd=size)# lwd=AIC.dat$samp.size/sizer, 
points(AICc.dat$num.partitions, AICc.dat$AICc, col="darkred", pch=AICc.dat$plot.pch, lwd=size)
points(BIC.dat$num.partitions-0.1, BIC.dat$AICc, col="darkblue", pch=BIC.dat$plot.pch, lwd=size)
points(jitter(noPart.dat$num.partitions, factor=3.5), noPart.dat$AICc, col="green3", pch=noPart.dat$plot.pch, lwd=size)
legend("topright", c("none", "BIC", "AIC", "AICc"), title="Partitioning criterion", col=c("green3", "darkblue", "red", "darkred"), pch=15, cex=0.95, bty="n", ncol=2, pt.lwd=3)
#legend(1,2.2e6, c("10", "20", "40", "80"), title="Num. genes sampled", col=c("black"), pch=1, cex=0.95, bty="n", ncol=2, pt.lwd=c(10/sizer, 20/sizer, 40/sizer, 80/sizer))
legend(1,2.25e6, c("10", "20", "40", "80"), pch=c(0,1,2,6), title="Num. genes sampled", bty="n", ncol=2, pt.lwd=2, cex=0.95)
dev.off()



pdf("likelihoods_treebuilding.pdf", width=5, height=4)
par(mar=c(4,4,0.85,0.5)); size <- 1.75
plot(stats.dat$num.partitions, stats.dat$max.log.lik, ylim=c(-150000, -80000), xlab="'best' number of partitions", ylab="max. log likelihood", col="white", xaxt="n")
axis(side=1, at=c(0:7), c(" ", "1", "2", "3", "4", "5", "6", "7"))
mtext(side=1, at=0, "no", line=0.5)
mtext(side=1, at=0, "partitioning", line=1.5)
points(jitter(AIC.dat$num.partitions+0.1, factor=0.1), as.numeric(AIC.dat$max.log.lik), col="red", pch=AIC.dat$plot.pch, lwd=size)# lwd=AIC.dat$samp.size/sizer, 
points(jitter(AICc.dat$num.partitions, factor=0.1), AICc.dat$max.log.lik, col="darkred", pch=AICc.dat$plot.pch, lwd=size)
points(jitter(BIC.dat$num.partitions-0.1, factor=0.1), BIC.dat$max.log.lik, col="darkblue", pch=BIC.dat$plot.pch, lwd=size)
points(jitter(noPart.dat$num.partitions, factor=3.5), noPart.dat$max.log.lik, col="green3", pch=noPart.dat$plot.pch, lwd=size)
legend("bottomright", c("none", "BIC", "AIC", "AICc"), title="Partitioning criterion", col=c("green3", "darkblue", "red", "darkred"), pch=15, cex=0.95, bty="n", ncol=2, pt.lwd=3)


plot(stats.dat$num.partitions, stats.dat$max.log.lik, ylim=c(-150000, -30000), xlab="'best' number of partitions", ylab="max. log likelihood", col="white", xaxt="n")
axis(side=1, at=c(0:7), c(" ", "1", "2", "3", "4", "5", "6", "7"))
mtext(side=1, at=0, "no", line=0.5)
mtext(side=1, at=0, "partitioning", line=1.5)
points(jitter(AIC.dat$num.partitions+0.1, factor=0.1), as.numeric(AIC.dat$max.log.lik), col="red", pch=AIC.dat$plot.pch, lwd=size)# lwd=AIC.dat$samp.size/sizer, 
points(jitter(AICc.dat$num.partitions, factor=0.1), AICc.dat$max.log.lik, col="darkred", pch=AICc.dat$plot.pch, lwd=size)
points(jitter(BIC.dat$num.partitions-0.1, factor=0.1), BIC.dat$max.log.lik, col="darkblue", pch=BIC.dat$plot.pch, lwd=size)
points(jitter(noPart.dat$num.partitions, factor=3.5), noPart.dat$max.log.lik, col="green3", pch=noPart.dat$plot.pch, lwd=size)
legend("bottomright", c("none", "BIC", "AIC", "AICc"), title="Partitioning criterion", col=c("green3", "darkblue", "red", "darkred"), pch=15, cex=0.95, bty="n", ncol=2, pt.lwd=3)


plot(stats.dat$num.partitions, stats.dat$max.log.lik, xlab="'best' number of partitions", ylab="max. log likelihood", col="white", xlim=c(0,7), ylim=c(-1e6, -40000), xaxt="n")
axis(side=1, at=c(0:7), c(" ", "1", "2", "3", "4", "5", "6", "7"))
mtext(side=1, at=0, "no", line=0.5)
mtext(side=1, at=0, "partitioning", line=1.5)
points(AIC.dat$num.partitions+0.1, as.numeric(AIC.dat$max.log.lik), col="red", pch=AIC.dat$plot.pch, lwd=size)# lwd=AIC.dat$samp.size/sizer, 
points(AICc.dat$num.partitions, AICc.dat$max.log.lik, col="darkred", pch=AICc.dat$plot.pch, lwd=size)
points(BIC.dat$num.partitions-0.1, BIC.dat$max.log.lik, col="darkblue", pch=BIC.dat$plot.pch, lwd=size)
points(jitter(noPart.dat$num.partitions, factor=3.5), noPart.dat$max.log.lik, col="green3", pch=noPart.dat$plot.pch, lwd=size)
legend("bottomright", c("none", "BIC", "AIC", "AICc"), title="Partitioning criterion", col=c("green3", "darkblue", "red", "darkred"), pch=15, cex=0.95, bty="n", ncol=2, pt.lwd=3)
#legend(1,2.2e6, c("10", "20", "40", "80"), title="Num. genes sampled", col=c("black"), pch=1, cex=0.95, bty="n", ncol=2, pt.lwd=c(10/sizer, 20/sizer, 40/sizer, 80/sizer))
##legend(1,-7.6e5, c("10", "20", "40", "80"), pch=c(0,1,2,6), title="Num. genes sampled", bty="n", ncol=2, pt.lwd=2, cex=0.95)
dev.off()





plot(stats.dat$num.partitions, stats.dat$AIC, xlab="number of partitions", ylab="treebuilding AIC", col="white")
points(AIC.dat$num.partitions, AIC.dat$AIC, col="red", lwd=2)
points(AICc.dat$num.partitions, AICc.dat$AIC, col="darkred", lwd=2)
points(BIC.dat$num.partitions, BIC.dat$AIC, col="darkblue", lwd=2)
points(noPart.dat$num.partitions, noPart.dat$AIC, col="green3", lwd=2)





pdf("distances.pdf", width=7, height=4)
par(mfrow=c(1,2), mar=c(3,3,1,0.5))
plot(0, type="n", xlim=c(0.5,4.5), ylim=c(0,0.25), xaxt="n", xlab="", ylab="Topological distance", mgp=c(2.15,1,0))
axis(side=1, at=c(1:4), c("", "AIC", "AICc", "BIC"))
mtext(side=1, at=c(1:4), c("no", "", "", ""), line=0.5)
mtext(side=1, at=c(1:4), c("partitioning", "", "", ""), line=1.5)
boxplot(noPart.dat$rf_dist, add=TRUE, at=1, pch=20, col="green3", border="black")
boxplot(AIC.dat$rf_dist, add=TRUE, at=2, pch=20, col="red", border="black")
boxplot(AIC.dat$rf_dist, add=TRUE, at=3, pch=20, col="darkred", border="black")
boxplot(AIC.dat$rf_dist, add=TRUE, at=4, pch=20, col="darkblue", border="black")
plot(0, type="n", xlim=c(0.5,4.5), ylim=c(1,3), xaxt="n", xlab="", ylab="Top. + branch length distance", mgp=c(2.15,1,0))
axis(side=1, at=c(1:4), c("", "AIC", "AICc", "BIC"))
mtext(side=1, at=c(1:4), c("no", "", "", ""), line=0.5)
mtext(side=1, at=c(1:4), c("partitioning", "", "", ""), line=1.5)
boxplot(noPart.dat$eucl_dist/noPart.dat$mut.rate, add=TRUE, at=1, pch=20, col="green3", border="black")
boxplot(AIC.dat$eucl_dist/AIC.dat$mut.rate, add=TRUE, at=2, pch=20, col="red", border="black")
boxplot(AICc.dat$eucl_dist/AICc.dat$mut.rate, add=TRUE, at=3, pch=20, col="darkred", border="black")
boxplot(BIC.dat$eucl_dist/BIC.dat$mut.rate, add=TRUE, at=4, pch=20, col="darkblue", border="black")
dev.off()









par(mar=c(4,4,0.85,0.5), mfrow=c(1,1)); sizer <- 6

plot(stats.dat$num.partitions, stats.dat$rf_dist, xlab="'best' number of partitions", ylab="Robinson Foulds distance", col="white")
points(AIC.dat$num.partitions+0.1, AIC.dat$rf_dist, col="red", pch=AIC.dat$plot.pch, lwd=size)# lwd=AIC.dat$samp.size/sizer, 
points(AICc.dat$num.partitions, AICc.dat$rf_dist, col="darkred", pch=AICc.dat$plot.pch, lwd=size)
points(BIC.dat$num.partitions-0.1, BIC.dat$rf_dist, col="darkblue", pch=BIC.dat$plot.pch, lwd=size)
points(jitter(noPart.dat$num.partitions, factor=3.5), noPart.dat$rf_dist, col="green3", pch=noPart.dat$plot.pch, lwd=size)
legend("topright", c("none", "BIC", "AIC", "AICc"), title="Partitioning criterion", col=c("green3", "darkblue", "red", "darkred"), pch=15, cex=0.95, bty="n", ncol=2, pt.lwd=3)
legend(1,2.25e6, c("10", "20", "40", "80"), pch=c(0,1,2,6), title="Num. genes sampled", bty="n", ncol=2, pt.lwd=2, cex=0.95)

plot(stats.dat$num.partitions, stats.dat$eucl_dist/stats.dat$mut.rate, xlab="best number of partitions", ylab="scaled Euclidean distance", col="white")
points(AIC.dat$num.partitions+0.1, AIC.dat$eucl_dist/AIC.dat$mut.rate, col="red", pch=AIC.dat$plot.pch, lwd=size)# lwd=AIC.dat$samp.size/sizer, 
points(AICc.dat$num.partitions, AICc.dat$eucl_dist/AICc.dat$mut.rate, col="darkred", pch=AICc.dat$plot.pch, lwd=size)
points(BIC.dat$num.partitions-0.1, BIC.dat$eucl_dist/BIC.dat$mut.rate, col="darkblue", pch=BIC.dat$plot.pch, lwd=size)
points(jitter(noPart.dat$num.partitions, factor=3.5), noPart.dat$eucl_dist/noPart.dat$mut.rate, col="green3", pch=noPart.dat$plot.pch, lwd=size)
dev.off()

















plot(0, type="n", xlim=c(0,14), xaxt="n", ylab="partitioning raw log likelihood", ylim=range(stats.dat$part.log.lik, na.rm=TRUE))
x.it <- 1
for(i in uniq.names){
	temp.dat <- stats.dat[stats.dat$uniqueID == i ,]
	if(unique(temp.dat$species) == "lopho") color <- "red"
	if(unique(temp.dat$species) == "myria") color <- "darkblue"
	if(unique(temp.dat$search.algo) == "rcluster") pt.pch <- 0
	if(unique(temp.dat$search.algo) == "greedy") pt.pch <- 1
	if(unique(temp.dat$search.algo) == "noPart") pt.pch <- 2
	
	points(rep(x.it, dim(temp.dat)[1]), temp.dat$part.log.lik, col=color, pch=pt.pch)
	x.it <- x.it + 1
}


plot(stats.dat$num.partition.steps, stats.dat$part.log.lik, xlab="number steps to partition", ylab="log likelihood partition model", col="white")
# sub1 <- which(stats.dat$species == "lopho" & stats.dat$search.algo == "rcluster")
# sub2 <- which(stats.dat$species == "lopho" & stats.dat$search.algo == "greedy")
# sub3 <- which(stats.dat$species == "myria" & stats.dat$search.algo == "rcluster")
# sub4 <- which(stats.dat$species == "myria" & stats.dat$search.algo == "greedy")
# points(stats.dat$num.partition.steps[sub1], stats.dat$part.log.lik[sub1], col="red", pch=1)
# points(stats.dat$num.partition.steps[sub2], stats.dat$part.log.lik[sub2], col="red", pch=2)
# points(stats.dat$num.partition.steps[sub3], stats.dat$part.log.lik[sub3], col="darkblue", pch=1)
# points(stats.dat$num.partition.steps[sub4], stats.dat$part.log.lik[sub4], col="darkblue", pch=2)

points(AIC.dat$num.partition.steps, AIC.dat$part.log.lik, col="red")
points(BIC.dat$num.partition.steps, BIC.dat$part.log.lik, col="blue")
points(AICc.dat$num.partition.steps, AICc.dat$part.log.lik, col="darkred")





















plot(stats.dat$num.partitions, stats.dat$BIC, xlab="number of partitions", ylab="treebuilding BIC", col="white")
points(AIC.dat$num.partitions, AIC.dat$BIC, col="red", lwd=2)
points(AICc.dat$num.partitions, AICc.dat$BIC, col="darkred", lwd=2)
points(BIC.dat$num.partitions, BIC.dat$BIC, col="darkblue", lwd=2)
points(noPart.dat$num.partitions, noPart.dat$BIC, col="green3", lwd=2)


sd=0.05
plot(0, type="n", xlim=c(0.5,4.5), xaxt="n", xlab="model", ylab="treebuilding AIC", ylim=range(as.numeric(stats.dat$AIC), na.rm=TRUE))
points(rnorm(n=dim(noPart.dat)[1], mean=1, sd=sd), noPart.dat$AIC)
points(rnorm(n=dim(AIC.dat)[1], mean=2, sd=sd), AIC.dat$AIC)
points(rnorm(n=dim(AICc.dat)[1], mean=3, sd=sd), AICc.dat$AIC)
points(rnorm(n=dim(BIC.dat)[1], mean=4, sd=sd), BIC.dat$AIC)

plot(0, type="n", xlim=c(0.5,4.5), xaxt="n", xlab="model", ylab="treebuilding AICc", ylim=range(as.numeric(stats.dat$AICc), na.rm=TRUE))
points(rnorm(n=dim(noPart.dat)[1], mean=1, sd=sd), noPart.dat$AICc)
points(rnorm(n=dim(AIC.dat)[1], mean=2, sd=sd), AIC.dat$AICc)
points(rnorm(n=dim(AICc.dat)[1], mean=3, sd=sd), AICc.dat$AICc)
points(rnorm(n=dim(BIC.dat)[1], mean=4, sd=sd), BIC.dat$AICc)

plot(0, type="n", xlim=c(0.5,4.5), xaxt="n", xlab="model", ylab="treebuilding BIC", ylim=range(as.numeric(stats.dat$BIC), na.rm=TRUE))
points(rnorm(n=dim(noPart.dat)[1], mean=1, sd=sd), noPart.dat$BIC)
points(rnorm(n=dim(AIC.dat)[1], mean=2, sd=sd), AIC.dat$BIC)
points(rnorm(n=dim(AICc.dat)[1], mean=3, sd=sd), AICc.dat$BIC)
points(rnorm(n=dim(BIC.dat)[1], mean=4, sd=sd), BIC.dat$BIC)







plot(0, type="n", ylab="RF distance", xlim=c(0.5,4.5), ylim=c(0,0.25), xaxt="n", xlab="sample size (genes)")
axis(side=1, at=1:4, c("10", "20", "40", "80"))
points(rep(1, 11), dat2$rf_dist[c(1:4,11:17)]) # samp size 10
points(rep(2, 8), dat2$rf_dist[c(5:8,18:21)]) # samp size 20
points(rep(3, 5), dat2$rf_dist[c(9,22:25)]) # samp size 40
points(rep(4, 2), dat2$rf_dist[c(10,26)]) # samp size 80

points(rep(1.1, 7), dat4$rf_dist[c(1:4,10:12)], col="red") # samp size 10
points(rep(2.1, 6), dat4$rf_dist[c(5:7,13:15)], col="red") # samp size 20
points(rep(3.1, 2), dat4$rf_dist[c(8,16)], col="red") # samp size 40
points(rep(4.1, 2), dat4$rf_dist[c(9,17)], col="red") # samp size 80

points(rep(0.9, 8), dat4$rf_dist[c(4,8:14)], col="blue") # samp size 10
points(rep(1.9, 5), dat4$rf_dist[c(5,15:18)], col="blue") # samp size 20
points(rep(2.9, 2), dat4$rf_dist[c(6,19)], col="blue") # samp size 40
points(rep(3.9, 2), dat4$rf_dist[c(7,20)], col="blue") # samp size 80



plot(0, type="n", ylab="Euclidean distance", xlim=c(0.5,4.5), ylim=c(400,1000))
points(rep(1, 11), dat2$eucl_dist[c(1:4,11:17)]) # samp size 10
points(rep(2, 8), dat2$eucl_dist[c(5:8,18:21)]) # samp size 20
points(rep(3, 5), dat2$eucl_dist[c(9,22:25)]) # samp size 40
points(rep(4, 2), dat2$eucl_dist[c(10,26)]) # samp size 80

points(rep(1.1, 7), dat4$eucl_dist[c(1:4,10:12)], col="red") # samp size 10
points(rep(2.1, 6), dat4$eucl_dist[c(5:7,13:15)], col="red") # samp size 20
points(rep(3.1, 2), dat4$eucl_dist[c(8,16)], col="red") # samp size 40
points(rep(4.1, 2), dat4$eucl_dist[c(9,17)], col="red") # samp size 80

points(rep(0.9, 8), dat4$eucl_dist[c(4,8:14)], col="blue") # samp size 10
points(rep(1.9, 5), dat4$eucl_dist[c(5,15:18)], col="blue") # samp size 20
points(rep(2.9, 2), dat4$eucl_dist[c(6,19)], col="blue") # samp size 40
points(rep(3.9, 2), dat4$eucl_dist[c(7,20)], col="blue") # samp size 80
