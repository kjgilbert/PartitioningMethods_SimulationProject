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
				"AllStatsOutput_Newest_Rep20_100-1000.csv"
				)
infiles.b <- c("RF_Eucl_results_New_Rep1_100-1000.txt",
				"RF_Eucl_results_New_Rep2_100-1000.txt",
				"RF_Eucl_results_New_Rep3_100-1000.txt",
				"RF_Eucl_results_New_Rep4_100-1000.txt",
				"RF_Eucl_results_New_Rep5_100-1000.txt",
				"RF_Eucl_results_New_Rep6_100-1000.txt",
				"RF_Eucl_results_New_Rep7_100-1000.txt",
				"RF_Eucl_results_New_Rep8_100-1000.txt",
				"RF_Eucl_results_New_Rep9_100-1000.txt",
				"RF_Eucl_results_New_Rep10_100-1000.txt",
				"RF_Eucl_results_New_Rep11_100-1000.txt",
				"RF_Eucl_results_New_Rep12_100-1000.txt",
				"RF_Eucl_results_New_Rep13_100-1000.txt",
				"RF_Eucl_results_New_Rep14_100-1000.txt",
				"RF_Eucl_results_New_Rep15_100-1000.txt",
				"RF_Eucl_results_New_Rep16_100-1000.txt",
				"RF_Eucl_results_New_Rep17_100-1000.txt",
				"RF_Eucl_results_New_Rep18_100-1000.txt",
				"RF_Eucl_results_New_Rep19_100-1000.txt",
				"RF_Eucl_results_New_Rep20_100-1000.txt"
				)

dat <- NULL
for(i in 1:15){#length(infiles.a)){ #**************UNCOMMENT ME FOR ALL REPS
	temp.dat1 <- read.csv(infiles.a[i], header=TRUE, stringsAsFactors=FALSE)
	temp.dat1$uniqueID <- paste(temp.dat1$species, temp.dat1$samp.size, temp.dat1$search.algo, temp.dat1$part.IC, sep="_")
	temp.dat2 <- read.table(infiles.b[i], header=TRUE, stringsAsFactors=FALSE)
	temp.dat2$uniqueID <- make.merge.id(temp.dat2)
	temp.dat <- merge(temp.dat1, temp.dat2, by="uniqueID")
	temp.dat$gene.size <- 100
	temp.dat$mut.rate <- 1000
	temp.dat$replicate <- i
	
	dat <- rbind(dat, temp.dat)
}

no.part.dat <- dat[dat$part.IC == "noPart" ,]


## IS there a statistically significant difference between the RF/eucl distance to true tree of no partitioning vs with partitioning??
# compare at each sample size, for each species, for both RF and euclidean distance

stats <- c("rf_dist", "eucl_dist")
spp <- c("myria", "lopho")
samps <- c(10,20,40,80)
stat.data.set <- data.frame(matrix(NA, ncol=8, nrow=8))
names(stat.data.set) <- c("species", "samp.size", "rf_twoside_p", "rf_less_p", "rf_greater_p", "eucl_twoside_p", "eucl_less_p", "eucl_greater_p")

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
		
		stat.data.set[iterate, ] <- c(j, k, rf_ts_p, rf_less_p, rf_greater_p, eucl_ts_p, eucl_less_p, eucl_greater_p)
		iterate <- iterate + 1
	}
}



####################################
samp.size.col <- c("steelblue4", "darkorange2", "green4", "red3")
outline.col <- "gray30"
samp.sizes <- c(10,20,40,80)
####################################
ylim1 <- c(-0.13, 0.13)
ylim2 <- c(-1.3, 1.3)
box.width <- 1
line.type <- 1
line.width <- 1.75
ytext1 <- 0.11
ytext2 <- 1.1



pdf("Results_Plots/Boxplots_withStats.pdf", width=6, height=8)
par(mfrow=c(2,2), mar=c(3.25,3.25,1.5,0.25))


part.crit <- "BIC"


plot(0, type="n", xlim=c(0.5,4.5), ylim=ylim1, xlab="", xaxt="n", ylab="diff in RF dist from true tree for 'no part' - 'part'", main="Lophotrochozoa", mgp=c(2.1,0.9,0))
axis(side=1, at=c(1:4), c(1:4))
abline(h=0, col="gray75", lty=3)
text(2.5,ytext1, expression(italic("partitioning better")), cex=1, col="green4")
text(2.5,-ytext1, expression(italic("no partitioning better")), cex=1, col="red3")
for(i in 1:4){
	plot.dat <- get.boxplot.dat(dat=dat, samp.size=samp.sizes[i], spp="lopho", partitioning.criterion=part.crit, distance="rf_dist")
	if(is.na(unique(plot.dat)) & length(unique(plot.dat)) == 1) next
	boxplot(plot.dat, add=TRUE, at=i, pch=23, col=samp.size.col[i], border=outline.col, width=box.width, lty=line.type, lwd=line.width, axes=FALSE)
}
# add p values from 2-sided wilcoxon test
text(x=c(0.65,1:4), y=rep(0.125, 5), cex=0.75, c(expression(italic("p")*" = "), 
	format(as.numeric(stat.data.set$rf_twoside_p[stat.data.set$species == "lopho" & stat.data.set$samp.size == 10][1]), digits=2), 
	format(as.numeric(stat.data.set$rf_twoside_p[stat.data.set$species == "lopho" & stat.data.set$samp.size == 20][1]), digits=2),
	format(as.numeric(stat.data.set$rf_twoside_p[stat.data.set$species == "lopho" & stat.data.set$samp.size == 40][1]), digits=2),
	format(as.numeric(stat.data.set$rf_twoside_p[stat.data.set$species == "lopho" & stat.data.set$samp.size == 80][1]), digits=2)))


plot(0, type="n", xlim=c(0.5,4.5), ylim=ylim1, xlab="", xaxt="n", ylab="", main="Myriapoda", mgp=c(2.1,0.9,0))
axis(side=1, at=c(1:4), c(1:4))
abline(h=0, col="gray75", lty=3)
text(2.5,ytext1, expression(italic("partitioning better")), cex=1, col="green4")
text(2.5,-ytext1, expression(italic("no partitioning better")), cex=1, col="red3")
for(i in 1:4){
	plot.dat <- get.boxplot.dat(dat=dat, samp.size=samp.sizes[i], spp="myria", partitioning.criterion=part.crit, distance="rf_dist")
	if(is.na(unique(plot.dat)) & length(unique(plot.dat)) == 1) next
	boxplot(plot.dat, add=TRUE, at=i, pch=23, col=samp.size.col[i], border=outline.col, width=box.width, lty=line.type, lwd=line.width, axes=FALSE)
}
text(x=c(0.65,1:4), y=rep(0.125, 5), cex=0.75, c(expression(italic("p")*" = "), 
	format(as.numeric(stat.data.set$rf_twoside_p[stat.data.set$species == "myria" & stat.data.set$samp.size == 10][1]), digits=2), 
	format(as.numeric(stat.data.set$rf_twoside_p[stat.data.set$species == "myria" & stat.data.set$samp.size == 20][1]), digits=2),
	format(as.numeric(stat.data.set$rf_twoside_p[stat.data.set$species == "myria" & stat.data.set$samp.size == 40][1]), digits=2),
	format(as.numeric(stat.data.set$rf_twoside_p[stat.data.set$species == "myria" & stat.data.set$samp.size == 80][1]), digits=2)))
	
	
plot(0, type="n", xlim=c(0.5,4.5), ylim=ylim2, xlab="num. genes", xaxt="n", ylab="diff in eucl dist from true tree for 'no part' - 'part'", main="", mgp=c(2.1,0.9,0))
axis(side=1, at=c(1:4), c(1:4))
abline(h=0, col="gray75", lty=3)
text(2.5,ytext2, expression(italic("partitioning better")), cex=1, col="green4")
text(2.5,-ytext2, expression(italic("no partitioning better")), cex=1, col="red3")
for(i in 1:4){
	plot.dat <- get.boxplot.dat(dat=dat, samp.size=samp.sizes[i], spp="lopho", partitioning.criterion=part.crit, distance="eucl_dist")
	if(is.na(unique(plot.dat)) & length(unique(plot.dat)) == 1) next
	boxplot(plot.dat, add=TRUE, at=i, pch=23, col=samp.size.col[i], border=outline.col, width=box.width, lty=line.type, lwd=line.width, axes=FALSE)
}
text(x=c(0.65,1:4), y=rep(1.25, 5), cex=0.75, c(expression(italic("p")*" = "), 
	format(as.numeric(stat.data.set$eucl_twoside_p[stat.data.set$species == "lopho" & stat.data.set$samp.size == 10][1]), digits=2), 
	format(as.numeric(stat.data.set$eucl_twoside_p[stat.data.set$species == "lopho" & stat.data.set$samp.size == 20][1]), digits=2),
	format(as.numeric(stat.data.set$eucl_twoside_p[stat.data.set$species == "lopho" & stat.data.set$samp.size == 40][1]), digits=2),
	format(as.numeric(stat.data.set$eucl_twoside_p[stat.data.set$species == "lopho" & stat.data.set$samp.size == 80][1]), digits=2)))


plot(0, type="n", xlim=c(0.5,4.5), ylim=ylim2, xlab="num. genes", xaxt="n", ylab="", main="", mgp=c(2.1,0.9,0))
axis(side=1, at=c(1:4), c(1:4))
abline(h=0, col="gray75", lty=3)
text(2.5,ytext2, expression(italic("partitioning better")), cex=1, col="green4")
text(2.5,-ytext2, expression(italic("no partitioning better")), cex=1, col="red3")
for(i in 1:4){
	plot.dat <- get.boxplot.dat(dat=dat, samp.size=samp.sizes[i], spp="myria", partitioning.criterion=part.crit, distance="eucl_dist")
	if(is.na(unique(plot.dat)) & length(unique(plot.dat)) == 1) next
	boxplot(plot.dat, add=TRUE, at=i, pch=23, col=samp.size.col[i], border=outline.col, width=box.width, lty=line.type, lwd=line.width, axes=FALSE)
}
text(x=c(0.65,1:4), y=rep(1.25, 5), cex=0.75, c(expression(italic("p")*" = "), 
	format(as.numeric(stat.data.set$eucl_twoside_p[stat.data.set$species == "myria" & stat.data.set$samp.size == 10][1]), digits=2), 
	format(as.numeric(stat.data.set$eucl_twoside_p[stat.data.set$species == "myria" & stat.data.set$samp.size == 20][1]), digits=2),
	format(as.numeric(stat.data.set$eucl_twoside_p[stat.data.set$species == "myria" & stat.data.set$samp.size == 40][1]), digits=2),
	format(as.numeric(stat.data.set$eucl_twoside_p[stat.data.set$species == "myria" & stat.data.set$samp.size == 80][1]), digits=2)))


dev.off()




