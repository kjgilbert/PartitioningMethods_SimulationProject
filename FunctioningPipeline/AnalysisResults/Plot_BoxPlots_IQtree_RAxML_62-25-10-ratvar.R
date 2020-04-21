setwd("~/Documents/UNIL/PartitioningMethods_SimulationProject/FunctioningPipeline/AnalysisResults")



##source('~/Documents/UNIL/PartitioningMethods_SimulationProject/FunctioningPipeline/AnalysisResults/ReadInData_Source.R', chdir = TRUE)











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
box.width <- 1
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

pdf(paste(c("Results_Plots/Boxplots_", part.crit, "_IQtree-RAxML.pdf"), collapse=""), width=12, height=10)
par(mar=c(2.5,4.25,1.75,0.25))
layout(matrix(c(1,1,3,3,1,1,3,3,1,1,3,3,2,2,4,4,2,2,4,4,2,2,4,4,5,5,5,5), ncol=4, byrow=TRUE))


#_________________________________________________________________________________________________________#

plot(0, type="n", xlim=c(0.5,18), ylim=ylim1, xlab="", xaxt="n", ylab="difference in distance from true tree", main="Robinson Fould distance - IQtree", mgp=c(3,0.9,0))
axis(side=1, at=c(2.5, 7, 11.5, 16), c("Lophotrochozoa (62 spp)", "Myriapoda (25 spp)", "10 spp", "10 spp w/ rate variation"))
mtext(side=2, "(no partitioning - partitioning)", cex=0.75, line=2)
abline(h=0, col="gray75", lty=3, lwd=1.5)
text(9,ytext1, expression(italic("partitioning better")), cex=1, col="green4")
text(9,-ytext1, expression(italic("no partitioning better")), cex=1, col="red3")

# lopho (62 spp)
for(i in 1:4){
	plot.dat <- get.boxplot.dat(dat=dat, samp.size=samp.sizes[i], spp="lopho", partitioning.criterion=part.crit, distance="rf_dist")
	boxplot(plot.dat, add=TRUE, at=i, pch=23, col=samp.size.col[i], border=outline.col, width=box.width, lty=line.type, lwd=line.width, axes=FALSE)
}
# add p values from 2-sided wilcoxon test
text(x=c(0.35,1:4), y=rep(p.height1, 5), cex=text.size, c(expression(italic("p")*" = "), 
	format(as.numeric(stat.data.set$rf_twoside_p[stat.data.set$species == "lopho" & stat.data.set$samp.size == 10][1]), digits=2), 
	format(as.numeric(stat.data.set$rf_twoside_p[stat.data.set$species == "lopho" & stat.data.set$samp.size == 20][1]), digits=2),
	format(as.numeric(stat.data.set$rf_twoside_p[stat.data.set$species == "lopho" & stat.data.set$samp.size == 40][1]), digits=2),
	format(as.numeric(stat.data.set$rf_twoside_p[stat.data.set$species == "lopho" & stat.data.set$samp.size == 80][1]), digits=2)))

# myria (25 spp)
for(i in 1:4){
	plot.dat <- get.boxplot.dat(dat=dat, samp.size=samp.sizes[i], spp="myria", partitioning.criterion=part.crit, distance="rf_dist")
	boxplot(plot.dat, add=TRUE, at=i+4.5, pch=23, col=samp.size.col[i], border=outline.col, width=box.width, lty=line.type, lwd=line.width, axes=FALSE)
}
text(x=1:4+4.5, y=rep(p.height1, 4), cex=text.size, c(
	format(as.numeric(stat.data.set$rf_twoside_p[stat.data.set$species == "myria" & stat.data.set$samp.size == 10][1]), digits=2), 
	format(as.numeric(stat.data.set$rf_twoside_p[stat.data.set$species == "myria" & stat.data.set$samp.size == 20][1]), digits=2),
	format(as.numeric(stat.data.set$rf_twoside_p[stat.data.set$species == "myria" & stat.data.set$samp.size == 40][1]), digits=2),
	format(as.numeric(stat.data.set$rf_twoside_p[stat.data.set$species == "myria" & stat.data.set$samp.size == 80][1]), digits=2)))

# 10 spp
for(i in 1:4){
	plot.dat <- get.boxplot.dat(dat=dat, samp.size=samp.sizes[i], spp="TenSpp", partitioning.criterion=part.crit, distance="rf_dist")
	boxplot(plot.dat, add=TRUE, at=i+9, pch=23, col=samp.size.col[i], border=outline.col, width=box.width, lty=line.type, lwd=line.width, axes=FALSE)
}
text(x=1:4+9, y=rep(p.height1, 4), cex=text.size, c(
	format(as.numeric(stat.data.set$rf_twoside_p[stat.data.set$species == "TenSpp" & stat.data.set$samp.size == 10][1]), digits=2), 
	format(as.numeric(stat.data.set$rf_twoside_p[stat.data.set$species == "TenSpp" & stat.data.set$samp.size == 20][1]), digits=2),
	format(as.numeric(stat.data.set$rf_twoside_p[stat.data.set$species == "TenSpp" & stat.data.set$samp.size == 40][1]), digits=2),
	format(as.numeric(stat.data.set$rf_twoside_p[stat.data.set$species == "TenSpp" & stat.data.set$samp.size == 80][1]), digits=2)))

# rate variation
for(i in 1:4){
	plot.dat <- get.boxplot.dat(dat=dat, samp.size=samp.sizes[i], spp="RateVar", partitioning.criterion=part.crit, distance="rf_dist")
	boxplot(plot.dat, add=TRUE, at=i++13.5, pch=23, col=samp.size.col[i], border=outline.col, width=box.width, lty=line.type, lwd=line.width, axes=FALSE)
}
text(x=1:4+13.5, y=rep(p.height1, 4), cex=text.size, c(
	format(as.numeric(stat.data.set$rf_twoside_p[stat.data.set$species == "RateVar" & stat.data.set$samp.size == 10][1]), digits=2), 
	format(as.numeric(stat.data.set$rf_twoside_p[stat.data.set$species == "RateVar" & stat.data.set$samp.size == 20][1]), digits=2),
	format(as.numeric(stat.data.set$rf_twoside_p[stat.data.set$species == "RateVar" & stat.data.set$samp.size == 40][1]), digits=2),
	format(as.numeric(stat.data.set$rf_twoside_p[stat.data.set$species == "RateVar" & stat.data.set$samp.size == 80][1]), digits=2)))

legend("bottomright", c("10", "20", "40", "80"), pch=15, col=samp.size.col, title = "num. genes", cex=0.875, pt.cex=1.5) # , bg="white"


#_________________________________________________________________________________________________________#


# eucl distance
plot(0, type="n", xlim=c(0.5,18), ylim=ylim2, xlab="", xaxt="n", ylab="difference in distance from true tree", main="Euclidean distance - IQtree", mgp=c(3,0.9,0))
mtext(side=2, "(no partitioning - partitioning)", cex=0.75, line=2)
axis(side=1, at=c(2.5, 7, 11.5, 16), c("Lophotrochozoa (62 spp)", "Myriapoda (25 spp)", "10 spp", "10 spp w/ rate variation"))
abline(h=0, col="gray75", lty=3, lwd=1.5)
text(9,ytext2, expression(italic("partitioning better")), cex=1, col="green4")
text(9,-ytext2, expression(italic("no partitioning better")), cex=1, col="red3")
for(i in 1:4){
	plot.dat <- get.boxplot.dat(dat=dat, samp.size=samp.sizes[i], spp="lopho", partitioning.criterion=part.crit, distance="eucl_dist")
	boxplot(plot.dat, add=TRUE, at=i, pch=23, col=samp.size.col[i], border=outline.col, width=box.width, lty=line.type, lwd=line.width, axes=FALSE)
}
text(x=c(0.35,1:4), y=rep(p.height2, 5), cex=text.size, c(expression(italic("p")*" = "), 
	format(as.numeric(stat.data.set$eucl_twoside_p[stat.data.set$species == "lopho" & stat.data.set$samp.size == 10][1]), digits=2), 
	format(as.numeric(stat.data.set$eucl_twoside_p[stat.data.set$species == "lopho" & stat.data.set$samp.size == 20][1]), digits=2),
	format(as.numeric(stat.data.set$eucl_twoside_p[stat.data.set$species == "lopho" & stat.data.set$samp.size == 40][1]), digits=2),
	format(as.numeric(stat.data.set$eucl_twoside_p[stat.data.set$species == "lopho" & stat.data.set$samp.size == 80][1]), digits=2)))

for(i in 1:4){
	plot.dat <- get.boxplot.dat(dat=dat, samp.size=samp.sizes[i], spp="myria", partitioning.criterion=part.crit, distance="eucl_dist")
	boxplot(plot.dat, add=TRUE, at=i+4.5, pch=23, col=samp.size.col[i], border=outline.col, width=box.width, lty=line.type, lwd=line.width, axes=FALSE)
}
text(x=1:4+4.5, y=rep(p.height2, 4), cex=text.size, c(
	format(as.numeric(stat.data.set$eucl_twoside_p[stat.data.set$species == "myria" & stat.data.set$samp.size == 10][1]), digits=2), 
	format(as.numeric(stat.data.set$eucl_twoside_p[stat.data.set$species == "myria" & stat.data.set$samp.size == 20][1]), digits=2),
	format(as.numeric(stat.data.set$eucl_twoside_p[stat.data.set$species == "myria" & stat.data.set$samp.size == 40][1]), digits=2),
	format(as.numeric(stat.data.set$eucl_twoside_p[stat.data.set$species == "myria" & stat.data.set$samp.size == 80][1]), digits=2)))

# 10 spp
for(i in 1:4){
	plot.dat <- get.boxplot.dat(dat=dat, samp.size=samp.sizes[i], spp="TenSpp", partitioning.criterion=part.crit, distance="eucl_dist")
	boxplot(plot.dat, add=TRUE, at=i+9, pch=23, col=samp.size.col[i], border=outline.col, width=box.width, lty=line.type, lwd=line.width, axes=FALSE)
}
text(x=1:4+9, y=rep(p.height2, 4), cex=text.size, c(
	format(as.numeric(stat.data.set$eucl_twoside_p[stat.data.set$species == "TenSpp" & stat.data.set$samp.size == 10][1]), digits=2), 
	format(as.numeric(stat.data.set$eucl_twoside_p[stat.data.set$species == "TenSpp" & stat.data.set$samp.size == 20][1]), digits=2),
	format(as.numeric(stat.data.set$eucl_twoside_p[stat.data.set$species == "TenSpp" & stat.data.set$samp.size == 40][1]), digits=2),
	format(as.numeric(stat.data.set$eucl_twoside_p[stat.data.set$species == "TenSpp" & stat.data.set$samp.size == 80][1]), digits=2)))

# rate variation
for(i in 1:4){
	plot.dat <- get.boxplot.dat(dat=dat, samp.size=samp.sizes[i], spp="RateVar", partitioning.criterion=part.crit, distance="eucl_dist")
	boxplot(plot.dat, add=TRUE, at=i++13.5, pch=23, col=samp.size.col[i], border=outline.col, width=box.width, lty=line.type, lwd=line.width, axes=FALSE)
}
text(x=1:4+13.5, y=rep(p.height2, 4), cex=text.size, c(
	format(as.numeric(stat.data.set$eucl_twoside_p[stat.data.set$species == "RateVar" & stat.data.set$samp.size == 10][1]), digits=2), 
	format(as.numeric(stat.data.set$eucl_twoside_p[stat.data.set$species == "RateVar" & stat.data.set$samp.size == 20][1]), digits=2),
	format(as.numeric(stat.data.set$eucl_twoside_p[stat.data.set$species == "RateVar" & stat.data.set$samp.size == 40][1]), digits=2),
	format(as.numeric(stat.data.set$eucl_twoside_p[stat.data.set$species == "RateVar" & stat.data.set$samp.size == 80][1]), digits=2)))


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
