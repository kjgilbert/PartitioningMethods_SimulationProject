setwd("~/Documents/UNIL/PartitioningMethods_SimulationProject/FunctioningPipeline/AnalysisResults")

# dat <- readLines("prelimResults.txt")

dat <- read.csv("prelim_lopho_N10.csv", na.strings="-")

## LOWER IC IS BETTER!

plot(1:3, c(dat$AIC[3], dat$AICc[2], dat$BIC[1]), xlim=c(0.5,3.5), ylim=c(103400,104700), pch=16, xlab="Criterion", ylab="Score", xaxt="n", main="Partitioning")
axis(side=1, at=1:3, labels=c("AIC", "AICc", "BIC"))
points(1.1:3.1, c(dat$AIC[4], dat$AICc[6], dat$BIC[5]), pch=17, col="red")
legend("topleft", c("rcluster", "greedy"), col=c("black", "red"), pch=16)

text(0.6, 103525, "Number")
text(0.6, 103475, "Partitions")
text(x=c(1,1.1,2,2.1,3,3.1), y=rep(103500,6), labels=c(dat$numParts[3:4], dat$numParts[2], dat$numParts[6], dat$numParts[1], dat$numParts[5]))



plot(1:3, c(dat$AIC[3], dat$AICc[2], dat$BIC[1]), xlim=c(0.5,3.5), ylim=c(103000,105000), pch=16, xlab="Criterion", ylab="Score", xaxt="n", main="Treebuilding")
axis(side=1, at=1:3, labels=c("AIC", "AICc", "BIC"))
points(1.1:3.1, c(dat$AIC[4], dat$AICc[6], dat$BIC[5]), pch=17, col="red")
legend("topleft", c("rcluster", "greedy"), col=c("black", "red"), pch=16)


plot(dat$sumBranchLengths, dat$sumInternalBranchLengths)
plot(dat$iqTreeAIC, dat$iqTreeAICc, ylim=c(103800,104600))
points(dat$iqTreeAIC, dat$iqTreeBIC, col="red")





plot(1:3, c(dat$AIC[3], dat$AICc[2], dat$BIC[1]), xlim=c(0.5,4.5), ylim=c(103850,104600), pch=16, xlab="Partitioning Criterion", ylab="Score", xaxt="n", main="Partitioning & Treebuilding")
axis(side=1, at=1:4, labels=c("AIC", "AICc", "BIC", "none"))
points(1.1:3.1, c(dat$AIC[4], dat$AICc[6], dat$BIC[5]), pch=17, col="red")
points(c(1,1.1,2,2.1,3,3.1,4), c(dat$iqTreeAIC[3], dat$iqTreeAIC[4], dat$iqTreeAIC[2], dat$iqTreeAIC[6], dat$iqTreeAIC[1], dat$iqTreeAIC[5], dat$iqTreeAIC[7]), col="blue", lwd=2)
##points(c(1,1.1,2,2.1,3,3.1,4), c(dat$iqTreeAICc[3], dat$iqTreeAICc[4], dat$iqTreeAICc[2], dat$iqTreeAICc[6], dat$iqTreeAICc[1], dat$iqTreeAICc[5], dat$iqTreeAICc[7]), col="green3", lwd=2)
points(c(1,1.1,2,2.1,3,3.1,4), c(dat$iqTreeBIC[3], dat$iqTreeBIC[4], dat$iqTreeBIC[2], dat$iqTreeBIC[6], dat$iqTreeBIC[1], dat$iqTreeBIC[5], dat$iqTreeBIC[7]), col="orange", lwd=2)
legend("left", c("rcluster", "greedy"), col=c("black", "red"), pch=16)
legend("right", c("IQTree AIC", "IQTree BIC"), col=c("blue", "orange"), pch=1, pt.lwd=2)

text(x=c(1,1.1,2,2.1,3,3.1,4), y=rep(103870,6), labels=c(dat$numParts[3:4], dat$numParts[2], dat$numParts[6], dat$numParts[1], dat$numParts[5], 0))





plot(1:3, c(dat$AIC[3], dat$AICc[2], dat$BIC[1]), xlim=c(0.5,4.5), ylim=c(103920,103980), pch=16, xlab="Partitioning Criterion", ylab="Score", xaxt="n", main="Partitioning & Treebuilding")
axis(side=1, at=1:4, labels=c("AIC", "AICc", "BIC", "none"))
points(1.1:3.1, c(dat$AIC[4], dat$AICc[6], dat$BIC[5]), pch=17, col="red")
points(c(1,1.1,2,2.1,3,3.1,4), c(dat$iqTreeAIC[3], dat$iqTreeAIC[4], dat$iqTreeAIC[2], dat$iqTreeAIC[6], dat$iqTreeAIC[1], dat$iqTreeAIC[5], dat$iqTreeAIC[7]), col="blue", lwd=2)
##points(c(1,1.1,2,2.1,3,3.1,4), c(dat$iqTreeAICc[3], dat$iqTreeAICc[4], dat$iqTreeAICc[2], dat$iqTreeAICc[6], dat$iqTreeAICc[1], dat$iqTreeAICc[5], dat$iqTreeAICc[7]), col="green3", lwd=2)
points(c(1,1.1,2,2.1,3,3.1,4), c(dat$iqTreeBIC[3], dat$iqTreeBIC[4], dat$iqTreeBIC[2], dat$iqTreeBIC[6], dat$iqTreeBIC[1], dat$iqTreeBIC[5], dat$iqTreeBIC[7]), col="orange", lwd=2)
#legend("left", c("rcluster", "greedy"), col=c("black", "red"), pch=16)
legend("bottomright", c("IQTree AIC"), col=c("blue"), pch=1, pt.lwd=2)

text(x=c(1,1.1,2,2.1,3,3.1,4), y=rep(103940,6), labels=c(dat$numParts[3:4], dat$numParts[2], dat$numParts[6], dat$numParts[1], dat$numParts[5], 0))


plot(1:3, c(dat$AIC[3], dat$AICc[2], dat$BIC[1]), xlim=c(0.5,4.5), ylim=c(104535,104565), pch=16, xlab="Partitioning Criterion", ylab="Score", xaxt="n", main="Partitioning & Treebuilding")
axis(side=1, at=1:4, labels=c("AIC", "AICc", "BIC", "none"))
points(1.1:3.1, c(dat$AIC[4], dat$AICc[6], dat$BIC[5]), pch=17, col="red")
points(c(1,1.1,2,2.1,3,3.1,4), c(dat$iqTreeAIC[3], dat$iqTreeAIC[4], dat$iqTreeAIC[2], dat$iqTreeAIC[6], dat$iqTreeAIC[1], dat$iqTreeAIC[5], dat$iqTreeAIC[7]), col="blue", lwd=2)
##points(c(1,1.1,2,2.1,3,3.1,4), c(dat$iqTreeAICc[3], dat$iqTreeAICc[4], dat$iqTreeAICc[2], dat$iqTreeAICc[6], dat$iqTreeAICc[1], dat$iqTreeAICc[5], dat$iqTreeAICc[7]), col="green3", lwd=2)
points(c(1,1.1,2,2.1,3,3.1,4), c(dat$iqTreeBIC[3], dat$iqTreeBIC[4], dat$iqTreeBIC[2], dat$iqTreeBIC[6], dat$iqTreeBIC[1], dat$iqTreeBIC[5], dat$iqTreeBIC[7]), col="orange", lwd=2)
#legend("left", c("rcluster", "greedy"), col=c("black", "red"), pch=16)
legend("topright", c("IQTree BIC"), col=c("orange"), pch=1, pt.lwd=2)

text(x=c(1,1.1,2,2.1,3,3.1,4), y=rep(104557,6), labels=c(dat$numParts[3:4], dat$numParts[2], dat$numParts[6], dat$numParts[1], dat$numParts[5], 0))


plot(c(1:4, 1.1:3.1), rep(0.0847457627118644, 7), xaxt="n", ylim=c(0.08, 0.09), xlab="Partitioning Criterion", ylab="Normalized Robinson-Foulds", pch=15)
axis(side=1, at=1:4, labels=c("AIC", "AICc", "BIC", "none"))






plot(c(1:4, 1.1:3.1), c(498.1729185398443, 498.17346684526217, 498.1721147408412,   498.1721147408412, 498.16853197394056, 498.1586883819632     , 498.1758378443058), xaxt="n", ylim=c(498.14,498.2), xlab="Partitioning Criterion", ylab="Euclidean Distance", pch=15)
axis(side=1, at=1:4, labels=c("AIC", "AICc", "BIC", "none"))
