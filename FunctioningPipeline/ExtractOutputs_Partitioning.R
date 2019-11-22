
setwd("~/Documents/UNIL/PartitioningMethods_SimulationProject/FunctioningPipeline")



dat <- "partitioning_treebuilding_output.txt"

search.algo <- matrix(unlist(strsplit(system(paste(c("grep -n search ", dat), collapse=""), intern=TRUE), split=":")), ncol=3, byrow=TRUE)[,3]
scheme.criteria <- matrix(unlist(strsplit(system(paste(c("grep -n Scheme ", dat), collapse=""), intern=TRUE), split="Scheme")), ncol=6, byrow=TRUE)[,3:4]
part.log.lik <- as.numeric(matrix(unlist(strsplit(scheme.criteria[,1], split=": ")), ncol=2, byrow=TRUE)[,2])
part.IC <- matrix(unlist(strsplit(scheme.criteria[,2], split=": ")), ncol=2, byrow=TRUE)[,1]
part.IC.score <- as.numeric(matrix(unlist(strsplit(scheme.criteria[,2], split=": ")), ncol=2, byrow=TRUE)[,2])
partitions <- matrix( unlist(strsplit(system(paste(c("grep -n Scheme ", dat), collapse=""), intern=TRUE), split="= ")), ncol=2, byrow=TRUE)[,2]

max.lik.lines <- as.numeric(matrix(unlist(strsplit(system(paste(c("grep -n MAXIMUM ", dat), collapse=""), intern=TRUE), split=":")), ncol=10, byrow=TRUE)[,1])+1
max.lik.stats <- matrix(unlist(strsplit(system(paste(c("grep -n MAXIMUM ", dat), collapse=""), intern=TRUE), split=" ")), ncol=62, byrow=TRUE)[,c(9,16,25,31,38,44,52,58)]
max.lik.tree <- dat[max.lik.lines,] ## somehow grep or echo or cat those lines!!

cons.lines <- as.numeric(matrix(unlist(strsplit(system(paste(c("grep -n CONSENSUS ", dat), collapse=""), intern=TRUE), split=":")), ncol=4, byrow=TRUE)[,1])+1
cons.stats <- matrix(unlist(strsplit(system(paste(c("grep -n CONSENSUS ", dat), collapse=""), intern=TRUE), split=" ")), ncol=24, byrow=TRUE)[,c(15,24)]
cons.tree <- 

out.dat <- data.frame(cbind(search.algo, part.log.lik, part.IC, part.IC.score, partitions, max.lik.stats, cons.stats))
names(out.dat) <- c("search.algo", "part.log.lik", "part.IC", "part.IC.score", "partitions", "max.log.lik", "uncon.log.lik", "num.free.params", "AIC", "AICc", "BIC", "sum.branch.lens", "sum.int.branch.lens", "cons.log.lik", "cons.max.RF.dist")

