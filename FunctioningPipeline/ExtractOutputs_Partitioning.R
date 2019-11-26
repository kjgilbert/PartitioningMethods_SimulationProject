
args<-commandArgs(trailingOnly=TRUE)

# args[1]: -- input directory
setwd(args[1])



dat <- "partitioning_treebuilding_output.txt"

species <- matrix(unlist(strsplit(system(paste(c("grep -n sampledMSAs ", dat), collapse=""), intern=TRUE), split="_")), ncol=3, byrow=TRUE)[,2]
# each output file should have 7 tree inferences, so repeat each species 7 times
species.col <- rep(species, each=7)

search.algo <- trimws(matrix(unlist(strsplit(system(paste(c("grep -n search ", dat), collapse=""), intern=TRUE), split=":")), ncol=3, byrow=TRUE)[,3])
scheme.criteria <- matrix(unlist(strsplit(system(paste(c("grep -n Scheme ", dat), collapse=""), intern=TRUE), split="Scheme")), ncol=6, byrow=TRUE)[,3:4]
part.log.lik <- as.numeric(matrix(unlist(strsplit(scheme.criteria[,1], split=": ")), ncol=2, byrow=TRUE)[,2])
part.IC <- trimws(matrix(unlist(strsplit(scheme.criteria[,2], split=": ")), ncol=2, byrow=TRUE)[,1])
part.IC.score <- as.numeric(matrix(unlist(strsplit(scheme.criteria[,2], split=": ")), ncol=2, byrow=TRUE)[,2])
partitions <- matrix( unlist(strsplit(system(paste(c("grep -n Scheme ", dat), collapse=""), intern=TRUE), split="= ")), ncol=2, byrow=TRUE)[,2]

max.lik.lines <- as.numeric(matrix(unlist(strsplit(system(paste(c("grep -n MAXIMUM ", dat), collapse=""), intern=TRUE), split=":")), ncol=10, byrow=TRUE)[,1])+1
max.lik.stats <- matrix(unlist(strsplit(system(paste(c("grep -n MAXIMUM ", dat), collapse=""), intern=TRUE), split=" ")), ncol=62, byrow=TRUE)[,c(9,16,25,31,38,44,52,58)]

cons.lines <- as.numeric(matrix(unlist(strsplit(system(paste(c("grep -n CONSENSUS ", dat), collapse=""), intern=TRUE), split=":")), ncol=4, byrow=TRUE)[,1])+1
cons.stats <- matrix(unlist(strsplit(system(paste(c("grep -n CONSENSUS ", dat), collapse=""), intern=TRUE), split=" ")), ncol=24, byrow=TRUE)[,c(15,24)]


# add on the no partitioning output:
search.algo <- c(search.algo, "none")
part.log.lik <- c(part.log.lik, NA)
part.IC <- c(part.IC, NA)
part.IC.score <- c(part.IC.score, NA)
partitions <- c(partitions, NA)


# should be seven lines for each species, so can iterate through
max.lik.trees <- rep(NA, length(max.lik.lines))
it <- 1
for(i in max.lik.lines){
	con <- file(dat, open = "r")
	line <- readLines(con, i)
	close(con)
	tree.line <- line[i]
	max.lik.trees[it] <- tree.line
	write(tree.line, file=paste(c("AnalysisResults/maxtree_", species[it], "_", search.algo[it], part.IC[it], ".nwk"), collapse=""))
	it <- it + 1
}


cons.trees <- rep(NA, length(cons.lines))
it <- 1
for(i in cons.lines){
	con <- file(dat, open = "r")
	line <- readLines(con, i)
	close(con)
	tree.line <- line[i]
	cons.trees[it] <- tree.line
	write(tree.line, file=paste(c("AnalysisResults/constree_", species[it], "_", search.algo[it], part.IC[it], ".nwk"), collapse=""))
	it <- it + 1
}



out.dat <- data.frame(cbind(species.col, search.algo, part.log.lik, part.IC, part.IC.score, partitions, max.lik.stats, cons.stats, max.lik.trees, cons.trees))
names(out.dat) <- c("species","search.algo", "part.log.lik", "part.IC", "part.IC.score", "partitions", "max.log.lik", "uncon.log.lik", "num.free.params", "AIC", "AICc", "BIC", "sum.branch.lens", "sum.int.branch.lens", "cons.log.lik", "cons.max.RF.dist", "max.lik.tree", "consen.tree")


write.csv(out.dat, file="AnalysisResults/AllStatsOutput.csv")
