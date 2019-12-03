
args<-commandArgs(trailingOnly=TRUE)

# args[1]: -- input directory
setwd(args[1])

system("mkdir AnalysisResults")

##setwd("~/Documents/UNIL/TEST_analyses")
outfiles <- system(paste("ls part_N*/partitioning_treebuilding_output.txt"), intern=TRUE)

num.files <- length(outfiles)
out.dat <- data.frame(matrix(NA, ncol=19, nrow=num.files))
names(out.dat) <- c("species","search.algo", "part.log.lik", "part.IC", "part.IC.score", "num.partition.steps", "partitions", "max.log.lik", "uncon.log.lik", "num.free.params", "AIC", "AICc", "BIC", "sum.branch.lens", "sum.int.branch.lens", "cons.log.lik", "cons.max.RF.dist", "max.lik.tree", "consen.tree")

iterate <- 1
for(i in outfiles){
	dat <- i
	
	output.exists <- system(paste(c("if grep -q MAXIMUM ", dat, "; then
    echo 1
else
    echo 0
fi"), collapse=""), intern=TRUE)

	if(output.exists == 0) next

	if(length(grep("noPart", dat)) == 0){		
		species <- matrix(unlist(strsplit(system(paste(c("head -n 1 ", dat), collapse=""), intern=TRUE), split="_")), ncol=3, byrow=TRUE)[,1]
		samp.size <- matrix(unlist(strsplit(system(paste(c("head -n 1 ", dat), collapse=""), intern=TRUE), split="_")), ncol=3, byrow=TRUE)[,2]

		search.algo <- trimws(matrix(unlist(strsplit(system(paste(c("grep -n search ", dat), collapse=""), intern=TRUE), split=":")), ncol=3, byrow=TRUE)[,3])
		scheme.criteria <- matrix(unlist(strsplit(system(paste(c("grep -n Scheme ", dat), collapse=""), intern=TRUE), split="Scheme")), ncol=6, byrow=TRUE)[,c(2:4,6)]
		part.log.lik <- as.numeric(matrix(unlist(strsplit(scheme.criteria[2], split=": ")), ncol=2, byrow=TRUE)[,2])
		part.IC <- trimws(matrix(unlist(strsplit(scheme.criteria[3], split=": ")), ncol=2, byrow=TRUE)[,1])
		part.IC.score <- as.numeric(matrix(unlist(strsplit(scheme.criteria[3], split=": ")), ncol=2, byrow=TRUE)[,2])
		partitioning.steps <- as.numeric(unlist(strsplit(matrix(unlist(strsplit(scheme.criteria[1], split=": ")), ncol=2, byrow=TRUE)[,2], split="_"))[2])
		partitions <- matrix(unlist(strsplit(scheme.criteria[4], split="= ")), ncol=2, byrow=TRUE)[,2]
			
		max.lik.lines <- as.numeric(matrix(unlist(strsplit(system(paste(c("grep -n MAXIMUM ", dat), collapse=""), intern=TRUE), split=":")), ncol=10, byrow=TRUE)[,1])+1
		max.lik.stats <- matrix(unlist(strsplit(system(paste(c("grep -n MAXIMUM ", dat), collapse=""), intern=TRUE), split=" ")), ncol=62, byrow=TRUE)[,c(9,16,25,31,38,44,52,58)]
		
		cons.lines <- as.numeric(matrix(unlist(strsplit(system(paste(c("grep -n CONSENSUS ", dat), collapse=""), intern=TRUE), split=":")), ncol=4, byrow=TRUE)[,1])+1
		cons.stats <- matrix(unlist(strsplit(system(paste(c("grep -n CONSENSUS ", dat), collapse=""), intern=TRUE), split=" ")), ncol=24, byrow=TRUE)[,c(15,24)]
		
		
		# and write the max likelihood tree to a file
		con <- file(dat, open = "r")
		line <- readLines(con, max.lik.lines)
		line2 <- readLines(con, cons.lines)
		close(con)
		max.tree <- line[max.lik.lines]
		cons.tree <- line2[2]
		write(max.tree, file=paste(c("AnalysisResults/maxtree_", species,"_", samp.size, "_", search.algo, part.IC, ".nwk"), collapse=""))
		write(cons.tree, file=paste(c("AnalysisResults/constree_", species,"_", samp.size, "_", search.algo, part.IC, ".nwk"), collapse=""))		

	}else{ # then the results are from IQtree without paritioning beforehand

		species <- unlist(strsplit(dat, split="_"))[3]
		samp.size <- unlist(strsplit(unlist(strsplit(dat, split="_"))[2], split="N"))[2]

		search.algo <- NA
		scheme.criteria <- NA
		part.log.lik <- NA
		part.IC <- NA
		part.IC.score <- NA
		partitioning.steps <- NA
		partitions <- NA
			
		max.lik.lines <- as.numeric(matrix(unlist(strsplit(system(paste(c("grep -n MAXIMUM ", dat), collapse=""), intern=TRUE), split=":")), ncol=10, byrow=TRUE)[,1])+1
		max.lik.stats <- matrix(unlist(strsplit(system(paste(c("grep -n MAXIMUM ", dat), collapse=""), intern=TRUE), split=" ")), ncol=62, byrow=TRUE)[,c(9,16,25,31,38,44,52,58)]
		
		cons.lines <- as.numeric(matrix(unlist(strsplit(system(paste(c("grep -n CONSENSUS ", dat), collapse=""), intern=TRUE), split=":")), ncol=4, byrow=TRUE)[,1])+1
		cons.stats <- matrix(unlist(strsplit(system(paste(c("grep -n CONSENSUS ", dat), collapse=""), intern=TRUE), split=" ")), ncol=24, byrow=TRUE)[,c(15,24)]
		
		
		# and write the max likelihood tree to a file
		con <- file(dat, open = "r")
		line <- readLines(con, max.lik.lines)
		line2 <- readLines(con, cons.lines)
		close(con)
		max.tree <- line[max.lik.lines]
		cons.tree <- line2[2]
		write(max.tree, file=paste(c("AnalysisResults/maxtree_", species,"_", samp.size, "_noPart.nwk"), collapse=""))
		write(cons.tree, file=paste(c("AnalysisResults/constree_", species,"_", samp.size, "_noPart.nwk"), collapse=""))		

	}
	
	out.dat[iterate ,] <- c(species, search.algo, part.log.lik, part.IC, part.IC.score, partitioning.steps, partitions, max.lik.stats, cons.stats, max.tree, cons.tree)	
	iterate <- iterate + 1
}







write.csv(out.dat, file=paste(c("AnalysisResults/AllStatsOutput_", tail(unlist(strsplit(getwd(), split="/")), n=1), ".csv"), collapse=""), quote=FALSE)
