Hi Kim, if you are interested, I’d like to involve you on a little subproject that could be 
part of a paper we would like to resubmit. The paper tests partitioning methods for phylogenetic 
inference. The idea of partitioning methods is that there are different subsets of sites, each 
evolving under different models of evolution, so they try hard to find partitions of the data 
that maximise some overall information criterion.

Like AIC or BIC. Perhaps you are familiar with these.

Anyway, the point is that we have empirical evidence that while these measures indeed 
increase, the trees don’t really get better

My conjecture is that this is kind of overfitting

When we repeatedly do statistical tests, it’s now common knowledge that we need some sort 
of multiple testing correction

Yet, as far as I can tell, for model testing, somehow people don’t seem to worry about 
overfitting and multiple testing corrections.

So to illustrate the problem, and convince the readers that this is potentially a problem, 
I think we could do a little simulation study

Basically generate synthetic datasets without partitioning, and then run these methods, 
which test a huge number of models, and show that they always end up selecting the wrong 
model if we try enough combinations.

That seems like a pretty straightforward conceptual argument but we can still learn from 
the simulation what sort of delta AIC we can expect etc.