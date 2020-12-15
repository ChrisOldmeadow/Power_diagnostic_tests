library(boot)
library(epitools) # for expanding a table
library(pbapply) # for progress bar



sim_confusion <- function(n, pd, se, sp) {
    # simulate data from the pre period
    # n = sample size
    # pd = prob of disease
    # se = sensitivty
    # sp = specificity

    nd <- rbinom(1, n, pd)
    tp <- rbinom(1, nd, se)
    tn <- rbinom(1, n-nd, sp)
    matrix(c(tp, nd-tp, n-nd-tn, tn), 2, 2)
}


lr_pos <- function(df,i){
    # relative comparison of two LR positives
    # df: data frame with the the Test result and Truth for two tests in wide
    #     format
    # i: indicies for selecting subsamples in the boot function
    tab.pre <- prop.table(table(df$TestPre[i],df$TruthPre[i]),margin = 2)
    sens.pre <- tab.pre[1,1]
    spec.pre <- tab.pre[2,2]
    lr.pre <- sens.pre/(1-spec.pre)
    tab.post <- prop.table(table(df$TestPost[i],df$TruthPost[i]),margin = 2)
    sens.post <- tab.post[1,1]
    spec.post <- tab.post[2,2]
    lr.post <- sens.post/(1-spec.post)
    lr <- lr.post/lr.pre
    return(lr)
}



lr_test_power <- function(){
    # simulate the data, bootstrap the CI for the relative difference
    # output is the lower bound of the CI
    simpre <- sim_confusion(n=200, pd=0.55, se=0.74, sp=0.65)
    simpost <- sim_confusion(n=200, pd=0.55, se=0.85, sp=0.80)

    dimnames(simpre) <- list(c("No","Yes"), c("No","Yes"))
    names(dimnames(simpre)) <- c("TestPre", "TruthPre")
    simdatpre <- expand.table(simpre)

    dimnames(simpost) <- list(c("No","Yes"), c("No","Yes"))
    names(dimnames(simpost)) <- c("TestPost", "TruthPost")
    simdatpost <- expand.table(simpost)

    simdat <- cbind(simdatpre,simdatpost)

    tmp <- boot(simdat,lr_pos,R=1000)
    res <- boot.ci(tmp,type = "norm")$norm[2] > 1
    return(res)
}


set.seed(123)
pow <- pbreplicate(1000,lr_test_power())
mean(pow)
