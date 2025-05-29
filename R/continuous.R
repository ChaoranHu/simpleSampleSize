##' Sample size calculation: continuous endpoint and Z/t statistic
##'
##' @description
##' Sample size calculation for continous endpoint with Z/t statistic.
##'
##' `getPower_Zstat()` and `getPower_Tstat()` provides power with given sample size and effect size.
##' 
##' `getEffectSize_Zstat()` provides minimal effect size can be detected (defined as statistical significant) with given sample size.
##' 
##' `getSampleSize_Zstat()` and `getSampleSize_Tstat` provide total sample size needed to reach specified power and effect size.
##'
##'
##' @param total_sample_size total sample size of two arms
##' @param ss_ratio randomization ratio between two arms
##' @param mean_AC treatment effect of active arm
##' @param mean_PLB PLB effect
##' @param sd standard deviation
##' @param drop_rate early termination rate
##' @param conf.level confidence level
##' @param power power
##'
##' @return see description part
##'
##' @examples
##' ## using Z-test as an example
##' ## get power for a 60 pts study with 10% drop out rate,
##' ## randomization ratio 1:1, 5 and 3 mean for
##' ## two arms separately with common sd 3, at 0.05 significance level
##'
##' getPower_Zstat(60, 1, 5, 3, 3, 0.1, 0.95)
##'
##' ## search sample size needed to reach 90% power with same setup
##'
##' getSampleSize_Zstat(0.9, 1, 5, 3, 3, 0.1, 0.95)
##'
##' ## calculate the minimal effect size can be detected
##' ## (that is statistical significant) with 60 pts study
##' getEffectSize_Zstat(60, 1, 3, 0.1, 0.95)
##' 
##' @export
getPower_Zstat <- function(total_sample_size, ss_ratio = 1,
                           mean_AC, mean_PLB, sd, drop_rate = 0.1,
                           conf.level = 0.95) {
    total_sample_size <- total_sample_size * (1-drop_rate)
    sample_size_AC <- total_sample_size*(ss_ratio/(ss_ratio+1))
    sample_size_PLB <- total_sample_size - sample_size_AC
    
    z_critical <- qnorm((1 - conf.level)/2 + conf.level)
    cartA <- (mean_AC - mean_PLB) / sd / sqrt(1/sample_size_AC + 1/sample_size_PLB)
    1 - pnorm(z_critical - cartA) + pnorm(-z_critical - cartA)
}


##' @rdname getPower_Zstat
##' @export
getEffectSize_Zstat <- function(total_sample_size, ss_ratio = 1,
                                sd, drop_rate = 0.1,
                                conf.level = 0.95){
    total_sample_size <- total_sample_size * (1-drop_rate)
    sample_size_AC <- total_sample_size*(ss_ratio/(ss_ratio+1))
    sample_size_PLB <- total_sample_size - sample_size_AC
    
    z_critical <- -qnorm((1 - conf.level)/2)

    z_critical * sd * sqrt(1/sample_size_AC + 1/sample_size_PLB)
}


##' @rdname getPower_Zstat
##' @export
getSampleSize_Zstat <- function(power, ss_ratio = 1,
                                mean_AC, mean_PLB,
                                sd, drop_rate = 0.1,
                                conf.level = 0.95){
    f <- function(total_sample_size) {
        getPower_Zstat(total_sample_size, ss_ratio,
                       mean_AC, mean_PLB,
                       sd, drop_rate,
                       conf.level) - power
    }
    round(uniroot(f, c(10, 5000))$root)
}


##' @rdname getPower_Zstat
##' @export
getPower_Tstat <- function(total_sample_size, ss_ratio = 1,
                           mean_AC, mean_PLB, sd, drop_rate = 0.1,
                           conf.level = 0.95) {
    total_sample_size <- total_sample_size * (1-drop_rate)
    sample_size_AC <- total_sample_size*(ss_ratio/(ss_ratio+1))
    sample_size_PLB <- total_sample_size - sample_size_AC
    
 
    active_ss <- round(sample_size_AC)
    plb_ss <- round(sample_size_PLB)

    do1 <- function(){
        cart_a <- rnorm(active_ss, mean_AC, sd)
        cart_c <- rnorm(plb_ss, mean_PLB, sd)
        t.test(cart_a, cart_c, var.equal = TRUE)$p.value
    }

    result <- unlist(replicate(10000, do1()))
    sum(result <= (1 - conf.level)) / 10000
}

##' @rdname getPower_Zstat
##' @export
getSampleSize_Tstat <- function(power, ss_ratio = 1,
                                mean_AC, mean_PLB,
                                sd, drop_rate = 0.1,
                                conf.level = 0.95){
    f <- function(total_sample_size) {
        getPower_Tstat(total_sample_size, ss_ratio,
                       mean_AC, mean_PLB,
                       sd, drop_rate,
                       conf.level) - power
    }
    round(uniroot(f, c(3, 5000))$root)
}

