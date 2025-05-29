## ## power function for binary outcome with normal approx
## power_binary_normalapprox <- function(active_p, plb_p,
##                                       active_ss, plb_ss,
##                                       alpha) {
##     se <- sqrt(active_p*(1-active_p)/active_ss + plb_p*(1-plb_p)/plb_ss)
##     z_crt <- -qnorm(alpha/2)
##     delta <- active_p - plb_p
##     1 - pnorm(z_crt-delta/se) + pnorm(-z_crt-delta/se)
## }


## ss_binary_normalapprox <- function(active_p, plb_p,
##                                    rand_ratio, power,
##                                    alpha) {
##     power_fnct <- function(sample_size) {
##         active_ss <- sample_size/(rand_ratio + 1)*rand_ratio
##         plb_ss <- sample_size - active_ss
##         cart <- power_binary_normalapprox(active_p, plb_p, active_ss, plb_ss, alpha)
##         cart - power
##     }
##     uniroot(power_fnct, lower = 5, upper = 5000)$root
## }



##' Sample size calculation: binary endpoint with Fisher test
##'
##' @description
##' Sample size calculation for binary endpoint with Fisher test.
##'
##' `getPower_binary_fisher()` provides power with given sample size and effect size.
##' 
##' `getSampleSize_Zstat()` provide total sample size needed to reach specified power and effect size.
##'
##' @param total_sample_size total sample size of two arms
##' @param ss_ratio randomization ratio between two arms
##' @param mean_AC treatment effect of active arm
##' @param mean_PLB PLB effect
##' @param drop_rate early termination rate
##' @param conf.level confidence level
##' @param power power
##'
##' @return see description part
##'
##' @examples
##' ## get power for a 240 pts study with 10% drop out rate,
##' ## randomization ratio 1:1, 50% and 30% response rate for
##' ## two arms separately, at 0.05 significance level
##'
##' getPower_binary_fisher(240, 1, 0.5, 0.3, 0.1, 0.95)
##'
##' ## search sample size needed to reach 90% power with same setup
##'
##' # comment out due to time-consuming code
##' #getSampleSize_binary_fisher(0.9, 1, 0.5, 0.3, 0.1, 0.95)
##'
##' @export
getPower_binary_fisher <- function(total_sample_size, ss_ratio = 1,
                                   mean_AC, mean_PLB, drop_rate = 0.1,
                                   conf.level = 0.95){
    total_sample_size <- total_sample_size * (1-drop_rate)
    sample_size_AC <- total_sample_size*(ss_ratio/(ss_ratio+1))
    sample_size_PLB <- total_sample_size - sample_size_AC

    active_p <- mean_AC
    plb_p <- mean_PLB
    active_ss <- round(sample_size_AC)
    plb_ss <- round(sample_size_PLB)

    do1 <- function(){
        cart_a <- rbinom(1, active_ss, active_p)
        cart_b <- active_ss - cart_a
        cart_c <- rbinom(1, plb_ss, plb_p)
        cart_d <- plb_ss - cart_c
        cart <- matrix(c(cart_a, cart_b, cart_c, cart_d), nrow = 2)
        fisher.test(cart)$p.value
    }

    result <- unlist(replicate(5000, do1()))
    sum(result <= (1 - conf.level)) / 5000
    
}

##' @rdname getPower_binary_fisher
##' @export
getSampleSize_binary_fisher <- function(power, ss_ratio = 1,
                                        mean_AC, mean_PLB,
                                        drop_rate = 0.1,
                                        conf.level = 0.95){
    f <- function(total_sample_size) {
        getPower_binary_fisher(total_sample_size, ss_ratio,
                               mean_AC, mean_PLB,
                               drop_rate,
                               conf.level) - power
    }

    if (getPower_binary_fisher(30,  ss_ratio,
                               mean_AC, mean_PLB,
                               drop_rate,
                               conf.level) >= power) return("<=30 sample size needed")
    
    round(uniroot(f, c(30, 5000))$root)
}
