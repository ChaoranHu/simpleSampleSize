
rqpois <- function(n, mu, theta) rnbinom(n=n, mu=mu, size=mu/(theta-1))


##' Sample size calculation: count endpoint with overdispersion
##'
##' @description
##' Sample size calculation for count endpoint with overdispersion.
##' In the simulation, overdispersion data is simulated with
##' negative binomial. Generated data is fitted with quasipoisson glm.
##'
##' `getPower_count_ovdisp()` provide power with given sample size
##' 
##' `getSampleSize_count_ovdisp` provide sample size with given power
##'
##' @param total_sample_size total sample size of two arms
##' @param ss_ratio randomization ratio between two arms
##' @param rate_AC event rate of active arm
##' @param rate_PLB event rate of PLB arm
##' @param over_disp_rate overdispersion rate, var/mean = 1+rate
##' @param drop_rate early termination rate
##' @param conf.level confidence level
##' @param n_sim number of simulations
##' @param power power
##' @param n_core number of cores used for parallel computation
##'
##' @return see description part
##'
##' @examples
##' ## comment out due to time consuming example
##' ## and parallel computation limitation in CRAN
##' 
##' ## get power for a 80 pts study with 10% drop out rate,
##' ## randomization ratio 1:1, 1 and 0.5 as the event rate for
##' ## two arms separately, 20% over-dispersive at 0.05 significance level
##'
##' #getPower_count_ovdisp(80, 1, 1, 0.5, 0.2, 0.1, 0.95)
##'
##' ## search sample size needed to reach 90% power with same setup
##' 
##' #getSampleSize_count_ovdisp(0.9, 1, 1, 0.5, 0.2, 0.1, 0.95)
##' 
##' @export
getPower_count_ovdisp <- function(total_sample_size, ss_ratio = 1,
                                  rate_AC, rate_PLB,
                                  over_disp_rate = 0.2, drop_rate = 0.1,
                                  conf.level = 0.95, n_sim = 10000,
                                  n_core = 6) {
    total_sample_size <- total_sample_size * (1-drop_rate)
    sample_size_AC <- total_sample_size*(ss_ratio/(ss_ratio+1))
    sample_size_PLB <- total_sample_size - sample_size_AC
    
    ## result <- numeric(n_sim)
    ## for (i in seq_len(n_sim)) {
    ##     dat_1 <- rqpois(sample_size_AC, rate_AC, over_disp_rate + 1)
    ##     dat_2 <- rqpois(sample_size_PLB, rate_PLB, over_disp_rate + 1)

    ##     dat <- rbind(cbind(dat_1, "AC"),
    ##                  cbind(dat_2, "PLB"))
    ##     dat <- as.data.frame(dat)
    ##     colnames(dat) <- c("event", "arm")
    ##     dat$event <- as.numeric(dat$event)

    ##     cart <- glm(event ~ arm, dat, family = quasipoisson)
    ##                                     # cart <- glm(event ~ arm, dat, family = poisson)

    ##     result[i] <- coef(summary(cart))[2, 4]
    ## }

    result <- unlist(parallel::mclapply(1:n_sim,
                                        function(x) {
                                            dat_1 <- rqpois(sample_size_AC, rate_AC, over_disp_rate + 1)
                                            dat_2 <- rqpois(sample_size_PLB, rate_PLB, over_disp_rate + 1)

                                            dat <- rbind(cbind(dat_1, "AC"),
                                                         cbind(dat_2, "PLB"))
                                            dat <- as.data.frame(dat)
                                            colnames(dat) <- c("event", "arm")
                                            dat$event <- as.numeric(dat$event)

                                            cart <- glm(event ~ arm, dat, family = quasipoisson)
                                        # cart <- glm(event ~ arm, dat, family = poisson)

                                            coef(summary(cart))[2, 4]
                                        },
                                        mc.cores = n_core))

    
    sum(result <= (1 - conf.level)) / n_sim
}



##' @rdname getPower_count_ovdisp
##' @export
getSampleSize_count_ovdisp <- function(power, ss_ratio = 1,
                                       rate_AC, rate_PLB,
                                       over_disp_rate = 0.2, drop_rate = 0.1,
                                       conf.level = 0.95, n_sim = 10000, n_core = 6){

    check <- getPower_count_ovdisp(5000, ss_ratio,
                              rate_AC, rate_PLB,
                              over_disp_rate, drop_rate,
                              conf.level, n_sim, n_core) < power
    if (check) {
        message("total sample size > 5000; return 9999 instead")
        return(9999)
    }
    
    f <- function(total_sample_size) {
        getPower_count_ovdisp(total_sample_size, ss_ratio,
                              rate_AC, rate_PLB,
                              over_disp_rate, drop_rate,
                              conf.level, n_sim, n_core) - power
    }
    round(uniroot(f, c(10, 5000))$root)
}
