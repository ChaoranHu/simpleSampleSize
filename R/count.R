
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
##' @param glm_model specify which glm model used in simulation.
##' It can be c("negativebinomial", "quasipoisson").
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
                                  n_core = 6,
                                  glm_model = "negativebinomial") {
    total_sample_size <- total_sample_size * (1-drop_rate)
    sample_size_AC <- total_sample_size*(ss_ratio/(ss_ratio+1))
    sample_size_PLB <- total_sample_size - sample_size_AC
    

    result <- parallel::mclapply(1:n_sim,
                                 function(x) {
                                     tryCatch({
                                         dat_1 <- rqpois(sample_size_AC, rate_AC, over_disp_rate + 1)
                                         dat_2 <- rqpois(sample_size_PLB, rate_PLB, over_disp_rate + 1)
                                         dat <- rbind(cbind(dat_1, "AC"),
                                                      cbind(dat_2, "PLB"))
                                         dat <- as.data.frame(dat)
                                         colnames(dat) <- c("event", "arm")
                                         dat$event <- as.numeric(dat$event)
                                         
                                         if (glm_model == "quasipoisson") {
                                             cart <- glm(event ~ arm, dat, family = quasipoisson)
                                         }
                                         if (glm_model == "negativebinomial") {
                                             cart <- MASS::glm.nb(event ~ arm, data = dat)
                                         }
                                         coef(summary(cart))[2, 4]
                                     },
                                     error = function(e) {
                                         # message("Sim ", x, " error: ", conditionMessage(e))
                                         NA_real_
                                     })
                                 },
                                 mc.cores = n_core
                                 )

    result <- unlist(result)

                                        # check how many failed
    # n_failed <- sum(is.na(result))
    # if (n_failed > 0) message(n_failed, " of ", n_sim, " simulations failed, power will be calculated with success fit only. The impact should be minimal.")

                                        # downstream analysis on successful runs only
    result_clean <- result[!is.na(result)]
    
    sum(result_clean <= (1 - conf.level)) / length(result_clean)
}


##' @rdname getPower_count_ovdisp
##' @export
getSampleSize_count_ovdisp <- function(power, ss_ratio = 1,
                                       rate_AC, rate_PLB,
                                       over_disp_rate = 0.2, drop_rate = 0.1,
                                       conf.level = 0.95, n_sim = 10000, n_core = 6,
                                       glm_model = "negativebinomial"){

    check <- getPower_count_ovdisp(5000, ss_ratio,
                                   rate_AC, rate_PLB,
                                   over_disp_rate, drop_rate,
                                   conf.level, n_sim, n_core, glm_model) < power
    if (check) {
        message("total sample size > 5000; return 9999 instead")
        return(9999)
    }
    
    f <- function(total_sample_size) {
        getPower_count_ovdisp(total_sample_size, ss_ratio,
                              rate_AC, rate_PLB,
                              over_disp_rate, drop_rate,
                              conf.level, n_sim, n_core, glm_model) - power
    }
    round(uniroot(f, c(50, 5000))$root)
}


