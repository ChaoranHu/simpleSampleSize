##' Basic Bayesian functions related to CSF simulation
##'
##' @description
##' Several basic working functions to generate posterior
##' sample for common data types. Currently, functions are
##' avaliable for data with normal, negative binomial,
##' binary distributions.
##'
##' `getPost_normal1`: normal likelihood with non-informative prior
##' (mean ~ normal(0, 10000), sd ~ unif(0, 8))
##'
##' `getPost_normal2`: normal likelihood with conjugated prior (normal inverse gamma)
##'
##' `getPost_ngbn`: negative binomial likelihood with non-informative prior
##' (rate ~ unif(0, 8), p ~ beta(1, 1))
##'
##' `getPost_binary`: Bernoulli likelihood with beta prior
##'
##' @param data input dataset
##' @param 

## normal likelihood with non-informative prior
## mean ~ normal(0, 10000)
## sd ~ unif(0, 8)

getPost_normal1 <- function(dat){
    N <- length(dat)

    model_string <- "model{
for(i in 1:N){
dat[i] ~ dnorm(mu, tau2)
}

tau2 <- 1/sigma2
sigma2 <- sigma*sigma
sigma ~ dunif(0, 8)

mu ~ dnorm(0, 1e-04)
}
"



    dataList <- list(dat = dat, N = N)
    model <- jags.model(file = textConnection(model_string),
                        data = dataList)

    update(model, n.iter = 1000)
    Nrep = 5000 # number of values to simulate

    posterior_sample <- coda.samples(model,
                                     variable.names = c("mu", "sigma"),
                                     n.iter = Nrep,
                                     thin = 5)

    posterior_sample[[1]][, 1]
}

## data <- rnorm(1000, 0, 1)
## mean(getPost_normal1(data))







## normal likelihood with conjugated prior (normal inverse gamma)
getPost_normal2 <- function(dat,
                            prior_mu = 0, prior_nu = 0.001,
                            prior_alpha = 0.001, prior_beta = 0.001){
    post_alpha <- prior_alpha + length(dat)/2
    post_beta <- prior_beta + 0.5 * sum((dat - mean(dat))^2) + (length(dat) * prior_nu / (prior_nu + length(dat))) * 0.5 * (mean(dat) - prior_mu)^2
    post_mu <- (prior_nu * prior_mu + length(dat) * mean(dat)) / (prior_nu + length(dat))
    post_nu <- prior_nu + length(dat)

    cart1 <- sqrt(invgamma::rinvgamma(1000, shape = post_alpha, scale = post_beta)/post_nu)
    if(any(is.nan(cart1))) warning("WARNING: NaN occur in sigma")
    sapply(1:1000, function(x) rnorm(1, post_mu, sd = cart1[x]))
}






## negative binomial likelihood with non-informative prior
## rate ~ unif(0, 8)
## p ~ beta(1, 1)
getPost_ngbn <- function(dat){
    N <- length(dat)

    model_string <- "model{
for(i in 1:N){
dat[i] ~ dnegbin(p, r)
}

r <- mu/((1-p)/p)

mu ~ dunif(0,8)
p ~ dbeta(1,1)

}"

    dataList <- list(dat = dat, N = N)
    model <- jags.model(file = textConnection(model_string),
                        data = dataList)

    update(model, n.iter = 1000)
    Nrep = 5000 # number of values to simulate

    posterior_sample <- coda.samples(model,
                                     variable.names = c("mu"),
                                     n.iter = Nrep,
                                     thin = 5)

    posterior_sample[[1]][, 1]
}

## data <- simpleSampleSize:::rqpois(1000, 3, 1.2)
## mean(getPost_ngbn(data))

## data <- simpleSampleSize:::rqpois(1000, 3*0.5, 1.2) * 2
## mean(getPost_ngbn(data))



## Bernoulli likelihood with beta prior
getPost_binary <- function(dat, prior_alpha = 1, prior_beta = 1) {
    post_alpha <- prior_alpha + sum(dat)
    post_beta <- prior_beta + length(dat) - sum(dat)
    rbeta(1000, post_alpha, post_beta)
}


## data <- rbinom(10000, size = 1, prob = 0.4)
## mean(getPost_binary(data))
