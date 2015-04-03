################################################################
#
# Model fitting code for paper "Using Multilevel Models to Estimate Variation in Foraging Returns"
# Richard McElreath and Jeremy Koster
#
# Nov 14 2012
#
# Additional R packages required: rstan, chron
#
################################################################

##########################
# load and prep data

d <- read.csv("hunting.csv")

# fix some day coding errors in raw data
d$day[ d$day==3157 ] <- 31 # guess Oct "3157" 2001 is supposed to be Oct "31" 2001
d$day[ d$day==-74 ] <- 15 # no idea, so guess middle of month (April 1997)

# make date variable, days since Jan 1 1970
library(chron)
d$date <- julian( x=d$month , d=d$day , y=d$year )
# Missing values! 6 of them. Drop those rows
d.orig <- d
d <- subset( d , !is.na(d$date) )
d$date.s <- (d$date - mean(d$date))/sd(d$date)
# make standardized age variables
d$age.s <- ( d$age - mean(d$age) ) / sd(d$age)
nid <- length(unique(d$id))
d$hunter.id <- as.integer( as.factor( d$id ) )
# outcome variables
d$iszero <- ifelse( d$kg.meat==0 , 1 , 0 )
d$nonzero <- ifelse( d$kg.meat > 0 , d$kg.meat , 1 )
# day ids
Ndays <- length( unique(d$date) )
d$dayid <- as.integer( as.factor(d$date) )

# load rstan
library(rstan)

#################################
# Fixed effects model

# model code
ache_code_Fix <- '
  data {
    int<lower=0> N; // number of trips
    int<lower=0> J; // number of hunters
    real y[N]; // hunting returns
    int<lower=0,upper=1> iszero[N]; // indicator of zero return
    real nonzero[N]; // positive returns
    int hid[N]; // hunter ids
    real ages[N]; // standardized age at time of hunt
    real ages2[N]; // square of age
    real ages3[N]; // cube of age
  } 
  parameters {
    real theta; // log rate of gamma
    real alpha[4];
    real beta[4];
  }
  model {
    real pi;
    real mugamma;
    for ( k in 1:4 ) {
        alpha[k] ~ normal(0,100);
        beta[k] ~ normal(0,100);
    }
    for( i in 1:N ) {
        pi <- alpha[1] + alpha[2]*ages[i] + alpha[3]*ages2[i] + alpha[4]*ages3[i];
        mugamma <- beta[1] + beta[2]*ages[i] + beta[3]*ages2[i] + beta[4]*ages3[i];
        iszero[i] ~ bernoulli( inv_logit(pi) );
        for ( g in 1:(1-iszero[i]) )
            nonzero[i] ~ gamma( exp(mugamma)*exp(theta) , exp(theta) );
        //lp__ <- lp__ + if_else( iszero[i] , log(inv_logit(pi)) , log1m(inv_logit(pi)) + gamma_log( nonzero[i] , exp(mugamma)*exp(theta) , exp(theta) ) );
    }
  }
  generated quantities {
        real dev;
        real pi;
        real mugamma;
        dev <- 0;
        for ( i in 1:N ) {
            pi <- alpha[1] + alpha[2]*ages[i] + alpha[3]*ages2[i] + alpha[4]*ages3[i];
            mugamma <- beta[1] + beta[2]*ages[i] + beta[3]*ages2[i] + beta[4]*ages3[i];
            dev <- dev + (-2) * if_else( iszero[i] , log(inv_logit(pi)) , log1m(inv_logit(pi)) + gamma_log( nonzero[i] , exp(mugamma)*exp(theta) , exp(theta) ) );
        }
    }
'

# prep data
ache_dat0 <- list(
    N = nrow(d),
    J = 147,
    y = d$kg.meat,
    iszero = d$iszero,
    nonzero = d$nonzero,
    hid = d$hunter.id,
    ages = d$age.s,
    ages2 = d$age.s^2,
    ages3 = d$age.s^3
)

# initialize
initlist <- list(list(
    theta = log(0.3) ,
    alpha = c(0,0,0,0),
    beta = c(log(7),0,0,0)
))

fitFix <- stan( model_code = ache_code_Fix , data = ache_dat0 , iter = 6000 , warmup=1000 , init=initlist , chains = 1 )

#################################
# Vary_2

# model code
ache_code_Vary2 <- '
  data {
    int<lower=0> N; // number of trips
    int<lower=0> J; // number of hunters
    int<lower=0> K; // number of varying effects
    real y[N]; // hunting returns
    int<lower=0,upper=1> iszero[N]; // indicator of zero return
    real nonzero[N]; // positive returns
    int hid[N]; // hunter ids
    real ages[N]; // standardized age at time of hunt
    real ages2[N]; // square of age
    real ages3[N]; // cube of age
    matrix[K,K] W;
  } 
  parameters {
    real theta; // log rate of gamma
    real alphaA1; // age effects
    real alphaA2; // quadratic age effect
    real alphaA3; // cubic age effect
    real betaA1; // age effects
    real betaA2; // quadratic age effect
    real betaA3; // cubic age effect
    vector[K] mu;
    cov_matrix[K] Sigma;
    vector[K] slopes[J];
  }
  model {
    real pi;
    real mugamma;
    alphaA1 ~ normal(0,100);
    alphaA2 ~ normal(0,100);
    alphaA3 ~ normal(0,100);
    betaA1 ~ normal(0,100);
    betaA2 ~ normal(0,100);
    betaA3 ~ normal(0,100);
    Sigma ~ inv_wishart( K+1 , W );
    for ( k in 1:K ) mu[k] ~ normal(0,100);
    for( j in 1:J ) slopes[j] ~ multi_normal( mu , Sigma );
    for( i in 1:N ) {
        pi <- slopes[hid[i],1] + alphaA1*ages[i] + alphaA2*ages2[i] + alphaA3*ages3[i];
        mugamma <- slopes[hid[i],2] + betaA1*ages[i] + betaA2*ages2[i] + betaA3*ages3[i];
        iszero[i] ~ bernoulli( inv_logit(pi) );
        for ( g in 1:(1-iszero[i]) )
            nonzero[i] ~ gamma( exp(mugamma)*exp(theta) , exp(theta) );
        //lp__ <- lp__ + if_else( iszero[i] , log(inv_logit(pi)) , log1m(inv_logit(pi)) + gamma_log( nonzero[i] , exp(mugamma)*exp(theta) , exp(theta) ) );
    }
  }
  generated quantities {
        real dev;
        real pi;
        real mugamma;
        dev <- 0;
        for ( i in 1:N ) {
            pi <- slopes[hid[i],1] + alphaA1*ages[i] + alphaA2*ages2[i] + alphaA3*ages3[i];
            mugamma <- slopes[hid[i],2] + betaA1*ages[i] + betaA2*ages2[i] + betaA3*ages3[i];
            dev <- dev + (-2) * if_else( iszero[i] , log(inv_logit(pi)) , log1m(inv_logit(pi)) + gamma_log( nonzero[i] , exp(mugamma)*exp(theta) , exp(theta) ) );
        }
    }
'

# prep data
K <- 2
ache_dat0 <- list(
    N = nrow(d),
    J = 147,
    K = K,
    y = d$kg.meat,
    iszero = d$iszero,
    nonzero = d$nonzero,
    hid = d$hunter.id,
    ages = d$age.s,
    ages2 = d$age.s^2,
    ages3 = d$age.s^3,
    W = diag(K)
)

# initialize
init_slopes = matrix( 0 , nrow=147 , ncol=2 )
init_slopes[ , 1 ] <- rep( 0 , 147 )
init_slopes[ , 2 ] <- rep( log(7) , 147 )
initlist <- list(list(
    theta = log(0.3) ,
    alphaA1 = 0 ,
    alphaA2 = 0 ,
    alphaA3 = 0 ,
    betaA1 = 0 ,
    betaA2 = 0 ,
    betaA3 = 0 ,
    mu = c( 0 , log(7) ) ,
    Sigma = diag(2) ,
    slopes = init_slopes
))

fitVary2 <- stan( model_code = ache_code_Vary2 , data = ache_dat0 , iter = 6000 , warmup=1000 , init=initlist , chains = 1 )

#################################
# Vary_4 (varying intercepts and first age terms)

# model code
ache_code_Vary4 <- '
  data {
    int<lower=0> N; // number of trips
    int<lower=0> J; // number of hunters
    int<lower=0> K; // number of varying effects
    real y[N]; // hunting returns
    int<lower=0,upper=1> iszero[N]; // indicator of zero return
    real nonzero[N]; // positive returns
    int hid[N]; // hunter ids
    real ages[N]; // standardized age at time of hunt
    real ages2[N]; // square of age
    real ages3[N]; // cube of age
    matrix[K,K] W;
  } 
  parameters {
    real theta; // log rate of gamma
    real alphaA2; // quadratic age effect
    real alphaA3; // cubic age effect
    real betaA2; // quadratic age effect
    real betaA3; // cubic age effect
    vector[K] mu;
    cov_matrix[K] Sigma;
    vector[K] slopes[J];
  }
  model {
    real pi;
    real mugamma;
    alphaA2 ~ normal(0,100);
    alphaA3 ~ normal(0,100);
    betaA2 ~ normal(0,100);
    betaA3 ~ normal(0,100);
    Sigma ~ inv_wishart( K+1 , W );
    for ( k in 1:K ) mu[k] ~ normal(0,100);
    for( j in 1:J ) slopes[j] ~ multi_normal( mu , Sigma );
    for( i in 1:N ) {
        pi <- slopes[hid[i],1] + slopes[hid[i],2]*ages[i] + alphaA2*ages2[i] + alphaA3*ages3[i];
        mugamma <- slopes[hid[i],3] + slopes[hid[i],4]*ages[i] + betaA2*ages2[i] + betaA3*ages3[i];
        lp__ <- lp__ + if_else( iszero[i] , log(inv_logit(pi)) , log1m(inv_logit(pi)) + gamma_log( nonzero[i] , exp(mugamma)*exp(theta) , exp(theta) ) );
    }
  }
  generated quantities {
        real dev;
        real pi;
        real mugamma;
        dev <- 0;
        for ( i in 1:N ) {
            pi <- slopes[hid[i],1] + slopes[hid[i],2]*ages[i] + alphaA2*ages2[i] + alphaA3*ages3[i];
            mugamma <- slopes[hid[i],3] + slopes[hid[i],4]*ages[i] + betaA2*ages2[i] + betaA3*ages3[i];
            dev <- dev + (-2) * if_else( iszero[i] , log(inv_logit(pi)) , log1m(inv_logit(pi)) + gamma_log( nonzero[i] , exp(mugamma)*exp(theta) , exp(theta) ) );
        }
    }
'

# prep data
K <- 4
ache_dat0 <- list(
    N = nrow(d),
    J = 147,
    K = K,
    y = d$kg.meat,
    iszero = d$iszero,
    nonzero = d$nonzero,
    hid = d$hunter.id,
    ages = d$age.s,
    ages2 = d$age.s^2,
    ages3 = d$age.s^3,
    W = diag(K)
)

# initialize
init_slopes <- matrix( 0 , nrow=147 , ncol=K )
vec <- c( 0,0, log(7) ,0 )
for ( j in 1:147 ) init_slopes[j,] <- vec
initlist <- list(list(
    mu = vec ,
    Sigma = diag(K) ,
    slopes = init_slopes,
    alphaA2 = 0,
    betaA2 = 0,
    alphaA3 = 0,
    betaA3 = 0,
    theta = log(0.3)
))

fitVary4 <- stan( model_code = ache_code_Vary4 , data = ache_dat0 , iter = 6000 , warmup=1000 , init=initlist , chains = 1 )

#################################
# Vary_6

# model code
ache_code_Vary6 <- '
  data {
    int<lower=0> N; // number of trips
    int<lower=0> J; // number of hunters
    int<lower=0> K; // number of varying effects
    real y[N]; // hunting returns
    int<lower=0,upper=1> iszero[N]; // indicator of zero return
    real nonzero[N]; // positive returns
    int hid[N]; // hunter ids
    real ages[N]; // standardized age at time of hunt
    real ages2[N]; // square of age
    real ages3[N]; // cube of age
    matrix[K,K] W;
  } 
  parameters {
    real theta; // log rate of gamma
    //real alphaA2; // quadratic age effect
    real alphaA3; // cubic age effect
    //real betaA2; // quadratic age effect
    real betaA3; // cubic age effect
    vector[K] mu;
    cov_matrix[K] Sigma;
    vector[K] slopes[J];
  }
  model {
    real pi;
    real mugamma;
    //alphaA2 ~ normal(0,100);
    alphaA3 ~ normal(0,100);
    //betaA2 ~ normal(0,100);
    betaA3 ~ normal(0,100);
    Sigma ~ inv_wishart( K+1 , W );
    for ( k in 1:K ) mu[k] ~ normal(0,100);
    for( j in 1:J ) slopes[j] ~ multi_normal( mu , Sigma );
    for( i in 1:N ) {
        pi <- slopes[hid[i],1] + slopes[hid[i],2]*ages[i] + slopes[hid[i],3]*ages2[i] + alphaA3*ages3[i];
        mugamma <- slopes[hid[i],4] + slopes[hid[i],5]*ages[i] + slopes[hid[i],6]*ages2[i] + betaA3*ages3[i];
        lp__ <- lp__ + if_else( iszero[i] , log(inv_logit(pi)) , log1m(inv_logit(pi)) + gamma_log( nonzero[i] , exp(mugamma)*exp(theta) + 0.01 , exp(theta) ) );
    }
  }
  generated quantities {
        real dev;
        real pi;
        real mugamma;
        dev <- 0;
        for ( i in 1:N ) {
            pi <- slopes[hid[i],1] + slopes[hid[i],2]*ages[i] + slopes[hid[i],3]*ages2[i] + alphaA3*ages3[i];
            mugamma <- slopes[hid[i],4] + slopes[hid[i],5]*ages[i] + slopes[hid[i],6]*ages2[i] + betaA3*ages3[i];
            dev <- dev + (-2) * if_else( iszero[i] , log(inv_logit(pi)) , log1m(inv_logit(pi)) + gamma_log( nonzero[i] , exp(mugamma)*exp(theta) + 0.01 , exp(theta) ) );
        }
    }
'

# prep data
K <- 6
ache_dat0 <- list(
    N = nrow(d),
    J = 147,
    K = K,
    y = d$kg.meat,
    iszero = d$iszero,
    nonzero = d$nonzero,
    hid = d$hunter.id,
    ages = d$age.s,
    ages2 = d$age.s^2,
    ages3 = d$age.s^3,
    W = diag(K)
)

# initialize
init_slopes <- matrix( 0 , nrow=147 , ncol=K )
vec <- c( 0,0,0, log(7) ,0,0 )
for ( j in 1:147 ) init_slopes[j,] <- vec
initlist <- list(list(
    mu = vec ,
    Sigma = diag(K) ,
    slopes = init_slopes,
    alphaA3 = 0,
    betaA3 = 0,
    theta = log(0.3)
))

fitVary6 <- stan( model_code = ache_code_Vary6 , data = ache_dat0 , iter = 6000 , warmup=1000 , init=initlist , chains = 1 )

#################################
# Vary_7 --- six random slopes for Age polynomials, 1 more for theta

# model code
ache_code_Vary7 <- '
  data {
    int<lower=0> N; // number of trips
    int<lower=0> J; // number of hunters
    int<lower=0> K; // number of varying effects
    real y[N]; // hunting returns
    int<lower=0,upper=1> iszero[N]; // indicator of zero return
    real nonzero[N]; // positive returns
    int hid[N]; // hunter ids
    real ages[N]; // standardized age at time of hunt
    real ages2[N]; // square of age
    real ages3[N]; // cube of age
    matrix[K,K] W;
  } 
  parameters {
    real alphaA3; // cubic age effect
    real betaA3; // cubic age effect
    vector[K] mu;
    cov_matrix[K] Sigma;
    vector[K] slopes[J];
  }
  model {
    real pi;
    real mugamma;
    real thetaJ;
    alphaA3 ~ normal(0,100);
    betaA3 ~ normal(0,100);
    Sigma ~ inv_wishart( K+1 , W );
    for ( k in 1:K ) mu[k] ~ normal(0,100);
    for( j in 1:J ) slopes[j] ~ multi_normal( mu , Sigma );
    for( i in 1:N ) {
        pi <- inv_logit( slopes[hid[i],1] + slopes[hid[i],2]*ages[i] + slopes[hid[i],3]*ages2[i] + alphaA3*ages3[i] );
        mugamma <- slopes[hid[i],4] + slopes[hid[i],5]*ages[i] + slopes[hid[i],6]*ages2[i] + betaA3*ages3[i];
        thetaJ <- exp( slopes[hid[i],7] );
        lp__ <- lp__ + if_else( iszero[i] , log(pi) , log1m(pi) + gamma_log( nonzero[i] , exp(mugamma)*thetaJ , thetaJ ) );
    }
  }
  generated quantities {
        real dev;
        real pi;
        real mugamma;
        real thetaJ;
        dev <- 0;
        for ( i in 1:N ) {
            pi <- inv_logit( slopes[hid[i],1] + slopes[hid[i],2]*ages[i] + slopes[hid[i],3]*ages2[i] + alphaA3*ages3[i] );
            mugamma <- slopes[hid[i],4] + slopes[hid[i],5]*ages[i] + slopes[hid[i],6]*ages2[i] + betaA3*ages3[i];
            thetaJ <- exp( slopes[hid[i],7] );
            dev <- dev + (-2) * if_else( iszero[i] , log(pi) , log1m(pi) + gamma_log( nonzero[i] , exp(mugamma)*thetaJ , thetaJ ) );
        }
    }
'

# prep data
K <- 7
ache_dat0 <- list(
    N = nrow(d),
    J = 147,
    K = K,
    y = d$kg.meat,
    iszero = d$iszero,
    nonzero = d$nonzero,
    hid = d$hunter.id,
    ages = d$age.s,
    ages2 = d$age.s^2,
    ages3 = d$age.s^3,
    W = diag(K)
)

# initialize
init_slopes <- matrix( 0 , nrow=147 , ncol=K )
vec <- c( 0,0,0, log(7),0,0, log(0.3) )
for ( j in 1:147 ) init_slopes[j,] <- vec
initlist <- list(list(
    mu = vec ,
    Sigma = diag(K) ,
    slopes = init_slopes,
    alphaA3 = 0,
    betaA3 = 0
))

fitVary7 <- stan( model_code = ache_code_Vary7 , data = ache_dat0 , iter = 6000 , warmup=1000 , init=initlist , chains = 1 )

#################################
# Vary_8

# model code
ache_code_Vary8 <- '
  data {
    int<lower=0> N; // number of trips
    int<lower=0> J; // number of hunters
    int<lower=0> K; // number of varying effects
    real y[N]; // hunting returns
    int<lower=0,upper=1> iszero[N]; // indicator of zero return
    real nonzero[N]; // positive returns
    int hid[N]; // hunter ids
    real ages[N]; // standardized age at time of hunt
    real ages2[N]; // square of age
    real ages3[N]; // cube of age
    cov_matrix[K] W;
  } 
  parameters {
    vector[K] mu;
    cov_matrix[K] Sigma;
    vector[K] slopes[J];
    real theta;
  }
  model {
    real pi;
    real mugamma;
    Sigma ~ inv_wishart( K+1 , W );
    for( k in 1:K ) mu[k] ~ normal( 0 , 100 );
    for( j in 1:J ) slopes[j] ~ multi_normal( mu , Sigma );
    for( i in 1:N ) {
        pi <- slopes[hid[i],1] + slopes[hid[i],2]*ages[i] + slopes[hid[i],3]*ages2[i] + slopes[hid[i],4]*ages3[i];
        mugamma <- slopes[hid[i],5] + slopes[hid[i],6]*ages[i] + slopes[hid[i],7]*ages2[i] + slopes[hid[i],8]*ages3[i];
        lp__ <- lp__ + if_else( iszero[i] , log(inv_logit(pi)) , log1m(inv_logit(pi)) + gamma_log( nonzero[i] , exp(mugamma)*exp( theta )+0.001 , exp( theta ) ) );
    }
  }
  generated quantities {
        real dev;
        real pi;
        real mugamma;
        dev <- 0;
        for ( i in 1:N ) {
            pi <- slopes[hid[i],1] + slopes[hid[i],2]*ages[i] + slopes[hid[i],3]*ages2[i] + slopes[hid[i],4]*ages3[i];
            mugamma <- slopes[hid[i],5] + slopes[hid[i],6]*ages[i] + slopes[hid[i],7]*ages2[i] + slopes[hid[i],8]*ages3[i];
            dev <- dev + (-2) * if_else( iszero[i] , log(inv_logit(pi)) , log1m(inv_logit(pi)) + gamma_log( nonzero[i] , exp(mugamma)*exp(theta)+0.001 , exp(theta) ) );
        }
    }
'

# prep data
K <- 8
ache_dat0 <- list(
    N = nrow(d),
    J = 147,
    K = K,
    y = d$kg.meat,
    iszero = d$iszero,
    nonzero = d$nonzero,
    hid = d$hunter.id,
    ages = d$age.s,
    ages2 = d$age.s^2,
    ages3 = d$age.s^3,
    W = diag(K)
)

# initialize
init_slopes <- matrix( 0 , nrow=147 , ncol=K )
vec <- c( 0,0.01,0.01,0.01, log(7) ,0.01,0.01,0.01 )
for ( j in 1:147 ) init_slopes[j,] <- vec
initlist <- list(list(
    mu = vec ,
    Sigma = diag(K) ,
    slopes = init_slopes ,
    theta = log(0.3)
))

fitVary8 <- stan( model_code = ache_code_Vary8 , data = ache_dat0 , iter = 12000 , warmup=2000 , init=initlist , chains = 1 )

#################################
# Vary_9

# model code
ache_code_Vary9 <- '
  data {
    int<lower=0> N; // number of trips
    int<lower=0> J; // number of hunters
    int<lower=0> K; // number of varying effects
    real y[N]; // hunting returns
    int<lower=0,upper=1> iszero[N]; // indicator of zero return
    real nonzero[N]; // positive returns
    int hid[N]; // hunter ids
    real ages[N]; // standardized age at time of hunt
    real ages2[N]; // square of age
    real ages3[N]; // cube of age
    matrix[K,K] W;
  } 
  parameters {
    vector[K] mu;
    cov_matrix[K] Sigma;
    vector[K] slopes[J];
  }
  model {
    real pi;
    real mugamma;
    real thetaJ;
    Sigma ~ inv_wishart( K+1 , W );
    for( k in 1:K ) mu[k] ~ normal( 0 , 100 );
    for( j in 1:J ) slopes[j] ~ multi_normal( mu , Sigma);
    for( i in 1:N ) {
        pi <- slopes[hid[i],1] + slopes[hid[i],2]*ages[i] + slopes[hid[i],3]*ages2[i] + slopes[hid[i],4]*ages3[i];
        mugamma <- slopes[hid[i],5] + slopes[hid[i],6]*ages[i] + slopes[hid[i],7]*ages2[i] + slopes[hid[i],8]*ages3[i];
        thetaJ <- exp( slopes[hid[i],9] );
        iszero[i] ~ bernoulli( inv_logit(pi) );
        for ( g in 1:(1-iszero[i]) )
            nonzero[i] ~ gamma( exp(mugamma)*thetaJ , thetaJ );
        //lp__ <- lp__ + if_else( iszero[i] , log(inv_logit(pi)) , log1m(inv_logit(pi)) + gamma_log( nonzero[i] , exp(mugamma)*thetaJ , thetaJ ) );
    }
  }
  generated quantities {
        real dev;
        real pi;
        real mugamma;
        dev <- 0;
        for ( i in 1:N ) {
            pi <- slopes[hid[i],1] + slopes[hid[i],2]*ages[i] + slopes[hid[i],3]*ages2[i] + slopes[hid[i],4]*ages3[i];
            mugamma <- slopes[hid[i],5] + slopes[hid[i],6]*ages[i] + slopes[hid[i],7]*ages2[i] + slopes[hid[i],8]*ages3[i];
            dev <- dev + (-2) * if_else( iszero[i] , log(inv_logit(pi)) , log1m(inv_logit(pi)) + gamma_log( nonzero[i] , exp(mugamma)*exp(slopes[hid[i],9]) , exp(slopes[hid[i],9]) ) );
        }
    }
'

# prep data
K <- 9
ache_dat0 <- list(
    N = nrow(d),
    J = 147,
    K = K,
    y = d$kg.meat,
    iszero = d$iszero,
    nonzero = d$nonzero,
    hid = d$hunter.id,
    ages = d$age.s,
    ages2 = d$age.s^2,
    ages3 = d$age.s^3,
    W = diag(K)
)

# initialize
init_slopes <- matrix( 0 , nrow=147 , ncol=9 )
vec <- c( 0,0.01,0.01,0.01, log(7) ,0.01,0.01,0.01, log(0.3) )
for ( j in 1:147 ) init_slopes[j,] <- vec
initlist <- list(list(
    mu = vec ,
    Sigma = diag(9) ,
    slopes = init_slopes
))

fitVary9 <- stan( model_code = ache_code_Vary9 , data = ache_dat0 , iter = 11000 , warmup=1000 , init=initlist , chains = 1 , pars=c("mu","slopes","Sigma","dev") )

#################################
# Vary_9D - includes Julian date

# model code
ache_code_Vary9D <- '
  data {
    int<lower=0> N; // number of trips
    int<lower=0> J; // number of hunters
    int<lower=0> K; // number of varying effects
    real y[N]; // hunting returns
    int<lower=0,upper=1> iszero[N]; // indicator of zero return
    real nonzero[N]; // positive returns
    int hid[N]; // hunter ids
    real ages[N]; // standardized age at time of hunt
    real ages2[N]; // square of age
    real ages3[N]; // cube of age
    matrix[K,K] W;
    real date[N];
  } 
  parameters {
    vector[K] mu;
    cov_matrix[K] Sigma;
    vector[K] slopes[J];
    real alpha_date;
    real beta_date;
  }
  model {
    real pi;
    real mugamma;
    real thetaJ;
    Sigma ~ inv_wishart( K+1 , W );
    alpha_date ~ normal( 0 , 100 );
    beta_date ~ normal( 0 , 100 );
    for( k in 1:K ) mu[k] ~ normal( 0 , 100 );
    for( j in 1:J ) slopes[j] ~ multi_normal( mu , Sigma );
    for( i in 1:N ) {
        pi <- slopes[hid[i],1] + slopes[hid[i],2]*ages[i] + slopes[hid[i],3]*ages2[i] + slopes[hid[i],4]*ages3[i] + alpha_date*date[i];
        mugamma <- slopes[hid[i],5] + slopes[hid[i],6]*ages[i] + slopes[hid[i],7]*ages2[i] + slopes[hid[i],8]*ages3[i] + beta_date*date[i];
        thetaJ <- exp( slopes[hid[i],9] );
        iszero[i] ~ bernoulli( inv_logit(pi) );
        for ( g in 1:(1-iszero[i]) )
            nonzero[i] ~ gamma( exp(mugamma)*thetaJ , thetaJ );
        //lp__ <- lp__ + if_else( iszero[i] , log(inv_logit(pi)) , log1m(inv_logit(pi)) + gamma_log( nonzero[i] , exp(mugamma)*thetaJ , thetaJ ) );
    }
  }
  generated quantities {
        real dev;
        real pi;
        real mugamma;
        dev <- 0;
        for ( i in 1:N ) {
            pi <- slopes[hid[i],1] + slopes[hid[i],2]*ages[i] + slopes[hid[i],3]*ages2[i] + slopes[hid[i],4]*ages3[i] + alpha_date*date[i];
            mugamma <- slopes[hid[i],5] + slopes[hid[i],6]*ages[i] + slopes[hid[i],7]*ages2[i] + slopes[hid[i],8]*ages3[i] + beta_date*date[i];
            dev <- dev + (-2) * if_else( iszero[i] , log(inv_logit(pi)) , log1m(inv_logit(pi)) + gamma_log( nonzero[i] , exp(mugamma)*exp(slopes[hid[i],9]) , exp(slopes[hid[i],9]) ) );
        }
    }
'

# prep data
K <- 9
ache_dat0 <- list(
    N = nrow(d),
    J = 147,
    K = K,
    y = d$kg.meat,
    iszero = d$iszero,
    nonzero = d$nonzero,
    hid = d$hunter.id,
    ages = d$age.s,
    ages2 = d$age.s^2,
    ages3 = d$age.s^3,
    W = diag(K),
    date = d$date.s
)

# initialize
init_slopes <- matrix( 0 , nrow=147 , ncol=9 )
vec <- c( 0,0.01,0.01,0.01, log(7) ,0.01,0.01,0.01, log(0.3) )
for ( j in 1:147 ) init_slopes[j,] <- vec
initlist <- list(list(
    mu = vec ,
    Sigma = diag(9) ,
    slopes = init_slopes,
    alpha_date = 0,
    beta_date = 0
))

fitVary9D <- stan( model_code = ache_code_Vary9D , data = ache_dat0 , iter = 12000 , warmup=2000 , init=initlist , chains = 1 , pars=c("mu","alpha_date","beta_date","slopes","Sigma","dev") )

#################################
# Vary_9h
# need to impute missing hours values

# split data into two parts: hours observed and hours missing
obs_idx <- which( !is.na( d$hours ) )
miss_idx <- which( is.na( d$hours ) )
N_obs <- length(obs_idx)
N_miss <- length(d$hours) - N_obs

# prep data
K <- 9
ache_dat_h <- list(
    N = nrow(d),
    N_obs = N_obs,
    N_miss = N_miss,
    J = 147,
    K = K,
    y_obs = d$kg.meat[ obs_idx ],
    y_miss = d$kg.meat[ miss_idx ],
    iszero_obs = d$iszero[ obs_idx ],
    iszero_miss = d$iszero[ miss_idx ],
    nonzero_obs = d$nonzero[ obs_idx ],
    nonzero_miss = d$nonzero[ miss_idx ],
    hid_obs = d$hunter.id[ obs_idx ],
    hid_miss = d$hunter.id[ miss_idx ],
    ages_obs = d$age.s[ obs_idx ],
    ages2_obs = (d$age.s^2)[ obs_idx ],
    ages3_obs = (d$age.s^3)[ obs_idx ],
    ages_miss = d$age.s[ miss_idx ],
    ages2_miss = (d$age.s^2)[ miss_idx ],
    ages3_miss = (d$age.s^3)[ miss_idx ],
    hours_obs = d$hours[ obs_idx ],
    W = diag(K)
)

# model code
ache_code_Vary9h <- '
  data {
    int<lower=0> N; // number of trips
    int<lower=0> N_obs; // number of trips with observed hours
    int<lower=0> N_miss; // number of trips with missing hours
    int<lower=0> J; // number of hunters
    int<lower=0> K; // number of varying effects
    real y_obs[N_obs]; // hunting returns
    real y_miss[N_miss];
    int<lower=0,upper=1> iszero_obs[N_obs]; // indicator of zero return
    int<lower=0,upper=1> iszero_miss[N_miss];
    real nonzero_obs[N_obs]; // positive returns
    real nonzero_miss[N_miss];
    int hid_obs[N_obs]; // hunter ids
    int hid_miss[N_miss];
    real ages_obs[N_obs]; // standardized age at time of hunt
    real ages2_obs[N_obs]; // square of age
    real ages3_obs[N_obs]; // cube of age
    real ages_miss[N_miss];
    real ages2_miss[N_miss];
    real ages3_miss[N_miss];
    matrix[K,K] W;
    real<lower=0> hours_obs[N_obs]; // observed hours hunted values
  } 
  parameters {
    vector[K] mu;
    cov_matrix[K] Sigma;
    vector[K] slopes[J];
    real<lower=0> hours_miss[N_miss]; // missing values to impute
    real alpha_hours;
    real beta_hours;
    real mu_hours; // imputation populaton
    real<lower=0> sigma_hours;
  }
  model {
    real pi;
    real mugamma;
    real thetaJ;
    Sigma ~ inv_wishart( K+1 , W );
    mu_hours ~ normal( 7 , 100 );
    sigma_hours ~ uniform( 0, 20 );
    for ( k in 1:K ) mu[k] ~ normal(0,100);
    for( j in 1:J ) 
        slopes[j] ~ multi_normal( mu , Sigma );
    for ( i in 1:N_miss ) {
        hours_miss[i] ~ normal( mu_hours , sigma_hours );
        pi <- slopes[hid_miss[i],1] + slopes[hid_miss[i],2]*ages_miss[i] + slopes[hid_miss[i],3]*ages2_miss[i] + slopes[hid_miss[i],4]*ages3_miss[i] + alpha_hours*hours_miss[i];
        mugamma <- exp( slopes[hid_miss[i],5] + slopes[hid_miss[i],6]*ages_miss[i] + slopes[hid_miss[i],7]*ages2_miss[i] + slopes[hid_miss[i],8]*ages3_miss[i] + beta_hours*hours_miss[i] );
        thetaJ <- exp( slopes[hid_miss[i],9] );
        iszero_miss[i] ~ bernoulli( inv_logit(pi) );
        for ( g in 1:(1-iszero_miss[i]) )
            nonzero_miss[i] ~ gamma( mugamma*thetaJ , thetaJ );
        //lp__ <- lp__ + if_else( iszero_miss[i] , log(inv_logit(pi)) , log1m(inv_logit(pi)) + gamma_log( nonzero_miss[i] , mugamma*thetaJ , thetaJ ) );
    }
    for( i in 1:N_obs ) {
        hours_obs[i] ~ normal( mu_hours , sigma_hours );
        pi <- slopes[hid_obs[i],1] + slopes[hid_obs[i],2]*ages_obs[i] + slopes[hid_obs[i],3]*ages2_obs[i] + slopes[hid_obs[i],4]*ages3_obs[i] + alpha_hours*hours_obs[i];
        mugamma <- exp( slopes[hid_obs[i],5] + slopes[hid_obs[i],6]*ages_obs[i] + slopes[hid_obs[i],7]*ages2_obs[i] + slopes[hid_obs[i],8]*ages3_obs[i] + beta_hours*hours_obs[i] );
        thetaJ <- exp( slopes[hid_obs[i],9] );
        iszero_obs[i] ~ bernoulli( inv_logit(pi) );
        for ( g in 1:(1-iszero_obs[i]) )
            nonzero_obs[i] ~ gamma( mugamma*thetaJ , thetaJ );
        //lp__ <- lp__ + if_else( iszero_obs[i] , log(inv_logit(pi)) , log1m(inv_logit(pi)) + gamma_log( nonzero_obs[i] , mugamma*thetaJ , thetaJ ) );
    }
  }

  generated quantities {
        real dev;
        real pi;
        real thetaJ;
        real mugamma;
        dev <- 0;
        for( i in 1:N_obs ) {
            pi <- slopes[hid_obs[i],1] + slopes[hid_obs[i],2]*ages_obs[i] + slopes[hid_obs[i],3]*ages2_obs[i] + slopes[hid_obs[i],4]*ages3_obs[i] + alpha_hours*hours_obs[i];
            mugamma <- exp( slopes[hid_obs[i],5] + slopes[hid_obs[i],6]*ages_obs[i] + slopes[hid_obs[i],7]*ages2_obs[i] + slopes[hid_obs[i],8]*ages3_obs[i] + beta_hours*hours_obs[i] );
            thetaJ <- exp( slopes[hid_obs[i],9] );
            dev <- dev + (-2) * if_else( iszero_obs[i] , log(inv_logit(pi)) , log1m(inv_logit(pi)) + gamma_log( nonzero_obs[i] , mugamma*thetaJ , thetaJ ) );
        }
        for ( i in 1:N_miss ) {
            pi <- slopes[hid_miss[i],1] + slopes[hid_miss[i],2]*ages_miss[i] + slopes[hid_miss[i],3]*ages2_miss[i] + slopes[hid_miss[i],4]*ages3_miss[i] + alpha_hours*hours_miss[i];
            mugamma <- exp( slopes[hid_miss[i],5] + slopes[hid_miss[i],6]*ages_miss[i] + slopes[hid_miss[i],7]*ages2_miss[i] + slopes[hid_miss[i],8]*ages3_miss[i] + beta_hours*hours_miss[i] );
            thetaJ <- exp( slopes[hid_miss[i],9] );
            dev <- dev + (-2) * if_else( iszero_miss[i] , log(inv_logit(pi)) , log1m(inv_logit(pi)) + gamma_log( nonzero_miss[i] , mugamma*thetaJ , thetaJ ) );
        }
    }
'

# initialize
init_slopes <- matrix( 0 , nrow=147 , ncol=9 )
vec <- c( 0,0,0,0, log(7) ,0,0,0, log(0.3) )
for ( j in 1:147 ) init_slopes[j,] <- vec
initlist <- list(list(
    mu = vec ,
    Sigma = diag(9) ,
    slopes = init_slopes,
    alpha_hours = 0,
    beta_hours = 0,
    mu_hours = mean(d$hours,na.rm=TRUE),
    sigma_hours = sd(d$hours,na.rm=TRUE),
    hours_miss = rep( mean(d$hours,na.rm=TRUE) , N_miss )
))

fitVary9h <- stan( model_code = ache_code_Vary9h , data = ache_dat_h , init=initlist , iter = 8000 , warmup=3000 , chains = 1 , pars=c("mu","slopes","Sigma","alpha_hours","beta_hours","mu_hours","sigma_hours","hours_miss","dev") )

#################################
# Vary_9hD
# hours and Julian date

# split data into two parts: hours observed and hours missing
obs_idx <- which( !is.na( d$hours ) )
miss_idx <- which( is.na( d$hours ) )
N_obs <- length(obs_idx)
N_miss <- length(d$hours) - N_obs

# prep data
K <- 9
ache_dat_h <- list(
    N = nrow(d),
    N_obs = N_obs,
    N_miss = N_miss,
    J = 147,
    K = K,
    y_obs = d$kg.meat[ obs_idx ],
    y_miss = d$kg.meat[ miss_idx ],
    iszero_obs = d$iszero[ obs_idx ],
    iszero_miss = d$iszero[ miss_idx ],
    nonzero_obs = d$nonzero[ obs_idx ],
    nonzero_miss = d$nonzero[ miss_idx ],
    hid_obs = d$hunter.id[ obs_idx ],
    hid_miss = d$hunter.id[ miss_idx ],
    ages_obs = d$age.s[ obs_idx ],
    ages2_obs = (d$age.s^2)[ obs_idx ],
    ages3_obs = (d$age.s^3)[ obs_idx ],
    ages_miss = d$age.s[ miss_idx ],
    ages2_miss = (d$age.s^2)[ miss_idx ],
    ages3_miss = (d$age.s^3)[ miss_idx ],
    hours_obs = d$hours[ obs_idx ],
    W = diag(K),
    date_obs = d$date.s[ obs_idx ],
    date_miss = d$date.s[ miss_idx ]
)

# model code
ache_code_Vary9hD <- '
  data {
    int<lower=0> N; // number of trips
    int<lower=0> N_obs; // number of trips with observed hours
    int<lower=0> N_miss; // number of trips with missing hours
    int<lower=0> J; // number of hunters
    int<lower=0> K; // number of varying effects
    real y_obs[N_obs]; // hunting returns
    real y_miss[N_miss];
    int<lower=0,upper=1> iszero_obs[N_obs]; // indicator of zero return
    int<lower=0,upper=1> iszero_miss[N_miss];
    real nonzero_obs[N_obs]; // positive returns
    real nonzero_miss[N_miss];
    int hid_obs[N_obs]; // hunter ids
    int hid_miss[N_miss];
    real ages_obs[N_obs]; // standardized age at time of hunt
    real ages2_obs[N_obs]; // square of age
    real ages3_obs[N_obs]; // cube of age
    real ages_miss[N_miss];
    real ages2_miss[N_miss];
    real ages3_miss[N_miss];
    matrix[K,K] W;
    real<lower=0> hours_obs[N_obs]; // observed hours hunted values
    real date_obs[N_obs];
    real date_miss[N_miss];
  } 
  parameters {
    vector[K] mu;
    cov_matrix[K] Sigma;
    vector[K] slopes[J];
    real<lower=0> hours_miss[N_miss]; // missing values to impute
    real alpha_hours;
    real beta_hours;
    real mu_hours; // imputation populaton
    real<lower=0> sigma_hours;
    real alpha_date;
    real beta_date;
  }
  model {
    real pi;
    real mugamma;
    real thetaJ;
    Sigma ~ inv_wishart( K+1 , W );
    mu_hours ~ normal( 7 , 100 );
    sigma_hours ~ uniform( 0, 20 );
    alpha_date ~ normal( 0 , 100 );
    beta_date ~ normal( 0 , 100 );
    for ( k in 1:K ) mu[k] ~ normal(0,100);
    for( j in 1:J ) slopes[j] ~ multi_normal( mu , Sigma );
    for ( i in 1:N_miss ) {
        hours_miss[i] ~ normal( mu_hours , sigma_hours );
        pi <- slopes[hid_miss[i],1] + slopes[hid_miss[i],2]*ages_miss[i] + slopes[hid_miss[i],3]*ages2_miss[i] + slopes[hid_miss[i],4]*ages3_miss[i] + alpha_hours*hours_miss[i] + alpha_date*date_miss[i];
        mugamma <- exp( slopes[hid_miss[i],5] + slopes[hid_miss[i],6]*ages_miss[i] + slopes[hid_miss[i],7]*ages2_miss[i] + slopes[hid_miss[i],8]*ages3_miss[i] + beta_hours*hours_miss[i] + beta_date*date_miss[i] );
        thetaJ <- exp( slopes[hid_miss[i],9] );
        iszero_miss[i] ~ bernoulli( inv_logit(pi) );
        for ( g in 1:(1-iszero_miss[i]) )
            nonzero_miss[i] ~ gamma( mugamma*thetaJ , thetaJ );
        //lp__ <- lp__ + if_else( iszero_miss[i] , log(inv_logit(pi)) , log1m(inv_logit(pi)) + gamma_log( nonzero_miss[i] , mugamma*thetaJ , thetaJ ) );
    }
    for( i in 1:N_obs ) {
        hours_obs[i] ~ normal( mu_hours , sigma_hours );
        pi <- slopes[hid_obs[i],1] + slopes[hid_obs[i],2]*ages_obs[i] + slopes[hid_obs[i],3]*ages2_obs[i] + slopes[hid_obs[i],4]*ages3_obs[i] + alpha_hours*hours_obs[i] + alpha_date*date_obs[i];
        mugamma <- exp( slopes[hid_obs[i],5] + slopes[hid_obs[i],6]*ages_obs[i] + slopes[hid_obs[i],7]*ages2_obs[i] + slopes[hid_obs[i],8]*ages3_obs[i] + beta_hours*hours_obs[i] + beta_date*date_obs[i] );
        thetaJ <- exp( slopes[hid_obs[i],9] );
        iszero_obs[i] ~ bernoulli( inv_logit(pi) );
        for ( g in 1:(1-iszero_obs[i]) )
            nonzero_obs[i] ~ gamma( mugamma*thetaJ , thetaJ );
        //lp__ <- lp__ + if_else( iszero_obs[i] , log(inv_logit(pi)) , log1m(inv_logit(pi)) + gamma_log( nonzero_obs[i] , mugamma*thetaJ , thetaJ ) );
    }
  }

  generated quantities {
        real dev;
        real pi;
        real thetaJ;
        real mugamma;
        dev <- 0;
        for( i in 1:N_obs ) {
            pi <- slopes[hid_obs[i],1] + slopes[hid_obs[i],2]*ages_obs[i] + slopes[hid_obs[i],3]*ages2_obs[i] + slopes[hid_obs[i],4]*ages3_obs[i] + alpha_hours*hours_obs[i] + alpha_date*date_obs[i];
            mugamma <- exp( slopes[hid_obs[i],5] + slopes[hid_obs[i],6]*ages_obs[i] + slopes[hid_obs[i],7]*ages2_obs[i] + slopes[hid_obs[i],8]*ages3_obs[i] + beta_hours*hours_obs[i] + beta_date*date_obs[i] );
            thetaJ <- exp( slopes[hid_obs[i],9] );
            dev <- dev + (-2) * if_else( iszero_obs[i] , log(inv_logit(pi)) , log1m(inv_logit(pi)) + gamma_log( nonzero_obs[i] , mugamma*thetaJ , thetaJ ) );
        }
        for ( i in 1:N_miss ) {
            pi <- slopes[hid_miss[i],1] + slopes[hid_miss[i],2]*ages_miss[i] + slopes[hid_miss[i],3]*ages2_miss[i] + slopes[hid_miss[i],4]*ages3_miss[i] + alpha_hours*hours_miss[i] + alpha_date*date_miss[i];
            mugamma <- exp( slopes[hid_miss[i],5] + slopes[hid_miss[i],6]*ages_miss[i] + slopes[hid_miss[i],7]*ages2_miss[i] + slopes[hid_miss[i],8]*ages3_miss[i] + beta_hours*hours_miss[i] + beta_date*date_miss[i] );
            thetaJ <- exp( slopes[hid_miss[i],9] );
            dev <- dev + (-2) * if_else( iszero_miss[i] , log(inv_logit(pi)) , log1m(inv_logit(pi)) + gamma_log( nonzero_miss[i] , mugamma*thetaJ , thetaJ ) );
        }
    }
'

# initialize
init_slopes <- matrix( 0 , nrow=147 , ncol=9 )
vec <- c( 0,0,0,0, log(7) ,0,0,0, log(0.3) )
for ( j in 1:147 ) init_slopes[j,] <- vec
initlist <- list(list(
    mu = vec ,
    Sigma = diag(9) ,
    slopes = init_slopes,
    alpha_hours = 0,
    beta_hours = 0,
    mu_hours = mean(d$hours,na.rm=TRUE),
    sigma_hours = sd(d$hours,na.rm=TRUE),
    hours_miss = rep( mean(d$hours,na.rm=TRUE) , N_miss ),
    alpha_date = 0,
    beta_date = 0
))

fitVary9hD <- stan( model_code = ache_code_Vary9hD , data = ache_dat_h , init=initlist , iter = 12000 , warmup=2000 , chains = 1 , pars=c("mu","slopes","Sigma","alpha_date","beta_date","alpha_hours","beta_hours","mu_hours","sigma_hours","hours_miss","dev") )

#################################
# Vary_9Dv - includes Julian date and varying effect of day (cluster by day)

# model code
ache_code_Vary9Dv <- '
  data {
    int<lower=0> N; // number of trips
    int<lower=0> J; // number of hunters
    int<lower=0> K; // number of varying effects
    real y[N]; // hunting returns
    int<lower=0,upper=1> iszero[N]; // indicator of zero return
    real nonzero[N]; // positive returns
    int hid[N]; // hunter ids
    real ages[N]; // standardized age at time of hunt
    real ages2[N]; // square of age
    real ages3[N]; // cube of age
    matrix[K,K] W;
    real date[N];
    int<lower=0> Ndays;                  // number of unique days
    int<lower=0> dayid[N];               // unique day id for each trip
    vector[2] mu_days;                   // constant c(0,0) location vector for day effects
  } 
  parameters {
    vector[K] mu;
    cov_matrix[K] Sigma;
    vector[K] slopes[J];
    real alpha_date;
    real beta_date;
    vector[2] eta_day[Ndays];   // day specific effects on zeros
    cov_matrix[2] Sigma_day;    // var-cov of varying effects on unique day (means are zero!)
  }
  model {
    real pi;
    real mugamma;
    real thetaJ;
    Sigma ~ inv_wishart( K+1 , W );
    alpha_date ~ normal( 0 , 100 );
    beta_date ~ normal( 0 , 100 );
    for( k in 1:K ) mu[k] ~ normal( 0 , 100 );
    for( j in 1:J ) slopes[j] ~ multi_normal( mu , Sigma );
    for( j in 1:Ndays ) eta_day[j] ~ multi_normal( mu_days , Sigma_day );
    for( i in 1:N ) {
        pi <- slopes[hid[i],1] + slopes[hid[i],2]*ages[i] + slopes[hid[i],3]*ages2[i] + slopes[hid[i],4]*ages3[i] + alpha_date*date[i] + eta_day[ dayid[i] , 1 ];
        mugamma <- slopes[hid[i],5] + slopes[hid[i],6]*ages[i] + slopes[hid[i],7]*ages2[i] + slopes[hid[i],8]*ages3[i] + beta_date*date[i] + eta_day[ dayid[i] , 2 ];
        thetaJ <- exp( slopes[hid[i],9] );
        iszero[i] ~ bernoulli( inv_logit(pi) );
        for ( g in 1:(1-iszero[i]) )
            nonzero[i] ~ gamma( exp(mugamma)*thetaJ , thetaJ );
        //lp__ <- lp__ + if_else( iszero[i] , log(inv_logit(pi)) , log1m(inv_logit(pi)) + gamma_log( nonzero[i] , exp(mugamma)*thetaJ , thetaJ ) );
    }
  }
  generated quantities {
        real dev;
        real pi;
        real mugamma;
        dev <- 0;
        for ( i in 1:N ) {
            pi <- slopes[hid[i],1] + slopes[hid[i],2]*ages[i] + slopes[hid[i],3]*ages2[i] + slopes[hid[i],4]*ages3[i] + alpha_date*date[i] + eta_day[ dayid[i] , 1 ];
            mugamma <- slopes[hid[i],5] + slopes[hid[i],6]*ages[i] + slopes[hid[i],7]*ages2[i] + slopes[hid[i],8]*ages3[i] + beta_date*date[i] + eta_day[ dayid[i] , 2 ];
            dev <- dev + (-2) * if_else( iszero[i] , log(inv_logit(pi)) , log1m(inv_logit(pi)) + gamma_log( nonzero[i] , exp(mugamma)*exp(slopes[hid[i],9]) , exp(slopes[hid[i],9]) ) );
        }
    }
'

# prep data
K <- 9
Ndays <- length( unique(d$date) )
dayid <- as.numeric( as.factor(d$date) )
ache_dat0 <- list(
    N = nrow(d),
    J = 147,
    K = K,
    y = d$kg.meat,
    iszero = d$iszero,
    nonzero = d$nonzero,
    hid = d$hunter.id,
    ages = d$age.s,
    ages2 = d$age.s^2,
    ages3 = d$age.s^3,
    W = diag(K),
    date = d$date.s,
    Ndays = Ndays,
    dayid = dayid,
    mu_days = c(0,0)
)

# initialize
init_slopes <- matrix( 0 , nrow=147 , ncol=9 )
vec <- c( 0,0.01,0.01,0.01, log(7) ,0.01,0.01,0.01, log(0.3) )
for ( j in 1:147 ) init_slopes[j,] <- vec
initlist <- list(list(
    mu = vec ,
    Sigma = diag(9) ,
    slopes = init_slopes,
    alpha_date = 0,
    beta_date = 0,
    eta_day = matrix( 0 , nrow=Ndays , ncol=2 ),
    Sigma_day = diag(2)
))

fitVary9Dv <- stan( model_code = ache_code_Vary9Dv , data = ache_dat0 , iter = 7000 , warmup=2000 , init=initlist , chains = 1 , pars=c("mu","alpha_date","beta_date","slopes","Sigma","eta_day","Sigma_day","dev") )

#################################
# Vary_9v - includes varying effect of day (cluster by day)

# model code
ache_code_Vary9v <- '
  data {
    int<lower=0> N; // number of trips
    int<lower=0> J; // number of hunters
    int<lower=0> K; // number of varying effects
    real y[N]; // hunting returns
    int<lower=0,upper=1> iszero[N]; // indicator of zero return
    real nonzero[N]; // positive returns
    int hid[N]; // hunter ids
    real ages[N]; // standardized age at time of hunt
    real ages2[N]; // square of age
    real ages3[N]; // cube of age
    matrix[K,K] W;
    real date[N];
    int<lower=0> Ndays;                  // number of unique days
    int<lower=0> dayid[N];               // unique day id for each trip
    vector[2] mu_days;                   // constant c(0,0) location vector for day effects
  } 
  parameters {
    vector[K] mu;
    cov_matrix[K] Sigma;
    vector[K] slopes[J];
    vector[2] eta_day[Ndays];   // day specific effects on zeros
    cov_matrix[2] Sigma_day;    // var-cov of varying effects on unique day (means are zero!)
  }
  model {
    real pi;
    real mugamma;
    real thetaJ;
    Sigma ~ inv_wishart( K+1 , W );
    for( k in 1:K ) mu[k] ~ normal( 0 , 100 );
    for( j in 1:J ) slopes[j] ~ multi_normal( mu , Sigma );
    for( j in 1:Ndays ) eta_day[j] ~ multi_normal( mu_days , Sigma_day );
    for( i in 1:N ) {
        pi <- slopes[hid[i],1] + slopes[hid[i],2]*ages[i] + slopes[hid[i],3]*ages2[i] + slopes[hid[i],4]*ages3[i] + eta_day[ dayid[i] , 1 ];
        mugamma <- slopes[hid[i],5] + slopes[hid[i],6]*ages[i] + slopes[hid[i],7]*ages2[i] + slopes[hid[i],8]*ages3[i] + eta_day[ dayid[i] , 2 ];
        thetaJ <- exp( slopes[hid[i],9] );
        iszero[i] ~ bernoulli( inv_logit(pi) );
        for ( g in 1:(1-iszero[i]) )
            nonzero[i] ~ gamma( exp(mugamma)*thetaJ , thetaJ );
        //lp__ <- lp__ + if_else( iszero[i] , log(inv_logit(pi)) , log1m(inv_logit(pi)) + gamma_log( nonzero[i] , exp(mugamma)*thetaJ , thetaJ ) );
    }
  }
  generated quantities {
        real dev;
        real pi;
        real mugamma;
        dev <- 0;
        for ( i in 1:N ) {
            pi <- slopes[hid[i],1] + slopes[hid[i],2]*ages[i] + slopes[hid[i],3]*ages2[i] + slopes[hid[i],4]*ages3[i] + eta_day[ dayid[i] , 1 ];
            mugamma <- slopes[hid[i],5] + slopes[hid[i],6]*ages[i] + slopes[hid[i],7]*ages2[i] + slopes[hid[i],8]*ages3[i] + eta_day[ dayid[i] , 2 ];
            dev <- dev + (-2) * if_else( iszero[i] , log(inv_logit(pi)) , log1m(inv_logit(pi)) + gamma_log( nonzero[i] , exp(mugamma)*exp(slopes[hid[i],9]) , exp(slopes[hid[i],9]) ) );
        }
    }
'

# prep data
K <- 9
Ndays <- length( unique(d$date) )
dayid <- as.numeric( as.factor(d$date) )
ache_dat0 <- list(
    N = nrow(d),
    J = 147,
    K = K,
    y = d$kg.meat,
    iszero = d$iszero,
    nonzero = d$nonzero,
    hid = d$hunter.id,
    ages = d$age.s,
    ages2 = d$age.s^2,
    ages3 = d$age.s^3,
    W = diag(K),
    date = d$date.s,
    Ndays = Ndays,
    dayid = dayid,
    mu_days = c(0,0)
)

# initialize
init_slopes <- matrix( 0 , nrow=147 , ncol=9 )
vec <- c( 0,0.01,0.01,0.01, log(7) ,0.01,0.01,0.01, log(0.3) )
for ( j in 1:147 ) init_slopes[j,] <- vec
initlist <- list(list(
    mu = vec ,
    Sigma = diag(9) ,
    slopes = init_slopes,
    eta_day = matrix( 0 , nrow=Ndays , ncol=2 ),
    Sigma_day = diag(2)
))

fitVary9v <- stan( model_code = ache_code_Vary9v , data = ache_dat0 , iter = 7000 , warmup=2000 , init=initlist , chains = 1 , pars=c("mu","slopes","Sigma","eta_day","Sigma_day","dev") )

#################################
# Vary_9hv
# hours and varying effects on day (NO Julian date secular trend)

# split data into two parts: hours observed and hours missing
obs_idx <- which( !is.na( d$hours ) )
miss_idx <- which( is.na( d$hours ) )
N_obs <- length(obs_idx)
N_miss <- length(d$hours) - N_obs

# prep data
K <- 9
Ndays <- length( unique(d$date) )
dayid <- as.numeric( as.factor(d$date) )
ache_dat_h <- list(
    N = nrow(d),
    N_obs = N_obs,
    N_miss = N_miss,
    J = 147,
    K = K,
    y_obs = d$kg.meat[ obs_idx ],
    y_miss = d$kg.meat[ miss_idx ],
    iszero_obs = d$iszero[ obs_idx ],
    iszero_miss = d$iszero[ miss_idx ],
    nonzero_obs = d$nonzero[ obs_idx ],
    nonzero_miss = d$nonzero[ miss_idx ],
    hid_obs = d$hunter.id[ obs_idx ],
    hid_miss = d$hunter.id[ miss_idx ],
    ages_obs = d$age.s[ obs_idx ],
    ages2_obs = (d$age.s^2)[ obs_idx ],
    ages3_obs = (d$age.s^3)[ obs_idx ],
    ages_miss = d$age.s[ miss_idx ],
    ages2_miss = (d$age.s^2)[ miss_idx ],
    ages3_miss = (d$age.s^3)[ miss_idx ],
    hours_obs = d$hours[ obs_idx ],
    W = diag(K),
    date_obs = d$date.s[ obs_idx ],
    date_miss = d$date.s[ miss_idx ],
    Ndays = Ndays,
    dayid_obs = dayid[ obs_idx ],
    dayid_miss = dayid[ miss_idx ],
    mu_days = c(0,0)
)

# model code
ache_code_Vary9hv <- '
  data {
    int<lower=0> N; // number of trips
    int<lower=0> N_obs; // number of trips with observed hours
    int<lower=0> N_miss; // number of trips with missing hours
    int<lower=0> J; // number of hunters
    int<lower=0> K; // number of varying effects
    real y_obs[N_obs]; // hunting returns
    real y_miss[N_miss];
    int<lower=0,upper=1> iszero_obs[N_obs]; // indicator of zero return
    int<lower=0,upper=1> iszero_miss[N_miss];
    real nonzero_obs[N_obs]; // positive returns
    real nonzero_miss[N_miss];
    int hid_obs[N_obs]; // hunter ids
    int hid_miss[N_miss];
    real ages_obs[N_obs]; // standardized age at time of hunt
    real ages2_obs[N_obs]; // square of age
    real ages3_obs[N_obs]; // cube of age
    real ages_miss[N_miss];
    real ages2_miss[N_miss];
    real ages3_miss[N_miss];
    matrix[K,K] W;
    real<lower=0> hours_obs[N_obs]; // observed hours hunted values
    real date_obs[N_obs];
    real date_miss[N_miss];
    int Ndays;                      // number of unique days in sample
    int dayid_obs[N_obs];           // day ids corresponding to cases with observed hours
    int dayid_miss[N_miss];         // day ids corresponding to cases with MISSING hours
    vector[2] mu_days;              // c(0,0) location vector for day varying effects
  } 
  parameters {
    vector[K] mu;
    cov_matrix[K] Sigma;
    vector[K] slopes[J];
    real<lower=0> hours_miss[N_miss]; // missing values to impute
    real alpha_hours;
    real beta_hours;
    real mu_hours; // imputation populaton
    real<lower=0> sigma_hours;
    vector[2] eta_day[Ndays];
    cov_matrix[2] Sigma_day;
  }
  model {
    real pi;
    real mugamma;
    real thetaJ;
    Sigma ~ inv_wishart( K+1 , W );
    mu_hours ~ normal( 7 , 100 );
    sigma_hours ~ uniform( 0, 20 );
    for ( k in 1:K ) mu[k] ~ normal(0,100);
    for( j in 1:J ) slopes[j] ~ multi_normal( mu , Sigma );
    for( j in 1:Ndays ) eta_day[j] ~ multi_normal( mu_days , Sigma_day );
    for ( i in 1:N_miss ) {
        hours_miss[i] ~ normal( mu_hours , sigma_hours );
        pi <- slopes[hid_miss[i],1] + slopes[hid_miss[i],2]*ages_miss[i] + slopes[hid_miss[i],3]*ages2_miss[i] + slopes[hid_miss[i],4]*ages3_miss[i] + alpha_hours*hours_miss[i] + eta_day[ dayid_miss[i] , 1 ];
        mugamma <- exp( slopes[hid_miss[i],5] + slopes[hid_miss[i],6]*ages_miss[i] + slopes[hid_miss[i],7]*ages2_miss[i] + slopes[hid_miss[i],8]*ages3_miss[i] + beta_hours*hours_miss[i] + eta_day[ dayid_miss[i] , 2 ] );
        thetaJ <- exp( slopes[hid_miss[i],9] );
        iszero_miss[i] ~ bernoulli( inv_logit(pi) );
        for ( g in 1:(1-iszero_miss[i]) )
            nonzero_miss[i] ~ gamma( mugamma*thetaJ , thetaJ );
        //lp__ <- lp__ + if_else( iszero_miss[i] , log(inv_logit(pi)) , log1m(inv_logit(pi)) + gamma_log( nonzero_miss[i] , mugamma*thetaJ , thetaJ ) );
    }
    for( i in 1:N_obs ) {
        hours_obs[i] ~ normal( mu_hours , sigma_hours );
        pi <- slopes[hid_obs[i],1] + slopes[hid_obs[i],2]*ages_obs[i] + slopes[hid_obs[i],3]*ages2_obs[i] + slopes[hid_obs[i],4]*ages3_obs[i] + alpha_hours*hours_obs[i] + eta_day[ dayid_obs[i] , 1 ];
        mugamma <- exp( slopes[hid_obs[i],5] + slopes[hid_obs[i],6]*ages_obs[i] + slopes[hid_obs[i],7]*ages2_obs[i] + slopes[hid_obs[i],8]*ages3_obs[i] + beta_hours*hours_obs[i] + eta_day[ dayid_obs[i] , 2 ] );
        thetaJ <- exp( slopes[hid_obs[i],9] );
        iszero_obs[i] ~ bernoulli( inv_logit(pi) );
        for ( g in 1:(1-iszero_obs[i]) )
            nonzero_obs[i] ~ gamma( mugamma*thetaJ , thetaJ );
        //lp__ <- lp__ + if_else( iszero_obs[i] , log(inv_logit(pi)) , log1m(inv_logit(pi)) + gamma_log( nonzero_obs[i] , mugamma*thetaJ , thetaJ ) );
    }
  }

  generated quantities {
        real dev;
        real pi;
        real thetaJ;
        real mugamma;
        dev <- 0;
        for( i in 1:N_obs ) {
            pi <- slopes[hid_obs[i],1] + slopes[hid_obs[i],2]*ages_obs[i] + slopes[hid_obs[i],3]*ages2_obs[i] + slopes[hid_obs[i],4]*ages3_obs[i] + alpha_hours*hours_obs[i] + eta_day[ dayid_obs[i] , 1 ];
            mugamma <- exp( slopes[hid_obs[i],5] + slopes[hid_obs[i],6]*ages_obs[i] + slopes[hid_obs[i],7]*ages2_obs[i] + slopes[hid_obs[i],8]*ages3_obs[i] + beta_hours*hours_obs[i] + eta_day[ dayid_obs[i] , 2 ] );
            thetaJ <- exp( slopes[hid_obs[i],9] );
            dev <- dev + (-2) * if_else( iszero_obs[i] , log(inv_logit(pi)) , log1m(inv_logit(pi)) + gamma_log( nonzero_obs[i] , mugamma*thetaJ , thetaJ ) );
        }
        for ( i in 1:N_miss ) {
            pi <- slopes[hid_miss[i],1] + slopes[hid_miss[i],2]*ages_miss[i] + slopes[hid_miss[i],3]*ages2_miss[i] + slopes[hid_miss[i],4]*ages3_miss[i] + alpha_hours*hours_miss[i] + eta_day[ dayid_miss[i] , 1 ];
            mugamma <- exp( slopes[hid_miss[i],5] + slopes[hid_miss[i],6]*ages_miss[i] + slopes[hid_miss[i],7]*ages2_miss[i] + slopes[hid_miss[i],8]*ages3_miss[i] + beta_hours*hours_miss[i] + eta_day[ dayid_miss[i] , 2 ] );
            thetaJ <- exp( slopes[hid_miss[i],9] );
            dev <- dev + (-2) * if_else( iszero_miss[i] , log(inv_logit(pi)) , log1m(inv_logit(pi)) + gamma_log( nonzero_miss[i] , mugamma*thetaJ , thetaJ ) );
        }
    }
'

# initialize
init_slopes <- matrix( 0 , nrow=147 , ncol=9 )
vec <- c( 0,0,0,0, log(7) ,0,0,0, log(0.3) )
for ( j in 1:147 ) init_slopes[j,] <- vec
initlist <- list(list(
    mu = vec ,
    Sigma = diag(9) ,
    slopes = init_slopes,
    alpha_hours = 0,
    beta_hours = 0,
    mu_hours = mean(d$hours,na.rm=TRUE),
    sigma_hours = sd(d$hours,na.rm=TRUE),
    hours_miss = rep( mean(d$hours,na.rm=TRUE) , N_miss ),
    eta_day = matrix( 0 , nrow=Ndays , ncol=2 ),
    Sigma_day = diag(2)
))

fitVary9hv <- stan( model_code = ache_code_Vary9hv , data = ache_dat_h , init=initlist , iter = 7000 , warmup=2000 , chains = 1 , pars=c("mu","slopes","Sigma","alpha_hours","beta_hours","mu_hours","sigma_hours","hours_miss","eta_day","Sigma_day","dev") )

#################################
# Vary_9hDv
# hours and Julian date and varying effects on day

# split data into two parts: hours observed and hours missing
obs_idx <- which( !is.na( d$hours ) )
miss_idx <- which( is.na( d$hours ) )
N_obs <- length(obs_idx)
N_miss <- length(d$hours) - N_obs

# prep data
K <- 9
Ndays <- length( unique(d$date) )
dayid <- as.numeric( as.factor(d$date) )
ache_dat_h <- list(
    N = nrow(d),
    N_obs = N_obs,
    N_miss = N_miss,
    J = 147,
    K = K,
    y_obs = d$kg.meat[ obs_idx ],
    y_miss = d$kg.meat[ miss_idx ],
    iszero_obs = d$iszero[ obs_idx ],
    iszero_miss = d$iszero[ miss_idx ],
    nonzero_obs = d$nonzero[ obs_idx ],
    nonzero_miss = d$nonzero[ miss_idx ],
    hid_obs = d$hunter.id[ obs_idx ],
    hid_miss = d$hunter.id[ miss_idx ],
    ages_obs = d$age.s[ obs_idx ],
    ages2_obs = (d$age.s^2)[ obs_idx ],
    ages3_obs = (d$age.s^3)[ obs_idx ],
    ages_miss = d$age.s[ miss_idx ],
    ages2_miss = (d$age.s^2)[ miss_idx ],
    ages3_miss = (d$age.s^3)[ miss_idx ],
    hours_obs = d$hours[ obs_idx ],
    W = diag(K),
    date_obs = d$date.s[ obs_idx ],
    date_miss = d$date.s[ miss_idx ],
    Ndays = Ndays,
    dayid_obs = dayid[ obs_idx ],
    dayid_miss = dayid[ miss_idx ],
    mu_days = c(0,0)
)

# model code
ache_code_Vary9hDv <- '
  data {
    int<lower=0> N; // number of trips
    int<lower=0> N_obs; // number of trips with observed hours
    int<lower=0> N_miss; // number of trips with missing hours
    int<lower=0> J; // number of hunters
    int<lower=0> K; // number of varying effects
    real y_obs[N_obs]; // hunting returns
    real y_miss[N_miss];
    int<lower=0,upper=1> iszero_obs[N_obs]; // indicator of zero return
    int<lower=0,upper=1> iszero_miss[N_miss];
    real nonzero_obs[N_obs]; // positive returns
    real nonzero_miss[N_miss];
    int hid_obs[N_obs]; // hunter ids
    int hid_miss[N_miss];
    real ages_obs[N_obs]; // standardized age at time of hunt
    real ages2_obs[N_obs]; // square of age
    real ages3_obs[N_obs]; // cube of age
    real ages_miss[N_miss];
    real ages2_miss[N_miss];
    real ages3_miss[N_miss];
    matrix[K,K] W;
    real<lower=0> hours_obs[N_obs]; // observed hours hunted values
    real date_obs[N_obs];
    real date_miss[N_miss];
    int Ndays;                      // number of unique days in sample
    int dayid_obs[N_obs];           // day ids corresponding to cases with observed hours
    int dayid_miss[N_miss];         // day ids corresponding to cases with MISSING hours
    vector[2] mu_days;              // c(0,0) location vector for day varying effects
  } 
  parameters {
    vector[K] mu;
    cov_matrix[K] Sigma;
    vector[K] slopes[J];
    real<lower=0> hours_miss[N_miss]; // missing values to impute
    real alpha_hours;
    real beta_hours;
    real mu_hours; // imputation populaton
    real<lower=0> sigma_hours;
    real alpha_date;
    real beta_date;
    vector[2] eta_day[Ndays];
    cov_matrix[2] Sigma_day;
  }
  model {
    real pi;
    real mugamma;
    real thetaJ;
    Sigma ~ inv_wishart( K+1 , W );
    mu_hours ~ normal( 7 , 100 );
    sigma_hours ~ uniform( 0, 20 );
    alpha_date ~ normal( 0 , 100 );
    beta_date ~ normal( 0 , 100 );
    for ( k in 1:K ) mu[k] ~ normal(0,100);
    for( j in 1:J ) slopes[j] ~ multi_normal( mu , Sigma );
    for( j in 1:Ndays ) eta_day[j] ~ multi_normal( mu_days , Sigma_day );
    for ( i in 1:N_miss ) {
        hours_miss[i] ~ normal( mu_hours , sigma_hours );
        pi <- slopes[hid_miss[i],1] + slopes[hid_miss[i],2]*ages_miss[i] + slopes[hid_miss[i],3]*ages2_miss[i] + slopes[hid_miss[i],4]*ages3_miss[i] + alpha_hours*hours_miss[i] + alpha_date*date_miss[i] + eta_day[ dayid_miss[i] , 1 ];
        mugamma <- exp( slopes[hid_miss[i],5] + slopes[hid_miss[i],6]*ages_miss[i] + slopes[hid_miss[i],7]*ages2_miss[i] + slopes[hid_miss[i],8]*ages3_miss[i] + beta_hours*hours_miss[i] + beta_date*date_miss[i] + eta_day[ dayid_miss[i] , 2 ] );
        thetaJ <- exp( slopes[hid_miss[i],9] );
        iszero_miss[i] ~ bernoulli( inv_logit(pi) );
        for ( g in 1:(1-iszero_miss[i]) )
            nonzero_miss[i] ~ gamma( mugamma*thetaJ , thetaJ );
        //lp__ <- lp__ + if_else( iszero_miss[i] , log(inv_logit(pi)) , log1m(inv_logit(pi)) + gamma_log( nonzero_miss[i] , mugamma*thetaJ , thetaJ ) );
    }
    for( i in 1:N_obs ) {
        hours_obs[i] ~ normal( mu_hours , sigma_hours );
        pi <- slopes[hid_obs[i],1] + slopes[hid_obs[i],2]*ages_obs[i] + slopes[hid_obs[i],3]*ages2_obs[i] + slopes[hid_obs[i],4]*ages3_obs[i] + alpha_hours*hours_obs[i] + alpha_date*date_obs[i] + eta_day[ dayid_obs[i] , 1 ];
        mugamma <- exp( slopes[hid_obs[i],5] + slopes[hid_obs[i],6]*ages_obs[i] + slopes[hid_obs[i],7]*ages2_obs[i] + slopes[hid_obs[i],8]*ages3_obs[i] + beta_hours*hours_obs[i] + beta_date*date_obs[i] + eta_day[ dayid_obs[i] , 2 ] );
        thetaJ <- exp( slopes[hid_obs[i],9] );
        iszero_obs[i] ~ bernoulli( inv_logit(pi) );
        for ( g in 1:(1-iszero_obs[i]) )
            nonzero_obs[i] ~ gamma( mugamma*thetaJ , thetaJ );
        //lp__ <- lp__ + if_else( iszero_obs[i] , log(inv_logit(pi)) , log1m(inv_logit(pi)) + gamma_log( nonzero_obs[i] , mugamma*thetaJ , thetaJ ) );
    }
  }

  generated quantities {
        real dev;
        real pi;
        real thetaJ;
        real mugamma;
        dev <- 0;
        for( i in 1:N_obs ) {
            pi <- slopes[hid_obs[i],1] + slopes[hid_obs[i],2]*ages_obs[i] + slopes[hid_obs[i],3]*ages2_obs[i] + slopes[hid_obs[i],4]*ages3_obs[i] + alpha_hours*hours_obs[i] + alpha_date*date_obs[i] + eta_day[ dayid_obs[i] , 1 ];
            mugamma <- exp( slopes[hid_obs[i],5] + slopes[hid_obs[i],6]*ages_obs[i] + slopes[hid_obs[i],7]*ages2_obs[i] + slopes[hid_obs[i],8]*ages3_obs[i] + beta_hours*hours_obs[i] + beta_date*date_obs[i] + eta_day[ dayid_obs[i] , 2 ] );
            thetaJ <- exp( slopes[hid_obs[i],9] );
            dev <- dev + (-2) * if_else( iszero_obs[i] , log(inv_logit(pi)) , log1m(inv_logit(pi)) + gamma_log( nonzero_obs[i] , mugamma*thetaJ , thetaJ ) );
        }
        for ( i in 1:N_miss ) {
            pi <- slopes[hid_miss[i],1] + slopes[hid_miss[i],2]*ages_miss[i] + slopes[hid_miss[i],3]*ages2_miss[i] + slopes[hid_miss[i],4]*ages3_miss[i] + alpha_hours*hours_miss[i] + alpha_date*date_miss[i] + eta_day[ dayid_miss[i] , 1 ];
            mugamma <- exp( slopes[hid_miss[i],5] + slopes[hid_miss[i],6]*ages_miss[i] + slopes[hid_miss[i],7]*ages2_miss[i] + slopes[hid_miss[i],8]*ages3_miss[i] + beta_hours*hours_miss[i] + beta_date*date_miss[i] + eta_day[ dayid_miss[i] , 2 ] );
            thetaJ <- exp( slopes[hid_miss[i],9] );
            dev <- dev + (-2) * if_else( iszero_miss[i] , log(inv_logit(pi)) , log1m(inv_logit(pi)) + gamma_log( nonzero_miss[i] , mugamma*thetaJ , thetaJ ) );
        }
    }
'

# initialize
init_slopes <- matrix( 0 , nrow=147 , ncol=9 )
vec <- c( 0,0,0,0, log(7) ,0,0,0, log(0.3) )
for ( j in 1:147 ) init_slopes[j,] <- vec
initlist <- list(list(
    mu = vec ,
    Sigma = diag(9) ,
    slopes = init_slopes,
    alpha_hours = 0,
    beta_hours = 0,
    mu_hours = mean(d$hours,na.rm=TRUE),
    sigma_hours = sd(d$hours,na.rm=TRUE),
    hours_miss = rep( mean(d$hours,na.rm=TRUE) , N_miss ),
    alpha_date = 0,
    beta_date = 0,
    eta_day = matrix( 0 , nrow=Ndays , ncol=2 ),
    Sigma_day = diag(2)
))

fitVary9hDv <- stan( model_code = ache_code_Vary9hDv , data = ache_dat_h , init=initlist , iter = 7000 , warmup=2000 , chains = 1 , pars=c("mu","slopes","Sigma","alpha_date","beta_date","alpha_hours","beta_hours","mu_hours","sigma_hours","hours_miss","eta_day","Sigma_day","dev") )

#################################
# Alt_9hDv
# alternative age model, with hours and Julian date and varying effects on day

# split data into two parts: hours observed and hours missing
obs_idx <- which( !is.na( d$hours ) )
miss_idx <- which( is.na( d$hours ) )
N_obs <- length(obs_idx)
N_miss <- length(d$hours) - N_obs

# prep data
K <- 9
Ndays <- length( unique(d$date) )
dayid <- as.numeric( as.factor(d$date) )
ache_dat_h <- list(
    N = nrow(d),
    N_obs = N_obs,
    N_miss = N_miss,
    J = 147,
    K = K,
    y_obs = d$kg.meat[ obs_idx ],
    y_miss = d$kg.meat[ miss_idx ],
    iszero_obs = d$iszero[ obs_idx ],
    iszero_miss = d$iszero[ miss_idx ],
    nonzero_obs = d$nonzero[ obs_idx ],
    nonzero_miss = d$nonzero[ miss_idx ],
    hid_obs = d$hunter.id[ obs_idx ],
    hid_miss = d$hunter.id[ miss_idx ],
    ages_obs = d$age.s[ obs_idx ],
    ages2_obs = (d$age.s^2)[ obs_idx ],
    ages3_obs = (d$age.s^3)[ obs_idx ],
    ages_miss = d$age.s[ miss_idx ],
    ages2_miss = (d$age.s^2)[ miss_idx ],
    ages3_miss = (d$age.s^3)[ miss_idx ],
    hours_obs = d$hours[ obs_idx ],
    W = diag(K),
    date_obs = d$date.s[ obs_idx ],
    date_miss = d$date.s[ miss_idx ],
    Ndays = Ndays,
    dayid_obs = dayid[ obs_idx ],
    dayid_miss = dayid[ miss_idx ],
    mu_days = c(0,0)
)

# model code
ache_code_Alt9hDv <- '
  data {
    int<lower=0> N;                 // number of trips
    int<lower=0> N_obs;             // number of trips with observed hours
    int<lower=0> N_miss;            // number of trips with missing hours
    int<lower=0> J;                 // number of hunters
    int<lower=0> K;                 // number of varying effects
    real y_obs[N_obs];              // hunting returns
    real y_miss[N_miss];
    int<lower=0,upper=1> iszero_obs[N_obs]; // indicator of zero return
    int<lower=0,upper=1> iszero_miss[N_miss];
    real nonzero_obs[N_obs];        // positive returns
    real nonzero_miss[N_miss];
    int hid_obs[N_obs];             // hunter ids
    int hid_miss[N_miss];
    real ages_obs[N_obs];           // standardized age at time of hunt
    real ages2_obs[N_obs];          // square of age
    real ages3_obs[N_obs];          // cube of age
    real ages_miss[N_miss];
    real ages2_miss[N_miss];
    real ages3_miss[N_miss];
    matrix[K,K] W;
    real<lower=0> hours_obs[N_obs]; // observed hours hunted values
    real date_obs[N_obs];
    real date_miss[N_miss];
    int Ndays;                      // number of unique days in sample
    int dayid_obs[N_obs];           // day ids corresponding to cases with observed hours
    int dayid_miss[N_miss];         // day ids corresponding to cases with MISSING hours
    vector[2] mu_days;              // c(0,0) location vector for day varying effects
  } 
  parameters {
    vector[K] mu;
    cov_matrix[K] Sigma;
    vector[K] slopes[J];
    real<lower=0> hours_miss[N_miss]; // missing values to impute
    real alpha_hours;
    real beta_hours;
    real mu_hours; // imputation populaton
    real<lower=0> sigma_hours;
    real alpha_date;
    real beta_date;
    vector[2] eta_day[Ndays];
    cov_matrix[2] Sigma_day;
  }
  model {
    real pi;
    real mugamma;
    real thetaJ;
    real coef; // temp for calculating age coefficient
    Sigma ~ inv_wishart( K+1 , W );
    mu_hours ~ normal( 7 , 100 );
    sigma_hours ~ uniform( 0, 20 );
    alpha_date ~ normal( 0 , 100 );
    beta_date ~ normal( 0 , 100 );
    for ( k in 1:K ) mu[k] ~ normal(0,100);
    for( j in 1:J ) slopes[j] ~ multi_normal( mu , Sigma );
    for( j in 1:Ndays ) eta_day[j] ~ multi_normal( mu_days , Sigma_day );
    for ( i in 1:N_miss ) {
        hours_miss[i] ~ normal( mu_hours , sigma_hours );
        // use int_step() to ask if age > midpoint
        coef <- if_else( int_step(ages_miss[i]-slopes[hid_miss[i],2]) , slopes[hid_miss[i],4] , slopes[hid_miss[i],3] );
        pi <- slopes[hid_miss[i],1] + coef*pow(ages_miss[i]-slopes[hid_miss[i],2],2) + alpha_hours*hours_miss[i] + alpha_date*date_miss[i] + eta_day[ dayid_miss[i] , 1 ];
        coef <- if_else( int_step(ages_miss[i]-slopes[hid_miss[i],6]) , slopes[hid_miss[i],8] , slopes[hid_miss[i],7] );
        mugamma <- exp( slopes[hid_miss[i],5] + coef*pow(ages_miss[i]-slopes[hid_miss[i],6],2) + beta_hours*hours_miss[i] + beta_date*date_miss[i] + eta_day[ dayid_miss[i] , 2 ] );
        thetaJ <- exp( slopes[hid_miss[i],9] );
        iszero_miss[i] ~ bernoulli( inv_logit(pi) );
        for ( g in 1:(1-iszero_miss[i]) )
            nonzero_miss[i] ~ gamma( mugamma*thetaJ , thetaJ );
        //lp__ <- lp__ + if_else( iszero_miss[i] , log(inv_logit(pi)) , log1m(inv_logit(pi)) + gamma_log( nonzero_miss[i] , mugamma*thetaJ , thetaJ ) );
    }
    for( i in 1:N_obs ) {
        hours_obs[i] ~ normal( mu_hours , sigma_hours );
        coef <- if_else( int_step(ages_obs[i]-slopes[hid_obs[i],2]) , slopes[hid_obs[i],4] , slopes[hid_obs[i],3] );
        pi <- slopes[hid_obs[i],1] + coef*pow(ages_obs[i]-slopes[hid_obs[i],2],2) + alpha_hours*hours_obs[i] + alpha_date*date_obs[i] + eta_day[ dayid_obs[i] , 1 ];
        coef <- if_else( int_step(ages_obs[i]-slopes[hid_obs[i],6]) , slopes[hid_obs[i],8] , slopes[hid_obs[i],7] );
        mugamma <- exp( slopes[hid_obs[i],5] + coef*pow(ages_obs[i]-slopes[hid_obs[i],6],2) + beta_hours*hours_obs[i] + beta_date*date_obs[i] + eta_day[ dayid_obs[i] , 2 ] );
        thetaJ <- exp( slopes[hid_obs[i],9] );
        iszero_obs[i] ~ bernoulli( inv_logit(pi) );
        for ( g in 1:(1-iszero_obs[i]) )
            nonzero_obs[i] ~ gamma( mugamma*thetaJ , thetaJ );
        //lp__ <- lp__ + if_else( iszero_obs[i] , log(inv_logit(pi)) , log1m(inv_logit(pi)) + gamma_log( nonzero_obs[i] , mugamma*thetaJ , thetaJ ) );
    }
  }

  generated quantities {
        real dev;
        real pi;
        real thetaJ;
        real mugamma;
        dev <- 0;
        for( i in 1:N_obs ) {
            real coef;
            coef <- if_else( int_step(ages_obs[i]-slopes[hid_obs[i],2]) , slopes[hid_obs[i],4] , slopes[hid_obs[i],3] );
            pi <- slopes[hid_obs[i],1] + coef*pow(ages_obs[i]-slopes[hid_obs[i],2],2) + alpha_hours*hours_obs[i] + alpha_date*date_obs[i] + eta_day[ dayid_obs[i] , 1 ];
            coef <- if_else( int_step(ages_obs[i]-slopes[hid_obs[i],6]) , slopes[hid_obs[i],8] , slopes[hid_obs[i],7] );
            mugamma <- exp( slopes[hid_obs[i],5] + coef*pow(ages_obs[i]-slopes[hid_obs[i],6],2) + beta_hours*hours_obs[i] + beta_date*date_obs[i] + eta_day[ dayid_obs[i] , 2 ] );
            thetaJ <- exp( slopes[hid_obs[i],9] );
            dev <- dev + (-2) * if_else( iszero_obs[i] , log(inv_logit(pi)) , log1m(inv_logit(pi)) + gamma_log( nonzero_obs[i] , mugamma*thetaJ , thetaJ ) );
        }
        for ( i in 1:N_miss ) {
            real coef;
            coef <- if_else( int_step(ages_miss[i]-slopes[hid_miss[i],2]) , slopes[hid_miss[i],4] , slopes[hid_miss[i],3] );
            pi <- slopes[hid_miss[i],1] + coef*pow(ages_miss[i]-slopes[hid_miss[i],2],2) + alpha_hours*hours_miss[i] + alpha_date*date_miss[i] + eta_day[ dayid_miss[i] , 1 ];
            coef <- if_else( int_step(ages_miss[i]-slopes[hid_miss[i],6]) , slopes[hid_miss[i],8] , slopes[hid_miss[i],7] );
            mugamma <- exp( slopes[hid_miss[i],5] + coef*pow(ages_miss[i]-slopes[hid_miss[i],6],2) + beta_hours*hours_miss[i] + beta_date*date_miss[i] + eta_day[ dayid_miss[i] , 2 ] );
            thetaJ <- exp( slopes[hid_miss[i],9] );
            dev <- dev + (-2) * if_else( iszero_miss[i] , log(inv_logit(pi)) , log1m(inv_logit(pi)) + gamma_log( nonzero_miss[i] , mugamma*thetaJ , thetaJ ) );
        }
    }
'

# initialize
init_slopes <- matrix( 0 , nrow=147 , ncol=9 )
vec <- c( 0,0,0,0, log(7) ,0,0,0, log(0.3) )
for ( j in 1:147 ) init_slopes[j,] <- vec
initlist <- list(list(
    mu = vec ,
    Sigma = diag(9) ,
    slopes = init_slopes,
    alpha_hours = 0,
    beta_hours = 0,
    mu_hours = mean(d$hours,na.rm=TRUE),
    sigma_hours = sd(d$hours,na.rm=TRUE),
    hours_miss = rep( mean(d$hours,na.rm=TRUE) , N_miss ),
    alpha_date = 0,
    beta_date = 0,
    eta_day = matrix( 0 , nrow=Ndays , ncol=2 ),
    Sigma_day = diag(2)
))

fitAlt9hDv <- stan( model_code = ache_code_Alt9hDv , data = ache_dat_h , init=initlist , iter = 7000 , warmup=2000 , chains = 1 , pars=c("mu","slopes","Sigma","alpha_date","beta_date","alpha_hours","beta_hours","mu_hours","sigma_hours","hours_miss","eta_day","Sigma_day","dev") )

