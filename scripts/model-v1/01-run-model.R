# Prep data for SAM model
library(tidyverse)
library(rjags)
library(mcmcplots)
library(postjags)

gpp <- read_csv("data/GPPtest.csv") %>%
  filter(PT != "S3") %>%
  mutate(treat = factor(PT, levels = c("S1", "S2","S4")),
         trt = as.numeric(treat),
         doy = as.numeric(Date - as.Date("2021-01-01") + 1))
vpd <- read_csv("data/VPDtest.csv")
# vwc <- read_csv("data/VWCtest.csv")

dat_list <- list(
  GPP = gpp$GPP,
  N = nrow(gpp),
  trt = gpp$trt,
  PAR = as.vector(scale(gpp$Par)),
  doy = gpp$doy,
  Dmean = as.vector(scale(vpd$VPD)),
  nlagA = 10,
  pA = 1,
  alphaA = rep(1, 10),
  NParam = 4,
  NTrt = length(unique(gpp$PT))
)

# Create initials 
# Function to generate initials
init <- function() {
  list(B = matrix(rnorm(dat_list$NParam * dat_list$NTrt, 0, 10), 
                  nrow = dat_list$NParam),
       tau = runif(1, 0, 1),
       deltaA = runif(dat_list$nlagA, 0, 1)
  )
}

inits_list <- list(init(), init(), init())

# Initialize model
jm <- jags.model("scripts/model-v1/PAR_Dant.jags",
                  data = dat_list,
                  inits = inits_list,
                  n.chains = 3)

update(jm, n.iter = 10000)
dic.samples(jm, n.iter = 3000, thin = 5)

# Monitor
params <- c("deviance", "Dsum", "R2_resid",
             "B",
             "tau", "sig",
             "wA",
             "deltaA")

coda.out <- coda.samples(jm,
                         variable.names = params,
                         n.iter = 3000,
                         n.thin = 5)

save(coda.out, file = "scripts/model-v1/coda/coda.Rdata")

# Inspect chains visually
mcmcplot(coda.out, parms = c("deviance", "Dsum", "R2_resid",
                             "B",
                             "tau", "sig",
                             "wA"))
caterplot(coda.out, regex = "^B\\[", reorder = FALSE)
caterplot(coda.out, regex = "^wA", reorder = FALSE)
# caterplot(coda.out, regex = "^wB\\[", reorder = FALSE)

# Save state
final <- initfind(coda.out, OpenBUGS = FALSE)
final[[1]]
saved_state <- removevars(final, variables = c(2:3, 5, 7))
saved_state[[1]]

save(saved_state, file = "scripts/model-v1/inits/saved_state.Rdata")

