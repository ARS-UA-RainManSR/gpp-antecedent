# Prep data for SAM model
library(tidyverse)
library(rjags)
library(mcmcplots)
library(postjags)

gpp <- read_csv("data/GPPtest.csv") %>%
  filter(PT != "S3") %>% # few S3s, not enough matching vwc
  mutate(treat = factor(PT, levels = c("S1", "S2","S4")),
         trt = as.numeric(treat),
         doy = as.numeric(Date - as.Date("2021-01-01") + 1))
vpd <- read_csv("data/VPDtest.csv")
vwc12 <- read_csv("data/VWC12.csv")
vwc25 <- read_csv("data/VWC25.csv")

dat_list <- list(
  GPP = gpp$GPP,
  N = nrow(gpp),
  trt = gpp$trt,
  PAR = as.vector(scale(gpp$Par)),
  doy = gpp$doy,
  Dmean = as.vector(scale(vpd$VPD)),
  W12 = matrix(as.vector(scale(vwc12[,2:4])), ncol = 3),
  W25 = matrix(as.vector(scale(vwc25[,2:4])), ncol = 3),
  nlagA = 7,
  nlagB = 5, 
  nlagC = 7, 
  pA = 1,
  pB = 1,
  pC = 1,
  alphaA = rep(1, 5),
  alphaB = rep(1, 10),
  alphaC = rep(1, 7),
  NParam = 9,
  NTrt = length(unique(gpp$PT))
)
str(dat_list)

# Create initials 
# Function to generate initials
init <- function() {
  list(B = matrix(rnorm(dat_list$NParam * dat_list$NTrt, 0, 10), 
                  nrow = dat_list$NParam),
       tau = runif(1, 0, 1),
       deltaA = runif(dat_list$nlagA, 0, 1),
       deltaB = runif(dat_list$nlagB, 0, 1),
       deltaC = runif(dat_list$nlagC, 0, 1)
  )
}

inits_list <- list(init(), init(), init())

# Initialize model
jm <- jags.model("scripts/model-v2/PAR_Dant_Want.jags",
                  data = dat_list,
                  inits = inits_list,
                  n.chains = 3)

update(jm, n.iter = 10000)
dic <- dic.samples(jm, n.iter = 3000, thin = 5)

# Monitor
params <- c("deviance", "Dsum", "R2_resid",
             "B",
             "tau", "sig",
             "wA", "wB", "wC",
             "deltaA", "deltaB", "deltaC")

coda.out <- coda.samples(jm,
                         variable.names = params,
                         n.iter = 3000,
                         n.thin = 5)

save(coda.out, file = "scripts/model-v2/coda/coda.Rdata")

# Inspect chains visually
mcmcplot(coda.out, parms = c("deviance", "Dsum", "R2_resid",
                             "B",
                             "tau", "sig",
                             "wA", "wB", "wC"))
caterplot(coda.out, regex = "^B\\[[1-9],1", reorder = FALSE)
caterplot(coda.out, regex = "^B\\[[1-9],2", reorder = FALSE)
caterplot(coda.out, regex = "^B\\[[1-9],3", reorder = FALSE)
caterplot(coda.out, regex = "^wA\\[", reorder = FALSE)
caterplot(coda.out, regex = "^wB\\[", reorder = FALSE)
caterplot(coda.out, regex = "^wC\\[", reorder = FALSE)

# Save state
final <- initfind(coda.out, OpenBUGS = FALSE)
final[[1]]
saved_state <- removevars(final, variables = c(2:3, 7, 9:11))
saved_state[[1]]

save(saved_state, file = "scripts/model-v2/inits/saved_state.Rdata")

