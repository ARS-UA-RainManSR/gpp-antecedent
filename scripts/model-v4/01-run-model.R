# Prep data for SAM model
library(tidyverse)
library(rjags)
library(mcmcplots)
library(postjags)
library(broom.mixed)

gpp <- read_csv("data/GPPtest.csv") %>%
  filter(PT != "S3") %>% # few S3s, not enough matching vwc
  mutate(treat = factor(PT, levels = c("S1", "S2","S4")),
         trt = as.numeric(treat),
         doy = as.numeric(Date - as.Date("2021-01-01") + 1),
         house = factor(block, levels = c("H1", "H2", "H3", "H4", "H5")))
vpd <- read_csv("data/VPDtest.csv")
vwc12 <- read_csv("data/VWC12.csv")
vwc25 <- read_csv("data/VWC25.csv")

dat_list <- list(
  GPP = gpp$GPP,
  N = nrow(gpp),
  trt = gpp$trt,
  block = gpp$house,
  PAR = as.vector(scale(gpp$Par)),
  doy = gpp$doy,
  Dmean = as.vector(scale(vpd$VPD)),
  W12 = matrix(as.vector(scale(vwc12[,2:4])), ncol = 3),
  W25 = matrix(as.vector(scale(vwc25[,2:4])), ncol = 3),
  nlagA = 7,
  nlagB = 7, 
  nlagC = 7, 
  pA = 1,
  pB = 1,
  pC = 1,
  alphaA = rep(1, 7),
  alphaB = matrix(rep(1, 3*7), nrow = 3),
  alphaC = matrix(rep(1, 3*7), nrow = 3),
  NParam = 7,
  NTrt = length(unique(gpp$PT)),
  NBlock = length(unique(gpp$house))
)
str(dat_list)

# Create initials 
# Function to generate initials
init <- function() {
  list(B = matrix(rnorm(dat_list$NParam * dat_list$NTrt, 0, 10), 
                  nrow = dat_list$NParam),
       tau = runif(1, 0, 1),
       deltaA = runif(dat_list$nlagA, 0, 1),
       deltaB = matrix(runif(dat_list$nlagB * dat_list$NTrt, 0, 1), nrow = dat_list$NTrt),
       deltaC = matrix(runif(dat_list$nlagC * dat_list$NTrt, 0, 1), nrow = dat_list$NTrt),
       sig.eps = runif(1, 0, 1)
  )
}

inits_list <- list(init(), init(), init())

# Initialize model
jm <- jags.model("scripts/model-v4/PAR_Dant_Want_bytrt_RE.jags",
                  data = dat_list,
                  inits = inits_list,
                  n.chains = 3)

update(jm, n.iter = 10000)
dic <- dic.samples(jm, n.iter = 3000, thin = 5)

# Monitor
params <- c("deviance", "Dsum", "R2_resid",
             "B", "Bstar", "Estar",
             "tau", "sig",
             "wA", "wB", "wC",
             "deltaA", "deltaB", "deltaC",
            "sig.eps")

coda.out <- coda.samples(jm,
                         variable.names = params,
                         n.iter = 3000,
                         n.thin = 5)

save(coda.out, file = "scripts/model-v4/coda/coda.Rdata")

# Inspect chains visually
mcmcplot(coda.out, parms = c("deviance", "Dsum", "R2_resid",
                             "Bstar", "Estar",
                             "tau", "sig",
                             "wA", "wB", "wC"))
caterplot(coda.out, regex = "^Bstar\\[[1-7],1", reorder = FALSE)
caterplot(coda.out, regex = "^Bstar\\[[1-7],2", reorder = FALSE)
caterplot(coda.out, regex = "^Bstar\\[[1-7],3", reorder = FALSE)
caterplot(coda.out, regex = "^Estar\\[[1-5],1", reorder = FALSE)
caterplot(coda.out, regex = "^Estar\\[[1-5],2", reorder = FALSE)
caterplot(coda.out, regex = "^Estar\\[[1-5],3", reorder = FALSE)
caterplot(coda.out, regex = "^wA\\[", reorder = FALSE)
caterplot(coda.out, regex = "^wB\\[1", reorder = FALSE)
caterplot(coda.out, regex = "^wB\\[2", reorder = FALSE)
caterplot(coda.out, regex = "^wB\\[3", reorder = FALSE)
caterplot(coda.out, regex = "^wC\\[1", reorder = FALSE)
caterplot(coda.out, regex = "^wC\\[2", reorder = FALSE)
caterplot(coda.out, regex = "^wC\\[3", reorder = FALSE)

# Check convergence diagnostic
gel <- gelman.diag(coda.out, multivariate = FALSE)
gel$psrf %>%
  data.frame() %>%
  tibble::rownames_to_column() %>%
  filter(grepl("Dsum", rowname) |grepl("deviance", rowname))

gel$psrf %>%
  data.frame() %>%
  tibble::rownames_to_column() %>%
  filter(grepl("sig", rowname))

gel$psrf %>%
  data.frame() %>%
  tibble::rownames_to_column() %>%
  filter(grepl("^Bstar", rowname))

gel$psrf %>%
  data.frame() %>%
  tibble::rownames_to_column() %>%
  filter(grepl("^Estar", rowname))

gel$psrf %>%
  data.frame() %>%
  tibble::rownames_to_column() %>%
  filter(grepl("^wA", rowname))

gel$psrf %>%
  data.frame() %>%
  tibble::rownames_to_column() %>%
  filter(grepl("^wB", rowname))

gel$psrf %>%
  data.frame() %>%
  tibble::rownames_to_column() %>%
  filter(grepl("^wC", rowname))

# Save state
final <- initfind(coda.out, OpenBUGS = FALSE)
final[[1]]
saved_state <- removevars(final, variables = c(2:5, 9, 12:14))
saved_state[[1]]

save(saved_state, file = "scripts/model-v4/inits/saved_state.Rdata")

# Run for replicated data

# Monitor
coda.rep <- coda.samples(jm,
                         variable.names = "GPP.rep",
                         n.iter = 3000,
                         n.thin = 5)

save(coda.rep, file = "scripts/model-v4/coda/coda.Rdata")

# Summarize replicated output
coda_sum <- tidyMCMC(coda.rep,
                     conf.int = TRUE,
                     conf.method = "HPDinterval") %>%
  rename(pred.mean = estimate,
         pred.lower = conf.low,
         pred.upper = conf.high)


# Check model fit
pred <- cbind.data.frame(gpp, coda_sum)

m_all <- lm(pred.mean ~ GPP, data = pred)
summary(m_all) # R2 = 0.6174

m1 <- lm(pred.mean ~ GPP, data = filter(pred, trt == 1))
summary(m1) # R2 = 0.5882

m2 <- lm(pred.mean ~ GPP, data = filter(pred, trt == 2))
summary(m2) # R2 = 0.6859

m3 <- lm(pred.mean ~ GPP, data = filter(pred, trt == 3))
summary(m3) # R2 = 0.5423

pred %>%
  ggplot(aes(x = GPP, y =pred.mean,
             color = factor(trt))) +
  geom_abline(intercept = 0, slope = 1, col = "red") +
  geom_errorbar(aes(ymin = pred.lower, ymax = pred.upper),
                alpha = 0.25) +
  geom_point() +
  scale_x_continuous("Observed") +
  scale_y_continuous("Predicted") +
  theme_bw(base_size = 12) +
  facet_wrap(~trt) +
  coord_equal()

