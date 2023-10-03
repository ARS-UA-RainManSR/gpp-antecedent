# Average vwc across houses
library(tidyverse)

vwc <- read_csv("data/VWCtest.csv")

VWC_mean <- vwc |> 
  group_by(Date, PT, Dep) |> 
  summarize(vwc_mean = mean(VWC))

# split into depths 1 and 2, make wide for each trt

vwc12 <- VWC_mean |> 
  filter(Dep == 1,
         PT != "S3") |> 
  select(-Dep) |> 
  pivot_wider(names_from = PT, values_from = vwc_mean)

vwc25 <- VWC_mean |> 
  filter(Dep == 2) |> 
  select(-Dep) |> 
  pivot_wider(names_from = PT, values_from = vwc_mean)

write_csv(vwc12, file = "data/VWC12.csv")
write_csv(vwc25, file = "data/VWC25.csv")
