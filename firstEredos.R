library(igraph)
library(tidyverse)

seed = 33100

irg <- erdos.renyi.game(n = 10, p.or.m = .2, type = "gnp", directed = TRUE)

plot(irg)

degree_distribution(irg)

bag <- sample_pa(n = 20, algorithm = "psumtree-multiple")

plot(bag)

get.edge.attribute(bag)

as_data_frame(bag)

bag_df <- igraph::as_data_frame(bag) %>% rename(to = "from", from = "to") %>% select(from, to) %>%
  mutate(lag = rbinom(length(from), 1, 0.5))

bag_02 <-graph_from_data_frame(bag_df) 

plot(bag_02)


## Lets say I want a data frame, with random, non-autocorrelated true abundances.
## Each variable has a mean which draws from a gamma distribution. 
# Then we add a normal distribution on top of that.

make_microbe_series <- function(n, base = rgamma(1, 5), base_weight = 1, ac_weight = 100, neighbor_weight = 100, neighbor_vec = NULL, lag = 0, sd = 0.1){
  mts = rep(NA, length.out = n)
  mts[1] = base + rnorm(1, 0, 1)
  if(mts[1] < 0){mts[1] = 0}
  for(iter in 2:n){
    mts[iter] = (100 * mts[iter - 1] + 1 * base)/101 + rnorm(1, 0, sd)
    if(mts[iter] < 0){mts[iter] = 0}
    }
  #mts[mts < 0] <- 0
  mts
}

make_microbe_series(100, 5) %>% plot(type = "l")
