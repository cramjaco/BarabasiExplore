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

library(lubridate)

dateVec <- seq(from = as_date("1983-March-20"), by = "month", length.out = 1000)

sin(c(1:3 * pi))

make_microbe_series <- function(n, base = rgamma(1, 5), base_weight = 1,
                                ac_weight = 100,
                                neighbor_vec = NULL, lag = 0, neighbor_weight = 0,
                                season_lag = runif(1, 0, 11), season_weight = 0,
                                sd = 0.1){
  mts = rep(NA, length.out = n)
  if(is.null(neighbor_vec)){
    neighbor_vec = rep(0, n)
    neighbor_weight = 0
  }
  
  season_vec = sin((season_lag + 0:n-1) * 2 * pi / 12)
  
  if(!(lag %in% c(0,1))){stop("lag must be 0 or 1")}
  ## first element only
  mts[1] = (
    base * base_weight + # base value
      neighbor_vec[1] * neighbor_weight # somewhere between the value of base and the value of the connected vector -- we ignore the time lag for this step
    ) /(base_weight + neighbor_weight) + # normalize to weights
    season_vec[1] * season_weight + 
    rnorm(1, 0, sd) # add random noise
  if(mts[1] < 0){mts[1] = 0} # if the answer is less than zero, make it zero
  
  ## second through last element
  for(iter in 2:n){
    mts[iter] = (
      ## each value is the sum of several components
      base * base_weight + #base is the starting value, but also a value towards wich things trend 
        mts[iter - 1] * ac_weight+ # autocorrelation influence, the values tend to be similar to the values before and after 
        neighbor_vec[iter - lag] * neighbor_weight + # the value of the neighboring vector
        0 # placeholder, can delete
    )/(base_weight + ac_weight + neighbor_weight) + # divide by weights
      season_vec[iter] * season_weight + # this one goes outside of the expression above because its added to everything else
      rnorm(1, 0, sd) # finally, add some random noise on top of everything with a specified sd = standard deviation
    if(mts[iter] < 0){mts[iter] = 0} # if the value is below zero, just make it zero
  }
  #mts[mts < 0] <- 0
  mts
}

## Second version

make_microbe_proto_series <- function(n,
                                ac_weight = .1,
                                neighbor_vec = NULL, lag = 0, neighbor_weight = 0,
                                season_lag = runif(1, 0, 11), season_weight = 0,
                                trend_weight = 0,
                                sd = 0.1){
  mts = rep(NA, length.out = n)
  if(is.null(neighbor_vec)){
    neighbor_vec = rep(0, n)
    neighbor_weight = 0
  }
  
  season_vec = sin((season_lag + 0:n-1) * 2 * pi / 12)
  trend_vec = seq(from = -1, to = 1, length.out = n)
  
  if(!(lag %in% c(0,1))){stop("lag must be 0 or 1")}
  ## first element only
  mts[1] = (
    0 +
      season_vec[1] * season_weight +
      trend_vec[1] * trend_weight +
      neighbor_vec[1] * neighbor_weight +
      rnorm(1, 0, sd)
)
  
  ## second through last element
  for(iter in 2:n){
    mts[iter] =
      0 +
      season_vec[iter] * season_weight + 
      trend_vec[iter] * trend_weight +
      neighbor_vec[iter] * neighbor_weight +
        mts[iter - 1] * ac_weight +
        rnorm(1, 0, sd)
  }
  #mts[mts < 0] <- 0
  tibble(step = 1:n, mts = mts)
}

protoseries_to_counts <- function(protoSeries, base = rgamma(1, shape = 1)){
  mts <- protoSeries$mts
  
  mts_rebase <- mts + base
  mts_exp <- exp(mts_rebase)
  
  ct = map_int(mts_exp, ~rpois(n = 1, lambda = .))
  
  tibble(step = protoSeries$step, ct = ct)
  
}

# I'm tempted to run this to make the pseudo network file. Then do some post processing where I add on a "base" value, exp() tranform everthing and then pull an rpois of that.
# The rpois being the random draw from the underlyind distribution.

# Step 1 -- Make Network
library(igraph)
seed(33100)

NNodes <- 20
theNetwork <- sample_pa(n = NNodes, algorithm = "psumtree-multiple")

theNetwork_df <- igraph::as_data_frame(bag) %>% rename(to = "from", from = "to") %>% select(from, to) %>%
  mutate(lag = rbinom(length(from), 1, 0))

theNetwork01 <-graph_from_data_frame(bag_df) 

plot(theNetwork01)


# Step 2, generate pseudo series for each node

pseudo_df <- tibble(node = 1:NNodes, series = vector(mode = "list", length = NNodes))

loc_network_df <- theNetwork_df

pseudo_df <- tibble(node = 1:NNodes, series = vector(mode = "list", length = NNodes))
count_df <- tibble(node = 1:NNodes, series = vector(mode = "list", length = NNodes))
for(nodeNum in 1:NNodes){
  neighborData <- NULL
  # Is node in the network table's "to" column?
  if(nodeNum %in% loc_network_df$to){
    # if so, save the neigbor's value to neighbor
    neighborNode <- loc_network_df$from[which(loc_network_df$to == nodeNum)]
    neighborData = loc_network_df[["series"]][[neighborNode]]
  }
  loc_protoseries <- make_microbe_proto_series(100,
                            ac_weight = rbeta(1, 1.5, 1.5) * .4 - .4/2,
                            season_weight = rbeta(1, 1.5, 1.5) * .4 - .4/2,
                            trend_weight = rbeta(1, 1.5, 1.5) * .2 - .2/2, sd = 1,
                            neighbor_vec = neighborData, neighbor_weight = 0.5) 
  loc_count <- protoseries_to_counts(loc_protoseries, base = rgamma(1, shape = 1))
  
  pseudo_df[["series"]][[nodeNum]] <- loc_protoseries
  count_df[["series"]][[nodeNum]] <- loc_count
}

## 2b cat to data frame

# Step 3, 

count_df_2 <- unnest(count_df)

count_df_wide <- pivot_wider(count_df_2, names_from = "step", values_from = "ct")



### Graveyard/Scratch


set.seed(seed)
everyonesNeighbor <- make_microbe_proto_series(100, ac_weight = .6, season_weight = .2, trend_weight = .1, sd = 1) 
everyonesNeighbor%>% plot(type = "l")

# additional challenges: unevenly spaced vectors

##Lets explore this parameter space...

paramSpace <- expand_grid(
  base = c(0, 0.1, 0.5, 1, 2, 5),
  ac_weight = c(0, 1, 10, 100),
  neighbor_weight = c(0, 1, 10, 100),
  season_weight = c(0, 0.001, 0.01, 0.1, 1),
  sd = c(0, 0.001, 0.1, 1, 5, 10),
  rep = 1:10
)


mms_wrap <- function(base, ac_weight, neighbor_weight, season_weight, sd, rep){
  series <- make_microbe_series(n = 100, base = base, base_weight = 1, 
                      ac_weight = ac_weight,
                      neighbor_vec = everyonesNeighbor, lag = 0, neighbor_weight = neighbor_weight,
                      season_lag = 0, season_weight = season_weight,
                      sd = sd)
  
  tibble(n = 1:100, series = series)
}


paramSeries <- paramSpace %>%
  mutate(series = pmap(., mms_wrap))

paramSeries <- paramSeries %>% unnest(cols = c(series))

paramSeries %>% filter(base %in% c(0, 1, 5),neighbor_weight == 0, season_weight == 0, sd == 1, rep <=5) %>%
  ggplot(aes(x = n, y = series, group = rep)) + geom_path() +
  facet_grid(base ~ac_weight , scales = "free")

paramSeries %>% filter(base %in% c(5), neighbor_weight == 0, sd == 1, rep <=5) %>%
  ggplot(aes(x = n, y = series, group = rep)) + geom_path() +
  facet_grid(ac_weight ~season_weight )

wnp <- paramSeries %>% filter(base %in% c(5), season_weight == 0.1, sd == .1, rep <=5) %>%
  ggplot(aes(x = n, y = series, group = rep)) + geom_path() +
  facet_grid(ac_weight~neighbor_weight, scales = "free") 

enp <- tibble(n = 1:100, series = everyonesNeighbor) %>% ggplot(aes(x = n, y = series)) + geom_path() + theme(plot.margin = unit(c(2, 0, 2, 0), "in"))

cowplot::plot_grid(enp, wnp, rel_widths = c(1,2))

