#' Compute probability of observations given linear change model
#' This is essentially stats::dbinom() but allowing for the probability of heads to linearly change from the starting value. By default it increases by 10% per flip, but this can be set to other values. 
#' 
#' The idea is that the coin before handling has some probability of heads, but each time it is picked up that probability could change (maybe it is bent by the statistician's mighty thumb). The slope gives the amount of change in this probability each flip: for example, a coin that starts fair and which has a slope of 0.01 has a probability of heads of 0.51 (0.50 + 0.01) on its first flip, 0.52 on its second, and so forth. If the act of flipping has absolutely no effect on the probability of heads, slope can be set to be zero, though using stats::dbinom() for this particular edge case should be faster.
#' 
#' Of course, if all we have is the total number of heads and total number of flips, we do not know if it was HTT, THT, or TTH. For the particular case of a slope set to exactly zero the order does not matter, but in the general case it will. For example, if the probability of heads increases with each flip, HTT is less likely than TTH even though each has one heads out of three flips. The current code looks at all possibilities exhaustively, but more efficient ways to calculate this undoubtedly exist. Pull requests are welcome. It also means this may be slow as the number of flips increases.
#' 
#' For some slopes and preflip probabilities,the probabilities of heads on a given flip may be outside the 0 to 1 bounds. By default, if this happens the function returns NA. If outside_bounds_is_NA is FALSE, it moves the probabilities to the nearer bound. 
#' @param nheads Number of heads
#' @inheritParams prob_heads_linear
#' @param log If TRUE return log transformed probabilities.
#' @param outside_bounds_is_NA If TRUE, if any probability of heads is outside the bounds of probability, the function returns NA. Otherwise, it sets the value to the nearer bound.
#' @return The likelihood of the data (or log likelihood if log=TRUE)
#' @export
dcoin_linear <- function(nheads, nflips, preflip_prob=0.5, slope=0.1, log=FALSE, outside_bounds_is_NA=TRUE) {
    probabilities <- prob_heads_linear (nflips=nflips, preflip_prob=preflip_prob, slope=slope)
  if(any(probabilities>1 | probabilities<0)) {
      if(outside_bounds_is_NA) {
          return(NA)
      } else {
          probabilities[probabilities>1] <- 1
          probabilities[probabilities<0] <- 0
      }
  }
  totalheadsprobs <- rep(0, nflips+1) #so it goes from 0 to nflips total heads
  for (i in seq_along(probabilities)){
    probhead <- probabilities[i]
    probtail <- 1-probhead
    if(i==1) {
      totalheadsprobs[0+1] <- probtail
      totalheadsprobs[1+1] <- probhead
    } else {
      totalheadsprobs <- totalheadsprobs*probtail + c(0,totalheadsprobs[1:(length(totalheadsprobs)-1)]*probhead)
    }
  }
  return(ifelse(log, log(totalheadsprobs[nheads+1]), totalheadsprobs[nheads+1]))
}

#' Compute the probability of heads with each flip given a linear change model.
#' @param nflips Total number of flips (heads and tails)
#' @param preflip_prob Probability of heads before the coin is handled
#' @param slope How much the probability changes each time the coin is flipped
#' @return Vector of probability of heads for the first flip, second flip, etc.
#' @export
prob_heads_linear <- function(nflips, preflip_prob=0.5, slope=0.1) {
  probabilities <- slope*sequence(nflips)/nflips + preflip_prob
  return(probabilities)
}

#' Compute the probability of heads with each flip given an exponential model
#' The model assumes 100% chance of heads before a coin is picked up and it drops exponentially each time the coin is handled.
#' @param nflips Total number of flips (heads and tails)
#' @param halflife How many flips to get to 50% heads
#' @return Vector of probability of heads for the first flip, second flip, etc.
#' @export
prob_heads_exponential <- function(nflips, halflife) {
  lambda <- log(2)/halflife
  probabilities <- 1*exp(-lambda*seq(from=1, to=nflips, by=1)) #by not starting at exp(-lambda * 0), it allows for something that has no heads on the first flip to have positive probability
  return(probabilities)
}

#' Exhauatively get all possible sets of outcomes that result in a specified number of heads out of a certain number of flips
#' 
#' This grows very large with the number of flips. It will throw an error if you try too many flips.
#' @param nheads Number of heads
#' @param nflips Total number of flips (heads and tails)
#' @return data.frame with each potential trial as a row. 1=heads, 0=tails.
#' @export
get_possibilities <- function(nheads, nflips) {
 if(nflips>25) {
    stop("Too many flips for the expand.grid silly approach to work")
 } 
  possibilities <- expand.grid(rep(list(c(0,1)),nflips))
  possibilities <- possibilities[rowSums(possibilities)==nheads,]
  return(possibilities)
}

#' Compute probability of observations given an exponential model
#' The idea is that the coin before handling has 100% chance of heads, but each time it is picked up that probability will decrease (maybe it is bent by the statistician's mighty thumb). After halflife times handling it, the probability of heads is 50%, and it keeps dropping from there.
#' @param nheads Number of heads
#' @param nflips Total number of flips (heads and tails)
#' @param halflife How many flips to get to 50% heads
#' @param log If TRUE return log transformed probabilities.
#' @param possibilities All possible sequences of flips that lead to the observed number of heads
#' @return The likelihood of the data (or log likelihood if log=TRUE)
#' @export
dcoin_exponential <- function(nheads, nflips, halflife, log=FALSE, possibilities=get_possibilities(nheads,nflips)) {
  pheads <- prob_heads_exponential(nflips, halflife)
  ptails <- 1-pheads
  prob_matrix <- rbind(ptails, pheads) # heads = 1, tails =0, so when we offset, tails is first row, heads second
  
  likelihood <- 0
  for (i in sequence(nrow(possibilities))) {
    positions <- as.numeric(possibilities[i,]+1)
    local_prob <- 1
    for (j in sequence(ncol(possibilities))) {
      local_prob <- local_prob * unname(prob_matrix[positions[j],j])
    }
    likelihood <- likelihood + local_prob
  }
  print(c(halflife, likelihood))
  return(ifelse(log, log(likelihood), likelihood))
}

