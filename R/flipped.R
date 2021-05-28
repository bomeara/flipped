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
dcoin_linear <- function(nheads, nflips, preflip_prob=0.5, slope=0.1, log=FALSE, outside_bounds_is_NA=FALSE) {
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
prob_heads_exponential_decay <- function(nflips, halflife) {
  lambda <- log(2)/halflife
  probabilities <- 1*exp(-lambda*seq(from=1, to=nflips, by=1)) #by not starting at exp(-lambda * 0), it allows for something that has no heads on the first flip to have positive probability
  return(probabilities)
}

#' Compute the probability of heads with each flip given a multiplier model
#' The model assumes 50% chance of heads before a coin is picked up and it changes as a percentage of the previous value each flip. i.e., the probability of heads is 101% of the probability the previous flip with a multiplier of 1.01. 
#' @param nflips Total number of flips (heads and tails)
#' @param multiplier Factor to multiply the previous probability by
#' @param outside_bounds_is_NA If TRUE, if any probability of heads is outside the bounds of probability, the function returns NA. Otherwise, it sets the value to the nearer bound.
#' @return Vector of probability of heads for the first flip, second flip, etc.
#' @export
prob_heads_multiplicative <- function(nflips, multiplier, outside_bounds_is_NA=FALSE) {
  probabilities <- 0.5*(multiplier^(sequence(nflips)))
	if(any(probabilities>1 | probabilities<0)) {
		if(outside_bounds_is_NA) {
			return(NA)
		} else {
			probabilities[probabilities>1] <- 1
			probabilities[probabilities<0] <- 0
		}
	}
  return(probabilities)
}



#' Compute probability of observations given an exponential decay model
#' @param nheads Number of heads
#' @param nflips Total number of flips (heads and tails)
#' @param multiplier How much to multiply by each flip
#' @param log If TRUE return log transformed probabilities.
#' @param possibilities All possible sequences of flips that lead to the observed number of heads
#' @param outside_bounds_is_NA If TRUE, if any probability of heads is outside the bounds of probability, the function returns NA. Otherwise, it sets the value to the nearer bound.
#' @return The likelihood of the data (or log likelihood if log=TRUE)
#' @export
d_coin_multiplicative <- function(nheads, nflips, multiplier, log=FALSE, possibilities=get_possibilities(nheads,nflips), outside_bounds_is_NA=FALSE) {
	pheads <- prob_heads_multiplicative(nflips=nflips, multiplier=multiplier, outside_bounds_is_NA=outside_bounds_is_NA)
	return(dcoin_from_probability(pheads=pheads, nheads=nheads, nflips=nflips, log=log, possibilities=possibilities))
}



#' Exhaustively get all possible sets of outcomes that result in a specified number of heads out of a certain number of flips
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

#' Compute probability of observations given an exponential decay model
#' The idea is that the coin before handling has 100% chance of heads, but each time it is picked up that probability will decrease (maybe it is bent by the statistician's mighty thumb). After halflife times handling it, the probability of heads is 50%, and it keeps dropping from there.
#' @param nheads Number of heads
#' @param nflips Total number of flips (heads and tails)
#' @param halflife How many flips to get to 50% heads
#' @param log If TRUE return log transformed probabilities.
#' @param possibilities All possible sequences of flips that lead to the observed number of heads
#' @return The likelihood of the data (or log likelihood if log=TRUE)
#' @export
dcoin_exponential_decay <- function(nheads, nflips, halflife, log=FALSE, possibilities=get_possibilities(nheads,nflips)) {
  pheads <- prob_heads_exponential_decay(nflips, halflife)
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
  #print(c(halflife, likelihood))
  return(ifelse(log, log(likelihood), likelihood))
}

#' Computes the likelihood for a range of values using a linear coin model 
#' @inheritParams prob_heads_linear
#' @param param_range Range of parameters to try
#' @param number_of_steps How many values of the parameter to try
#' @export
#' @return vector of likelihoods
#' @examples
#' nheads <- 8
#' nflips <- 10
#' linear_results <- profile_linear_model(nheads, nflips)
#' plot(x=linear_results$preflip_prob, y=linear_results$likelihood, type="l")
#' dbinom_proportions <- seq(from=0, to=1, length.out=1000)
#' lines(dbinom_proportions, dbinom(nheads, nflips, dbinom_proportions), col="red")
#' best_param <- linear_results$preflip_prob[which.max(linear_results$likelihood, na.rm=TRUE)]
#' print(best_param)
profile_linear_model <- function(nheads, nflips, param_range=c(0,1), slope=0.1, number_of_steps=1000,log=FALSE, outside_bounds_is_NA=FALSE) {
	starting_prob_heads <- seq(from=min(param_range), to=max(param_range), length.out=number_of_steps)
	likelihoods <- rep(NA, number_of_steps)
	for (i in sequence(number_of_steps)) {
		likelihoods[i] <- dcoin_linear(nheads=nheads, nflips=nflips, preflip_prob=starting_prob_heads[i], slope=slope, log=log, outside_bounds_is_NA=outside_bounds_is_NA)
	}
	return(data.frame(preflip_prob=starting_prob_heads, likelihood=likelihoods))
}

#' Computes the likelihood for a range of values using an exponential coin model 
#' @inheritParams prob_heads_linear
#' @param param_range Range of parameters to try
#' @param number_of_steps How many values of the parameter to try
#' @export
#' @return vector of likelihoods
#' @examples
#' nheads <- 8
#' nflips <- 10
#' exp_results <- profile_exponential_decay_model(nheads, nflips)
#' plot(x=exp_results$preflip_prob, y=exp_results$likelihood, type="l")
#' best_param <- exp_results$halflife[which.max(exp_results$likelihood, na.rm=TRUE)]
#' print(best_param)
profile_exponential_decay_model <- function(nheads, nflips, param_range=c(0,nflips*10), number_of_steps=1000,log=FALSE) {
	halflives <- seq(from=min(param_range), to=max(param_range), length.out=number_of_steps)
	likelihoods <- rep(NA, number_of_steps)
	for (i in sequence(number_of_steps)) {
		likelihoods[i] <- dcoin_exponential_decay(nheads=nheads, nflips=nflips, halflife=halflives[i], log=log)
	}
	return(data.frame(halflife=halflives, likelihood=likelihoods))
}



#' Compute probability of observations given a vector of probability of heads
#' @param pheads Vector with the probability of a heads on flip 1, 2, etc.
#' @param nheads Number of heads
#' @param nflips Total number of flips (heads and tails)
#' @param log If TRUE return log transformed probabilities.
#' @param possibilities All possible sequences of flips that lead to the observed number of heads
#' @param diff_value If not NULL, the final likelihood will be abs(likelihood - diff_value) for minimizing a function
#' @return The likelihood of the data (or log likelihood if log=TRUE)
#' @export
dcoin_from_probability <- function(pheads, nheads, nflips, log=FALSE, possibilities=get_possibilities(nheads,nflips), diff_value=NULL) {
  if(length(pheads)!=nflips) {
	stop("Length of pheads should equal nflips")  
  }
  if(any(pheads<0) | any(pheads>1)) {
	if(!is.null(diff_value)) {
		return(1e6) # b/c we're optimizing, most likely
	} else {
		return(NA)	
	}
  }
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
  #print(c(halflife, likelihood))
  final_value <- ifelse(log, log(likelihood), likelihood)
  if(!is.null(diff_value)) {
	final_value <- abs(final_value-diff_value)  
  }
  return(final_value)
}

#' Compute probability of observations across many potential vectors
#' This will try (1/stepsize)^nflips possible vectors, computing the probability of the observation for each
#' @param nheads Number of heads
#' @param nflips Total number of flips (heads and tails)
#' @param number_samples How many vectors to sample
#' @param stopval How large a difference in probability is considered close enough between the flat model and others
#' @param log If TRUE return log transformed probabilities.
#' @param possibilities All possible sequences of flips that lead to the observed number of heads
#' @return The likelihood of the data (or log likelihood if log=TRUE)
#' @export
try_many_vectors <- function(nheads, nflips, number_samples=1000, stopval=0.00001, log=FALSE, possibilities=get_possibilities(nheads,nflips)) {
	old_probs <- rep(nheads/nflips, nflips)
	old_lik <- stats::dbinom(nheads, nflips, nheads/nflips)
	starting_lik <- old_lik
	results <- data.frame(matrix(nrow=number_samples+1, ncol=1+nflips))
	results[1,] <- c(old_lik, old_probs)
	print(results[1,])
	# for (i in sequence(number_samples)) {
	# 	new_probs <- old_probs * runif(nflips, min=0.9, max=1.1)
	# 	new_probs[new_probs>1] <- 1
	# 	new_lik <- dcoin_from_probability(nheads, nflips, pheads=new_probs, possibilities=possibilities)
	# 	results[i+1,] <- c(new_lik, new_probs)
	# 	if((new_lik-starting_lik)>0) {
	# 		old_probs <- new_probs	
	# 	}
	# 	if((starting_lik-new_lik)<0.01 & runif(1)<0.3) {
	# 		old_probs <- new_probs	
	# 	}
	# 	print(results[i+1,])
	# }
	for (i in sequence(number_samples)) {
		difference <- Inf
		while(difference>stopval) {
			start_prob <- sort(runif(nflips, min=0, max=1), decreasing=(runif(1,0,1)<0.5))
			optimized <- nloptr::nloptr(start_prob, dcoin_from_probability, lb=rep(0, nflips), ub=rep(1, nflips), nheads=nheads, nflips=nflips, possibilities=possibilities, diff_value=starting_lik, log=FALSE, opts=list("algorithm"="NLOPT_LN_NEWUOA", "stopval"=stopval))
			difference <- optimized$objective
		}
		likelihood <- dcoin_from_probability(optimized$solution, nheads=nheads, nflips=nflips, possibilities=possibilities, diff_value=NULL, log=FALSE)
		results[i+1,] <- c(likelihood, optimized$solution)
		print(results[i+1,])
	}
	colnames(results) <- c("likelihood_of_data", paste0("flip_", sequence(nflips)))
	return(results)
}


#' Find congruent models to a simple binomial model
#' This will find the parameter values for other models that equal the likelihood for a simple binomial model. This may not be the MLE for these other models
#' @param nheads Number of heads
#' @param nflips Total number of flips (heads and tails)
#' @param slopes Vector of slopes to use
#' @param stopval How large a difference in probability is considered close enough between the flat model and others
#' @return A list containing the parameter estimates with likelihoods for each model and the probabilities for heads at each model
#' @export
find_congruent_models <- function(nheads, nflips, slopes=c(0, 0.1, -0.05), stopval=0.0001) {
	old_probs <- rep(nheads/nflips, nflips)
	old_lik <- stats::dbinom(nheads, nflips, nheads/nflips)
	starting_lik <- old_lik
	probabilities_per_flip <- data.frame(matrix(nrow=0, ncol=3))
	colnames(probabilities_per_flip) <- c("Flip", "Model", "Probability_of_heads_this_flip")
	multiplicative_model <- function(x, nheads, nflips,  diff_value) {
		return(abs(diff_value - d_coin_multiplicative(nheads=nheads, nflips=nflips, multiplier=x)))
	}
	
	slope_model <- function(x, nheads, nflips,  diff_value, slope) {
		return(abs(diff_value - dcoin_linear(nheads=nheads, nflips=nflips, preflip_prob=x, slope=slope)))
	}
	
	exponential_model <- function(x, nheads, nflips,  diff_value) {
		return(abs(diff_value - dcoin_exponential_decay(nheads=nheads, nflips=nflips, halflife=x)))
	}	
	
	formatted_results <- data.frame(matrix(nrow=0, ncol=4))
	colnames(formatted_results) <- c("likelihood_of_data", "model", "param_name", "param_value")
	
	difference <- Inf
	while(difference>stopval) {
		print("Optimizing multiplicative model")
		optimized <- nloptr::nloptr(.99, multiplicative_model, nheads=nheads, nflips=nflips, diff_value=starting_lik, opts=list("algorithm"="NLOPT_LN_NELDERMEAD", "stopval"=stopval))
		difference <- optimized$objective
		print(paste0("Got within ", difference, " of the target"))
	}
	formatted_results[1,] <- c(multiplicative_model(optimized$solution, nheads=nheads, nflips=nflips, diff_value=0), "Multiplicative", "Multiplier", optimized$solution)
	probabilities_per_flip <- rbind(probabilities_per_flip, data.frame(Flips=seq(from=0, to=nflips, by=1), Model=rep("Multiplier", nflips+1), Probability_of_heads_this_flip=c(0.5,prob_heads_multiplicative(nflips=nflips, multiplier=optimized$solution ))))
	
	difference <- Inf
	while(difference>stopval) {
		print("Optimizing exponential decay model")
		optimized <- nloptr::nloptr(nflips/2, exponential_model, nheads=nheads, nflips=nflips, diff_value=starting_lik, opts=list("algorithm"="NLOPT_LN_NELDERMEAD", "stopval"=stopval))
		difference <- optimized$objective
		print(paste0("Got within ", difference, " of the target"))

	}
	formatted_results[2,] <- c(exponential_model(optimized$solution, nheads=nheads, nflips=nflips, diff_value=0), "Exponential Decay", "Halflife", optimized$solution)
	probabilities_per_flip <- rbind(probabilities_per_flip, data.frame(Flips=seq(from=0, to=nflips, by=1), Model=rep("Exponential Decay", nflips+1), Probability_of_heads_this_flip=c(1,prob_heads_exponential_decay(nflips=nflips, halflife=optimized$solution ))))

	for(i in seq_along(slopes)) {
		difference <- Inf
		while(difference>stopval) {
			print(paste0("Optimizing slope ", slopes[i], " model"))
			optimized <- nloptr::nloptr(nheads/nflips, slope_model, nheads=nheads, nflips=nflips, diff_value=starting_lik, slope=slopes[i], opts=list("algorithm"="NLOPT_LN_NELDERMEAD", "stopval"=stopval))
			difference <- optimized$objective
			print(paste0("Got within ", difference, " of the target"))
		}
		formatted_results[nrow(formatted_results)+1,] <- c(slope_model(optimized$solution, nheads=nheads, nflips=nflips, diff_value=0, slope=slopes[i]), paste0("Slope of ", slopes[i]), "Preflip Probability of Heads", optimized$solution)
		probabilities_per_flip <- rbind(probabilities_per_flip, data.frame(Flips=seq(from=0, to=nflips, by=1), Model=rep(paste0("Slope of ", slopes[i]), nflips+1), Probability_of_heads_this_flip=c(optimized$solution,prob_heads_linear(nflips=nflips, preflip_prob=optimized$solution, slope=slopes[i] ))))
	}
	return(list(estimates=formatted_results, probabilities_per_flip=probabilities_per_flip))
}