
# Source files from:
# https://github.com/avehtari/rhat_ess/blob/master/code/monitornew.R
# https://github.com/avehtari/rhat_ess/blob/master/code/monitorplot.R

fft_next_good_size <- function(N) {
	# Find the optimal next size for the FFT so that
	# a minimum number of zeros are padded.
	if (N <= 2)
		return(2)
	while (TRUE) {
		m = N
		while ((m %% 2) == 0) m = m / 2
		while ((m %% 3) == 0) m = m / 3
		while ((m %% 5) == 0) m = m / 5
		if (m <= 1)
			return(N)
		N = N + 1
	}
}

autocovariance <- function(y) {
	# Compute autocovariance estimates for every lag for the specified
	# input sequence using a fast Fourier transform approach.
	N <- length(y)
	M <- fft_next_good_size(N)
	Mt2 <- 2 * M
	yc <- y - mean(y)
	yc <- c(yc, rep.int(0, Mt2 - N))
	transform <- fft(yc)
	ac <- fft(Conj(transform) * transform, inverse = TRUE)
	ac <- Re(ac)[1:N] / (N * 2 * seq(N, 1, by = -1))
	ac
}

autocorrelation <- function(y) {
	# Compute autocorrelation estimates for every lag for the specified
	# input sequence using a fast Fourier transform approach.
	ac <- autocovariance(y)
	ac <- ac / ac[1]
}

z_scale <- function(x) {
	S <- length(x)
	r <- rank(x, ties.method = 'average')
	z <- qnorm((r - 1 / 2) / S)
	if (!is.null(dim(x))) {
		# output should have the input dimension
		z <- array(z, dim = dim(x), dimnames = dimnames(x))
	}
	z
}

u_scale <- function(x) {
	S <- length(x)
	r <- rank(x, ties.method = 'average')
	u <- (r - 1 / 2) / S
	if (!is.null(dim(x))) {
		# output should have the input dimension
		u <- array(u, dim = dim(x), dimnames = dimnames(x))
	}
	u
}

r_scale <- function(x) {
	S <- length(x)
	r <- rank(x, ties.method = 'average')
	if (!is.null(dim(x))) {
		# output should have the input dimension
		r <- array(r, dim = dim(x), dimnames = dimnames(x))
	}
	r
}

ess_rfun <- function(sims) {
	# Compute the effective sample size for samples of several chains 
	# for one parameter; see the C++ code of function  
	# effective_sample_size in chains.cpp 
	# Args:
	#   sims: a 2D array _without_ warmup samples (# iter * # chains) 
	# Returns:
	#   A single numeric value
	if (is.vector(sims)) {
		dim(sims) <- c(length(sims), 1)
	}
	chains <- ncol(sims)
	n_samples <- nrow(sims)
	
	acov <- lapply(seq_len(chains), function(i) autocovariance(sims[, i])) 
	acov <- do.call(cbind, acov)
	chain_mean <- apply(sims, 2, mean)
	mean_var <- mean(acov[1, ]) * n_samples / (n_samples - 1) 
	var_plus <- mean_var * (n_samples - 1) / n_samples
	if (chains > 1) 
		var_plus <- var_plus + var(chain_mean)
	
	# Geyer's initial positive sequence
	rho_hat_t <- rep.int(0, n_samples)
	t <- 0
	rho_hat_even <- 1;
	rho_hat_t[t + 1] <- rho_hat_even
	rho_hat_odd <- 1 - (mean_var - mean(acov[t + 2, ])) / var_plus
	rho_hat_t[t + 2] <- rho_hat_odd
	t <- 2  
	while (t < nrow(acov) - 1 && !is.nan(rho_hat_even + rho_hat_odd) && 
		   (rho_hat_even + rho_hat_odd > 0)) {
		rho_hat_even = 1 - (mean_var - mean(acov[t + 1, ])) / var_plus
		rho_hat_odd = 1 - (mean_var - mean(acov[t + 2, ])) / var_plus
		if ((rho_hat_even + rho_hat_odd) >= 0) {
			rho_hat_t[t + 1] <- rho_hat_even
			rho_hat_t[t + 2] <- rho_hat_odd
		}
		t <- t + 2
	}
	max_t <- t
	# Geyer's initial monotone sequence
	t <- 2
	while (t <= max_t - 2) {  
		if (rho_hat_t[t + 1] + rho_hat_t[t + 2] >
			rho_hat_t[t - 1] + rho_hat_t[t]) {
			rho_hat_t[t + 1] = (rho_hat_t[t - 1] + rho_hat_t[t]) / 2;
			rho_hat_t[t + 2] = rho_hat_t[t + 1];
		}
		t <- t + 2
	}
	ess <- chains * n_samples
	ess <- ess / (-1 + 2 * sum(rho_hat_t[1:max_t]))
	ess 
} 

tail_ess <- function(sims) {
	I05 <- sims <= quantile(sims, 0.05)
	q05_ess <- ess_rfun(z_scale(split_chains(I05)))
	I95 <- sims <= quantile(sims, 0.95)
	q95_ess <- ess_rfun(z_scale(split_chains(I95)))
	min(q05_ess, q95_ess)
}

rhat_rfun <- function(sims) {
	# Compute the rhat convergence diagnostic for a single parameter
	# For split-rhat, just call this with splitted chains
	# Args:
	#   sims: a 2D array _without_ warmup samples (# iter * # chains) 
	# Returns:
	#   A single numeric value
	if (is.vector(sims)) {
		dim(sims) <- c(length(sims), 1)
	}
	chains <- ncol(sims)
	n_samples <- nrow(sims)
	chain_mean <- numeric(chains)
	chain_var <- numeric(chains)
	for (i in seq_len(chains)) {
		chain_mean[i] <- mean(sims[, i])
		chain_var[i] <- var(sims[, i])
	} 
	var_between <- n_samples * var(chain_mean)
	var_within <- mean(chain_var) 
	sqrt((var_between/var_within + n_samples - 1) / n_samples)
} 

quantile_mcse <- function(sims, prob) {
	# compute Markov-chain SE of quantiles for a single parameter
	# prob must be a single probability for the quantile of interest
	if (is.vector(sims)) {
		dim(sims) <- c(length(sims), 1)
	}
	I <- sims <= quantile(sims, prob)
	Seff <- ess_rfun(z_scale(split_chains(I)))
	p <- c(0.1586553, 0.8413447, 0.05, 0.95)
	a <- qbeta(p, Seff * prob + 1, Seff * (1 - prob) + 1)
	ssims <- sort(sims)
	S <- length(ssims)
	th1 <- ssims[max(round(a[1] * S), 1)]
	th2 <- ssims[min(round(a[2] * S), S)]
	mcse <- (th2 - th1) / 2
	th1 <- ssims[max(round(a[3] * S), 1)]
	th2 <- ssims[min(round(a[4] * S), S)]
	data.frame(mcse = mcse, Q05 = th1, Q95 = th2, Seff = Seff)
}

split_chains <- function(sims) {
	# split Markov chains
	# Args:
	#   sims: a 2D array of samples (# iter * # chains) 
	if (is.vector(sims)) {
		dim(sims) <- c(length(sims), 1)
	}
	niter <- dim(sims)[1]
	half <- niter / 2
	cbind(sims[1:floor(half), ], sims[ceiling(half + 1):niter, ])
}

is_constant <- function(x, tol = .Machine$double.eps) {
	abs(max(x) - min(x)) < tol
}

monitor <- function(sims, warmup = 0, probs = c(0.05, 0.50, 0.95)) { 
	# print a summary for general simulation results 
	# of 3D array: # iter * # chains * # parameters 
	# Args:
	#   sims: a 3D array described above 
	#   warmup: the number of iterations used for warmup 
	#   probs: probs of summarizing quantiles 
	# Return: 
	#   A summary data.frame of class 'simsummary'
	if (inherits(sims, "stanfit")) {
		chains <- sims@sim$chains
		iter <- sims@sim$iter
		warmup <- sims@sim$warmup
		parnames <- names(sims)
		sims <- as.array(sims)
	} else {
		dim_sims <- dim(sims)
		if (is.null(dim_sims)) {
			dim(sims) <- c(length(sims), 1, 1) 
		} else if (length(dim_sims) == 2) {
			dim(sims) <- c(dim_sims, 1)
		} else if (length(dim_sims) > 3) {
			stop("'sims' has more than 3 dimensions") 
		}
		parnames <- dimnames(sims)[[3]]
		if (is.null(parnames)) {
			parnames <- paste0("V", seq_len(dim(sims)[3]))
		}
		iter <- dim(sims)[1]
		chains <- dim(sims)[2]
		if (warmup > dim(sims)[1]) {
			stop("warmup is larger than the total number of iterations")
		}
		if (warmup >= 1) {
			sims <- sims[-seq_len(warmup), , , drop = FALSE] 
		}
	}
	
	mcse_quan_fun <- function(p, sims) quantile_mcse(sims, p)$mcse
	out <- vector("list", length(parnames))
	out <- setNames(out, parnames)
	for (i in seq_along(out)) {
		sims_i <- sims[, , i]
		valid <- all(is.finite(sims_i))
		quan <- unname(quantile(sims_i, probs = probs))
		mcse_quan <- sapply(probs, mcse_quan_fun, sims_i)
		
		zsims_split <- z_scale(split_chains(sims_i))
		zsplit_rhat <- rhat_rfun(zsims_split)
		bulk_ess <- round(ess_rfun(zsims_split))
		tail_ess <- round(tail_ess(sims_i))
		
		sims_folded <- abs(sims_i - median(sims_i))
		zsims_folded_split <- z_scale(split_chains(sims_folded))
		zfsplit_rhat <- rhat_rfun(zsims_folded_split)
		rhat <- max(zsplit_rhat, zfsplit_rhat)
		
		ess <- ess_rfun(sims_i)
		mean <- mean(sims_i)
		sd <- sd(sims_i)
		mcse_mean <- sd / sqrt(ess)
		# mcse_sd assumes normality of sims and uses Stirling's approximation
		# min of ess for sims and sims^2 is used instead of ess
		ess2 <- ess_rfun(sims_i^2)
		essmin <- min(c(ess,ess2))
		fac_mcse_sd <- sqrt(exp(1) * (1 - 1 / essmin)^(essmin - 1) - 1)
		mcse_sd <- sd * fac_mcse_sd
		
		out[[i]] <- c(
			valid, quan, mean, sd, mcse_quan, mcse_mean,
			mcse_sd, rhat, bulk_ess, tail_ess
		)
	}
	
	out <- as.data.frame(do.call(rbind, out))
	str_quan <- paste0("Q", probs * 100)
	str_mcse_quan <- paste0("MCSE_", str_quan)
	colnames(out) <- c(
		"valid", str_quan, "Mean", "SD", str_mcse_quan, 
		"MCSE_Mean", "MCSE_SD", "Rhat", "Bulk_ESS", "Tail_ESS"
	)
	rownames(out) <- parnames
	
	# replace NAs with appropriate values if draws are valid
	S <- prod(dim(sims)[1:2]) 
	out$Rhat[out$valid & !is.finite(out$Rhat)] <- 1
	out$Bulk_ESS[out$valid & !is.finite(out$Bulk_ESS)] <- S
	out$Tail_ESS[out$valid & !is.finite(out$Tail_ESS)] <- S
	SE_vars <- colnames(out)[grepl("^SE_", colnames(out))]
	for (v in SE_vars) {
		out[[v]][out$valid & !is.finite(out[[v]])] <- 0
	}
	out$valid <- NULL
	
	structure(
		out,
		chains = chains,
		iter = iter,
		warmup = warmup,
		class = c("simsummary", "data.frame") 
	)
} 

monitor_extra <- function(sims, warmup = 0, probs = c(0.05, 0.50, 0.95)) { 
	# print an extended summary for general simulation results 
	# of 3D array: # iter * # chains * # parameters 
	# Args:
	#   sims: a 3D array described above 
	#   warmup: the number of iterations used for warmup 
	#   probs: probs of summarizing quantiles 
	# Return: 
	#   A summary data.frame of class 'simsummary'
	if (inherits(sims, "stanfit")) {
		chains <- sims@sim$chains
		iter <- sims@sim$iter
		warmup <- sims@sim$warmup
		parnames <- names(sims)
		sims <- as.array(sims)
	} else {
		dim_sims <- dim(sims)
		if (is.null(dim_sims)) {
			dim(sims) <- c(length(sims), 1, 1) 
		} else if (length(dim_sims) == 2) {
			dim(sims) <- c(dim_sims, 1)
		} else if (length(dim_sims) > 3) {
			stop("'sims' has more than 3 dimensions") 
		}
		parnames <- dimnames(sims)[[3]]
		if (is.null(parnames)) {
			parnames <- paste0("V", seq_len(dim(sims)[3]))
		}
		iter <- dim(sims)[1]
		chains <- dim(sims)[2]
		if (warmup > dim(sims)[1]) {
			stop("warmup is larger than the total number of iterations")
		}
		if (warmup >= 1) {
			sims <- sims[-seq_len(warmup), , , drop = FALSE] 
		}
	}
	
	mcse_fun <- function(p, sims) quantile_mcse(sims, p)$mcse
	out <- vector("list", length(parnames))
	out <- setNames(out, parnames)
	for (i in seq_along(out)) {
		sims_i <- sims[, , i]
		nsamples <- prod(dim(sims_i))
		mean <- mean(sims_i)
		sd <- sd(sims_i)
		quan <- unname(quantile(sims_i, probs = probs))
		ess <- round(ess_rfun(sims_i))
		ress <- ess / nsamples
		sem <- sd / sqrt(ess)
		rhat <- rhat_rfun(sims_i)
		
		split_ess <- round(ess_rfun(split_chains(sims_i)))
		split_rhat <- rhat_rfun(split_chains(sims_i))
		zess <- round(ess_rfun(z_scale(sims_i)))
		zrhat <- rhat_rfun(z_scale(sims_i))
		
		zsims_split <- z_scale(split_chains(sims_i))
		zsplit_rhat <- rhat_rfun(zsims_split)
		zsplit_ess <- round(ess_rfun(zsims_split))
		zsplit_ress <- zsplit_ess / nsamples
		
		sims_centered <- sims_i - median(sims_i)
		sims_folded <- abs(sims_centered)
		zsims_folded_split <- z_scale(split_chains(sims_folded))
		zfsplit_rhat <- rhat_rfun(zsims_folded_split)
		zfsplit_ess <- round(ess_rfun(zsims_folded_split))
		zfsplit_ress <- zfsplit_ess / nsamples
		
		tail_ess <- round(tail_ess(sims_i))
		tail_ress <- tail_ess / nsamples 
		
		sims_med <- (sims_centered <= 0) * 1
		sims_mad <- ((sims_folded - median(sims_folded)) <= 0) * 1
		medsplit_ess <- round(ess_rfun(z_scale(split_chains(sims_med))))
		medsplit_ress <- medsplit_ess / nsamples 
		madsplit_ess <- round(ess_rfun(z_scale(split_chains(sims_mad))))
		madsplit_ress <- madsplit_ess / nsamples 
		
		out[[i]] <- c(
			mean, sem, sd, quan, ess, ress, split_ess, zess, zsplit_ess, zsplit_ress, 
			rhat, split_rhat, zrhat, zsplit_rhat, zfsplit_rhat, zfsplit_ess, 
			zfsplit_ress, tail_ess, tail_ress, medsplit_ess, medsplit_ress, 
			madsplit_ess, madsplit_ress
		)
	}
	
	out <- as.data.frame(do.call(rbind, out))
	probs_str <- paste0("Q", probs * 100)
	colnames(out) <- c(
		"mean", "se_mean", "sd", probs_str, "seff", "reff", "sseff", "zseff", 
		"zsseff", "zsreff", "Rhat", "sRhat", "zRhat", "zsRhat", "zfsRhat", 
		"zfsseff", "zfsreff", "tailseff", "tailreff", "medsseff", "medsreff", 
		"madsseff", "madsreff"
	)
	rownames(out) <- parnames
	structure(
		out,
		chains = chains,
		iter = iter,
		warmup = warmup,
		extra = TRUE,
		class = c("simsummary", "data.frame") 
	)
}

print.simsummary <- function(x, digits = 3, se = FALSE, ...) {
	atts <- attributes(x)
	px <- x
	if (!se) {
		px <- px[, !grepl("^MCSE_", colnames(px))]
	}
	class(px) <- "data.frame"
	decimal_places <- max(1, digits - 1)
	px$Rhat <- round(px$Rhat, digits = max(2, decimal_places))
	estimates <- setdiff(names(px), c("Rhat", "Bulk_ESS", "Tail_ESS"))
	px[, estimates] <- round(px[, estimates], digits = decimal_places)
	#px[, estimates] <- signif(px[, estimates], digits = digits)
	# add a space between summary and convergence estimates
	names(px)[names(px) %in% "Rhat"] <- " Rhat"
	cat(
		"Inference for the input samples (", atts$chains, 
		" chains: each with iter = ", atts$iter, 
		"; warmup = ", atts$warmup, "):\n\n", sep = ""
	)
	print(px, ...)
	if (!isTRUE(atts$extra)) {
		cat(
			"\nFor each parameter, Bulk_ESS and Tail_ESS are crude measures of \n",
			"effective sample size for bulk and tail quantities respectively (good values is \n",
			"ESS > 400), and Rhat is the potential scale reduction factor on rank normalized\n",
			"split chains (at convergence, Rhat = 1).\n", sep = ""
		)
	}
	invisible(x)
}


fft_next_good_size <- function(N) {
	# Find the optimal next size for the FFT so that
	# a minimum number of zeros are padded.
	if (N <= 2)
		return(2)
	while (TRUE) {
		m = N
		while ((m %% 2) == 0) m = m / 2
		while ((m %% 3) == 0) m = m / 3
		while ((m %% 5) == 0) m = m / 5
		if (m <= 1)
			return(N)
		N = N + 1
	}
}

autocovariance <- function(y) {
	# Compute autocovariance estimates for every lag for the specified
	# input sequence using a fast Fourier transform approach.
	N <- length(y)
	M <- fft_next_good_size(N)
	Mt2 <- 2 * M
	yc <- y - mean(y)
	yc <- c(yc, rep.int(0, Mt2 - N))
	transform <- fft(yc)
	ac <- fft(Conj(transform) * transform, inverse = TRUE)
	ac <- Re(ac)[1:N] / (N * 2 * seq(N, 1, by = -1))
	ac
}

autocorrelation <- function(y) {
	# Compute autocorrelation estimates for every lag for the specified
	# input sequence using a fast Fourier transform approach.
	ac <- autocovariance(y)
	ac <- ac / ac[1]
}

z_scale <- function(x) {
	S <- length(x)
	r <- rank(x, ties.method = 'average')
	z <- qnorm((r - 1 / 2) / S)
	if (!is.null(dim(x))) {
		# output should have the input dimension
		z <- array(z, dim = dim(x), dimnames = dimnames(x))
	}
	z
}

u_scale <- function(x) {
	S <- length(x)
	r <- rank(x, ties.method = 'average')
	u <- (r - 1 / 2) / S
	if (!is.null(dim(x))) {
		# output should have the input dimension
		u <- array(u, dim = dim(x), dimnames = dimnames(x))
	}
	u
}

r_scale <- function(x) {
	S <- length(x)
	r <- rank(x, ties.method = 'average')
	if (!is.null(dim(x))) {
		# output should have the input dimension
		r <- array(r, dim = dim(x), dimnames = dimnames(x))
	}
	r
}

ess_rfun <- function(sims) {
	# Compute the effective sample size for samples of several chains 
	# for one parameter; see the C++ code of function  
	# effective_sample_size in chains.cpp 
	# Args:
	#   sims: a 2D array _without_ warmup samples (# iter * # chains) 
	# Returns:
	#   A single numeric value
	if (is.vector(sims)) {
		dim(sims) <- c(length(sims), 1)
	}
	chains <- ncol(sims)
	n_samples <- nrow(sims)
	
	acov <- lapply(seq_len(chains), function(i) autocovariance(sims[, i])) 
	acov <- do.call(cbind, acov)
	chain_mean <- apply(sims, 2, mean)
	mean_var <- mean(acov[1, ]) * n_samples / (n_samples - 1) 
	var_plus <- mean_var * (n_samples - 1) / n_samples
	if (chains > 1) 
		var_plus <- var_plus + var(chain_mean)
	
	# Geyer's initial positive sequence
	rho_hat_t <- rep.int(0, n_samples)
	t <- 0
	rho_hat_even <- 1;
	rho_hat_t[t + 1] <- rho_hat_even
	rho_hat_odd <- 1 - (mean_var - mean(acov[t + 2, ])) / var_plus
	rho_hat_t[t + 2] <- rho_hat_odd
	t <- 2  
	while (t < nrow(acov) - 1 && !is.nan(rho_hat_even + rho_hat_odd) && 
		   (rho_hat_even + rho_hat_odd > 0)) {
		rho_hat_even = 1 - (mean_var - mean(acov[t + 1, ])) / var_plus
		rho_hat_odd = 1 - (mean_var - mean(acov[t + 2, ])) / var_plus
		if ((rho_hat_even + rho_hat_odd) >= 0) {
			rho_hat_t[t + 1] <- rho_hat_even
			rho_hat_t[t + 2] <- rho_hat_odd
		}
		t <- t + 2
	}
	max_t <- t
	# Geyer's initial monotone sequence
	t <- 2
	while (t <= max_t - 2) {  
		if (rho_hat_t[t + 1] + rho_hat_t[t + 2] >
			rho_hat_t[t - 1] + rho_hat_t[t]) {
			rho_hat_t[t + 1] = (rho_hat_t[t - 1] + rho_hat_t[t]) / 2;
			rho_hat_t[t + 2] = rho_hat_t[t + 1];
		}
		t <- t + 2
	}
	ess <- chains * n_samples
	ess <- ess / (-1 + 2 * sum(rho_hat_t[1:max_t]))
	ess 
} 

tail_ess <- function(sims) {
	I05 <- sims <= quantile(sims, 0.05)
	q05_ess <- ess_rfun(z_scale(split_chains(I05)))
	I95 <- sims <= quantile(sims, 0.95)
	q95_ess <- ess_rfun(z_scale(split_chains(I95)))
	min(q05_ess, q95_ess)
}

rhat_rfun <- function(sims) {
	# Compute the rhat convergence diagnostic for a single parameter
	# For split-rhat, just call this with splitted chains
	# Args:
	#   sims: a 2D array _without_ warmup samples (# iter * # chains) 
	# Returns:
	#   A single numeric value
	if (is.vector(sims)) {
		dim(sims) <- c(length(sims), 1)
	}
	chains <- ncol(sims)
	n_samples <- nrow(sims)
	chain_mean <- numeric(chains)
	chain_var <- numeric(chains)
	for (i in seq_len(chains)) {
		chain_mean[i] <- mean(sims[, i])
		chain_var[i] <- var(sims[, i])
	} 
	var_between <- n_samples * var(chain_mean)
	var_within <- mean(chain_var) 
	sqrt((var_between/var_within + n_samples - 1) / n_samples)
} 

quantile_mcse <- function(sims, prob) {
	# compute Markov-chain SE of quantiles for a single parameter
	# prob must be a single probability for the quantile of interest
	if (is.vector(sims)) {
		dim(sims) <- c(length(sims), 1)
	}
	I <- sims <= quantile(sims, prob)
	Seff <- ess_rfun(z_scale(split_chains(I)))
	p <- c(0.1586553, 0.8413447, 0.05, 0.95)
	a <- qbeta(p, Seff * prob + 1, Seff * (1 - prob) + 1)
	ssims <- sort(sims)
	S <- length(ssims)
	th1 <- ssims[max(round(a[1] * S), 1)]
	th2 <- ssims[min(round(a[2] * S), S)]
	mcse <- (th2 - th1) / 2
	th1 <- ssims[max(round(a[3] * S), 1)]
	th2 <- ssims[min(round(a[4] * S), S)]
	data.frame(mcse = mcse, Q05 = th1, Q95 = th2, Seff = Seff)
}

split_chains <- function(sims) {
	# split Markov chains
	# Args:
	#   sims: a 2D array of samples (# iter * # chains) 
	if (is.vector(sims)) {
		dim(sims) <- c(length(sims), 1)
	}
	niter <- dim(sims)[1]
	half <- niter / 2
	cbind(sims[1:floor(half), ], sims[ceiling(half + 1):niter, ])
}

is_constant <- function(x, tol = .Machine$double.eps) {
	abs(max(x) - min(x)) < tol
}

monitor <- function(sims, warmup = 0, probs = c(0.05, 0.50, 0.95)) { 
	# print a summary for general simulation results 
	# of 3D array: # iter * # chains * # parameters 
	# Args:
	#   sims: a 3D array described above 
	#   warmup: the number of iterations used for warmup 
	#   probs: probs of summarizing quantiles 
	# Return: 
	#   A summary data.frame of class 'simsummary'
	if (inherits(sims, "stanfit")) {
		chains <- sims@sim$chains
		iter <- sims@sim$iter
		warmup <- sims@sim$warmup
		parnames <- names(sims)
		sims <- as.array(sims)
	} else {
		dim_sims <- dim(sims)
		if (is.null(dim_sims)) {
			dim(sims) <- c(length(sims), 1, 1) 
		} else if (length(dim_sims) == 2) {
			dim(sims) <- c(dim_sims, 1)
		} else if (length(dim_sims) > 3) {
			stop("'sims' has more than 3 dimensions") 
		}
		parnames <- dimnames(sims)[[3]]
		if (is.null(parnames)) {
			parnames <- paste0("V", seq_len(dim(sims)[3]))
		}
		iter <- dim(sims)[1]
		chains <- dim(sims)[2]
		if (warmup > dim(sims)[1]) {
			stop("warmup is larger than the total number of iterations")
		}
		if (warmup >= 1) {
			sims <- sims[-seq_len(warmup), , , drop = FALSE] 
		}
	}
	
	mcse_quan_fun <- function(p, sims) quantile_mcse(sims, p)$mcse
	out <- vector("list", length(parnames))
	out <- setNames(out, parnames)
	for (i in seq_along(out)) {
		sims_i <- sims[, , i]
		valid <- all(is.finite(sims_i))
		quan <- unname(quantile(sims_i, probs = probs))
		mcse_quan <- sapply(probs, mcse_quan_fun, sims_i)
		
		zsims_split <- z_scale(split_chains(sims_i))
		zsplit_rhat <- rhat_rfun(zsims_split)
		bulk_ess <- round(ess_rfun(zsims_split))
		tail_ess <- round(tail_ess(sims_i))
		
		sims_folded <- abs(sims_i - median(sims_i))
		zsims_folded_split <- z_scale(split_chains(sims_folded))
		zfsplit_rhat <- rhat_rfun(zsims_folded_split)
		rhat <- max(zsplit_rhat, zfsplit_rhat)
		
		ess <- ess_rfun(sims_i)
		mean <- mean(sims_i)
		sd <- sd(sims_i)
		mcse_mean <- sd / sqrt(ess)
		# mcse_sd assumes normality of sims and uses Stirling's approximation
		# min of ess for sims and sims^2 is used instead of ess
		ess2 <- ess_rfun(sims_i^2)
		essmin <- min(c(ess,ess2))
		fac_mcse_sd <- sqrt(exp(1) * (1 - 1 / essmin)^(essmin - 1) - 1)
		mcse_sd <- sd * fac_mcse_sd
		
		out[[i]] <- c(
			valid, quan, mean, sd, mcse_quan, mcse_mean,
			mcse_sd, rhat, bulk_ess, tail_ess
		)
	}
	
	out <- as.data.frame(do.call(rbind, out))
	str_quan <- paste0("Q", probs * 100)
	str_mcse_quan <- paste0("MCSE_", str_quan)
	colnames(out) <- c(
		"valid", str_quan, "Mean", "SD", str_mcse_quan, 
		"MCSE_Mean", "MCSE_SD", "Rhat", "Bulk_ESS", "Tail_ESS"
	)
	rownames(out) <- parnames
	
	# replace NAs with appropriate values if draws are valid
	S <- prod(dim(sims)[1:2]) 
	out$Rhat[out$valid & !is.finite(out$Rhat)] <- 1
	out$Bulk_ESS[out$valid & !is.finite(out$Bulk_ESS)] <- S
	out$Tail_ESS[out$valid & !is.finite(out$Tail_ESS)] <- S
	SE_vars <- colnames(out)[grepl("^SE_", colnames(out))]
	for (v in SE_vars) {
		out[[v]][out$valid & !is.finite(out[[v]])] <- 0
	}
	out$valid <- NULL
	
	structure(
		out,
		chains = chains,
		iter = iter,
		warmup = warmup,
		class = c("simsummary", "data.frame") 
	)
} 

monitor_extra <- function(sims, warmup = 0, probs = c(0.05, 0.50, 0.95)) { 
	# print an extended summary for general simulation results 
	# of 3D array: # iter * # chains * # parameters 
	# Args:
	#   sims: a 3D array described above 
	#   warmup: the number of iterations used for warmup 
	#   probs: probs of summarizing quantiles 
	# Return: 
	#   A summary data.frame of class 'simsummary'
	if (inherits(sims, "stanfit")) {
		chains <- sims@sim$chains
		iter <- sims@sim$iter
		warmup <- sims@sim$warmup
		parnames <- names(sims)
		sims <- as.array(sims)
	} else {
		dim_sims <- dim(sims)
		if (is.null(dim_sims)) {
			dim(sims) <- c(length(sims), 1, 1) 
		} else if (length(dim_sims) == 2) {
			dim(sims) <- c(dim_sims, 1)
		} else if (length(dim_sims) > 3) {
			stop("'sims' has more than 3 dimensions") 
		}
		parnames <- dimnames(sims)[[3]]
		if (is.null(parnames)) {
			parnames <- paste0("V", seq_len(dim(sims)[3]))
		}
		iter <- dim(sims)[1]
		chains <- dim(sims)[2]
		if (warmup > dim(sims)[1]) {
			stop("warmup is larger than the total number of iterations")
		}
		if (warmup >= 1) {
			sims <- sims[-seq_len(warmup), , , drop = FALSE] 
		}
	}
	
	mcse_fun <- function(p, sims) quantile_mcse(sims, p)$mcse
	out <- vector("list", length(parnames))
	out <- setNames(out, parnames)
	for (i in seq_along(out)) {
		sims_i <- sims[, , i]
		nsamples <- prod(dim(sims_i))
		mean <- mean(sims_i)
		sd <- sd(sims_i)
		quan <- unname(quantile(sims_i, probs = probs))
		ess <- round(ess_rfun(sims_i))
		ress <- ess / nsamples
		sem <- sd / sqrt(ess)
		rhat <- rhat_rfun(sims_i)
		
		split_ess <- round(ess_rfun(split_chains(sims_i)))
		split_rhat <- rhat_rfun(split_chains(sims_i))
		zess <- round(ess_rfun(z_scale(sims_i)))
		zrhat <- rhat_rfun(z_scale(sims_i))
		
		zsims_split <- z_scale(split_chains(sims_i))
		zsplit_rhat <- rhat_rfun(zsims_split)
		zsplit_ess <- round(ess_rfun(zsims_split))
		zsplit_ress <- zsplit_ess / nsamples
		
		sims_centered <- sims_i - median(sims_i)
		sims_folded <- abs(sims_centered)
		zsims_folded_split <- z_scale(split_chains(sims_folded))
		zfsplit_rhat <- rhat_rfun(zsims_folded_split)
		zfsplit_ess <- round(ess_rfun(zsims_folded_split))
		zfsplit_ress <- zfsplit_ess / nsamples
		
		tail_ess <- round(tail_ess(sims_i))
		tail_ress <- tail_ess / nsamples 
		
		sims_med <- (sims_centered <= 0) * 1
		sims_mad <- ((sims_folded - median(sims_folded)) <= 0) * 1
		medsplit_ess <- round(ess_rfun(z_scale(split_chains(sims_med))))
		medsplit_ress <- medsplit_ess / nsamples 
		madsplit_ess <- round(ess_rfun(z_scale(split_chains(sims_mad))))
		madsplit_ress <- madsplit_ess / nsamples 
		
		out[[i]] <- c(
			mean, sem, sd, quan, ess, ress, split_ess, zess, zsplit_ess, zsplit_ress, 
			rhat, split_rhat, zrhat, zsplit_rhat, zfsplit_rhat, zfsplit_ess, 
			zfsplit_ress, tail_ess, tail_ress, medsplit_ess, medsplit_ress, 
			madsplit_ess, madsplit_ress
		)
	}
	
	out <- as.data.frame(do.call(rbind, out))
	probs_str <- paste0("Q", probs * 100)
	colnames(out) <- c(
		"mean", "se_mean", "sd", probs_str, "seff", "reff", "sseff", "zseff", 
		"zsseff", "zsreff", "Rhat", "sRhat", "zRhat", "zsRhat", "zfsRhat", 
		"zfsseff", "zfsreff", "tailseff", "tailreff", "medsseff", "medsreff", 
		"madsseff", "madsreff"
	)
	rownames(out) <- parnames
	structure(
		out,
		chains = chains,
		iter = iter,
		warmup = warmup,
		extra = TRUE,
		class = c("simsummary", "data.frame") 
	)
}

print.simsummary <- function(x, digits = 3, se = FALSE, ...) {
	atts <- attributes(x)
	px <- x
	if (!se) {
		px <- px[, !grepl("^MCSE_", colnames(px))]
	}
	class(px) <- "data.frame"
	decimal_places <- max(1, digits - 1)
	px$Rhat <- round(px$Rhat, digits = max(2, decimal_places))
	estimates <- setdiff(names(px), c("Rhat", "Bulk_ESS", "Tail_ESS"))
	px[, estimates] <- round(px[, estimates], digits = decimal_places)
	#px[, estimates] <- signif(px[, estimates], digits = digits)
	# add a space between summary and convergence estimates
	names(px)[names(px) %in% "Rhat"] <- " Rhat"
	cat(
		"Inference for the input samples (", atts$chains, 
		" chains: each with iter = ", atts$iter, 
		"; warmup = ", atts$warmup, "):\n\n", sep = ""
	)
	print(px, ...)
	if (!isTRUE(atts$extra)) {
		cat(
			"\nFor each parameter, Bulk_ESS and Tail_ESS are crude measures of \n",
			"effective sample size for bulk and tail quantities respectively (good values is \n",
			"ESS > 400), and Rhat is the potential scale reduction factor on rank normalized\n",
			"split chains (at convergence, Rhat = 1).\n", sep = ""
		)
	}
	invisible(x)
}


plot_ranknorm <- function(theta, n, m = 1, interval = FALSE) {
	df <- data.frame(theta = theta) %>%
		mutate(
			gid = gl(m, n),
			r = r_scale(theta),
			u = u_scale(theta),
			z = z_scale(theta)
		)
	size <- 1.1
	alpha <- if (m > 1) 0.5 else 1
	blue <- color_scheme_get(scheme = 'blue', i = 4)[[1]]
	p2 <- ggplot(df, aes(x = theta, grp = gid)) +
		stat_ecdf(color = blue, size = size, alpha = alpha, pad = FALSE) +
		labs(x = 'theta', y = 'ECDF')
	p3 <- ggplot(df, aes(x = r / (n * m), grp = gid)) +
		stat_ecdf(color = blue, size = size, alpha = alpha, pad = FALSE) +
		labs(x = 'Scaled rank', y = 'ECDF')
	
	if (interval) {
		df <- df %>% mutate(
			psd = sqrt((r + 1) * (n - r + 1) / (n + 2)^2 / (n + 3)),
			pm2sd = (r + 1) / (n + 2) - 1.96 * psd,
			pp2sd = (r + 1) / (n + 2) + 1.96 * psd,
			p975 = qbeta(0.975, r + 1, n - r + 1),
			p025 = qbeta(0.025, r + 1, n - r + 1)
		)
		p2b <- p2 + 
			geom_line(data = df, aes(y = p025), color = blue) +
			geom_line(data = df, aes(y = p975), color = blue) +
			ylab('ECDF + beta 95% interval')
		p3b <- p3 + 
			geom_line(data = df, aes(y = p025), color = blue) +
			geom_line(data = df, aes(y = p975), color = blue) + 
			ylab('ECDF + beta 95% interval')
		p3c <- p3 + 
			geom_line(data = df, aes(y = pm2sd), color = blue) +
			geom_line(data = df, aes(y = pp2sd), color = blue) +
			ylab('ECDF + normal approximated 95% interval')
		out <- gridExtra::grid.arrange(p2b, p3b, p3c, nrow = 1)
	} else {
		p1 <- bayesplot::mcmc_hist(as.data.frame(theta)) +
			xlab('theta')
		p4 <- ggplot(df, aes(x = z, grp = gid)) +
			stat_ecdf(color = blue, size = size, alpha = alpha, pad = FALSE) +
			labs(x = 'z', y = 'ECDF')
		out <- gridExtra::grid.arrange(p1, p2, p3, p4, nrow = 1)
	}
	invisible(out)
}

mcmc_hist_r_scale <- function(x, nbreaks = 50, ...) {
	max <- prod(dim(x)[1:2])
	bayesplot::mcmc_hist(
		r_scale(x), 
		breaks = seq(0, max, by = max / nbreaks) + 0.5,
		...
	) +
		theme(axis.line.y = element_blank())
}

plot_rhat <- function(res) {
	res$par <- rownames(res)
	p1 <- ggplot(res, aes(x = par, y = sRhat)) + 
		geom_point() + 
		ggtitle('Classic split-Rhat') + 
		geom_hline(yintercept = 1.005, linetype = 'dashed') + 
		geom_hline(yintercept = 1) + 
		ylim(c(.99 ,1.26))
	
	p2 <- ggplot(res, aes(x = par, y = zsRhat)) + 
		geom_point() + 
		ggtitle('Rank normalized split-Rhat') + 
		geom_hline(yintercept = 1.005, linetype = 'dashed') + 
		geom_hline(yintercept = 1) + 
		ylim(c(.99,1.26))
	
	p3 <- ggplot(res, aes(x = par, y = zfsRhat)) + 
		geom_point() + 
		ggtitle('Folded rank norm. split-Rhat') + 
		geom_hline(yintercept = 1.005, linetype = 'dashed') + 
		geom_hline(yintercept = 1) + 
		ylim(c(.99,1.26))
	
	gridExtra::grid.arrange(p1, p2, p3, nrow = 1)
}

plot_ess <- function(res) {
	att <- attributes(res)
	max_seff <- with(res,
					 max(c(seff, zsseff, zfsseff, medsseff, madsseff), na.rm = TRUE)
	)
	S <- att$iter * att$chains
	if (!length(S)) S <- max_seff
	ymax <- round(max(S, max_seff * 1.15))
	ylimits <- c(0, ymax)
	res$par <- rownames(res)
	
	p1 <- ggplot(res, aes(x = par, y = seff)) + 
		geom_point() + 
		ggtitle('Classic ESS') + 
		geom_hline(yintercept = 0) +
		geom_hline(yintercept = 400, linetype = 'dashed') + 
		scale_y_continuous(limits = ylimits) 
	
	p2 <- ggplot(res, aes(x = par, y = zsseff)) + 
		geom_point() + 
		ggtitle('Bulk-ESS') + 
		geom_hline(yintercept = 0) + 
		geom_hline(yintercept = 400, linetype = 'dashed') + 
		scale_y_continuous(limits = ylimits) 
	
	p3 <- ggplot(res, aes(x = par, y = tailseff)) + 
		geom_point() + 
		ggtitle('Tail-ESS') + 
		geom_hline(yintercept = 0) + 
		geom_hline(yintercept = 400, linetype = 'dashed') + 
		scale_y_continuous(limits = ylimits) 
	
	p4 <- ggplot(res, aes(x = par, y = medsseff)) + 
		geom_point() + 
		ggtitle('Median-ESS') + 
		geom_hline(yintercept = 0) + 
		geom_hline(yintercept = 400, linetype = 'dashed') + 
		scale_y_continuous(limits = ylimits) 
	
	p5 <- ggplot(res, aes(x = par, y = madsseff)) + 
		geom_point() + 
		ggtitle('MAD-ESS') + 
		geom_hline(yintercept = 0) + 
		geom_hline(yintercept = 400, linetype = 'dashed') + 
		scale_y_continuous(limits = ylimits) 
	
	blank <- grid::grid.rect(gp = grid::gpar(col = "white"), draw = FALSE)
	gridExtra::grid.arrange(p1, p2, p3, blank, p4, p5, nrow = 2)
}

plot_local_ess <- function(fit, par, nalpha = 20, rank = TRUE) {
	if (length(par) != 1L) {
		stop("'par' should be of length 1.")
	}
	if (inherits(fit, "stanfit")) {
		if (!is.character(par)) {
			par <- names(fit)[par]
		}
		sims <- as.array(fit, pars = par)[, , 1]
		params <- as.data.frame(fit, pars = par)
		sampler_params <- get_sampler_params(fit, inc_warmup = FALSE)
		divergent <- do.call(rbind, sampler_params)[, 'divergent__']
		max_depth <- attr(fit@sim$samples[[1]], "args")$control$max_treedepth
		treedepths <- do.call(rbind, sampler_params)[, 'treedepth__']
		params$divergent <- divergent
		params$max_depth <- (treedepths == max_depth) * 1
		params$urank <- u_scale(params[, par])
		params$value <- params[, par]
	} else {
		if (!is.character(par)) {
			par <- dimnames(fit)[[3]][par]
		}
		sims <- fit[, , par]
		params <- data.frame(value = as.vector(sims))
		params$divergent <- 0
		params$max_depth <- 0
		params$urank <- u_scale(params$value)
	}   
	
	# compute local Seff
	delta <- 1 / nalpha
	alphas <- seq(0, 1 - delta, by = delta)
	zsseffs <- rep(NA, length(alphas))
	for (i in seq_along(alphas)) {
		alpha <- alphas[i]
		I <- sims > quantile(sims, alpha) & sims <= quantile(sims, alpha + delta)
		zsseffs[i] <- ess_rfun(z_scale(split_chains(I)))
	}
	S <- prod(dim(I))
	
	# create the plot
	df <- data.frame(
		quantile = seq(0, 1, by = delta),
		value = quantile(params$value, seq(0, 1, by = delta)),
		zsseff = c(zsseffs, zsseffs[nalpha])
	)
	ymax <- max(S, round(max(zsseffs, na.rm = TRUE) * 1.15, 1))
	xname <- if (rank) "quantile" else "value"
	xrug <- if (rank) "urank" else "value"
	out <- ggplot(data = df, aes_string(x = xname, y = "zsseff")) +
		geom_step() + 
		geom_hline(yintercept = c(0, 1)) + 
		geom_hline(yintercept = 400, linetype = 'dashed') + 
		scale_y_continuous(
			breaks = seq(0, ymax, by = round(0.25*S)), 
			limits = c(0, ymax)
		) +
		geom_rug(
			data = params[params$divergent == 1, ], 
			aes_string(x = xrug, y = NULL), sides = "b", color = "red"
		) +
		geom_rug(
			data = params[params$max_depth == 1, ], 
			aes_string(x = xrug, y = NULL), sides = "b", color = "orange"
		) +
		ylab('ESS for small intervals')
	if (rank) {
		out <- out +
			scale_x_continuous(breaks = seq(0, 1, by = 0.1)) + 
			xlab('Quantile')
	} else {
		out <- out + xlab(par)
	}
	out
}

plot_quantile_ess <- function(fit, par, nalpha = 20, rank = TRUE) {
	if (length(par) != 1L) {
		stop("'par' should be of length 1.")
	}
	if (inherits(fit, "stanfit")) {
		if (!is.character(par)) {
			par <- names(fit)[par]
		}
		sims <- as.array(fit, pars = par)[, , 1]
		params <- as.data.frame(fit, pars = par)
		sampler_params <- get_sampler_params(fit, inc_warmup = FALSE)
		divergent <- do.call(rbind, sampler_params)[, 'divergent__']
		max_depth <- attr(fit@sim$samples[[1]], "args")$control$max_treedepth
		treedepths <- do.call(rbind, sampler_params)[, 'treedepth__']
		params$divergent <- divergent
		params$max_depth <- (treedepths == max_depth) * 1
		params$urank <- u_scale(params[, par])
		params$value <- params[, par]
	} else {
		if (!is.character(par)) {
			par <- dimnames(fit)[[3]][par]
		}
		sims <- fit[, , par]
		params <- data.frame(value = as.vector(sims))
		params$divergent <- 0
		params$max_depth <- 0
		params$urank <- u_scale(params$value)
	}
	
	# compute quantile Seff
	delta <- 1 / nalpha
	alphas <- seq(delta, 1 - delta, by = delta)
	zsseffs <- rep(NA, length(alphas))
	for (i in seq_along(alphas)) {
		alpha <- alphas[i]
		I <- sims <= quantile(sims, alpha)
		zsseffs[i] <- ess_rfun(z_scale(split_chains(I)))
	}
	S <- prod(dim(I))
	
	# create the plot
	df <- data.frame(
		quantile = seq(delta, 1 - delta, by = delta), 
		value = quantile(params$value, seq(delta, 1 - delta, by = delta)),
		zsseff = zsseffs
	)
	ymax <- max(S, round(max(zsseffs, na.rm = TRUE) * 1.15, 1))
	xname <- if (rank) "quantile" else "value"
	xrug <- if (rank) "urank" else "value"
	out <- ggplot(data = df, aes_string(x = xname, y = "zsseff")) +
		geom_point() + 
		geom_hline(yintercept = c(0, 1)) + 
		geom_hline(yintercept = 400, linetype = 'dashed') + 
		scale_y_continuous(
			breaks = seq(0, ymax, by = round(0.25*S)), 
			limits = c(0, ymax)
		) +
		geom_rug(
			data = params[params$divergent == 1, ], 
			aes_string(x = xrug, y = NULL), sides = "b", color = "red"
		) +
		geom_rug(
			data = params[params$max_depth == 1, ], 
			aes_string(x = xrug, y = NULL), sides = "b", color = "orange"
		) +
		ylab("ESS for quantiles")
	if (rank) {
		out <- out +
			scale_x_continuous(breaks = seq(0, 1, by = 0.1)) + 
			xlab('Quantile')
	} else {
		out <- out + xlab(par)
	}
	out
}

plot_change_ess <- function(fit, par, breaks = seq(0.1, 1, 0.05), 
							yaxis = c("absolute", "relative")) {
	if (length(par) != 1L) {
		stop("'par' should be of length 1.")
	}
	yaxis <- match.arg(yaxis)
	if (inherits(fit, "stanfit")) {
		if (!is.character(par)) {
			par <- names(fit)[par]
		}
		sims <- as.array(fit, pars = par)[, , 1]
	} else {
		if (!is.character(par)) {
			par <- dimnames(fit)[[3]][par]
		}
		sims <- fit[, , par]
	}
	
	iter_breaks <- round(breaks * NROW(sims))
	nbreaks <- length(iter_breaks)
	bulk_seff <- tail_seff <- bulk_reff <- 
		tail_reff <- rep(NA, length(nbreaks))
	for (i in seq_along(iter_breaks)) {
		sims_i <- sims[seq_len(iter_breaks[i]), ]
		nsamples <- prod(dim(sims_i))
		bulk_seff[i] <- ess_rfun(z_scale(split_chains(sims_i)))
		tail_seff[i] <- tail_ess(sims_i)
		bulk_reff[i] <- bulk_seff[i] / nsamples
		tail_reff[i] <- tail_seff[i] / nsamples
	}
	df <- data.frame(
		breaks = breaks,
		ndraws = iter_breaks * NCOL(sims),
		seff = c(bulk_seff, tail_seff),
		reff = c(bulk_reff, tail_reff), 
		type = rep(c("bulk", "tail"), each = nbreaks)
	)
	blues <- bayesplot::color_scheme_get(scheme = "blue", i = c(4, 2))
	blues <- unname(unlist(blues))
	if (yaxis == "absolute") {
		out <- ggplot(df, aes(ndraws, seff, color = type)) +
			ylab("ESS") +
			geom_hline(yintercept = 0, linetype = 1) +
			geom_hline(yintercept = 400, linetype = 2)
	} else if (yaxis == "relative") {
		out <- ggplot(df, aes(ndraws, reff, color = type)) +
			ylab("Relative efficiency") +
			geom_hline(yintercept = 0, linetype = 2)
	}
	out +  
		geom_line() +
		geom_point() +
		xlab("Total number of draws") +
		scale_colour_manual(values = blues)
}