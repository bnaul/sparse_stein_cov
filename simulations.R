# NOTE use 1/rho instead of rho (original parametrization)
# NOTE rho==Inf --> soft thresholding? is this used in simulations?

#install.packages(c('mvtnorm', 'ggplot2', 'reshape', 'svd', 'spcov'))
#source("http://bioconductor.org/biocLite.R"); biocLite(); biocLite("genefilter")
library(mvtnorm)
library(ggplot2)
library(reshape)
library(svd)
library(spcov)
library(Biobase)
library(genefilter)
source('GGDescent_patch.R')
source('cglasso.R')
source('sparse_stein.R')
source('sparse_mle.R')

theme_set(theme_bw(base_size=16))

# Helper functions
stein_loss <- function(Sigma_hat, Sigma_inv) {
  if (any(is.na(Sigma_hat))) return(NA)
  return(sum(diag(Sigma_hat%*%Sigma_inv)) - log(det(Sigma_hat%*%Sigma_inv)) - nrow(Sigma_inv))
}
nz_vals <- function(A,tol=1e-8) {A[upper.tri(A) & abs(A) > tol]}
nz_count <- function(A,tol=1e-8) {length(nz_vals(A,tol))}
lseq <- function(from, to, length.out) exp(seq(log(from), log(to), length.out=length.out))
# Cross-validation functions; adapted from code provided by Sang Oh, UC Berkeley
getbeta <- function(omega) {
  omegaD <- diag(diag(omega))
  return(-omegaD %*% (omega - omegaD))
}
getresid <- function(X, omega){
  if (any(is.na(omega))) return(NA)
  norm(X %*% (diag(ncol(omega)) - getbeta(omega)), 'F')
}
# Func takes inputs X, lambda and returns list(Omega=Sigma^-1 est.)
# Returns matrix of errors; colMeans(cv(...)) gives avg. error for each lambda
cv <- function(lvec, X, func, fold=min(nrow(X),10), type='regression') {
  n <- nrow(X)
  p <- ncol(X)
  g <- length(lvec)
  S <- cov(X)

  sections <- cut(1:n, fold, labels=1:fold)

  negll <- matrix(0, fold, g)
  for (i in c(1:fold)) {
    tsindx <- which(sections==i)
    trindx <- which(sections!=i)

    Xtr <- X[trindx,,drop=FALSE]; ntr <- nrow(Xtr)
    Xts <- X[tsindx,,drop=FALSE]; nts <- nrow(Xts)

    meanresid <- lapply(lvec,
      function(lambda, Xtr, Xts){
        sol <- func(Xtr,lambda)$Omega
        if (type=='regression') {
          return(getresid(Xts,sol)/nrow(Xts))
        } else if (type=='likelihood') {
          return(-log(det(sol)) + sum(sol*S))
        } else stop('Unknown cv type')
      }, Xtr=Xtr, Xts=Xts)
      negll[i,] <- unlist(meanresid)
    }

  return(negll)
}

# Find penalty parameter that yields a specific number of edges (for plotting)
# Inefficient, but the resulting plots are more uniform than if we just iterate over a
# grid of penalty parameters
lambda_search <- function(func, S, n, num_edges, lambda_min=0, lambda_max=5, N_lambda=6) {
  while (TRUE) {
    lambda_list <- seq(lambda_min, lambda_max, len=N_lambda)
    nnz_list <- numeric(N_lambda)
    for (i in 1:N_lambda) {
      out <- func(S, n, lambda_list[i])
      nnz_list[i] <- length(nz_vals(out))
    }
    if (is.unsorted(rev(nnz_list))) {warning('Non-monotonic sparsity')}
    if (any(nnz_list == num_edges)) {break}
    else if (num_edges < min(nnz_list) || num_edges > max(nnz_list)) {
      warning('lambda search failed')
      return(1e3)
    }
    else {
      lambda_min <- max(lambda_list[nnz_list>num_edges])
      lambda_max <- min(lambda_list[nnz_list<num_edges])
      #			cat(sprintf('New lambda range: %f (%i)-%f (%i), delta=%e\n', lambda_min, nnz_list[1], lambda_max, nnz_list[length(nnz_list)], lambda_max-lambda_min))
    }
  }
  return(min(lambda_list[whichClosest(nnz_list, num_edges)]))
}

# Iterate over grid of penalty parameters and return cv scores / sparsity levels
cv_paths <- function(estimator, S, X, n, lambda_range, rho_range) {
  cv_table <- t(sapply(lambda_range, function(lambda) {
    colMeans(cv(rho_range, X, function(X_i,rho) {
      list(Omega=solve(estimator(cov(X_i), nrow(X_i), lambda, rho)))}
    ))}
  ))
  nnz_table <- outer(lambda_range, rho_range,
                     Vectorize(function(lambda, rho) {
                       nz_count(estimator(S,n,lambda,rho))
                     }))
  return(list(cv_table=cv_table, nnz_table=nnz_table))
}

# Combine multiple simulation outputs into single data frame for plotting
aggregate_data <- function(full_data, x_var, y_var, span, max_points, separate_rho=FALSE) {
  if (length(unique(full_data$k)) != 1) warning('Multiple test cases Sigma_k combined')
  plot_data <- full_data[NULL,c('k','algorithm','rho',x_var,y_var)]
  if (!separate_rho) full_data$rho <- Inf
  for (algorithm in unique(full_data$algorithm)) {
    for (rho in unique(full_data$rho)) { # for spcov, rho=Inf, so this still applies
      to_smooth <- full_data[full_data$algorithm==algorithm & full_data$rho==rho,]
      if (nrow(to_smooth) > 0) {
        if (all(is.infinite(to_smooth[,y_var]))) {
          smoothed <- list(x=to_smooth[x_var][,1], y=to_smooth[y_var][,1])
          inds <- rep(TRUE,length(smoothed$x))
        }
        else {
          if (x_var == 'sensitivity' && any(to_smooth[,x_var]==1)) {
            to_smooth <- to_smooth[!(to_smooth[,x_var]==1 & to_smooth[,y_var]!=max(to_smooth[to_smooth[,x_var]==1,y_var])),]
          }
          if (algorithm %in% c('Xue','Liu')) { # smooth Xue, Liu bc of extra noise
            bass <- 1
          } else {
            bass <- 0
          }
          smoothed <- supsmu(to_smooth[,x_var], to_smooth[,y_var], span='cv', bass=bass)
          width <- diff(range(smoothed$x))/max_points# maximum ~150 points
          inds <- c(TRUE,diff(round(smoothed$x/width)*width)!=0)
          if (x_var == 'sensitivity' && any(to_smooth[,x_var]==1)) {
            smoothed$x <- c(smoothed$x, 1)
            smoothed$y <- c(smoothed$y, 0)
            inds <- c(inds, TRUE)
          }
        }
        plot_data <- rbind(plot_data, data.frame(full_data$k[1], algorithm, rho, x=smoothed$x[inds], y=smoothed$y[inds]))
      }
    }
  }
  colnames(plot_data)[(ncol(plot_data)-1):(ncol(plot_data))] <- c(x_var, y_var)
  plot_data$Type <- paste(plot_data$algorithm, ifelse(plot_data$rho!=Inf & length(unique(plot_data$rho))>2 & separate_rho, paste(", rho=",plot_data$rho,sep=""), ''), sep='')
  # have to manually set levels so that they're ordered properly
  plot_data$Type <- factor(plot_data$Type, levels=unique(plot_data$Type))
  if (x_var %in% c('sensitivity','specificity')) {
    plot_data$sensitivity <- pmax(pmin(plot_data$sensitivity, 1), 0)
    plot_data$specificity <- pmax(pmin(plot_data$specificity, 1), 0)
  }
  return(plot_data)
}

# Compare estimated matrix to ground truth and return summary statistics
check_Sigma <- function(Sigma_hat, Sigma, tol=1e-4) {
    if (!all(dim(Sigma_hat) == dim(Sigma))) stop('Incompatible dimensions')
    # Sigma=NA if estimation failed, just return NAs
    if (all(is.na(Sigma))) {result[1:length(result)] <- NA; return(result)}
    result <- data.frame(tp=integer(1), fp=integer(1), tn=integer(1), fn=integer(1), rmse=numeric(1), entropy=numeric(1), min_eig=numeric(1))

    Sigma_hat_offdiag <- Sigma_hat[upper.tri(Sigma_hat, diag=FALSE)]
    Sigma_offdiag <- Sigma[upper.tri(Sigma, diag=FALSE)]
    result$tp <- as.numeric(sum(abs(Sigma_hat_offdiag) >= tol & abs(Sigma_offdiag) >= tol))
    result$fp <- as.numeric(sum(abs(Sigma_hat_offdiag) >= tol & abs(Sigma_offdiag) <  tol))
    result$tn <- as.numeric(sum(abs(Sigma_hat_offdiag) <  tol & abs(Sigma_offdiag) <  tol))
    result$fn <- as.numeric(sum(abs(Sigma_hat_offdiag) <  tol & abs(Sigma_offdiag) >= tol))

    result$rmse <- norm(Sigma_hat-Sigma, 'F')/nrow(Sigma)
    result$entropy <- stein_loss(Sigma_hat, solve(Sigma))
    result$min_eig <- min(eigen(Sigma_hat, symmetric=TRUE)$val)
    return(result)
}

# Produce plots like those from Bien & Tibshirani 2011, Section 5.1 for a given
# set of simulation outputs
bien_tibs_plots <- function(plot_data, span=0.02, max_points=75, separate_rho=FALSE) {
  plots <- list()
  label_size <- 28
  plots[[1]] <- ggplot(aggregate_data(plot_data, 'sensitivity', 'specificity', span=span, max_points=max_points, separate_rho=separate_rho), aes(x=sensitivity, y=specificity, color=Type, shape=Type)) + geom_line() + geom_point() + expand_limits(y=c(0,1)) + theme(legend.position="right") + theme(text=element_text(size=label_size))# + theme(axis.title.x = element_text(size=label_size), axis.title.y = element_text(size=label_size))
  plots[[2]] <- ggplot(aggregate_data(plot_data, 'nnz', 'rmse', span=span, max_points=max_points, separate_rho=separate_rho), aes(x=nnz, y=rmse, color=Type, shape=Type)) + geom_line() + geom_point() + expand_limits(y=0) + theme(legend.position="right") + xlab('number of non-zeros') + ylab('rmse') + theme(text=element_text(size=label_size))
  plots[[3]] <- ggplot(aggregate_data(plot_data, 'nnz', 'entropy', span=span, max_points=max_points, separate_rho=separate_rho), aes(x=nnz, y=entropy, color=Type, shape=Type)) + geom_line() + geom_point() + expand_limits(y=0) + theme(legend.position="right") + xlab('number of non-zeros') + ylab('entropy loss') + theme(text=element_text(size=label_size))
  plots[[4]] <- ggplot(aggregate_data(plot_data, 'nnz', 'min_eig', span=span, max_points=max_points, separate_rho=separate_rho), aes(x=nnz, y=min_eig, color=Type, shape=Type)) + geom_line() + geom_point() + expand_limits(y=0 + theme(legend.position="right")) + xlab('number of non-zeros') + ylab('minimum eigenvalue') + theme(text=element_text(size=label_size))
  names(plots) <- paste(c('roc','rmse','entropy','min_eig'), '_Sigma', plot_data$k[1], sep='')
  return(plots)
}

rand_sign_symm <- function(p, neg_prob=0.5) {
  A <- matrix(sign(runif(p^2)-neg_prob),p,p)
  A[lower.tri(A)] <- 0
  A <- A + t(A)
  diag(A) <- 1
  return(A)
}
scale_cond <- function(A, k_star) {
  lambda <- eigen(A, symmetric=TRUE)$val
  d_star <- (max(lambda) - k_star*min(lambda))/(k_star - 1)
  return(A + d_star*diag(nrow(A)))
}
# Test covariance matrices from Bien & Tibshirani 2011
bien_tibs_sigmas <- function(p, num_blocks, nz_prob, MA_const, kappa_star) {
  block_size <- p / num_blocks
  i <- rep(1:p, p)                        # row indices for all elements of Sigma
  j <- rep(1:p, each=p)                   # col indices for all elements of Sigma
  block_inds_i <- ceiling(i / block_size) # block containing row i
  block_inds_j <- ceiling(j / block_size) # block containing col j
  Sigma_list <- list()

  Sigma_list[[1]] <- rand_sign_symm(p)
  Sigma_list[[1]][block_inds_i != block_inds_j] <- 0
  Sigma_list[[1]] <- scale_cond(Sigma_list[[1]], kappa_star)

  Sigma_list[[2]] <- rand_sign_symm(p)
  Sigma_list[[2]][block_inds_i != block_inds_j] <- 0
  Sigma_list[[2]][((i %% block_size) != 0) & ((j %% block_size) != 0) & (i != j)] <- 0
  Sigma_list[[2]] <- scale_cond(Sigma_list[[2]], kappa_star)

  Sigma_list[[3]] <- rand_sign_symm(p)
  zero_inds <- matrix(runif(p^2) > nz_prob, p, p)
  zero_inds[lower.tri(zero_inds, diag=TRUE)] <- 0
  Sigma_list[[3]][zero_inds | t(zero_inds)] <- 0
  Sigma_list[[3]] <- scale_cond(Sigma_list[[3]], kappa_star)

  Sigma_list[[4]] <- diag(p)
  Sigma_list[[4]][cbind(1:(p-1),2:p)] <- MA_const
  Sigma_list[[4]][cbind(2:p,1:(p-1))] <- MA_const

  return(Sigma_list)
}

cov_glasso_sim <- function(estimators, p, n, N, num_sparsity_levels=401, save_data=TRUE) {
  set.seed(3)
  num_blocks <- 5
  nz_prob <- 0.02
  MA_const <- 0.4
  Sigma_list <- bien_tibs_sigmas(p, num_blocks, nz_prob, MA_const, p)

  lambda_range <- rev(lseq(1e-5, 30, num_sparsity_levels)[-1])
  rho_range <- 1/c(lseq(2, 1000, len=20),Inf)
  tau_range <- seq(0, 2, len=21)
  num_edges_range <- round(seq(1, choose(p,2)*0.9, len=num_sparsity_levels))[-1]

  full_data <- data.frame(k=integer(), i=integer(),
                          algorithm=factor(levels=c('Cov Glasso','Soft','Stein Soft',
                                                    'Sparse Stein','Stein Hard','ECM',
                                                    'Adaptive Sparse Stein','Xue','Liu')),
                          lambda=numeric(), rho=numeric(), tp=integer(),
                          fp=integer(), tn=integer(), fn=integer(),
                          rmse=numeric(), entropy=numeric(), min_eig=numeric(),
                          cv_score=numeric())
  all_plots <- list()
  for (k in seq_along(Sigma_list)) {
    X_list <- replicate(N, rmvnorm(n, rep(0,p), Sigma_list[[k]]), simplify=FALSE)
    S_list <- lapply(X_list, function(X_i) t(X_i) %*% X_i)
    for (i in 1:N) {
      cat(sprintf('k=%i, i=%i (%s)\n', k, i, Sys.time()))
      set.seed(i)
      if ('Cov Glasso' %in% estimators) {
        for (lambda in lambda_range) {
#          cat(sprintf('k=%i, i=%i, lambda=%f (%s)\n', k, i, lambda, Sys.time()))
          spcov_out <- tryCatch(spcov(S_list[[i]]/n+1e-3*diag(p), S_list[[i]]/n+1e-3*diag(p), lambda, 100), error=function(e) list(Sigma=NA))
          spcov_stats <- check_Sigma(spcov_out$Sig, Sigma_list[[k]])
          full_data[nrow(full_data)+1,] <- c(k, i, 'Cov Glasso', lambda, Inf, spcov_stats)
        }
      }

      # Used for test case where a single estimator is selected via cross-validation
      if ('Cov Glasso CV' %in% estimators) {
        cv_out <- cv_paths(function(S,n,lambda,rho) tryCatch(spcov(S_list[[i]]/n+1e-3*diag(p), S_list[[i]]/n+1e-3*diag(p), lambda, 100)$Sig, error=function(e) NA), S_list[[i]]/n, X_list[[i]], n, lambda_range, Inf)
        for (j in seq_along(lambda_range)) {
          lambda <- lambda_range[j]
          cv_score <- cv_out$cv_table[j]
          spcov_out <- tryCatch(spcov(S_list[[i]]/n+1e-3*diag(p), S_list[[i]]/n+1e-3*diag(p), lambda, 100), error=function(e) list(Sigma=NA))
          spcov_stats <- check_Sigma(spcov_out$Sig, Sigma_list[[k]])
          full_data[nrow(full_data)+1,] <- c(k, i, 'Cov Glasso CV', lambda, Inf, spcov_stats, cv_score)
        }
      }

      if ('Sparse Stein' %in% estimators) {
        cv_out <- cv_paths(function(S,n,lambda,rho) sparse_stein_cov(n*S,n,lambda,rho)$Theta, S_list[[i]]/n, X_list[[i]], n, lambda_range, rho_range)
        for (j in seq_along(num_edges_range)[-1]) {
          matching_inds <- which(cv_out$nnz_table >= num_edges_range[j-1] & cv_out$nnz_table <= num_edges_range[j], arr.ind=TRUE)
          if (nrow(matching_inds) > 0) {
            inds_star <- matching_inds[which.min(cv_out$cv_table[matching_inds]),]
            lambda_star <- lambda_range[inds_star[1]]
            rho_star <- rho_range[inds_star[2]]
            cv_score <- cv_out$cv_table[inds_star[1],inds_star[2]]
            if (nrow(full_data)==0 || !(lambda_star==full_data$lambda[nrow(full_data)] && rho_star==full_data$rho[nrow(full_data)])) {
              stein_out <- sparse_stein_cov(S_list[[i]], n, lambda_star, rho_star, delta=0)
              stein_stats <- check_Sigma(stein_out$Theta, Sigma_list[[k]])
              full_data[nrow(full_data)+1,] <- c(k, i, 'Sparse Stein', lambda_star, rho_star, stein_stats, cv_score)
            }
          }
        }
      }

      if ('Adaptive Sparse Stein' %in% estimators) {
        cv_out <- cv_paths(function(S,n,lambda,rho) sparse_stein_cov(n*S,n,lambda,rho,delta=2)$Theta, S_list[[i]]/n, X_list[[i]], n, lambda_range, rho_range)
        for (j in seq_along(num_edges_range)[-1]) {
          matching_inds <- which(cv_out$nnz_table >= num_edges_range[j-1] & cv_out$nnz_table <= num_edges_range[j], arr.ind=TRUE)
          if (nrow(matching_inds) > 0) {
            inds_star <- matching_inds[which.min(cv_out$cv_table[matching_inds]),]
            lambda_star <- lambda_range[inds_star[1]]
            rho_star <- rho_range[inds_star[2]]
            cv_score <- cv_out$cv_table[inds_star[1],inds_star[2]]
            if (nrow(full_data)==0 || !(lambda_star==full_data$lambda[nrow(full_data)] && rho_star==full_data$rho[nrow(full_data)])) {
              stein_out <- sparse_stein_cov(S_list[[i]], n, lambda_star, rho_star, delta=2)
              stein_stats <- check_Sigma(stein_out$Theta, Sigma_list[[k]])
              full_data[nrow(full_data)+1,] <- c(k, i, 'Adaptive Sparse Stein', lambda_star, rho_star, stein_stats, cv_score)
            }
          }
        }
      }

      if ('Xue' %in% estimators) {
        cv_out <- cv_paths(function(S,n,lambda,tau) xue_cov(S,lambda,tau,tol=1e-3,max_it=100), S_list[[i]]/n, X_list[[i]], n, lambda_range, tau_range)
        for (j in seq_along(num_edges_range)[-1]) {
          matching_inds <- which(cv_out$nnz_table >= num_edges_range[j-1] & cv_out$nnz_table <= num_edges_range[j], arr.ind=TRUE)
          if (nrow(matching_inds) > 0) {
            inds_star <- matching_inds[which.min(cv_out$cv_table[matching_inds]),]
            lambda_star <- lambda_range[inds_star[1]]
            tau_star <- tau_range[inds_star[2]]
            cv_score <- cv_out$cv_table[inds_star[1],inds_star[2]]
            if (nrow(full_data)==0 || !(lambda_star==full_data$lambda[nrow(full_data)] && tau_star==full_data$rho[nrow(full_data)])) {
              xue_out <- xue_cov(S_list[[i]]/n, lambda_star, tau_star)
              xue_stats <- check_Sigma(xue_out, Sigma_list[[k]])
              full_data[nrow(full_data)+1,] <- c(k, i, 'Xue', lambda_star, tau_star, xue_stats, cv_score)
            }
          }
        }
      }

      if ('Liu' %in% estimators) {
        cv_out <- cv_paths(function(S,n,lambda,tau) liu_cov(S,lambda,tau,tol=1e-3,max_it=100), S_list[[i]]/n, X_list[[i]], n, lambda_range, tau_range)
        for (j in seq_along(num_edges_range)[-1]) {
          matching_inds <- which(cv_out$nnz_table >= num_edges_range[j-1] & cv_out$nnz_table <= num_edges_range[j], arr.ind=TRUE)
          if (nrow(matching_inds) > 0) {
            inds_star <- matching_inds[which.min(cv_out$cv_table[matching_inds]),]
            lambda_star <- lambda_range[inds_star[1]]
            tau_star <- tau_range[inds_star[2]]
            cv_score <- cv_out$cv_table[inds_star[1],inds_star[2]]
            if (nrow(full_data)==0 || !(lambda_star==full_data$lambda[nrow(full_data)] && tau_star==full_data$rho[nrow(full_data)])) {
              liu_out <- liu_cov(S_list[[i]]/n, lambda_star, tau_star)
              liu_stats <- check_Sigma(liu_out, Sigma_list[[k]])
              full_data[nrow(full_data)+1,] <- c(k, i, 'Liu', lambda_star, tau_star, liu_stats, cv_score)
            }
          }
        }
      }

      full_data$sensitivity <- full_data$tp/(full_data$tp + full_data$fn)
      full_data$specificity <- full_data$tn/(full_data$fp + full_data$tn)
      full_data$nnz <- full_data$tp + full_data$fp

      if (!exists('filename')) {
        if (length(unique(full_data$i)) > 1) {
          filename <- sprintf('%s_p%i_k%i.RData', paste(gsub(' ', '', tolower(unique(full_data$algorithm))), collapse='_'), p, k)
        } else {
          filename <- sprintf('%s_p%i_k%i_i%i.RData', paste(gsub(' ', '', tolower(unique(full_data$algorithm))), collapse='_'), p, k, i)
        }
      }
      if (save_data) save.image(filename)  # save before plotting since ggplot can crash
      all_plots <- c(all_plots, unlist(lapply(split(full_data, full_data$k), bien_tibs_plots), recursive=FALSE))
      if (save_data) save.image(filename)  # re-save with plots
    }
  }
  return(list(data=full_data, plots=all_plots))
}

# Should run in < 1hr; can reduce # of penalty parameters tested for sparser/faster plot
times_small <- function() {
  set.seed(1)
  load('nkitiny.RData')
  n <- ncol(exprs(nkitiny))
  p <- nrow(exprs(nkitiny))
  S <- cov(t(exprs(nkitiny)))

  lambda_range <- lseq(0.2, 15, len=101)
  spcov_times <- numeric(length(lambda_range))
  spcov_nnz_list <- integer(length(lambda_range))
  cd_times <- numeric(length(lambda_range))
  cd_nnz_list <- integer(length(lambda_range))
  for (i in seq_along(lambda_range)) {
    spcov_times[i] <- as.numeric(system.time({
      spcov_out <- spcov(S, S, lambda_range[i], step.size=100)
      spcov_nnz_list[i] <- nz_count(spcov_out$Sigma)
    }))[1]
    cat(sprintf('cgl: %i, %fs\n', spcov_nnz_list[i], spcov_times[i]))
    cd_times[i] <- as.numeric(system.time({
      cd_out <- CovGlassoCD(S, matrix(lambda_range[i],p,p), S+1e-3*diag(p), 1e-4, 1e-2, 1e4, 1e4)
      cd_nnz_list[i] <- nz_count(cd_out$Sig)
    }))[1]
    cat(sprintf('cd: %i, %fs\n', cd_nnz_list[i], cd_times[i]))
  }

  rho_range <- 1/c(lseq(2,1000,len=20),Inf)
  lambda_range <- lseq(0.01, 15.0, len=101)
  stein_times <- numeric(length(lambda_range))
  stein_nnz_list <- integer(length(lambda_range))
  for (i in seq_along(lambda_range)) {
    stein_times[i] <- as.numeric(system.time({
      for (rho in rho_range) {
        stein_out <- sparse_stein_cov(n*S, n, lambda_range[i], rho)
      }
      stein_nnz_list[i] <- nz_count(stein_out$Theta)
    }))[1]
    cat(sprintf('stein: %i, %fs\n', cd_nnz_list[i], cd_times[i]))
  }

  tau_range <- rev(seq(0,2,len=21))
  lambda_range <- lseq(0.01, 15.0, len=101)
  xue_times <- numeric(length(lambda_range))
  xue_nnz_list <- integer(length(lambda_range))
  for (i in seq_along(lambda_range)) {
    xue_times[i] <- as.numeric(system.time({
      for (tau in tau_range) {
        xue_out <- xue_cov(S, lambda_range[i], tau)
      }
      xue_nnz_list[i] <- nz_count(xue_out)
    }))[1]
    cat(sprintf('xue: %i, %fs\n', xue_nnz_list[i], xue_times[i]))
  }

  size <- 6
  time_data <- data.frame(algorithm=factor(levels=c('Cov Glasso','Coordinate Descent','Sparse Stein','Xue et al.')), nnz=integer(), time=numeric())
  time_data[1:length(spcov_times),] <- data.frame(algorithm='Cov Glasso', nnz=spcov_nnz_list, time=spcov_times)
  time_data <- rbind(time_data, data.frame(algorithm='Sparse Stein', nnz=stein_nnz_list, time=stein_times))
  time_data <- rbind(time_data, data.frame(algorithm='Coordinate Descent', nnz=cd_nnz_list, time=cd_times))
  time_data <- rbind(time_data, data.frame(algorithm='Xue et al.', nnz=xue_nnz_list, time=xue_times))
  gp <- ggplot(time_data, aes(x=nnz, y=time, color=algorithm, shape=algorithm)) + geom_point() + ylab('time (sec)') + xlab('number of edges') + theme(legend.title=element_blank()) + theme(text=element_text(size=16)) + scale_y_log10(breaks=10^((-2):2))
  return(list(times=time_data, plot=gp))
}

# Values are hard-coded here since they were computed on many servers in parallel
# Running everything below on a single machine will take days/weeks
times_reduced <- function() {
  load('nkireduced.RData')
  n <- ncol(exprs(nkireduced))
  p <- nrow(exprs(nkireduced))
  S <- cov(t(exprs(nkireduced)))

# Uncomment below to re-run:
#  lambda_range <- lseq(0.01, 2.0, len=21)
#  rho_range <- 1/c(lseq(50,1000,len=20),Inf)
#  stein_times <- numeric(length(lambda_range))
#  stein_nnz_list <- integer(length(lambda_range))
#  for (i in seq_along(lambda_range)) {
#    stein_times[i] <- as.numeric(system.time({
#      for (rho in rho_range) {
#        stein_out <- sparse_stein_cov(n*S, n, lambda_range[i], rho)
#      }
#      stein_nnz_list[i] <- nz_count(stein_out$Theta)
#    }))[1]
#    print(c(stein_nnz_list[i], stein_times[i]))
#  }
#
#  lambda_range <- rev(lseq(1e-3, 15, len=101))
#  cd_times <- numeric(length(lambda_range))
#  cd_nnz_list <- integer(length(lambda_range))
#  for (i in seq_along(lambda_range)) {
#    cd_times[i] <- as.numeric(system.time({
#      cd_out <- CovGlassoCD(S, matrix(lambda_range[i],p,p), diag(diag(S)), 1e-4, 1e-2, 1e4, 1e4)
#      cd_nnz_list[i] <- nz_count(cd_out$Sig)
#    }))[1]
#    print(c(cd_nnz_list[i], cd_times[i]))
#  }
#
#  tau_range <- rev(seq(0, 2,len=21))
#  lambda_range <- lseq(0.01, 2.0, len=101)
#  xue_times <- numeric(length(lambda_range))
#  xue_nnz_list <- integer(length(lambda_range))
#  for (i in seq_along(lambda_range)) {
#    xue_times[i] <- as.numeric(system.time({
#      for (tau in tau_range) {
#        xue_out <- xue_cov(S, lambda_range[i], tau)
#      }
#      xue_nnz_list[i] <- nz_count(xue_out)
#    }))[1]
#    print(c(xue_nnz_list[i], xue_times[i]))
#  }

  size <- 6
  time_data <- data.frame(algorithm=factor(levels=c('Cov Glasso','Coordinate Descent','Xue', 'Sparse Stein')), nnz=integer(), time=numeric())

  # Comment out if re-running:
  spcov_nnz_list <- c(2488, 3968, 13586, 21631, 38634, 66546, 78329, 88153, 95818, 264294)
  spcov_times <- c(84600, 324972.684258576, 165499.694826901, 133689.043754563, 
    266017.106939536, 435230.816833321, 311793.861372694, 454998.208480887, 
    362271.350709312, 464040)
  stein_times <- c(99.144, 92.28, 87.744, 82.584, 82.776, 78.504, 75.336, 71.04, 70.392, 67.488, 63.504, 62.832, 58.296, 55.344, 55.344, 51.072, 50.52, 46.848, 46.536, 42.528, 42.84, 42.552, 39.192, 38.088, 36.024, 34.656, 34.2, 34.992, 31.344, 31.056, 30.768, 28.32, 27.984, 27.96, 27.816, 27.864, 25.728, 25.032, 24.576, 25.176, 24.384, 25.056, 24.24, 24.216, 24.864, 21.84, 21.864, 22.584, 21.816, 21.72, 22.464, 21.696, 21.576, 22.608, 21.672, 21.816, 22.464, 21.72, 21.768, 22.416, 21.792, 19.536, 19.416, 20.304, 19.416, 19.344, 19.44, 19.56, 20.232, 19.392, 19.488, 19.56, 19.392, 19.464, 20.256, 19.512)
  stein_nnz_list <- c(546597, 536860, 526407, 515002, 502970, 490246, 476798, 462160, 446735, 430466, 413255, 395538, 376790, 357152, 336975, 316246, 294914, 273383, 251707, 230087, 208721, 187859, 167826, 148578, 130459, 113137, 97540, 83185, 70402, 58759, 48757, 39922, 32492, 26015, 20748, 16173, 12681, 9850, 7666, 5963, 4704, 3749, 3015, 2420, 1939, 1601, 1317, 1103, 913, 754, 623, 530, 423, 358, 296, 235, 200, 172, 143, 128, 98, 79, 61, 52, 45, 38, 28, 23, 18, 12, 8, 6, 4, 3, 2, 1)
	xue_times <- c(198.102 ,198.494 ,197.898 ,199.965 ,196.549 ,196.546 ,211.799 ,206.183 ,199.562 ,192.621 ,187.215 ,182.442 ,183.776 ,179.591 ,177.18 ,177.897 ,169.029 ,165.345 ,159.476 ,167.41 ,165.941 ,157.104 ,151.232 ,150.227 ,183.844 ,175.98 ,166.585 ,171.981 ,162.249 ,157.505 ,153.898 ,149.881 ,144.942 ,147.33 ,147.932 ,147.326 ,144.437 ,142.064 ,140.716 ,139.372 ,140.179 ,145.632 ,148.077 ,150.015 ,150.848 ,153.382 ,158.345 ,160.928 ,169.938 ,171.832 ,178.726 ,187.341 ,200.663 ,208.22 ,224.133 ,250.315 ,266.239 ,288.847 ,305.641 ,342.168 ,371.504 ,447.324 ,550.622 ,583.397 ,614.822 ,734.59 ,831.176 ,910.093 ,916.704 ,939.14 ,924.657 ,946.192 ,942.634 ,935.992 ,936.296 ,928.33 ,934.523 ,936.565 ,939.8 ,927.803 ,930.229 ,936.881 ,928.513 ,922.2 ,918.641 ,919.09 ,922.497 ,924.554 ,918.86 ,922.185 ,921.899 ,840.706 ,535.088 ,328.967 ,241.298 ,209.082 ,159.44 ,134.194 ,117.052 ,104.412 ,1.306)*4
	xue_nnz_list <- c(664789, 664627, 664442, 664230, 664020, 663773, 663476, 663158, 662821, 662477, 662128, 661766, 661313, 660799, 660323.00, 659771, 659198, 658558, 657909, 657195.00, 656410, 655620, 654785, 653901, 652847, 651758.00, 650621, 649445, 648166, 646767, 645282, 643743, 641943, 640120.00, 638155, 635898, 633560, 631135, 628524, 625635, 622548, 619164, 615525, 611566, 607226, 602535, 597336, 591673, 585413, 578686, 571605, 563490, 554866, 545451.00, 535105, 523752, 511298, 497733, 482987, 466799, 449231, 430281, 409635, 387626, 364122, 339311.00, 313308, 286557, 258924, 231492.00, 204272, 177580, 152407, 128701, 107107, 87963.00, 71400, 57143, 45169.0, 35212, 26960, 20353, 15241, 11392.0, 8350, 6153.00, 4578, 3446, 2627.00, 2034, 1583, 1223, 963, 748, 587, 448, 347.00, 281, 208, 163, 132)
	cd_nnz_list <- c(NA)
	cd_times <- c(NA)

  time_data[1:length(spcov_times),] <- data.frame(algorithm='Cov Glasso', nnz=spcov_nnz_list, time=spcov_times)
  time_data <- rbind(time_data, data.frame(algorithm='Sparse Stein', nnz=stein_nnz_list, time=stein_times))
	time_data <- rbind(time_data, data.frame(algorithm='Xue et al.', nnz=xue_nnz_list, time=xue_times))
	time_data <- rbind(time_data, data.frame(algorithm='Coordinate Descent', nnz=cd_nnz_list, time=cd_times))
	time_data <- time_data[time_data$nnz <= max(stein_nnz_list),]
	gp <- ggplot(time_data, aes(x=nnz, y=time, color=algorithm, shape=algorithm)) + geom_point() + ylab('time (sec)') + xlab('number of edges') + theme(legend.title=element_blank()) + theme(text=element_text(size=16)) + scale_y_log10(breaks=10^(1:5))
  return(list(times=time_data, plot=gp))
}
