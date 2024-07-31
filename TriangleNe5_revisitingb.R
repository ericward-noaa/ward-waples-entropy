library(EigenR)
library(ggplot2)
library(dplyr)
library(RSpectra)
library(Rcpp)
library(RcppArmadillo)

set.seed(123)
# ### this uses squared correlations and actual matings
combos = expand.grid(
  Ne = c(10,20,40,80,160,640,1280), # fine
  NOffspring = c(50, 100), # fine
  NGens1 = 6, # consistent through all simulations
  NReps = 50, # could bump to 20
  NLoci = c(25, 50, 100, 200), # limiting, (25, 50, 100, 200),
  max_loci = 500
)

# # Define the Rcpp function to calculate TC using the full matrix
# cppFunction('double calculateeigval(const arma::mat& covMatrix) {
#   arma::vec eigval = arma::eig_sym(covMatrix);
#   double sum_logeigval = sum(arma::log(eigval));
#   return sum_logeigval;
# }', depends = "RcppArmadillo")
#
# # Return log determinant
# cppFunction('double logdet(const arma::mat& covMatrix) {
#   double log_det_val;
#   double sign;
#   arma::log_det(log_det_val, sign, covMatrix);
#   return log_det_val;
# }', depends = "RcppArmadillo")
#
# # Return positive eigenvalues
# cppFunction('arma::vec calculateEigval(const arma::mat& covMatrix) {
#   arma::vec eigval = arma::eig_sym(covMatrix);
#   eigval.elem(arma::find(eigval < 0)).zeros(); // Set negative eigenvalues to zero
#   return eigval;
# }', depends = "RcppArmadillo")
#
cppFunction('
arma::mat calculateCovmat(const arma::mat& data) {
    // Centering the data (subtracting the mean)
    arma::mat centered = data.each_row() - arma::mean(data, 0);

    // Calculate the covariance matrix
    arma::mat cov = (centered.t() * centered) / (data.n_rows - 1);

    return cov;
}', depends = "RcppArmadillo")

# returns a full positive definite symmetric matrix
full_matrix = function(x) {
  x[lower.tri(x)]  <- t(x)[lower.tri(x)]
  return(x)
}

calc_det = function(x) {
  # involves 2 terms -- the first is a trace that is 0 for correlation matrix
  # this is equivalent to entropy
  x[lower.tri(x)]  <- t(x)[lower.tri(x)]
  #diag(x) <- 1
  z = as.numeric(determinant(x, logarithm = TRUE)$modulus)
  #z = log(abs(EigenR::Eigen_det(x)))
  return(z)
}

###########
reproduce <- function(parents) {
  # N = number of inds (rows)
  n_parents = dim(parents)[[1]]
  # L = number of loci (columns)
  L = dim(parents)[[2]]

  # construct matrices of the parental genotypes
  parents_geno1 = parents[parents1,]
  parents_geno2 = parents[parents2,]

  # offspring are a combination of the two parents
  offspring <-
    rbinom(n = length(parents_geno1),
           size = 1,
           p = parents_geno1 / 2) +
    rbinom(n = length(parents_geno2),
           size = 1,
           p = parents_geno2 / 2)
  # convert offspring vector back to a matrix
  offspring <-
    matrix(data = offspring,
           nrow = length(parents1),
           ncol = L)

  return(offspring)
} # end function
#############

####### #######
GetLD <- function(Geno)     {
  mat = data.matrix(Geno)
  mat = t(mat)
  # Center each variable
  mat = mat - rowMeans(mat)
  # Standardize each variable
  mat = mat / sqrt(rowSums(mat ^ 2))
  # Calculate correlations
  r = tcrossprod(mat)
  rsq = r ^ 2
  return(rsq)
} # end function

##############

#### Main program
## We are simulating a monoecious population with random selfing = original Wright-Fisher ideal population
out_df <- NULL

for (ii in 1:nrow(combos)) {
  print(ii)
  NPairs = combos$NLoci[ii] * (combos$NLoci[ii] - 1) / 2
  PairPairs = NPairs * (NPairs - 1) / 2
  bottom = combos$NLoci[ii] * (combos$NLoci[ii] - 1) / 2 - 1
  top = (combos$NLoci[ii] - 2) * (combos$NLoci[ii] - 3) / 2

  PairLoci = matrix(NA, 2, NPairs)

  counter = 0
  for (j in 1:(combos$NLoci[ii] - 1))  {
    for (k in (j + 1):combos$NLoci[ii])  {
      counter = counter + 1
      PairLoci[1, counter] = j
      PairLoci[2, counter] = k
    }
  }  ##2 end for j,k

  BigRSQ = array(NA, dim = c(combos$NLoci[ii], combos$NLoci[ii], combos$NReps[ii]))

    ##Initialize population as 100% heterozygotes so starting P = 0.5 at each locus
    Geno0 = matrix(1, combos$Ne[ii], combos$max_loci[ii])

    ##Burnin with NGens1 generations of random mating and genetic drift to equilibrate LD
    for (kk in 1:(combos$NGens1[ii] - 1))    {
      parents1 <- sample.int(n = combos$Ne[ii],
                             size = combos$Ne[ii],
                             replace = TRUE)
      parents2 <- sample.int(n = combos$Ne[ii],
                             size = combos$Ne[ii],
                             replace = TRUE)
      ## create Ne offspring in generation k, each with NLoci genotypes
      Geno1 = reproduce(Geno0)

      ## offpsring become parents of next generation
      Geno0 = Geno1
    } # end for kk

    for (jj in 1:combos$NReps[ii])  {

    parents1 <-
      sample.int(n = combos$Ne[ii],
                 size = combos$NOffspring[ii],
                 replace = TRUE)
    parents2 <-
      sample.int(n = combos$Ne[ii],
                 size = combos$NOffspring[ii],
                 replace = TRUE)
    GenoLast = reproduce(Geno0)

    if(jj == 1) {
      ## For the first replicate, sample loci that are variable. These will be
      ## used in the other replicates too

      cut = 0.05 ## Screen loci and eliminate any with MAF < 0.05
      P = colMeans(GenoLast) / 2
      loci_to_cut <- rep(0, ncol(GenoLast))
      loci_to_cut[which(P < cut)] = 1
      loci_to_cut[which(P > (1-cut))] = 1
      #GenoFinal = GenoLast[, P > cut & P < (1 - cut), drop = F]
      P = apply(GenoLast, 2, sd)
      loci_to_cut[which(P == 0)] = 1
      #GenoFinal <- GenoFinal[,which(P > 0)]
      loci_to_save <- sample(seq_len(ncol(GenoLast))[which(loci_to_cut==0)], size = combos$NLoci[ii], replace=FALSE)
      #sample.int(ncol(GenoFinal), size = combos$NLoci[ii], replace=F)
    }
    # Apply NLoci here -- so each set of replicates has a constant
    GenoFinal <- GenoLast[ , loci_to_save]
    BigRSQ[, , jj] = GetLD(GenoFinal)
  }

  rho_mat = calculateCovmat(t(matrix(BigRSQ, dim(BigRSQ)[1]*dim(BigRSQ)[1], dim(BigRSQ)[3])))
  rm(BigRSQ)

  e_cov <- eigs_sym(rho_mat, k =500, lower=TRUE)
  n_pos_cov <- length(which(e_cov$values > 0))
  log_det <- sum(log(e_cov$values[which(e_cov$values > 0)]))

  det_cov <- log_det
  misra_cov = (combos$NLoci[ii]/2) * (1+log(2*pi)) + 0.5 * det_cov
  det_cor <- NA
  misra_cor <- NA
  tc_cor <- NA
  if(combos$NLoci[ii] < 200) {
    cor_mat <- cov2cor(rho_mat)
    e_cor <- eigs_sym(cor_mat, k =500, lower=TRUE)
    n_pos_cor <- length(which(e_cor$values > 0))
    log_det <- sum(log(e_cor$values[which(e_cor$values > 0)]))
    det_cor <- log_det
    misra_cor = (combos$NLoci[ii]/2) * (1+log(2*pi)) + 0.5 * det_cor
    log_det_estimate <- det_cor * nrow(rho_mat) / n_pos_cor
    tc_cor <- det_cor - log_det_estimate
  }

    # Estimate the log determinant (approximation)
    log_det_estimate <- det_cov * nrow(rho_mat) / n_pos_cov
    tc_cov <- det_cov - log_det_estimate

  output = data.frame(
    det_covmat = det_cov,
    det_cor = det_cor,
    misra_cov = misra_cov,
    misra_cor = misra_cor,
    tc_cov = tc_cov,
    tc_cor = tc_cor,
    loci = combos$NLoci[ii],
    Ne = combos$Ne[ii],
    NOffspring = combos$NOffspring[ii]
  )

  if (is.null(out_df)) {
    out_df <- output
  } else {
    out_df <- rbind(out_df, output)
    saveRDS(out_df, file = paste0("output5_revisitedb.rds"))
  }

}


###############################################################
# This loop is for the InfRSq scenarios, and only adds a small
# chunk of code to the above loop. This is only run for the large
# Ne scenario
###############################################################
combos <- dplyr::filter(combos, Ne == max(combos$Ne))

out_df <- NULL
for (ii in 1:nrow(combos)) {
  print(ii)
  NPairs = combos$NLoci[ii] * (combos$NLoci[ii] - 1) / 2

  # Same as from Triangle NE program
  PairLoci = matrix(NA, 2, NPairs)
  counter = 0
  for (j in 1:(combos$NLoci[ii] - 1))  {
    for (k in (j + 1):combos$NLoci[ii])  {
      counter = counter + 1
      PairLoci[1, counter] = j
      PairLoci[2, counter] = k
    }
  }  ##2 end for j,k

  BigRSQ = array(NA, dim = c(combos$NLoci[ii], combos$NLoci[ii], combos$NReps[ii]))

    ##Initialize population as 100% heterozygotes so starting P = 0.5 at each locus
    Geno0 = matrix(1, combos$Ne[ii], combos$max_loci[ii])

    ##Burnin with NGens1 generations of random mating and genetic drift to equilibrate LD
    for (kk in 1:(combos$NGens1[ii] - 1))    {
      parents1 <- sample.int(n = combos$Ne[ii],
                             size = combos$Ne[ii],
                             replace = TRUE)
      parents2 <- sample.int(n = combos$Ne[ii],
                             size = combos$Ne[ii],
                             replace = TRUE)
      ## create Ne offspring in generation k, each with NLoci genotypes
      Geno1 = reproduce(Geno0)

      ## offpsring become parents of next generation
      Geno0 = Geno1
    } # end for kk

    parents1 <-
      sample.int(n = combos$Ne[ii],
                 size = combos$NOffspring[ii] * 10,
                 replace = TRUE)
    parents2 <-
      sample.int(n = combos$Ne[ii],
                 size = combos$NOffspring[ii] * 10,
                 replace = TRUE)
    GenoLast = reproduce(Geno0)

    ## Now screen loci and eliminate any with MAF < 0.05
    cut = 0.05 ## removes loci with low frequencies, or monomorphic
    P = colMeans(GenoLast) / 2
    GenoFinal = GenoLast[, P > cut & P < (1 - cut), drop = F]
    P = apply(GenoFinal, 2, sd)
    GenoFinal <- GenoFinal[,which(P > 0)]

    # Apply NLoci here -- so each set of replicates has a constant
    GenoFinal <- GenoFinal[, sample.int(ncol(GenoFinal), size = combos$NLoci[ii], replace=F)]

    GenoInf = matrix(NA, nrow(GenoFinal), ncol(GenoFinal))
    PP = colMeans(GenoFinal) / 2  ## allele frequencies in trimmed dataset
    X = c(2, 1, 0)  ## possible genotypes
    ## generate HW expected genotype frequencies
    A = PP ^ 2
    B = 2 * PP * (1 - PP)
    C = (1 - PP) ^ 2
    ## pick each genotype at each locus randomly for NOfffspring individuals
    while(1) {
      for (j in 1:ncol(GenoInf))  {
        GenoInf[, j] = sample(X,
                              combos$NOffspring[ii],
                              replace = T,
                              prob = c(A[j], B[j], C[j]))
      }  # end for j
      # keep a sample that is not monomorphic
      if(length(which(apply(GenoInf,2,sd) == 0)) == 0) break
    }

  for (jj in 1:combos$NReps[ii])  {
    BigRSQ[,,jj] = GetLD(GenoInf[sample(1:nrow(GenoInf), size=combos$NOffspring[ii], replace=T),])
  }  ## end for jj

  rho_mat = cov(t(matrix(BigRSQ, dim(BigRSQ)[1]*dim(BigRSQ)[1], dim(BigRSQ)[3])), use = "pairwise.complete.obs")
  rm(BigRSQ)

  e_cov <- eigs_sym(rho_mat, k =500, lower=TRUE)
  n_pos_cov <- length(which(e_cov$values > 0))
  log_det <- sum(log(e_cov$values[which(e_cov$values > 0)]))

  det_cov <- log_det
  misra_cov = (combos$NLoci[ii]/2) * (1+log(2*pi)) + 0.5 * det_cov
  det_cor <- NA
  misra_cor <- NA
  tc_cor <- NA
  if(combos$NLoci[ii] < 200) {
    cor_mat <- cov2cor(rho_mat)
    e_cor <- eigs_sym(cor_mat, k =500, lower=TRUE)
    n_pos_cor <- length(which(e_cor$values > 0))
    log_det <- sum(log(e_cor$values[which(e_cor$values > 0)]))
    det_cor <- log_det
    misra_cor = (combos$NLoci[ii]/2) * (1+log(2*pi)) + 0.5 * det_cor
    log_det_estimate <- det_cor * nrow(rho_mat) / n_pos_cor
    tc_cor <- det_cor - log_det_estimate
  }

  # Estimate the log determinant (approximation)
  log_det_estimate <- det_cov * nrow(rho_mat) / n_pos_cov
  tc_cov <- det_cov - log_det_estimate

  output = data.frame(
    det_covmat = det_cov,
    det_cor = det_cor,
    misra_cov = misra_cov,
    misra_cor = misra_cor,
    tc_cov = tc_cov,
    tc_cor = tc_cor,
    loci = combos$NLoci[ii],
    Ne = combos$Ne[ii],
    NOffspring = combos$NOffspring[ii]
  )

  if (is.null(out_df)) {
    out_df <- output
  } else {
    out_df <- rbind(out_df, output)
    saveRDS(out_df, file = paste0("output_inf_v5_revistedb.rds"))
  }

}


# bind data together
out <- readRDS("output5_revisitedb.rds")
inf <- readRDS("output_inf_v5_revistedb.rds")
out <- dplyr::select(out, -det_cor, -misra_cor, -tc_cor)
inf <- dplyr::select(inf, -det_cor, -misra_cor, -tc_cor)

names(inf)[1:3] = paste0("inf_", names(inf)[1:3])
out <- dplyr::left_join(out, dplyr::select(inf, -Ne)) # add entropy stats
write.csv(out, "all_data_combined_Dec13.csv", row.names = FALSE)

out$Ne <- as.factor(out$Ne)
p1 <- out %>%
ggplot(aes(loci, tc_cov, col = Ne)) +
  facet_wrap(~NOffspring) +
  geom_point(alpha=0.4) +
  geom_line(alpha=0.4) +
  theme_bw() +
  scale_x_log10() +  # Apply log scale to x-axis
  scale_y_log10() +
  xlab("") + ylab("Total correlation (TC)") +
  scale_color_viridis_d(end=0.8, option="magma") +
  theme(strip.background =element_rect(fill="white"))
p2 <- out %>%
  ggplot(aes(loci, -misra_cov, col = Ne)) +
  facet_wrap(~NOffspring) +
  geom_point(alpha=0.4) +
  geom_line(alpha=0.4) +
  theme_bw() +
  scale_x_log10() +  # Apply log scale to x-axis
  scale_y_log10() +
  xlab("Loci (L)") + ylab("-H(x)") +
  scale_color_viridis_d(end=0.8, option="magma") +
  theme(strip.background =element_rect(fill="white"))

library(patchwork)
# Combine the plots
combined_plot <- p1 + p2 +
  plot_layout(ncol = 1, guides = "collect") &
  theme(legend.position = "right")

ggsave("Fig_1_loci_vs_entropy.png")

p1 <- out %>%
  ggplot(aes(loci, tc_cov/inf_tc_cov, col = Ne)) +
  facet_wrap(~NOffspring) +
  geom_point(alpha=0.4) +
  geom_line(alpha=0.4) +
  scale_x_log10() +
  theme_bw() +
  xlab("") + ylab(expression(paste("TC / TC"[infinity]))) +
  scale_color_viridis_d(end=0.8, option="magma") +
  theme(strip.background =element_rect(fill="white"))
p2 <- out %>%
  ggplot(aes(loci, misra_cov/inf_misra_cov, col = Ne)) +
  facet_wrap(~NOffspring) +
  geom_point(alpha=0.4) +
  geom_line(alpha=0.4) +
  scale_x_log10() +
  theme_bw() +
  xlab("Loci (L)") + ylab(expression(paste("H(x) / H(x)"[infinity]))) +
  scale_color_viridis_d(end=0.8, option="magma") +
  theme(strip.background =element_rect(fill="white"))

combined_plot <- p1 + p2 +
  plot_layout(ncol = 1, guides = "collect") &
  theme(legend.position = "right")
ggsave("Fig_2_entropy_inf.png")

out$Ne <- as.numeric(as.character(out$Ne))


df <- expand.grid(dim = seq(10,100,by=10), i = 1:10)
df$tc_cov <- 0
df$tc_cor <- 0
df$misra_cov <- 0
df$misra_cor <- 0
for(i in 1:nrow(df)) {
  set.seed(i)
  m = diag(runif(df$dim[i], 1, 3))
  det_cov <- calc_det(m)
  det_cor <- calc_det(cov2cor(m))
  df$misra_cov[i] = (df$dim[i]/2) * (1+log(2*pi)) + 0.5 * det_cov
  df$misra_cor[i] = (df$dim[i]/2) * (1+log(2*pi)) + 0.5 * det_cor
  df$tc_cov[i] <- sum(log(eigen(m)$values)) - det_cov
  df$tc_cor[i] <- sum(log(eigen(cov2cor(m))$values)) - det_cor
}

# Explore relationship with entropy and random numbers
df <- data.frame(dim = seq(10,70,by=5),
                 entropy = NA,
                 misra_cov = NA, tc_cov = NA)
for(i in 1:nrow(df)) {
  # simulate a random array
  arr <- array(rnorm(df$dim[i]*df$dim[i]*20,0,1), dim = c(df$dim[i], df$dim[i], 20))
  #
  covMat <- cov(t(matrix(arr, dim(arr)[1]*dim(arr)[1], dim(arr)[3])), use = "pairwise.complete.obs")
  det_cov <- logdet(covMat)
  df$misra_cov[i] = (df$dim[i]/2) * (1+log(2*pi)) + 0.5 * det_cov
  eigval <- calculateEigval(covMat)
  df$tc_cov[i] <- sum(log(eigval[which(eigval>0)])) - det_cov
}

# linear models
out <- readRDS("output5_revisitedb.rds")
out$y <- log(out$tc_cov)
library(glmmTMB)
fit <- glmmTMB(y ~ log10(loci) + log10(Ne) + log10(NOffspring), data = out)
lm_fit <- lm(y ~ log(loci) + log(Ne) + log(NOffspring), data = out)

out$pred <- predict(lm_fit)
out$Ne <- as.factor(out$Ne)
out %>%
  ggplot(aes(log(loci), pred, col = Ne)) +
  facet_wrap(~NOffspring) +
  geom_point(alpha=0.4) +
  geom_line(alpha=0.4) +
  theme_bw() +
  xlab("Ln loci") + ylab("Ln Entropy") +
  scale_color_viridis_d(end=0.8, option="magma") +
  theme(strip.background =element_rect(fill="white"))
ggsave("Figure_S1.png")



out %>%
  ggplot(aes(log(loci), pred, col = Ne)) +
  facet_wrap(~NOffspring) +
  geom_point(alpha=0.4) +
  geom_line(alpha=0.4) +
  theme_bw() +
  xlab("Ln loci") + ylab("Ln Entropy") +
  scale_color_viridis_d(end=0.8, option="magma") +
  theme(strip.background =element_rect(fill="white"))
