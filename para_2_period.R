library(tidyverse)
library(data.table)
library(lfe)
library(mvtnorm)
library(survival)
library(evd)
library(furrr)
library(stargazer)
library(parallel)
plan(multiprocess(workers = detectCores() - 1))

N <- 1000
f <- 20
t <- 2
b_1 <- 0.5
P <- 1
C <- 1
option_price <- 0.5

eta_set <- rgumbel(N, 0, 1)
omega_set <- rnorm(t*f, 0, 1)
epsilon_set <- rnorm(N, 0, 1)
X_set <- rnorm(N, 0, 1)
cov_mix <- matrix(c(1, 0, 0,1), ncol = 2) 
X_matrix <- rmvnorm(N, c(0, 0), cov_mix)

data <- as.data.table(seq(1, N, 1))
names(data)[1] <- "location"

omega_table <-
    as.data.table(omega_set)[, `:=` (time = ifelse(seq_len(.N) %% 2 == 1, 1, 2),
    	                             firm = ceiling(seq_len(.N)/2))]

data <-
    data[, `:=` (eta_i = eta_set[location],
				 X_i = X_matrix[location, 1],
				 epsilon_i = epsilon_set[location],
				 firm = ceiling(seq_len(.N)/ (N/f)))][
		 , well_quality := eta_i + b_1 * X_i][order(firm, -well_quality),][
		 , drill_first := ifelse(P * exp(well_quality +
		 	                     omega_set[2 * (firm - 1) + 1]) >=
		                         C + log(1 + option_price), 1, 0)][
		 , drill_second := ifelse((P * exp(well_quality +
		 	                       omega_set[2 * (firm - 1) + 2]) >=
		 	                       C) & (drill_first == 0), 1, 0)][
		 , time := ifelse(drill_first == 1,
		 	              1, ifelse(drill_second == 1, 2, 3))][
		 , status := ifelse(time <= 2, 1, 0)][
		 , y_i_t := ifelse(status == 1,
		 	               omega_set[2 * firm + time - 2] +
		 	               well_quality + epsilon_i, 0)][
		 , drilling_order := seq_len(.N) * (time <= 2), by = c("firm", "time")][
		 , firm_time := firm + time/10][
		 , firm_order := seq_len(.N) * status, by = "firm"][
		 , firm_max_order := max(firm_order), by = "firm"][
		 , firm_order := ifelse(firm_order != 0, firm_order, firm_max_order)]

ols_outcome <- felm(y_i_t ~ X_i | time, data = data, subset=(status == 1))

est_b <- coxph(Surv(firm_order, status) ~ X_i + strata(firm),
	           data = data)

print(summary(est_b))


b_hat <- as.numeric(est_b$coeff)

for (T in seq(1, t, 1)) {

	for(F in seq(1, f, 1)) {

	    sub_data <- data[time == T & firm == F]

	    total_well <- nrow(sub_data)

	    if (total_well >= 1) {

	    	for (j in seq(1, total_well, 1)) {

	    		location_id <- sub_data$location[j]

	    		good_eta <- c()

	    		good_eta <- future_map(1:2500, function(x) {
	    				eta <- rgumbel(total_well, 0, 1)

	                    sub_data$eta_trial <- eta

	                    sub_data <- sub_data[, eta_x := b_hat * X_i + eta_trial][
					                         , ordering := order(-eta_x)]

	                    if (sub_data$ordering[j] == sub_data$drilling_order[j]) {

		                    return(sub_data$eta_trial[j])
	            }
	    	})

		        data$order_stat[data$location == location_id] <-
		          mean(unlist(good_eta[lengths(good_eta) != 0]))
	        }
        }
    }
}

outcome <- felm(y_i_t ~ X_i + order_stat | firm_time, data = data,
				subset = (status == 1))



