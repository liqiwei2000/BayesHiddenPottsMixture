# ***README***
# The following script is used to load the marked spatial point pattern data analyzed in
# the submitted manuscript titled "A Bayesian Hidden Potts Mixture Model for Analyzing 
# Lung Cancer Pathology Images."

# Please only run the code of one block each time. The outputs are a data.frame object
# with three columns: x, y, and z, which denote the x, y coordinates, and marks of each 
# point. Marks should be labelled as natural numbers 1, 2, 3...
# ***END***



# ========================================================================================
# ========================================================================================
# Load simulated data, where the points are generated from a homogeneous Poisson process
# and Delta shown in Figure 2(a)
# ========================================================================================
# ========================================================================================
data <- read.table("simulated_data_homogeneous_poisson_scenario_1.txt", header = TRUE);
par(pty = "s", mar=c(4.5, 2.5, 1, 1));
plot(data$x, data$y, xlim = c(1, 50), ylim = c(1, 50), col = data$z, pch = data$z, 
     main = "", xlab = "x", ylab = "y", cex.lab = 1.5, cex.axis = 1.5, cex = 0.5);



# ========================================================================================
# ========================================================================================
# Load simulated data, where the points are generated from a homogeneous Poisson process
# and Delta shown in Figure 2(b)
# ========================================================================================
# ========================================================================================
data <- read.table("simulated_data_homogeneous_poisson_scenario_2.txt", header = TRUE);
par(pty = "s", mar=c(2.5, 2.5, 1, 1));
plot(data$x, data$y, xlim = c(1, 50), ylim = c(1, 50), col = data$z, pch = data$z, 
     main = "", xlab = "x", ylab = "y", cex.lab = 1.5, cex.axis = 1.5, cex = 0.5);



# ========================================================================================
# ========================================================================================
# Load simulated data, where the points are generated from a log Gaussian Cox process and
# Delta shown in Figure 2(a)
# ========================================================================================
# ========================================================================================
data <- read.table("simulated_data_log_gaussian_cox_scenario_1.txt", header = TRUE);
par(pty = "s", mar=c(2.5, 2.5, 1, 1));
plot(data$x, data$y, xlim = c(1, 50), ylim = c(1, 50), col = data$z, pch = data$z, 
     main = "", xlab = "x", ylab = "y", cex.lab = 1.5, cex.axis = 1.5, cex = 0.5);



# ========================================================================================
# ========================================================================================
# Load simulated data, where the points are generated from a log Gaussian Cox process and
# Delta shown in Figure 2(b)
# ========================================================================================
# ========================================================================================
data <- read.table("simulated_data_log_gaussian_cox_scenario_2.txt", header = TRUE);
par(pty = "s", mar=c(2.5, 2.5, 1, 1));
plot(data$x, data$y, xlim = c(1, 50), ylim = c(1, 50), col = data$z, pch = data$z, 
     main = "", xlab = "x", ylab = "y", cex.lab = 1.5, cex.axis = 1.5, cex = 0.5);



# ========================================================================================
# ========================================================================================
# Load NLST data, as shown in Figure 4(a)
# ========================================================================================
# ========================================================================================
data <- read.table("real_data_fig_4a.txt", header = TRUE); # Q = 3, n = 10533
par(pty = "s", mar=c(2.5, 2.5, 1, 1));
plot(data$x, data$y, xlim = c(0, 1), ylim = c(0, 1), col = data$z, pch = data$z, 
     main = "", xlab = "x", ylab = "y", cex.lab = 1.5, cex.axis = 1.5, cex = 0.5);



# ========================================================================================
# ========================================================================================
# Load NLST data, as shown in Figure 4(b)
# ========================================================================================
# ========================================================================================
data <- read.table("real_data_fig_4b.txt", header = TRUE); # Q = 3, n = 13874
par(pty = "s", mar=c(2.5, 2.5, 1, 1));
plot(data$x, data$y, xlim = c(0, 1), ylim = c(0, 1), col = data$z, pch = data$z, 
     main = "", xlab = "x", ylab = "y", cex.lab = 1.5, cex.axis = 1.5, cex = 0.5);


