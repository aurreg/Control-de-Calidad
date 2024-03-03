winzordraze <- function(data, all_sample, q) {
  num_values <- length(data)
  num_replace <- round(q * num_values)
  
  start_replace <- min(num_replace, num_values)
  end_replace <- max(num_values - num_replace + 1, 1)
  
  data[1:start_replace] <- mean(all_sample)
  data[end_replace:num_values] <- mean(all_sample)
  
  return(data)
}

mixture <- function(n, phi) {
  probabilities <- runif(n)
  x <- numeric(n)
  for (i in 1:n) {
    x[i] <- ifelse(probabilities[i] > phi, rnorm(1), rnorm(1) + rchisq(1, df = 1))
  }
  return(x)
}



n<-10

m<-200
matriz_sigma <- matrix(nrow = m, ncol = n)
Ls <- seq(5.2, 5.3, by = 0.01)

results_Winsorize <- data.frame(L = numeric(length(Ls)),
                                Mean_ARL = numeric(length(Ls)),
                                SD_ARL = numeric(length(Ls)),
                                Quantile_10 = numeric(length(Ls)),
                                Quantile_25 = numeric(length(Ls)),
                                Quantile_50 = numeric(length(Ls)),
                                Quantile_75 = numeric(length(Ls)),
                                Quantile_90 = numeric(length(Ls)),
                                ARL_d_Risk = numeric(length(Ls)))

VLR<-NULL

repeticiones<-10000

for (i in 1:length(Ls)){
  L<-Ls[i]
  
  print(L)  
  
  for (j in 1:repeticiones) {
    
    phi<-0
    
    
    for (k in 1:m) {
      matriz_sigma[k, ] <- mixture(n, phi)
    }
    
    sigma_barra_W <- mean(apply(matriz_sigma, 1, function(x) {
      filtered_data <- winzordraze(x, q=0.1, all_sample = as.vector(matriz_sigma))
      var(filtered_data)
    }))
    
    
    
    VLR[j] <- 1/(1-(pchisq((((n-1)*sigma_barra_W)* (1 + (L * (sqrt(2 / (n - 1)))))),n-1)))
  }
  
  ARL_l <- ARL_0 - (E * SDRL_0)
  ARL_u <- ARL_0 + (E * SDRL_0)
  ARL_d_risk_winsorize <- mean(ARL_l <= VLR & VLR <= ARL_u)
  
  results_Winsorize[i, "L"] <- L
  results_Winsorize[i, "Mean_ARL"] <- mean(VLR)
  results_Winsorize[i, "SD_ARL"] <- sd(VLR)
  results_Winsorize[i, "Quantile_10"] <- quantile(VLR, probs = 0.1)
  results_Winsorize[i, "Quantile_25"] <- quantile(VLR, probs = 0.25)
  results_Winsorize[i, "Quantile_50"] <- quantile(VLR, probs = 0.5)
  results_Winsorize[i, "Quantile_75"] <- quantile(VLR, probs = 0.75)
  results_Winsorize[i, "Quantile_90"] <- quantile(VLR, probs = 0.9)
  results_Winsorize[i, "ARL_d_Risk"] <- 1-ARL_d_risk_winsorize
  
} 
print(results_Winsorize)   

set.seed(13)

n<-10

m<-200

L=3.799594

UCL <- function(sigma, L, n) {
  UCL <- sigma * (1 + (L * (sqrt(2 / (n - 1)))))
  return(UCL)
}


library(DescTools)


f_Tukey <- function(data, eta, all_sample = NULL) {
  if (is.null(all_sample)) {
    Q1 <- quantile(data, 0.25)
    Q3 <- quantile(data, 0.75)
  } else {
    Q1 <- quantile(all_sample, 0.25)
    Q3 <- quantile(all_sample, 0.75)
  }
  
  LDLTukey <- Q1 - eta * (Q3 - Q1)
  UDLTukey <- Q3 + eta * (Q3 - Q1)
  
  filtered_data <- data[data >= LDLTukey & data <= UDLTukey]
  
  return(filtered_data)
}
f_MAD <- function(data, eta, all_sample) {
  median_val <- median(all_sample)
  MAD <- median(abs(all_sample - median_val))
  
  LDLMAD <- median_val - eta * MAD/0.6745
  UDLMAD <- median_val + eta * MAD/0.6745
  
  filtered_data <- data[data >= LDLMAD & data <= UDLMAD]
  
  return(filtered_data)
}

f_Zscore <- function(data, confidence_level, all_sample) {
  mean_val <- mean(all_sample)
  sd_val <- sd(all_sample)
  
  z_score_threshold <- qnorm(1 - (1 - confidence_level) / 2)
  
  lower_limit <- mean_val - z_score_threshold * sd_val
  upper_limit <- mean_val + z_score_threshold * sd_val
  
  filtered_data <- data[data >= lower_limit & data <= upper_limit]
  
  return(filtered_data)
}

winzordraze <- function(data, all_sample, q) {
  num_values <- length(data)
  num_replace <- round(q * num_values)
  
  start_replace <- min(num_replace, num_values)
  end_replace <- max(num_values - num_replace + 1, 1)
  
  data[1:start_replace] <- mean(all_sample)
  data[end_replace:num_values] <- mean(all_sample)
  
  return(data)
}

mixture <- function(n, phi) {
  probabilities <- runif(n)
  x <- numeric(n)
  for (i in 1:n) {
    x[i] <- ifelse(probabilities[i] > phi, rnorm(1), rnorm(1) + rchisq(1, df = 1))
  }
  return(x)
}



etaM=3.642245

etaT=2.2

sigma_barra<-NULL

sigma<-0


LR<-0

VLR<-NULL
VLR1<-NULL
VLR2<-NULL
VLR3<-NULL
VLR4<-NULL

phi_values <- seq(0, 0.1, by = 0.01)

results <- data.frame(Phi = numeric(length(phi_values)),
                      Mean_ARL = numeric(length(phi_values)),
                      SD_ARL = numeric(length(phi_values)),
                      Quantile_10 = numeric(length(phi_values)),
                      Quantile_25 = numeric(length(phi_values)),
                      Quantile_50 = numeric(length(phi_values)),
                      Quantile_75 = numeric(length(phi_values)),
                      Quantile_90 = numeric(length(phi_values)),
                      ARL_d_Risk = numeric(length(phi_values)))

results_Tukey <- data.frame(Phi = numeric(length(phi_values)),
                            Mean_ARL = numeric(length(phi_values)),
                            SD_ARL = numeric(length(phi_values)),
                            Quantile_10 = numeric(length(phi_values)),
                            Quantile_25 = numeric(length(phi_values)),
                            Quantile_50 = numeric(length(phi_values)),
                            Quantile_75 = numeric(length(phi_values)),
                            Quantile_90 = numeric(length(phi_values)),
                            ARL_d_Risk = numeric(length(phi_values)))

results_MAD <- data.frame(Phi = numeric(length(phi_values)),
                          Mean_ARL = numeric(length(phi_values)),
                          SD_ARL = numeric(length(phi_values)),
                          Quantile_10 = numeric(length(phi_values)),
                          Quantile_25 = numeric(length(phi_values)),
                          Quantile_50 = numeric(length(phi_values)),
                          Quantile_75 = numeric(length(phi_values)),
                          Quantile_90 = numeric(length(phi_values)),
                          ARL_d_Risk = numeric(length(phi_values)))

results_zscore <- data.frame(Phi = numeric(length(phi_values)),
                             Mean_ARL = numeric(length(phi_values)),
                             SD_ARL = numeric(length(phi_values)),
                             Quantile_10 = numeric(length(phi_values)),
                             Quantile_25 = numeric(length(phi_values)),
                             Quantile_50 = numeric(length(phi_values)),
                             Quantile_75 = numeric(length(phi_values)),
                             Quantile_90 = numeric(length(phi_values)),
                             ARL_d_Risk = numeric(length(phi_values)))


results_Winsorize <- data.frame(Phi = numeric(length(phi_values)),
                                Mean_ARL = numeric(length(phi_values)),
                                SD_ARL = numeric(length(phi_values)),
                                Quantile_10 = numeric(length(phi_values)),
                                Quantile_25 = numeric(length(phi_values)),
                                Quantile_50 = numeric(length(phi_values)),
                                Quantile_75 = numeric(length(phi_values)),
                                Quantile_90 = numeric(length(phi_values)),
                                ARL_d_Risk = numeric(length(phi_values)))

repeticiones <- 100000

matriz_sigma <- matrix(nrow = m, ncol = n)

for (i in 1:length(phi_values)){
  
  phi<-phi_values[i]
  
  print(phi)  
  
  for (j in 1:repeticiones) {
    
    
    
    
    for (k in 1:m) {
      matriz_sigma[k, ] <- mixture(n, phi)
    }
    
    
    sigma_barra <- mean(apply(matriz_sigma,1,var))
    
    sigma_barra_T <- mean(apply(matriz_sigma, 1, function(x) {
      filtered_data <- f_Tukey(x, etaT, all_sample = as.vector(matriz_sigma))
      var(filtered_data)
    }))
    
    sigma_barra_M <- mean(apply(matriz_sigma, 1, function(x) {
      filtered_data <- f_MAD(x, etaM, all_sample = as.vector(matriz_sigma))
      var(filtered_data)
    }))
    
    sigma_barra_Z <- mean(apply(matriz_sigma, 1, function(x) {
      filtered_data <- f_Zscore(x, 0.9999, all_sample = as.vector(matriz_sigma))
      var(filtered_data)
    }))
    sigma_barra_W <- mean(apply(matriz_sigma, 1, function(x) {
      filtered_data <- winzordraze(x, q=0.1, all_sample = as.vector(matriz_sigma))
      var(filtered_data)
    }))
    
    
    
    VLR[j] <- 1/(1-(pchisq((((n-1)*sigma_barra)* (1 + (L * (sqrt(2 / (n - 1)))))),n-1)))
    VLR1[j] <- 1/(1-(pchisq((((n-1)*sigma_barra_T)* (1 + (L * (sqrt(2 / (n - 1)))))),n-1)))
    VLR2[j] <- 1/(1-(pchisq((((n-1)*sigma_barra_M)* (1 + (L * (sqrt(2 / (n - 1)))))),n-1)))
    VLR3[j] <- 1/(1-(pchisq((((n-1)*sigma_barra_Z)* (1 + (L * (sqrt(2 / (n - 1)))))),n-1)))
    VLR4[j] <- 1/(1-(pchisq((((n-1)*sigma_barra_W)* (1 + (5.27 * (sqrt(2 / (n - 1)))))),n-1)))
  }
  
  
  
  ARL_0 <- 370.37     
  SDRL_0 <- 370.37     
  E <- 1/4            # Constante E
  ARL_l <- ARL_0 - (E * SDRL_0)
  ARL_u <- ARL_0 + (E * SDRL_0)
  ARL_d_risk <- mean(ARL_l <= VLR & VLR <= ARL_u)
  
  results[i, "Phi"] <- phi
  results[i, "Mean_ARL"] <- mean(VLR)
  results[i, "SD_ARL"] <- sd(VLR)
  results[i, "Quantile_10"] <- quantile(VLR, probs = 0.1)
  results[i, "Quantile_25"] <- quantile(VLR, probs = 0.25)
  results[i, "Quantile_50"] <- quantile(VLR, probs = 0.5)
  results[i, "Quantile_75"] <- quantile(VLR, probs = 0.75)
  results[i, "Quantile_90"] <- quantile(VLR, probs = 0.9)
  results[i, "ARL_d_Risk"] <- 1-ARL_d_risk
  
  
  ARL_d_risk_Tukey <- mean(ARL_l <= VLR1 & VLR1 <= ARL_u)
  
  results_Tukey[i, "Phi"] <- phi
  results_Tukey[i, "Mean_ARL"] <- mean(VLR1)
  results_Tukey[i, "SD_ARL"] <- sd(VLR1)
  results_Tukey[i, "Quantile_10"] <- quantile(VLR1, probs = 0.1)
  results_Tukey[i, "Quantile_25"] <- quantile(VLR1, probs = 0.25)
  results_Tukey[i, "Quantile_50"] <- quantile(VLR1, probs = 0.5)
  results_Tukey[i, "Quantile_75"] <- quantile(VLR1, probs = 0.75)
  results_Tukey[i, "Quantile_90"] <- quantile(VLR1, probs = 0.9)
  results_Tukey[i, "ARL_d_Risk"] <- 1-ARL_d_risk_Tukey
  
  ARL_d_risk_MAD <- mean(ARL_l <= VLR2 & VLR2 <= ARL_u)
  
  results_MAD[i, "Phi"] <- phi
  results_MAD[i, "Mean_ARL"] <- mean(VLR2)
  results_MAD[i, "SD_ARL"] <- sd(VLR2)
  results_MAD[i, "Quantile_10"] <- quantile(VLR2, probs = 0.1)
  results_MAD[i, "Quantile_25"] <- quantile(VLR2, probs = 0.25)
  results_MAD[i, "Quantile_50"] <- quantile(VLR2, probs = 0.5)
  results_MAD[i, "Quantile_75"] <- quantile(VLR2, probs = 0.75)
  results_MAD[i, "Quantile_90"] <- quantile(VLR2, probs = 0.9)
  results_MAD[i, "ARL_d_Risk"] <- 1-ARL_d_risk_MAD
  
  ARL_d_risk_zscore <- mean(ARL_l <= VLR3 & VLR3 <= ARL_u)
  
  results_zscore[i, "Phi"] <- phi
  results_zscore[i, "Mean_ARL"] <- mean(VLR3)
  results_zscore[i, "SD_ARL"] <- sd(VLR3)
  results_zscore[i, "Quantile_10"] <- quantile(VLR3, probs = 0.1)
  results_zscore[i, "Quantile_25"] <- quantile(VLR3, probs = 0.25)
  results_zscore[i, "Quantile_50"] <- quantile(VLR3, probs = 0.5)
  results_zscore[i, "Quantile_75"] <- quantile(VLR3, probs = 0.75)
  results_zscore[i, "Quantile_90"] <- quantile(VLR3, probs = 0.9)
  results_zscore[i, "ARL_d_Risk"] <-1- ARL_d_risk_zscore
  
  
  ARL_d_risk_winsorize <- mean(ARL_l <= VLR4 & VLR4 <= ARL_u)
  
  results_Winsorize[i, "Phi"] <- phi
  results_Winsorize[i, "Mean_ARL"] <- mean(VLR4)
  results_Winsorize[i, "SD_ARL"] <- sd(VLR4)
  results_Winsorize[i, "Quantile_10"] <- quantile(VLR4, probs = 0.1)
  results_Winsorize[i, "Quantile_25"] <- quantile(VLR4, probs = 0.25)
  results_Winsorize[i, "Quantile_50"] <- quantile(VLR4, probs = 0.5)
  results_Winsorize[i, "Quantile_75"] <- quantile(VLR4, probs = 0.75)
  results_Winsorize[i, "Quantile_90"] <- quantile(VLR4, probs = 0.9)
  results_Winsorize[i, "ARL_d_Risk"] <- 1-ARL_d_risk_winsorize
  
}

print(results)
print(results_Tukey)
print(results_MAD)
print(results_zscore)
print(results_Winsorize)
