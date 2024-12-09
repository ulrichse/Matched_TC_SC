
library(data.table)
library(dplyr)
library(arrow)
library(lme4)
library(sjPlot)
setwd("C:/Users/ulric/OneDrive - University of North Carolina at Chapel Hill/RStudio/Matched-Analysis/TC-Exposure/Hurricane_Maternal_Health/")

dat <- read.csv("Hurricane_Maternal_Health/final_merged_data.csv") %>%
  mutate(date=as.Date(date))%>%
  mutate(year=year(date),
         doy=yday(date),
         dow=wday(date))%>%
  setorder(county, date)

dat[is.na(dat)]<-0

matched_df <- read_parquet("Matched_DF_SC_TC_30day_SU.parquet")

# Function to create crossbasis for lags 0 to 7
matched_cb_lag30 <- function(data, n_controls = 3, lag_range = 0:30, control_doy_range = -3:3){
  
  county_list <- unique(dat$county)
  setorder(dat, county, date)
  
  # Use "dlnm" package to generate the distributed lag function for "disaster"
  for (i in 1:length(county_list)) {
    orig_dat <- subset(dat, county == county_list[i])
    match_dat <- subset(matched_df, county == county_list[i])
    
    orig_cb <- dlnm::crossbasis(orig_dat$disaster, lag = c(0, 30),
                                argvar = list(fun = "lin"),
                                arglag = list(fun = "integer"))
    obs_n <- nrow(orig_dat)
    orig_cb_matr <- as.data.frame(subset(orig_cb, nrow = obs_n))
    orig_cb_matr$date <- orig_dat$date
    matched_date <- match_dat %>% dplyr::select(date)
    matched_cb_matr <- orig_cb_matr %>%
      dplyr::right_join(matched_date, by = "date") %>%
      dplyr::select(-date) %>% as.matrix()
    
    if (i == 1) {
      matched_cb_matrix <- matched_cb_matr
    } else {
      matched_cb_matrix <- rbind(matched_cb_matrix, matched_cb_matr)
    }
    # Add attributes
    matched_dim <- dim(matched_cb_matrix)
    attr <- attributes(orig_cb)
    attr$dim <- matched_dim
    matched_cb <- matched_cb_matrix
    attributes(matched_cb) <- attr
  }
  
  return(matched_cb)
  
  gc()
  
}

# Create crossbasis
matched_cb <- matched_cb_lag30(dat)

# alternative model 1
fit_am1 <- glm(smm20  ~ matched_cb + factor(dow) + factor(year) + factor(county)+factor(stratum),
               data = matched_df, 
               family = poisson(link = "log"))
pred_am1 <- dlnm::crosspred(matched_cb, fit_am1, at = 1)
plot(pred_am1 , "slices", var=1, ci="bars")

# we used "delta method" to calculate confidence interval for the overall 
# RR for the entire heatwave period

pred_am1$matRRfit
pred_am1$matRRlow; pred_am1$matRRhigh

# estimate of overall RR for the entire heatwave-exposure period (i.e., ten days)
over_rr <- sum(pred_am1$matRRfit)/8

# we used "delta method" to calculate confidence interval for the overall 
# RR for the entire heatwave period
library(msm)
estvar <- pred_am1$vcov
estmean <- c(pred_am1$coefficients)

over_rr_se <- msm::deltamethod(~ (exp(x1) + exp(x2) + exp(x3) + exp(x4) + 
                                    exp(x5) + exp(x6) + exp(x7) + 
                                    exp(x8))/8, 
                               estmean, estvar)
over_rr_low <- over_rr / exp(1.96*over_rr_se)
over_rr_high <- over_rr * exp(1.96*over_rr_se)

# get model results

outcomes <- c("smm20", "smm21", "dep", "preterm", "lbw", "gest_htn", "htn", "mental_disorders","gest_diabetes", "gen_anxiety")

get_daily_lags <- function(matched_df, matched_cb) {
  # Initialize lists to store results
  datasets <- list()
  combined_data_list <- list()
  
  for (outcome in outcomes) {
    # Build formula dynamically
    formula <- as.formula(
      paste0(outcome, " ~ matched_cb + factor(dow) + factor(year) + factor(county) + factor(stratum)")
    )
    
    # Fit GLM and predict using DLNM
    fit_model <- glm(formula, data = matched_df, family = poisson(link = "log"))
    predictions <- dlnm::crosspred(matched_cb, fit_model, at = 1)
    
    # Store predictions in datasets
    datasets[[outcome]] <- predictions
    
    # Print prediction matrices (optional, consider removing if not needed)
    print(predictions$matRRfit)
    print(predictions$matRRlow)
    print(predictions$matRRhigh)
  }
  
  # Determine max number of columns across all matrices
  max_columns <- max(sapply(datasets, function(pred) {
    max(ncol(pred$matRRfit), ncol(pred$matRRlow), ncol(pred$matRRhigh))
  }))
  
  # Create combined data frame
  for (outcome in names(datasets)) {
    pred <- datasets[[outcome]]
    
    lag_names <- paste0("lag", 0:(ncol(pred$matRRfit) - 1))
    # Consolidate prediction matrices into a long-format data frame
    combined_data <- rbind(
      cbind(pred$matRRfit, Type = "matRRfit", Pred_Model = outcome, Lag = lag_names),
      cbind(pred$matRRlow, Type = "matRRlow", Pred_Model = outcome, Lag = lag_names),
      cbind(pred$matRRhigh, Type = "matRRhigh", Pred_Model = outcome, Lag = lag_names)
    )
    
    # Assign to combined_data_list with outcome-specific name
    combined_data_list[[paste0(outcome)]] <- combined_data
  }
  
  # Return the list of combined data
  return(combined_data_list)
}

cumulative_lag0_lag7 <- function(matched_df, matched_cb){
  
  outcome_data <- list()
  
  for (outcome in outcomes) {
    # Build formula dynamically
    formula <- as.formula(
      paste0(outcome, " ~ matched_cb + factor(dow) + factor(year) + factor(county) + factor(stratum)")
    )
    
    # Fit GLM and predict using DLNM
    fit_model <- glm(formula, data = matched_df, family = poisson(link = "log"))
    pred <- dlnm::crosspred(matched_cb, fit_model, at = 1)
    
    over_rr <- sum(pred$matRRfit) / 8
    
    library(msm)
    estvar <- pred$vcov
    estmean <- c(pred$coefficients)
    
    over_rr_se <- msm::deltamethod(~ (exp(x1) + exp(x2) + exp(x3) + exp(x4) + exp(x5) + exp(x6) + exp(x7) + exp(x8)) / 8, 
                                   estmean, estvar)
    over_rr_low <- over_rr / exp(1.96 * over_rr_se)
    over_rr_high <- over_rr * exp(1.96 * over_rr_se)
    
    # Create a data frame for the current outcome
    outcome_df <- data.frame(outcome = outcome,
                             over_rr = over_rr,
                             over_rr_se = over_rr_se,
                             over_rr_low = over_rr_low,
                             over_rr_high = over_rr_high,
                             estvar = estvar,
                             estmean = estmean)
    
    outcome_data[[paste0(outcome, "_cumlag")]] <- outcome_df
  }
  
  data_frames <- unname(outcome_data)
  merged_df <- do.call(rbind, data_frames)
  return(merged_df)
  
  gc()
  
}

#combined_data <- get_daily_lags(matched_df, matched_cb)
#write.csv(combined_data, file = "daily_lag_all_outcomes_30d.csv", row.names = FALSE)

merged_df <- cumulative_lag0_lag7(matched_df, matched_cb)
merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
write.csv(merged_df, "lag0_lag7_cumulative_all_outcomes_30d.csv")
gc()

# Email from Burrows: 

#Without looking at the data or results, my guess is that the crossbasis is changing the magnitude of the effect because 
#you’re accounting for lagged effects, not just looking at an immediate effect.
#You may want to consider looking at non-linear effects (changing argvar) as a sensitivity anlaysis
#In previous studies, my team has looked separately at single counties or zips and then combined them to estimate a relative risk. 
# I think the way you’ve set it up is fine, but if you want to capture local variations or allow for location-specific flexibility 
#that could be an option. However, depending on the sample size in each zip you may need to run the estimates in the pooled way that 
#you have them set up now.
#You may find it helpful to look at plot(pred_am5, "slices", var=1, ci="bars") to visualize the effects over the lagged period.
#Depending on how you did the matching, you will probably need to include a strata term to account for the matched groups in your model.
#The eliminate = factor(zip) term in your gnm function is effectively controlling for zip code-specific baseline differences. I think it 
#might be redundant to have that eliminate term in there as well as the interaction between zip and year. I’d consider removing.





