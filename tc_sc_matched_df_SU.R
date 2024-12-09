
library(dplyr)
library(data.table)
library(arrow)
getwd()
setwd("C:/Users/ulric/OneDrive - University of North Carolina at Chapel Hill/RStudio/Matched-Analysis/TC-Exposure/Hurricane_Maternal_Health/")

dat <- read.csv("Hurricane_Maternal_Health/final_merged_data.csv") %>%
  mutate(date=as.Date(date))%>%
  mutate(year=year(date),
         doy=yday(date),
         dow=wday(date))%>%
  setorder(county, date)

dat[is.na(dat)]<-0

no_tc <- dat %>%
  group_by(county)%>%
  summarize(tc=sum(disaster))%>% # Create list of counties without heatwaves
  filter(tc==0)%>%
  dplyr::select(county)

dat <- dat %>% # Remove counties without heatwaves from the data
  filter(!county %in% no_tc$county)%>%
  setorder(county, date)

gc()



##### Match TC-exposed days to comparable unexposed days ######

#### For lag0 to lag7 ####

set.seed(123)
county_list <- unique(dat$county)

for(i in 1:length(county_list)){
  df <- subset(dat, county == county_list[i])
  
  df$time <- 1:nrow(df)
  disaster_idx <- which(df$disaster == 1)
  
  # Handle edge cases for index boundaries
  cand_control_idx <- unique(c(disaster_idx, disaster_idx + 1, disaster_idx - 1))
  cand_control_idx <- cand_control_idx[cand_control_idx > 0 & cand_control_idx <= nrow(df)] # Ensure valid indices
  
  df$cand_control <- TRUE
  df$cand_control[cand_control_idx] <- FALSE # Mark disaster days and adjacent as non-controls
  
  case_dates <- subset(df, disaster == 1)
  control_dates <- subset(df, cand_control == TRUE) 
  
  for(j in 1:nrow(case_dates)){
    # choose lags (lagged 0 to lagged 7)
    lag_dates <- case_dates[j, ]$date + 0:30
    lag_case <- subset(df, date %in% lag_dates)
    
    # choose 10 comparable unexposed days for each disaster-exposed day
    control_range <- case_dates[j, ]$doy + -3:3
    control_subset <- subset(control_dates,
                             control_dates$year != case_dates[j, ]$year &
                               doy %in% control_range)
    
    # Check if there are enough controls
    if(nrow(control_subset) < 3){
      next # Skip this iteration if there are not enough controls
    }
    
    controls <- dplyr::sample_n(control_subset, 3)
    
    # choose lagged days for selected unexposed days
    la_con <- c(1:30)
    lag_control <- NULL
    for(p in 1:length(la_con)){
      lag_control_dates <- controls$date + la_con[p]
      lag_control_each <- subset(df, date %in% lag_control_dates & cand_control == TRUE) # Modified
      
      if(nrow(lag_control_each) == 0){
        next # Skip this lag if no valid control days
      }
      
      if(is.null(lag_control)){
        lag_control <- lag_control_each
      }else{
        lag_control <- rbind(lag_control, lag_control_each)
      }
    }
    
    # Combine lag_case, controls, and lagged controls
    j_stratum <- rbind(lag_case, controls, lag_control)
    
    # Dynamically set the 'status' and 'lag' vectors based on actual data size
    num_lag_case <- nrow(lag_case)
    num_controls <- nrow(controls)
    num_lag_control <- nrow(lag_control)
    
    stratum <- paste("stratum", j, sep = ".")
    j_stratum$stratum <- stratum
    
    status <- c(rep("case", num_lag_case), rep("control", num_controls + num_lag_control))
    j_stratum$status <- status
    
    lag <- c(0:30, rep(0, num_controls), rep(c(1:30), each = num_controls))
    j_stratum$lag <- lag[1:nrow(j_stratum)]  # Ensure lag vector matches row count
    
    if(j == 1){
      new_df <- j_stratum
    }else{
      new_df <- rbind(new_df, j_stratum)
    }
  }
  
  if(i == 1){
    matched_df <- new_df
  }else{
    matched_df <- rbind(matched_df, new_df)
  }
}

# matched_df is a matched multi-county data set

matched_df_check <- matched_df %>%
  filter(cand_control == "FALSE" & status == "control")

matched_control <- matched_df %>%
  filter(status == 'control')
table(matched_control$disaster)

write_parquet(matched_df, "Matched_DF_SC_TC_30day_SU.parquet")


