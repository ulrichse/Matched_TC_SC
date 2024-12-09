library(dplyr)
library(data.table)
library(arrow)
library(forestploter)
getwd()
setwd("C:/Users/ulric/OneDrive - University of North Carolina at Chapel Hill/RStudio/Matched-Analysis/TC-Exposure/Hurricane_Maternal_Health/")

dt <- read.csv("daily_lag_all_outcomes.csv")

#forest_data <- pivot_longer(dt, cols = c("lag0", "lag1", "lag2", "lag3", "lag4", "lag5", "lag6", "lag7"), names_to = "Lag")
forest_data <- pivot_wider(dt, names_from = c("type"), values_from = c(starts_with("lag")))
dt <- forest_data

dt$`Lag0 RR` <- paste(rep(" ", 10), collapse = " ")
dt$`Lag1 RR` <- paste(rep(" ", 10), collapse = " ")
dt$`Lag2 RR` <- paste(rep(" ", 10), collapse = " ")
dt$`Lag3 RR` <- paste(rep(" ", 10), collapse = " ")
dt$`Lag4 RR` <- paste(rep(" ", 10), collapse = " ")
dt$`Lag5 RR` <- paste(rep(" ", 10), collapse = " ")
dt$`Lag6 RR` <- paste(rep(" ", 10), collapse = " ")
dt$`Lag7 RR` <- paste(rep(" ", 10), collapse = " ")

dt$'Lag0 95% CI' <- paste(sprintf("%.2f (%.2f-%.2f)", dt$'lag0_matRRfit', dt$'lag0_matRRlow', dt$'lag0_matRRhigh'))
dt$'Lag1 95% CI' <- paste(sprintf("%.2f (%.2f-%.2f)", dt$'lag1_matRRfit', dt$'lag1_matRRlow', dt$'lag1_matRRhigh'))
dt$'Lag2 95% CI' <- paste(sprintf("%.2f (%.2f-%.2f)", dt$'lag2_matRRfit', dt$'lag2_matRRlow', dt$'lag2_matRRhigh'))
dt$'Lag3 95% CI' <- paste(sprintf("%.2f (%.2f-%.2f)", dt$'lag3_matRRfit', dt$'lag3_matRRlow', dt$'lag3_matRRhigh'))
dt$'Lag4 95% CI' <- paste(sprintf("%.2f (%.2f-%.2f)", dt$'lag4_matRRfit', dt$'lag4_matRRlow', dt$'lag4_matRRhigh'))
dt$'Lag5 95% CI' <- paste(sprintf("%.2f (%.2f-%.2f)", dt$'lag5_matRRfit', dt$'lag5_matRRlow', dt$'lag5_matRRhigh'))
dt$'Lag6 95% CI' <- paste(sprintf("%.2f (%.2f-%.2f)", dt$'lag6_matRRfit', dt$'lag6_matRRlow', dt$'lag6_matRRhigh'))
dt$'Lag7 95% CI' <- paste(sprintf("%.2f (%.2f-%.2f)", dt$'lag7_matRRfit', dt$'lag7_matRRlow', dt$'lag7_matRRhigh'))


est = list(dt$'lag0_matRRfit',
           dt$'lag1_matRRfit',
           dt$'lag2_matRRfit',
           dt$'lag3_matRRfit',
           dt$'lag4_matRRfit',
           dt$'lag5_matRRfit',
           dt$'lag6_matRRfit',
           dt$'lag7_matRRfit')
lower =list(dt$'lag0_matRRlow',
            dt$'lag1_matRRlow',
            dt$'lag2_matRRlow',
            dt$'lag3_matRRlow',
            dt$'lag4_matRRlow',
            dt$'lag5_matRRlow',
            dt$'lag6_matRRlow',
            dt$'lag7_matRRlow')
upper = list(dt$'lag0_matRRhigh',
             dt$'lag1_matRRhigh',
             dt$'lag2_matRRhigh',
             dt$'lag3_matRRhigh',
             dt$'lag4_matRRhigh',
             dt$'lag5_matRRhigh',
             dt$'lag6_matRRhigh',
             dt$'lag7_matRRhigh')

dt <- dt %>%
  dplyr::select('outcome',
                'Lag0 RR',
                'Lag0 95% CI',
                'Lag1 RR',
                'Lag1 95% CI',
                'Lag2 RR',
                'Lag2 95% CI',
                'Lag3 RR',
                'Lag3 95% CI',
                'Lag4 RR',
                'Lag4 95% CI',
                'Lag5 RR',
                'Lag5 95% CI',
                'Lag6 RR',
                'Lag6 95% CI',
                'Lag7 RR',
                'Lag7 95% CI',
                )

tm <- forest_theme(base_size = 8,
                   core=list(bg_params=list(fill = "white")))
p <- forest(dt,
            est = est,
            lower = lower, 
            upper = upper,
            ci_column = c(2, 4, 6, 8, 10, 12, 14, 16),
            xlim = c(0.5, 2),
            ticks_at = c(0.5, 1, 2),
            ref_line = 1,
            x_trans = "log",
            theme = tm)


plot(p) 

# lag0 lag7

dt <- read.csv("lag0_lag7_cumulative_all_outcomes_30d.csv")

dt <- dt %>%
  dplyr::select(-X, -estmean)%>%
  dplyr::distinct(over_rr, over_rr_se, over_rr_high, over_rr_low, outcome)%>%
  dplyr::rename(Outcome=outcome)

dt$`Lag0-Lag7 RR` <- paste(rep(" ", 20), collapse = " ")

dt$'95% CI' <- paste(sprintf("%.2f (%.2f-%.2f)", dt$'over_rr', dt$'over_rr_low', dt$'over_rr_high'))

est = list(dt$'over_rr')

lower =list(dt$'over_rr_low')

upper = list(dt$'over_rr_high')

dt <- dt %>%
  dplyr::select(Outcome,
                'Lag0-Lag7 RR',
                '95% CI')
tm <- forest_theme(base_size = 8,
                   core=list(bg_params=list(fill = "white")))
p <- forest(dt,
            est = est,
            lower = lower, 
            upper = upper,
            ci_column = c(2),
            ref_line = 1,
            x_trans = "log",
            theme = tm
)

plot(p)
