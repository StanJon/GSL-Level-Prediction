library(tidyhydat)
library(tidyverse)
library(ggplot2)
library(zoo)
library(patchwork)

#*****************************INPUT SECTION*************************************
#Declare start date
Start_date <- Sys.Date()

#Declare duration of prediction
duration <- 182

#Declare value to lag prediction from sys.date (max 27)
prediction_lag <- 0

#Declare prediction quantile and limits of plot ribbon
prediction_percentile <- 0.10
upper_limit <- 0.15
lower_limit <- 0.05

#Declare lake area variable in m^2
Larea <- 28568 * 1000 * 1000

#Declare conversion factor from cms to metres of lake level per day
cms <- (60 * 60 * 24) / (28568 * 1000 * 1000)

#*****************************NTS CALCULATION***********************************

#Set the sequence of dates for the Duration
dates <- seq(Start_date-prediction_lag, by = "day", length.out = duration + prediction_lag)

#Exclude the additional day if the year is leap year
dates <- dates[format(dates, "%m-%d") != "02-29"]

#Define date range for historic NTS calculation
start_d <- "1995-10-01"
end_d <- "2021-09-30"

#Calculate number of years for matrix sizing
num_years <- as.numeric(substr(end_d, 1, 4)) - as.numeric(substr(start_d, 1, 4))

#Pull historic lake level and lake outlet data for defined time period
GSL_levels <- subset(hy_daily_levels("07SB001", start_date = start_d, end_date = end_d), select = c(Date, Value))
GSL_outlet <- subset(hy_daily_flows("10FB006",start_date=start_d,end_date=end_d),select=c(Date,Value))

#Check the number of missing values in each time series
num_na_levels <- sum(is.na(GSL_levels))
num_na_outlet <- sum(is.na(GSL_outlet))

#Merge historic data into a single data frame
merged <-merge(GSL_levels,GSL_outlet, by="Date", all=TRUE)
merged <- merged[format(merged$Date, "%m-%d") != "02-29", ]
colnames(merged) <- c("Date","Level","Outlet")
num_na_merged <- sum(is.na(merged))


#Lag level data and calculate daily change in water level
merged$DeltaLevel <- (merged$Level-lag(merged$Level,n=1))
merged$DeltaLevel[1] <- merged$DeltaLevel[2]

#Convert daily change in water level to flow in CMS 
merged$DeltaFlow <- merged$DeltaLevel * Larea / (60 * 60 * 24)

#Calculate daily NTS
merged$NTS <- merged$DeltaFlow + merged$Outlet

#Add year and month-day columns
merged$Year <- format(merged$Date, "%Y")
merged$Month_day <- format(merged$Date, "%m-%d")

#Arrange historic time series data for the duration
#Create a matrix for historic water levels
matched_level <- data.frame(matrix(ncol = length(dates), nrow = num_years))
colnames(matched_level) <- paste0("Day_", 1:length(dates))

#Create a matrix for historic daily NTS
matched_NTS <- data.frame(matrix(ncol = length(dates), nrow = num_years))
colnames(matched_NTS) <- paste0("Day_", 1:length(dates))

#Arrange NTS and level data in the two matrices
#Each column in matrix is a day in a time sequence of prediction duration
#Each row of matrix is a separate year
for (i in 1:length(dates)) {
  current_month_day <- format(dates[i], "%m-%d")
  date_index <- which(merged$Month_day == current_month_day)
  
  if (length(date_index) > 0) {
    matched_level[, i] <- merged[date_index, "Level"]
    matched_NTS[, i] <- merged[date_index, "NTS"]
  } else {
    if (i > 1) {
      matched_level[, i] <- matched_level[, i - 1]
      matched_NTS[, i] <- matched_NTS[, i - 1]
    } else {
      matched_level[, i] <- NA
      matched_NTS[, i] <- NA}}}

#*****************************REALTIME DATA*************************************
#Pull realtime lake water level data
#Downloads hourly measurement values for the last 30 days
realtime <- na.omit(subset(realtime_dd("07SB001"), select = c(Date, Value)))
realtime$Date <- as.POSIXct(realtime$Date)

#Calculate daily average water levels
daily_avg <- realtime %>%
  group_by(Date = as.Date(Date)) %>%
  summarize(DailyAverage = mean(Value))

#Separate daily average values for the last 28 days
last_28_days <- daily_avg[daily_avg$Date >= max(daily_avg$Date) - 27, ]

#Fit linear model to determine the trend
linear_model <- lm(DailyAverage ~ Date, data = last_28_days)

#Generate statistics of data over the last 28 days
#Report min, max, mean, standard deviation and linear trendline slope
summary_current <- data.frame(
  min <- min(last_28_days$DailyAverage),
  max <- max(last_28_days$DailyAverage),
  mean <- mean(last_28_days$DailyAverage),
  sd <- sd(last_28_days$DailyAverage),
  slope <- coef(linear_model)[2])
colnames(summary_current) <- c("min", "max", "mean", "standard_deviation", "Slope")

#*************************FITTING MONTHLY POLYNOMIAL MODELS**********************

#Create a list to store equations for each month
month_eq_list <- list()

#Fit equations for each month
for (month_idx in 1:12) {
  #Filter data for the current month
  month_data <- merged[month(merged$Date) == month_idx, ]
  
  #Fit polynomial equation to the current month
  month_eq <- lm(Outlet ~ poly(Level, 2), data = month_data)
  
  #Store equation for the current month
  month_eq_list[[month_idx]] <- month_eq
}

#**Evaluation of models**
#Initialize vectors to store R-squared values for each month
month_r_squared <- numeric(12)

#Calculate R-squared values for each month
for (month_idx in 1:12) {
  #Filter data for the current month
  month_data <- merged[month(merged$Date) == month_idx, ]
  
  #Predict values for the current month
  predicted_values <- predict(month_eq_list[[month_idx]], newdata = month_data)
  
  #Calculate R-squared for the current month
  month_r_squared[month_idx] <- summary(month_eq_list[[month_idx]])$r.squared
}

#Function to calculate RMSE
rmse <- function(actual, predicted) {
  sqrt(mean((predicted - actual)^2))
}

#Function to calculate R-squared
calculate_r_squared <- function(observed, predicted) {
  obs_mean <- mean(observed)
  ss_total <- sum((observed - obs_mean)^2)
  ss_residual <- sum((observed - predicted)^2)
  return(1 - ss_residual / ss_total)
}

#Predictions from polynomial models
predicted_values_list <- list()
for (month_idx in 1:12) {
  predicted_values_list[[month_idx]] <- predict(month_eq_list[[month_idx]])
}

#Calculate R-squared, RMSE, and NSE for each month
r_squared_list <- numeric(12)
rmse_list <- numeric(12)
nse_list <- numeric(12)
for (month_idx in 1:12) {
  r_squared_list[month_idx] <- calculate_r_squared(merged[month(merged$Date) == month_idx, ]$Outlet, predicted_values_list[[month_idx]])
  rmse_list[month_idx] <- rmse(merged[month(merged$Date) == month_idx, ]$Outlet, predicted_values_list[[month_idx]])
}

#Print the calculated values for each month
for (month_idx in 1:12) {
  cat(month.name[month_idx], "model: R-squared:", r_squared_list[month_idx], "RMSE:", rmse_list[month_idx], "\n\n")
}

# Check the range of level values
level_min <- min(merged$Level)
level_max <- max(merged$Level)

#set the x-axis range
x_range <- range(156:157.5)
x_ticks <- seq(156, 157.5, by = 0.5)

#Check the range of outlet values
outlet_min <- min(merged$Outlet)
outlet_max <- max(merged$Outlet)

#set the y-axis range
y_range <- range(1000:13000)
y_ticks <- seq(1000, 13000, by = 4000)

#Plot data for each month with fitted models
plot_list <- list()
for (month_idx in 1:12) {
  r_squared_label <- substitute(R^2 == r_squared_value, list(r_squared_value = round(r_squared_list[month_idx], 2)))

  plot_month <- ggplot(merged[month(merged$Date) == month_idx, ], aes(x = Level, y = Outlet)) +
    geom_point(size = 0.5) +
    geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = FALSE, color = "blue") +
    labs(title = paste(month.name[month_idx], "Model"), x = "Level (m)", y = expression("Flow" ~ (m^3/s))) +
    theme_minimal() +
    scale_x_continuous(limits = x_range, breaks = x_ticks) +
    scale_y_continuous(limits = y_range, breaks = y_ticks) +
    annotate("text", x = Inf, y = -Inf, label = r_squared_label, hjust = 2, vjust = -7, size = 3, color = "black")  # Add R-squared annotation
  plot_list[[month_idx]] <- plot_month
}


# Arrange all monthly plots into one image
consolidated_plot <- wrap_plots(
  plot_list[[1]], plot_list[[2]], plot_list[[3]],
  plot_list[[4]], plot_list[[5]], plot_list[[6]],
  plot_list[[7]], plot_list[[8]], plot_list[[9]],
  plot_list[[10]], plot_list[[11]], plot_list[[12]],
  nrow = 4
)

# Print the consolidated plot
# print(consolidated_plot)

#*****************************LEVEL PREDICTION**********************************
#Create a matrix for predicted level values
result_matrix <- matrix(NA, nrow = nrow(matched_NTS), ncol = ncol(matched_NTS))
m_loc <- 28 - prediction_lag
result_matrix[] <- last_28_days$DailyAverage[m_loc]

#Create a matrix to store assigned daily NTS, outlet and intermediate level values (DIAGNOSTICS)
nts_matrix <- matrix(NA, nrow = nrow(matched_NTS), ncol = ncol(matched_NTS))
intermediate_matrix <- matrix(NA, nrow = nrow(matched_NTS), ncol = ncol(matched_NTS))
outlet_matrix <- matrix(NA, nrow = nrow(matched_NTS), ncol = ncol(matched_NTS))

#Loop over each row in matched_NTS
for (i in 1:nrow(matched_NTS)) {
  #Initialize current level with the initial value
  current_level <- last_28_days$DailyAverage[m_loc]
  
  #Loop over each column (day) in matched_NTS
  #Calculate water level using the level pool routing method:
  #(Inlet1+Inlet2) / 2 * dt - (Outlet1+Outlet2) / 2 * dt = dS
    for (j in 2:ncol(matched_NTS)) { 
    current_month <- month(dates[j])
    current_nts <- matched_NTS[i, j]
    
    #Get equation for the current month
    current_month_eq <- month_eq_list[[current_month]]
    
    #Calculate intermediate level
    intermediate_level <- current_level + matched_NTS[i, j] * cms
    
    #Predict outlet using intermediate level
    predicted_outlet <- predict(current_month_eq, newdata = data.frame(Level = intermediate_level))
    
    #Calculate delta level
    delta_level <- (((matched_NTS[i, j] + matched_NTS[i, j-1]) / 2) - ((predicted_outlet + predict(current_month_eq, newdata = data.frame(Level = current_level)))) / 2) * cms
    
    #Update current level using the calculated result
    current_level <- current_level + delta_level
    
    #Ensure the matrix index is within bounds before assignment
    if (j <= ncol(result_matrix)) {
      result_matrix[i, j] <- current_level
      nts_matrix[i, j] <- current_nts
      intermediate_matrix[i, j] <- intermediate_level
      outlet_matrix[i, j] <- predicted_outlet
        
      #Print relevant information
      cat("Scenario:", i, ", Day:", j, ", Equation:", current_month, ", Water level:", current_level, ", NTS:", current_nts,  "\n")
    } else {
      break  #Stop the loop if index exceeds matrix bounds
    }
  }
}

#Move results to a data frame
df_results <- as.data.frame(result_matrix)

#Mutate prediction results for plotting
df_long <- df_results %>%
  mutate(row = row_number()) %>%
  pivot_longer(-row, names_to = "Interval", values_to = "Value") %>%
  mutate(Interval = gsub("V", "", Interval),
         Interval = as.numeric(Interval))

#*****************************PLOTTING******************************************

#Create plot showing all possible time series
plot_TS <- ggplot(df_long, aes(x = (Start_date + Interval - prediction_lag - 1), y = Value, color = as.factor(row))) +
  geom_line() +
  geom_point(data = last_28_days[nrow(last_28_days), ], aes(x = Date, y = DailyAverage), color = "black", size = 1.5) +
  labs(title = "Predicted water level change", x = "Month", y = "Value") +
  geom_line(data = last_28_days, aes(x = Date, y = DailyAverage), color = "blue") +
  scale_x_date(date_labels = "%b", date_breaks = "1 month") +
  theme_minimal()+
  theme(legend.position = "none")

#Plot all possible timeseries
print(plot_TS)

# Calculate the 10th and 90th percentiles based on the median
perc_low_vals <- apply(result_matrix, 2, quantile, prob = lower_limit, na.rm = TRUE)
predict_vals <- apply(result_matrix, 2, quantile, prob = prediction_percentile, na.rm = TRUE)
mean_vals <- apply(result_matrix, 2, mean, na.rm = TRUE)
perc_high_vals <- apply(result_matrix, 2, quantile, prob = upper_limit, na.rm = TRUE)

# Create a data frame for plotting
plot_data <- data.frame(
  Interval = seq_along(predict_vals),
  perc_l = perc_low_vals,
  predict = predict_vals,
  mean = mean_vals,
  perc_h = perc_high_vals
)

#Calculate min and max values for each row in matched_level_sub
#To find historic low and high water levels to be used in plotting
min_vals <- apply(matched_level, 2, min, na.rm = TRUE)
max_vals <- apply(matched_level, 2, max, na.rm = TRUE)

#Add new columns historic_low and historic_high to plot_data
plot_data$historic_low <- min_vals
plot_data$historic_high <- max_vals

#Plotting using ggplot
plot_Q <- ggplot(plot_data, aes(x = (Start_date + Interval - prediction_lag - 1))) +
  geom_line(aes(y = predict, color = "Predicted"), show.legend = TRUE) +
  #geom_line(aes(y = mean, color = "Mean"), show.legend = TRUE)
  geom_ribbon(aes(ymin = perc_l, ymax = perc_h, fill = "Prediction Interval"), alpha = 0.5, show.legend = TRUE) +
  geom_point(data = last_28_days[nrow(last_28_days), ], aes(x = Date, y = DailyAverage), color = "black", size = 1.0) +
  labs(title = "Predicted water level change", x = "Month", y = "Level (m)") +
  geom_line(data = last_28_days, aes(x = Date, y = DailyAverage, color = "Observed"), show.legend = TRUE) +
  geom_ribbon(aes(ymin = historic_low, ymax = historic_high, fill = "Historic Range"), alpha = 0.2, show.legend = TRUE)+
  scale_x_date(date_labels = "%b", date_breaks = "1 month") +
  scale_fill_manual(name = "Interval Type", values = c("Prediction Interval" = "skyblue", "Historic Range" = "grey")) +
  scale_color_manual(name = "Data", values = c("Observed" = "blue", "Predicted" = "red", "Mean" = "green")) +
  theme_minimal()


#Plot all possible timeseries
print(plot_TS)

#Plot the selected quantiles
print(plot_Q)
