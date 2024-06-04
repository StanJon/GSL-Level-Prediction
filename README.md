# GSL_hydraulic_model

The hydraulic model model is intended for Great Slave Lake water level prediction. Historic water level and discharge observations are used to estimate Residual Net Total Supply (NTS) time-series. Residual NTS is applied to lake storage. Monthly stage-discharge relationships are fitted to historic observations and used to predict the discharge. Pool Level Routing method is used to calculate water level on a daily time step. The model is divided into the following sections:

# Input Section
Contains definition of model parameters. The following parameters can be defined by the user:

Duration - Sets the duration of simulation

Prediction_lag - Lag of simulation start date in days. Setting at 0 starts simulation from current date. Can be set up to 27, which will start simulation 27 days before current date and overlay real-time observations on the final plot.

Prediction_percentile - Sets the prediction percentile for final plot. Value range between 0 and 1.

upper_limit & lower_limit - Set the percentiles of prediction envelope. Value range between 0 and 1.

# NTS Calculation Section
In this section, the date range of historic data is set. Discharge and water level historic observations are downloaded from Hydat. Daily water level change is calculated. Conversion to flow is performed and outlet is added to calculate Residual Net Total Supply (NTS). Daily NTS values are arranged in a matrix with each row representing a year from the historic time-series and columns containing NTS values for the relevant day of prediction period.

# Raltime Data Section
Realtime water level for the past 30 days is downloaded and daily average values are calculated. Linear trendline is fitted and statistics are calculated and stored in summary_current data frame.

# Fitting Monthly Polynomial Models Section
Monthly stage-discharge relationships are established by fitting second-order polynomial models for each month. Coefficient of Determination statistics are calculated for each model. An optional plot showing the fitted models can be enabled by printing consolidated_plots at the end of this section. 

# Level Prediction Section
Matrices for predicted outlet and water level are created. Intermediate water level including NTS is calculated and discharge is predicted using stage-discharge relationship. Pool Level Routing method is applied and results for the entire ensemble are stored in result_matrix.

# Plotting Section
Two plots are produced. The first plot contains the entire ensemble of predicted water levels. The second plot shows the prediction_percentile time series defined in input section and an envelope of uncertainty defined by upper_limit & lower_limit.

