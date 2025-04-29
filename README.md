# Code for R 4.4.2
# Code for variable (mean temperature)
library(openair)
library(dlnm)
library(splines)
library(tsModel)
library(ggplot2)
bc <- read.csv("C:\\Users\\hudada12\\Documents\\R-Studio\\data.csv")
head (bc)
bc$date <- as.Date(bc$date)
bc$dow <- weekdays(bc$date)
# Time-series plot
png("Temperature.png", width=3000, height=2000, res=300)
timePlot(mydata=bc, pollutant=c("Temperature"), y.relation="free", xlab="Time (year)", ylab=" Mean temperature (°c)", col="#BD6263")
dev.off()
# Bubble chart code
scatter_plot <- function(data, input1, piece) {
  input <- data[[input1]]
  m <- max(input, na.rm = TRUE)
  n <- min(input, na.rm = TRUE)
  output <- data.frame(value = numeric(), case = numeric(), number = integer())
  for (i in ceiling(piece * n / m):(piece - 1)) {
    condition <- (m * (i) / piece < input) & (input < m * (i + 1) / piece)
    y_rate <- mean(data[condition, 'case'], na.rm = TRUE) 
    number <- sum(condition, na.rm = TRUE)
    output <- rbind(output, data.frame(value = (i / piece) * m + m / (2 * piece), case = y_rate, number = number))
  }
  ggplot(output, aes(x = value, y = case, size = number)) +
    geom_point(alpha = 0.7, color = "#08306b") +
    theme(legend.position = "bottom") +
    xlab("Mean temperature (℃)") + 
    ylab("Case / day") +            theme_bw() +
    scale_size_continuous(range = c(0.5, 5))
}
png(file = "Mean temperatuure Bubble chart.png", width = 1500, height =1500, res = 300)
scatter_plot(bc, 'Temperature', 145) 
dev.off() 
# GAM
library(foreign)
library(mgcv)
bc <- read.csv("C:\\Users\\hudada12\\Documents\\R-Studio\\data.csv")
head (bc)
mgam1 <- gam(case ~ s(Temperature),
             data = bc,
             family = poisson(),
             method = "GCV.Cp")
summary(mgam1)
gam.check(mgam1)
disp_ratio1 <- deviance(mgam1)/df.residual(mgam1)
tiff("GAM_-temperatue.tiff", width = 20, height = 16, units = "cm", res = 300,
     compression = "lzw", family = "Times")
par(mfrow = c(1, 1),
    oma = c(1, 1, 1, 1),  
    mar = c(3.5, 3.5, 1.5, 1.5), 
    mgp = c(2.2, 0.6, 0),  
    tcl = -0.4)  
plot_effect <- function(model, term, xlabel)
{
  plot(model, select = term,
       shade = TRUE,
       shade.col = rgb(0.96, 0.5, 0.5, 0.25),  
       col = "#8B0000",  
       lwd = 1.2,
       xlab = xlabel,
       cex.axis = 0.85,
       cex.lab = 0.95,
       rug = FALSE,
       scheme = 1,
       axes = FALSE)
  axis(1, cex.axis = 0.85, col.axis = "gray30") 
  axis(2, cex.axis = 0.85, col.axis = "gray30")
  box(col = "gray60")  
}
plot_effect(mgam, 1, "Temperature (℃)")
dev.off()
# VIF collinearity diagnosis
lm_for_vif <- glm(case ~ Temperature + Atmospheric.pressure + Relative.humidity + Wind.speed + Sunshine.duration,
                  data = bc,
                  family = poisson())
library(car)
vif_values <- vif(lm_for_vif)
print(vif_values)
# Multi-GAM
library(foreign)
library(mgcv)
Mgam1 <- gam(case ~ s(Temperature, k = 15) + s(Relative.humidity) + s(Wind.speed) + s(Sunshine.duration) + s(date.month, bs = "cc", k = 12) + s(date.year, k = 10),
             data = bc,
             family = poisson(),
             knots = list(date.month = c(1, 12)),
             method = "GCV.Cp")
summary(Mgam1)
gam.check(Mgam1)
disp_ratio <- deviance(Mgam1)/df.residual(Mgam1)
tiff("GAM_more.tiff", width = 20, height = 16, units = "cm", res = 300,
     compression = "lzw", family = "Times")
par(mfrow = c(2, 2),
    oma = c(1, 1, 1, 1),  
    mar = c(3.5, 3.5, 1.5, 1.5), 
    mgp = c(2.2, 0.6, 0),  
    tcl = -0.4)  
plot_effect <- function(model, term, xlabel) {
  plot(model, select = term,
       shade = TRUE,
       shade.col = rgb(0.96, 0.5, 0.5, 0.25), 
       col = "#8B0000",  
       lwd = 1.2,
       xlab = xlabel,
       cex.axis = 0.85,
       cex.lab = 0.95,
       rug = FALSE,
       scheme = 1,
       axes = FALSE)
  axis(1, cex.axis = 0.85, col.axis = "gray30")  
  axis(2, cex.axis = 0.85, col.axis = "gray30")  
  box(col = "gray60")  
}
plot_effect(Mgam1, 1, "Temperature (°C)")
plot_effect(Mgam1, 2, "Relative humidity (%)")
plot_effect(Mgam1, 3, "Wind speed (m/s)")
plot_effect(Mgam1, 4, "Sunshine duration (h)")
dev.off()

# Interaction Variable: Temperature and “Relative.humidity”
Mgam2 <- gam(case ~ s(Temperature,Relative.humidity, bs="ad")+ s(date.month, bs = "cc", k = 12) + s(date.year, k = 10),
             data = bc,
             family = poisson(),
             knots = list(date.month = c(1, 12)),
             method = "GCV.Cp")
gam.check(Mgam2)
disp_ratio2 <- deviance(Mgam2)/df.residual(Mgam2)
tiff("GAM-TR1.tiff",
     width = 16,              
     height = 12,
     units = "cm",
     res = 600,               
     compression = "lzw")       
par(
  mar = c(3.2, 3.2, 1.5, 1.5),  
  mgp = c(2, 0.6, 0),          
  cex.axis = 0.75,             
  cex.lab = 0.85,             
  font.lab = 1                )
vis.gam(Mgam2,
        view = c("Temperature", "Relative.humidity"),
        theta = 30,
        ticktype = "detailed",
        xlab = "Temperature (°C)",
        ylab = "Relative humidity (%)",
        axis.args = list(cex.axis = 0.7),   
        legend.args = list(text = "Risk", cex = 0.8)  
)
dev.off()

# DLNM model eg: Variable Temperature
library(dlnm)
library(splines)
library(tsModel)
library(ggplot2)
bc <- read.csv("C:\\Users\\hudada12\\Documents\\R-Studio\\data.csv")
head (bc)
bc$date <- as.Date(bc$date)
bc$dow <- weekdays(bc$date)
cb1.temp <- crossbasis (bc$Temperature, lag=14, argvar=list(fun="ns",df=4), arglag=list(fun="ns", df=3))
model <- glm(case ~ cb1.temp + ns(time, 8*7) + dow, family=quasipoisson(), data=bc)
pred1.temp = crosspred(cb1.temp, model, cen=round(median(bc$Temperature)), bylag=0.2)
#Figure 5
png(file = "single_temp_3D.png", width = 3000, height = 2500, res = 300)
plot(pred1.temp,ticktype='detailed',border='#3366FF',xlab=" Mean temperature (°C)",ylab="Lag (days)",zlab="RR",col='#99FFCC',shade = 0.1,cex.lab=1.3,cex.axis=1.3,lwd=1,theta = 20, phi = 25,ltheta = -35)
dev.off()
#Figure 6
png(file = "Temperature 0-14days.png", width = 1600, height = 1600, res = 300)
par(mfrow=c(1,1))
crall <- crossreduce(cb1.temp,model,cen= round(median(bc$Temperature)),type="overall",lag=c(0,14))
plot(crall,xlab="Mean temperature (°C)",ylab="RR", col=2,lwd=2,cex.lab=1.2,cex.axis=1.2,mar=c(1,2,0,1))
mtext(text="Overall cumulative association 0-14 days",cex=0.89)
dev.off()
#Figure 7
percentile_5 <- quantile(bc$Temperature, 0.05, na.rm = TRUE)
percentile_95 <- quantile(bc$Temperature, 0.95, na.rm = TRUE)
pred <- crosspred(cb1.temp, model, bylag = 1,cen= round(median(bc$Temperature)),cumul=TRUE)
png(file = "lag-P95.png", width = 3000, height = 2500, res = 300)
plot(pred, var= percentile_95, cumul=TRUE, ci="bars", type="p", xlab=" Lag (days)",
     ylab="RR", col="#FF0000", main= "Mean temperature-P95", lwd=2,cex.lab=1.2,cex.axis=1.2,mar=c(1,2,0,1), pch=20, bg="#FF0000")
dev.off()




