# Tue Jun 15 10:29:42 2021 ------------------------------


#'@T.Alonaizi
#'Header# Bacillus subtilis 11715

#'Meta data
#'Raw data file= 'b.sub.may21'
#'OD: @550nm
#'Time: 24 hour sampling (1 hour intervals up to 24h hours)
#'Rep: three replicates
#'CFU/ml: colony count was adjusted to cfu/ml in raw data
#'
#':::::::::::PACKAGES USED:::::::::::'#

#############LOAD PACKAGES###############
library(Hmisc)
library(patchwork)
library(ggplot2)
library(growthcurver)



#Import data#
b.sub <- read.table('Raw data/b.sub.may21.txt', header = TRUE, 
                    stringsAsFactors = TRUE,na.strings = " NA") #Create dataframe obj



str(b.sub) #display the structure of your dataframe
#'----results below---'#
# 'data.frame':	75 obs. of  4 variables:
#   $ rep : int  1 2 3 1 2 3 1 2 3 1 ...
# $ time: int  0 0 0 1 1 1 2 2 2 3 ...
# $ CFU : num  293333 263333 246667 350000 273333 ...
# $ OD  : num  0.012 0.012 0.013 0.012 0.014 0.014 0.025 0.017 0.02 0.025 ...'
summary(b.sub) #display the summary of your dataframe
# rep         time         CFU                  OD       
# Min.   :1   Min.   : 0   Min.   :   246667   Min.   :0.012  
# 1st Qu.:1   1st Qu.: 6   1st Qu.: 16666667   1st Qu.:0.258  
# Median :2   Median :12   Median :290000000   Median :2.620  
# Mean   :2   Mean   :12   Mean   :219177556   Mean   :2.441  
# 3rd Qu.:3   3rd Qu.:18   3rd Qu.:356666667   3rd Qu.:4.380  
# Max.   :3   Max.   :24   Max.   :453333333   Max.   :4.850 

#####Data exploration#####
#convert time and Rep into factors
b.sub$Ftime <- as.factor(b.sub$time)
b.sub$FRep <- as.factor(b.sub$rep)

##plot the OD values##

obj1<- ggplot(b.sub, aes(x=Ftime, y=OD, group=1)) +  #Normal values
  stat_summary(
    fun = mean,
    geom='line') +
  stat_summary(
    fun=mean,
    geom='point') +
  stat_summary(
    fun.data=mean_cl_boot,
    geom='errorbar',
    width=0.5) +
  labs(x = "Time (hours)", y = "OD @550nm", 
       title = "B.subtilis 11715 OD growth assay") +
  theme_bw()

obj2 <- ggplot(b.sub, aes(x=Ftime, y=log(OD), group=1)) +  #Natural log
  stat_summary(
    fun = mean,
    geom='line') +
  stat_summary(
    fun=mean,
    geom='point') +
  stat_summary(
    fun.data=mean_cl_boot,
    geom='errorbar',
    width=0.5) +
  labs(x = "Time (hours)", y = "OD @550nm", 
       title = "[Natural log] B.subtilis 11715 OD growth assay") +
  theme_bw()

##Cutommise graphs graphs##

n.obj1<- obj1+theme(
  plot.title = element_text(color="red", size=8, face="bold.italic"),
  axis.title.x = element_text(color="blue", size=9, face="bold"),
  axis.title.y = element_text(color="#993333", size=9, face="bold")
)
n.obj2  <- obj2+theme(
  plot.title = element_text(color="red", size=8, face="bold.italic"),
  axis.title.x = element_text(color="blue", size=9, face="bold"),
  axis.title.y = element_text(color="#993333", size=9, face="bold")
)

###Print graphs###
n.obj1+n.obj2

##Print CFU curve##
obj3<- ggplot(b.sub, aes(x=Ftime, y=CFU, group=1)) +  #Normal values
  stat_summary(
    fun = mean,
    geom='line') +
  stat_summary(
    fun=mean,
    geom='point') +
  stat_summary(
    fun.data=mean_se,
    geom='errorbar',
    width=0.5) +
  labs(x = "Time (hours)", y = "CFU/ml", 
       title = "B.subtilis 11715 growth assay") +
  theme_bw()

obj4 <- ggplot(b.sub, aes(x=Ftime, y=log(CFU), group=1)) +  #Normal values
  stat_summary(
    fun = mean,
    geom='line') +
  stat_summary(
    fun=mean,
    geom='point') +
  stat_summary(
    fun.data=mean_se,
    geom='errorbar',
    width=0.5) +
  labs(x = "Time (hours)", y = "LOG CFU/ml", 
       title = "[LOG] B.subtilis 11715 growth assay") +
  theme_bw()

#customise CFU graphs
n.obj3<- obj3+theme(
  plot.title = element_text(color="red", size=8, face="bold.italic"),
  axis.title.x = element_text(color="blue", size=9, face="bold"),
  axis.title.y = element_text(color="#993333", size=9, face="bold")
)
n.obj4  <- obj4+theme(
  plot.title = element_text(color="red", size=8, face="bold.italic"),
  axis.title.x = element_text(color="blue", size=9, face="bold"),
  axis.title.y = element_text(color="#993333", size=9, face="bold")
)

#print CFU graphs
n.obj3+n.obj4
obj3+obj4


#Plot OD values based on Growthcurver package
gc_bsub <- SummarizeGrowth(b.sub$time, b.sub$OD)
gc_bsub
plot(gc_bsub)
str(gc_bsub$vals)


###Standard line ODxCFU###
#import standard data
OD.CFU_std <- read.table('Raw data/bsubmeanstd.txt', header = TRUE, 
                          stringsAsFactors = TRUE )
#
ODxCFU <- lm(data = OD.CFU_std,OD_mean~ CFU_mean)
par(mfrow=c(2,2))
plot(ODxCFU)

#75% of the residuals are explained from the linear model
#Taking the mean value will increase the value
plot(ODxCFU, which = 4)
##The results show that the residuals are normal in the
#Q-Q plot however there are some residuals that have influnce of the lm
# residuals number 35,36 and 38

#Checking and Plotting the mean values
odcfu.lm <- lm(data = OD.CFU_std,CFU_mean~ OD_mean)
par(mfrow=c(2,2))
plot(odcfu.lm)
#mean values look more presentable
summary(odcfu.lm)
#98% of the residuals fit the model
anova(odcfu.lm)

OD.CFU_std$ODCFULM <- predict(odcfu.lm)

my.data <- data.frame(OD_mean = seq(from = min(OD.CFU_std$OD_mean),
                               to = max(OD.CFU_std$OD_mean),
                               length = 28))
pred.vals <- predict(odcfu.lm, newdata = my.data)

par(mfrow=c(1,1))
plot(OD.CFU_std$OD_mean, OD.CFU_std$CFU_mean, xlab = "OD", ylab = "CFU", main = 'CFU/ml|OD standard line')

lines(my.data$OD, pred.vals, lty = 1,col = 2)

pred.vals.se <- predict(odcfu.lm, newdata = my.data, se.fit = TRUE)
lines(my.data$OD, pred.vals.se$fit, lty = 1,col = 1)


# add the upper 95% confidence interval
lines(my.data$OD, pred.vals.se$fit + (1.96 * pred.vals.se$se.fit), lty = 2, col = 2)

# add the lower 95% confidence interval
lines(my.data$OD, pred.vals.se$fit - (1.96 * pred.vals.se$se.fit), lty = 2, col = 2)