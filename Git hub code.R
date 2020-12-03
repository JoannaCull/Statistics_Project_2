#Author: Joanna Cull
#Version: Edit 3
#R version: Version 1.3.1093
#Contact: joanna.cull@dtc.ox.ac.uk


####LOOKING AT DATA SET ####
install.packages('lmerTest')
install.packages('emmeans')
install.packages('languageR')
install.packages('multcomp')
library(tidyverse)
library(cowplot)
library(reshape2)
library(lmerTest)
library(lme4)
library(emmeans)
require(car)
library(languageR)
library(multcomp)
library(ggpubr)


#Getting Data into R #
mice_info<-read.table('Blood_Pressure_Data_edit.txt', header = TRUE)
mice<-as_tibble(mice_info)
View(mice)

#Get basal data together: 
mice_bp_basal_info<-mice %>% dplyr::select(Mouse_Number, Basal_1_SystolicBloodPressure, Basal_2_SystolicBloodPressure, Basal_3_SystolicBloodPressure, Basal_4_SystolicBloodPressure)
View(mice_bp_basal_info)
mice_hr_basal_info<-mice %>% dplyr:: select(Mouse_Number, Basal_1_Pulse, Basal_2_Pulse, Basal_3_Pulse, Basal_4_Pulse)
View(mice_hr_basal_info)

#Take means for blood pressure and heart rate basal information: 
basal_bp_mean<- mice_bp_basal_info %>%dplyr::transmute(Mouse_Number,
                                                       Basal_Systolic_BP_Mean = rowMeans(dplyr::select(., 2:5)))

#View(basal_bp_mean)

basal_hr_mean<-mice_hr_basal_info %>% dplyr::transmute(Mouse_Number,
                                                       Basal_HR_Mean = rowMeans(dplyr::select(. , 2:5)))
#View(basal_hr_mean)

#Add means of basal bp and hr into main mice table as new column:
mice<- cbind(mice, basal_bp_mean[,2])
mice<-cbind(mice,basal_hr_mean[,2])
View(mice)


#### QUESTION 1 ####
#Does genetic background influence systolic blood pressure and heart rate:
#KO and WT are not related - independent get into their subset groups
wildtype<-mice[mice$Genotype=='WT', ]
View(wildtype)
knockout<-mice[mice$Genotype=='KO', ]
### 1) Genetics vs Basal Systolic BP ###
t.test(wildtype$Basal_Systolic_BP_Mean, knockout$Basal_Systolic_BP_Mean) 
#P value >0.05 so not significant (0.5522)
### 2) Genetics vs Basal HR ###
t.test(wildtype$Basal_HR_Mean, knockout$Basal_HR_Mean)
#P value >0.05 so not significant (0.5072)

###Box plot of Basal Systolic BP in WT vs KO Visualization ###
boxplotbp<-ggplot(mice, aes(x= Genotype, y= Basal_Systolic_BP_Mean, color=Genotype))+
  stat_boxplot(geom='errorbar', width=.5)+
  geom_boxplot()+
  geom_point(pch=18, color='black', size=2.5)+
  theme_bw()+
  labs(y= 'Mean Basal Systolic Blood Pressure (mmHg)', 
       title= 'Mean Basal Systolic Blood Pressure in Wild Type and Knockout')
boxplothr<-ggplot(mice, aes(x= Genotype, y= Basal_HR_Mean, color=Genotype))+
  stat_boxplot(geom='errorbar', width=.5)+
  geom_boxplot()+
  geom_point(pch=18, color='black', size=2.5)+
  theme_bw()+
  labs(y= 'Mean Basal Pulse (bpm)', 
       title= 'Mean Basal Pulse in Wild Type and Knockout')

plot_grid(boxplotbp, boxplothr, labels='AUTO')


#### QUESTION 2 ####
#Does Genetics influence BP or HR in pregnancy?
#Making and separating blood pressure data into long form:
mice_bp_data_wide <-mice %>% dplyr::select(Mouse_Number, Genotype, Basal_Systolic_BP_Mean, E4.5_SystolicBloodPressure, E6.5_SystolicBloodPressure, 
                                           E8.5_SystolicBloodPressure, E10.5_SystolicBloodPressure,E12.5_SystolicBloodPressure, 
                                           E14.5_SystolicBloodPressure, E16.5_SystolicBloodPressure, E18.5_SystolicBloodPressure) 
mice_bp_data<- melt(mice_bp_data_wide, id.vars=c("Mouse_Number","Genotype"),
                    variable.name = 'Time_Point', 
                    value.name='Systolic_Blood_Pressure')
#Making and separating HR data into long form:
mice_hr_data_wide<-mice %>% dplyr::select(Mouse_Number,Genotype, Basal_HR_Mean, E4.5_Pulse, E6.5_Pulse, E8.5_Pulse, E10.5_Pulse, E12.5_Pulse,
                                          E14.5_Pulse, E16.5_Pulse, E18.5_Pulse)
mice_hr_data<- melt(mice_hr_data_wide, id.vars=c("Mouse_Number","Genotype"),
                    variable.name = 'Time_Point',  
                    value.name = 'Heart_Rate')
#THEN MADE THE TIME POINTS INTO FACTORS - got time points instead of titles:
mice_bp_data$Time_Point<-factor(mice_bp_data$Time_Point, 
                                levels= c('Basal_Systolic_BP_Mean','E4.5_SystolicBloodPressure','E6.5_SystolicBloodPressure',
                                          'E8.5_SystolicBloodPressure','E10.5_SystolicBloodPressure', 'E12.5_SystolicBloodPressure',
                                          'E14.5_SystolicBloodPressure','E16.5_SystolicBloodPressure','E18.5_SystolicBloodPressure'), 
                                labels= c(0,4.5,6.5,8.5,10.5,12.5,14.5,16.5,18.5))
View(mice_bp_data)

mice_hr_data$Time_Point<-factor(mice_hr_data$Time_Point, 
                                levels= c('Basal_HR_Mean','E4.5_Pulse','E6.5_Pulse',
                                          'E8.5_Pulse','E10.5_Pulse', 'E12.5_Pulse',
                                          'E14.5_Pulse','E16.5_Pulse','E18.5_Pulse'), 
                                labels= c(0,4.5,6.5,8.5,10.5,12.5,14.5,16.5,18.5))
View(mice_hr_data)


overal_hr_bp<- cbind(mice_bp_data,mice_hr_data[,4])
names(overal_hr_bp)<-c('Mouse_Number', 'Genotype', 'Time_Point', 'Systolic_Blood_Pressure', 'Heart_Rate')
View(overal_hr_bp)

####JUST FOR DATA 1) HEART RATE AND 2) BLOOD PRESSURE VISUALISATION FOR MYSELF ####
require(ggpubr)
d1<-ggplot(mice_bp_data %>%filter(Genotype =='WT'), aes(x=Systolic_Blood_Pressure))+
  ggtitle("Wild Type Systolic Blood Pressure")+
  xlab("Systolic Blood Pressure (mmHg)")+
  ylab("Frequency")+
  theme_bw()+
  geom_histogram(colour='black', fill='grey')
d2<-ggplot(mice_bp_data %>%filter(Genotype == 'KO'), aes(x=Systolic_Blood_Pressure))+
  ggtitle("Knock Out Systolic Blood Pressure")+
  xlab("Systolic Blood Pressure (mmHg)")+
  ylab("Frequency")+
  theme_bw()+
  geom_histogram(colour='black', fill='grey')
d3<-ggplot(mice_hr_data %>%filter(Genotype =='WT'), aes(x=Heart_Rate))+
  ggtitle("Wild Type Heart Rate")+
  xlab("Heart Rate (bpm) ")+
  ylab("Frequency")+
  theme_bw()+
  geom_histogram(colour='black', fill='grey')
d4<-ggplot(mice_hr_data %>%filter(Genotype=='KO'), aes(x=Heart_Rate))+
  ggtitle("Knock Out Heart Rate")+
  xlab("Heart Rate (bpm) ")+
  ylab("Frequency")+
  theme_bw()+
  geom_histogram(colour='black', fill='grey')

ggarrange(d1,d2,d3,d4, ncol=2,nrow=2)

#### INITIAL DATA VISUALS FOR DATA VISUALISATION SLIDE ####
# Scatter plot BP vs HR 
# Histogram of BP 
# Histogram of HR 
a<-ggplot(overal_hr_bp, aes(x=Systolic_Blood_Pressure, y= Heart_Rate, color=Genotype))+
  xlab('Systolic Blood Pressure (mmHg)')+
  ylab('Heart Rate (bpm)')+
  geom_point()
b<-ggplot(mice_bp_data, aes(x=Systolic_Blood_Pressure, color=Genotype))+
  labs(title="Systolic Blood Pressure histogram plot",x="Systolic Blood Pressure (mmHg)", y = "Frequency")+
  geom_histogram()
c<-ggplot(mice_hr_data, aes(x= Heart_Rate, color=Genotype))+
  labs(title="Heart Rate histogram plot",x="Heart Rate (bpm)", y = "Frequency")+
  geom_histogram()

ggarrange(a, 
          ggarrange(b,c, ncol=2,nrow=1, labels=c('B','C')),
          nrow=2,
          labels='A')




#### QUESTION 2 ANCOVA ####
#Check assumptions
#1) Normal distribution using Shapiro-Wilk Test:
shapiro.test(mice_bp_data$Systolic_Blood_Pressure) #P>0.05 - is normal p=0.708
shapiro.test(mice_hr_data$Heart_Rate) #P>0.05 - is normal (only just! p=0.05426)
#2) Homoscedasticity using Levenes Test:
leveneTest(mice_bp_data$Systolic_Blood_Pressure~mice_bp_data$Genotype) #P=0.2924 can go ahead!
leveneTest(mice_hr_data$Heart_Rate~mice_hr_data$Genotype)#P=0.6326 so can go ahead! 

View(mice_bp_data)
#Does genetic background influence systolic blood pressure and heart rate in pregnancy 
#ANCOVA - Systolic Blood Pressure: 
#I looked at effect of both genotype and time point:
q2bp<-aov(mice_bp_data$Systolic_Blood_Pressure~mice_bp_data$Genotype+mice_bp_data$Time_Point)
summary(q2bp)
TukeyHSD(q2bp)#post hoc 
#ANCOVA -Heart Rate
#I looked at effect of both genotype and time point --> significant genotype effect when don't account for basal like aksel
q2hr<-aov(mice_hr_data$Heart_Rate~mice_hr_data$Genotype+mice_hr_data$Time_Point)
summary(q2hr)
TukeyHSD(q2hr)#post hoc 


#### QUESTION 3 ####
#Does pregnancy influence systolic blood pressure and heart rate?
#Pregnant- TRUE/FALSE column added into data frame:
mice_pregnant_data<- merge(mice_bp_data,mice_hr_data)
mice_pregnant_data$Pregnant<- c(ifelse(mice_pregnant_data$Time_Point == 0, 'No', 'Yes'))
View(mice_pregnant_data)
#ANCOVA of pregnancy and blood pressure
#Same results as aksel
q3bp<- aov(Systolic_Blood_Pressure~Genotype+Pregnant, data=mice_pregnant_data)
summary(q3bp)
#POST HOC -TUKEY:
post_hoc_bp_1<-glht(q3bp, linfct=mcp(Genotype='Tukey'))
summary(post_hoc_bp_1)

#ANCOVA of pregnancy and heart rate 
q3hr<- aov(Heart_Rate~Genotype+Pregnant, data=mice_pregnant_data)
summary(q3hr)
#POST HOC - TUKEY:
post_hoc_hr_1<-glht(q3hr, linfct=mcp(Genotype='Tukey'))
summary(post_hoc_hr_1)



####QUESTION 4 GLM ####
#Are there interactions between systolic blood pressure and :  > Pregnant state > Gestation Day  > Genetic Background  
#Assumptions are not met for this as not linear data :( 
#This says we want pregnancy, time point and genotype to explain systolic bp 
linearmodel1<-glm(Systolic_Blood_Pressure ~ Pregnant+Time_Point+Genotype, data= mice_pregnant_data)
summary(linearmodel1)
#plot(linearmodel1)
#this says we want pregnancy, time and genotype to explain heart rate 
linearmodel2<-glm(Heart_Rate~ Pregnant+Time_Point+Genotype, data=mice_pregnant_data)
summary(linearmodel2)









#### POWER TESTING ####
#N CALCULATION FOR 0.8 POWER TEST FOR T TEST WE WANT #
t_test_power<-pwr::pwr.t.test(d=1.18, sig.level = 0.05, 
                              power = 0.8,type= 'two.sample', 
                              alternative = 'two.sided')

t_test_power
plot(t_test_power)

#POWER OF THE DATA WE HAVE:
#n is n per group
t_test_actual_power<-pwr::pwr.t.test(n=6, d=0.1,
                                     sig.level = 0.05,
                                     type='two.sample')
t_test_actual_power

#N CALCULATION FOR 0.8 POWER TEST FOR ANOVA WE WANT #
#k=groups = 2 genotypes
#pwr package can only do 1 way anovas 
#f = 0.1 = small effect size 
anova_power<-pwr::pwr.anova.test(k = 2,
                                 f = 0.1, 
                                 sig.level =0.05, 
                                 power =0.8 )
anova_power
plot(anova_power)


#POWER OF THE DATA WE HAVE:

#Pregnancy state effect bit more but small 
#n is number in each group 
anova_actual_power<-pwr::pwr.anova.test(k=2,
                                        f=0.1,
                                        sig.level=0.05,
                                        n=6)
anova_actual_power
plot(anova_actual_power)


#### Working out effect size using MEANS AND SD  ####
#Means of time 0 and time 18.5 
#Difference between them / SD of total data 

#WILDTYPE - MEANS 
#MEAN OF TIME 0 FOR BP - WT 
power_wildtype<-first_last_data[first_last_data$Genotype=='WT', ]
View(power_wildtype)
time_zero_wt_mean<-mean(power_wildtype[power_wildtype$Time_Point==0, ]$Systolic_Blood_Pressure, na.rm= TRUE)
time_zero_wt_mean
#MEAN OF TIME 18.5 BP -WT
time_18.5_wt_mean<-mean(power_wildtype[power_wildtype$Time_Point==18.5, ]$Systolic_Blood_Pressure, na.rm=TRUE)
time_18.5_wt_mean
#MEAN OF TIME 0 FOR HR - WT
time_zero_wt_mean_hr<-mean(power_wildtype[power_wildtype$Time_Point==0, ]$Heart_Rate, na.rm= TRUE)
time_zero_wt_mean_hr
#MEAN OF TIME 18.5 BP -WT
time_18.5_wt_mean_hr<-mean(power_wildtype[power_wildtype$Time_Point==18.5, ]$Heart_Rate, na.rm=TRUE)
time_18.5_wt_mean_hr

#KNOCKOUT MEANS
#MEAN OF TIME 0 - BP FOR KO 
power_knock<-first_last_data[first_last_data$Genotype=='KO', ]
View(power_knock)
time_zero_ko_mean<-mean(power_knock[power_knock$Time_Point==0, ]$Systolic_Blood_Pressure, na.rm= TRUE)
time_zero_ko_mean
#MEAN OF TIME 18.5 BP 
time_18.5_ko_mean<-mean(power_knock[power_knock$Time_Point==18.5, ]$Systolic_Blood_Pressure, na.rm=TRUE)
time_18.5_ko_mean
#MEAN OF TIME 0 - HR
time_zero_ko_mean_hr<-mean(power_knock[power_knock$Time_Point==0, ]$Heart_Rate, na.rm= TRUE)
time_zero_ko_mean_hr
#MEAN OF TIME 18.5 HR 
time_18.5_ko_mean_hr<-mean(power_knock[power_knock$Time_Point==18.5, ]$Heart_Rate, na.rm=TRUE)
time_18.5_ko_mean_hr

#SD OF TOTAL DATA - BP + HR 
sd(mice_pregnant_data$Systolic_Blood_Pressure, na.rm=TRUE)
sd(mice_pregnant_data$Heart_Rate, na.rm=TRUE)












#### ALTERNATIVE APPROACH FROM HERE ON - LOOKING AT TIME 0 vs TIME 18.5 ####




####BASAL TIME POINT AND LAST TIME E18.5 - DIFFERENCE DATA - PREGNANT VS NOT####
#ALTERNATIVE ROUTE - TAKING TWO TIME POINTS IN DATA AND COMPARING  
#Making data set:
first_last_data<-mice_pregnant_data[ mice_pregnant_data$Time_Point == 0 | mice_pregnant_data$Time_Point == 18.5, ]
View(first_last_data)

#Box plots of data from baseline and 18.5 Gestation to see what the data is like 
first_last_bp_plot<-ggplot(first_last_data, aes(x =Genotype, y = Systolic_Blood_Pressure, fill = Genotype))+
  stat_boxplot(geom='errorbar', width=.5)+
  geom_boxplot() +
  scale_colour_grey(start = 0.2, end = 0.5)+
  geom_point(shape = 18,
             col=ifelse(first_last_data$Time_Point==0,"black","grey30"),
             size=2.5) +
  labs(y= 'Systolic Blood Pressure (mmHg)', title='Day 0 vs Day 18.5')+
  theme_bw()
first_last_hr_plot<-ggplot(first_last_data, aes(x =Genotype, y = Heart_Rate, fill = Genotype))+
  stat_boxplot(geom='errorbar', width=.5)+
  geom_boxplot() +
  scale_colour_grey(start = 0.2, end = 0.5)+
  geom_point(shape = 18,
             col=ifelse(first_last_data$Time_Point==0,"black","grey30"),
             size=2.5) +
  labs(y= 'Heart Rate (bpm)', title='Day 0 vs Day 18.5')+
  theme_bw()
plot_grid(first_last_bp_plot,first_last_hr_plot, labels='AUTO')


####BASELINE AND E18.5 DATA POINTS ####
#Looking at data 
#What are the BLOOD PRESSURE means and SD for the data
summary(first_last_data)
mean(first_last_data[first_last_data$Time_Point==0, ]$Systolic_Blood_Pressure, na.rm=TRUE)
sd(first_last_data[first_last_data$Time_Point==0, ]$Systolic_Blood_Pressure, na.rm=TRUE)
mean(first_last_data[first_last_data$Time_Point==18.5, ]$Systolic_Blood_Pressure, na.rm=TRUE)
sd(first_last_data[first_last_data$Time_Point==18.5, ]$Systolic_Blood_Pressure, na.rm=TRUE)

#What are the HEART RATE means and SD for the data 
mean(first_last_data[first_last_data$Time_Point==0, ]$Heart_Rate, na.rm=TRUE)
sd(first_last_data[first_last_data$Time_Point==0, ]$Heart_Rate, na.rm=TRUE)
mean(first_last_data[first_last_data$Time_Point==18.5, ]$Heart_Rate, na.rm=TRUE)
sd(first_last_data[first_last_data$Time_Point==18.5, ]$Heart_Rate, na.rm=TRUE)


#### Normality and Homoscedasticity Test for both genotypes WT and KO in first to last time point data ####
#Blood Pressure - ALL OK
shapiro.test(first_last_data[first_last_data$Genotype=='WT',]$Systolic_Blood_Pressure) #p=0.3454
shapiro.test(first_last_data[first_last_data$Genotype=='KO',]$Systolic_Blood_Pressure) #p=0.4427
shapiro.test(first_last_data$Systolic_Blood_Pressure)#p=0.2302

#Heart Rate - ALL OK 
shapiro.test(first_last_data[first_last_data$Genotype=='WT',]$Heart_Rate)#p=0.2003
shapiro.test(first_last_data[first_last_data$Genotype=='KO',]$Heart_Rate)#p=0.6417
shapiro.test(first_last_data$Heart_Rate)#p=0.1637

#Levene Test 
#ALL >0.05 so good to go 
leveneTest(first_last_data$Systolic_Blood_Pressure~first_last_data$Genotype)
leveneTest(first_last_data$Heart_Rate~first_last_data$Genotype)



#### FIRST-> LAST ANOVA QUESTION 2 ####
#Does genetic background influence systolic blood pressure and heart rate in pregnancy 
#NOT SIGNIFICANT FOR EITHER -- this confirms aksels data for qu 2 (accounting for baseline) matching non significance 
first_last_bp_anova<-aov(first_last_data$Systolic_Blood_Pressure~first_last_data$Genotype+first_last_data$Time_Point)
summary(first_last_bp_anova)
first_last_hr_anova<-aov(first_last_data$Heart_Rate~first_last_data$Genotype+first_last_data$Time_Point)
summary(first_last_hr_anova)


#### FIRST-> LAST ANCOVA BP QUESTION 3 ####
#Does pregnancy influence systolic blood pressure and heart rate?
interactionAssumption<-lm(first_last_data$Systolic_Blood_Pressure~first_last_data$Genotype*first_last_data$Time_Point)
summary(interactionAssumption)

first_last_ancova<-lm(first_last_data$Systolic_Blood_Pressure~first_last_data$Genotype+first_last_data$Time_Point)
summary(interactionAssumption)

#NO SIGNIFICANCE 
first_last_AOVobject<-aov(Systolic_Blood_Pressure ~ Genotype+Time_Point, data=first_last_data)
summary(first_last_AOVobject)

#plot(first_last_data)

#### FIRST-> LAST ANCOVA HR QUESTION 3 ####
first_last_hr_AOVobject<-aov(Heart_Rate~ Genotype+Time_Point, data=first_last_data)
summary(first_last_AOVobject)

#### FIRST -> LAST QUESTION 4 - GLM ####
linearmodel1<-glm(Systolic_Blood_Pressure ~ Pregnant+Time_Point+Genotype, data= first_last_data)
summary(linearmodel1)
#plot(linearmodel1)
#this says we want pregnancy, time and genotype to explain heart rate 
linearmodel2<-glm(Heart_Rate~ Pregnant+Time_Point+Genotype, data=first_last_data)
summary(linearmodel2)


#### WILDTYPE AND KNOCKOUT BLOOD PRESSURE/ HEART RATE INFO BOXPLOT FOR MYSELF  ####
wildtype_bp<-mice_bp_data[mice_bp_data$Genotype=='WT', ]
#View(wildtype_bp)
wt1<-ggplot(wildtype_bp, aes(x=Time_Point, y=Systolic_Blood_Pressure))+
  geom_boxplot()+
  labs(y= 'Systolic Blood Pressure (mmHg)', x='Time Point', title= 'WT Systolic Blood Pressure (mmHg)')+
  theme_bw()

knockout_bp<-mice_bp_data[mice_bp_data$Genotype=='KO',]
#View(knockout_bp)
ko1<-ggplot(knockout_bp, aes(x=Time_Point, y=Systolic_Blood_Pressure))+
  geom_boxplot()+
  labs(y= 'Systolic Blood Pressure (mmHg)', x='Time Point', title= 'KO Systolic Blood Pressure (mmHg)')+
  theme_bw()

wildtype_hr<-mice_hr_data[mice_hr_data$Genotype=='WT', ]
#View(wildtype_hr)
wt2<-ggplot(wildtype_hr, aes(x=Time_Point, y=Heart_Rate))+
  geom_boxplot()+
  labs(y= 'Heart Rate (bpm)', x='Time Point', title= 'WT Heart Rate (bpm)')+
  theme_bw()

knockout_hr<-mice_hr_data[mice_hr_data$Genotype=='KO',]
#View(knockout_bp)
ko2<-ggplot(knockout_hr, aes(x=Time_Point, y=Heart_Rate))+
  geom_boxplot()+
  labs(y= 'Heart Rate (bpm)', x='Time Point', title= 'KO Heart Rate (bpm)')+
  theme_bw()

plot_grid(wt1,wt2, ko1, ko2, labels='AUTO')


####PLOTS OF DAY 0 Vs DAY 18.5 ####
#Scatter BP vs HR 
#Histogram of BP 
#Histogram of HR
p1<-ggplot(mice_pregnant_data, aes(x=Systolic_Blood_Pressure, y= Heart_Rate, color=Genotype))+
  geom_point()
p2<-ggplot(mice_pregnant_data, aes(x=Systolic_Blood_Pressure, color=Genotype))+
  geom_histogram()
p3<-ggplot(mice_pregnant_data, aes(x= Heart_Rate, color=Genotype))+
  geom_histogram()

ggarrange(p1,
          ggarrange(p2, p3, ncol=2, nrow=1, labels=c('B','C')),
          nrow=2,
          labels='A')





