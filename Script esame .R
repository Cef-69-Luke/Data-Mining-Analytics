#Librerie utilizzate------ 
#HMM#
library(mhsmm) #Hidden Markov Models 
library(ggplot2) #Grafici      
library(GGally) #Grafici correlazioni 
library(FourWayHMM)#Hidden Markov Models  
library(plotly)#plot
library(TSstudio)
library(xts)
library(MSwM)
#PCA#
library(FactoMineR)
library(FactoInvestigate)
library(Factoshiny)
library(factoextra)
library("corrplot") #esegue plot nelle dimensioni 
#REGRESSIONE MULTIVARIATA#
library(lme4) # for the analysis
library(haven) # to load the SPSS .sav file
library(tidyverse) # needed for data manipulation.
library(RColorBrewer) # needed for some extra colours in one of the graphs
library(lmerTest)# to get p-value estimations that are not part of the standard lme4 packages
#LIBRERIA SU INQUINAMENTO#
library(openair)

#Importazione dataset-----

uk_2015_20 <- importAURN(site = "kc1", year = 2015:2020, meta = TRUE) #Scaricamento Dati da libreria 

save(uk_2015_20, file = "Emisisoni UK 2015-2020.RData")#salvataggio in RData
load("~/Desktop/Esame Maurotti DataMining/Emisisoni UK 2015-2020.RData")#Dcaricamento dati 
uk_15_20 <- uk_2015_20[,-c(12,13,14,15,21)]#Tolgo variabili che non utilizzo 

#dati per HMM
logT_dta <- log(uk_15_20[,4:11]+1)#log trasformation x hmm ì
logT_dta$pm2.5[logT_dta$pm2.5 == -Inf] <- 0
logT_dta$pm10[logT_dta$pm10 == -Inf] <- 0


#Analisi preliminari -----
ggpairs(uk_15_20[,-c(1,2,3,15,16)]) #Correlazioni tra Inquinanti 
ggpairs(uk_15_20[,c(5:7)])# correlazioni ossidi 
ggpairs(uk_15_20[,c(4:7)])#corr variabili ntressanti 
ggpairs(uk_15_20[,c(8:11)])#corr variabili supp.

#SUMMARYPLOT
library(openair)
summaryPlot(uk_15_20)#plot sommario del data

#WINDROSE??
windRose(uk_15_20, type = "year", layout = c(3, 2))

#TIMEVARIATION
timeVariation(uk_15_20, pollutant = c("co","nox","no2",
                                        "no","so2","o3",
                                        "pm2.5","pm10"), type = "year")
#CALENDAR PLOT
#SO2
calendarPlot(uk_15_20, pollutant = "so2", year = 2015)
calendarPlot(uk_15_20, pollutant = "so2", year = 2016)
calendarPlot(uk_15_20, pollutant = "so2", year = 2017)
calendarPlot(uk_15_20, pollutant = "so2", year = 2018)
calendarPlot(uk_15_20, pollutant = "so2", year = 2019)
calendarPlot(uk_15_20, pollutant = "so2", year = 2020)

#o3
calendarPlot(uk_15_20, pollutant = "o3", year = 2015)
calendarPlot(uk_15_20, pollutant = "o3", year = 2016)
calendarPlot(uk_15_20, pollutant = "o3", year = 2017)
calendarPlot(uk_15_20, pollutant = "o3", year = 2018)
calendarPlot(uk_15_20, pollutant = "o3", year = 2019)
calendarPlot(uk_15_20, pollutant = "o3", year = 2020)

#PCA-----
library(FactoMineR)
library(FactoInvestigate)
library(Factoshiny)
library(factoextra)
library("corrplot")
#Prove PCA 
PCA(uk_15_20)
pca_uk1 <- PCA(uk_15_20, quali.sup = c(1,2,3), quanti.sup = (4:11))#prova 1 no 
pca_uk2 <- PCA(uk_15_20, quali.sup = c(1,2,3), quanti.sup = (8:11))#prova 2 da valutare 
pca_uk3 <- PCA(uk_15_20[,-c(1,2,3,8,9,10,11)])#prova 3 da valutare 
pca_uk4 <- PCA(uk_15_20, quali.sup = c(1,2,3), quanti.sup = c(8,9,10,11,12,13,15,16))#prova 4  Migliore di tutti


pca_uk5 <- PCA(uk_15_20[, -c(15,16)], quali.sup = c(1,2,3), quanti.sup = c(12,13,14)) #quali sup vetno e temperatura migliora 
res.PCA<-PCA(uk_15_20[, -c(15, 16)],quali.sup=c(1,2,3),quanti.sup=c(12,13,14),graph=FALSE)
plot.PCA(res.PCA,choix='var',col.quanti.sup='#0000FF')
plot.PCA(res.PCA,invisible=c('quali','ind.sup'),label =c('ind'))


# Rapp. dimensioni dati 
fviz_eig(pca_uk5, addlabels = TRUE)#Rapp. Dimensioni 

#Qualità dalle rappresentazioni 
fviz_cos2(pca_uk5, choice = "var") #Qualita delle rappresentazioni 

#ordine delle variabili per dimensioni 
corrplot(pca_uk5$var$cos2, is.corr=FALSE)#ordine delle dimensioni 

fviz_pca_var(pca_uk5, col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE # Avoid text overlappin
             )#Importanza aper cos2 

#Clustering delle varibili 
# Create a grouping variable using kmeans
# Create 3 groups of variables (centers = 3)
set.seed(123)
res.km <- kmeans(pca_uk5$var$coord, centers = 3, nstart = 25)
grp <- as.factor(res.km$cluster)
# Color variables by groups
fviz_pca_var(pca_uk5, col.var = grp, 
             palette = c("#0073C2FF", "#EFC000FF", "#868686FF"),
             legend.title = "Cluster")
#

fviz_pca_biplot(pca_uk5, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969"  # Individuals color
)
Investigate(pca_uk4)
PCAshiny(pca_uk4)
summary(pca_uk)
summary(pca_uk2)

pca_uk1 <- PCA(uk_15_20, quali.sup = c(1,2,3), quanti.sup = (12:16))



#PCR #principal component Regression------
#Installare pacchetti 
#install.packages("pls")
#install.packages("ISLR")
library(ISLR)
library(dplyr)
library(tidyr)
library(pls) #pcr
#GRAFICI FREQUENZE 
hist(uk_15_20$co) #p
hist(uk_15_20$nox) #p
hist(uk_15_20$no2) #g
hist(uk_15_20$no) #g
hist(uk_15_20$o3) #g
hist(uk_15_20$so2) #g
hist(uk_15_20$pm10) #g
hist(uk_15_20$pm2.5) #g

#link dell'esempio 
https://www.r-bloggers.com/2016/07/performing-principal-components-regression-pcr-in-r/
  

#Prova PCR----
set.seed (1234)
pcr_model <- pcr(nox ~ co + no + no2, data = uk_15_20, scale = TRUE, validation = "CV")
summary(pcr_model)
validationplot(pcr_model)
validationplot(pcr_model, val.type = "MSEP" )#convalida incrociata MSEP
validationplot(pcr_model, val.type = "R2")
predplot(pcr_model)#predizione modello 
coefplot(pcr_model)#coeff 



set.seed (56789)
pcr_model <- pcr(pm10 +pm2.5 ~ air_temp + wd + ws , data = uk_15_20, scale = TRUE, validation = "CV")
summary(pcr_model)
validationplot(pcr_model)
validationplot(pcr_model, val.type = "MSEP" )#convalida incrociata MSEP
validationplot(pcr_model, val.type = "R2")
predplot(pcr_model)#predizione modello 
coefplot(pcr_model)#coeff 

na_count <-sapply(uk_15_20, function(y) sum(length(which(is.na(y)))))
na_count <- data.frame(na_count)

