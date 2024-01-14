train<-read.csv("C:/Users/ugneo/OneDrive/Dokumentai/train_tiriamasis_cenz_2.csv") 
test<-read.csv("C:/Users/ugneo/OneDrive/Dokumentai/test_tiriamasis_cenz_2.csv")

# duomenu sutvarkymas -----------------------------------------------------

train$Age<-as.factor(train$Age)
test$Age<-as.factor(test$Age)

train$Tumour_Stage<-as.factor(train$Tumour_Stage)
test$Tumour_Stage<-as.factor(test$Tumour_Stage)

train$Histology<-as.factor(train$Histology)
test$Histology<-as.factor(test$Histology)

train$ER.status<-as.factor(train$ER.status)
test$ER.status<-as.factor(test$ER.status)

train$PR.status<-as.factor(train$PR.status)
test$PR.status<-as.factor(test$PR.status)

train$HER2.status<-as.factor(train$HER2.status)
test$HER2.status<-as.factor(test$HER2.status)

train$Surgery_type<-as.factor(train$Surgery_type)
test$Surgery_type<-as.factor(test$Surgery_type)

train$Patient_Status<-as.factor(train$Patient_Status)
test$Patient_Status<-as.factor(test$Patient_Status)

library(ggplot2)
library(ggpubr)
library(dplyr)
library(corrplot)
library(survival)
library(survminer)
library(TSHRC)
library(cmprsk)
library(mstate)

df <- read.csv("C:/Users/Austejos/Downloads/breast_sutvarkytas_2.csv",sep=",",header=TRUE)

View(df)


# PRADINE ANALIZE ---------------------------------------------------------


# dazniu lenteles

# ivykiai, cenzuravimas
table(df$Patient_Status)

# Age
table(df$Patient_Status, df$Age)

# Tumour_Stage
table(df$Patient_Status, df$Tumour_Stage)

# Histology
table(df$Patient_Status, df$Histology)

# ER.status
table(df$Patient_Status, df$ER.status)

# PR.status
table(df$Patient_Status, df$PR.status)

# HER2.status
table(df$Patient_Status, df$HER2.status)

# Surgery_type
table(df$Patient_Status, df$Surgery_type)

# staciakampes diagramos
mano_tema<-theme(
  axis.title.x = element_text(size = 14),
  axis.text.x = element_text(size = 12),
  axis.title.y = element_text(size = 14),
  axis.text.y = element_text(size = 12),
  title = element_text(size=14),
  plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), 
  axis.text = element_text(colour='black'))

df$Patient_Status <- as.factor(df$Patient_Status)
grupuoti<-df%>%group_by(Patient_Status)

# isskirtys
isskirtis <- function(x) {
  return(x <= quantile(x, .25) - 3*IQR(x) | x >= quantile(x, .75) + 3*IQR(x))
}


salygine_isskirtis<-function(x){
  return(((x <= quantile(x, .25) - 1.5*IQR(x))&(x > quantile(x, .25) - 3*IQR(x))) | (x >= quantile(x, .75) + 1.5*IQR(x))&
           (x < quantile(x, .75) + 3*IQR(x)))
}

isskirtys <- grupuoti %>% 
  select(Protein1, Protein2, Protein3, Protein4) %>% 
  mutate(outlier_Protein1=ifelse(isskirtis(Protein1), Protein1, NA))%>%
  mutate(salygine_outlier_Protein1=ifelse(salygine_isskirtis(Protein1), Protein1, NA)) %>% 
  mutate(outlier_Protein2=ifelse(isskirtis(Protein2), Protein2, NA))%>%
  mutate(salygine_outlier_Protein2=ifelse(salygine_isskirtis(Protein2), Protein2, NA)) %>% 
  mutate(outlier_Protein3=ifelse(isskirtis(Protein3), Protein3, NA))%>%
  mutate(salygine_outlier_Protein3=ifelse(salygine_isskirtis(Protein3), Protein3, NA)) %>%
  mutate(outlier_Protein4=ifelse(isskirtis(Protein4), Protein4, NA))%>%
  mutate(salygine_outlier_Protein4=ifelse(salygine_isskirtis(Protein4), Protein4, NA))

Protein1<-ggplot(grupuoti,
                 aes(x=Patient_Status, y=Protein1))+
  geom_boxplot()+ylim(c(-3, 3))+
  geom_point(data = subset(isskirtys, outlier_Protein1 != "NA"),
             aes(x = Patient_Status, y=Protein1), size = 2, color = "red") +
  ggtitle("Sta�iakamp�s pacien�i� I proteino lygio diagramos \npagal i�gyvenamum�") +
  labs(subtitle = "Pa�ym�ti kvantiliai")+
  scale_x_discrete(labels=c("0"="Cenz�ruotas","1"="Mirtis nuo v��io", "2"="Mirtis ne nuo v��io")) + mano_tema+
  xlab("I�gyvenamumo stadija") + ylab("I proteino lygis")+
  stat_summary(fun = "quantile", size = 4,
               geom = "text", aes(label = round(after_stat(y),2)),
               position = position_nudge(x = 0.5))

Protein2<-ggplot(grupuoti,
                 aes(x=Patient_Status, y=Protein2))+
  geom_boxplot()+ylim(c(-1, 4))+
  geom_point(data = subset(isskirtys, outlier_Protein2 != "NA"),
             aes(x = Patient_Status, y=Protein2), size = 2, color = "red") +
  ggtitle("Sta�iakamp�s pacien�i� II proteino lygio diagramos \npagal i�gyvenamum�") +
  labs(subtitle = "Pa�ym�ti kvantiliai")+
  scale_x_discrete(labels=c("0"="Cenz�ruotas","1"="Mirtis nuo v��io", "2"="Mirtis ne nuo v��io")) + mano_tema+
  xlab("I�gyvenamumo stadija") + ylab("II proteino lygis")+
  stat_summary(fun = "quantile", size = 4,
               geom = "text", aes(label = round(after_stat(y),2)),
               position = position_nudge(x = 0.5))

Protein3<-ggplot(grupuoti,
                 aes(x=Patient_Status, y=Protein3))+
  geom_boxplot()+ylim(c(-2, 3))+
  geom_point(data = subset(isskirtys, outlier_Protein3 != "NA"),
             aes(x = Patient_Status, y=Protein3), size = 2, color = "red") +
  ggtitle("Sta�iakamp�s pacien�i� III proteino lygio diagramos \npagal i�gyvenamum�") +
  labs(subtitle = "Pa�ym�ti kvantiliai")+
  scale_x_discrete(labels=c("0"="Cenz�ruotas","1"="Mirtis nuo v��io", "2"="Mirtis ne nuo v��io")) + mano_tema+
  xlab("I�gyvenamumo stadija") + ylab("III proteino lygis")+
  stat_summary(fun = "quantile", size = 4,
               geom = "text", aes(label = round(after_stat(y),2)),
               position = position_nudge(x = 0.5))

Protein4<-ggplot(grupuoti,
                 aes(x=Patient_Status, y=Protein4))+
  geom_boxplot()+ylim(c(-2.5, 2))+
  geom_point(data = subset(isskirtys, outlier_Protein4 != "NA"),
             aes(x = Patient_Status, y=Protein4), size = 2, color = "red") +
  ggtitle("Sta�iakamp�s pacien�i� IV proteino lygio diagramos \npagal i�gyvenamum�") +
  labs(subtitle = "Pa�ym�ti kvantiliai")+
  scale_x_discrete(labels=c("0"="Cenz�ruotas","1"="Mirtis nuo v��io", "2"="Mirtis ne nuo v��io")) + mano_tema+
  xlab("I�gyvenamumo stadija") + ylab("IV proteino lygis")+
  stat_summary(fun = "quantile", size = 4,
               geom = "text", aes(label = round(after_stat(y),2)),
               position = position_nudge(x = 0.5))

ggarrange(Protein1, Protein2, Protein3, Protein4, ncol=2, nrow = 2)

# koreliacijos
cor_df <- df %>% 
  select(Protein1, Protein2, Protein3, Protein4)
corelation <- cor(cor_df)
corrplot(corelation, method = 'color', order = 'alphabet', tl.col = "black")


# KM KREIVES --------------------------------------------------------------

train <- read.csv("C:/Users/Austejos/Downloads/train_tiriamasis_cenz_2.csv",sep=",",header=TRUE)

str(train)
train$Age<-as.factor(train$Age)
train$Protein1<-as.numeric(train$Protein1)
train$Protein2<-as.numeric(train$Protein2)
train$Protein3<-as.numeric(train$Protein3)
train$Protein4<-as.numeric(train$Protein4)
train$Tumour_Stage<-as.factor(train$Tumour_Stage)
train$Histology<-as.factor(train$Histology)
train$ER.status<-as.factor(train$ER.status)
train$PR.status<-as.factor(train$PR.status)
train$HER2.status<-as.factor(train$HER2.status)
train$Surgery_type<-as.factor(train$Surgery_type)
train$Patient_Status<-as.factor(train$Patient_Status)
train$days_of_last_follow_up<-as.numeric(train$days_of_last_follow_up)
train$days_of_last_follow_up<-ifelse(train$days_of_last_follow_up <0, 0, train$days_of_last_follow_up)

mano_tema<-theme(
  axis.title.x = element_text(size = 12),
  axis.text.x = element_text(size = 12),
  axis.title.y = element_text(size = 12),
  axis.text.y = element_text(size = 12),
  title = element_text(size=13),
  plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), 
  axis.text = element_text(colour='black'))

# neskisrtome pagal joki pozymi

KM_1 <- survfit(Surv(time = days_of_last_follow_up, 
                     event = Patient_Status == 1 ) ~ 1, 
                data = train)
summary(KM_1)
print(KM_1, print.rmean = T)

g1<-ggsurvplot(KM_1,
               conf.int=TRUE, # add confidence intervals
               risk.table=c( "nrisk_cumcensor"), # show a risk table below the plot
               palette=c("dodgerblue4"),
               title="Kaplan-Meier kreiv? i?gyvenamumui", # add title to plot
               risk.table.height=.2,
               legend="none",
               legend.labs="Visi",
               xlab="Laikas (dienomis)",
               ylab="I?gyvenamumo tikimyb?",
               legend.title="",
               risk.table.title = "Pacien?i?, esan?i? rizikos grup?je, skai?ius \n(cenz?ruot? pacien?i? sukauptinis da?nis)",
               ggtheme = mano_tema,surv.median.line = "hv")  
g1

# skirstome pagal amziu
unique(train$Age) # 6 grupes

KM_age <- survfit(Surv(time = days_of_last_follow_up, 
                       event = Patient_Status == 1 ) ~ Age, 
                  data = train)
summary(KM_age)
print(KM_age, print.rmean = T)


g_age<-ggsurvplot(KM_age,
                  conf.int=TRUE, # add confidence intervals
                  risk.table=c( "nrisk_cumcensor"), # show a risk table below the plot
                  palette=c("navy", "deeppink", "lightgreen", "orangered", "turquoise", "darkred"),
                  title="Kaplan-Meier kreiv? i?gyvenamumui pagal am?iaus kategorij?", # add title to plot
                  risk.table.height=.2,
                  legend.labs=c("(28,40]", "(40,50]", "(50,60]", "(60,70]", "(70,80]", "(80,91]"),
                  xlab="Laikas (dienomis)",
                  ylab="I?gyvenamumo tikimyb?",
                  legend.title="Am?iaus kategorija",
                  risk.table.title = "Pacien?i?, esan?i? rizikos grup?je, skai?ius \n(cenz?ruot? pacien?i? sukauptinis da?nis)",
                  ggtheme = mano_tema,surv.median.line = "hv")  
g_age

# skirstome pagal stadija
unique(train$Tumour_Stage) # 3 grupes

KM_Tumour_Stage <- survfit(Surv(time = days_of_last_follow_up, 
                                event = Patient_Status == 1 ) ~ Tumour_Stage, 
                           data = train)
summary(KM_Tumour_Stage)
print(KM_Tumour_Stage, print.rmean = T)


g_Tumour_Stage<-ggsurvplot(KM_Tumour_Stage,
                           conf.int=TRUE, # add confidence intervals
                           risk.table=c( "nrisk_cumcensor"), # show a risk table below the plot
                           palette=c("navy", "deeppink", "lightgreen"),
                           title="Kaplan-Meier kreiv? i?gyvenamumui pagal stadij?", # add title to plot
                           risk.table.height=.3,
                           legend.labs=c("I", "II", "III"),
                           xlab="Laikas (dienomis)",
                           ylab="I?gyvenamumo tikimyb?",
                           legend.title="Stadija",
                           risk.table.title = "Pacien?i?, esan?i? rizikos grup?je, skai?ius \n(cenz?ruot? pacien?i? sukauptinis da?nis)",
                           ggtheme = mano_tema,surv.median.line = "hv")  
g_Tumour_Stage

# grupes kertasi

# skirstome pagal histologija
unique(train$Histology) # 3 grupes

KM_Histology <- survfit(Surv(time = days_of_last_follow_up, 
                             event = Patient_Status == 1 ) ~ Histology, 
                        data = train)
summary(KM_Histology)
print(KM_Histology, print.rmean = T)


g_Histology<-ggsurvplot(KM_Histology,
                        conf.int=TRUE, # add confidence intervals
                        risk.table=c( "nrisk_cumcensor"), # show a risk table below the plot
                        palette=c("navy", "deeppink", "lightgreen"),
                        title="Kaplan-Meier kreiv? i?gyvenamumui pagal histologij? (v??io tip?)", # add title to plot
                        risk.table.height=.25,
                        legend.labs=c("Infiltruojanti duktalin? karcinoma", "Infiltruojanti lobulin? karcinoma", "Mucinozin? karcinoma"),
                        xlab="Laikas (dienomis)",
                        ylab="I?gyvenamumo tikimyb?",
                        legend.title="Histologija",
                        risk.table.title = "Pacien?i?, esan?i? rizikos grup?je, skai?ius \n(cenz?ruot? pacien?i? sukauptinis da?nis)",
                        ggtheme = mano_tema,surv.median.line = "hv")  
g_Histology

# grupes kertasi

# skirstome pagal HER2 statusa
unique(train$HER2.status) # 2 grupes

KM_HER2.status <- survfit(Surv(time = days_of_last_follow_up, 
                               event = Patient_Status == 1 ) ~ HER2.status, 
                          data = train)
summary(KM_HER2.status)
print(KM_HER2.status, print.rmean = T)


g_HER2.status<-ggsurvplot(KM_HER2.status,
                          conf.int=TRUE, # add confidence intervals
                          risk.table=c( "nrisk_cumcensor"), # show a risk table below the plot
                          palette=c("navy", "deeppink", "lightgreen"),
                          title="Kaplan-Meier kreiv? i?gyvenamumui pagal HER 2 status?", # add title to plot
                          risk.table.height=.25,
                          legend.labs=c("Neigiamas", "Teigiamas"),
                          xlab="Laikas (dienomis)",
                          ylab="I?gyvenamumo tikimyb?",
                          legend.title=" HER 2 statusas",
                          risk.table.title = "Pacien?i?, esan?i? rizikos grup?je, skai?ius \n(cenz?ruot? pacien?i? sukauptinis da?nis)",
                          ggtheme = mano_tema,surv.median.line = "hv")  
g_HER2.status

her2 <- train
her2$Patient_Status<-as.numeric(her2$Patient_Status)
her2$Patient_Status<-ifelse(her2$Patient_Status == 2, 1, 0)
her2$HER2.status<-as.numeric(her2$HER2.status)
her2$HER2.status<-ifelse(her2$HER2.status == 2, 1, 0)
her2$days_of_last_follow_up<-as.numeric(her2$days_of_last_follow_up)
twostage(her2$days_of_last_follow_up, her2$Patient_Status, her2$HER2.status, nboot = 1000)

# grupes kertasi

# skirstome pagal operacijos tipa
unique(train$Surgery_type) # 4 grupes

KM_Surgery_type <- survfit(Surv(time = days_of_last_follow_up, 
                                event = Patient_Status == 1 ) ~ Surgery_type, 
                           data = train)
summary(KM_Surgery_type)
print(KM_Surgery_type, print.rmean = T)


g_Surgery_type<-ggsurvplot(KM_Surgery_type,
                           conf.int=TRUE, # add confidence intervals
                           risk.table=c( "nrisk_cumcensor"), # show a risk table below the plot
                           palette=c("navy", "deeppink", "darkred",  "turquoise"),
                           title="Kaplan-Meier kreiv? i?gyvenamumui pagal operacijos tip?", # add title to plot
                           risk.table.height=.25,
                           legend.labs=c("Kita", " Lumpektomija", "Modifikuota radikali mastektomija", "Paprasta mastektomija"),
                           xlab="Laikas (dienomis)",
                           ylab="I?gyvenamumo tikimyb?",
                           legend.title="Operacijos tipas",
                           risk.table.title = "Pacien?i?, esan?i? rizikos grup?je, skai?ius \n(cenz?ruot? pacien?i? sukauptinis da?nis)",
                           ggtheme = mano_tema,surv.median.line = "hv")  
g_Surgery_type

# grupes kertasi


# 1 - KM ------------------------------------------------------------------

train_new <- within(train, {
  status.cancer <- as.numeric({Patient_Status == 1})
  status.other <- as.numeric({Patient_Status == 2})})


result.cancer.km <- survfit(Surv(time = days_of_last_follow_up, 
                                 event = status.cancer) ~ 1, 
                            data = train_new)

result.other.km <- survfit(Surv(time = days_of_last_follow_up,
                                event = status.other) ~ 1,
                           data = train_new)

surv.other.km <- result.other.km$surv
time.km <- result.other.km$time/365
# CIF arba CDF, t.y. 1-Survival:
surv.cancer.km <- result.cancer.km$surv
cumDist.cancer.km <- 1 - surv.cancer.km

##########################################

windows(height=5, width=7)
layout(matrix(c(2,1,3), ncol=3),
       heights=c(5,5,5), widths=c(0.5, 5, 0.5))
layout.show(3)
par(mar=c(5,2,2,2), cex.axis=1.65, cex.lab=1.65)
plot(cumDist.cancer.km ~ time.km, type="s", ylim=c(0,1), xlim = c(0,8), lwd=2, xlab="Metai nuo kr?ties v??io diagnoz?s nustatymo",
     axes=F, col="blue")
lines(surv.other.km ~ time.km, type="s", col="green", lwd=2)
axis(1)
axis(2, at=seq(0,1,0.2), las=1)
#axis(3)
axis(4, at=seq(0,1,0.2), labels=seq(1,0,-0.2), las=1)
box()
text(0.75, 0.25, label="Mirtis nuo\nkr?ties v??io", cex=1.65)
text(0.75, 0.73, label="Mirtis ne nuo\nkr?ties v??io", cex=1.65)

x <- c(0,5)
y <- c(0,0.25)
par(mar=c(0,0,0,0))
plot(x ~ y,type="n",axes=F)
text(x=0.1,y=2.5,"Tikimyb? numirti nuo kr?ties v??io",srt=90,cex=1.65)
plot(x ~ y,type="n",axes=F)
text(x=0.15,y=2.5,"Tikimyb? numirti ne nuo kr?ties v??io",srt=-90,cex=1.65)


##############################################
ci.cancer <- Cuminc(time=train$days_of_last_follow_up, 
                    status=train$Patient_Status)
head(ci.cancer)


ci1 <- ci.cancer$CI.1  # CI.1 is for prostate cancer 
ci2 <- ci.cancer$CI.2  # CI.2 is for other causes 
times <- ci.cancer$time/365  # convert months to years 
Rci2 <- 1 - ci2

###############################################

windows(height=5, width=7)
layout(matrix(c(2,1,3), ncol=3),
       heights=c(5,5,5), widths=c(0.5, 5, 0.5))
layout.show(3)
par(mar=c(5,2,2,2), cex.axis=1.65, cex.lab=1.65)

#par(mar=c(5, 4, 4, 4) + 0.1)
plot(cumDist.cancer.km ~ time.km, type="s", ylim=c(0,1), lwd=1, xlab="Years from prostate cancer diagnosis",
     axes=F, col="lightblue", xlim=c(0,10))
lines(surv.other.km ~ time.km, type="s", col="lightgreen", lwd=1)
axis(1)
axis(2, at=seq(0,1,0.2), las=1)
#axis(3)
axis(4, at=seq(0,1,0.2), labels=seq(1,0,-0.2), las=1)
box()
text(20/12, 0.27, label="Death from\nprostate cancer", cex=1.65)
text(20/12, 0.73, label="Death from\nother causes", cex=1.65)

lines(Rci2 ~ times, type="s", ylim=c(0,1), lwd=2, col="green")
lines(ci1 ~ times, type="s", lwd=2, col="blue")

x <- c(0,5)
y <- c(0,0.25)
par(mar=c(0,0,0,0))
plot(x ~ y,type="n",axes=F)
text(x=0.1,y=2.5,"Probability of death from prostate cancer",srt=90,cex=1.65)
plot(x ~ y,type="n",axes=F)
text(x=0.15,y=2.5,"Probability of death from other causes",srt=-90,cex=1.65)


##############################################################

windows(height=5, width=7)

par(mar=c(5, 4, 2, 2) + 0.1, cex.lab=1.1, cex.axis=1.1)  # default

plot(ci1 ~ times, type="s", ylim=c(0,1), lwd=2, xlab="Metai nuo kr?ties v??io diagnoz?s nustatymo",
     ylab="Tikimyb?, jog pacient? mir?", axes=F, col="blue", xlim=c(0,8))
ci1.ci2 <- ci1 + ci2
lines(ci1.ci2 ~ times, type="s", lwd=2, col="green")

axis(1)
axis(2, at=seq(0,1,0.2), labels=seq(0,1,0.2), las=1)
#axis(3)
#axis(4)
box()
text(7, 0.18, "Mirtis nuo\nkr?ties v??io", cex=1.1)
text(7, 0.78, "Mirtis nuo\nkitos prie?asties", cex=1.1)



# CIF ---------------------------------------------------------------------

# Age
cuminc(Surv(days_of_last_follow_up, Patient_Status) ~ Age, train) %>%
  ggcuminc(outcome = "1") +
  labs(
    x = "Days"
  ) +
  scale_ggsurvfit()

# Tumour_Stage
cuminc(Surv(days_of_last_follow_up, Patient_Status) ~ Tumour_Stage, train) %>%
  ggcuminc(outcome = c("1","2"),theme = mano_tema, linewidth = 0.8) +
  labs(
    x = "Laikas (dienomis)",
    y = "CIF",
    title = "CIF funkcija pagal stadijos tip?"
  ) +
  scale_ggsurvfit()

cuminc(Surv(days_of_last_follow_up, Patient_Status) ~ Tumour_Stage, train) %>% 
  tbl_cuminc(
    times = 1095.75, 
    label_header = "**{time/365.25}-year cuminc**") %>% 
  add_p()

# Histology
cuminc(Surv(days_of_last_follow_up, Patient_Status) ~ Histology, train) %>%
  ggcuminc(outcome = c("1","2"),theme = mano_tema, linewidth = 0.8) +
  labs(
    x = "Laikas (dienomis)",
    y = "CIF",
    title = "CIF funkcija pagal histologij?"
  ) +
  scale_ggsurvfit()

cuminc(Surv(days_of_last_follow_up, Patient_Status) ~ Histology, train) %>% 
  tbl_cuminc(
    times = 1095.75, 
    label_header = "**{time/365.25}-year cuminc**") %>% 
  add_p()

# HER2.status
cuminc(Surv(days_of_last_follow_up, Patient_Status) ~ HER2.status, train) %>%
  ggcuminc(outcome = c("1","2"),theme = mano_tema, linewidth = 0.8) +
  labs(
    x = "Laikas (dienomis)",
    y = "CIF",
    title = "CIF funkcija pagal HER 2 status?"
  ) +
  scale_ggsurvfit()

cuminc(Surv(days_of_last_follow_up, Patient_Status) ~ HER2.status, train) %>% 
  tbl_cuminc(
    times = 1095.75, 
    label_header = "**{time/365.25}-year cuminc**") %>% 
  add_p()

# Surgery_type
cuminc(Surv(days_of_last_follow_up, Patient_Status) ~ Surgery_type, train) %>%
  ggcuminc(outcome = c("1","2"),theme = mano_tema, linewidth = 0.8) +
  labs(
    x = "Laikas (dienomis)",
    y = "CIF",
    title = "CIF funkcija pagal operacijos tip?"
  ) +
  scale_ggsurvfit()

cuminc(Surv(days_of_last_follow_up, Patient_Status) ~ Surgery_type, train) %>% 
  tbl_cuminc(
    times = 1095.75, 
    label_header = "**{time/365.25}-year cuminc**") %>% 
  add_p()

cuminc(Surv(days_of_last_follow_up, Patient_Status) ~ 1, train) %>%
  ggcuminc(outcome = c("1","2"),theme = mano_tema, linewidth = 0.8) +
  labs(
    x = "Laikas (dienomis)",
    y = "CIF",
    title = "CIF funkcija"
  ) +
  scale_ggsurvfit()


# MODELIS SU cmprsk PAKETU ------------------------------------------------

library(cmprsk)

library(survival)
m1 <- coxph(Surv(days_of_last_follow_up, Patient_Status==1) ~ Age+Protein1+Protein2+Protein3+Protein4+
              Tumour_Stage+Histology+HER2.status+Surgery_type, data = train)
m1_test<-coxph(Surv(days_of_last_follow_up, Patient_Status==1) ~ Age+Protein1+Protein2+Protein3+Protein4+
                 Tumour_Stage+Histology+HER2.status+Surgery_type, data = test)
modelis1 <- with(train, crr(days_of_last_follow_up, Patient_Status, cov1 = model.matrix(m1), failcode = 1, cencode = 0))
summary(modelis1, Exp = T)

m2 <- coxph(Surv(days_of_last_follow_up, Patient_Status==2) ~ Age+Protein1+Protein2+Protein3+Protein4+
              Tumour_Stage+Histology+HER2.status+Surgery_type, data = train)
m2_test<- coxph(Surv(days_of_last_follow_up, Patient_Status==2) ~ Age+Protein1+Protein2+Protein3+Protein4+
                  Tumour_Stage+Histology+HER2.status+Surgery_type, data = test)
modelis2 <- with(train, crr(days_of_last_follow_up, Patient_Status, cov1 = model.matrix(m2), failcode = 2, cencode = 0))
summary(modelis2, Exp = T)



# multikolinearumas -------------------------------------------------------

summary(with(train, crr(days_of_last_follow_up, Patient_Status, cov1 = model.matrix(m1)[,6], failcode = 1, cencode = 0)))#ok
summary(with(train, crr(days_of_last_follow_up, Patient_Status, cov1 = model.matrix(m1)[,7], failcode = 1, cencode = 0)))#ok
summary(with(train, crr(days_of_last_follow_up, Patient_Status, cov1 = model.matrix(m1)[,8], failcode = 1, cencode = 0)))#ok
summary(with(train, crr(days_of_last_follow_up, Patient_Status, cov1 = model.matrix(m1)[,9], failcode = 1, cencode = 0)))#ok

summary(with(train, crr(days_of_last_follow_up, Patient_Status, cov1 = model.matrix(m2)[,6], failcode = 2, cencode = 0)))#ok
summary(with(train, crr(days_of_last_follow_up, Patient_Status, cov1 = model.matrix(m2)[,7], failcode = 2, cencode = 0)))#ok
summary(with(train, crr(days_of_last_follow_up, Patient_Status, cov1 = model.matrix(m2)[,8], failcode = 2, cencode = 0)))#ok
summary(with(train, crr(days_of_last_follow_up, Patient_Status, cov1 = model.matrix(m2)[,9], failcode = 2, cencode = 0)))#netinka



# PH ----------------------------------------------------------------------
library(goftte)

train_1<-train[train$Patient_Status==1,]

par(mfrow=c(2, 2))
for (j in 1:17) {
  plot(modelis1$res[, j],
       ylab=colnames(model.matrix(m1))[j],
       xlab="Event = 1 time index")
  abline(h=0, lty=2)
}

train_2<-train[train$Patient_Status==2,]

par(mfrow=c(2, 2))
for (j in 1:17) {
  plot(modelis2$res[, j],
       ylab=colnames(model.matrix(m2))[j],
       xlab="Event = 2 time index")
  abline(h=0, lty=2)
}


# pazingsnine -----------------------------------------------------------------
fgr1<-FGR(Hist(days_of_last_follow_up, Patient_Status)~ Age+Protein1+Protein2+Protein3+Protein4+
            Tumour_Stage+Histology+HER2.status+Surgery_type,data=train, cause=1)#ismetame protein2

summary(FGR(Hist(days_of_last_follow_up, Patient_Status)~ Age+Protein1+Protein3+Protein4+
              Tumour_Stage+Histology+HER2.status+Surgery_type,data=train, cause=1))#ismetame age

summary(FGR(Hist(days_of_last_follow_up, Patient_Status)~ Protein1+Protein3+Protein4+
              Tumour_Stage+Histology+HER2.status+Surgery_type,data=train, cause=1))#ismetame surgery

summary(FGR(Hist(days_of_last_follow_up, Patient_Status)~ Protein1+Protein3+Protein4+
              Tumour_Stage+Histology+HER2.status,data=train, cause=1))#ismetame protein3

summary(FGR(Hist(days_of_last_follow_up, Patient_Status)~ Protein1+Protein4+
              Tumour_Stage+Histology+HER2.status,data=train, cause=1))#ismetame protein1

summary(FGR(Hist(days_of_last_follow_up, Patient_Status)~ Protein4+
              Tumour_Stage+Histology+HER2.status,data=train, cause=1))#ismetame histology

summary(FGR(Hist(days_of_last_follow_up, Patient_Status)~ Protein4+
              Tumour_Stage+HER2.status,data=train, cause=1))#ismetame her2 status

summary(FGR(Hist(days_of_last_follow_up, Patient_Status)~ Protein4+
              Tumour_Stage,data=train, cause=1))#visi reiksmingi
modelis1_reiksmingas<-FGR(Hist(days_of_last_follow_up, Patient_Status)~ Protein4+
                            Tumour_Stage,data=train, cause=1)

#2 event
summary(modelis2, Exp = T)#ismetame protein4

summary( with(train, crr(days_of_last_follow_up, Patient_Status, cov1 = model.matrix(m2)[,c(-9)], 
                         failcode = 2, cencode = 0)), Exp=T)#trumor stage

summary( with(train, crr(days_of_last_follow_up, Patient_Status, cov1 = model.matrix(m2)[,c(-9,-10,-11)], 
                         failcode = 2, cencode = 0)), Exp=T)#protein1

summary( with(train, crr(days_of_last_follow_up, Patient_Status, cov1 = model.matrix(m2)[,c(-6,-9,-10,-11)], 
                         failcode = 2, cencode = 0)), Exp=T)#age

summary( with(train, crr(days_of_last_follow_up, Patient_Status, cov1 = model.matrix(m2)[,c(-1,-2,-3,-4,-5,-6,-9,-10,-11)], 
                         failcode = 2, cencode = 0)), Exp=T)#histology


summary( with(train, crr(days_of_last_follow_up, Patient_Status, cov1 = model.matrix(m2)[,c(7,8,14,15,16,17)], 
                         failcode = 2, cencode = 0)), Exp=T)#protein2

summary( with(train, crr(days_of_last_follow_up, Patient_Status, cov1 = model.matrix(m2)[,c(8,14,15,16,17)], 
                         failcode = 2, cencode = 0)), Exp=T)#prtotein3

crr2<-with(train, crr(days_of_last_follow_up, Patient_Status, cov1 = model.matrix(m2)[,c(14,15,16,17)], 
                      failcode = 2, cencode = 0))
summary( with(train, crr(days_of_last_follow_up, Patient_Status, cov1 = model.matrix(m2)[,c(14,15,16,17)], 
                         failcode = 2, cencode = 0)), Exp=T)#HER2 ir surgery kai alfa0,15


# HARRELIO ----------------------------------------------------------------
library(riskRegression)
library(pec)
library(cmprsk)

fgr1<-FGR(Hist(days_of_last_follow_up, Patient_Status)~ Protein4+
            Tumour_Stage,data=train, cause=1)

fgr2<-FGR(Hist(days_of_last_follow_up, Patient_Status)~HER2.status+
            Surgery_type,data=train, cause=2)

pec::cindex(list(fgr1), formula = Hist(days_of_last_follow_up, Patient_Status)~ 1,
            data = test, cause = 1, 
            cens.model = "marginal")
206/381#0.54

pec::cindex(fgr1, formula = Hist(days_of_last_follow_up, Patient_Status)~ 1,
            data = train, cause = 1, 
            cens.model = "marginal")
4077/6242#0.65

pec::cindex(fgr2, formula = Hist(days_of_last_follow_up, Patient_Status)~ 1,
            data = test, cause = 2)
173/255#0.68
pec::cindex(fgr2, formula = Hist(days_of_last_follow_up, Patient_Status)~ 1,
            data = train, cause = 2)
4286/6842#0.63


# predict risk ------------------------------------------------------------

p1<-predictRisk(fgr1,times=test$days_of_last_follow_up,newdata=test)
p2<-predictRisk(fgr2,times=test$days_of_last_follow_up,newdata=test)


# calplot -----------------------------------------------------------------
par(mfrow=c(1, 1))
c<-calPlot(fgr1,cause=1, bars=T, q=5,hanging=F,
           formula=Hist(days_of_last_follow_up, Patient_Status)~ 1, data=train, xlim=c(0, 0.25), ylim=c(0,0.25))
calPlot(fgr1,cause=1, bars=T, q=5,
        formula=Hist(days_of_last_follow_up, Patient_Status)~ 1, data=test, xlim=c(0, 0.5), ylim=c(0,0.5))

calPlot(fgr1,cause=1,legend=F,
        formula=Hist(days_of_last_follow_up, Patient_Status)~ 1, data=train, xlim=c(0, 0.25), ylim=c(0,0.25))
calPlot(fgr1,cause=1, legend=F,
        formula=Hist(days_of_last_follow_up, Patient_Status)~ 1, data=test, xlim=c(0, 0.25), ylim=c(0,0.25))

calPlot(fgr2,cause=2,legend = F,
        formula=Hist(days_of_last_follow_up, Patient_Status)~ 1, data=train, xlim=c(0, 0.3), ylim=c(0,0.3))

calPlot(fgr2,cause=2, legend = F,
        formula=Hist(days_of_last_follow_up, Patient_Status)~ 1, data=test, xlim=c(0, 0.3), ylim=c(0,0.3))


library(splitstackshape)
library(openxlsx)
library(dynpred)
library(Epi)
library(rms)
library(MASS)
library(survminer)
library(survival)
library(ggplot2)
library(survminer)
library(car)
library(QHScrnomo)
library(pROC)


# IVYKIS 1 ----------------------------------------------------------------
# modeliu sukurimai -------------------------------------------------------
cox_univ_amzius<-coxph(Surv(time = days_of_last_follow_up ,event = Patient_Status == 1)~ Age,
                       data =  train) 

cox_univ_protein1<-coxph(Surv(time = days_of_last_follow_up ,event = Patient_Status == 1)~ Protein1,
                         data =  train)  

cox_univ_protein2<-coxph(Surv(time = days_of_last_follow_up ,event = Patient_Status == 1)~ Protein2,
                         data =  train)  

cox_univ_protein3<-coxph(Surv(time = days_of_last_follow_up ,event = Patient_Status == 1)~ Protein3,
                         data =  train)  

cox_univ_protein4<-coxph(Surv(time = days_of_last_follow_up ,event = Patient_Status == 1)~ Protein4,
                         data =  train)  

cox_univ_tumour<-coxph(Surv(time = days_of_last_follow_up ,event = Patient_Status == 1)~ Tumour_Stage,
                       data =  train) 

cox_univ_histology<-coxph(Surv(time = days_of_last_follow_up ,event = Patient_Status == 1)~ Histology,
                          data =  train)   

cox_univ_her<-coxph(Surv(time = days_of_last_follow_up ,event = Patient_Status == 1)~ HER2.status,
                    data =  train)  

cox_univ_sur<-coxph(Surv(time = days_of_last_follow_up ,event = Patient_Status == 1)~ Surgery_type,
                    data =  train)  

# PRIELAIDU PATIKRINIMAS --------------------------------------------------
# isskirtys ---------------------------------------------------------------
hold<-2/sqrt(nrow(train))
hold

dfbeta_amzius<-residuals(cox_univ_amzius, type="dfbeta")
dfbeta_protein1<-residuals(cox_univ_protein1, type="dfbeta")
dfbeta_protein2<-residuals(cox_univ_protein2, type="dfbeta")
dfbeta_protein3<-residuals(cox_univ_protein3, type="dfbeta")
dfbeta_protein4<-residuals(cox_univ_protein4, type="dfbeta")
dfbeta_tumour<-residuals(cox_univ_tumour, type="dfbeta")
dfbeta_histology<-residuals(cox_univ_histology, type="dfbeta")
dfbeta_her<-residuals(cox_univ_her, type="dfbeta")
dfbeta_sur<-residuals(cox_univ_sur, type="dfbeta")

vardai<-c("Amžius", "I proteino lygis","II proteino lygis", 
          "III proteino lygis", "IV proteino lygis", 
          "Stadija","Histologija",
          "HER2 statusas","Operacija")
dfbeta<-data.frame(dfbeta_amzius, dfbeta_protein1, dfbeta_protein2, dfbeta_protein3, dfbeta_protein4,
                   dfbeta_tumour, dfbeta_histology, dfbeta_her,  dfbeta_sur)

par(mfrow=c(2, 3))
for (j in 7:9) { #1:4 - amzius, protein1, protein2, protein3; 3:8 - protein4; 4:9 - operacija
  plot(dfbeta[, j],
       ylab=vardai[j], xlab="")
  title(xlab="Indeksas")
  abline(h=0, lty=2)
  abline(h = hold, lty = 2, col="red")
  abline(h = -hold, lty=2, col="red")
}#yra isskirciu

#kurie indeksai yra iskirtis
which(abs(dfbeta[, 1]) > hold)
which(abs(dfbeta[, 2]) > hold)
which(abs(dfbeta[, 3]) > hold)
which(abs(dfbeta[, 4]) > hold)
which(abs(dfbeta[, 5]) > hold)
which(abs(dfbeta[, 9]) > hold)

#pasalinti tas isskirtis
train <- train[-c(89, 137, 193, 223, 252, 17, 55, 140, 162, 26), ]

# kov. reiksmingumas ir ph prielaida --------------------------------------

cox_univ <- coxph(Surv(time = days_of_last_follow_up ,event = Patient_Status == 1)~ Age+Protein1+Protein2+Protein3+Protein4+
                    Tumour_Stage+Histology+HER2.status+Surgery_type,
                  data =  train) 
summary(cox_univ)
cox.zph(cox_univ)

#multikolinerumas 
cox_kiek <- coxph(Surv(time = days_of_last_follow_up ,event = Patient_Status == 2)~ Protein1+Protein2+Protein3+Protein4,
                  data =  train)
vif(cox_kiek)

# tiesiskumo prielaida ----------------------------------------------------
#protein4
par(mfrow=c(1, 2))
res_univ_protein4 <- residuals(cox_univ_protein4, type="martingale")
plot(train$Protein4, res_univ_protein4, xlab="IV proteino lygis",
     ylab="Martingaliosios liekanos")
abline(h=0, lty=2)
lines(lowess(train$Protein4, res_univ_protein4, iter=0), col="red")
b <- coef(cox_univ_protein4) # regression coefficients
plot(train$Protein4, b*train$Protein4 + res_univ_protein4, xlab="Protein4",
     ylab="Komponentė + martingaliosios liekanos")
abline(lm(b*train$Protein4+ res_univ_protein4 ~ train$Protein4), lty=2)
lines(lowess(train$Protein4, b*train$Protein4 + res_univ_protein4, iter=0), col="red")#tenkina,
#tiesiskumo prielaida tenkinama

# modelis -----------------------------------------------------------------
cox<-coxph(Surv(time = days_of_last_follow_up ,event = Patient_Status == 1 )~  Protein4+Tumour_Stage,
           data =  train)

summary(cox)
#konkordacija 0.679 
AIC(cox)
#AIC - 361.9418

stepAIC(cox_univ, method="both")

# kalibravimas ------------------------------------------------------------
fs1=coxph(Surv(time = days_of_last_follow_up ,event = Patient_Status == 1 )~Protein4+Tumour_Stage,data=train,x=1)
xs=Score(list(Cox1=fs1),Surv(time = days_of_last_follow_up ,event = Patient_Status == 1 )~1,data=train,
         plots="cal",metrics=NULL)
plotCalibration(xs, xlim = c(0, 0.4))
#plotCalibration(xs,cens.method="local",pseudo=1)
#plotCalibration(xs, method = "quantile", cens.method = "jackknife", breaks = seq(0, 1, 0.1))

fs1=coxph(Surv(time = days_of_last_follow_up ,event = Patient_Status == 1 )~Protein4+Tumour_Stage,data=test,x=1)
xs=Score(list(Cox1=fs1),Surv(time = days_of_last_follow_up ,event = Patient_Status == 1 )~1,data=test,
         plots="cal",metrics=NULL)
plotCalibration(xs, xlim = c(0, 0.4))

# IVYKIS 2 ----------------------------------------------------------------
# modeliu sukurimai -------------------------------------------------------
cox_univ_amzius<-coxph(Surv(time = days_of_last_follow_up ,event = Patient_Status == 2)~ Age,
                       data =  train) 

cox_univ_protein1<-coxph(Surv(time = days_of_last_follow_up ,event = Patient_Status == 2)~ Protein1,
                         data =  train)  

cox_univ_protein2<-coxph(Surv(time = days_of_last_follow_up ,event = Patient_Status == 2)~ Protein2,
                         data =  train)  

cox_univ_protein3<-coxph(Surv(time = days_of_last_follow_up ,event = Patient_Status == 2)~ Protein3,
                         data =  train)  

cox_univ_protein4<-coxph(Surv(time = days_of_last_follow_up ,event = Patient_Status == 2)~ Protein4,
                         data =  train)  

cox_univ_tumour<-coxph(Surv(time = days_of_last_follow_up ,event = Patient_Status == 2)~ Tumour_Stage,
                       data =  train) 

cox_univ_histology<-coxph(Surv(time = days_of_last_follow_up ,event = Patient_Status == 2)~ Histology,
                          data =  train)   

cox_univ_her<-coxph(Surv(time = days_of_last_follow_up ,event = Patient_Status == 2)~ HER2.status,
                    data =  train)  

cox_univ_sur<-coxph(Surv(time = days_of_last_follow_up ,event = Patient_Status == 2)~ Surgery_type,
                    data =  train)  

# PRIELAIDU PATIKRINIMAS --------------------------------------------------
# isskirtys ---------------------------------------------------------------
hold<-2/sqrt(nrow(train))
hold

dfbeta_amzius<-residuals(cox_univ_amzius, type="dfbeta")
dfbeta_protein1<-residuals(cox_univ_protein1, type="dfbeta")
dfbeta_protein2<-residuals(cox_univ_protein2, type="dfbeta")
dfbeta_protein3<-residuals(cox_univ_protein3, type="dfbeta")
dfbeta_protein4<-residuals(cox_univ_protein4, type="dfbeta")
dfbeta_tumour<-residuals(cox_univ_tumour, type="dfbeta")
dfbeta_histology<-residuals(cox_univ_histology, type="dfbeta")
dfbeta_her<-residuals(cox_univ_her, type="dfbeta")
dfbeta_sur<-residuals(cox_univ_sur, type="dfbeta")

vardai<-c("Amžius", "I proteino lygis","II proteino lygis", 
          "III proteino lygis", "IV proteino lygis", 
          "Stadija","Histologija",
          "HER2 statusas","Operacija")
dfbeta<-data.frame(dfbeta_amzius, dfbeta_protein1, dfbeta_protein2, dfbeta_protein3, dfbeta_protein4,
                   dfbeta_tumour, dfbeta_histology, dfbeta_her,  dfbeta_sur)

par(mfrow=c(2, 3))
for (j in 7:9) { #1:4 - amzius, protein1, protein2, protein3; 3:8 - protein4
  plot(dfbeta[, j],
       ylab=vardai[j], xlab="")
  title(xlab="Indeksas")
  abline(h=0, lty=2)
  abline(h = hold, lty = 2, col="red")
  abline(h = -hold, lty=2, col="red")
}#yra isskirciu

#kurie indeksai yra iskirtis
which(abs(dfbeta[, 1]) > hold)
which(abs(dfbeta[, 2]) > hold)
which(abs(dfbeta[, 3]) > hold)
which(abs(dfbeta[, 4]) > hold)
which(abs(dfbeta[, 5]) > hold)


#pasalinti tas isskirtis
train <- train[-c(5, 8, 241, 248, 110, 17, 41, 55, 140, 169), ]


# multikolinerumas --------------------------------------------------------
cox_kiek <- coxph(Surv(time = days_of_last_follow_up ,event = Patient_Status == 2)~ Protein1+Protein2+Protein3+Protein4,
                  data =  train)

vif(cox_kiek)

# kov. reiksmingumas ir ph prielaida --------------------------------------
#reiksmingumo lygmuo 0.15

cox_univ <- coxph(Surv(time = days_of_last_follow_up ,event = Patient_Status == 2)~ Age+Protein1+Protein2+Protein3+Protein4+
                    Tumour_Stage+Histology+HER2.status+Surgery_type,
                  data =  train)

summary(cox_univ) 
cox.zph(cox_univ) 

stepAIC(cox_univ, method="both")

# modelis -----------------------------------------------------------------
cox<-coxph(Surv(time = days_of_last_follow_up ,event = Patient_Status == 2 )~ Surgery_type,
           data =  train)
summary(cox)
#konkordacija 0.619
AIC(cox)
#AIC - 414.0342


# kalibravimas ------------------------------------------------------------
fs2=coxph(Surv(time = days_of_last_follow_up ,event = Patient_Status == 2 )~Surgery_type,data=train,x=1)
xs=Score(list(Cox1=fs2),Surv(time = days_of_last_follow_up ,event = Patient_Status == 2 )~1,data=train,
         plots="cal",metrics=NULL)
plotCalibration(xs, xlim = c(0, 0.4))
#plotCalibration(xs,cens.method="local",pseudo=1)
#plotCalibration(xs, method = "quantile", cens.method = "jackknife", breaks = seq(0, 1, 0.1))

fs2=coxph(Surv(time = days_of_last_follow_up ,event = Patient_Status == 2 )~Surgery_type,data=test,x=1)
xs=Score(list(Cox1=fs2),Surv(time = days_of_last_follow_up ,event = Patient_Status == 2 )~1,data=test,
         plots="cal",metrics=NULL)
plotCalibration(xs, xlim = c(0, 0.4))

# cox nomograma ---------------------------------------------------------------
library(QHScrnomo)
cox2<-coxph(Surv(time = days_of_last_follow_up ,event = Patient_Status == 2 )~  HER2.status+Surgery_type,
            data =  train)

cox1<-coxph(Surv(time = days_of_last_follow_up ,event = Patient_Status == 1 )~  Protein4+Tumour_Stage,
            data =  train)


dd <- datadist(train)
options(datadist = "dd")

cox2.f <- cph(Surv(time = days_of_last_follow_up ,event = Patient_Status == 2 )~  Surgery_type,
              data =  train,x = TRUE, y= TRUE, surv=TRUE)
cox2.crr <- crr.fit(cox2.f,cencode = 0,failcode = 2)
nomogram.crr(cox2.crr,failtime = 1080,lp=FALSE,
             funlabel = "Predicted 3-year risk for non-cancer death")


cox1.f <- cph(Surv(time = days_of_last_follow_up ,event = Patient_Status == 1 )~  Protein4+Tumour_Stage,
              data =  train,x = TRUE, y= TRUE, surv=TRUE)
cox1.crr <- crr.fit(cox1.f,cencode = 0,failcode = 1)
nomogram.crr(cox1.crr,failtime = 1080,lp=FALSE,
             funlabel = "Predicted 3-year risk for cancer death")

