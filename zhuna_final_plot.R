
 library(verification)
 library(pROC)    # ROC
 library(rms)     # nomogram
 library(glmnet)  # cv.glmnet


 .
 #--------------- Environment settings -----------------
 rm(list=ls())  # 清空
 setwd("C:/Users/212481425/Desktop/0114_zhuna/6_for_replot_gy")
 getwd()

#---------------------------------- Nomogram ------------------------------------------
 df_train<- read.csv("For_nomogram.csv",header = T)
 df_test <- read.csv("For_nomogram.csv",header = T)
 response_train<-as.factor(df_train$label) 
 response_test<-as.factor(df_test$label) 
 names(df_train)
 train<-df_train[-1]
 test<-df_test[-1]
 names(train)

train_for_nomogram<-train
names(train_for_nomogram)


tiff(filename = "C:/Users/212481425/Desktop/0114_zhuna/8_figures/nomogram_600dpi.tiff",width=5000,height = 3000,res=600)

# nomogram_train
{
   dt <- data.frame(train_for_nomogram)
   ddist    <- datadist(dt)
   options(datadist="ddist")
   flrm     <- lrm(response_train ~ .,data=dt,x=T,y=T)
   
   # va <- validate(flrm,method="boot",B=150,dxy=T,pr=T)
   # cal  <- calibrate(flrm,method="boot",B=150)
   nomogplot <- nomogram(flrm,
                         fun=plogis,
                         fun.at=c(.1, .5,.9),
                         lp=F,
                         funlabel="Risk of Malignance")
   par(mfrow=c(1,1))
   plot(nomogplot)
}
dev.off()

#---------------------------------------------  ROC  -----------------------------------------------------


df_1_train<- read.csv("Conventional_SCORE_seed_5_train.csv",header = T)
df_1_test<- read.csv("Conventional_SCORE_seed_5_test.csv",header = T)

df_2_train<- read.csv("Radiomics_SCORE_Seed_70_train.csv",header = T)
df_2_test<- read.csv("Radiomics_SCORE_Seed_70_test.csv",header = T)

df_3_train<- read.csv("Integrated_score_seed_56_train.csv",header = T)
df_3_test<- read.csv("Integrated_score_seed_56_test.csv",header = T)


response_train_1<-as.factor(df_1_train$label)
response_test_1<-as.factor(df_1_test$label)

response_train_2<-as.factor(df_2_train$label)
response_test_2<-as.factor(df_2_test$label)

response_train_3<-as.factor(df_3_train$label)
response_test_3<-as.factor(df_3_test$label)



#------------------- Delong -------------------------------------

roc_train_conv<-roc(response_train_1,df_1_train$Conv_Score,ci=T)
roc_train_rad<-roc(response_train_2,df_2_train$Rad_Score,ci=T)
roc_train_int<-roc(response_train_3,df_3_train$Int_Score,ci=T)

roc_test_conv<-roc(response_test_1,df_1_test$Conv_Score,ci=T)
roc_test_rad<-roc(response_test_2,df_2_test$Rad_Score,ci=T)
roc_test_int<-roc(response_test_3,df_3_test$Int_Score,ci=T)


roc.test(roc_train_conv,roc_train_rad)$p.value
roc.test(roc_train_conv,roc_train_int)$p.value
roc.test(roc_train_int,roc_train_rad)$p.value

roc.test(roc_test_conv,roc_test_rad)$p.value
roc.test(roc_test_conv,roc_test_int)$p.value
roc.test(roc_test_int,roc_test_rad)$p.value


roc.test(roc_train_conv,roc_test_conv)$p.value
roc.test(roc_train_int,roc_test_int)$p.value
roc.test(roc_train_rad,roc_test_rad)$p.value

#-------------------------- ROC Curve ------------------------------------


tiff(file = "C:/Users/212481425/Desktop/0114_zhuna/8_figures/ROC_train_600.tiff", res = 600,width = 4500, height = 4500)


plot(roc_train_conv, col=1, main="The ROC curves in the training dataset")
plot(roc_train_rad, col=4,add=T)
plot(roc_train_int, col=2,add=T)

legend(0.85,0.15,
       c("Conventional Score: AUC=0.882,95% CI=0.823-0.941",
         "Radiomics Score:     AUC=0.961,95% CI=0.925-0.996",
         "Integrated Score:      AUC=0.969,95% CI=0.940-0.997"
         
       ),
       
       border=3,
       cex=1.1,
       text.width = 0.75,
       col=c(1,4,2),
       lty= 1,
       lwd= 3)

dev.off() 

tiff(file = "C:/Users/212481425/Desktop/0114_zhuna/8_figures/ROC_test_600.tiff", res = 600,width = 4500, height = 4500)

plot(roc_test_conv, col=1, main="The ROC curves in the testing dataset")
plot(roc_test_rad, col=4,add=T)
plot(roc_test_int, col=2,add=T)

legend(0.85,0.15,
       c("Conventional Score: AUC=0.857,95% CI=0.767-0.948",
         "Radiomics Score:     AUC=0.955,95% CI=0.914-0.996",
         "Integrated Score:      AUC=0.966,95% CI=0.928-1.000"
         
       ),
       
       border=3,
       cex=1.1,
       text.width = 0.75,
       col=c(1,4,2),
       lty= 1,
       lwd= 3)

dev.off() 


plot(roc_train_conv,print.thres=T,print.auc=T)
plot(roc_test_conv,print.thres=0.597,print.auc=T)
plot(roc_train_rad,print.thres=T,print.auc=T)
plot(roc_test_rad,print.thres=0.905,print.auc=T)
plot(roc_train_int,print.thres=1.014,print.auc=T)
plot(roc_test_int,print.thres=1.014,print.auc=T)
#--------------------------- Correlation map ------------------------------------------
library(corrplot)
setwd("C:/Users/212481425/Desktop")
mt<- read.csv("radiomics_90_60_icc_0.8_left_train.csv")

cor_matr<- cor(mt)

tiff(file = "C:/Users/212481425/Desktop/correlation_41_train_600.tiff", res = 600,width = 3500, height = 3500)

corrplot(cor_matr, tl.col="black", tl.srt=45)

dev.off() 


mt_2<- read.csv("for_correlation_map_102_train.csv")

cor_matr_2<- cor(mt_2)

tiff(file = "C:/Users/212481425/Desktop/correlation_102_train_600.tiff", res = 600,width = 7000, height = 7000)
corrplot(cor_matr_2, tl.col="black", tl.srt=45)
dev.off() 














#______________________________ Calibration Curve _________________________________


library(ModelGood)

flrm_conv_train<-lrm(label~.,data=df_1_train,x=T,y=T )
flrm_conv_test<-lrm(label~.,data=df_1_test,x=T,y=T )

flrm_rad_train<-lrm(label~.,data=df_2_train,x=T,y=T )
flrm_rad_test<-lrm(label~.,data=df_2_test,x=T,y=T )

flrm_int_train<-lrm(label~.,data=df_3_train,x=T,y=T )
flrm_int_test<-lrm(label~.,data=df_3_test,x=T,y=T )

?calPlot2

library(DescTools)
HosmerLemeshowTest(cal_1_train$Frame$lrm,cal_1_train$Frame$jack)
HosmerLemeshowTest(cal_2_train$Frame$lrm,cal_2_train$Frame$jack)
HosmerLemeshowTest(cal_3_train$Frame$lrm,cal_3_train$Frame$jack)



tiff(filename = "C:/Users/212481425/Desktop/0114_zhuna/8_figures/Calibration_train_600dpi.tiff",width=4500,height = 4500,res=600)   


cal_1_train<-calPlot2( list(flrm_conv_train),col=1,lty=1,legend = F,showY =F)
cal_2_train<-calPlot2( list(flrm_rad_train),col=4,lty=1,showY =F,legend = F,add=T)
cal_3_train<-calPlot2( list(flrm_int_train),col=2,lty=1,showY =F,legend = F,add=T)


         
 
legend(0.35,0.2,
       c("Conventional Score: P value=0.364",
         "Radiomics Score:     P value=0.521 ",
         "Integrated Score:      P value=0.564"
         
       ),
       
       border=3,
       cex=1.1,
       text.width = 0.5,
       col=c(1,4,2),
       lty= 1,
       lwd= 3)
dev.off()
       
tiff(filename = "C:/Users/212481425/Desktop/0114_zhuna/8_figures/Calibration_test_600dpi.tiff",width=4500,height = 4500,res=600)   

HosmerLemeshowTest(cal_1_test$Frame$lrm,cal_1_test$Frame$jack)
HosmerLemeshowTest(cal_2_test$Frame$lrm,cal_2_test$Frame$jack)
HosmerLemeshowTest(cal_3_test$Frame$lrm,cal_3_test$Frame$jack)



cal_1_test<-calPlot2( list(flrm_conv_test),col=1,lty=1,legend = F,showY =F)
cal_2_test<-calPlot2( list(flrm_rad_test),col=4,lty=1,showY =F,legend = F,add=T)
cal_3_test<-calPlot2( list(flrm_int_test),col=2,lty=1,showY =F,legend = F,add=T)




legend(0.35,0.2,
       c("Conventional Score: P value=0.390",
         "Radiomics Score:     P value=0.838 ",
         "Integrated Score:      P value=0.859"
         
       ),
       
       border=3,
       cex=1.1,
       text.width = 0.5,
       col=c(1,4,2),
       lty= 1,
       lwd= 3)
dev.off()

#________________________________ DCA _______________________________________________
library(rmda) 
names(df_2_train)

dca_conv_train <- decision_curve(label~Conv_Score,data=df_1_train)
dca_rad_train <- decision_curve(label~Rad_Score,data=df_2_train)
dca_int_train<-decision_curve(label~Int_Score, data=df_3_train)

dca_conv_test <- decision_curve(label~Conv_Score,data=df_1_test)
dca_rad_test <- decision_curve(label~Rad_Score,data=df_2_test)
dca_int_test<-decision_curve(label~Int_Score, data=df_3_test)


tiff(filename = "C:/Users/212481425/Desktop/0114_zhuna/8_figures/DCA_train_600dpi.tiff",width=4500,height = 3500,res=600)

plot_decision_curve(list( dca_conv_train, dca_rad_train,dca_int_train), 
                    confidence.intervals = F, 
                    col = c(1,4,2), 
                    curve.names = c('Conventional Score','Radiomics Score','Integrated Score'),
                    lty=c(1,1,1),
                    lwd=3,
                    xlab="Threshold Probability",
                    ylab="Net Benefit",
                    legend.position = 'left', 
                    cost.benefits = T
)
dev.off()


?plot_decision_curve

tiff(filename = "C:/Users/212481425/Desktop/0114_zhuna/8_figures/DCA_test_600dpi.tiff",width=4500,height = 3500,res=600)

plot_decision_curve(list( dca_conv_test, dca_rad_test,dca_int_test), 
                    confidence.intervals = F, 
                    col = c(1,4,2), 
                    curve.names = c('Conventional Score','Radiomics Score','Integrated Score'),
                    lty=c(1,1,1),
                    lwd=3,
                    xlab="Threshold Probability",
                    ylab="Net Benefit",
                    legend.position = 'left', 
                    cost.benefits = T
)
dev.off()










#---------------------------------- vio plot ----------------------
names(df_1_train)

library(ggplot2)
library(ggpubr)
tiff(filename = "C:/Users/212481425/Desktop/0114_zhuna/8_figures/Conv_Score_vioplot_train_600dpi.tiff",width=3500,height = 3500,res=600)

vioplot::vioplot(df_1_train$Conv_Score~df_1_train$label,
                 main="Distribution of Conventional Score in Training dataset ",
                 col=c("deepskyblue3","orange2"),
                 ylim=c(-5,20),
                 xlab="Label",ylab="Conventional Score")
abline(h = 0.597, v = 0, col = 1,lty=2)
legend(1,6, "Cutoff=0.597",bty="n")

dev.off()
tiff(filename = "C:/Users/212481425/Desktop/0114_zhuna/8_figures/Conv_Score_vioplot_test_600dpi.tiff",width=3500,height = 3500,res=600)
vioplot::vioplot(df_1_test$Conv_Score~df_1_test$label,
                 main="Distribution of Conventional Score in Testing dataset ",
                 col=c("deepskyblue3","orange2"),ylim=c(-5,20),
                 xlab="Label",ylab="Conventional Score")
abline(h = 0.597, v = 0, col = 1,lty=2)
legend(1,10, "Cutoff=0.597",bty="n")
dev.off()








tiff(filename = "C:/Users/212481425/Desktop/0114_zhuna/8_figures/Rad_Score_vioplot_train_600dpi.tiff",width=3500,height = 3500,res=600)

vioplot::vioplot(df_2_train$Rad_Score~df_2_train$label,
                 main="Distribution of Radiomics Score in Training dataset ",
                 col=c("deepskyblue3","orange2"),
                 ylim=c(-22,12),
                 xlab="Label",ylab="Radiomics Score")
abline(h = 0.905, v = 0, col = 1,lty=2)
legend(1,8, "Cutoff=0.905",bty="n")

dev.off()
tiff(filename = "C:/Users/212481425/Desktop/0114_zhuna/8_figures/Rad_Score_vioplot_test_600dpi.tiff",width=3500,height = 3500,res=600)
vioplot::vioplot(df_2_test$Rad_Score~df_2_test$label,
                 main="Distribution of Radiomics Score in Testing dataset ",
                 col=c("deepskyblue3","orange2"),
                 ylim=c(-22,12),
                 xlab="Label",ylab="Radiomics Score")
abline(h = 0.905, v = 0, col = 1,lty=2)
legend(1,5, "Cutoff=0.905",bty="n")
dev.off()






tiff(filename = "C:/Users/212481425/Desktop/0114_zhuna/8_figures/Int_Score_vioplot_train_600dpi.tiff",width=3500,height = 3500,res=600)

vioplot::vioplot(df_3_train$Int_Score~df_3_train$label,
                 main="Distribution of Integrated Score in Training dataset ",
                 col=c("deepskyblue3","orange2"),
                 ylim=c(-25,12),
                 xlab="Label",ylab="Integrated Score")
abline(h = 1.014, v = 0, col = 1,lty=2)
legend(1,8, "Cutoff=1.014",bty="n")

dev.off()
tiff(filename = "C:/Users/212481425/Desktop/0114_zhuna/8_figures/Int_Score_vioplot_test_600dpi.tiff",width=3500,height = 3500,res=600)
vioplot::vioplot(df_3_test$Int_Score~df_3_test$label,
                 main="Distribution of Integrated Score in Testing dataset ",
                 col=c("deepskyblue3","orange2"),
                 ylim=c(-25,12),
                 xlab="Label",ylab="Integrated Score")
abline(h = 1.014, v = 0, col = 1,lty=2)
legend(1,7, "Cutoff=1.014",bty="n")
dev.off()


#------------------ delong --------------------



setwd("C:/Users/212481425/Desktop/0114_zhuna/5_reanalysis")

df_1<- read.csv("A_65_SCORE_Seed_10_train.csv",header = T)
df_2 <- read.csv("V_65_SCORE_Seed_55_train.csv",header = T)

df_3<- read.csv("A_65_SCORE_Seed_10_test.csv",header = T)
df_4 <- read.csv("V_65_SCORE_Seed_55_test.csv",header = T)


roc_1<-roc(df_1$label,df_1$A65_Score,ci=T)
roc_2<-roc(df_2$label,df_2$V65_Score,ci=T)
roc_3<-roc(df_3$label,df_3$A65_Score,ci=T)
roc_4<-roc(df_4$label,df_4$V65_Score,ci=T)

roc.test(roc_1,roc_2)
roc.test(roc_3,roc_4)



