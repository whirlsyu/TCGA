args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).\n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  dis_name = args[1]
}
cat('######################',dis_name,'#######################','\n',sep = " ")
source('function.R')
require(survival)
require(EnMCB)
library(readxl)
library(data.table)
library(readr)

mmc1 <- readxl::read_excel("mmc1.xlsx", sheet = "TCGA-CDR")

library(RTCGA)
downloadTCGA(cancerTypes = dis_name, '.Merge_methylation__humanmethylation450__jhu_usc_edu__Level_3', destDir ='.', date = "2016-01-28")

#This url may need further fixed when using different OS.

url_tcga= paste( "./", 
                 dis_name, "/gdac.broadinstitute.org_", 
                 dis_name,".Merge_methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.Level_3.2016012800.0.0/",
                 dis_name,".methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt", sep="")

cat('pre process the data for',dis_name,'\n')

#prepare the methylation matrix
cat('downloading the files for ',dis_name,'\n')
met_prepare<-as.matrix(fread(input = url_tcga
                             ,sep = "\t",header = T))

cat('Checking the methylation files. If error occurs, please change the script manually, since the format of files may further changed.',cancer_name,'\n')
rownames(met_prepare)<-met_prepare[,1]
met_prepare<-met_prepare[,grep("Beta_value",met_prepare[1,])]

met_prepare<-met_prepare[,grep("-01A-",colnames(met_prepare))]
pdata<-data.frame(mmc1[mmc1$type %in% dis_name,])
rownames(pdata)<-pdata$bcr_patient_barcode

remaining<-intersect(substr(colnames(met_prepare),1,12),pdata$bcr_patient_barcode)

colnames(met_prepare)<-substr(colnames(met_prepare),1,12)

met_prepare<-met_prepare[-1,remaining]
pdata<-pdata[remaining,]

met_final<-pre_process_methylation(as.numeric_matrix(met_prepare),Mvalue=FALSE)


dir.create(dis_name)
setwd(dis_name)
MCB_INFO<-mcb$MCBinformation
pdata<-data.frame(mmc1[mmc1$type %in% dis_name,])
rownames(pdata)<-pdata$bcr_patient_barcode
remaining<-intersect(colnames(met_final),pdata$bcr_patient_barcode)
pdata<-pdata[remaining,]
survival_data<-Surv(time = pdata$PFI.time/365.25,event = pdata$PFI)
set.seed(1)
tr_rz<-colnames(met_final) %in% sample(colnames(met_final),0.7*length(colnames(met_final)))
#tr_rz<-colnames(met_final) %in% colnames(models$MCB_cox_matrix_training) 
te_rz<-!tr_rz
trainingset<-met_final[,tr_rz]
testingset <-met_final[,te_rz]
y_surv_train <- survival_data[tr_rz]
y_surv_test <- survival_data[te_rz]
total_res_select<-MCB_INFO[as.numeric(MCB_INFO[,'CpGs_num'])>=5,]
#feature selection based on the mean value of CpG sites in MCB using L1 penalty 
res_mean_all<-apply(total_res_select,1, function(x,data_set=trainingset){
   CpGs<-strsplit(x['CpGs']," ")[[1]]
   colMeans(data_set[CpGs,])
 })
rz=!(is.na(y_surv_train)|y_surv_train[,1]==0)
library(glmnet)

cvg<-cv.glmnet(res_mean_all[rz,],y_surv_train[rz],family='cox',type.measure = 'C')
plot(cvg)


lambda <- cvg$lambda.1se
repeat{
  #at least 10 MCBs are included
  if (sum(coef(cvg, s=lambda)!=0) >=10 & sum(coef(cvg, s=lambda)!=0) <=20) break
  if (sum(coef(cvg, s=lambda)!=0) < 10) lambda = lambda - 0.0001
  if (sum(coef(cvg, s=lambda)!=0) > 20) lambda = lambda + 0.0001
}
  
total_res_select_filtered<- total_res_select[which(coef(cvg, s=lambda)!=0),]

# #S.t.x method results are not as good as cv.glmnet, therefore here we only use cv.glmnet
# stx_selection<-apply(res_mean_all[rz,1:10], 2, 
#       function(x,...) {
#         s_t_x(Stime=y_surv_train[rz,1],
#               status = y_surv_train[rz,2],
#                 marker = x, predict.time = 5)$ks_res$p.value
#       })
# stx_selection<-unlist(stx_selection)
# stx_selection<-p.adjust(stx_selection)
# #### end selections ####

new_res_cv<-NULL
for (i in seq(nrow(total_res_select_filtered))) {
  single_in_mcb_new_res_select_filtered<-total_res_select_filtered[i,]
  try(new_res_cv<-rbind(new_res_cv,auc_cal_cv_tt(single_in_mcb_new_res_select_filtered,
                                                 trainingset,
                                                 testingset,
                                                 y_surv_train,
                                                 y_surv_test,
                                                 nfold = 10,seed = 6)))
}
write.csv(new_res_cv,file = paste('results',dis_name,'.csv',sep = ""))

RP_data<-new_res_cv[,c('COX_AUC_CV.AUC','SVM_AUC_CV.AUC','eNet_AUC_CV.AUC','mboost_AUC_CV.AUC','EnMCB_AUC_CV.AUC')]
RP_data<-as.numeric_matrix(RP_data)

#boxplot
large_data<-data.frame(rbind(
  cbind(RP_data[,'COX_AUC_CV.AUC'],rep('Cox',nrow(RP_data))),
  cbind(RP_data[,'SVM_AUC_CV.AUC'],rep('eNet',nrow(RP_data))),
  cbind(RP_data[,'eNet_AUC_CV.AUC'],rep('SVR',nrow(RP_data))),
  cbind(RP_data[,'mboost_AUC_CV.AUC'],rep('mboost',nrow(RP_data))),
  cbind(RP_data[,'EnMCB_AUC_CV.AUC'],rep('Ensemble',nrow(RP_data)))
))
colnames(large_data)<-c("value", "group")
large_data$value<-as.numeric(as.character(large_data$value))
boxplot_large_data(large_data)

#calculate the rank
RP_data_cal<-matrix(data=apply(RP_data, 2, function(x){(length(x)-rank(x,ties.method='first')+1)}),nrow = nrow(RP_data),ncol = ncol(RP_data))
RP_data_cal<-data.frame(RP_data_cal,RP=apply(RP_data_cal, 1, function(x)sum(log(x))))
rownames(RP_data_cal)<-new_res_cv[,'MCB_no']
res_cv_RP<-data.frame(total_res_select_filtered[total_res_select_filtered[,'MCB_no'] %in% rownames(RP_data_cal),],new_res_cv,RP=RP_data_cal$RP)
write.csv(res_cv_RP,file = paste('res_cv_RP_',dis_name,'.csv',sep = ""))

#get the best one
MCBblocks_selected=new_res_cv[which.min(res_cv_RP$RP),'MCB_no']

single_res<-t(as.matrix(total_res_select_filtered[total_res_select_filtered[,'MCB_no'] %in% MCBblocks_selected,]))

em<-ensemble_model(single_res = single_res,training_set = trainingset,Surv_training = y_surv_train)

#use the model to predict the responses in training set
ep<-ensemble_prediction(ensemble_model = em,prediction_data = trainingset,mutiple_results = T)
#ep<-ensemble_prediction_lp(ensemble_model = em,prediction_data = trainingset)
#calculation of AUC, plot and save the results in pdf files
ROC_multiple(test_frame = data.frame(cox=ep['cox',],
                                     svm=ep['svm',],
                                     eNet=ep['enet',],
                                     mboost=ep['mboost',],
                                     ensemble=ep['ensemble',]
),
y = y_surv_train,file_name = paste(dis_name,"multiple model training")
)

#use the model to predict the response in testing set
ep_t<-ensemble_prediction(ensemble_model = em,prediction_data = testingset,mutiple_results = T)
#ep_t<-ensemble_prediction_lp(ensemble_model = em,prediction_data = testingset)
#calculation of AUC, plot and save the results in pdf files
ROC_multiple(test_frame = data.frame(cox=ep_t['cox',],
                                     svm=ep_t['svm',],
                                     eNet=ep_t['enet',],
                                     mboost=ep_t['mboost',],
                                     ensemble=ep_t['ensemble',]
),
y = y_surv_test,file_name = paste(dis_name,"multiple model test") 
)

#draw survival curve for training set
EnMCB::draw_survival_curve(ep['ensemble',],
                           living_days = y_surv_train[,1],
                           living_events = y_surv_train[,2],
                           write_name = paste(dis_name,"ensemble model training"),file = T)

#draw survival curve for testing set
EnMCB::draw_survival_curve(ep_t['ensemble',],
                           living_days = y_surv_test[,1],
                           living_events = y_surv_test[,2],
                           write_name = paste(dis_name,"ensemble model testing"),threshold = median(ep['ensemble',]),file = T)

cat('######################',dis_name,' finished! #######################','\n',sep = " ")
