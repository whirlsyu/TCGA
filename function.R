# tested environment:
# those functions are compatible with 3.6.0 or higher

require('survivalROC')
require('prognosticROC')

as.numeric_matrix<-function(mat){
  matf<-matrix(as.numeric(mat),
               nrow = NROW(mat),
               ncol = NCOL(mat))
  rownames(matf)<-rownames(mat)
  colnames(matf)<-colnames(mat)
  matf
}

auc_cal_cv_tt <- function(single_one,all_data_train,data_test,y_surv_train,y_surv_test,nfold=10,predict_time=3,seed=NULL) {
  
  require(EnMCB)
  #remove all the NA and zero survival data
  na_or_zero_data<-(is.na(y_surv_train[,1])|y_surv_train[,1]==0)
  all_data_train = all_data_train[,!na_or_zero_data]
  y_surv_train = y_surv_train[!na_or_zero_data]
  
  na_or_zero_data<-(is.na(y_surv_test[,1])|y_surv_test[,1]==0)
  data_test = data_test[,!na_or_zero_data]
  y_surv_test = y_surv_test[!na_or_zero_data]
  
  set.seed(seed)
  sp<-sample(1:length(y_surv_train),replace = F)
  order_sp<-order(sp)
  y_surv_train = y_surv_train[sp]
  all_data_train = all_data_train[,sp]
  
  train<-1:length(y_surv_train)
  folds <- cut(seq(1,length(train)),breaks=nfold,labels=FALSE)
  survival_predictions<-NULL
  for (i in seq(unique(folds))) {
    rz<- which(folds==i,arr.ind=TRUE)
    em<-ensemble_model(t(single_one),all_data_train[,train[-rz]],y_surv_train[-rz])
    survival_prediction<-ensemble_prediction(em,all_data_train[,train[rz]],mutiple_results = T)
    survival_predictions<-cbind(survival_predictions,survival_prediction)
  }
  
  em_all<-ensemble_model(t(single_one),all_data_train,y_surv_train)
  pre_test<-ensemble_prediction(em_all,data_test,mutiple_results = T)
  
  auc_COX_cv<-calculate_auc_ci(survival = y_surv_train,survival_predictions['cox',])
  auc_SVR_cv<-calculate_auc_ci(survival = y_surv_train,survival_predictions['svm',])
  auc_eNet_cv<-calculate_auc_ci(survival = y_surv_train,survival_predictions['enet',])
  auc_mboost_cv<-calculate_auc_ci(survival = y_surv_train,survival_predictions['mboost',])
  auc_em_cv<-calculate_auc_ci(survival = y_surv_train,survival_predictions['ensemble',])
  
  auc_COX_test<-calculate_auc_ci(survival = y_surv_test,pre_test['cox',],ref = survival_predictions['cox',])
  auc_SVR_test<-calculate_auc_ci(survival = y_surv_test,pre_test['svm',],ref = survival_predictions['svm',])
  auc_eNet_test<-calculate_auc_ci(survival = y_surv_test,pre_test['enet',],ref = survival_predictions['enet',])
  auc_mboost_test<-calculate_auc_ci(survival = y_surv_test,pre_test['mboost',],ref = survival_predictions['mboost',])
  auc_em_test<-calculate_auc_ci(survival = y_surv_test,pre_test['ensemble',],ref = survival_predictions['ensemble',])
  
  return(c(single_one,
           COX_AUC_CV = auc_COX_cv['AUC'],COX_95CI_CV = paste("(",auc_COX_cv['CI95'],")",sep = " "),COX_brier_score_CV = auc_COX_cv['brier_score'],
           SVM_AUC_CV = auc_SVR_cv['AUC'], SVM_95CI_CV = paste("(",auc_SVR_cv['CI95'],")",sep = " "),SVM_brier_score_CV = auc_SVR_cv['brier_score'],
           eNet_AUC_CV = auc_eNet_cv['AUC'],eNet_95CI_CV = paste("(",auc_eNet_cv['CI95'],")",sep = " "),eNet_brier_score_CV = auc_eNet_cv['brier_score'],
           mboost_AUC_CV = auc_mboost_cv['AUC'],mboost_95CI_CV= paste("(",auc_mboost_cv['CI95'],")",sep = " "),mboost_brier_score_CV = auc_mboost_cv['brier_score'],
           EnMCB_AUC_CV = auc_em_cv['AUC'],EnMCB_95CI_CV = paste("(",auc_em_cv['CI95'],")",sep = " "),EnMCB_brier_score_CV = auc_em_cv['brier_score'],
           
           COX_AUC_test =auc_COX_test['AUC'], COX_95CI_test = paste("(",auc_COX_test['CI95'],")",sep = " "),COX_brier_score_test = auc_COX_test['brier_score'],
           SVM_AUC_test =auc_SVR_test['AUC'], SVM_95CI_test = paste("(",auc_SVR_test['CI95'],")",sep = " "),SVM_brier_score_test = auc_SVR_test['brier_score'],
           eNet_AUC_test =auc_eNet_test['AUC'], eNet_95CI_test = paste("(",auc_eNet_test['CI95'],")",sep = " "),eNet_brier_score_test = auc_eNet_test['brier_score'],
           mboost_AUC_test =auc_mboost_test['AUC'], mboost_95CI_test= paste("(",auc_mboost_test['CI95'],")",sep = " "),mboost_brier_score_test = auc_mboost_test['brier_score'],
           EnMCB_AUC_test =auc_em_test['AUC'], EnMCB_95CI_test = paste("(",auc_em_test['CI95'],")",sep = " "),EnMCB_brier_score_test = auc_em_test['brier_score']))
}


ROC_multiple<-function(test_frame,y_surv,file_name="title",ntime=3){
  require("prognosticROC") #version 0.7
  require("ggplot2")
  require("plotROC")
  sroclong_all<-NULL
  na_or_0 <- is.na(y_surv)|y_surv[,1]==0
  y_surv = y_surv[!na_or_0]
  test_frame = test_frame[!na_or_0,]
  res_list<-list()
  for (n in 1:ncol(test_frame)) {
    ROC_res= survivalROC::survivalROC(Stime=y_surv[,1],
                                      status=y_surv[,2],
                                      marker =as.numeric(test_frame[,n]),
                                      predict.time = ntime,span =0.25*length(y_surv)^(-0.20),window = "asymmetric" )#
    set.seed(6)
    ci_res = boot::boot(data=data.frame(survival = y_surv[,1],
                                        survival_status = y_surv[,2],
                                        marker = as.numeric(test_frame[,n])),
                        statistic=bs_ci, R=1000, predict.time = ntime)
    res = boot::boot.ci(ci_res,type="perc")
    AUC_95CI<-paste(format(ROC_res$AUC,digits = 3),'(' ,format(res$percent[,4], digits = 3)," - ", format(res$percent[,5], digits = 3),')',sep = "")
    res_list$AUC<-cbind(res_list$AUC,ROC_res$AUC)
    res_list$AUC_95CI<-cbind(res_list$AUC_95C,AUC_95CI)
    sroclong_all<-ROCdata_save(sroclong_all,ROC_res,mark = paste(ntime,"year AUC at",colnames(test_frame)[n],AUC_95CI,collapse = " "))
  }
  colnames(res_list$AUC)<-colnames(test_frame)
  colnames(res_list$AUC_95CI)<-colnames(test_frame)
  # The palette with grey:
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  # The palette with black:
  # cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  # To use for fills, add
  #scale_fill_manual(values=cbPalette)
  # To use for line and point colors, add
  gg<-ggplot2::ggplot(sroclong_all, aes(x = FPF, y = TPF, label = c, color = group))+
    coord_cartesian(ylim=c(0, 1.05),xlim=c(0, 1.05),expand=FALSE)+
    geom_roc(labels = FALSE, stat = "identity") + style_roc(theme = theme_light())+
    theme(legend.position=c(1,0.365), legend.justification=c(1,1),legend.title = element_blank(),legend.text = element_text(size = 8))+
    theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),axis.line = element_line(colour = "black"))+
    theme(axis.text.x = element_text(color = "grey20", size = 10),axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          axis.text.y = element_text(color = "grey20", size = 10))+
    geom_abline(slope=1, colour="black")+scale_colour_manual(values=cbPalette)
  ggplot2::ggsave(filename = paste("Time ROC of ",file_name,".jpeg",sep=""),plot = gg,device ="jpeg" ,
                  path = getwd(),dpi = 300,units = "in",width = 5, height = 4.5,
                  limitsize=F)
  return(res_list)
}

ROCdata_save<-function(origin=NULL,perf,mark="none"){
  sroclong<-data.frame(TPF = perf$TP, FPF =perf$FP, 
                       c = perf$cut.values, 
                       group = rep(mark, length(perf$TP)) )
  sroclong<-rbind(origin,sroclong)
  sroclong
}

bs_ci <- function(data, indices, predict.time = 3) {
  d <- data[indices, ]
  surv.res = survivalROC::survivalROC(Stime = d$survival, status = d$survival_status, 
                                      marker = d$marker, predict.time = predict.time,window = "asymmetric" ,
                                      span = 0.25 * NROW(d)^(-0.2))
  return(surv.res$AUC)
}
calculate_auc_ci <- function(survival, marker, predict_time = 3,select_time = F,ci = T,threshold = NA, ref = NULL,seed=6){
  na_or_0 <- is.na(survival)|survival[,1]==0
  survival = survival[!na_or_0]
  marker = marker[!na_or_0]  
  if (ci){
    # The lower the Brier score is for a set of predictions, the better the predictions are calibrated. 
    # Note that the Brier score, in its most common formulation, takes on a value between zero and one, 
    # since this is the square of the largest possible difference between a predicted probability 
    # (which must be between zero and one) and the actual outcome (which can take on values of only 0 or 1). 
    # In the original (1950) formulation of the Brier score, the range is double, from zero to two.
    # 
    # The Brier score is appropriate for binary and categorical outcomes that can be structured as true 
    # or false, but it is inappropriate for ordinal variables which can take on three or more values.
    if (select_time){
      predict_time = as.numeric(which.min(unlist(lapply(X = seq(floor(max(y_surv_train[,1]))), 
                                                        FUN = function(survival,marker,X){ipred::sbrier(survival,marker,X)},
                                                        survival=y_surv_train,marker = trainingset[2,]))))
    }
    require(scales)
    if (!is.null(ref)){
      sd_ref = sd(ref)
      threshold = median(ref)
    } else{
      sd_ref = sd(marker)
      threshold = median(marker)
    } 
    marker_z = (marker - threshold)/sd_ref
    marker_bs = rescale(marker_z,to =c(1,0))
    brier_score = as.numeric(ipred::sbrier(survival,marker_bs,predict_time))
    #   cat(brier_score,'\n')
    if (brier_score>1){
      cat('some problems may occurs in the data, see\n',survival,'\n')
      cat(marker,'\n')
    }
    set.seed(seed)
    ci_res = boot::boot(data=data.frame(survival = survival[,1],
                                        survival_status = survival[,2],
                                        marker = marker),
                        statistic=bs_ci, R=1000, predict.time = predict_time)
    res = boot::boot.ci(ci_res,type="perc")
    AUC = survivalROC::survivalROC(Stime = survival[,1], 
                                   status = survival[,2], 
                                   marker = marker, 
                                   predict.time = predict_time,
                                   window = "asymmetric",
                                   span = 0.25*length(survival)^(-0.20))$AUC
    return(c( AUC = format(AUC,digits = 4),
              CI95 = paste(format(res$percent[,4], digits = 4),"-", format(res$percent[,5], digits = 4)),
              brier_score = brier_score, predict_time = predict_time))
  }else{
    res = survivalROC::survivalROC(Stime = survival[,1], 
                                   status = survival[,2], 
                                   marker = marker, 
                                   predict.time = predict_time,
                                   window = "asymmetric", 
                                   span = 0.25*length(survival)^(-0.20))
    return(res)
  }
}

boxplot_large_data<-function(large_data){
  require(ggpubr)
  require(ggsignif)
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  GGM<-ggboxplot(large_data, #数据对象
                 x = 'group',# 选择x轴用那一列数据
                 y = 'value', #选择y轴用什么数据
                 fill = 'group', #颜色根据哪一列决定
                 bxp.errorbar = T, #是否添加 error bar
                 bxp.errorbar.width = 0.2, #error bar的长度
                 #palette = 'npg' #颜色风格
                 #add = 'point' # 是否添加boxplot上面的点点
  ) + labs(title = 'AUCs for different models in top ranked MCBs using cross validation', # 添加主标题
           subtitle = 'Plot of value by group', # 添加次标记
           caption = 'Data source: AUC', #添加脚注
           x = 'Groups (different models)', # x轴的名字
           y = 'AUCs' # y轴的名字
  ) + 
    # scale_y_continuous(limits = c(0, max(large_data$value)+0.05))+
    scale_fill_manual(values=cbPalette)+
    # geom_signif(comparisons = list(c('CpG Beta value', 'CpG M value'), c('MCB Beta value','MCB M value'), c('CpG Beta value','MCB Beta value'),c('CpG M value','MCB M value')), # 设置要对比的组
    #             y_position = c(34,36,38), #设置3个显著性标记的高度
    #             tip_length = c(0), #设置显著性那条横线两头向下的长度
    #             map_signif_level = T, #设置是否标记显著性的*号，还是直接标记数值
    #             test = t.test #设置显著性计算方式
    #             ) +
    # min.mean.sd.max <- function(x) {
    #   r <- c(min(x), mean(x) - sd(x), mean(x), mean(x) + sd(x), max(x))
    #   names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
    #   r
    # }
  # 
  # p1 <- ggplot(aes(y = value, x = factor(group)), data = large_data)
  # p1 <- p1 + stat_summary(fun.data = min.mean.sd.max, geom = "boxplot") +
  #   geom_jitter(position=position_jitter(width=.2), size=3) + 
  #   ggtitle("Boxplot con media, 95%CI, valore min. e max.") + 
  #   xlab("group") + ylab("Beta value")+
  theme( plot.title = element_text(color = 'black', size = 16, hjust = 0.5), 
         plot.subtitle = element_text(color = 'black', size = 16,hjust = 0.5), 
         plot.caption = element_text(color = 'black', size = 16,face = 'italic', hjust = 1), 
         axis.text.x = element_text(color = 'black', size = 16, angle = 0), 
         axis.text.y = element_text(color = 'black', size = 16, angle = 0), 
         axis.title.x = element_text(color = 'black', size = 16, angle = 0), 
         axis.title.y = element_text(color = 'black', size = 16, angle = 90), 
         legend.title = element_text(color = 'black', size = 16), 
         legend.text = element_text(color = 'black', size = 16), 
         axis.line.y = element_line(color = 'black', linetype = 'solid'), # y轴线特征 
         axis.line.x = element_line (color = 'black',linetype = 'solid'), # x轴线特征 
         panel.border = element_rect(linetype = 'solid', size = 1.2,fill = NA) # 图四周框起来 
  )
  ggsave(filename = "boxplot_AUC_topcv.jpg",plot = GGM,device ="jpeg" ,
         path = getwd(),dpi = 300,units = "in",width = 9, height = 8,
         limitsize=F)
}

ensemble_prediction_lp<-function (ensemble_model, prediction_data, mutiple_results = FALSE) 
{
  if (mutiple_results) {
    return(ensemble_prediction.m(ensemble_model, prediction_data))
  }
  prediction_data <- prediction_data[ensemble_model$cox$cox_model$CpGs, 
                                     ]
  if (nrow(prediction_data) != length(rownames(prediction_data))) {
    stop("ERROR: The predition data and the model have wrong dimensions!")
  }
  rank_svm <- stats::predict(ensemble_model$svm$svm_model, data.frame(t(prediction_data)))$predicted[1,]
  svm <- predict(ensemble_model$svm$hr_model,data.frame(rank_svm = rank_svm),type='lp')
  cox <- stats::predict(ensemble_model$cox$cox_model, data.frame(t(prediction_data)))
  enet <- stats::predict(ensemble_model$enet$enet_model, t(prediction_data), 
                         s = ensemble_model$enet$`corrected_lambda(min)`)
  mboost <- stats::predict(ensemble_model$mboost$mboost_model, 
                           t(prediction_data))[, 1]
  data <- rbind(cox, svm, t(enet), mboost)
  rownames(data) <- c("cox", "svm", "enet", "mboost")
  data <- t(data)
  data_f <- as.data.frame(data)
  if (class(ensemble_model$stacking)[1] == "cv.glmnet") 
    ensemble = as.numeric(stats::predict(ensemble_model$stacking, as.matrix(data_f)))
  else if (class(ensemble_model$stacking)[1] == "cph") 
    ensemble = stats::predict(ensemble_model$stacking, as.matrix(data_f),type="lp")
  else stop("ERROR: The ensemble predition model is invaild !")
  return(t(cbind(data,ensemble)))
}


survival_curve_p<-function (exp, living_days, living_events, write_name, title_name = "", 
                            threshold = NA, file = FALSE) 
{
  if (is.na(threshold)) threshold = median(exp)
  group_sur <- factor(1 * (exp > threshold), levels = 0:1, 
                      labels = c("low risk", "high risk"))
  pdata_gene <- list(time = living_days, event = living_events, 
                     gene = exp, group_sur = group_sur)
  fit <- survival::survfit(survival::Surv(time, event) ~ group_sur, 
                           data = pdata_gene)
  data.survdiff <- survival::survdiff(survival::Surv(time, 
                                                     event) ~ group_sur, data = pdata_gene)
  p.val = 1 - stats::pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  p.val
}


s_t_x<-function (Stime, status, marker, entry = NULL, predict.time) 
{
  span = 0.25*length(marker)^(-0.20)
  lambda = NULL
  times = Stime
  x <- marker
  if (is.null(entry)) 
    entry <- rep(0, length(times))
  bad <- is.na(times) | is.na(status) | is.na(x) | is.na(entry)
  entry <- entry[!bad]
  times <- times[!bad]
  status <- status[!bad]
  x <- x[!bad]
  #if (sum(bad) > 0) 
  #cat(paste("\n", sum(bad), "records with missing values dropped. \n"))
  cut.values <- unique(x)
  cut.values <- cut.values[order(cut.values)]
  ncuts <- length(cut.values)
  ooo <- order(times)
  times <- times[ooo]
  status <- status[ooo]
  x <- x[ooo]
  s0 <- 1
  unique.t0 <- unique(times)
  unique.t0 <- unique.t0[order(unique.t0)]
  n.times <- sum(unique.t0 <= predict.time)
  for (j in 1:n.times) {
    n <- sum(entry <= unique.t0[j] & times >= unique.t0[j])
    d <- sum((entry <= unique.t0[j]) & (times == unique.t0[j]) & 
               (status == 1))
    if (n > 0) 
      s0 <- s0 * (1 - d/n)
  }
  s.pooled <- s0
  roc.matrix <- matrix(NA, ncuts, 2)
  roc.matrix[ncuts, 1] <- 0
  roc.matrix[ncuts, 2] <- 1
  x.unique <- unique(x)
  x.unique <- x.unique[order(x.unique)]
  S.t.x <- rep(0, length(x.unique))
  t.evaluate <- unique(times[status == 1])
  t.evaluate <- t.evaluate[order(t.evaluate)]
  t.evaluate <- t.evaluate[t.evaluate <= predict.time]
  for (j in 1:length(x.unique)) {
    if (!is.null(span)) {
      ddd <- (x - x.unique[j])
      n <- length(x)
      ddd <- ddd[order(ddd)]
      index0 <- sum(ddd < 0) + 1
      index1 <- index0 + trunc(n * span + 0.5)
      if (index1 > n) 
        index1 <- n
      lambda <- ddd[index1]
      wt <- as.integer(((x - x.unique[j]) <= lambda) & 
                         ((x - x.unique[j]) >= 0))
      index0 <- sum(ddd <= 0)
      index2 <- index0 - trunc(n * span/2)
      if (index2 < 1) 
        index2 <- 1
      lambda <- abs(ddd[index2])
      set.index <- ((x - x.unique[j]) >= -lambda) & 
        ((x - x.unique[j]) <= 0)
      wt[set.index] <- 1
      
    }
    else {
      cat("requires a valid span! \n")
      stop(0)
    }
    s0 <- 1
    for (k in 1:length(t.evaluate)) {
      n <- sum(wt * (entry <= t.evaluate[k]) & (times >= 
                                                  t.evaluate[k]))
      d <- sum(wt * (entry <= t.evaluate[k]) & (times == 
                                                  t.evaluate[k]) * (status == 1))
      if (n > 0) 
        s0 <- s0 * (1 - d/n)
    }
    S.t.x[j] <- s0
  }
  S.all.x <- S.t.x[match(x, x.unique)]
  ks_res<-suppressWarnings(ks.test(x=S.all.x,y="exp",alternative = "t"))
  list(S.all.x = S.all.x,ks_res = ks_res)
}

multi_coxph <- function(dataframe,y_surv,digits=4,asnumeric=TRUE) {
  dataframe <- as.data.frame(dataframe)
  for (i in seq(ncol(dataframe))) {
    if (!is.numeric(dataframe[1,i]))  dataframe[,i] <-as.numeric(as.factor(dataframe[,i] ))
  }
  #print(head(dataframe))
  covariates <- paste(colnames(dataframe),collapse = ' + ')
  formula_str<- paste('y_surv',covariates,sep = ' ~ ')
  #cat(formula_str)
  multi_models<-survival::coxph(as.formula(formula_str),data = dataframe)
  multi_models %>% finalfit::fit2df()
}
