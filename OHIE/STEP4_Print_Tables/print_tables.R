## load notrim estimates

rm(list=ls())
library(xtable)
my_path<-"/net/holyparkesec/data/tata/leebounds/"

### read in results
## notrim results are in single table
table_notrim<-read.csv(paste0(my_path,"/OHIE/STEP2_Standard_Trim/csv/estimates_no_trim.csv"),row.names =1)
sd_notrim<-table_notrim[,4:6]


table_basic_trim<-read.csv(paste0(my_path,"/OHIE/STEP2_Standard_Trim/csv/estimates_basic_trim.csv"),row.names =1)
sd_basic_trim<-read.csv(paste0(my_path,"/OHIE/STEP2_Standard_Trim/csv/sd_basic_trim_clustered.csv"),row.names = 1)
sd_basic_trim<-apply(sd_basic_trim,2,round,4)


#table_trim_ml<-read.csv(paste0(my_path,"/OHIE/STEP3_ML/csv/estimates_ml_trim.csv"),row.names=1)

estimates_trim_discrete<-read.csv(paste0(my_path,"/OHIE/STEP3_ML/csv/estimates_ml_trim_discrete.csv"),row.names=1)
estimates_trim_numeric_health<-read.csv(paste0(my_path,"/OHIE/STEP3_ML/csv/estimates_ml_trim_health.csv"),row.names=1)
estimates_trim_utilization<-read.csv(paste0(my_path,"/OHIE/STEP3_ML/csv/estimates_ml_trim_utilization.csv"),row.names=1)

sd_trim_numeric_health<-read.csv(paste0(my_path,"/OHIE/STEP3_ML/csv/sd_ml_trim_numeric_health.csv"),row.names=1)
sd_trim_discrete<-read.csv(paste0(my_path,"/OHIE/STEP3_ML/csv/sd_ml_trim_discrete.csv"),row.names=1)
sd_trim_utilization<-read.csv(paste0(my_path,"/OHIE/STEP3_ML/csv/sd_ml_trim_utilization.csv"),row.names=1)

table_trim_ml<-rbind(estimates_trim_utilization,estimates_trim_numeric_health,estimates_trim_discrete)
sd_trim_ml<-rbind(sd_trim_utilization,sd_trim_numeric_health,sd_trim_discrete)


print_table<-function(names,digs) {

  estimates<-cbind(
                   table_notrim[match(names,rownames(table_notrim)),2], table_basic_trim[match(names,rownames(table_basic_trim)),2],table_trim_ml[match(names,rownames(table_trim_ml)),2],
                   table_notrim[match(names,rownames(table_notrim)),3], table_basic_trim[match(names,rownames(table_basic_trim)),3],table_trim_ml[match(names,rownames(table_trim_ml)),3])
  
  sd<-cbind( sd_notrim[match(names,rownames(sd_notrim)),2], sd_basic_trim[match(names,rownames(sd_basic_trim)),2],sd_trim_ml[match(names,rownames(sd_trim_ml)),2],
            sd_notrim[match(names,rownames(sd_notrim)),3], sd_basic_trim[match(names,rownames(sd_basic_trim)),3],sd_trim_ml[match(names,rownames(sd_trim_ml)),3])
  
  sd<-apply(sd,2,round,digs)
  estimates<-apply(estimates,2,round,digs)
  M<-matrix(NA,8,6)
  M[2*(1:dim(estimates)[1])-1,]<-apply(estimates,2,as.character)
  M[2*(1:dim(estimates)[1]),]<-apply(sd,2,function (x) paste0("(",x,")"))
  
  
  return(M)
}
print_table(rownames(table_notrim)[8:11],digs=3)
table<-print(xtable(print_table(rownames(table_notrim)[8:11],digs=3)),type="latex",include.rownames =FALSE )
write.table("Panel_A_Extensive_Margin",paste0(my_path,"/OHIE/STEP4_Print_Tables/utilization.txt"))
write.table(table,paste0(my_path,"/OHIE/STEP4_Print_Tables/utilization.txt"),append=FALSE)

table<-print(xtable(print_table(rownames(table_notrim)[1:4],digs=3)),type="latex",include.rownames =FALSE )
write.table("Panel_B_Total_Utilization",paste0(my_path,"/OHIE/STEP4_Print_Tables/utilization.txt"),append=TRUE)
write.table(table,paste0(my_path,"/OHIE/STEP4_Print_Tables/utilization.txt"),append=TRUE)



table<-print(xtable(print_table(rownames(table_notrim)[12:15],digs=3)),type="latex",include.rownames =FALSE )
write.table("Panel_A_Binary_Outcomes",paste0(my_path,"/OHIE/STEP4_Print_Tables/health.txt"))
write.table(table,paste0(my_path,"/OHIE/STEP4_Print_Tables/health.txt"),append=FALSE)

table<-print(xtable(print_table(rownames(table_notrim)[5:7],digs=3)),type="latex",include.rownames =FALSE )
write.table("Panel_B_Continuous_Outcomes",paste0(my_path,"/OHIE/STEP4_Print_Tables/health.txt"),append=TRUE)
write.table(table,paste0(my_path,"/OHIE/STEP4_Print_Tables/health.txt"),append=TRUE)

