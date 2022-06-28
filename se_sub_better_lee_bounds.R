#library(hdm)
#library(devtools)
#library(leebounds)

data <- read.csv("test_data.csv", header=TRUE, stringsAsFactors=FALSE)
names(data)[1] <- 'treatment'

leedata <- data.frame(treat=data$treatment, selection=data$selection, outcome=data$outcome)
covars <- data.frame(straumxa1=data$stratumxa1, stratumxa2=data$stratumxa2, stratumxa3=data$stratumxa3)

#Function as it is defined on Github:
#estimate_selection<-function(leedata,form,selection_function,selection_function_name,variables_for_selection=NULL,names_to_include=c(),treat_name="treat+",yname="selection",myweights=NULL,...) 

#Trying two different versions, differences are on how form argument and names_to_include argument are defined
s_hat <- estimate_selection(leedata,outcome~treat,selection_function,"rlassologit",covars,names_to_include=c( ),treat_name="treat+",y_name="selection",myweights=NULL)
s_hat <- estimate_selection(leedata,yname~treat_name,selection_function,"rlassologit",covars,names_to_include=(stratumxa1 / stratumxa2 / stratumxa3),treat_name="treat+",yname="selection",myweights=NULL)

bounds <- leebounds_wout_monotonicity(leedata,s_hat)