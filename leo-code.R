
# Cargar librerias --------------------------------------------------------

install.packages("devtools")
library(devtools) 
install_github("vsemenova/leebounds")
library(leebounds)

# install.packages("rqPen") Lasso Logit
#install.packages("reldist")
# install.packages("nleqslv") 
# library(hdm)


# Load data ---------------------------------------------------------------

JobCorps_baseline <- leebounds::JobCorps_data_baseline
JobCorps_employment <- leebounds::JobCorps_data_employment
JobCorps_wages<- leebounds::JobCorps_data_wages

# compute basic Lee (2009) bounds for ATE in week 208
leedata_rep=data.frame(treat=JobCorps_baseline$TREATMNT.y,selection=JobCorps_employment$week_208,outcome=JobCorps_wages$week_208)
GetBounds(leebounds(leedata_rep))

# Estimate selection
estimate_selection(leedata_rep, form = leedata_rep$selection~leedata_rep$treat + leedata_rep$outcome, selection_function_name = "rlassologit",
                   variables_for_selection = c("outcome", "treat"))
