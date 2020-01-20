library(data.table)
library(survival)
library(flexsurv)
library(dplyr)

d= fread('/mnt/hdd/data/mobaqs/p540/mfr_entire.csv')

d= filter(d, FLERFODSEL==0, 
		(DODKAT<6 | DODKAT> 10),
		SVLEN_UL_DG<308,
		SVLEN_UL_DG> 154,
		IVF== 0,
		ABRUPTIOP== 0,
		PLACENTA_PREVIA== 0,
		FOSTERV_POLYHYDRAMNION== 0,
		C00_MALF_ALL== 0)


ids= fread('/mnt/hdd/data/mobaqs/p540/mother_ids.csv')
ids= ids[sample(1:nrow(ids)),]

ids= ids[!duplicated(ids$M_ID),]

d= filter(d, PREG_ID_540 %in% ids$PREG_ID)

d= mutate(d, spont= as.numeric(FSTART==1 & (is.na(KSNITT) | KSNITT>1 | KSNITT == 0) &
                INDUKSJON_PROSTAGLANDIN==0 & 
		INDUKSJON_ANNET==0 &
                INDUKSJON_OXYTOCIN==0 & 
		INDUKSJON_AMNIOTOMI==0))


dists <- c("exp", "weibull", "gamma", "lognormal", "llogis")
dists_long= c('Exponential', 'Weibull', 'Gamma', 'Log-normal', 'Log-logistic')
parametric_haz <- vector(mode = "list", length = length(dists))
for (i in 1:length(dists)){
  fit <- flexsurvreg(Surv(SVLEN_UL_DG, spont) ~ 1, data = d, dist = dists[i]) 
  parametric_haz[[i]] <- summary(fit, type = "hazard", ci = FALSE, tidy = TRUE)
  parametric_haz[[i]]$method <- dists_long[i]
  print(fit$AIC)
}

parametric_haz <- rbindlist(parametric_haz)


