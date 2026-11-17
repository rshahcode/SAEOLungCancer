# Exploratory regression inference; no GBD posterior uncertainty propagation.
# Reference specification fixed before the first fit. AR1 is a post-diagnostic sensitivity.
if(!requireNamespace('sandwich',quietly=TRUE)||!requireNamespace('nlme',quietly=TRUE))stop('Install sandwich and nlme')
NeweyWest<-sandwich::NeweyWest
fitone<-function(z,start=1990,end=2023) {
 z<-z[z$year>=start & z$year<=end,];z<-z[order(z$year),]
 stopifnot(identical(as.integer(z$year),seq.int(start,end)))
 x<-z$year-mean(z$year);y<-log(z$asr);f<-lm(y~x);n<-nrow(z);k<-2
 V<-NeweyWest(f,lag=3,prewhite=FALSE,adjust=TRUE)
 b<-unname(coef(f)[2]);se<-sqrt(V[2,2]);crit<-qt(.975,n-k)
 olsse<-sqrt(vcov(f)[2,2]);ep<-function(v)100*expm1(v)
 pval<-2*pt(-abs(b/se),n-k)
 quad<-lm(y~x+I(x^2));qV<-NeweyWest(quad,lag=3,prewhite=FALSE,adjust=TRUE)
 qp<-2*pt(-abs(coef(quad)[3]/sqrt(qV[3,3])),n-3)
 e<-resid(f);acf1<-unname(acf(e,lag.max=1)$acf[2]);dw<-sum(diff(e)^2)/sum(e^2)
 lb<-Box.test(e,lag=5,type='Ljung-Box',fitdf=0)$p.value
 out<-data.frame(location=z$location[1],sex=z$sex[1],measure=z$measure[1],age_range=z$age_range[1],start=start,end=end,n=n,
  intercept=unname(coef(f)[1]),slope=b,slope_se_hac3=se,slope_var_hac3=V[2,2],df=n-k,t_critical=crit,
  eapc=ep(b),ci_lower=ep(b-crit*se),ci_upper=ep(b+crit*se),p_hac3=pval,
  ols_ci_lower=ep(b-crit*olsse),ols_ci_upper=ep(b+crit*olsse),ols_p=2*pt(-abs(b/olsse),n-k),
  r_squared=summary(f)$r.squared,residual_acf1=acf1,durbin_watson=dw,ljung_box_lag5_p=lb,
  quadratic_hac3_p=unname(qp),linear_rmse=sqrt(mean(e^2)),quadratic_rmse=sqrt(mean(resid(quad)^2)),
  leave_last_year_out_eapc=ep(unname(coef(lm(y[-n]~x[-n]))[2])),
  interval_interpretation='Approximate regression CI conditional on annual point estimates; not propagated GBD UI')
 out$trend<-if(out$ci_lower>0)'Positive CI excludes zero' else if(out$ci_upper<0)'Negative CI excludes zero' else 'CI includes zero'
 for(L in c(1,2,4,5)) {
  ss<-sqrt(NeweyWest(f,lag=L,prewhite=FALSE,adjust=TRUE)[2,2])
  out[[paste0('hac',L,'_ci_lower')]]<-ep(b-crit*ss);out[[paste0('hac',L,'_ci_upper')]]<-ep(b+crit*ss)
 }
 out$curvature_flag<-out$quadratic_hac3_p<.05
 out$high_residual_persistence_flag<-abs(acf1)>.7
 out$lag_sensitivity_same_zero_classification<-all(sapply(c(1,2,4,5),function(L) {
  lo<-out[[paste0('hac',L,'_ci_lower')]];hi<-out[[paste0('hac',L,'_ci_upper')]]
  (lo>0)==(out$ci_lower>0) && (hi<0)==(out$ci_upper<0)
 }))

 out
}
run_trends<-function(d) {
 group<-function(z)interaction(z$location,z$sex,z$measure,z$age_range,drop=TRUE,lex.order=TRUE)
 primary<-d[d$age_range=='15-49',]
 results<-do.call(rbind,lapply(split(primary,group(primary)),fitone))
 stopifnot(nrow(results)==78)
 results$bh_p_78_primary_tests<-p.adjust(results$p_hac3,method='BH')
 ar1<-lapply(split(primary,group(primary)),function(z) {
  z<-z[order(z$year),];z$x<-z$year-mean(z$year);z$y<-log(z$asr)
  out<-data.frame(location=z$location[1],sex=z$sex[1],measure=z$measure[1],age_range='15-49',ar1_eapc=NA_real_,ar1_ci_lower=NA_real_,ar1_ci_upper=NA_real_,ar1_rho=NA_real_,ar1_p=NA_real_,status='Not fitted')
  f<-try(nlme::gls(y~x,data=z,correlation=nlme::corAR1(form=~year),method='REML',control=nlme::glsControl(msMaxIter=200)),silent=TRUE)
  if(inherits(f,'try-error')){out$status<-'GLS failed';return(out)}
  b<-coef(f)[2];se<-sqrt(vcov(f)[2,2]);crit<-qt(.975,nrow(z)-2)
  out$ar1_eapc<-100*expm1(b);out$ar1_ci_lower<-100*expm1(b-crit*se);out$ar1_ci_upper<-100*expm1(b+crit*se)
  out$ar1_rho<-as.numeric(coef(f$modelStruct$corStruct,unconstrained=FALSE));out$ar1_p<-2*pt(-abs(b/se),nrow(z)-2)
  out$status<-if(abs(out$ar1_rho)>.98)'Near-boundary rho; unstable model sensitivity' else 'Converged; conditional AR1 model sensitivity'
  out
 })
 sa<-d[d$location=='South Asia',];sens<-list();j<-1
 for(z in split(sa,group(sa)))for(end in c(2019,2023)) {
  if(z$age_range[1]=='15-49'&&end==2023)next
  sens[[j]]<-fitone(z,end=end);j<-j+1
 }
 sens<-do.call(rbind,sens);stopifnot(nrow(sens)==30)
 results<-results[order(results$location,results$sex,results$measure),]
 list(primary=results,sensitivity=sens,ar1=do.call(rbind,ar1))
}
