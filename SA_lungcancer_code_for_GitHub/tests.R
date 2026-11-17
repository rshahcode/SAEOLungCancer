script<-sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])
root<-dirname(normalizePath(script))
source(file.path(root,'standardization.R'))
source(file.path(root,'trends.R'))
# All test observations are synthetic and exist only in memory.
w<-data.frame(age_range='15-49',age=c('15-19 years','20-24 years'),restricted_weight=c(.4,.6))
d<-expand.grid(location='Synthetic',sex='Both',measure='Incidence',year=c(1990,2019,2023),age=w$age,stringsAsFactors=FALSE)
d$val<-ifelse(d$age==w$age[1],2,4)*(1+(d$year-1990)/100)
d$lower<-d$val*.8;d$upper<-d$val*1.2
a<-standardize_rates(d,w)
stopifnot(abs(a$asr[a$year==1990]-3.2)<1e-12)
stopifnot(all(abs(a$approximate_lower-.8*a$asr)<1e-12))
ch<-endpoint_changes(a)
stopifnot(abs(ch$percent_change[ch$year_start==1990&ch$year_end==2023]-33)<1e-10)
# Slight deterministic irregularity avoids a degenerate zero-residual covariance.
z<-data.frame(location='Synthetic',sex='Both',measure='Incidence',age_range='15-49',year=1990:2023)
z$asr<-exp(-1+.01*(z$year-1990)+.01*sin(seq_len(nrow(z))))
f<-fitone(z)
stopifnot(f$eapc>0,f$ci_lower<f$eapc,f$eapc<f$ci_upper)
# Independent Bartlett lag 3 sandwich calculation for the regression slope.
X<-cbind(1,z$year-mean(z$year));y<-log(z$asr)
beta<-solve(crossprod(X),crossprod(X,y));e<-as.numeric(y-X%*%beta);u<-X*e
S<-crossprod(u)
for(L in 1:3){G<-crossprod(u[(L+1):nrow(u),,drop=FALSE],u[1:(nrow(u)-L),,drop=FALSE]);S<-S+(1-L/4)*(G+t(G))}
B<-solve(crossprod(X));V<-B%*%S%*%B*nrow(X)/(nrow(X)-2)
stopifnot(abs(V[2,2]-f$slope_var_hac3)<1e-14)
cat('Synthetic standardization, descriptive change, and independent HAC checks passed.\n')
