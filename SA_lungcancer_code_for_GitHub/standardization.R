AGE_GROUPS <- paste0(seq(15,45,5),'-',seq(19,49,5),' years')
LOCATIONS <- c('Global','South Asia','Bangladesh','Bhutan','India','Nepal','Pakistan',
 'Central Europe, Eastern Europe, and Central Asia','High-income',
 'Latin America and Caribbean','North Africa and Middle East',
 'Southeast Asia, East Asia, and Oceania','Sub-Saharan Africa')
SEXES <- c('Both','Male','Female')
key <- function(d,fields) do.call(paste,c(d[fields],sep='\r'))
normalize_rates <- function(d) {
 for(nm in c('location','sex','measure','age','cause','metric')) {
  a<-paste0(nm,'_name'); if(!nm %in% names(d)&&a %in% names(d))names(d)[names(d)==a]<-nm
  if(!nm %in% names(d))stop('Missing field ',nm)
  d[[nm]]<-trimws(d[[nm]])
 }
 d$age<-gsub('\u2013','-',d$age,fixed=TRUE)
 d$location[d$location=='High-Income']<-'High-income'
 d$location[d$location=='Latin America and the Caribbean']<-'Latin America and Caribbean'
 d$location[d$location=='North Africa and the Middle East']<-'North Africa and Middle East'
 for(nm in c('year','val','lower','upper')) {
  d[[nm]]<-as.numeric(d[[nm]]);if(any(!is.finite(d[[nm]])))stop('Invalid ',nm)
 }
 d
}
validate_rates <- function(d) {
 fields<-c('location','sex','measure','year','age')
 e<-expand.grid(location=LOCATIONS,sex=SEXES,measure=c('Incidence','Deaths'),year=1990:2023,age=AGE_GROUPS,stringsAsFactors=FALSE)
 if(anyNA(d[fields])||anyDuplicated(d[fields])||!setequal(key(d,fields),key(e,fields)))stop('Incomplete, duplicate, or unexpected rate grid')
 if(any(d$cause!='Tracheal, bronchus, and lung cancer')||any(d$metric!='Rate'))stop('Wrong cause or metric')
 if(any(d$lower<0|d$lower>d$val|d$val>d$upper))stop('Invalid rate bounds')
 invisible(TRUE)
}
read_standard <- function(path) {
 if(!requireNamespace('readxl',quietly=TRUE))stop('Install readxl')
 x<-readxl::read_excel(path,col_names=FALSE,.name_repair='minimal')
 labels<-paste(seq(15,45,5),'to',seq(19,49,5))
 ix<-match(labels,trimws(as.character(x[[1]])))
 if(anyNA(ix)||anyDuplicated(ix))stop('Cannot find standard age rows in Appendix S14')
 values<-as.numeric(x[[2]][ix])
 if(any(!is.finite(values)|values<=0))stop('Invalid world standard weights')
 out<-list();j<-1
 for(rng in c('15-49','15-39','40-49')) {
  keep<-switch(rng,'15-49'=1:7,'15-39'=1:5,'40-49'=6:7)
  out[[j]]<-data.frame(age_range=rng,age=AGE_GROUPS[keep],standard_percent=values[keep],restricted_weight=values[keep]/sum(values[keep]))
  j<-j+1
 }
 do.call(rbind,out)
}
standardize_rates <- function(d,w) {
 out<-list();i<-1
 for(rng in unique(w$age_range)) {
  ww<-w[w$age_range==rng,]
  stopifnot(abs(sum(ww$restricted_weight)-1)<1e-12)
  z<-d[d$age %in% ww$age,]
  for(b in split(z,key(z,c('location','sex','measure','year')))) {
   b<-b[match(ww$age,b$age),];stopifnot(nrow(b)==nrow(ww),!anyNA(b$val))
   r<-b[1,c('location','sex','measure','year')];r$age_range<-rng
   r$asr<-sum(b$val*ww$restricted_weight)
   r$approximate_lower<-sum(b$lower*ww$restricted_weight)
   r$approximate_upper<-sum(b$upper*ww$restricted_weight)
   r$interval_method<-'Approximate weighted marginal bounds; not formal 95% UIs'
   out[[i]]<-r;i<-i+1
  }
 }
 do.call(rbind,out)
}
endpoint_changes <- function(a) {
 out<-list();i<-1
 for(z in split(a,key(a,c('location','sex','measure','age_range'))))for(p in list(c(1990,2023),c(1990,2019),c(2019,2023))) {
  x<-z[z$year==p[1],];y<-z[z$year==p[2],];stopifnot(nrow(x)==1,nrow(y)==1,x$asr>0)
  r<-x[1,c('location','sex','measure','age_range')];r$year_start<-p[1];r$year_end<-p[2]
  r$start_asr<-x$asr;r$end_asr<-y$asr;r$percent_change<-100*(y$asr/x$asr-1)
  r$method<-'Descriptive ratio of standardized point estimates; no change interval or test'
  out[[i]]<-r;i<-i+1
 }
 do.call(rbind,out)
}

