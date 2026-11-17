validate_numbers <- function(d) {
 measures<-c('Incidence','Deaths','YLLs (Years of Life Lost)','DALYs (Disability-Adjusted Life Years)')
 locs<-c('South Asia','Bangladesh','Bhutan','India','Nepal','Pakistan','Global','High-income')
 fields<-c('location_name','sex_name','year','measure_name')
 e<-expand.grid(location_name=locs,sex_name=SEXES,year=c(1990,2019,2023),measure_name=measures,stringsAsFactors=FALSE)
 if(anyDuplicated(d[fields])||!setequal(key(d,fields),key(e,fields)))stop('Incomplete or duplicate absolute burden grid')
 if(any(d$metric_name!='Number'|d$age_name!='15-49 years'|d$cause_name!='Tracheal, bronchus, and lung cancer'))stop('Wrong absolute burden extraction')
 if(any(!is.finite(as.matrix(d[c('val','lower','upper')]))))stop('Nonfinite numbers')
 if(any(d$lower<0|d$lower>d$val|d$val>d$upper))stop('Invalid number intervals')
 for(loc in locs)for(yr in c(1990,2019,2023))for(m in measures) {
  b<-d[d$location_name==loc & d$year==yr & d$measure_name==m,]
  total<-b$val[b$sex_name=='Both'];parts<-sum(b$val[b$sex_name!='Both'])
  stopifnot(abs(total-parts)<1e-7*max(1,total))
 }
 for(s in SEXES)for(yr in c(1990,2019,2023))for(m in measures) {
  b<-d[d$sex_name==s & d$year==yr & d$measure_name==m,]
  total<-b$val[b$location_name=='South Asia'];parts<-sum(b$val[b$location_name %in% c('Bangladesh','Bhutan','India','Nepal','Pakistan')])
  stopifnot(abs(total-parts)<1e-7*max(1,total))
 }
 for(loc in locs)for(s in SEXES)for(yr in c(1990,2019,2023)) {
  b<-d[d$location_name==loc & d$sex_name==s & d$year==yr,]
  stopifnot(b$val[b$measure_name==measures[4]]>=b$val[b$measure_name==measures[3]])
 }
 d
}
number_changes <- function(d) {
 out<-list();i<-1
 for(z in split(d,key(d,c('location_name','sex_name','measure_name'))))for(p in list(c(1990,2023),c(1990,2019),c(2019,2023))) {
  x<-z[z$year==p[1],];y<-z[z$year==p[2],]
  r<-x[1,c('location_name','sex_name','measure_name')];r$year_start<-p[1];r$year_end<-p[2]
  r$start_number<-x$val;r$end_number<-y$val;r$absolute_difference<-y$val-x$val
  r$percent_change_ratio_of_means<-100*(y$val/x$val-1)
  out[[i]]<-r;i<-i+1
 }
 do.call(rbind,out)
}
sex_contributions <- function(d) {
 out<-list();i<-1
 for(z in split(d,key(d,c('location_name','measure_name'))))for(p in list(c(1990,2023),c(1990,2019),c(2019,2023))) {
  total<-z$val[z$sex_name=='Both'&z$year==p[2]]-z$val[z$sex_name=='Both'&z$year==p[1]]
  for(s in c('Male','Female')) {
   delta<-z$val[z$sex_name==s&z$year==p[2]]-z$val[z$sex_name==s&z$year==p[1]]
   r<-z[1,c('location_name','measure_name')];r$sex<-s;r$year_start<-p[1];r$year_end<-p[2]
   r$absolute_number_difference<-delta;r$both_sex_number_difference<-total
   r$share_of_net_number_difference_percent<-if(total!=0)100*delta/total else NA_real_
   out[[i]]<-r;i<-i+1
  }
 }
 do.call(rbind,out)
}
age_specific_endpoints <- function(d,c) {
 if(any(c$year_start!=1990|c$year_end!=2023))stop('Wrong source change period')
 c$year<-c$year_end
 c<-normalize_rates(c);fields<-c('location','sex','measure','age')
 if(anyDuplicated(c[fields])||nrow(c)!=546)stop('Expected 546 unique source age-specific changes')
 e<-expand.grid(location=LOCATIONS,sex=SEXES,measure=c('Incidence','Deaths'),age=AGE_GROUPS,stringsAsFactors=FALSE)
 if(!setequal(key(c,fields),key(e,fields))||any(c$metric!='Rate'|c$cause!='Tracheal, bronchus, and lung cancer')||any(c$lower>c$val|c$val>c$upper))stop('Invalid source age-specific changes')
 out<-c[fields]
 for(yr in c(1990,2023)) {
  z<-d[d$year==yr,];ix<-match(key(c,fields),key(z,fields));if(anyNA(ix))stop('Missing age endpoint')
  for(f in c('val','lower','upper'))out[[paste0('rate_',yr,'_',f)]]<-z[[f]][ix]
 }
 for(f in c('val','lower','upper'))out[[paste0('GBD_source_relative_change_percent_',f)]]<-100*c[[f]]
 out$source_change_UI_excludes_zero<-c$lower>0|c$upper<0
 out$change_definition<-'GBD-provided age-specific relative change; not the custom age-standardized change'
 out
}
risk_tables <- function(v,c) {
 fields<-c('location','sex','rei')
 if(nrow(v)!=468||anyDuplicated(v[fields])||anyDuplicated(c[fields]))stop('Wrong risk grid')
 if(any(!is.finite(as.matrix(v[c('val','lower','upper')]))))stop('Nonfinite risk fractions')
 if(any(!is.finite(as.matrix(c[c('val','lower','upper')]))))stop('Nonfinite source risk changes')
 if(any(c$metric!='Percent'|c$measure!='Deaths'|c$age!='15-49 years'|c$year_start!=1990|c$year_end!=2023|c$cause!='Tracheal, bronchus, and lung cancer'))stop('Wrong risk change extraction')
 if(any(c$lower>c$val|c$val>c$upper))stop('Invalid source risk change intervals')
 if(any(v$metric!='Percent'|v$measure!='Deaths'|v$age!='15-49 years'|v$year!=2023|v$cause!='Tracheal, bronchus, and lung cancer'))stop('Wrong risk extraction')
 if(any(v$val<0|v$val>1|v$lower>v$val|v$val>v$upper))stop('Invalid risk summaries')
 e<-expand.grid(location=c('South Asia','Bangladesh','Bhutan','India','Nepal','Pakistan'),sex=SEXES,rei=unique(v$rei),stringsAsFactors=FALSE)
 if(!setequal(key(v,fields),key(e,fields))||length(unique(v$rei))!=26)stop('Incomplete risk categories')
 ix<-match(key(v,fields),key(c,fields));z<-c[ix,]
 if(nrow(c)!=458||sum(is.na(ix))!=10||any(v$rei[is.na(ix)]!='Occupational exposure to asbestos'))stop('Unexpected missing risk change records')
 out<-data.frame(location=v$location,sex=v$sex,risk=v$rei,age='15-49 years',fraction_2023_percent=100*v$val,fraction_2023_lower=100*v$lower,fraction_2023_upper=100*v$upper,fraction_2023_UI_includes_zero=v$lower<=0&v$upper>=0)
 for(f in c('val','lower','upper'))out[[paste0('relative_fraction_change_1990_2023_percent_',f)]]<-100*z[[f]]
 out$change_UI_relation_to_zero<-ifelse(is.na(ix),'No source change estimate',ifelse(z$lower>0,'Positive, UI excludes zero',ifelse(z$upper<0,'Negative, UI excludes zero','Includes zero')))
 out$note<-'Overlapping GBD modelled attributable fractions; not attributable death counts or age-standardized mortality rates'
 out
}
