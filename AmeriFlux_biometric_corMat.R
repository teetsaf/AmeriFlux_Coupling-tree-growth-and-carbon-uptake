

#####
##### Monthly correlation matrix  
##### Annual variable (wood biomass increment) vs. monthly variables (NEP, GPP)
##### 
##### Code written by Aaron Teets with help from Bijan Seyednasrollah
##### May 2021


#libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)
library(lubridate)
library(boot)


### set working drive - locate data files in the same folder where script is located 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

### import biometric data
biom=read.csv('wood_biomass_inc.csv')
head(biom)
biom_dt=data.table(biom[,1:3])

ggplot(biom_dt,aes(x=year,y=inc))+
  geom_point()+  geom_line()+ theme_bw()+ 
  ylab("Annual biomass increment in grams of carbon per meter squared")+
  facet_wrap(~site)

### import daily eddy covariance data
dayflux=read.csv('daily_carbon_uptake.csv')
dayflux$date <- as.Date(with(dayflux, paste(year,doy,sep="-")), "%Y-%j")
dayflux$month = as.integer(format(dayflux$date,"%m"))
dayflux=subset(dayflux,dayflux$doy<366 & dayflux$year<2018)
head(dayflux)

ggplot(dayflux,aes(x=date,y=NEE_U50))+
  geom_point()+ theme_bw()+ 
  ylab("Daily net ecosystem exchange in grams of carbon per meter squared")+
  facet_wrap(~site)

ggplot(dayflux,aes(x=date,y=GPP_U50))+
  geom_point()+ theme_bw()+ 
  ylab("Daily gross primary productivity in grams of carbon per meter squared")+
  facet_wrap(~site)


### integrate daily eddy covariance data to monthly data
df=data.table(dayflux)
monthflux=df[,.(NEP=sum(-NEE_U50),GPP=sum(GPP_U50),Reco=sum(Reco_U50)),.(site,year,month)]
monthflux


### import monthly climate data
clim=read.csv('monthly_climate.csv')
clim$date <- as.Date(with(clim, paste(year,month,"01",sep="-")), "%Y-%m-%d")
head(clim)

### merge eddy covariance data and climate drivers
drivers=merge(clim,monthflux)
drivers=drivers[order(drivers$site,drivers$date),]
head(drivers)

dt0=data.table(drivers)
dt0


#### Set up correlation, p-value, uncertainty, and slope matrices

summ_dt_orig2=rbind(data.frame(site='US.Bar',year=2003),unique(dt0[,c(1,2)]))
summ_dt_orig=as.data.frame(summ_dt_orig2 %>%
                             group_by(site) %>%
                             filter(year > min(year)))
summ_dt_orig
Sites=factor(unique(summ_dt_orig$site))

corMat1 <- matrix(NA, nrow = 24, ncol = 24)
corMat2 <- matrix(NA, nrow = 24, ncol = 24)
corMat3 <- matrix(NA, nrow = 24, ncol = 24)
corMat4 <- matrix(NA, nrow = 24, ncol = 24)
corMat5 <- matrix(NA, nrow = 24, ncol = 24)
corMat6 <- matrix(NA, nrow = 24, ncol = 24)
corMat11 <- matrix(NA, nrow = 24, ncol = 24)
corMat12 <- matrix(NA, nrow = 24, ncol = 24)
corMat13 <- matrix(NA, nrow = 24, ncol = 24)
corMat14 <- matrix(NA, nrow = 24, ncol = 24)
corMat15 <- matrix(NA, nrow = 24, ncol = 24)
corMat16 <- matrix(NA, nrow = 24, ncol = 24)
corMat21 <- matrix(NA, nrow = 24, ncol = 24)
corMat22 <- matrix(NA, nrow = 24, ncol = 24)
corMat23 <- matrix(NA, nrow = 24, ncol = 24)
corMat24 <- matrix(NA, nrow = 24, ncol = 24)
corMat25 <- matrix(NA, nrow = 24, ncol = 24)
corMat26 <- matrix(NA, nrow = 24, ncol = 24)
corMat31 <- matrix(NA, nrow = 24, ncol = 24)
corMat32 <- matrix(NA, nrow = 24, ncol = 24)
corMat33 <- matrix(NA, nrow = 24, ncol = 24)
corMat34 <- matrix(NA, nrow = 24, ncol = 24)
corMat35 <- matrix(NA, nrow = 24, ncol = 24)
corMat36 <- matrix(NA, nrow = 24, ncol = 24)

pb <- txtProgressBar(min = 1, max = 24, style = 3)

### load functions 
## assembles the datasets to run correlations
f <- function(s, y, i, j){
  refmo <- as.Date(paste0(y - 1 , '-01-01'))
  as.numeric(dt0[site == s & date >= (refmo %m+% months(i)) & date <= (refmo %m+% months(j)),
                .(sum(GPP),sum(NEP),sum(Reco),sum(precip),mean(temp),mean(Rg))])}

## bootstrap uncertainty
cor.mu <- function(df, n) {
  df = df[n,]
  x <- df$inc
  y <- df$V1
  return(cor(x, y))
}


## Correlation Matrix for loop
## runs correlations, p-values, bootstrapped uncertainty, and slope 
## for all combinations of months for 6 different variables

Sites

for(h in levels(Sites)){

summ_dt_orig0=subset(summ_dt_orig,summ_dt_orig$site==paste(h))  
  
#test different combinations of months
for(i in 1:24){
    for(j in i:24){
     summ_dt <- data.table(summ_dt_orig0,t(mapply(FUN = f,
                                    y = summ_dt_orig0$year,
                                    s = summ_dt_orig0$site,
                                    i = i,
                                    j = j)))
    if(nrow(summ_dt)==0) next()
    dt <- merge(biom_dt[,.(inc = unique(inc)),.(year, site)],summ_dt)
    #correlations
    corMat1[i, j] = dt[,as.numeric(cor(inc, V1))]
    corMat2[i, j] = dt[,as.numeric(cor(inc, V2))]
    corMat3[i, j] = dt[,as.numeric(cor(inc, V3))]
    corMat4[i, j] = dt[,as.numeric(cor(inc, V4))]
    corMat5[i, j] = dt[,as.numeric(cor(inc, V5))]
    corMat6[i, j] = dt[,as.numeric(cor(inc, V6))]
    #p-values
    corMat11[i, j] = dt[,cor.test(inc,V1)$p.value]
    corMat12[i, j] = dt[,cor.test(inc,V2)$p.value]
    corMat13[i, j] = dt[,cor.test(inc,V3)$p.value]
    corMat14[i, j] = dt[,cor.test(inc,V4)$p.value]
    corMat15[i, j] = dt[,cor.test(inc,V5)$p.value]
    corMat16[i, j] = dt[,cor.test(inc,V6)$p.value]
    #uncertainty
    corMat21[i, j] = sd(boot(dt, cor.mu, R = 500)$t)
    corMat22[i, j] = sd(boot(dt, cor.mu, R = 500)$t)
    corMat23[i, j] = sd(boot(dt, cor.mu, R = 500)$t)
    corMat24[i, j] = sd(boot(dt, cor.mu, R = 500)$t)
    corMat25[i, j] = sd(boot(dt, cor.mu, R = 500)$t)
    corMat26[i, j] = sd(boot(dt, cor.mu, R = 500)$t)
    #slopes
    corMat31[i, j] = dt[,coef(lm(inc~V1))[2]]
    corMat32[i, j] = dt[,coef(lm(inc~V2))[2]]
    corMat33[i, j] = dt[,coef(lm(inc~V3))[2]]
    corMat34[i, j] = dt[,coef(lm(inc~V4))[2]]
    corMat35[i, j] = dt[,coef(lm(inc~V5))[2]]
    corMat36[i, j] = dt[,coef(lm(inc~V6))[2]]
    }
  setTxtProgressBar(pb, value = i)
}


#### Format correlation matrices

cmGPP=data.frame(corMat1)
cmGPP$start=1:24
cmGPP=gather(cmGPP,end2,r,-start,na.rm=T)
cmGPP$end=as.numeric(substring(cmGPP$end2,2))
cmGPP2=data.frame(start=cmGPP$start,end=cmGPP$end,r=cmGPP$r)
cmGPP2=cmGPP2[order(cmGPP2$start,cmGPP2$end),]
cmGPP2$length=cmGPP2$end-cmGPP2$start
cmGPP2$meas='GPP'

cmNEP=data.frame(corMat2)
cmNEP$start=1:24
cmNEP=gather(cmNEP,end2,r,-start,na.rm=T)
cmNEP$end=as.numeric(substring(cmNEP$end2,2))
cmNEP2=data.frame(start=cmNEP$start,end=cmNEP$end,r=cmNEP$r)
cmNEP2=cmNEP2[order(cmNEP2$start,cmNEP2$end),]
cmNEP2$length=cmNEP2$end-cmNEP2$start
cmNEP2$meas='NEP'

cmReco=data.frame(corMat3)
cmReco$start=1:24
cmReco=gather(cmReco,end2,r,-start,na.rm=T)
cmReco$end=as.numeric(substring(cmReco$end2,2))
cmReco2=data.frame(start=cmReco$start,end=cmReco$end,r=cmReco$r)
cmReco2=cmReco2[order(cmReco2$start,cmReco2$end),]
cmReco2$length=cmReco2$end-cmReco2$start
cmReco2$meas='Reco'

cmprecip=data.frame(corMat4)
cmprecip$start=1:24
cmprecip=gather(cmprecip,end2,r,-start,na.rm=T)
cmprecip$end=as.numeric(substring(cmprecip$end2,2))
cmprecip2=data.frame(start=cmprecip$start,end=cmprecip$end,r=cmprecip$r)
cmprecip2=cmprecip2[order(cmprecip2$start,cmprecip2$end),]
cmprecip2$length=cmprecip2$end-cmprecip2$start
cmprecip2$meas='Precipitation'

cmtemp=data.frame(corMat5)
cmtemp$start=1:24
cmtemp=gather(cmtemp,end2,r,-start,na.rm=T)
cmtemp$end=as.numeric(substring(cmtemp$end2,2))
cmtemp2=data.frame(start=cmtemp$start,end=cmtemp$end,r=cmtemp$r)
cmtemp2=cmtemp2[order(cmtemp2$start,cmtemp2$end),]
cmtemp2$length=cmtemp2$end-cmtemp2$start
cmtemp2$meas='Temperature'

cmRg=data.frame(corMat6)
cmRg$start=1:24
cmRg=gather(cmRg,end2,r,-start,na.rm=T)
cmRg$end=as.numeric(substring(cmRg$end2,2))
cmRg2=data.frame(start=cmRg$start,end=cmRg$end,r=cmRg$r)
cmRg2=cmRg2[order(cmRg2$start,cmRg2$end),]
cmRg2$length=cmRg2$end-cmRg2$start
cmRg2$meas='Rg'

cm=rbind(cmGPP2,cmNEP2,cmReco2,cmprecip2,cmtemp2,cmRg2)
cm$meas=factor(cm$meas,levels=c('GPP','NEP','Reco','Precipitation','Temperature','Rg'))


##### Format p-value matrices

pmGPP=data.frame(corMat11)
pmGPP$start=1:24
pmGPP=gather(pmGPP,end2,p.value,-start,na.rm=T)
pmGPP$end=as.numeric(substring(pmGPP$end2,2))
pmGPP2=data.frame(start=pmGPP$start,end=pmGPP$end,p.value=pmGPP$p.value)
pmGPP2=pmGPP2[order(pmGPP2$start,pmGPP2$end),]
pmGPP2$length=pmGPP2$end-pmGPP2$start
pmGPP2$meas='GPP'

pmNEP=data.frame(corMat12)
pmNEP$start=1:24
pmNEP=gather(pmNEP,end2,p.value,-start,na.rm=T)
pmNEP$end=as.numeric(substring(pmNEP$end2,2))
pmNEP2=data.frame(start=pmNEP$start,end=pmNEP$end,p.value=pmNEP$p.value)
pmNEP2=pmNEP2[order(pmNEP2$start,pmNEP2$end),]
pmNEP2$length=pmNEP2$end-pmNEP2$start
pmNEP2$meas='NEP'

pmReco=data.frame(corMat13)
pmReco$start=1:24
pmReco=gather(pmReco,end2,p.value,-start,na.rm=T)
pmReco$end=as.numeric(substring(pmReco$end2,2))
pmReco2=data.frame(start=pmReco$start,end=pmReco$end,p.value=pmReco$p.value)
pmReco2=pmReco2[order(pmReco2$start,pmReco2$end),]
pmReco2$length=pmReco2$end-pmReco2$start
pmReco2$meas='Reco'

pmprecip=data.frame(corMat14)
pmprecip$start=1:24
pmprecip=gather(pmprecip,end2,p.value,-start,na.rm=T)
pmprecip$end=as.numeric(substring(pmprecip$end2,2))
pmprecip2=data.frame(start=pmprecip$start,end=pmprecip$end,p.value=pmprecip$p.value)
pmprecip2=pmprecip2[order(pmprecip2$start,pmprecip2$end),]
pmprecip2$length=pmprecip2$end-pmprecip2$start
pmprecip2$meas='Precipitation'

pmtemp=data.frame(corMat15)
pmtemp$start=1:24
pmtemp=gather(pmtemp,end2,p.value,-start,na.rm=T)
pmtemp$end=as.numeric(substring(pmtemp$end2,2))
pmtemp2=data.frame(start=pmtemp$start,end=pmtemp$end,p.value=pmtemp$p.value)
pmtemp2=pmtemp2[order(pmtemp2$start,pmtemp2$end),]
pmtemp2$length=pmtemp2$end-pmtemp2$start
pmtemp2$meas='Temperature'

pmRg=data.frame(corMat16)
pmRg$start=1:24
pmRg=gather(pmRg,end2,p.value,-start,na.rm=T)
pmRg$end=as.numeric(substring(pmRg$end2,2))
pmRg2=data.frame(start=pmRg$start,end=pmRg$end,p.value=pmRg$p.value)
pmRg2=pmRg2[order(pmRg2$start,pmRg2$end),]
pmRg2$length=pmRg2$end-pmRg2$start
pmRg2$meas='Rg'

pm=rbind(pmGPP2,pmNEP2,pmReco2,pmprecip2,pmtemp2,pmRg2)
pm$meas=factor(pm$meas,levels=c('GPP','NEP','Reco','Precipitation','Temperature','Rg'))


#### Format uncertainty matrices

umGPP=data.frame(corMat21)
umGPP$start=1:24
umGPP=gather(umGPP,end2,uncertainty,-start,na.rm=T)
umGPP$end=as.numeric(substring(umGPP$end2,2))
umGPP2=data.frame(start=umGPP$start,end=umGPP$end,uncertainty=umGPP$uncertainty)
umGPP2=umGPP2[order(umGPP2$start,umGPP2$end),]
umGPP2$length=umGPP2$end-umGPP2$start
umGPP2$meas='GPP'

umNEP=data.frame(corMat22)
umNEP$start=1:24
umNEP=gather(umNEP,end2,uncertainty,-start,na.rm=T)
umNEP$end=as.numeric(substring(umNEP$end2,2))
umNEP2=data.frame(start=umNEP$start,end=umNEP$end,uncertainty=umNEP$uncertainty)
umNEP2=umNEP2[order(umNEP2$start,umNEP2$end),]
umNEP2$length=umNEP2$end-umNEP2$start
umNEP2$meas='NEP'

umReco=data.frame(corMat23)
umReco$start=1:24
umReco=gather(umReco,end2,uncertainty,-start,na.rm=T)
umReco$end=as.numeric(substring(umReco$end2,2))
umReco2=data.frame(start=umReco$start,end=umReco$end,uncertainty=umReco$uncertainty)
umReco2=umReco2[order(umReco2$start,umReco2$end),]
umReco2$length=umReco2$end-umReco2$start
umReco2$meas='Reco'

umprecip=data.frame(corMat24)
umprecip$start=1:24
umprecip=gather(umprecip,end2,uncertainty,-start,na.rm=T)
umprecip$end=as.numeric(substring(umprecip$end2,2))
umprecip2=data.frame(start=umprecip$start,end=umprecip$end,uncertainty=umprecip$uncertainty)
umprecip2=umprecip2[order(umprecip2$start,umprecip2$end),]
umprecip2$length=umprecip2$end-umprecip2$start
umprecip2$meas='Precipitation'

umtemp=data.frame(corMat25)
umtemp$start=1:24
umtemp=gather(umtemp,end2,uncertainty,-start,na.rm=T)
umtemp$end=as.numeric(substring(umtemp$end2,2))
umtemp2=data.frame(start=umtemp$start,end=umtemp$end,uncertainty=umtemp$uncertainty)
umtemp2=umtemp2[order(umtemp2$start,umtemp2$end),]
umtemp2$length=umtemp2$end-umtemp2$start
umtemp2$meas='Temperature'

umRg=data.frame(corMat26)
umRg$start=1:24
umRg=gather(umRg,end2,uncertainty,-start,na.rm=T)
umRg$end=as.numeric(substring(umRg$end2,2))
umRg2=data.frame(start=umRg$start,end=umRg$end,uncertainty=umRg$uncertainty)
umRg2=umRg2[order(umRg2$start,umRg2$end),]
umRg2$length=umRg2$end-umRg2$start
umRg2$meas='Rg'

um=rbind(umGPP2,umNEP2,umReco2,umprecip2,umtemp2,umRg2)
um$meas=factor(um$meas,levels=c('GPP','NEP','Reco','Precipitation','Temperature','Rg'))


##### Format slope matrices

smGPP=data.frame(corMat31)
smGPP$start=1:24
smGPP=gather(smGPP,end2,slope,-start,na.rm=T)
smGPP$end=as.numeric(substring(smGPP$end2,2))
smGPP2=data.frame(start=smGPP$start,end=smGPP$end,slope=smGPP$slope)
smGPP2=smGPP2[order(smGPP2$start,smGPP2$end),]
smGPP2$length=smGPP2$end-smGPP2$start
smGPP2$meas='GPP'

smNEP=data.frame(corMat32)
smNEP$start=1:24
smNEP=gather(smNEP,end2,slope,-start,na.rm=T)
smNEP$end=as.numeric(substring(smNEP$end2,2))
smNEP2=data.frame(start=smNEP$start,end=smNEP$end,slope=smNEP$slope)
smNEP2=smNEP2[order(smNEP2$start,smNEP2$end),]
smNEP2$length=smNEP2$end-smNEP2$start
smNEP2$meas='NEP'

smReco=data.frame(corMat33)
smReco$start=1:24
smReco=gather(smReco,end2,slope,-start,na.rm=T)
smReco$end=as.numeric(substring(smReco$end2,2))
smReco2=data.frame(start=smReco$start,end=smReco$end,slope=smReco$slope)
smReco2=smReco2[order(smReco2$start,smReco2$end),]
smReco2$length=smReco2$end-smReco2$start
smReco2$meas='Reco'

smprecip=data.frame(corMat34)
smprecip$start=1:24
smprecip=gather(smprecip,end2,slope,-start,na.rm=T)
smprecip$end=as.numeric(substring(smprecip$end2,2))
smprecip2=data.frame(start=smprecip$start,end=smprecip$end,slope=smprecip$slope)
smprecip2=smprecip2[order(smprecip2$start,smprecip2$end),]
smprecip2$length=smprecip2$end-smprecip2$start
smprecip2$meas='Precipitation'

smtemp=data.frame(corMat35)
smtemp$start=1:24
smtemp=gather(smtemp,end2,slope,-start,na.rm=T)
smtemp$end=as.numeric(substring(smtemp$end2,2))
smtemp2=data.frame(start=smtemp$start,end=smtemp$end,slope=smtemp$slope)
smtemp2=smtemp2[order(smtemp2$start,smtemp2$end),]
smtemp2$length=smtemp2$end-smtemp2$start
smtemp2$meas='Temperature'

smRg=data.frame(corMat36)
smRg$start=1:24
smRg=gather(smRg,end2,slope,-start,na.rm=T)
smRg$end=as.numeric(substring(smRg$end2,2))
smRg2=data.frame(start=smRg$start,end=smRg$end,slope=smRg$slope)
smRg2=smRg2[order(smRg2$start,smRg2$end),]
smRg2$length=smRg2$end-smRg2$start
smRg2$meas='Rg'

sm=rbind(smGPP2,smNEP2,smReco2,smprecip2,smtemp2,smRg2)
sm$meas=factor(sm$meas,levels=c('GPP','NEP','Reco','Precipitation','Temperature','Rg'))

### Combine correlation, p-value, uncertainty, slope data into one data frame per site
bm=merge(cm,pm)
nm=merge(bm,um)
mm=merge(nm,sm)
mm$site=1
head(mm)
mm$site=paste(substr(h,4,6))
assign(paste(substr(h,4,6),'_mm',sep=''),mm)

}



### Re-name sites
Bar_mm$site='US-Bar'
Ha1_mm$site='US-Ha1'
Ho1_mm$site='US-Ho1'
NR1_mm$site='US-NR1'
UMB_mm$site='US-UMB'
MMS_mm$site='US-MMS'

### Combine into one data frame 
CM=rbind(NR1_mm,Ho1_mm,Bar_mm,UMB_mm,Ha1_mm,MMS_mm)
CM$site=factor(CM$site,levels=c('US-NR1','US-Ho1','US-Bar','US-UMB','US-Ha1','US-MMS'))




### save for loop output
save(CM,file='AmeriFlux_biometric_corr_mat.Rdata')

### or

### load for loop output

load('AmeriFlux_biometric_corr_mat.Rdata')





### plot correlation matrices

### NEP and GPP plots

data_text <- data.frame(
  site = as.factor(c('US-NR1','US-NR1', 'US-Ho1', 'US-Ho1', 'US-Bar', 'US-Bar','US-UMB', 'US-UMB','US-Ha1','US-Ha1', 'US-MMS', 'US-MMS')),
  meas = as.factor(c('NEP','GPP','NEP','GPP','NEP','GPP','NEP','GPP','NEP','GPP','NEP','GPP')),
  start = 21,
  end = 3,
  label = c('(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)', '(i)', '(j)', '(k)', '(l)'),
  r=0)

CM2=subset(CM,CM$meas=='GPP' | CM$meas=='NEP')
CM2$meas=factor(CM2$meas,levels=c('NEP','GPP'))
CM2$site=factor(CM2$site,levels=c('US-NR1','US-Ho1','US-Bar','US-UMB','US-Ha1','US-MMS'))

flux_p=ggplot(CM2, aes(x = start, y = end, z = r,label = ifelse(p.value < 0.05, "+", "")))+
  geom_raster(aes(fill = slope))+
  scale_fill_gradientn(colours=c("blue3","lightsteelblue2","orange2","brown2"),
  breaks = c(-8,-4,-2,0.0,2,4,10),
                     limits=c(-1,1))+
  #scale_fill_gradientn(colours=c("blue3","lightsteelblue2","white","white","orange2","brown2"),
  #                     breaks = c(-1,-0.75,-0.50,-0.25, 0.0,0.25,0.5,0.75,1),
  #                     limits=c(-1,1))+
  geom_abline(slope=1,intercept=11,size=1.15,alpha=0.6,linetype=2)+
  geom_abline(slope=1,intercept=0,size=1.15,alpha=0.9)+
  geom_vline(xintercept=12.5,linetype=3)+
  geom_hline(yintercept=12.5,linetype=3)+
  #geom_point(aes(x=max_x,y=max_y),shape=23,fill='white',color='red',size=3.5)+
  geom_point(x=6,y=8,shape=21,size=5)+
  geom_point(x=9,y=11,shape=22,size=4)+
  geom_point(x=15,y=17,shape=23,size=4)+
  geom_point(x=18,y=20,shape=24,size=4)+
  #geom_text(size=5)+
  geom_text(data=data_text,aes(x=start, y=end, label=label), fontface='bold', size=6.5)+
  theme_classic(base_size=20)+
  theme(legend.position='bottom',legend.title = element_text(size = 18),
        legend.text = element_text(size = 15))+
  guides(label.position="bottom",fill = guide_colourbar(barwidth = 24, barheight = 1.15))+
  labs(x="Start Date", y="End Date", fill ="Correlation\nCoefficient (r)    ")+
  scale_x_continuous(expand = c(0, 0),
                     breaks = c(0, 3, 6, 9, 12, 15, 18, 21,24),
                     label = c("pJAN", "pMAR", "pJUN","pSEP","pDEC","MAR","JUN","SEP","DEC"))+
  scale_y_continuous(expand = c(0, 0),
                     breaks = c(0, 3, 6, 9, 12, 15, 18, 21,24),
                     label = c("pJAN", "pMAR", "pJUN","pSEP","pDEC","MAR","JUN","SEP","DEC"))+
  theme(axis.text.x = element_text(angle=45,hjust=1),
        legend.background = element_rect(fill = NA),
        panel.border=element_rect(fill=NA,colour="black",size=1),
        strip.background = element_blank(),
        strip.text = element_text(size = 20))+
  guides(color = guide_legend(override.aes = list(size = 0.85)))+
  guides(color = guide_legend(label.position = "left", label.hjust = 1))+
  facet_grid(rows=vars(site),cols=vars(meas))
flux_p


#jpeg("Fig.4.jpeg", width = 8, height = 20, units = 'in',res=1000)
#flux_p
#dev.off()


### Precipitation, Temperature, and Solar Radiation plots

data_text2 <- data.frame(
  site = as.factor(c('US-NR1','US-NR1', 'US-NR1','US-Ho1', 'US-Ho1', 'US-Ho1', 'US-Bar', 'US-Bar','US-Bar', 'US-UMB', 'US-UMB','US-UMB','US-Ha1','US-Ha1','US-Ha1', 'US-MMS', 'US-MMS', 'US-MMS')),
  meas = as.factor(c('Precipitation','Temperature','Solar radiation','Precipitation','Temperature','Solar radiation','Precipitation','Temperature','Solar radiation','Precipitation','Temperature','Solar radiation','Precipitation','Temperature','Solar radiation','Precipitation','Temperature','Solar radiation')),
  start = 21,
  end = 3,
  label = c('(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)', '(i)', '(j)', '(k)', '(l)', '(m)', '(n)', '(o)', '(p)', '(q)', '(r)'),
  r=0)

CM3=subset(CM,CM$meas=='Precipitation' | CM$meas=='Temperature' | CM$meas=='Rg')
CM3$meas=ifelse(CM3$meas=='Rg','Solar radiation',as.character(CM3$meas))
CM3$meas=factor(CM3$meas,levels=c('Precipitation','Temperature','Solar radiation'))

clim_p=ggplot(CM3, aes(x = start, y = end, z = r,label = ifelse(p.value < 0.05, "+", "")))+
  geom_raster(aes(fill = r))+
  scale_fill_gradientn(colours=c("blue3","lightsteelblue2","white","white","orange2","brown2"),
                       breaks = c(-1,-0.75,-0.50,-0.25, 0.0,0.25,0.5,0.75,1),
                       limits=c(-1,1))+
  geom_abline(slope=1,intercept=11,size=1.15,alpha=0.6,linetype=2)+
  geom_abline(slope=1,intercept=0,size=1.15,alpha=0.9)+
  geom_vline(xintercept=12.5,linetype=3)+
  geom_hline(yintercept=12.5,linetype=3)+
  #geom_point(aes(x=max_x,y=max_y),shape=23,fill='white',color='red',size=3.5)+
  geom_point(x=6,y=8,shape=21,size=5)+
  geom_point(x=9,y=11,shape=22,size=4)+
  geom_point(x=15,y=17,shape=23,size=4)+
  geom_point(x=18,y=20,shape=24,size=4)+
  geom_text(size=5)+
  geom_text(data=data_text2,aes(x=start, y=end, label=label), fontface='bold', size=6.5)+
  theme_classic(base_size=20)+
  theme(legend.position='bottom',legend.title = element_text(size = 18),
        legend.text = element_text(size = 15))+
  guides(label.position="bottom",fill = guide_colourbar(barwidth = 24, barheight = 1.15))+
  labs(x="Start Date", y="End Date", fill ="Correlation\nCoefficient (r)    ")+
  scale_x_continuous(expand = c(0, 0),
                     breaks = c(0, 3, 6, 9, 12, 15, 18, 21,24),
                     label = c("pJAN", "pMAR", "pJUN","pSEP","pDEC","MAR","JUN","SEP","DEC"))+
  scale_y_continuous(expand = c(0, 0),
                     breaks = c(0, 3, 6, 9, 12, 15, 18, 21,24),
                     label = c("pJAN", "pMAR", "pJUN","pSEP","pDEC","MAR","JUN","SEP","DEC"))+
  theme(axis.text.x = element_text(angle=45,hjust=1),
        legend.background = element_rect(fill = NA),
        panel.border=element_rect(fill=NA,colour="black",size=1),
        strip.background = element_blank(),
        strip.text = element_text(size = 20))+
  guides(color = guide_legend(override.aes = list(size = 0.85)))+
  guides(color = guide_legend(label.position = "left", label.hjust = 1))+
  facet_grid(rows=vars(site),cols=vars(meas))
clim_p


#jpeg("Fig.5.jpeg", width = 11, height = 20, units = 'in',res=400)
#clim_p
#dev.off()
