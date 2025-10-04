library(MASS)
library(boot)

estim_betaks<-read.csv2('betaks.csv',h=T)
varcovks<-read.csv2('varcovks.csv',h=F)

estim_beta<-read.csv2('betagg.csv',h=T)
varcov<-read.csv2('varcovgg.csv',h=F)

numbetaks<-c(4:18)
lbetaks<-length(numbetaks)
Muks<-estim_betaks[numbetaks,2]
Sigmaks<-varcovks[numbetaks,numbetaks]
n<-10000
rpks<-mvrnorm(n, Muks, Sigmaks)

numbeta<-c(9:188)
lbeta<-length(numbeta)
Mu<-estim_beta[numbeta,2]
Sigma<-varcov[numbeta,numbeta]
rp<-mvrnorm(n, Mu, Sigma)

# ADDED "O" MEANS OPERATIONAL - WITH ONLY REPRODUCTIVELY ACTIVE FEMALES

#GOLEM GRAD

fp09<-matrix(NA,nrow=n,ncol=1)
fp10<-matrix(NA,nrow=n,ncol=1)
fp11<-matrix(NA,nrow=n,ncol=1)
fp12<-matrix(NA,nrow=n,ncol=1)
fp13<-matrix(NA,nrow=n,ncol=1)
fp14<-matrix(NA,nrow=n,ncol=1)
fp15<-matrix(NA,nrow=n,ncol=1)
fp16<-matrix(NA,nrow=n,ncol=1)
fp17<-matrix(NA,nrow=n,ncol=1)
fp18<-matrix(NA,nrow=n,ncol=1)
fp19<-matrix(NA,nrow=n,ncol=1)
fp20<-matrix(NA,nrow=n,ncol=1)
fp21<-matrix(NA,nrow=n,ncol=1)
fp22<-matrix(NA,nrow=n,ncol=1)
fp23<-matrix(NA,nrow=n,ncol=1)

fb09<-matrix(NA,nrow=n,ncol=1)
fb10<-matrix(NA,nrow=n,ncol=1)
fb11<-matrix(NA,nrow=n,ncol=1)
fb12<-matrix(NA,nrow=n,ncol=1)
fb13<-matrix(NA,nrow=n,ncol=1)
fb14<-matrix(NA,nrow=n,ncol=1)
fb15<-matrix(NA,nrow=n,ncol=1)
fb16<-matrix(NA,nrow=n,ncol=1)
fb17<-matrix(NA,nrow=n,ncol=1)
fb18<-matrix(NA,nrow=n,ncol=1)
fb19<-matrix(NA,nrow=n,ncol=1)
fb20<-matrix(NA,nrow=n,ncol=1)
fb21<-matrix(NA,nrow=n,ncol=1)
fb22<-matrix(NA,nrow=n,ncol=1)
fb23<-matrix(NA,nrow=n,ncol=1)

ofp09<-matrix(NA,nrow=n,ncol=1)
ofp10<-matrix(NA,nrow=n,ncol=1)
ofp11<-matrix(NA,nrow=n,ncol=1)
ofp12<-matrix(NA,nrow=n,ncol=1)
ofp13<-matrix(NA,nrow=n,ncol=1)
ofp14<-matrix(NA,nrow=n,ncol=1)
ofp15<-matrix(NA,nrow=n,ncol=1)
ofp16<-matrix(NA,nrow=n,ncol=1)
ofp17<-matrix(NA,nrow=n,ncol=1)
ofp18<-matrix(NA,nrow=n,ncol=1)
ofp19<-matrix(NA,nrow=n,ncol=1)
ofp20<-matrix(NA,nrow=n,ncol=1)
ofp21<-matrix(NA,nrow=n,ncol=1)
ofp22<-matrix(NA,nrow=n,ncol=1)
ofp23<-matrix(NA,nrow=n,ncol=1)

ofb09<-matrix(NA,nrow=n,ncol=1)
ofb10<-matrix(NA,nrow=n,ncol=1)
ofb11<-matrix(NA,nrow=n,ncol=1)
ofb12<-matrix(NA,nrow=n,ncol=1)
ofb13<-matrix(NA,nrow=n,ncol=1)
ofb14<-matrix(NA,nrow=n,ncol=1)
ofb15<-matrix(NA,nrow=n,ncol=1)
ofb16<-matrix(NA,nrow=n,ncol=1)
ofb17<-matrix(NA,nrow=n,ncol=1)
ofb18<-matrix(NA,nrow=n,ncol=1)
ofb19<-matrix(NA,nrow=n,ncol=1)
ofb20<-matrix(NA,nrow=n,ncol=1)
ofb21<-matrix(NA,nrow=n,ncol=1)
ofb22<-matrix(NA,nrow=n,ncol=1)
ofb23<-matrix(NA,nrow=n,ncol=1)

mp09<-matrix(NA,nrow=n,ncol=1)
mp10<-matrix(NA,nrow=n,ncol=1)
mp11<-matrix(NA,nrow=n,ncol=1)
mp12<-matrix(NA,nrow=n,ncol=1)
mp13<-matrix(NA,nrow=n,ncol=1)
mp14<-matrix(NA,nrow=n,ncol=1)
mp15<-matrix(NA,nrow=n,ncol=1)
mp16<-matrix(NA,nrow=n,ncol=1)
mp17<-matrix(NA,nrow=n,ncol=1)
mp18<-matrix(NA,nrow=n,ncol=1)
mp19<-matrix(NA,nrow=n,ncol=1)
mp20<-matrix(NA,nrow=n,ncol=1)
mp21<-matrix(NA,nrow=n,ncol=1)
mp22<-matrix(NA,nrow=n,ncol=1)
mp23<-matrix(NA,nrow=n,ncol=1)

mb09<-matrix(NA,nrow=n,ncol=1)
mb10<-matrix(NA,nrow=n,ncol=1)
mb11<-matrix(NA,nrow=n,ncol=1)
mb12<-matrix(NA,nrow=n,ncol=1)
mb13<-matrix(NA,nrow=n,ncol=1)
mb14<-matrix(NA,nrow=n,ncol=1)
mb15<-matrix(NA,nrow=n,ncol=1)
mb16<-matrix(NA,nrow=n,ncol=1)
mb17<-matrix(NA,nrow=n,ncol=1)
mb18<-matrix(NA,nrow=n,ncol=1)
mb19<-matrix(NA,nrow=n,ncol=1)
mb20<-matrix(NA,nrow=n,ncol=1)
mb21<-matrix(NA,nrow=n,ncol=1)
mb22<-matrix(NA,nrow=n,ncol=1)
mb23<-matrix(NA,nrow=n,ncol=1)

psr09<-matrix(NA,nrow=n,ncol=1)
psr10<-matrix(NA,nrow=n,ncol=1)
psr11<-matrix(NA,nrow=n,ncol=1)
psr12<-matrix(NA,nrow=n,ncol=1)
psr13<-matrix(NA,nrow=n,ncol=1)
psr14<-matrix(NA,nrow=n,ncol=1)
psr15<-matrix(NA,nrow=n,ncol=1)
psr16<-matrix(NA,nrow=n,ncol=1)
psr17<-matrix(NA,nrow=n,ncol=1)
psr18<-matrix(NA,nrow=n,ncol=1)
psr19<-matrix(NA,nrow=n,ncol=1)
psr20<-matrix(NA,nrow=n,ncol=1)
psr21<-matrix(NA,nrow=n,ncol=1)
psr22<-matrix(NA,nrow=n,ncol=1)
psr23<-matrix(NA,nrow=n,ncol=1)

bsr09<-matrix(NA,nrow=n,ncol=1)
bsr10<-matrix(NA,nrow=n,ncol=1)
bsr11<-matrix(NA,nrow=n,ncol=1)
bsr12<-matrix(NA,nrow=n,ncol=1)
bsr13<-matrix(NA,nrow=n,ncol=1)
bsr14<-matrix(NA,nrow=n,ncol=1)
bsr15<-matrix(NA,nrow=n,ncol=1)
bsr16<-matrix(NA,nrow=n,ncol=1)
bsr17<-matrix(NA,nrow=n,ncol=1)
bsr18<-matrix(NA,nrow=n,ncol=1)
bsr19<-matrix(NA,nrow=n,ncol=1)
bsr20<-matrix(NA,nrow=n,ncol=1)
bsr21<-matrix(NA,nrow=n,ncol=1)
bsr22<-matrix(NA,nrow=n,ncol=1)
bsr23<-matrix(NA,nrow=n,ncol=1)

posr09<-matrix(NA,nrow=n,ncol=1)
posr10<-matrix(NA,nrow=n,ncol=1)
posr10<-matrix(NA,nrow=n,ncol=1)
posr11<-matrix(NA,nrow=n,ncol=1)
posr12<-matrix(NA,nrow=n,ncol=1)
posr13<-matrix(NA,nrow=n,ncol=1)
posr14<-matrix(NA,nrow=n,ncol=1)
posr15<-matrix(NA,nrow=n,ncol=1)
posr16<-matrix(NA,nrow=n,ncol=1)
posr17<-matrix(NA,nrow=n,ncol=1)
posr18<-matrix(NA,nrow=n,ncol=1)
posr19<-matrix(NA,nrow=n,ncol=1)
posr20<-matrix(NA,nrow=n,ncol=1)
posr21<-matrix(NA,nrow=n,ncol=1)
posr22<-matrix(NA,nrow=n,ncol=1)
posr23<-matrix(NA,nrow=n,ncol=1)

bosr09<-matrix(NA,nrow=n,ncol=1)
bosr10<-matrix(NA,nrow=n,ncol=1)
bosr11<-matrix(NA,nrow=n,ncol=1)
bosr12<-matrix(NA,nrow=n,ncol=1)
bosr13<-matrix(NA,nrow=n,ncol=1)
bosr14<-matrix(NA,nrow=n,ncol=1)
bosr15<-matrix(NA,nrow=n,ncol=1)
bosr16<-matrix(NA,nrow=n,ncol=1)
bosr17<-matrix(NA,nrow=n,ncol=1)
bosr18<-matrix(NA,nrow=n,ncol=1)
bosr19<-matrix(NA,nrow=n,ncol=1)
bosr20<-matrix(NA,nrow=n,ncol=1)
bosr21<-matrix(NA,nrow=n,ncol=1)
bosr22<-matrix(NA,nrow=n,ncol=1)
bosr23<-matrix(NA,nrow=n,ncol=1)

#KONJSKO 

mks11<-matrix(NA,nrow=n,ncol=1)
mks12<-matrix(NA,nrow=n,ncol=1)
mks13<-matrix(NA,nrow=n,ncol=1)
mks14<-matrix(NA,nrow=n,ncol=1)
mks15<-matrix(NA,nrow=n,ncol=1)
mks16<-matrix(NA,nrow=n,ncol=1)
mks17<-matrix(NA,nrow=n,ncol=1)
mks18<-matrix(NA,nrow=n,ncol=1)
mks19<-matrix(NA,nrow=n,ncol=1)
mks20<-matrix(NA,nrow=n,ncol=1)
mks21<-matrix(NA,nrow=n,ncol=1)
mks22<-matrix(NA,nrow=n,ncol=1)
mks23<-matrix(NA,nrow=n,ncol=1)

fks11<-matrix(NA,nrow=n,ncol=1)
fks12<-matrix(NA,nrow=n,ncol=1)
fks13<-matrix(NA,nrow=n,ncol=1)
fks14<-matrix(NA,nrow=n,ncol=1)
fks15<-matrix(NA,nrow=n,ncol=1)
fks16<-matrix(NA,nrow=n,ncol=1)
fks17<-matrix(NA,nrow=n,ncol=1)
fks18<-matrix(NA,nrow=n,ncol=1)
fks19<-matrix(NA,nrow=n,ncol=1)
fks20<-matrix(NA,nrow=n,ncol=1)
fks21<-matrix(NA,nrow=n,ncol=1)
fks22<-matrix(NA,nrow=n,ncol=1)
fks23<-matrix(NA,nrow=n,ncol=1)

ofks11<-matrix(NA,nrow=n,ncol=1)
ofks12<-matrix(NA,nrow=n,ncol=1)
ofks13<-matrix(NA,nrow=n,ncol=1)
ofks14<-matrix(NA,nrow=n,ncol=1)
ofks15<-matrix(NA,nrow=n,ncol=1)
ofks16<-matrix(NA,nrow=n,ncol=1)
ofks17<-matrix(NA,nrow=n,ncol=1)
ofks18<-matrix(NA,nrow=n,ncol=1)
ofks19<-matrix(NA,nrow=n,ncol=1)
ofks20<-matrix(NA,nrow=n,ncol=1)
ofks21<-matrix(NA,nrow=n,ncol=1)
ofks22<-matrix(NA,nrow=n,ncol=1)
ofks23<-matrix(NA,nrow=n,ncol=1)

srks11<-matrix(NA,nrow=n,ncol=1)
srks12<-matrix(NA,nrow=n,ncol=1)
srks13<-matrix(NA,nrow=n,ncol=1)
srks14<-matrix(NA,nrow=n,ncol=1)
srks15<-matrix(NA,nrow=n,ncol=1)
srks16<-matrix(NA,nrow=n,ncol=1)
srks17<-matrix(NA,nrow=n,ncol=1)
srks18<-matrix(NA,nrow=n,ncol=1)
srks19<-matrix(NA,nrow=n,ncol=1)
srks20<-matrix(NA,nrow=n,ncol=1)
srks21<-matrix(NA,nrow=n,ncol=1)
srks22<-matrix(NA,nrow=n,ncol=1)
srks23<-matrix(NA,nrow=n,ncol=1)

osrks11<-matrix(NA,nrow=n,ncol=1)
osrks12<-matrix(NA,nrow=n,ncol=1)
osrks13<-matrix(NA,nrow=n,ncol=1)
osrks14<-matrix(NA,nrow=n,ncol=1)
osrks15<-matrix(NA,nrow=n,ncol=1)
osrks16<-matrix(NA,nrow=n,ncol=1)
osrks17<-matrix(NA,nrow=n,ncol=1)
osrks18<-matrix(NA,nrow=n,ncol=1)
osrks19<-matrix(NA,nrow=n,ncol=1)
osrks20<-matrix(NA,nrow=n,ncol=1)
osrks21<-matrix(NA,nrow=n,ncol=1)
osrks22<-matrix(NA,nrow=n,ncol=1)
osrks23<-matrix(NA,nrow=n,ncol=1)

for(i in 1:n)
{
  fks11[i,1]<-87/(inv.logit(rpks[i,1]+rpks[i,4]))
  fks12[i,1]<-49/(inv.logit(rpks[i,1]+rpks[i,5]))
  fks13[i,1]<-36/(inv.logit(rpks[i,1]+rpks[i,6]))
  fks14[i,1]<-8/(inv.logit(rpks[i,1]+rpks[i,7]))
  fks15[i,1]<-15/(inv.logit(rpks[i,1]+rpks[i,8]))
  fks16[i,1]<-169/(inv.logit(rpks[i,1]+rpks[i,9]))
  fks17[i,1]<-176/(inv.logit(rpks[i,1]+rpks[i,10]))
  fks18[i,1]<-71/(inv.logit(rpks[i,1]+rpks[i,11]))
  fks19[i,1]<-98/(inv.logit(rpks[i,1]+rpks[i,12]))
  fks20[i,1]<-10/(inv.logit(rpks[i,1]+rpks[i,13]))
  fks21[i,1]<-79/(inv.logit(rpks[i,1]+rpks[i,14]))
  fks22[i,1]<-77/(inv.logit(rpks[i,1]+rpks[i,15]))
  fks23[i,1]<-55/(inv.logit(rpks[i,1]))
  
  ofks11[i,1]<-87/(inv.logit(rpks[i,1]+rpks[i,4]))*0.94
  ofks12[i,1]<-49/(inv.logit(rpks[i,1]+rpks[i,5]))*0.94
  ofks13[i,1]<-36/(inv.logit(rpks[i,1]+rpks[i,6]))*0.94
  ofks14[i,1]<-8/(inv.logit(rpks[i,1]+rpks[i,7]))*0.94
  ofks15[i,1]<-15/(inv.logit(rpks[i,1]+rpks[i,8]))*0.94
  ofks16[i,1]<-169/(inv.logit(rpks[i,1]+rpks[i,9]))*0.94
  ofks17[i,1]<-176/(inv.logit(rpks[i,1]+rpks[i,10]))*0.94
  ofks18[i,1]<-71/(inv.logit(rpks[i,1]+rpks[i,11]))*0.94
  ofks19[i,1]<-98/(inv.logit(rpks[i,1]+rpks[i,12]))*0.94
  ofks20[i,1]<-10/(inv.logit(rpks[i,1]+rpks[i,13]))*0.94
  ofks21[i,1]<-79/(inv.logit(rpks[i,1]+rpks[i,14]))*0.94
  ofks22[i,1]<-77/(inv.logit(rpks[i,1]+rpks[i,15]))*0.94
  ofks23[i,1]<-55/(inv.logit(rpks[i,1]))*0.94
  
  mks11[i,1]<-81/(inv.logit(rpks[i,1]+rpks[i,4]))
  mks12[i,1]<-46/(inv.logit(rpks[i,1]+rpks[i,5]))
  mks13[i,1]<-14/(inv.logit(rpks[i,1]+rpks[i,6]))
  mks14[i,1]<-10/(inv.logit(rpks[i,1]+rpks[i,7]))
  mks15[i,1]<-13/(inv.logit(rpks[i,1]+rpks[i,8]))
  mks16[i,1]<-97/(inv.logit(rpks[i,1]+rpks[i,9]))
  mks17[i,1]<-110/(inv.logit(rpks[i,1]+rpks[i,10]))
  mks18[i,1]<-32/(inv.logit(rpks[i,1]+rpks[i,11]))
  mks19[i,1]<-60/(inv.logit(rpks[i,1]+rpks[i,12]))
  mks20[i,1]<-22/(inv.logit(rpks[i,1]+rpks[i,13]))
  mks21[i,1]<-60/(inv.logit(rpks[i,1]+rpks[i,14]))
  mks22[i,1]<-24/(inv.logit(rpks[i,1]+rpks[i,15]))
  mks23[i,1]<-41/(inv.logit(rpks[i,1]))
  
  fp09[i,1]<-38/(inv.logit(rp[i,1]))
  fp10[i,1]<-26/(inv.logit(rp[i,13]))
  fp11[i,1]<-17/(inv.logit(rp[i,25]))
  fp12[i,1]<-15/(inv.logit(rp[i,37]))
  fp13[i,1]<-4/(inv.logit(rp[i,49]))
  fp14[i,1]<-10/(inv.logit(rp[i,61]))
  fp15[i,1]<-9/(inv.logit(rp[i,73]))
  fp16[i,1]<-27/(inv.logit(rp[i,85]))
  fp17[i,1]<-32/(inv.logit(rp[i,97]))
  fp18[i,1]<-19/(inv.logit(rp[i,109]))
  fp19[i,1]<-24/(inv.logit(rp[i,121]))
  fp20[i,1]<-12/(inv.logit(rp[i,133]))
  fp21[i,1]<-16/(inv.logit(rp[i,145]))
  fp22[i,1]<-8/(inv.logit(rp[i,157]))
  fp23[i,1]<-17/(inv.logit(rp[i,169]))
  
  ofp09[i,1]<-38/(inv.logit(rp[i,1]))*0.15
  ofp10[i,1]<-26/(inv.logit(rp[i,13]))*0.15
  ofp11[i,1]<-17/(inv.logit(rp[i,25]))*0.15
  ofp12[i,1]<-15/(inv.logit(rp[i,37]))*0.15
  ofp13[i,1]<-4/(inv.logit(rp[i,49]))*0.15
  ofp14[i,1]<-10/(inv.logit(rp[i,61]))*0.15
  ofp15[i,1]<-9/(inv.logit(rp[i,73]))*0.15
  ofp16[i,1]<-27/(inv.logit(rp[i,85]))*0.15
  ofp17[i,1]<-32/(inv.logit(rp[i,97]))*0.15
  ofp18[i,1]<-19/(inv.logit(rp[i,109]))*0.15
  ofp19[i,1]<-24/(inv.logit(rp[i,121]))*0.15
  ofp20[i,1]<-12/(inv.logit(rp[i,133]))*0.15
  ofp21[i,1]<-16/(inv.logit(rp[i,145]))*0.15
  ofp22[i,1]<-8/(inv.logit(rp[i,157]))*0.15
  ofp23[i,1]<-17/(inv.logit(rp[i,169]))*0.15
  
  fb09[i,1]<-20/(inv.logit(rp[i,3]))
  fb10[i,1]<-9/(inv.logit(rp[i,15]))
  fb11[i,1]<-15/(inv.logit(rp[i,27]))
  fb12[i,1]<-9/(inv.logit(rp[i,39]))
  fb13[i,1]<-8/(inv.logit(rp[i,51]))
  fb14[i,1]<-7/(inv.logit(rp[i,63]))
  fb15[i,1]<-5/(inv.logit(rp[i,75]))
  fb16[i,1]<-19/(inv.logit(rp[i,87]))
  fb17[i,1]<-8/(inv.logit(rp[i,99]))
  fb18[i,1]<-14/(inv.logit(rp[i,111]))
  fb19[i,1]<-15/(inv.logit(rp[i,123]))
  fb20[i,1]<-11/(inv.logit(rp[i,135]))
  fb21[i,1]<-9/(inv.logit(rp[i,147]))
  fb22[i,1]<-4/(inv.logit(rp[i,159]))
  fb23[i,1]<-15/(inv.logit(rp[i,171]))
  
  ofb09[i,1]<-20/(inv.logit(rp[i,3]))*0.71
  ofb10[i,1]<-9/(inv.logit(rp[i,15]))*0.71
  ofb11[i,1]<-15/(inv.logit(rp[i,27]))*0.71
  ofb12[i,1]<-9/(inv.logit(rp[i,39]))*0.71
  ofb13[i,1]<-8/(inv.logit(rp[i,51]))*0.71
  ofb14[i,1]<-7/(inv.logit(rp[i,63]))*0.71
  ofb15[i,1]<-5/(inv.logit(rp[i,75]))*0.71
  ofb16[i,1]<-19/(inv.logit(rp[i,87]))*0.71
  ofb17[i,1]<-8/(inv.logit(rp[i,99]))*0.71
  ofb18[i,1]<-14/(inv.logit(rp[i,111]))*0.71
  ofb19[i,1]<-15/(inv.logit(rp[i,123]))*0.71
  ofb20[i,1]<-11/(inv.logit(rp[i,135]))*0.71
  ofb21[i,1]<-9/(inv.logit(rp[i,147]))*0.71
  ofb22[i,1]<-4/(inv.logit(rp[i,159]))*0.71
  ofb23[i,1]<-15/(inv.logit(rp[i,171]))*0.71
  
  mp09[i,1]<-677/(inv.logit(rp[i,7]))
  mp10[i,1]<-797/(inv.logit(rp[i,19]))
  mp11[i,1]<-698/(inv.logit(rp[i,31]))
  mp12[i,1]<-479/(inv.logit(rp[i,43]))
  mp13[i,1]<-265/(inv.logit(rp[i,55]))
  mp14[i,1]<-180/(inv.logit(rp[i,67]))
  mp15[i,1]<-318/(inv.logit(rp[i,79]))
  mp16[i,1]<-739/(inv.logit(rp[i,91]))
  mp17[i,1]<-680/(inv.logit(rp[i,103]))
  mp18[i,1]<-502/(inv.logit(rp[i,115]))
  mp19[i,1]<-532/(inv.logit(rp[i,127]))
  mp20[i,1]<-237/(inv.logit(rp[i,139]))
  mp21[i,1]<-347/(inv.logit(rp[i,151]))
  mp22[i,1]<-222/(inv.logit(rp[i,163]))
  mp23[i,1]<-358/(inv.logit(rp[i,175]))
  
  mb09[i,1]<-47/(inv.logit(rp[i,9]))
  mb10[i,1]<-46/(inv.logit(rp[i,21]))
  mb11[i,1]<-45/(inv.logit(rp[i,33]))
  mb12[i,1]<-33/(inv.logit(rp[i,45]))
  mb13[i,1]<-19/(inv.logit(rp[i,57]))
  mb14[i,1]<-28/(inv.logit(rp[i,69]))
  mb15[i,1]<-19/(inv.logit(rp[i,81]))
  mb16[i,1]<-64/(inv.logit(rp[i,93]))
  mb17[i,1]<-30/(inv.logit(rp[i,105]))
  mb18[i,1]<-34/(inv.logit(rp[i,117]))
  mb19[i,1]<-41/(inv.logit(rp[i,129]))
  mb20[i,1]<-29/(inv.logit(rp[i,141]))
  mb21[i,1]<-22/(inv.logit(rp[i,153]))
  mb22[i,1]<-8/(inv.logit(rp[i,165]))
  mb23[i,1]<-22/(inv.logit(rp[i,177]))

#CALCULATE POP AND SR FOR COHORTS FROM MATRICES

  
  psr09[i,1]<-mp09[i,1]/fp09[i,1]
  psr10[i,1]<-mp10[i,1]/fp10[i,1]
  psr11[i,1]<-mp11[i,1]/fp11[i,1]
  psr12[i,1]<-mp12[i,1]/fp12[i,1]
  psr13[i,1]<-mp13[i,1]/fp13[i,1]
  psr14[i,1]<-mp14[i,1]/fp14[i,1]
  psr15[i,1]<-mp15[i,1]/fp15[i,1]
  psr16[i,1]<-mp16[i,1]/fp16[i,1]
  psr17[i,1]<-mp17[i,1]/fp17[i,1]
  psr18[i,1]<-mp18[i,1]/fp18[i,1]
  psr19[i,1]<-mp19[i,1]/fp19[i,1]
  psr20[i,1]<-mp20[i,1]/fp20[i,1]
  psr21[i,1]<-mp21[i,1]/fp21[i,1]
  psr22[i,1]<-mp22[i,1]/fp22[i,1]
  psr23[i,1]<-mp23[i,1]/fp23[i,1]
  
  bsr09[i,1]<-mb09[i,1]/fb09[i,1]
  bsr10[i,1]<-mb10[i,1]/fb10[i,1]
  bsr11[i,1]<-mb11[i,1]/fb11[i,1]
  bsr12[i,1]<-mb12[i,1]/fb12[i,1]
  bsr13[i,1]<-mb13[i,1]/fb13[i,1]
  bsr14[i,1]<-mb14[i,1]/fb14[i,1]
  bsr15[i,1]<-mb15[i,1]/fb15[i,1]
  bsr16[i,1]<-mb16[i,1]/fb16[i,1]
  bsr17[i,1]<-mb17[i,1]/fb17[i,1]
  bsr18[i,1]<-mb18[i,1]/fb18[i,1]
  bsr19[i,1]<-mb19[i,1]/fb19[i,1]
  bsr20[i,1]<-mb20[i,1]/fb20[i,1]
  bsr21[i,1]<-mb21[i,1]/fb21[i,1]
  bsr22[i,1]<-mb22[i,1]/fb22[i,1]
  bsr23[i,1]<-mb23[i,1]/fb23[i,1]
  
  posr09[i,1]<-mp09[i,1]/ofp09[i,1]
  posr10[i,1]<-mp10[i,1]/ofp10[i,1]
  posr11[i,1]<-mp11[i,1]/ofp11[i,1]
  posr12[i,1]<-mp12[i,1]/ofp12[i,1]
  posr13[i,1]<-mp13[i,1]/ofp13[i,1]
  posr14[i,1]<-mp14[i,1]/ofp14[i,1]
  posr15[i,1]<-mp15[i,1]/ofp15[i,1]
  posr16[i,1]<-mp16[i,1]/ofp16[i,1]
  posr17[i,1]<-mp17[i,1]/ofp17[i,1]
  posr18[i,1]<-mp18[i,1]/ofp18[i,1]
  posr19[i,1]<-mp19[i,1]/ofp19[i,1]
  posr20[i,1]<-mp20[i,1]/ofp20[i,1]
  posr21[i,1]<-mp21[i,1]/ofp21[i,1]
  posr22[i,1]<-mp22[i,1]/ofp22[i,1]
  posr23[i,1]<-mp23[i,1]/ofp23[i,1]
  
  bosr09[i,1]<-mb09[i,1]/ofb09[i,1]
  bosr10[i,1]<-mb10[i,1]/ofb10[i,1]
  bosr11[i,1]<-mb11[i,1]/ofb11[i,1]
  bosr12[i,1]<-mb12[i,1]/ofb12[i,1]
  bosr13[i,1]<-mb13[i,1]/ofb13[i,1]
  bosr14[i,1]<-mb14[i,1]/ofb14[i,1]
  bosr15[i,1]<-mb15[i,1]/ofb15[i,1]
  bosr16[i,1]<-mb16[i,1]/ofb16[i,1]
  bosr17[i,1]<-mb17[i,1]/ofb17[i,1]
  bosr18[i,1]<-mb18[i,1]/ofb18[i,1]
  bosr19[i,1]<-mb19[i,1]/ofb19[i,1]
  bosr20[i,1]<-mb20[i,1]/ofb20[i,1]
  bosr21[i,1]<-mb21[i,1]/ofb21[i,1]
  bosr22[i,1]<-mb22[i,1]/ofb22[i,1]
  bosr23[i,1]<-mb23[i,1]/ofb23[i,1]
  
  srks11[i,1]<-mks11[i,1]/fks11[i,1]
  srks12[i,1]<-mks12[i,1]/fks12[i,1]
  srks13[i,1]<-mks13[i,1]/fks13[i,1]
  srks14[i,1]<-mks14[i,1]/fks14[i,1]
  srks15[i,1]<-mks15[i,1]/fks15[i,1]
  srks16[i,1]<-mks16[i,1]/fks16[i,1]
  srks17[i,1]<-mks17[i,1]/fks17[i,1]
  srks18[i,1]<-mks18[i,1]/fks18[i,1]
  srks19[i,1]<-mks19[i,1]/fks19[i,1]
  srks20[i,1]<-mks20[i,1]/fks20[i,1]
  srks21[i,1]<-mks21[i,1]/fks21[i,1]
  srks22[i,1]<-mks22[i,1]/fks22[i,1]
  srks23[i,1]<-mks23[i,1]/fks23[i,1]
  
  osrks11[i,1]<-mks11[i,1]/ofks11[i,1]
  osrks12[i,1]<-mks12[i,1]/ofks12[i,1]
  osrks13[i,1]<-mks13[i,1]/ofks13[i,1]
  osrks14[i,1]<-mks14[i,1]/ofks14[i,1]
  osrks15[i,1]<-mks15[i,1]/ofks15[i,1]
  osrks16[i,1]<-mks16[i,1]/ofks16[i,1]
  osrks17[i,1]<-mks17[i,1]/ofks17[i,1]
  osrks18[i,1]<-mks18[i,1]/ofks18[i,1]
  osrks19[i,1]<-mks19[i,1]/ofks19[i,1]
  osrks20[i,1]<-mks20[i,1]/ofks20[i,1]
  osrks21[i,1]<-mks21[i,1]/ofks21[i,1]
  osrks22[i,1]<-mks22[i,1]/ofks22[i,1]
  osrks23[i,1]<-mks23[i,1]/ofks23[i,1]
}
 
 #### GET QUANTILES FOR ALL MATRICES 
#GOLEM GRAD

# calculate quantiles for each matrix
qfp1 <- quantile(fp09,probs=c(0.025,0.5,0.975))
qfp2 <- quantile(fp10,probs=c(0.025,0.5,0.975))
qfp3 <- quantile(fp11,probs=c(0.025,0.5,0.975))
qfp4 <- quantile(fp12,probs=c(0.025,0.5,0.975))
qfp5 <- quantile(fp13,probs=c(0.025,0.5,0.975))
qfp6 <- quantile(fp14,probs=c(0.025,0.5,0.975))
qfp7 <- quantile(fp15,probs=c(0.025,0.5,0.975))
qfp8 <- quantile(fp16,probs=c(0.025,0.5,0.975))
qfp9 <- quantile(fp17,probs=c(0.025,0.5,0.975))
qfp10 <- quantile(fp18,probs=c(0.025,0.5,0.975))
qfp11 <- quantile(fp19,probs=c(0.025,0.5,0.975))
qfp12 <- quantile(fp20,probs=c(0.025,0.5,0.975))
qfp13 <- quantile(fp21,probs=c(0.025,0.5,0.975))
qfp14 <- quantile(fp22,probs=c(0.025,0.5,0.975))
qfp15 <- quantile(fp23,probs=c(0.025,0.5,0.975))

# Create a data frame with quantiles

fp <- data.frame(
    Variable = as.numeric(c("2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016", "2017", "2018", "2019", "2020", "2021", "2022", "2023")),
     Q_0.025 = as.numeric(c(qfp1[1], qfp2[1], qfp3[1], qfp4[1], qfp5[1], qfp6[1], qfp7[1], qfp8[1], qfp9[1], qfp10[1], qfp11[1], qfp12[1], qfp13[1], qfp14[1], qfp15[1])),
     Q_0.5   = as.numeric(c(qfp1[2], qfp2[2], qfp3[2], qfp4[2], qfp5[2], qfp6[2], qfp7[2], qfp8[2], qfp9[2], qfp10[2], qfp11[2], qfp12[2], qfp13[2], qfp14[2], qfp15[2])),
     Q_0.975 = as.numeric(c(qfp1[3], qfp2[3], qfp3[3], qfp4[3], qfp5[3], qfp6[3], qfp7[3], qfp8[3], qfp9[3], qfp10[3], qfp11[3], qfp12[3], qfp13[3], qfp14[3], qfp15[3]))
 )
 
# Write the quantiles data frame to a CSV file

write.csv(fp, file="qfp.csv", row.names = FALSE)


#Do for all groups

qfb1 <- quantile(fb09,probs=c(0.025,0.5,0.975))
qfb2 <- quantile(fb10,probs=c(0.025,0.5,0.975))
qfb3 <- quantile(fb11,probs=c(0.025,0.5,0.975))
qfb4 <- quantile(fb12,probs=c(0.025,0.5,0.975))
qfb5 <- quantile(fb13,probs=c(0.025,0.5,0.975))
qfb6 <- quantile(fb14,probs=c(0.025,0.5,0.975))
qfb7 <- quantile(fb15,probs=c(0.025,0.5,0.975))
qfb8 <- quantile(fb16,probs=c(0.025,0.5,0.975))
qfb9 <- quantile(fb17,probs=c(0.025,0.5,0.975))
qfb10 <- quantile(fb18,probs=c(0.025,0.5,0.975))
qfb11 <- quantile(fb19,probs=c(0.025,0.5,0.975))
qfb12 <- quantile(fb20,probs=c(0.025,0.5,0.975))
qfb13 <- quantile(fb21,probs=c(0.025,0.5,0.975))
qfb14 <- quantile(fb22,probs=c(0.025,0.5,0.975))
qfb15 <- quantile(fb23,probs=c(0.025,0.5,0.975))

fb <- data.frame(
    Variable = as.numeric(c("2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016", "2017", "2018", "2019", "2020", "2021", "2022", "2023")),
     Q_0.025 = as.numeric(c(qfb1[1], qfb2[1], qfb3[1], qfb4[1], qfb5[1], qfb6[1], qfb7[1], qfb8[1], qfb9[1], qfb10[1], qfb11[1], qfb12[1], qfb13[1], qfb14[1], qfb15[1])),
     Q_0.5   = as.numeric(c(qfb1[2], qfb2[2], qfb3[2], qfb4[2], qfb5[2], qfb6[2], qfb7[2], qfb8[2], qfb9[2], qfb10[2], qfb11[2], qfb12[2], qfb13[2], qfb14[2], qfb15[2])),
     Q_0.975 = as.numeric(c(qfb1[3], qfb2[3], qfb3[3], qfb4[3], qfb5[3], qfb6[3], qfb7[3], qfb8[3], qfb9[3], qfb10[3], qfb11[3], qfb12[3], qfb13[3], qfb14[3], qfb15[3]))
 )

write.csv(fb, file="qfb.csv", row.names = FALSE)

# calculate quantiles for each matrix
qofp1 <- quantile(ofp09,probs=c(0.025,0.5,0.975))
qofp2 <- quantile(ofp10,probs=c(0.025,0.5,0.975))
qofp3 <- quantile(ofp11,probs=c(0.025,0.5,0.975))
qofp4 <- quantile(ofp12,probs=c(0.025,0.5,0.975))
qofp5 <- quantile(ofp13,probs=c(0.025,0.5,0.975))
qofp6 <- quantile(ofp14,probs=c(0.025,0.5,0.975))
qofp7 <- quantile(ofp15,probs=c(0.025,0.5,0.975))
qofp8 <- quantile(ofp16,probs=c(0.025,0.5,0.975))
qofp9 <- quantile(ofp17,probs=c(0.025,0.5,0.975))
qofp10 <- quantile(ofp18,probs=c(0.025,0.5,0.975))
qofp11 <- quantile(ofp19,probs=c(0.025,0.5,0.975))
qofp12 <- quantile(ofp20,probs=c(0.025,0.5,0.975))
qofp13 <- quantile(ofp21,probs=c(0.025,0.5,0.975))
qofp14 <- quantile(ofp22,probs=c(0.025,0.5,0.975))
qofp15 <- quantile(ofp23,probs=c(0.025,0.5,0.975))

# Create a data frame with quantiles

ofp <- data.frame(
    Variable = as.numeric(c("2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016", "2017", "2018", "2019", "2020", "2021", "2022", "2023")),
     Q_0.025 = as.numeric(c(qofp1[1], qofp2[1], qofp3[1], qofp4[1], qofp5[1], qofp6[1], qofp7[1], qofp8[1], qofp9[1], qofp10[1], qofp11[1], qofp12[1], qofp13[1], qofp14[1], qofp15[1])),
     Q_0.5   = as.numeric(c(qofp1[2], qofp2[2], qofp3[2], qofp4[2], qofp5[2], qofp6[2], qofp7[2], qofp8[2], qofp9[2], qofp10[2], qofp11[2], qofp12[2], qofp13[2], qofp14[2], qofp15[2])),
     Q_0.975 = as.numeric(c(qofp1[3], qofp2[3], qofp3[3], qofp4[3], qofp5[3], qofp6[3], qofp7[3], qofp8[3], qofp9[3], qofp10[3], qofp11[3], qofp12[3], qofp13[3], qofp14[3], qofp15[3]))
 )
 
# Write the quantiles data frame to a CSV file

write.csv(ofp, file="qofp.csv", row.names = FALSE)


#Do for all groups

qofb1 <- quantile(ofb09,probs=c(0.025,0.5,0.975))
qofb2 <- quantile(ofb10,probs=c(0.025,0.5,0.975))
qofb3 <- quantile(ofb11,probs=c(0.025,0.5,0.975))
qofb4 <- quantile(ofb12,probs=c(0.025,0.5,0.975))
qofb5 <- quantile(ofb13,probs=c(0.025,0.5,0.975))
qofb6 <- quantile(ofb14,probs=c(0.025,0.5,0.975))
qofb7 <- quantile(ofb15,probs=c(0.025,0.5,0.975))
qofb8 <- quantile(ofb16,probs=c(0.025,0.5,0.975))
qofb9 <- quantile(ofb17,probs=c(0.025,0.5,0.975))
qofb10 <- quantile(ofb18,probs=c(0.025,0.5,0.975))
qofb11 <- quantile(ofb19,probs=c(0.025,0.5,0.975))
qofb12 <- quantile(ofb20,probs=c(0.025,0.5,0.975))
qofb13 <- quantile(ofb21,probs=c(0.025,0.5,0.975))
qofb14 <- quantile(ofb22,probs=c(0.025,0.5,0.975))
qofb15 <- quantile(ofb23,probs=c(0.025,0.5,0.975))

ofb <- data.frame(
    Variable = as.numeric(c("2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016", "2017", "2018", "2019", "2020", "2021", "2022", "2023")),
     Q_0.025 = as.numeric(c(qofb1[1], qofb2[1], qofb3[1], qofb4[1], qofb5[1], qofb6[1], qofb7[1], qofb8[1], qofb9[1], qofb10[1], qofb11[1], qofb12[1], qofb13[1], qofb14[1], qofb15[1])),
     Q_0.5   = as.numeric(c(qofb1[2], qofb2[2], qofb3[2], qofb4[2], qofb5[2], qofb6[2], qofb7[2], qofb8[2], qofb9[2], qofb10[2], qofb11[2], qofb12[2], qofb13[2], qofb14[2], qofb15[2])),
     Q_0.975 = as.numeric(c(qofb1[3], qofb2[3], qofb3[3], qofb4[3], qofb5[3], qofb6[3], qofb7[3], qofb8[3], qofb9[3], qofb10[3], qofb11[3], qofb12[3], qofb13[3], qofb14[3], qofb15[3]))
 )

write.csv(ofb, file="qofb.csv", row.names = FALSE)

qmp1 <- quantile(mp09,probs=c(0.025,0.5,0.975))
qmp2 <- quantile(mp10,probs=c(0.025,0.5,0.975))
qmp3 <- quantile(mp11,probs=c(0.025,0.5,0.975))
qmp4 <- quantile(mp12,probs=c(0.025,0.5,0.975))
qmp5 <- quantile(mp13,probs=c(0.025,0.5,0.975))
qmp6 <- quantile(mp14,probs=c(0.025,0.5,0.975))
qmp7 <- quantile(mp15,probs=c(0.025,0.5,0.975))
qmp8 <- quantile(mp16,probs=c(0.025,0.5,0.975))
qmp9 <- quantile(mp17,probs=c(0.025,0.5,0.975))
qmp10 <- quantile(mp18,probs=c(0.025,0.5,0.975))
qmp11 <- quantile(mp19,probs=c(0.025,0.5,0.975))
qmp12 <- quantile(mp20,probs=c(0.025,0.5,0.975))
qmp13 <- quantile(mp21,probs=c(0.025,0.5,0.975))
qmp14 <- quantile(mp22,probs=c(0.025,0.5,0.975))
qmp15 <- quantile(mp23,probs=c(0.025,0.5,0.975))

mp <- data.frame(
    Variable = as.numeric(c("2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016", "2017", "2018", "2019", "2020", "2021", "2022", "2023")),
     Q_0.025 = as.numeric(c(qmp1[1], qmp2[1], qmp3[1], qmp4[1], qmp5[1], qmp6[1], qmp7[1], qmp8[1], qmp9[1], qmp10[1], qmp11[1], qmp12[1], qmp13[1], qmp14[1], qmp15[1])),
     Q_0.5   = as.numeric(c(qmp1[2], qmp2[2], qmp3[2], qmp4[2], qmp5[2], qmp6[2], qmp7[2], qmp8[2], qmp9[2], qmp10[2], qmp11[2], qmp12[2], qmp13[2], qmp14[2], qmp15[2])),
     Q_0.975 = as.numeric(c(qmp1[3], qmp2[3], qmp3[3], qmp4[3], qmp5[3], qmp6[3], qmp7[3], qmp8[3], qmp9[3], qmp10[3], qmp11[3], qmp12[3], qmp13[3], qmp14[3], qmp15[3]))
 )

write.csv(mp, file="qmp.csv", row.names = FALSE)

qmb1 <- quantile(mb09,probs=c(0.025,0.5,0.975))
qmb2 <- quantile(mb10,probs=c(0.025,0.5,0.975))
qmb3 <- quantile(mb11,probs=c(0.025,0.5,0.975))
qmb4 <- quantile(mb12,probs=c(0.025,0.5,0.975))
qmb5 <- quantile(mb13,probs=c(0.025,0.5,0.975))
qmb6 <- quantile(mb14,probs=c(0.025,0.5,0.975))
qmb7 <- quantile(mb15,probs=c(0.025,0.5,0.975))
qmb8 <- quantile(mb16,probs=c(0.025,0.5,0.975))
qmb9 <- quantile(mb17,probs=c(0.025,0.5,0.975))
qmb10 <- quantile(mb18,probs=c(0.025,0.5,0.975))
qmb11 <- quantile(mb19,probs=c(0.025,0.5,0.975))
qmb12 <- quantile(mb20,probs=c(0.025,0.5,0.975))
qmb13 <- quantile(mb21,probs=c(0.025,0.5,0.975))
qmb14 <- quantile(mb22,probs=c(0.025,0.5,0.975))
qmb15 <- quantile(mb23,probs=c(0.025,0.5,0.975))

mb <- data.frame(
    Variable = as.numeric(c("2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016", "2017", "2018", "2019", "2020", "2021", "2022", "2023")),
     Q_0.025 = as.numeric(c(qmb1[1], qmb2[1], qmb3[1], qmb4[1], qmb5[1], qmb6[1], qmb7[1], qmb8[1], qmb9[1], qmb10[1], qmb11[1], qmb12[1], qmb13[1], qmb14[1], qmb15[1])),
     Q_0.5   = as.numeric(c(qmb1[2], qmb2[2], qmb3[2], qmb4[2], qmb5[2], qmb6[2], qmb7[2], qmb8[2], qmb9[2], qmb10[2], qmb11[2], qmb12[2], qmb13[2], qmb14[2], qmb15[2])),
     Q_0.975 = as.numeric(c(qmb1[3], qmb2[3], qmb3[3], qmb4[3], qmb5[3], qmb6[3], qmb7[3], qmb8[3], qmb9[3], qmb10[3], qmb11[3], qmb12[3], qmb13[3], qmb14[3], qmb15[3]))
 )

write.csv(mb, file="qmb.csv", row.names = FALSE)

qpsr1 <- quantile(psr09,probs=c(0.025,0.5,0.975))
qpsr2 <- quantile(psr10,probs=c(0.025,0.5,0.975))
qpsr3 <- quantile(psr11,probs=c(0.025,0.5,0.975))
qpsr4 <- quantile(psr12,probs=c(0.025,0.5,0.975))
qpsr5 <- quantile(psr13,probs=c(0.025,0.5,0.975))
qpsr6 <- quantile(psr14,probs=c(0.025,0.5,0.975))
qpsr7 <- quantile(psr15,probs=c(0.025,0.5,0.975))
qpsr8 <- quantile(psr16,probs=c(0.025,0.5,0.975))
qpsr9 <- quantile(psr17,probs=c(0.025,0.5,0.975))
qpsr10 <- quantile(psr18,probs=c(0.025,0.5,0.975))
qpsr11 <- quantile(psr19,probs=c(0.025,0.5,0.975))
qpsr12 <- quantile(psr20,probs=c(0.025,0.5,0.975))
qpsr13 <- quantile(psr21,probs=c(0.025,0.5,0.975))
qpsr14 <- quantile(psr22,probs=c(0.025,0.5,0.975))
qpsr15 <- quantile(psr23,probs=c(0.025,0.5,0.975))

psr <- data.frame(
    Variable = as.numeric(c("2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016", "2017", "2018", "2019", "2020", "2021", "2022", "2023")),
     Q_0.025 = as.numeric(c(qpsr1[1], qpsr2[1], qpsr3[1], qpsr4[1], qpsr5[1], qpsr6[1], qpsr7[1], qpsr8[1], qpsr9[1], qpsr10[1], qpsr11[1], qpsr12[1], qpsr13[1], qpsr14[1], qpsr15[1])),
     Q_0.5   = as.numeric(c(qpsr1[2], qpsr2[2], qpsr3[2], qpsr4[2], qpsr5[2], qpsr6[2], qpsr7[2], qpsr8[2], qpsr9[2], qpsr10[2], qpsr11[2], qpsr12[2], qpsr13[2], qpsr14[2], qpsr15[2])),
     Q_0.975 = as.numeric(c(qpsr1[3], qpsr2[3], qpsr3[3], qpsr4[3], qpsr5[3], qpsr6[3], qpsr7[3], qpsr8[3], qpsr9[3], qpsr10[3], qpsr11[3], qpsr12[3], qpsr13[3], qpsr14[3], qpsr15[3]))
 )

write.csv(psr, file="qpsr.csv", row.names = FALSE)

qbsr1 <- quantile(bsr09,probs=c(0.025,0.5,0.975))
qbsr2 <- quantile(bsr10,probs=c(0.025,0.5,0.975))
qbsr3 <- quantile(bsr11,probs=c(0.025,0.5,0.975))
qbsr4 <- quantile(bsr12,probs=c(0.025,0.5,0.975))
qbsr5 <- quantile(bsr13,probs=c(0.025,0.5,0.975))
qbsr6 <- quantile(bsr14,probs=c(0.025,0.5,0.975))
qbsr7 <- quantile(bsr15,probs=c(0.025,0.5,0.975))
qbsr8 <- quantile(bsr16,probs=c(0.025,0.5,0.975))
qbsr9 <- quantile(bsr17,probs=c(0.025,0.5,0.975))
qbsr10 <- quantile(bsr18,probs=c(0.025,0.5,0.975))
qbsr11 <- quantile(bsr19,probs=c(0.025,0.5,0.975))
qbsr12 <- quantile(bsr20,probs=c(0.025,0.5,0.975))
qbsr13 <- quantile(bsr21,probs=c(0.025,0.5,0.975))
qbsr14 <- quantile(bsr22,probs=c(0.025,0.5,0.975))
qbsr15 <- quantile(bsr23,probs=c(0.025,0.5,0.975))

bsr <- data.frame(
    Variable = as.numeric(c("2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016", "2017", "2018", "2019", "2020", "2021", "2022", "2023")),
     Q_0.025 = as.numeric(c(qbsr1[1], qbsr2[1], qbsr3[1], qbsr4[1], qbsr5[1], qbsr6[1], qbsr7[1], qbsr8[1], qbsr9[1], qbsr10[1], qbsr11[1], qbsr12[1], qbsr13[1], qbsr14[1], qbsr15[1])),
     Q_0.5   = as.numeric(c(qbsr1[2], qbsr2[2], qbsr3[2], qbsr4[2], qbsr5[2], qbsr6[2], qbsr7[2], qbsr8[2], qbsr9[2], qbsr10[2], qbsr11[2], qbsr12[2], qbsr13[2], qbsr14[2], qbsr15[2])),
     Q_0.975 = as.numeric(c(qbsr1[3], qbsr2[3], qbsr3[3], qbsr4[3], qbsr5[3], qbsr6[3], qbsr7[3], qbsr8[3], qbsr9[3], qbsr10[3], qbsr11[3], qbsr12[3], qbsr13[3], qbsr14[3], qbsr15[3]))
 )

write.csv(bsr, file="qbsr.csv", row.names = FALSE)
qposr1 <- quantile(posr09,probs=c(0.025,0.5,0.975))
qposr2 <- quantile(posr10,probs=c(0.025,0.5,0.975))
qposr3 <- quantile(posr11,probs=c(0.025,0.5,0.975))
qposr4 <- quantile(posr12,probs=c(0.025,0.5,0.975))
qposr5 <- quantile(posr13,probs=c(0.025,0.5,0.975))
qposr6 <- quantile(posr14,probs=c(0.025,0.5,0.975))
qposr7 <- quantile(posr15,probs=c(0.025,0.5,0.975))
qposr8 <- quantile(posr16,probs=c(0.025,0.5,0.975))
qposr9 <- quantile(posr17,probs=c(0.025,0.5,0.975))
qposr10 <- quantile(posr18,probs=c(0.025,0.5,0.975))
qposr11 <- quantile(posr19,probs=c(0.025,0.5,0.975))
qposr12 <- quantile(posr20,probs=c(0.025,0.5,0.975))
qposr13 <- quantile(posr21,probs=c(0.025,0.5,0.975))
qposr14 <- quantile(posr22,probs=c(0.025,0.5,0.975))
qposr15 <- quantile(posr23,probs=c(0.025,0.5,0.975))

posr <- data.frame(
    Variable = as.numeric(c("2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016", "2017", "2018", "2019", "2020", "2021", "2022", "2023")),
     Q_0.025 = as.numeric(c(qposr1[1], qposr2[1], qposr3[1], qposr4[1], qposr5[1], qposr6[1], qposr7[1], qposr8[1], qposr9[1], qposr10[1], qposr11[1], qposr12[1], qposr13[1], qposr14[1], qposr15[1])),
     Q_0.5   = as.numeric(c(qposr1[2], qposr2[2], qposr3[2], qposr4[2], qposr5[2], qposr6[2], qposr7[2], qposr8[2], qposr9[2], qposr10[2], qposr11[2], qposr12[2], qposr13[2], qposr14[2], qposr15[2])),
     Q_0.975 = as.numeric(c(qposr1[3], qposr2[3], qposr3[3], qposr4[3], qposr5[3], qposr6[3], qposr7[3], qposr8[3], qposr9[3], qposr10[3], qposr11[3], qposr12[3], qposr13[3], qposr14[3], qposr15[3]))
 )

write.csv(posr, file="qposr.csv", row.names = FALSE)

qbosr1 <- quantile(bosr09,probs=c(0.025,0.5,0.975))
qbosr2 <- quantile(bosr10,probs=c(0.025,0.5,0.975))
qbosr3 <- quantile(bosr11,probs=c(0.025,0.5,0.975))
qbosr4 <- quantile(bosr12,probs=c(0.025,0.5,0.975))
qbosr5 <- quantile(bosr13,probs=c(0.025,0.5,0.975))
qbosr6 <- quantile(bosr14,probs=c(0.025,0.5,0.975))
qbosr7 <- quantile(bosr15,probs=c(0.025,0.5,0.975))
qbosr8 <- quantile(bosr16,probs=c(0.025,0.5,0.975))
qbosr9 <- quantile(bosr17,probs=c(0.025,0.5,0.975))
qbosr10 <- quantile(bosr18,probs=c(0.025,0.5,0.975))
qbosr11 <- quantile(bosr19,probs=c(0.025,0.5,0.975))
qbosr12 <- quantile(bosr20,probs=c(0.025,0.5,0.975))
qbosr13 <- quantile(bosr21,probs=c(0.025,0.5,0.975))
qbosr14 <- quantile(bosr22,probs=c(0.025,0.5,0.975))
qbosr15 <- quantile(bosr23,probs=c(0.025,0.5,0.975))

bosr <- data.frame(
    Variable = as.numeric(c("2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016", "2017", "2018", "2019", "2020", "2021", "2022", "2023")),
     Q_0.025 = as.numeric(c(qbosr1[1], qbosr2[1], qbosr3[1], qbosr4[1], qbosr5[1], qbosr6[1], qbosr7[1], qbosr8[1], qbosr9[1], qbosr10[1], qbosr11[1], qbosr12[1], qbosr13[1], qbosr14[1], qbosr15[1])),
     Q_0.5   = as.numeric(c(qbosr1[2], qbosr2[2], qbosr3[2], qbosr4[2], qbosr5[2], qbosr6[2], qbosr7[2], qbosr8[2], qbosr9[2], qbosr10[2], qbosr11[2], qbosr12[2], qbosr13[2], qbosr14[2], qbosr15[2])),
     Q_0.975 = as.numeric(c(qbosr1[3], qbosr2[3], qbosr3[3], qbosr4[3], qbosr5[3], qbosr6[3], qbosr7[3], qbosr8[3], qbosr9[3], qbosr10[3], qbosr11[3], qbosr12[3], qbosr13[3], qbosr14[3], qbosr15[3]))
 )

write.csv(bosr, file="qbosr.csv", row.names = FALSE)

#KONJSKO

qmks3 <- quantile(mks11,probs=c(0.025,0.5,0.975))
qmks4 <- quantile(mks12,probs=c(0.025,0.5,0.975))
qmks5 <- quantile(mks13,probs=c(0.025,0.5,0.975))
qmks6 <- quantile(mks14,probs=c(0.025,0.5,0.975))
qmks7 <- quantile(mks15,probs=c(0.025,0.5,0.975))
qmks8 <- quantile(mks16,probs=c(0.025,0.5,0.975))
qmks9 <- quantile(mks17,probs=c(0.025,0.5,0.975))
qmks10 <- quantile(mks18,probs=c(0.025,0.5,0.975))
qmks11 <- quantile(mks19,probs=c(0.025,0.5,0.975))
qmks12 <- quantile(mks20,probs=c(0.025,0.5,0.975))
qmks13 <- quantile(mks21,probs=c(0.025,0.5,0.975))
qmks14 <- quantile(mks22,probs=c(0.025,0.5,0.975))
qmks15 <- quantile(mks23,probs=c(0.025,0.5,0.975))

mks <- data.frame(
    Variable = as.numeric(c("2011", "2012", "2013", "2014", "2015", "2016", "2017", "2018", "2019", "2020", "2021", "2022", "2023")),
     Q_0.025 = as.numeric(c(qmks3[1], qmks4[1], qmks5[1], qmks6[1], qmks7[1], qmks8[1], qmks9[1], qmks10[1], qmks11[1], qmks12[1], qmks13[1], qmks14[1], qmks15[1])),
     Q_0.5   = as.numeric(c(qmks3[2], qmks4[2], qmks5[2], qmks6[2], qmks7[2], qmks8[2], qmks9[2], qmks10[2], qmks11[2], qmks12[2], qmks13[2], qmks14[2], qmks15[2])),
     Q_0.975 = as.numeric(c(qmks3[3], qmks4[3], qmks5[3], qmks6[3], qmks7[3], qmks8[3], qmks9[3], qmks10[3], qmks11[3], qmks12[3], qmks13[3], qmks14[3], qmks15[3]))
 )


write.csv(mks, file="qmks.csv", row.names = FALSE)

qfks3 <- quantile(fks11,probs=c(0.025,0.5,0.975))
qfks4 <- quantile(fks12,probs=c(0.025,0.5,0.975))
qfks5 <- quantile(fks13,probs=c(0.025,0.5,0.975))
qfks6 <- quantile(fks14,probs=c(0.025,0.5,0.975))
qfks7 <- quantile(fks15,probs=c(0.025,0.5,0.975))
qfks8 <- quantile(fks16,probs=c(0.025,0.5,0.975))
qfks9 <- quantile(fks17,probs=c(0.025,0.5,0.975))
qfks10 <- quantile(fks18,probs=c(0.025,0.5,0.975))
qfks11 <- quantile(fks19,probs=c(0.025,0.5,0.975))
qfks12 <- quantile(fks20,probs=c(0.025,0.5,0.975))
qfks13 <- quantile(fks21,probs=c(0.025,0.5,0.975))
qfks14 <- quantile(fks22,probs=c(0.025,0.5,0.975))
qfks15 <- quantile(fks23,probs=c(0.025,0.5,0.975))

fks <- data.frame(
    Variable = as.numeric(c("2011", "2012", "2013", "2014", "2015", "2016", "2017", "2018", "2019", "2020", "2021", "2022", "2023")),
     Q_0.025 = as.numeric(c(qfks3[1], qfks4[1], qfks5[1], qfks6[1], qfks7[1], qfks8[1], qfks9[1], qfks10[1], qfks11[1], qfks12[1], qfks13[1], qfks14[1], qfks15[1])),
     Q_0.5   = as.numeric(c(qfks3[2], qfks4[2], qfks5[2], qfks6[2], qfks7[2], qfks8[2], qfks9[2], qfks10[2], qfks11[2], qfks12[2], qfks13[2], qfks14[2], qfks15[2])),
     Q_0.975 = as.numeric(c(qfks3[3], qfks4[3], qfks5[3], qfks6[3], qfks7[3], qfks8[3], qfks9[3], qfks10[3], qfks11[3], qfks12[3], qfks13[3], qfks14[3], qfks15[3]))
 )

write.csv(fks, file="qfks.csv", row.names = FALSE)
 
qofks3 <- quantile(ofks11,probs=c(0.025,0.5,0.975))
qofks4 <- quantile(ofks12,probs=c(0.025,0.5,0.975))
qofks5 <- quantile(ofks13,probs=c(0.025,0.5,0.975))
qofks6 <- quantile(ofks14,probs=c(0.025,0.5,0.975))
qofks7 <- quantile(ofks15,probs=c(0.025,0.5,0.975))
qofks8 <- quantile(ofks16,probs=c(0.025,0.5,0.975))
qofks9 <- quantile(ofks17,probs=c(0.025,0.5,0.975))
qofks10 <- quantile(ofks18,probs=c(0.025,0.5,0.975))
qofks11 <- quantile(ofks19,probs=c(0.025,0.5,0.975))
qofks12 <- quantile(ofks20,probs=c(0.025,0.5,0.975))
qofks13 <- quantile(ofks21,probs=c(0.025,0.5,0.975))
qofks14 <- quantile(ofks22,probs=c(0.025,0.5,0.975))
qofks15 <- quantile(ofks23,probs=c(0.025,0.5,0.975))

ofks <- data.frame(
    Variable = as.numeric(c("2011", "2012", "2013", "2014", "2015", "2016", "2017", "2018", "2019", "2020", "2021", "2022", "2023")),
     Q_0.025 = as.numeric(c(qofks3[1], qofks4[1], qofks5[1], qofks6[1], qofks7[1], qofks8[1], qofks9[1], qofks10[1], qofks11[1], qofks12[1], qofks13[1], qofks14[1], qofks15[1])),
     Q_0.5   = as.numeric(c(qofks3[2], qofks4[2], qofks5[2], qofks6[2], qofks7[2], qofks8[2], qofks9[2], qofks10[2], qofks11[2], qofks12[2], qofks13[2], qofks14[2], qofks15[2])),
     Q_0.975 = as.numeric(c(qofks3[3], qofks4[3], qofks5[3], qofks6[3], qofks7[3], qofks8[3], qofks9[3], qofks10[3], qofks11[3], qofks12[3], qofks13[3], qofks14[3], qofks15[3]))
 )
 
 write.csv(ofks, file="qofks.csv", row.names = FALSE)
 
qsrks3 <- quantile(srks11,probs=c(0.025,0.5,0.975))
qsrks4 <- quantile(srks12,probs=c(0.025,0.5,0.975))
qsrks5 <- quantile(srks13,probs=c(0.025,0.5,0.975))
qsrks6 <- quantile(srks14,probs=c(0.025,0.5,0.975))
qsrks7 <- quantile(srks15,probs=c(0.025,0.5,0.975))
qsrks8 <- quantile(srks16,probs=c(0.025,0.5,0.975))
qsrks9 <- quantile(srks17,probs=c(0.025,0.5,0.975))
qsrks10 <- quantile(srks18,probs=c(0.025,0.5,0.975))
qsrks11 <- quantile(srks19,probs=c(0.025,0.5,0.975))
qsrks12 <- quantile(srks20,probs=c(0.025,0.5,0.975))
qsrks13 <- quantile(srks21,probs=c(0.025,0.5,0.975))
qsrks14 <- quantile(srks22,probs=c(0.025,0.5,0.975))
qsrks15 <- quantile(srks23,probs=c(0.025,0.5,0.975))

srks <- data.frame(
    Variable = as.numeric(c("2011", "2012", "2013", "2014", "2015", "2016", "2017", "2018", "2019", "2020", "2021", "2022", "2023")),
     Q_0.025 = as.numeric(c(qsrks3[1], qsrks4[1], qsrks5[1], qsrks6[1], qsrks7[1], qsrks8[1], qsrks9[1], qsrks10[1], qsrks11[1], qsrks12[1], qsrks13[1], qsrks14[1], qsrks15[1])),
     Q_0.5   = as.numeric(c(qsrks3[2], qsrks4[2], qsrks5[2], qsrks6[2], qsrks7[2], qsrks8[2], qsrks9[2], qsrks10[2], qsrks11[2], qsrks12[2], qsrks13[2], qsrks14[2], qsrks15[2])),
     Q_0.975 = as.numeric(c(qsrks3[3], qsrks4[3], qsrks5[3], qsrks6[3], qsrks7[3], qsrks8[3], qsrks9[3], qsrks10[3], qsrks11[3], qsrks12[3], qsrks13[3], qsrks14[3], qsrks15[3]))
 )

write.csv(srks, file="srks.csv", row.names = FALSE)

qosrks3 <- quantile(osrks11,probs=c(0.025,0.5,0.975))
qosrks4 <- quantile(osrks12,probs=c(0.025,0.5,0.975))
qosrks5 <- quantile(osrks13,probs=c(0.025,0.5,0.975))
qosrks6 <- quantile(osrks14,probs=c(0.025,0.5,0.975))
qosrks7 <- quantile(osrks15,probs=c(0.025,0.5,0.975))
qosrks8 <- quantile(osrks16,probs=c(0.025,0.5,0.975))
qosrks9 <- quantile(osrks17,probs=c(0.025,0.5,0.975))
qosrks10 <- quantile(osrks18,probs=c(0.025,0.5,0.975))
qosrks11 <- quantile(osrks19,probs=c(0.025,0.5,0.975))
qosrks12 <- quantile(osrks20,probs=c(0.025,0.5,0.975))
qosrks13 <- quantile(osrks21,probs=c(0.025,0.5,0.975))
qosrks14 <- quantile(osrks22,probs=c(0.025,0.5,0.975))
qosrks15 <- quantile(osrks23,probs=c(0.025,0.5,0.975))

osrks <- data.frame(
    Variable = as.numeric(c("2011", "2012", "2013", "2014", "2015", "2016", "2017", "2018", "2019", "2020", "2021", "2022", "2023")),
     Q_0.025 = as.numeric(c(qosrks3[1], qosrks4[1], qosrks5[1], qosrks6[1], qosrks7[1], qosrks8[1], qosrks9[1], qosrks10[1], qosrks11[1], qosrks12[1], qosrks13[1], qosrks14[1], qosrks15[1])),
     Q_0.5   = as.numeric(c(qosrks3[2], qosrks4[2], qosrks5[2], qosrks6[2], qosrks7[2], qosrks8[2], qosrks9[2], qosrks10[2], qosrks11[2], qosrks12[2], qosrks13[2], qosrks14[2], qosrks15[2])),
     Q_0.975 = as.numeric(c(qosrks3[3], qosrks4[3], qosrks5[3], qosrks6[3], qosrks7[3], qosrks8[3], qosrks9[3], qosrks10[3], qosrks11[3], qosrks12[3], qosrks13[3], qosrks14[3], qosrks15[3]))
 )

write.csv(osrks, file="qosrks.csv", row.names = FALSE)


#PLOT

library(ggplot2)
library(ggpattern) 


#Konjsko

MFKS<-ggplot(data=mks, aes(x = Variable, y = Q_0.5)) +
    geom_ribbon(data=mks,  aes(ymin= Q_0.025, ymax = Q_0.975), fill = "white", color="black")+
    geom_line(data=mks, aes(x=Variable, y=Q_0.5), colour="black", linewidth=1) +
    labs(x = "Year", y = "Sex-specific populations size", title="Konjsko", caption="adult males = white fill with black line; adult females = grey fill with white line; reproductively active females = transparent fill with dashed line") +
    scale_x_continuous(breaks = seq(from = 2011, to = 2023, by = 2)) +
    scale_y_continuous(limits = c(0, 1134)) +
    theme_minimal() +
    geom_ribbon(data=fks, aes(ymin= Q_0.025, ymax = Q_0.975), 
                fill = alpha("grey", 0.8))+
    geom_line(data=fks, aes(x=Variable, y=Q_0.5), linewidth=1, colour="white")+
    geom_ribbon(data=ofks, aes(ymin= Q_0.025, ymax = Q_0.975), fill = alpha("white", 0.2), color="black")+
    geom_line(data=ofks, aes(x=Variable, y=Q_0.5, ), colour=alpha("black"), linewidth=1, linetype="twodash")+
    theme(text = element_text(size = 8), plot.title = element_text(face = "bold", hjust = 0.5, vjust = 0))

     
MFKScol<-ggplot(data=mks, aes(x = Variable, y = Q_0.5)) +
  geom_ribbon(data=mks,  aes(ymin= Q_0.025, ymax = Q_0.975), fill = "#56b4e9")+
  geom_line(data=mks, aes(x=Variable, y=Q_0.5), colour="white", linewidth=0.5) +
  labs(x = "Year", y = "Sex-specific populations size") +
  scale_x_continuous(breaks = seq(from = 2009, to = 2023, by = 2)) +
  scale_y_continuous(limits = c(0, 1134)) +
  theme_minimal() +
  geom_ribbon(data=fks, aes(ymin= Q_0.025, ymax = Q_0.975), 
              fill = alpha("#e69f00", 0.7))+
  geom_line(data=fks, aes(x=Variable, y=Q_0.5), linewidth=0.5, colour="white")+
  geom_ribbon(data=ofks, aes(ymin= Q_0.025, ymax = Q_0.975), fill = alpha("#cc79a7", 0.7), colour="#cc79a7")+
  geom_line(data=ofks, aes(x=Variable, y=Q_0.5, ), colour="white", linewidth=0.5)+
  theme(text = element_text(size = 10), axis.title = element_text(size=12), plot.title = element_text(face = "bold", hjust = 0.5, vjust = 0), panel.grid.major = element_blank(), panel.grid.minor = element_blank())       

     
     
     
#Plateau

MFP<-ggplot(data=mp, aes(x = Variable, y = Q_0.5)) +
    geom_ribbon(data=mp,  aes(ymin= Q_0.025, ymax = Q_0.975), fill = "white", color="black")+
    geom_line(data=mp, aes(x=Variable, y=Q_0.5), colour="black", linewidth=1) +
    labs(x = "Year", y = "Sex-specific populations size", title="Golem Grad Plateau", caption="adult males = white fill with black line; adult females = grey fill with white line; reproductively active females = transparent fill with dashed line") +
    scale_x_continuous(breaks = seq(from = 2009, to = 2023, by = 2)) +
    scale_y_continuous(limits = c(0, 1145)) +
    theme_minimal() +
    geom_ribbon(data=fp, aes(ymin= Q_0.025, ymax = Q_0.975), 
                fill = alpha("grey", 0.8))+
    geom_line(data=fp, aes(x=Variable, y=Q_0.5), linewidth=1, colour="white")+
    geom_ribbon(data=ofp, aes(ymin= Q_0.025, ymax = Q_0.975), fill = alpha("white", 0.2), color="black")+
    geom_line(data=ofp, aes(x=Variable, y=Q_0.5, ), colour=alpha("black"), linewidth=1, linetype="twodash")+
    theme(text = element_text(size = 8), plot.title = element_text(face = "bold", hjust = 0.5, vjust = 0))
     
     
MFPcol<-ggplot(data=mp, aes(x = Variable, y = Q_0.5)) +
  geom_ribbon(data=mp,  aes(ymin= Q_0.025, ymax = Q_0.975), fill = "#56b4e9")+
  geom_line(data=mp, aes(x=Variable, y=Q_0.5), colour="white", linewidth=0.5) +
  labs(x = "Year", y = "Sex-specific populations size") +
  scale_x_continuous(breaks = seq(from = 2009, to = 2023, by = 2)) +
  scale_y_continuous(limits = c(0, 1145)) +
  theme_minimal() +
  geom_ribbon(data=fp, aes(ymin= Q_0.025, ymax = Q_0.975), 
              fill = alpha("#e69f00", 0.8))+
  geom_line(data=fp, aes(x=Variable, y=Q_0.5), linewidth=0.5, colour="white")+
  geom_ribbon(data=ofp, aes(ymin= Q_0.025, ymax = Q_0.975), fill = "#cc79a7", colour="#cc79a7")+
  geom_line(data=ofp, aes(x=Variable, y=Q_0.5, ), colour="white", linewidth=0.5)+
  theme(text = element_text(size = 10), axis.title = element_text(size=12), plot.title = element_text(face = "bold", hjust = 0.5, vjust = 0), panel.grid.major = element_blank(), panel.grid.minor = element_blank())      
     
#Beaches
MFB<-ggplot(data=mb, aes(x = Variable, y = Q_0.5)) +
    geom_ribbon(data=mb,  aes(ymin= Q_0.025, ymax = Q_0.975), fill = "white", color="black")+
    geom_line(data=mb, aes(x=Variable, y=Q_0.5), colour="black", linewidth=1) +
    labs(x = "Year", y = "Sex-specific populations size", title="Golem Grad Beaches", caption="adult males = white fill with black line; adult females = grey fill with white line; reproductively active females = transparent fill with dashed line") +
    scale_x_continuous(breaks = seq(from = 2009, to = 2023, by = 2)) +
    scale_y_continuous(limits = c(0, 148)) +
    theme_minimal() +
    geom_ribbon(data=fb, aes(ymin= Q_0.025, ymax = Q_0.975), 
                fill = alpha("grey", 0.8))+
    geom_line(data=fb, aes(x=Variable, y=Q_0.5), linewidth=1, colour="white")+
    geom_ribbon(data=ofb, aes(ymin= Q_0.025, ymax = Q_0.975), fill = alpha("white", 0.2), color="black")+
    geom_line(data=ofb, aes(x=Variable, y=Q_0.5, ), colour=alpha("black"), linewidth=1, linetype="twodash")+
    theme(text = element_text(size = 8), plot.title = element_text(face = "bold", hjust = 0.5, vjust = 0))
    
    
MFBcol<-ggplot(data=mb, aes(x = Variable, y = Q_0.5)) +
  geom_ribbon(data=mb,  aes(ymin= Q_0.025, ymax = Q_0.975), fill = "#56b4e9")+
  geom_line(data=mb, aes(x=Variable, y=Q_0.5), colour="white", linewidth=0.5) +
  labs(x = "Year", y = "Sex-specific populations size") +
  scale_x_continuous(breaks = seq(from = 2009, to = 2023, by = 2)) +
  scale_y_continuous(limits = c(0, 148)) +
  theme_minimal() +
  geom_ribbon(data=fb, aes(ymin= Q_0.025, ymax = Q_0.975), 
              fill = alpha("#e69f00", 0.8))+
  geom_line(data=fb, aes(x=Variable, y=Q_0.5), linewidth=0.5, colour="white")+
  geom_ribbon(data=ofb, aes(ymin= Q_0.025, ymax = Q_0.975), fill = alpha("#cc79a7", 0.8), colour="#cc79a7")+
  geom_line(data=ofb, aes(x=Variable, y=Q_0.5, ), colour="white", linewidth=0.5)+
  theme(text = element_text(size = 10), axis.title = element_text(size=12), plot.title = element_text(face = "bold", hjust = 0.5, vjust = 0), panel.grid.major = element_blank(), panel.grid.minor = element_blank())     
     
     
#plot sex ratios

SRPBKScolCB<-ggplot(bsr, aes(x = Variable, y = Q_0.5)) +
  geom_ribbon(data=psr, aes(ymin= Q_0.025, ymax = Q_0.975), 
              fill = alpha("#56b4e9",0.7))+
  geom_line(data=psr, aes(x=Variable, y=Q_0.5), linewidth=0.5, colour="white") +
  labs(x = "Year", y = "Sex ratio") +
  scale_x_continuous(breaks = seq(from = 2009, to = 2023, by = 1), expand = c(0, 0)) +
  scale_y_continuous(breaks = c(1, 10, 20, 30, 40, 50), expand = c(0, 0)) +
  theme_minimal() +
  geom_ribbon(data=bsr, 
              aes(ymin= Q_0.025, ymax = Q_0.975),
              fill = alpha("#e69f00", 0.7))+
  geom_line(data=bsr, aes(x=Variable, y=Q_0.5), linewidth=0.5, colour="white") +
  geom_hline(yintercept = 1, color = alpha("lightgrey", 0.9) ,linewidth = 5) +
  geom_line(data=srks, aes(x=Variable, y=Q_0.5), linewidth=2, colour=alpha("#cc79a7")) +
  geom_line(data=srks, aes(x=Variable, y=Q_0.5), linewidth=2, colour="black", linetype="dotted") +
  theme(text = element_text(size = 10), axis.title = element_text(size=12), plot.title = element_text(face = "bold", hjust = 0.5, vjust = 0), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

MFP
MFPcol
MFKS
MFKScol
MFB
MFBcol
SRPBKScolCB

#printplot
library(Cairo)
ggsave("plot.png", plot = SRPBKScolCB, width = 6, height = 4, dpi = 300, type = "cairo")
ggsave("plot.pdf", plot = SRPBKScolCB, width = 6, height = 4)
