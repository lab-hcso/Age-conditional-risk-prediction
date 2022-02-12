#**********************************************************************************
# This is a function to convert prevalence odds ratio (POR) to 
# lifetime relative risk, with cosnideration of competing risks. 
# Effect size measures are assumed to be compared to the baseline group.
# 
#  PA: allele frequency of risk allele 
#  OR2: Prevalence odds ratio of the risk allele (assuming OR are multiplicative with addition of one risk allele)
#  incidFile: A file with 4 columns. First column is the start age of the age group, 
#         2nd column is the end age of the age group, the 3rd column is the incidence rate, the 4th column is net mortality rate
#  a: current age
#  s: period of follow-up in years
#  headerPres: whether a header is present for incidFile
#  Note: requires the function "ORcompareAvg"
#***********************************************************************************


ORtoLifetimeRR <- function(PA,OR2,incidFile,PrD=NA,a=10,s=100,headerPres=T) {

file = read.table(incidFile,header=headerPres)

mid.age= (file[,1]+file[,2])/2
incidence = file[,3]
net.mortality =file[,4]

netmor.func <- approxfun(x= mid.age, y=net.mortality,rule=2)   ##linear interpolation of net mortality/incidence
incidence.func <- approxfun(x=mid.age, y= incidence,rule=2) 

#netmor.func <- splinefun(x= mid.age, y=net.mortality)   ##spline interpolation of net mortality/incidence
#incidence.func <- splinefun(x=mid.age, y= incidence) 

survival.func<- function(a,RR,t) {
F.integrand <- function(x) {
RR*incidence.func(x) + netmor.func(x) }
F.res = integrate(F.integrand ,lower=a, upper =t)$value 

return( exp(-F.res) ) 
}

absrisk.func <- function(s,a,RR) {

p.integrand<- function(t){
incidence.func(t) * survival.func(a=a, RR=RR, t) }

absrisk = RR* integrate(Vectorize(p.integrand),lower =a, upper= a+s )$value
return(absrisk)
}

ORstar= ORcompareAvg(PrD=PrD, PA=PA, OR2=OR2)

##a is the current age; s is the length of follow-up (in yrs)
absrisk.avg = absrisk.func(s=s,a=a,RR=1) 
absrisk.aa = absrisk.func(s=s,a=a,RR=ORstar$ORstar.aa) 
absrisk.Aa = absrisk.func(s=s,a=a,RR=ORstar$ORstar.Aa) 
absrisk.AA = absrisk.func(s=s,a=a,RR=ORstar$ORstar.AA) 

lifeRR2 = absrisk.Aa/absrisk.aa
lifeRR3 = absrisk.AA/absrisk.aa

return( list(absrisk.avg=absrisk.avg, absrisk.aa=absrisk.aa, absrisk.Aa=absrisk.Aa ,absrisk.AA=absrisk.AA,
lifeRR2=lifeRR2, lifeRR3=lifeRR3)   ) 

}   



