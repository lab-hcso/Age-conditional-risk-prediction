#**********************************************************************************
# This function calculates the odds ratio as compared to the general population
# PrD : probability of disease in general population
# PA: allele frequency of risk allele 
# OR2: Prevalence odds ratio of the risk allele (assuming OR are multiplicative with addition of one risk allele)
# **********************************************************************************

ORcompareAvg <- function(PrD, PA, OR2) {

OR3 = OR2^2    

Paa = (1-PA)^2      ##genotype freq
PAa = 2*PA*(1-PA)
PAA = PA^2 

func <- function(p) {
LHS = PrD 
RHS = p/(1+p)*Paa + OR2*p/(1+OR2*p)*PAa  + OR3*p/(1+OR3*p)*PAA
(LHS-RHS)^2 }


p = optimize(func, c(0,10) )$minimum

PrD.aa = p/(1+p) 
PrD.Aa = OR2*p/(1+OR2*p)
PrD.AA = OR3*p/(1+OR3*p) 

ORstar.aa =  PrD.aa/(1- PrD.aa) / (  PrD/(1-PrD) ) 
ORstar.Aa =  PrD.Aa/(1- PrD.Aa) / (  PrD/(1-PrD) ) 
ORstar.AA =  PrD.AA/(1- PrD.AA) / (  PrD/(1-PrD) ) 

return( list(ORstar.aa=ORstar.aa,ORstar.Aa=ORstar.Aa, ORstar.AA=ORstar.AA) )          ##OR as compared to the general population
}



PrD = 0.01 ## prob of disease
OR2 = 1.3
PA = 0.3  #risk allele freq 
ORcompareAvg (PrD,PA,OR2)
