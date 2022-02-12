##this function returns the age-conditional risk/lifetime risk from incidence rate ratio, given current age and length of follow-up

#********************************************************************************
# a= start age
# s = number of years of follow-up
# incidRR = incidence rate ratio compared to the general population
# incidFile:  four-column file specifying the relevant paramters for risk estimation
# The fisrt 2 columns specify the start age and end age of each age group (like from 10 to 15, 15 to 20)
# The 3rd column is the incidence rate
# The 4th column is the net mortality rate (ie mortality rate in disease-free population)
#
#     If incidFile is not given, one can also provide in the next 3 arguments the mid-age in each age group (mid.age)
#   , incidence rates (incidence) and net mortality (net.mortality)

# headerPres : Whether header is present in incidFile
#********************************************************************************


LifeRisk <- function(a = 10, s = 100, incidRR, incidFile = NA,
    mid.age = NA, incidence = NA, net.mortality = NA, headerPres = T) {

    RR = incidRR

    if (!is.na(incidFile)) {
        file = read.table(incidFile, header = headerPres)
        #colnames(file) <-  c('age.start','age.end', 'incidence','net.mortality' )
        mid.age = (file[, 1] + file[, 2])/2
        incidence = file[, 3]
        net.mortality = file[, 4]
    }

    netmor.func <- approxfun(x = mid.age, y = net.mortality,
        rule = 2)
    incidence.func <- approxfun(x = mid.age, y = incidence, rule = 2)


    survival.func <- function(a, RR, t) {
        F.integrand <- function(x) {
            RR * incidence.func(x) + netmor.func(x)
        }
        F.res = integrate(F.integrand, lower = a, upper = t)$value

        return(exp(-F.res))
    }


    p.integrand <- function(t) {
        incidence.func(t) * survival.func(a = a, RR = RR, t)
    }

    absrisk = RR * integrate(Vectorize(p.integrand), lower = a,
        upper = a + s)$value
    return(absrisk)
}



