# mortonsAPETcalc_V1.1.r
#  Script calculates Morton's aPET on the fly for input into
#  the CWYET rainfall-runoff modelling tool via the extractDataForCWYET_V2.1.r. 
#  This script must be compiled prior to running extractDataForCWYET_V2.1.r.
#  This script has been tested with CSIRO Land and Water data and is operational.
#  Data will be output for into the extractDataForCWYET_V2.1.r script.
#  
# Authors:  Alison Oke, Andrew Frost, Derek Bacon.
# v1.1 28/07/2009 separate subroutine for calculating Morton's aPET within AWAP extract script.
######################################################################
mortonsAET <- function(tmin, tmax, radin, eact,lat,height,date){
#############################################
#Start Morton's calculations
## READ IN DATASETS THAT DETAIL THE LATITUDE AND ELEVATION OF EACH GRID CELL  


tavg = (tmax + tmin) / 2.0
eact<-eact/10

##Calculate esat
                    ## formula 35
                    esatMax = 0.6108 * exp(17.27 * tmax / (tmax + 237.3))
                    esatMin = 0.6108 * exp(17.27 * tmin / (tmin + 237.3))
                    esatMean = 0.6108 * exp(17.27 * tavg / (tavg + 237.3))

                    ## formula 43
                    ea = 0.25 * esatMax + 0.5 * esatMean + 0.25 * esatMin
                    esat = ea

                    if (esat <= 0.0) esat = 0.0001
                    if (eact <= 0.0) eact = 0.0001
                    if (eact > esat) eact = esat

                    RH = eact / esat

## Calculate ratio of atmospheric at the station to that at the sea level (Pstn/Psea)
                    pratio = ((293.0 - 0.0065 * height) / 293.0)^5.26

## Calculate psychrometric constant (kPa/C)
                    if (tavg >= 0.0){
                        gammap = 0.066 * pratio
                    } else {
                        gammap = 0.0574 * pratio
                    }
## Calculate slope of saturation vapour pressure/temperature curve (kPa/C)

                    if (tavg >= 0.0) {
                        pdelta = (4098.0 * esat) / ((tavg + 237.3)^2)
                    } else {
                        pdelta = (5809.0 * esat) / ((tavg + 265.5)^2)
                    }
## Calculate extraterrestrial radiation (Ra=radextra)   (MJ/m2/day)

                    psai = ((lat / 180.0) * pi)

# **will need input  julday = today.DayOfYear
                    julday<-as.numeric(paste(format(date, '%j')))
                    dr = (1.0 + 0.033 * cos(0.0172 * julday))
                    delta = 0.409 * sin(0.0172 * julday - 1.39)
                    omega = acos(-1.0 * tan(psai) * tan(delta))
                    radextra = (118.1 / pi) * dr * (omega * sin(psai) * sin(delta) + cos(psai) * cos(delta) * sin(omega))


## Calculate incoming solar radiation (Rs=radin) using Angstorm formula
                 ####solar radiation calc not require, have data from awap
                    #daylight = (24.0 / pi) * omega
                    #nNratio = sunhrs / daylight
                    #radin = (as + bs * nNratio) * radextra

## Calculate nNratio based on radin
                    as <- 0.25 #Angstorm formula, regression constant
                    bs <- 0.5 #Angstorm formula, regression constant
                    nNratio = ((radin / radextra - as) / bs)

## Calculate NET incoming solar radiation (Rns=radinnet)

                    albedo <- 0.23 # copied from CSIRO L&W script
                    sigma = 4.903E-9   # Stefan-Boltzmann constant
                    radinnet = (1.0 - albedo) * radin

## Calculate net outgoing longwave radiation (Rnl=radout)

                    term1 = ((((tmax + 273.16)^4) + (tmin + 273.16)^4) / 2.0)
                    term2 = (0.34 - 0.14 * sqrt(eact)) * (0.10 + 0.9 * nNratio)
                    radout = sigma * term1 * term2
                    if (radout < 0.0) radout = 0.0

##  Calculate net radiation (Rn) : if negative set it to zero

                    radnet = radinnet - radout
                    if (radnet < 0.0) radnet = 0.0

## Calculate stability factor(stabfac), vapour pressure transfer
## coefficient(fa=vptc)and heat transfer coeffieient(lamda=htc)

                    if (tavg >= 0.0) {
                        fz = 24.19
                    } else {
                        fz = 27.82
                        }

                    ediff = esat - eact
                    if (ediff <= 0.0) ediff = 0.0001
                    term3 = gammap * (sqrt(1.0 / pratio)) * fz * ediff   ## term3=gammap*((1.0/pratio)**0.5)*fz*ediff !!!can be wrong because of Jai
                    stabfac = 1.0 / (0.28 * (1.0 + eact / esat) + pdelta * radnet / term3)
                    if (stabfac < 1.0) stabfac = 1.0                     ##!!!!! MODIFICATION

                    vptc = (sqrt(1.0 / pratio)) * fz / stabfac           ##vptc=((1.0/pratio)**0.5)*fz/stabfac
                    htc = gammap + (1.804E-8 * (tavg + 273.0)^3) / vptc  ##htc=gammap+(1.804E-8*(tavg+273.0)**3)/vptc

## Carryout iterative procedure to satisfy the energy balance and obtain equlibrium quantities

                    xesat = esat
                    xtemp = tavg
                    xdelta = pdelta

                    for(i in 1 : 100000) ##** Bit unsure about this section it was a bit dfferent in the C# code
                    {
                        tempinc = (radnet / vptc + eact + htc * (tavg - xtemp) - xesat) / (xdelta + htc)
                        tdiff = abs(tempinc)
                        if (tdiff < 0.01) {
                            break
                        } else {
                        
                            xtemp = xtemp + tempinc
                            if (xtemp >= 0.0) {
                                xesat = 0.6108 * exp((17.27 * xtemp) / (xtemp + 237.3))
                            } else {
                                xesat = 0.6108 * exp((21.88 * xtemp) / (xtemp + 265.5))
                            }

                            if (xtemp >= 0.0) {
                                xdelta = (4098.0 * xesat) / (xtemp + 237.3)^2
                            } else {
                                xdelta = (5809.0 * xesat) / (xtemp + 265.5)^2
                            }
                        }

                    } 


                    ETppx = radnet - htc * vptc * (xtemp - tavg)
                    ETpp = ETppx * 0.408
                    if (ETpp < 0.0) ETpp = 0.0

## Calculate Morton Wet Environment Areal Potential Evapotranspiration. ETwp

                    radnetx = ETppx + gammap * vptc * (xtemp - tavg)
                    Dpfact = xdelta / (gammap + xdelta)
                    ETwp = 0.408 * (1.2096 + 1.2 * Dpfact * radnetx)

## Constraints applied such that calculated areal potential value is unchanged

                    if (ETwp > ETpp) ETpp = ETwp
                    if (ETwp < 0.5 * ETpp) ETpp = 2.0 * ETwp

## Calculate actual evapotranspiration from Morton's complementary relationship

                    ETa = 2.0 * ETwp - ETpp
                    if (ETa < 0.0) ETa = 0.0                ##!!!!Modification

## Calculate Reference crop evapotranspiration (ET0) using Penmann-Monteith Equation

                    #if (winsp <= -99.0)
                    #    ET0 = -99.0
                    #else
                    #{
                    #    term4 = 0.408 * pdelta * radnet + gammap * (900.0 / (tavg + 273.0)) * winsp * (esat - eact)
                    #    term5 = pdelta + gammap * (1.0 + 0.34 * winsp)
                    #    ET0 = term4 / term5
                    #}
                
 mort<-c(ETwp,RH)
 names(mort)<-c('mortAPET','mortRH')
 mort
#### ETwp and RH is what we want to output and we want to round them to 4 decimal places ( round(ETwp, 4) )
} #End function aetCalc
##End of Morton's Subscript################################
########################################################### 

