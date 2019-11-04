# AWAPer
This R package builds netCDF files of the Bureau of Meteorology Australian Water Availability Project daily national climate grids and allows efficient extraction of daily catchment average precipitation, Tmin, Tmax, vapour pressure, solar radiation and then estimation of areal potential evaporation (Morton's). 

Using the default compression settings, each meteorlogical variable requires ~5GB of hard-drive storage for the full record (1900 to 2019). Additionally, the netCDF files should be stored locally, and not over a network, to minimise the time for data extraction. For package details see the PDF manual https://github.com/peterson-tim-j/AWAPer/blob/master/AWAPer.pdf.

Below are details of how to install AWAPer followed by four examples of how to use it.

# Getting Started
To get started using the package, do the following steps. Below this are two examples for building and extracting data:

1. Download the latest release of the package from https://github.com/peterson-tim-j/AWAPer/releases.
1. Unzip the package on your local machine.
1. Rename the unzipped files, say "AWAPer-1.0.zip" folder to "AWAPer".
1. Open R. 
1. Install the netCDF package using the following command: `install.packages("ncdf4")` . **Importantly** this step may require installation of netCDF software outside of R. Please read the output R console messages carefully.
1. Install the remaining required packages using the following R command:
`install.packages(c("R.utils", "sp", "raster", "chron", "maptools", "Evapotranspiration","devtools"))`
1. Install the AWAPer package using the following example R command (NOTE: use the full file path to the AWAPer folder):For PC `install.packages("C:\MY_FOLDER\AWAPer\", repos = NULL, type = "source")` and for Mac `install.packages(“~/Users/MyFolder/AWAPer/", repos = NULL, type = "source")`

# Example 1. Build netCDF files

This example shows the steps required to build the two netCDF files, each containing 1 years of data.

```R
# Set working directory.
setwd('~/')`

# Make two netCDF files of AWAP data.
makeNetCDF_file(ncdfFilename='AWAP_demo.nc', ncdfSolarFilename='AWAP_solar_demo.nc', 
                updateFrom='2010-1-1',updateTo='2011-12-1')
```
# Example 2. Extract point precip. data and check with osberved data.
This example was developed by Ms Xinyang Fan (Uni. Melbourne) and shows how to extract point estimates of daily precipitation at four goundwater bore locations and at one rainfall gauge. The extracted data is then plotted. The rain gauge is also compared against the observed rain gauge. The latter shows that the results are unbiased, but minor differences do exist due to AWAP data having a 5x5 km grid-cell resolution. The plots below show (1) the locattions of the five sites (2) bar graphs of the daily precip. and (3) plots of the observed vs AWAPer estimated precip. at the rainfall gauge.

Importantly, this example downloads the Australia 9-second DEM. This is ~3GB.

![Example 2. Site locations](https://user-images.githubusercontent.com/8623994/68077541-ee782700-fe19-11e9-93fa-c0eae4606af1.png)

![Example 2. Daily precip at sites](https://user-images.githubusercontent.com/8623994/68077652-9d693280-fe1b-11e9-9c2e-fbebce8013cd.png)

![Example 3. Obs. vs Est Precip](https://user-images.githubusercontent.com/8623994/68077658-c7baf000-fe1b-11e9-9074-10d22b1bf8cd.png)

```R
# load all necessary packages
library(AWAPer)

# Make the netCDF files of only AWAP precipitation data.
fnames = makeNetCDF_file(ncdfFilename ='AWAP_demo.nc',
              updateFrom=as.Date("2010-08-01","%Y-%m-%d"),
              updateTo=as.Date("2010-10-01","%Y-%m-%d"),
              urlTmin=NA, urlTmax=NA, urlVprp=NA, urlSolarrad=NA)

# Download and import the DEM if it's not in the working directory
if (!file.exists('DEM.RData')) {
  DEM_9s = getDEM()
  save(DEM_9s,file="DEM.RData" )
} else {
  load('DEM.RData')
}

# Set coordinates to four bores locations and one rainfall gauge.
coordinates.data = data.frame( ID =c('Bore-70015656','Bore-50038039','Bore-20057861','Bore-10084446','Rain-63005'),
                   Longitude = c(131.33588, 113.066933, 143.118263, 153.551875, 149.5559),
                   Latitude =  c(-12.660622, -25.860046, -38.668506,-28.517974,-33.4289))

# Convert the points to a spatial object
sp::coordinates(coordinates.data) <- ~Longitude + Latitude

# Set projection to GDA94
sp::proj4string(coordinates.data) = '+proj=longlat +ellps=GRS80'

# Plot coordinates on top of DEM 9s
raster::plot(DEM_9s)
sp::plot(coordinates.data, add =T)
with(coordinates.data, text(sp::coordinates(coordinates.data)[,1],sp::coordinates(coordinates.data)[,2],
                            labels = coordinates.data$ID, pos = 1))

# Extract time-series of daily precip data at all five sites
climateData.data = extractCatchmentData(ncdfFilename='AWAP_demo.nc',
                   extractFrom=as.Date("2010-08-01","%Y-%m-%d"),
                   extractTo=as.Date("2010-10-01","%Y-%m-%d"),
                   catchments=coordinates.data,
                   getTmin=F, getTmax=F, getVprp=F, getSolarrad=F, getET=F)


# Plot the daily rainfall at each site
par.default = par()
par(mfrow=c(5,1), mar = c(4,5,3,0))
for (i in 1:nrow(coordinates.data)){
    filt = climateData.data$CatchmentID.ID == coordinates.data$ID[i]

    data2plot = climateData.data$precip_mm[filt]
    names(data2plot) = paste(climateData.data$day[filt],'/',climateData.data$month[filt],sep='')
    barplot(data2plot, main=coordinates.data$ID[i], ylab='Precip [mm/d]', xlab='Date [day/month]')

}

# The following is hard-coded observed rainfall for gauge  63005
obsPrecip <- data.frame(
         year= c(2010, 2010, 2010, 2010, 2010, 2010, 2010, 2010, 2010, 2010, 2010, 2010, 2010, 2010, 2010, 2010, 2010, 2010, 2010, 2010, 2010,
                 2010, 2010, 2010, 2010, 2010, 2010, 2010, 2010, 2010, 2010, 2010, 2010, 2010, 2010, 2010, 2010, 2010, 2010, 2010, 2010, 2010,
                 2010, 2010, 2010, 2010, 2010, 2010, 2010, 2010, 2010, 2010, 2010, 2010, 2010, 2010, 2010, 2010, 2010, 2010, 2010, 2010),
         month = c(8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9,
                 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 10),
         day = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 1, 2, 3, 4, 5,
                 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 1),
         precip_mm = c(0.6, 5.2, 0.8, 0.4, 0.0, 0.0, 0.0, 0.0, 0.0, 15.8, 15.6, 7.6, 0.7, 0.4, 1.4, 1.0, 0.0, 0.0, 30.4, 1.0, 0.0,
                 0.0, 5.0, 2.2, 0.3, 0.8, 13.8, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 13.8, 15.2, 0.4, 0.2, 0.0, 0.0, 12.4, 0.9,
                 0.0, 0.0, 0.1, 13.0, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.3, 0.0, 0.0, 0.0, 0.0))

# Plot the observed vs AWAPer rainfall at gauge ID 63005
par(mfrow=c(1,2), mar = c(4,5,3,2))
filt2 = climateData.data$CatchmentID.ID=='Rain-63005'
plot(obsPrecip$precip_mm,climateData.data$precip_mm[filt2],
     xlim = c(0,35),ylim = c(0,35),
     main='Obs. vs. AWAPer precip. at 63005',
     xlab='Obs. precip [mm/d]', ylab='AWAPer Precip [mm/d]')
abline(0,1, col='grey', lty=2)

# Plot the cumulative observed vs AWAPer rainfall at gauge ID 63005
plot(cumsum(obsPrecip$precip_mm),cumsum(climateData.data$precip_mm[filt2]),
     xlim = c(0,175),ylim = c(0,175),
     main='Cumulative obs. vs. AWAPer precip. at 63005',
     xlab='Obs. precip. [mm]', ylab='AWAPer precip. [mm]', type='l')
abline(0,1, col='grey', lty=2)

# Rest plotting parameters
par(par.default)
```

# Example 3. Calculate precip and evapotranspiration
This example calculates the catchment weighted precipitation at Bet Bet Creek (Victoria, Australia), the spatial standard deviation in precipitation and two measures of potential evapotranspiration. The example was developed by Dr Conrad Wasko. Below is a plot of the output.

![Example 3. Calculate precip and evapotranspiration](https://user-images.githubusercontent.com/39328041/66539735-e14a7e00-eb74-11e9-9890-dc3bce4148f2.png)

```R
# Make the netCDF files of AWAP data
makeNetCDF_file(ncdfFilename ='AWAP_demo.nc',ncdfSolarFilename='AWAP_solar_demo.nc', 
                updateFrom=as.Date("2010-1-1","%Y-%m-%d"),updateTo=as.Date("2011-12-1","%Y-%m-%d"))

# Download and import the DEM
DEM_9s = getDEM()

# Load example catchment boundaries.
data("catchments")

# Load the ET constants
data(constants)

# Extract catchment average data for Bet Bet Creek with 
# the Jensen Haise estimate of potential ET.
climateData.ET.JensenHaise.var = extractCatchmentData(ncdfFilename='AWAP_demo.nc',   
                                                      ncdfSolarFilename='AWAP_solar_demo.nc', extractFrom=as.Date("2010-1-1","%Y-%m-%d"),
                                                      extractTo=as.Date("2010-12-31","%Y-%m-%d"), catchments=catchments[2,], 
                                                      DEM=DEM_9s, ET.function='ET.JensenHaise',
                                                      ET.timestep='daily', ET.constants=constants);

# Extract catchment average data for Bet Bet Creek with 
# the Mortons CRAE estimate of potential ET.
climateData.ET.MortonCRAE.var = extractCatchmentData(ncdfFilename='AWAP_demo.nc',   
                                                     ncdfSolarFilename='AWAP_solar_demo.nc', extractFrom=as.Date("2010-1-1","%Y-%m-%d"),
                                                     extractTo=as.Date("2010-12-31","%Y-%m-%d"), catchments=catchments[2,], 
                                                     DEM=DEM_9s, ET.function='ET.MortonCRAE',
                                                     ET.timestep='monthly', ET.constants=constants);

# Set up start and end date indices for plotting
srt.date = as.Date("2010-06-01","%Y-%m-%d")
end.date = as.Date("2010-10-01","%Y-%m-%d")

# Convert year, month and day columns from extractions to a date.
climateData.ET.JensenHaise.var.date = as.Date(paste0(climateData.ET.JensenHaise.var$catchmentAvg$year, "-", 
                                                     climateData.ET.JensenHaise.var$catchmentAvg$month, "-",
                                                     climateData.ET.JensenHaise.var$catchmentAvg$day))

climateData.ET.MortonCRAE.var.date = as.Date(paste0(climateData.ET.MortonCRAE.var$catchmentAvg$year, "-", 
                                                    climateData.ET.MortonCRAE.var$catchmentAvg$month, "-",
                                                    climateData.ET.MortonCRAE.var$catchmentAvg$day))

i1.s = match(srt.date, climateData.ET.JensenHaise.var.date)
i1.e = match(end.date, climateData.ET.JensenHaise.var.date)
i2.s = match(srt.date, climateData.ET.MortonCRAE.var.date)
i2.e = match(end.date, climateData.ET.MortonCRAE.var.date)

# Plot rainfall and standard deviation against observations
# ---------------------------------------------------------
max.y = max(climateData.ET.JensenHaise.var$catchmentAvg$precip_mm[i1.s:i1.e] + 
              sqrt(climateData.ET.JensenHaise.var$catchmentvar$precip_mm[i1.s:i1.e]))

# Change the plot margins
par(mar =  c(5, 7.5, 4, 2.7) + 0.1)

# Rainfall
plot(climateData.ET.JensenHaise.var.date[i1.s:i1.e],
     climateData.ET.JensenHaise.var$catchmentAvg$precip_mm[i1.s:i1.e],
     type = "h", col = "#e31a1c", lwd = 3, mgp = c(2, 0.5, 0), ylim = c(0, 80),
     ylab = "", xlab = "2010", xaxs = "i", yaxt = "n", bty = "l", yaxs = "i")
axis(side = 2, mgp = c(2, 0.5, 0), line = 0.5, at = seq(from = 0, to = 80, by = 20),
     labels = c("0", "20", "40", "60", "80mm"), col = "#e31a1c", col.axis = "#e31a1c")

# Standard deviation
for (i in 1:length(climateData.ET.JensenHaise.var.date[i1.s:i1.e])) {
  x.plot = rep(climateData.ET.JensenHaise.var.date[i1.s:i1.e][i], 2)
  y.plot = c(climateData.ET.JensenHaise.var$catchmentAvg$precip_mm[i1.s:i1.e][i] + 
               sqrt(climateData.ET.JensenHaise.var$catchmentvar$precip_mm[i1.s:i1.e][i]),
             climateData.ET.JensenHaise.var$catchmentAvg$precip_mm[i1.s:i1.e][i] - 
               sqrt(climateData.ET.JensenHaise.var$catchmentvar$precip_mm[i1.s:i1.e][i]))
  lines(x.plot, y.plot, col = "black", lwd = 1.2)
  
}

# Plot evap data.
par(new = TRUE)
plot(climateData.ET.JensenHaise.var.date[i1.s:i1.e], climateData.ET.JensenHaise.var$catchmentAvg$ET_mm[i1.s:i1.e], col = "#bc80bd", lwd = 2, ylab = "",
     ylim = c(0, 4), lty = 1, xlab = "", xaxs = "i", yaxt = "n", xaxt = "n", type = "l", bty = "n", yaxs = "i")
axis(side = 2, line = 2.3, mgp = c(2, 0.5, 0), labels = c("0", "1", "2", "3", "4mm"), at = seq(from = 0, to = 4, by = 1), col = "#bc80bd", col.axis = "#bc80bd")
lines(climateData.ET.MortonCRAE.var.date[i2.s:i2.e], climateData.ET.MortonCRAE.var$catchmentAvg$ET_mm[i2.s:i2.e], col = "#bc80bd", lwd = 2, lty = 2)

# Legend
legend("topleft", cex = 0.8, lwd = 2, bty = "n", inset = c(0.01, -0.01),
       lty = c(1, 1, 2), pch = c(NA, NA, NA), 
       col = c("#e31a1c",  "#bc80bd", "#bc80bd"),
       legend = c("Precipitation (bars +/- one standard dev.)",  "Jensen-Haise PET", "Morton CRAE PET"), xpd = NA)
```

# Example 4. Calculate evapotranspiration

This example calculates and plot various estimates of evaportrnspiration using the netCDF data form the previous example. To do this, the Austrlia 9 second DEM is downloaded and two catchment boundaries are loaded from the package. The figure below shows the output plot from the example. It shows 10 different esitmates of area weighted evapotranspiration from 1/1/2010 to 31/12/2010 at catchment 407214 (Victoria, Australia).

![Example 4. ET plot](https://user-images.githubusercontent.com/8623994/64103416-8d9a8700-cdb5-11e9-977a-ea8fabdcfcf5.png)

```R
# Set working directory.
setwd('~/')`

# Make two netCDF files of AWAP data.
makeNetCDF_file(ncdfFilename='AWAP_demo.nc', ncdfSolarFilename='AWAP_solar_demo.nc', 
                updateFrom='2010-1-1',updateTo='2011-12-1')
                
# Download and import the DEM
DEM_9s = getDEM()

# Load example cacthment boundaries.
data("catchments")

# get ET constants
data(constants)

# Extract data and esitmate various types of ET and estimate the spatial variance 
# using the interquartile range.
#----------------------------------------------
climateData.ET.HargreavesSamani = extractCatchmentData(ncdfFilename='AWAP_demo.nc',   
      ncdfSolarFilename='AWAP_solar_demo.nc', extractFrom=as.Date("2009-1-1","%Y-%m-%d"),
      extractTo=as.Date("2010-12-1","%Y-%m-%d"), catchments=catchments, 
      spatial.function.name='IQR',DEM=DEM, ET.function='ET.HargreavesSamani',
      ET.timestep = 'daily', ET.constants= constants);

climateData.ET.JensenHaise = extractCatchmentData(ncdfFilename='AWAP_demo.nc', 
      ncdfSolarFilename='AWAP_solar_demo.nc', extractFrom=as.Date("2009-1-1","%Y-%m-%d"),
      extractTo=as.Date("2010-12-1","%Y-%m-%d"), catchments=catchments, 
      spatial.function.name='IQR',DEM=DEM, ET.function='ET.JensenHaise',
      ET.timestep = 'daily', ET.constants= constants);

climateData.ET.Makkink = extractCatchmentData(ncdfFilename='AWAP_demo.nc',
      ncdfSolarFilename='AWAP_solar_demo.nc',extractFrom=as.Date("2009-1-1","%Y-%m-%d"),
      extractTo=as.Date("2010-12-1","%Y-%m-%d"),catchments=catchments, 
      spatial.function.name='IQR',DEM=DEM, ET.function='ET.Makkink'
      ET.timestep = 'daily', ET.constants= constants);

climateData.ET.McGuinnessBordne = extractCatchmentData(ncdfFilename='AWAP_demo.nc',
      ncdfSolarFilename='AWAP_solar_demo.nc', extractFrom=as.Date("2009-1-1","%Y-%m-%d"),
      extractTo=as.Date("2010-12-1","%Y-%m-%d"),catchments=catchments, 
      spatial.function.name='IQR',DEM=DEM, ET.function='ET.McGuinnessBordne',
      ET.timestep = 'daily', ET.constants= constants);

climateData.ET.MortonCRAE = extractCatchmentData(ncdfFilename='AWAP_demo.nc', 
      ncdfSolarFilename='AWAP_solar_demo.nc', extractFrom=as.Date("2009-1-1","%Y-%m-%d"),
      extractTo=as.Date("2010-12-1","%Y-%m-%d"), catchments=catchments, 
      spatial.function.name='IQR',DEM=DEM, ET.function='ET.MortonCRAE',
      ET.timestep = 'monthly', ET.constants= constants);

climateData.ET.MortonCRAE.potentialET = extractCatchmentData(ncdfFilename='AWAP_demo.nc',
      ncdfSolarFilename='AWAP_solar_demo.nc', extractFrom=as.Date("2009-1-1","%Y-%m-%d"),
      extractTo=as.Date("2010-12-1","%Y-%m-%d"),catchments=catchments, 
      spatial.function.name='IQR',DEM=DEM, ET.function='ET.MortonCRAE',
      ET.timestep = 'monthly', ET.Mortons.est='potential ET', ET.constants= constants);

climateData.ET.MortonCRAE.wetarealET = extractCatchmentData(ncdfFilename='AWAP_demo.nc',
      ncdfSolarFilename='AWAP_solar_demo.nc', extractFrom=as.Date("2009-1-1","%Y-%m-%d"),
      extractTo=as.Date("2010-12-1","%Y-%m-%d"), catchments=catchments, 
      spatial.function.name='IQR',DEM=DEM, ET.function='ET.MortonCRAE',
      ET.timestep = 'monthly', ET.Mortons.est='wet areal ET', ET.constants= constants);

climateData.ET.MortonCRAE.actualarealET = extractCatchmentData(ncdfFilename='AWAP_demo.nc',
      ncdfSolarFilename='AWAP_solar_demo.nc',extractFrom=as.Date("2009-1-1","%Y-%m-%d"), 
      extractTo=as.Date("2010-12-1","%Y-%m-%d"), catchments=catchments, 
      spatial.function.name='IQR',DEM=DEM, ET.function='ET.MortonCRAE',
      ET.timestep = 'monthly', ET.Mortons.est='actual areal ET', ET.constants= constants);

climateData.ET.MortonCRWE = extractCatchmentData(ncdfFilename='AWAP_demo.nc',
      ncdfSolarFilename='AWAP_solar_demo.nc', extractFrom=as.Date("2009-1-1","%Y-%m-%d"),
      extractTo=as.Date("2010-12-1","%Y-%m-%d"),catchments=catchments, 
      spatial.function.name='IQR',DEM=DEM, ET.function='ET.MortonCRWE',
      ET.timestep = 'monthly', ET.constants= constants);

climateData.ET.MortonCRWE.shallowLake = extractCatchmentData(ncdfFilename='AWAP_demo.nc',
      ncdfSolarFilename='AWAP_solar_demo.nc', extractFrom=as.Date("2009-1-1","%Y-%m-%d"),
      extractTo=as.Date("2010-12-1","%Y-%m-%d"), catchments=catchments, 
      spatial.function.name='IQR',DEM=DEM, ET.function='ET.MortonCRWE',
      ET.timestep = 'monthly', ET.Mortons.est = 'shallow lake ET', ET.constants= constants);

climateData.ET.Turc = extractCatchmentData(ncdfFilename='AWAP_demo.nc',
      ncdfSolarFilename='AWAP_solar_demo.nc', extractFrom=as.Date("2009-1-1","%Y-%m-%d"),
      extractTo=as.Date("2010-12-1","%Y-%m-%d"),catchments=catchments, 
      spatial.function.name='IQR',DEM=DEM, ET.function='ET.Turc',
      ET.timestep = 'daily', ET.constants= constants);

# Plot the ET estimates for one of the catchmnts.
#----------------------------------------
filt = climateData.ET.HargreavesSamani$catchmentAvg$CatchID==407214
d = ISOdate(climateData.ET.HargreavesSamani$catchmentAvg$year, 
      climateData.ET.HargreavesSamani$catchmentAvg$month, 
      climateData.ET.HargreavesSamani$catchmentAvg$day)
plot(d[filt], climateData.ET.HargreavesSamani$catchmentAvg$ET_mm[filt], 
      col='black',lty=1, xlim = c(ISOdate(2010,1,1), ISOdate(2010,12,1)), ylim=c(0, 30),type='l', ylab='ET [mm/d]',xlab='Date')

filt = climateData.ET.JensenHaise$catchmentAvg$CatchID==407214
d = ISOdate(climateData.ET.JensenHaise$catchmentAvg$year, 
      climateData.ET.JensenHaise$catchmentAvg$month, climateData.ET.JensenHaise$catchmentAvg$day)
lines(d[filt], climateData.ET.JensenHaise$catchmentAvg$ET_mm[filt], col='red',lty=1)

filt = climateData.ET.Makkink$catchmentAvg$CatchID==407214
d = ISOdate(climateData.ET.Makkink$catchmentAvg$year, climateData.ET.Makkink$catchmentAvg$month, 
      climateData.ET.Makkink$catchmentAvg$day)
lines(d[filt], climateData.ET.Makkink$catchmentAvg$ET_mm[filt], col='green',lty=1)

filt = climateData.ET.McGuinnessBordne$catchmentAvg$CatchID==407214
d = ISOdate(climateData.ET.McGuinnessBordne$catchmentAvg$year, climateData.ET.McGuinnessBordne$catchmentAvg$month,      
      climateData.ET.McGuinnessBordne$catchmentAvg$day)
lines(d[filt], climateData.ET.McGuinnessBordne$catchmentAvg$ET_mm[filt], col='blue',lty=1)

filt = climateData.ET.MortonCRAE.potentialET$catchmentAvg$CatchID==407214
d = ISOdate(climateData.ET.MortonCRAE.potentialET$catchmentAvg$year, 
      climateData.ET.MortonCRAE.potentialET$catchmentAvg$month, climateData.ET.MortonCRAE.potentialET$catchmentAvg$day)
lines(d[filt], climateData.ET.MortonCRAE.potentialET$catchmentAvg$ET_mm[filt], col='black',lty=2)

filt = climateData.ET.MortonCRAE.wetarealET$catchmentAvg$CatchID==407214
d = ISOdate(climateData.ET.MortonCRAE.wetarealET$catchmentAvg$year,     
      climateData.ET.MortonCRAE.wetarealET$catchmentAvg$month, climateData.ET.MortonCRAE.wetarealET$catchmentAvg$day)
lines(d[filt], climateData.ET.MortonCRAE.wetarealET$catchmentAvg$ET_mm[filt], col='red',lty=2)

filt = climateData.ET.MortonCRAE.actualarealET$catchmentAvg$CatchID==407214
d = ISOdate(climateData.ET.MortonCRAE.actualarealET$catchmentAvg$year, 
      climateData.ET.MortonCRAE.actualarealET$catchmentAvg$month, 
      climateData.ET.MortonCRAE.actualarealET$catchmentAvg$day)
lines(d[filt], climateData.ET.MortonCRAE.actualarealET$catchmentAvg$ET_mm[filt], col='green',lty=2)

filt = climateData.ET.MortonCRWE$catchmentAvg$CatchID==407214
d = ISOdate(climateData.ET.MortonCRWE$catchmentAvg$year, climateData.ET.MortonCRWE$catchmentAvg$month, 
      climateData.ET.MortonCRWE$catchmentAvg$day)
lines(d[filt], climateData.ET.MortonCRWE$catchmentAvg$ET_mm[filt], col='blue',lty=2)

filt = climateData.ET.MortonCRWE.shallowLake$catchmentAvg$CatchID==407214
d = ISOdate(climateData.ET.MortonCRWE.shallowLake$catchmentAvg$year, 
      climateData.ET.MortonCRWE.shallowLake$catchmentAvg$month, climateData.ET.MortonCRWE.shallowLake$catchmentAvg$day)
lines(d[filt], climateData.ET.MortonCRWE.shallowLake$catchmentAvg$ET_mm[filt], col='black',lty=3)

filt = climateData.ET.Turc$catchmentAvg$CatchID==407214
d = ISOdate(climateData.ET.Turc$catchmentAvg$year, climateData.ET.Turc$catchmentAvg$month, 
      climateData.ET.Turc$catchmentAvg$day)
lines(d[filt], climateData.ET.Turc$catchmentAvg$ET_mm[filt], col='red',lty=3)

legend(x='topright', legend=c(
  'Hargreaves Samani (ref. crop)', 'Jensen Haise (PET)', 'Makkink (ref. crop)', 'McGuinness Bordne (PET)', 'Morton CRAE (PET)',
  'Morton CRAE (wet areal ET)', 'Morton CRAE (actual areal ET)', 'Morton CRWE (PET)', 'Morton CRWE (shallowLake)', 'Turc (ref. crop, non-humid'),
  lty = c(1,1,1,1,2,2,2,2,3,3), col=c('black','red','green','blue','black','red','green','blue','black','red')
  )
  ```
