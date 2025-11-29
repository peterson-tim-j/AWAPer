<!-- badges: start -->

[![CRANversion](http://www.r-pkg.org/badges/version/AWAPer)](https://cran.r-project.org/package=AWAPer) 
![Windows](https://github.com/peterson-tim-j/AWAPer/actions/workflows/windows-latest.yml/badge.svg)
![Ubuntu](https://github.com/peterson-tim-j/AWAPer/actions/workflows/ubuntu-latest.yml/badge.svg)
![Macos](https://github.com/peterson-tim-j/AWAPer/actions/workflows/macos-latest.yml/badge.svg)
[![Codecov](https://img.shields.io/codecov/c/github/peterson-tim-j/AWAPer?logo=CODECOV)](https://app.codecov.io/github/peterson-tim-j/AWAPer) 
[![CRAN RStudio mirror downloads](http://cranlogs.r-pkg.org/badges/grand-total/AWAPer)](https://cran.r-project.org/package=AWAPer)  [![CRAN RStudio mirror downloads](http://cranlogs.r-pkg.org/badges/AWAPer)](https://cran.r-project.org/package=AWAPer)
<!-- badges: end -->

# _AWAPer_ - an R-package for catchment-weighted climate data anywhere in Australia
This R package builds netCDF files of the Bureau of Meteorology Australian Water Availability Project daily national climate grids and allows efficient extraction of daily, weekly, monthly, quarterly or annual catchment average precipitation, Tmin, Tmax, vapour pressure, solar radiation and then estimation of various measures of potential evaporation. 

The package development was funded by the Victorian Government The Department of Environment, Land, Water and Planning _Climate and Water Initiate_ (https://www.water.vic.gov.au/climate-change/research/vicwaci). 

For details of the appraoch see the paper "Peterson, Tim J, Wasko, C, Saft, M, Peel, MC. AWAPer: An R package for area weighted catchment daily meteorological data anywhere within Australia. _Hydrological Processes_. 2020; 34: 1301– 1306. https://doi.org/10.1002/hyp.13637". 

For package details see the PDF manual https://github.com/peterson-tim-j/AWAPer/blob/master/AWAPer.pdf. Alternatively see the lastest CRAN published version at https://CRAN.R-project.org/package=AWAPer 

Package vignettes can also be found using the R-command:
```R
browseVignettes("AWAPer")
```

Using the default compression settings, each meteorlogical variable requires ~5GB of hard-drive storage for the full record (1900 to 2019). Additionally, the netCDF files should be stored locally, and not over a network, to minimise the time for data extraction. Below are details of the system requirements, how to install AWAPer and the following examples:
1. Building and updating the required netCDF files (![see here](https://github.com/peterson-tim-j/AWAPer#example-1-build-netcdf-files))
1. Extract and check daily point precipitation data against rain gauge observations (![see here](https://github.com/peterson-tim-j/AWAPer#example-2-extract-point-precip-data-and-check-with-osberved-data))
1. Extract daily areal weighted precipitation and calculate two measures of ET (![see here](https://github.com/peterson-tim-j/AWAPer#example-3-calculate-precip-and-evapotranspiration))
1. Calculate all measures of ET possible with AWAPer (![see here](https://github.com/peterson-tim-j/AWAPer#example-4-calculate-evapotranspiration))
1. Extract and map maximum daily temperature for all of Australia (![see here](https://github.com/peterson-tim-j/AWAPer#example-5-map-of-australia))
1. Extract monthly areal weighted total precipitation at two catchments and plot the totals and the spatial variability (![see here](https://github.com/peterson-tim-j/AWAPer/blob/master/README.md#example-6-extract-monthly-precipitation))
1. Extract and map the monthly total precipitation at two catchments (![see here](https://github.com/peterson-tim-j/AWAPer/blob/master/README.md#example-7-map-monthly-precipitation))

# System Requiements
On Windows OS only the program "7z" is required to uzip the ".Z" compressed grid files. Follow the steps below to download and install 7z.

1. Download and intall 7z from https://www.7-zip.org/download.html
1. Click "Search Windows", search "Edit environmental variables for your account" and click on it.
1. In the "User variables" window, select the "Path", and click "Edit...".
1. In the "Edit environmental variable" window, click "New".
1. Paste the path to the 7zip application folder, and click OK.
1. Restart Windows.
1. Check the setup by opening the "Command Prompt" and enter the command "7z". If 7z is correctly setup, output details such as the version, descriptions of commands, etc should be shown.')

# Getting Started
The package is available from the R library (i.e. CRAN) at https://cran.r-project.org/web/packages/AWAPer/index.html. To install the package from CRAN use the R command:

```R
install.packages('AWAPer')
```

Alternatively, to install the latest version:

1. Download the latest "tar.gz" (e.g. AWAPer_0.1.4.tar.gz) file from https://github.com/peterson-tim-j/AWAPer/releases.
1. Open R. 
1. Install the netCDF package using the following command: `install.packages("ncdf4")` . **Importantly** this step may require installation of netCDF software outside of R. Please read the output R console messages carefully.
1. Install the remaining required packages using the R command:
`install.packages(c("utils", "sp", "raster", "chron", "maptools", "Evapotranspiration","devtools","zoo", "methods", "xts"))`
1. Load the required packages using the R command: `library(c('Evapotranspiration', 'ncdf4', 'utils', 'raster', 'chron', 'maptools', 'sp', 'zoo', 'methods', 'xts')
1. Install the AWAPer package using the following example R command (NOTE: use the full file path to the AWAPer folder). For example on a PC `install.packages("C:\MY_FOLDER\AWAPer\AWAPer_0.1.4.tar.gz", repos = NULL, type = "source")` and for Mac `install.packages(“~/Users/MyFolder/AWAPer/AWAPer_0.1.4.tar.gz", repos = NULL, type = "source")`


# Example 4. Calculate evapotranspiration

This example calculates and plot various estimates of evaportranspiration using the netCDF data form the previous example. To do this, two catchment boundaries are loaded from the package. The figure below shows the output plot from the example. It shows 10 different esitmates of area weighted evapotranspiration from 1/1/2010 to 31/12/2010 at catchment 407214 (Victoria, Australia).

![Example 4. ET plot](https://user-images.githubusercontent.com/8623994/64103416-8d9a8700-cdb5-11e9-977a-ea8fabdcfcf5.png)

```R
# Set working directory.
setwd('~/')`

# Make two netCDF files of AWAP data.
makeNetCDF_file(ncdfFilename='AWAP_demo.nc', ncdfSolarFilename='AWAP_solar_demo.nc', 
                updateFrom='2010-1-1',updateTo='2011-12-1')

# Load example cacthment boundaries.
data("catchments")

# get ET constants
data(constants,package='Evapotranspiration')

# Extract data and esitmate various types of ET and estimate the spatial variance 
# using the interquartile range.
#----------------------------------------------
climateData.ET.HargreavesSamani = extractCatchmentData(ncdfFilename='AWAP_demo.nc',   
      ncdfSolarFilename='AWAP_solar_demo.nc', extractFrom=as.Date("2009-1-1","%Y-%m-%d"),
      extractTo=as.Date("2010-12-1","%Y-%m-%d"), locations=catchments, spatial.function.name='IQR',
      ET.function='ET.HargreavesSamani', ET.timestep = 'daily', ET.constants= constants);

climateData.ET.JensenHaise = extractCatchmentData(ncdfFilename='AWAP_demo.nc', 
      ncdfSolarFilename='AWAP_solar_demo.nc', extractFrom=as.Date("2009-1-1","%Y-%m-%d"),
      extractTo=as.Date("2010-12-1","%Y-%m-%d"), locations=catchments,  spatial.function.name='IQR', 
      ET.function='ET.JensenHaise', ET.timestep = 'daily', ET.constants= constants);

climateData.ET.Makkink = extractCatchmentData(ncdfFilename='AWAP_demo.nc',
      ncdfSolarFilename='AWAP_solar_demo.nc',extractFrom=as.Date("2009-1-1","%Y-%m-%d"),
      extractTo=as.Date("2010-12-1","%Y-%m-%d"),locations=catchments,  spatial.function.name='IQR', 
      ET.function='ET.Makkink', ET.timestep = 'daily', ET.constants= constants);

climateData.ET.McGuinnessBordne = extractCatchmentData(ncdfFilename='AWAP_demo.nc',
      ncdfSolarFilename='AWAP_solar_demo.nc', extractFrom=as.Date("2009-1-1","%Y-%m-%d"),
      extractTo=as.Date("2010-12-1","%Y-%m-%d"),locations=catchments, spatial.function.name='IQR',
      ET.function='ET.McGuinnessBordne', ET.timestep = 'daily', ET.constants= constants);

climateData.ET.MortonCRAE = extractCatchmentData(ncdfFilename='AWAP_demo.nc', 
      ncdfSolarFilename='AWAP_solar_demo.nc', extractFrom=as.Date("2009-1-1","%Y-%m-%d"),
      extractTo=as.Date("2010-12-1","%Y-%m-%d"), locations=catchments,  spatial.function.name='IQR', 
      ET.function='ET.MortonCRAE', ET.timestep = 'monthly', ET.constants= constants);

climateData.ET.MortonCRAE.potentialET = extractCatchmentData(ncdfFilename='AWAP_demo.nc',
      ncdfSolarFilename='AWAP_solar_demo.nc', extractFrom=as.Date("2009-1-1","%Y-%m-%d"),
      extractTo=as.Date("2010-12-1","%Y-%m-%d"),locations=catchments, spatial.function.name='IQR', 
      ET.function='ET.MortonCRAE', ET.timestep = 'monthly', ET.Mortons.est='potential ET', 
      ET.constants= constants);

climateData.ET.MortonCRAE.wetarealET = extractCatchmentData(ncdfFilename='AWAP_demo.nc',
      ncdfSolarFilename='AWAP_solar_demo.nc', extractFrom=as.Date("2009-1-1","%Y-%m-%d"),
      extractTo=as.Date("2010-12-1","%Y-%m-%d"), locations=catchments, spatial.function.name='IQR',
      ET.function='ET.MortonCRAE', ET.timestep = 'monthly', ET.Mortons.est='wet areal ET', 
      ET.constants= constants);

climateData.ET.MortonCRAE.actualarealET = extractCatchmentData(ncdfFilename='AWAP_demo.nc',
      ncdfSolarFilename='AWAP_solar_demo.nc',extractFrom=as.Date("2009-1-1","%Y-%m-%d"), 
      extractTo=as.Date("2010-12-1","%Y-%m-%d"), locations=catchments, spatial.function.name='IQR',
      ET.function='ET.MortonCRAE', ET.timestep = 'monthly', ET.Mortons.est='actual areal ET', 
      ET.constants= constants);

climateData.ET.MortonCRWE = extractCatchmentData(ncdfFilename='AWAP_demo.nc',
      ncdfSolarFilename='AWAP_solar_demo.nc', extractFrom=as.Date("2009-1-1","%Y-%m-%d"),
      extractTo=as.Date("2010-12-1","%Y-%m-%d"),locations=catchments, spatial.function.name='IQR',
      ET.function='ET.MortonCRWE', ET.timestep = 'monthly', ET.constants= constants);

climateData.ET.MortonCRWE.shallowLake = extractCatchmentData(ncdfFilename='AWAP_demo.nc',
      ncdfSolarFilename='AWAP_solar_demo.nc', extractFrom=as.Date("2009-1-1","%Y-%m-%d"),
      extractTo=as.Date("2010-12-1","%Y-%m-%d"), locations=catchments, spatial.function.name='IQR', 
      ET.function='ET.MortonCRWE', ET.timestep = 'monthly', ET.Mortons.est = 'shallow lake ET', 
      ET.constants= constants);

climateData.ET.Turc = extractCatchmentData(ncdfFilename='AWAP_demo.nc',
      ncdfSolarFilename='AWAP_solar_demo.nc', extractFrom=as.Date("2009-1-1","%Y-%m-%d"),
      extractTo=as.Date("2010-12-1","%Y-%m-%d"),locations=catchments, spatial.function.name='IQR',
      ET.function='ET.Turc', ET.timestep = 'daily', ET.constants= constants);

# Plot the ET estimates for one of the catchmnts.
#----------------------------------------
filt = climateData.ET.HargreavesSamani$catchmentTemporal.mean$CatchID==407214
d = ISOdate(climateData.ET.HargreavesSamani$catchmentTemporal.mean$year, 
      climateData.ET.HargreavesSamani$catchmentTemporal.mean$month, 
      climateData.ET.HargreavesSamani$catchmentTemporal.mean$day)
plot(d[filt], climateData.ET.HargreavesSamani$catchmentTemporal.mean$ET_mm[filt], 
      col='black',lty=1, xlim = c(ISOdate(2010,1,1), ISOdate(2010,12,1)), ylim=c(0, 30),type='l', ylab='ET [mm/d]',xlab='Date')

filt = climateData.ET.JensenHaise$catchmentTemporal.mean$CatchID==407214
d = ISOdate(climateData.ET.JensenHaise$catchmentTemporal.mean$year, 
      climateData.ET.JensenHaise$catchmentTemporal.mean$month, climateData.ET.JensenHaise$catchmentTemporal.mean$day)
lines(d[filt], climateData.ET.JensenHaise$catchmentTemporal.mean$ET_mm[filt], col='red',lty=1)

filt = climateData.ET.Makkink$catchmentTemporal.mean$CatchID==407214
d = ISOdate(climateData.ET.Makkink$catchmentTemporal.mean$year, climateData.ET.Makkink$catchmentTemporal.mean$month, 
      climateData.ET.Makkink$catchmentTemporal.mean$day)
lines(d[filt], climateData.ET.Makkink$catchmentTemporal.mean$ET_mm[filt], col='green',lty=1)

filt = climateData.ET.McGuinnessBordne$catchmentTemporal.mean$CatchID==407214
d = ISOdate(climateData.ET.McGuinnessBordne$catchmentTemporal.mean$year, climateData.ET.McGuinnessBordne$catchmentTemporal.mean$month,      
      climateData.ET.McGuinnessBordne$catchmentTemporal.mean$day)
lines(d[filt], climateData.ET.McGuinnessBordne$catchmentTemporal.mean$ET_mm[filt], col='blue',lty=1)

filt = climateData.ET.MortonCRAE.potentialET$catchmentTemporal.mean$CatchID==407214
d = ISOdate(climateData.ET.MortonCRAE.potentialET$catchmentTemporal.mean$year, 
      climateData.ET.MortonCRAE.potentialET$catchmentTemporal.mean$month, climateData.ET.MortonCRAE.potentialET$catchmentTemporal.mean$day)
lines(d[filt], climateData.ET.MortonCRAE.potentialET$catchmentTemporal.mean$ET_mm[filt], col='black',lty=2)

filt = climateData.ET.MortonCRAE.wetarealET$catchmentTemporal.mean$CatchID==407214
d = ISOdate(climateData.ET.MortonCRAE.wetarealET$catchmentTemporal.mean$year,     
      climateData.ET.MortonCRAE.wetarealET$catchmentTemporal.mean$month, climateData.ET.MortonCRAE.wetarealET$catchmentTemporal.mean$day)
lines(d[filt], climateData.ET.MortonCRAE.wetarealET$catchmentTemporal.mean$ET_mm[filt], col='red',lty=2)

filt = climateData.ET.MortonCRAE.actualarealET$catchmentTemporal.mean$CatchID==407214
d = ISOdate(climateData.ET.MortonCRAE.actualarealET$catchmentTemporal.mean$year, 
      climateData.ET.MortonCRAE.actualarealET$catchmentTemporal.mean$month, 
      climateData.ET.MortonCRAE.actualarealET$catchmentTemporal.mean$day)
lines(d[filt], climateData.ET.MortonCRAE.actualarealET$catchmentTemporal.mean$ET_mm[filt], col='green',lty=2)

filt = climateData.ET.MortonCRWE$catchmentTemporal.mean$CatchID==407214
d = ISOdate(climateData.ET.MortonCRWE$catchmentTemporal.mean$year, climateData.ET.MortonCRWE$catchmentTemporal.mean$month, 
      climateData.ET.MortonCRWE$catchmentTemporal.mean$day)
lines(d[filt], climateData.ET.MortonCRWE$catchmentTemporal.mean$ET_mm[filt], col='blue',lty=2)

filt = climateData.ET.MortonCRWE.shallowLake$catchmentTemporal.mean$CatchID==407214
d = ISOdate(climateData.ET.MortonCRWE.shallowLake$catchmentTemporal.mean$year, 
      climateData.ET.MortonCRWE.shallowLake$catchmentTemporal.mean$month, climateData.ET.MortonCRWE.shallowLake$catchmentTemporal.mean$day)
lines(d[filt], climateData.ET.MortonCRWE.shallowLake$catchmentTemporal.mean$ET_mm[filt], col='black',lty=3)

filt = climateData.ET.Turc$catchmentTemporal.mean$CatchID==407214
d = ISOdate(climateData.ET.Turc$catchmentTemporal.mean$year, climateData.ET.Turc$catchmentTemporal.mean$month, 
      climateData.ET.Turc$catchmentTemporal.mean$day)
lines(d[filt], climateData.ET.Turc$catchmentTemporal.mean$ET_mm[filt], col='red',lty=3)

legend(x='topright', legend=c(
  'Hargreaves Samani (ref. crop)', 'Jensen Haise (PET)', 'Makkink (ref. crop)', 'McGuinness Bordne (PET)', 'Morton CRAE (PET)',
  'Morton CRAE (wet areal ET)', 'Morton CRAE (actual areal ET)', 'Morton CRWE (PET)', 'Morton CRWE (shallowLake)', 'Turc (ref. crop, non-humid'),
  lty = c(1,1,1,1,2,2,2,2,3,3), col=c('black','red','green','blue','black','red','green','blue','black','red')
  )
  ```
# Example 5. Map of Australia
This example helps to illustrate how the netcdf file created can be interrogated to produce a map for daily maximum tempeature. 

![Example 5. Map of Australia](https://user-images.githubusercontent.com/39328041/80854309-6d8bd200-8c7a-11ea-8cac-96d14f4927c6.png)

```R
# Load the package
library(AWAPer)
library(ncdf4)
library(maps)
library(mapdata)

# Download the AWAP data to an NC file
makeNetCDF_file(ncdfFilename = "AWAP_demo.nc", ncdfSolarFilename="AWAP_solar_demo.nc", 
                updateFrom = as.Date("2017-12-14","%Y-%m-%d"),updateTo = as.Date("2017-12-15","%Y-%m-%d"))

# Connect to the NC file and list what is inside
ncfile = nc_open("AWAP_demo.nc")
names(ncfile$dim)  #  "Long" "Lat"  "time"
names(ncfile$var)  #  "tmin" "tmax" "vprp" "precip"

# Get netcdf geometry
lon = ncvar_get(ncfile, "Long")
lat = ncvar_get(ncfile, "Lat")

# Get times
print(ncfile$dim$time$units)
time = ncvar_get(ncfile, "time")
date = as.Date("1900-01-01") + time

# Time point to extract
tar.date = as.Date("2017-12-15")
tar.indx = match(tar.date, date)

# Extract an image/slice
im.tmax = ncvar_get(ncfile, varid = "tmax",   start = c(1, 1, tar.indx), count = c(-1, -1, 1)) # matrix

# Close the connection
nc_close(ncfile)

# Plot Tmax
col.vec = rev(heat.colors(10))
plot(NULL, xlim = c(110, 160), ylim = c(-45, -5), xlab = expression(degree*Longitude), ylab = expression(degree*Latitude),
     xaxs = "i", yaxs = "i", mgp = c(2, 0.6, 0))
image(lon, lat, im.tmax, zlim = c(0, 50), col = col.vec, add = TRUE)

# Remove data outisde coastal boundaries
outline = map("worldHires", regions = c("Australia"), exact = TRUE, plot = FALSE) # returns a list of x/y coords
xbox    = c(111, 156.2)
ybox    = c(-39.5, -9)
subset  = !is.na(outline$x)
polypath(c(outline$x[subset], NA, c(xbox, rev(xbox))),
         c(outline$y[subset], NA, rep(ybox, each = 2)),
         col = "white", rule = "evenodd", border = "white")

outline = map("worldHires", regions = c("Australia:Tasmania"), exact = TRUE, plot = FALSE) # returns a list of x/y coords
xbox    = c(111, 156.2)
ybox    = c(-44.9, -39.5)
subset  = !is.na(outline$x)
polypath(c(outline$x[subset], NA, c(xbox, rev(xbox))),
         c(outline$y[subset], NA, rep(ybox, each = 2)),
         col = "white", rule = "evenodd", border = "white")

# Add legend
for (i in 1:length(col.vec)) {
  rect(115+(i-1)*2, -42, 117+(i-1)*2, -42+1.4, col = col.vec[i], border = NULL)
}

text(118.0, -38.5, expression("Tmax ("*degree*"C)"), pos = 1, cex = 0.7)
text(115, -41.3, "0", pos = 1, cex = 0.7)
text(125, -41.3, "25", pos = 1, cex = 0.7)
text(135, -41.3, "50", pos = 1, cex = 0.7)

# Add bounding boxes to show AWAP data boundary
rect(min(lon), min(lat), max(lon), max(lat))
text(118.5, -7.5, "AWAP boundary", pos = 1, cex = 0.8)
  ```
# Example 6. Extract monthly precipitation
This example illustrates how to extract data at a monthly time step. Note, extracting data other than at a daily time step requires version 0.1.4 of AWAPer, which is available from https://github.com/peterson-tim-j/AWAPer/releases/tag/1.4. The image below shows the derived monthly total precipitation and the monthly spatial standard deviation.

![example6](https://user-images.githubusercontent.com/8623994/83482391-87bc0880-a4e3-11ea-97a2-d174cd68446c.png)

```R
# Load required librraies
x = c('Evapotranspiration', 'ncdf4', 'utils', 'raster', 'chron', 'maptools', 'sp', 'zoo', 'methods', 'xts')
lapply(x, library, character.only = TRUE)
        
# Set dates for building netCDFs and extracting data.
startDate = as.Date("2000-01-01","%Y-%m-%d")
endDate = as.Date("2000-12-31","%Y-%m-%d")

# Set names for netCDF files.
ncdfFilename = 'AWAPer_demo.nc'

# Build netCDF grid of ONLY precipitation and over a defined time period.
file.names = makeNetCDF_file(ncdfFilename=ncdfFilename,
                             updateFrom=startDate, updateTo=endDate,
                             urlTmin=NA, urlTmax=NA, urlVprp = NA, urlSolarrad = NA)

# Load the two example cacthment boundaries.
data("catchments")

# Extract the monthly total areal average precipitation for the two catchments..
monthlySumPrecip = extractCatchmentData(ncdfFilename=file.names$ncdfFilename,
                                         ncdfSolarFilename=file.names$ncdfSolarFilename,
                                         extractFrom=startDate, extractTo=endDate,
                                         locations=catchments,
                                         getTmin = F, getTmax = F, getVprp = F, getSolarrad = F, getET = F,
                                         temporal.timestep = 'monthly', temporal.function.name = 'sum',
                                         spatial.function.name='var')


# Extract the monthly precip. sum data and spatial variance for each of the catchments.
filt = monthlySumPrecip$catchmentTemporal.sum$CatchID == 407214
monthlySumPrecip.sum.407214 = monthlySumPrecip$catchmentTemporal.sum[filt,]
monthlySumPrecip.var.407214 = monthlySumPrecip$catchmentSpatial[filt,]
filt = monthlySumPrecip$catchmentTemporal.sum$CatchID == 407220
monthlySumPrecip.sum.407220 = monthlySumPrecip$catchmentTemporal.sum[filt,]
monthlySumPrecip.var.407220 = monthlySumPrecip$catchmentSpatial[filt,]

# Create data time series
dates2plot = ISOdate(monthlySumPrecip.sum.407220$year, monthlySumPrecip.sum.407220$month, monthlySumPrecip.sum.407220$day)

# Plot monthly sum and spatial variance in the monthly sum.
par(mfrow=c(1,2))
plot(dates2plot, monthlySumPrecip.sum.407214$precip_mm, type='l', col='red', xlab='Month', ylab='Monthly precip. [mm/month]')
lines(dates2plot, monthlySumPrecip.sum.407220$precip_mm, col='blue')
legend('topleft',legend=c('Cathcment 407214','Cathcment 407220'),col=c('red','blue'), lty=c(1,1))
plot(dates2plot, sqrt(monthlySumPrecip.var.407214$precip_mm), type='l', col='red', xlab='Month', ylab='Monthly precip. spatial standard deviation [mm/month]')
lines(dates2plot, sqrt(monthlySumPrecip.var.407220$precip_mm), col='blue')
```
# Example 7: Map monthly precipitation
This example illustrates how to map the monthly total precipitation across two catchments. Note, mapping data other than at a daily time step requires version 0.1.41 of AWAPer, which is available from https://github.com/peterson-tim-j/AWAPer/releases/tag/1.41. The image below shows maps of the monthly total precipitation..

![Monthly total precip.](https://user-images.githubusercontent.com/8623994/88503595-8dcaf300-d015-11ea-9d4a-0af5da2cba17.png)

```R
# Set dates for building netCDFs and extracting data.
startDate = as.Date("2000-01-01","%Y-%m-%d")
endDate = as.Date("2000-02-28","%Y-%m-%d")

# Set names for netCDF files.
ncdfFilename = 'AWAPer_demo.nc'
ncdfSolarFilename = 'AWAPer_solar_demo.nc'

# Remove any existing netCDF demo files.
if (file.exists(ncdfFilename))
 is.removed = file.remove(ncdfFilename)
if (file.exists(ncdfSolarFilename))
 is.removed = file.remove(ncdfSolarFilename)

# Build netCDF grids and over a defined time period.
file.names = makeNetCDF_file(ncdfFilename=ncdfFilename,
                            ncdfSolarFilename=ncdfSolarFilename,
                            updateFrom=startDate, updateTo=endDate)

# Load example cacthment boundaries.
data("catchments")

# Extract the monthly total precipitation.
monthlyPrecipData = extractCatchmentData(ncdfFilename=ncdfFilename,
                                         ncdfSolarFilename=ncdfSolarFilename,
                                         extractFrom=startDate, extractTo=endDate,
                                         locations=catchments,
                                         getTmin = F, getTmax = F, getVprp = F, getSolarrad = F, getET = F,spatial.function.name = '',
                                         temporal.timestep = 'monthly', temporal.function.name = 'sum')

# Map the monthly total repcip and overlay the catchment boundary.
v = list("sp.polygons", catchments, col = "red",first=FALSE)
spplot(monthlyPrecipData,2:ncol(monthlyPrecipData), sp.layout = list(v))
```
