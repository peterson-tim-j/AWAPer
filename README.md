# AWAPer
This R package builds netCDF files of the Bureau of Meteorology Australian Water Availability Project daily national climate grids and allows efficient extraction of daily catchment average precipitation, Tmin, Tmax, vapour pressure, solar radiation and then estimation of areal potential evaporation (Morton's). For package details see the PDF manual https://github.com/peterson-tim-j/AWAPer/blob/master/AWAPer.pdf.

# Getting Started
To get started using the package, do the following steps. Below this are two examples for building and extracting data:

1. Download the latest release of the package from https://github.com/peterson-tim-j/AWAPer/releases.
1. Unzip the package on your local machine.
1. Rename the unzipped files, say "AWAPer-1.0.zip" folder to "AWAPer".
1. Open R. 
1. Install the netCDF package using the following command: `install.packages("ncdf4")` . **Importantly** this step may require installation of netCDF software outside of R. Please read the output R console messages carefully.
1. Install the remaining required packages using the following R command:
`install.packages(c("R.utils", "sp", "raster", "chron", "maptools", "Evapotranspiration","devtools"))`
1. Install the AWAPer package using the following example R command (NOTE: use the full file path to the AWAPer folder): `install("C:\MY_FOLDER\AWAPer")` or `install.packages("C:\MY_FOLDER\AWAPer\", repos = NULL, type = "source")`

# Example 1. Build netCDF files

This example shows the steps required to build the two netCDF files, each containing 1 years of data.

```R
# Set working directory.
setwd('~/')`

# Make two netCDF files of AWAP data.
makeNetCDF_file(ncdfFilename='AWAP_demo.nc', ncdfSolarFilename='AWAP_solar_demo.nc', 
                updateFrom='2010-1-1',updateTo='2011-12-1')
```

# Example 2. Calculate evapotranspiration

This example calculates and plot various estimates of evaportrnspiration using the netCDF data form the previous example. To do this, the Austrlia 9 second DEM is downloaded and two catchment boundaries are loaded from the package. The figure below shows the output plot from the example. It shows 10 different esitmates of area weighted evapotranspiration from 1/1/2010 to 31/12/2010 at catchment 407214 (Victoria, Australia).

![image](https://user-images.githubusercontent.com/8623994/64103416-8d9a8700-cdb5-11e9-977a-ea8fabdcfcf5.png)

```R
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
