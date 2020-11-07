# This example shows ho to build and extract daily precipitation and PET at points (not catchment boundaries)
#-----------------------------------------------------------------------------------------------------------

# load all necessary packages
library(AWAPer)

# Make the netCDF files of only AWAP precipitation data.
fnames = makeNetCDF_file(ncdfFilename ='AWAP_demo.nc',
                         ncdfSolarFilename ='AWAP_solar_demo.nc',
                         updateFrom=as.Date("1989-01-01","%Y-%m-%d"),
                         updateTo=as.Date("1991-01-01","%Y-%m-%d"))
                         # updateFrom=as.Date("2010-08-01","%Y-%m-%d"),
                         # updateTo=as.Date("2010-10-01","%Y-%m-%d"))

# Set coordinates to four bores locations and one rainfall gauge.
coordinates.data = data.frame( ID =c('Bore-70015656','Bore-50038039','Bore-20057861','Bore-10084446','Rain-63005'),
                               Longitude = c(131.33588, 113.066933, 143.118263, 153.551875, 149.5559),
                               Latitude =  c(-12.660622, -25.860046, -38.668506,-28.517974,-33.4289))

# Convert the points to a spatial object
sp::coordinates(coordinates.data) <- ~Longitude + Latitude

# Set projection to GDA94
sp::proj4string(coordinates.data) = '+proj=longlat +ellps=GRS80'

# Get the constants required for ET estimation.
data(constants,package='Evapotranspiration')

# Download and import the Australian 9 second DEM.
# Note, the DEM only needs be downloaded if ET is be estimated.
#DEM = getDEM()
load('DEM.RDATA')

# Extract time-series of daily precip data at all five sites
climateData.data = extractCatchmentData(ncdfFilename='AWAP_demo.nc',
                                        ncdfSolarFilename ='AWAP_solar_demo.nc',
                                        extractFrom=as.Date("1989-01-01","%Y-%m-%d"),
                                        extractTo=as.Date("1991-01-01","%Y-%m-%d"),
                                        catchments=coordinates.data,
                                        getTmin=T, getTmax=T, getVprp=T, getSolarrad=T, getET=T,
                                        DEM = DEM,
                                        ET.constants = constants)
#extractFrom=as.Date("2010-08-01","%Y-%m-%d"),
#extractTo=as.Date("2010-10-01","%Y-%m-%d"),

# Plot the daily rainfall at each site
par.default = par()
par(mfrow=c(5,1), mar = c(4,5,3,0))
for (i in 1:nrow(coordinates.data)){
  filt = climateData.data$CatchmentID.ID == coordinates.data$ID[i]

  data2plot = climateData.data$precip_mm[filt]
  names(data2plot) = paste(climateData.data$day[filt],'/',climateData.data$month[filt],sep='')
  barplot(data2plot, main=coordinates.data$ID[i], ylab='Precip [mm/d]', xlab='Date [day/month]')

}

# Plot the daily ET at each site
par.default = par()
par(mfrow=c(5,1), mar = c(4,5,3,0))
for (i in 1:nrow(coordinates.data)){
  filt = climateData.data$CatchmentID.ID == coordinates.data$ID[i]

  data2plot = climateData.data$ET_mm[filt]
  names(data2plot) = paste(climateData.data$day[filt],'/',climateData.data$month[filt],sep='')
  barplot(data2plot, main=coordinates.data$ID[i], ylab='PET [mm/d]', xlab='Date [day/month]')

}

# Rest plotting parameters
par(par.default)
