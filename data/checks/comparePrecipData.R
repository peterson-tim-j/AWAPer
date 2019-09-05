# Filter for Campaspe River at Ashborne (406208)
filt=climateDataAvg$CatchID==406208
data=climateDataAvg[filt,]

# Plot comparison against Woodend
#----------------------------
# Read in guaged rainfall at Woodend.
guage.data = read.csv('/home/WATER/staff/timjp/Research/SRC/GIT-sandbox/AWAPer/data/Woodend_088061.csv')

# Filter to concurrent days
data.datetime = ISOdate(data$year, data$month, data$day)
guage.data.datetime = ISOdate(guage.data$Year, guage.data$Month, guage.data$Day)

startDate = max(min(data.datetime),min(guage.data.datetime))
endDate = min(max(data.datetime),max(guage.data.datetime))

filt = data.datetime>=startDate & data.datetime<=endDate
data = data[filt,]

filt = guage.data.datetime>=startDate & guage.data.datetime<=endDate
guage.data = guage.data[filt,]

# Check they are the same length.
message(paste('AWAP data length =',nrow(data)))
message(paste('Gauge data length =',nrow(guage.data)))
message(paste('AWAP No. NAs =',sum(is.na(data[,5])) ))
message(paste('Gauge data length =',sum(is.na(guage.data[,6])) ))

# Plot gauge vs AWAP rainfall
png('/home/WATER/staff/timjp/Research/SRC/GIT-sandbox/AWAPer/data/Woodend_088061_precip.png')
plot(guage.data[,6], data[,5],xlab='Rainfall @ Wooded (88061) [mm/d]',ylab='Catchment average rainfall @ Ashborne (406208) [mm/d]')
abline(0,1,col='grey')
dev.off()

df = data.frame(gauge=guage.data[,6], awap=data[,5])
f = lm(gauge ~ awap, df)

# Plot gauge vs AWAP rainfall
png('/home/WATER/staff/timjp/Research/SRC/GIT-sandbox/AWAPer/data/Woodend_088061_cumPrecip.png')
filt = !is.na(data[,5]) & !is.na(guage.data[,6])
plot(cumsum(guage.data[filt,6]), cumsum(data[filt,5]),xlab='Cum. rainfall @ Wooded (88061) [mm]',ylab='Cum. catchment average rainfall @ Ashborne (406208) [mm]', type='l',xlim=c(0,45000),ylim=c(0,45000))
abline(0,1,col='grey')
dev.off();
#----------------------------

# Plot comparison against Trentham
#----------------------------
# Read in guaged rainfall at Woodend.
guage.data = read.csv('/home/WATER/staff/timjp/Research/SRC/GIT-sandbox/AWAPer/data/Trentham_088059.csv')

# Filter to concurrent days
data.datetime = ISOdate(data$year, data$month, data$day)
guage.data.datetime = ISOdate(guage.data$Year, guage.data$Month, guage.data$Day)

startDate = max(min(data.datetime),min(guage.data.datetime))
endDate = min(max(data.datetime),max(guage.data.datetime))

filt = data.datetime>=startDate & data.datetime<=endDate
data = data[filt,]

filt = guage.data.datetime>=startDate & guage.data.datetime<=endDate
guage.data = guage.data[filt,]

# Check they are the same length.
message(paste('AWAP data length =',nrow(data)))
message(paste('Gauge data length =',nrow(guage.data)))
message(paste('AWAP No. NAs =',sum(is.na(data[,5])) ))
message(paste('Gauge data length =',sum(is.na(guage.data[,6])) ))

# Plot gauge vs AWAP rainfall
png('/home/WATER/staff/timjp/Research/SRC/GIT-sandbox/AWAPer/data/Trentham_088059_precip.png')
plot(guage.data[,6], data[,5],xlab='Rainfall @ Trentham (088059) [mm/d]',ylab='Catchment average rainfall @ Ashborne (406208) [mm/d]')
abline(0,1,col='grey')
dev.off()

df = data.frame(gauge=guage.data[,6], awap=data[,5])
f = lm(gauge ~ awap, df)

# Plot gauge vs AWAP rainfall
png('/home/WATER/staff/timjp/Research/SRC/GIT-sandbox/AWAPer/data/Trentham_088059_cumPrecip.png')
filt = !is.na(data[,5]) & !is.na(guage.data[,6])
plot(cumsum(guage.data[filt,6]), cumsum(data[filt,5]),xlab='Cum. rainfall @ Trentham (088059) [mm]',ylab='Cum. catchment average rainfall @ Ashborne (406208) [mm]', type='l',xlim=c(0,45000),ylim=c(0,45000))
abline(0,1,col='grey')
dev.off();
#----------------------------


# Plot comparison against SILO coentroid for HRS boundary
#----------------------------
# Read in SILO centroid rainfall
guage.data = read.table('/home/WATER/staff/timjp/Research/SRC/GIT-sandbox/AWAPer/data/checks/406208_SILO_MortonPET_-374_14445.txt',header=T)

# Filter to concurrent days
data.datetime = as.Date(ISOdate(data$year, data$month, data$day))
guage.data.datetime = as.Date(guage.data$Date2,'%d-%m-%Y')

startDate = max(min(data.datetime),min(guage.data.datetime))
endDate = min(max(data.datetime),max(guage.data.datetime))

filt = data.datetime>=startDate & data.datetime<=endDate
data = data[filt,]

filt = guage.data.datetime>=startDate & guage.data.datetime<=endDate
guage.data = guage.data[filt,]

# Check they are the same length.
message(paste('AWAP data length =',nrow(data)))
message(paste('SILO data length =',nrow(guage.data)))
message(paste('AWAP No. NAs =',sum(is.na(data[,5])) ))
message(paste('Gauge data length =',sum(is.na(guage.data[,8])) ))

# Plot gauge vs AWAP rainfall
png('/home/WATER/staff/timjp/Research/SRC/GIT-sandbox/AWAPer/data/SILO_precip.png')
plot(guage.data[,8], data[,5],xlab='Rainfall @ Trentham (088059) [mm/d]',ylab='Catchment average rainfall @ Ashborne (406208) [mm/d]')
abline(0,1,col='grey')
dev.off()

df = data.frame(gauge=guage.data[,8], awap=data[,5])
f = lm(gauge ~ awap, df)

# Plot gauge vs AWAP rainfall
png('/home/WATER/staff/timjp/Research/SRC/GIT-sandbox/AWAPer/data/SILO_cumPrecip.png')
filt = !is.na(data[,5]) & !is.na(guage.data[,8])
plot(cumsum(guage.data[filt,8]), cumsum(data[filt,5]),xlab='Cum. rainfall @ Trentham (088059) [mm]',ylab='Cum. catchment average rainfall @ Ashborne (406208) [mm]', type='l',xlim=c(0,45000),ylim=c(0,45000))
abline(0,1,col='grey')
dev.off();

# PET
#---------

data.datetime = as.Date(ISOdate(data$year, data$month, data$day))
guage.data.datetime = as.Date(guage.data$Date2,'%d-%m-%Y')

startDate = as.Date("1990-1-1")
endDate = min(max(data.datetime),max(guage.data.datetime))

filt = data.datetime>=startDate & data.datetime<=endDate
data = data[filt,]

filt = guage.data.datetime>=startDate & guage.data.datetime<=endDate
guage.data = guage.data[filt,]


# Plot gauge vs AWAP PET
png('/home/WATER/staff/timjp/Research/SRC/GIT-sandbox/AWAPer/data/SILO_PET.png')
plot(guage.data[,21], data[,11],xlab='SILO Mortons PET @ Trentham >=1990 (088059) [mm/d]',ylab='Catchment average Mortons PET >=1990 @ Ashborne (406208) [mm/d]')
abline(0,1,col='grey')
dev.off()

df = data.frame(gauge=guage.data[,21], awap=data[,11])
f = lm(gauge ~ awap, df)

# Plot gauge vs AWAP PET
png('/home/WATER/staff/timjp/Research/SRC/GIT-sandbox/AWAPer/data/SILO_cumPET.png')
filt = !is.na(data[,11]) & !is.na(guage.data[,21])
plot(cumsum(guage.data[filt,21]), cumsum(data[filt,11]),xlab='Cum. SILO Mortons PET @ Trentham >=1990 (088059) [mm]',ylab='Cum. catchment average Mortons PET >=1990 @ Ashborne (406208) [mm]', type='l',xlim=c(0,45000),ylim=c(0,45000))
abline(0,1,col='grey')
dev.off();

#----------------------------


