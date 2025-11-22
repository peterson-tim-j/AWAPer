#Test the create of netCDF files
test_that("netCDF grid can be created",
    {

      expect_no_error(
        {
          # Set dates for building netCDFs and extracting data from yesterday to one week ago.
          startDate = Sys.Date()-9
          endDate = Sys.Date()-2

          # define temp direcory for netCDF files
          fdir = tempdir()
          setwd(fdir)

          # Set names for netCDF files (in the system temp. directory).
          ncdfFilename = file.path(fdir,'data.nc')
          ncdfSolarFilename = file.path(fdir, 'solar.nc')

          # Build netCDF grids for all data but only over the defined time period.
          file.names = makeNetCDF_file(ncdfFilename=ncdfFilename,
                                       ncdfSolarFilename=ncdfSolarFilename,
                                       updateFrom=startDate, updateTo=endDate)
        },
        message='Testing creaion of netCDF grids.'
      )


      # Test the files were created
      expect_true(file.exists(ncdfFilename),'Testing creation of non-solar netCDF file for data from one week ago')
      expect_true(file.exists(ncdfSolarFilename),'Testing creation of solar netCDF file for data from one week ago')

      # test file size is reasonably large
      expect_gt(file.info(ncdfFilename)$size,1E6,'Testing non-solar netCDF file size is > 1MB')
      expect_gt(file.info(ncdfSolarFilename)$size,1E6,'Testing solar netCDF file size is > 1MB')

      # Test the netcDF files can be opened
      expect_s3_class( ncdf4::nc_open(ncdfFilename, write=T), 'ncdf4')
      expect_s3_class( ncdf4::nc_open(ncdfSolarFilename, write=T), 'ncdf4')

      # Open the netcdf files and check the dimensions
      ncout <- ncdf4::nc_open(ncdfFilename, write=T)
      ncout.solar <- ncdf4::nc_open(ncdfSolarFilename, write=T)

      # check the dimension names
      expect_true(ncdf4::ncatt_get(ncout, "time")$long_name == "time", "Testing time dimension exists in netCDF file")
      expect_true(ncdf4::ncatt_get(ncout, "Lat")$long_name == "Lat", "Testing Lat dimension exists in netCDF file")
      expect_true(ncdf4::ncatt_get(ncout, "Lat")$long_name == "Lat", "Testing Lat dimension exists in netCDF file")
      expect_true(ncdf4::ncatt_get(ncout.solar, "time")$long_name == "time", "Testing time dimension exists in solar netCDF file")
      expect_true(ncdf4::ncatt_get(ncout.solar, "Lat")$long_name == "Lat", "Testing Lat dimension exists in solar netCDF file")
      expect_true(ncdf4::ncatt_get(ncout.solar, "Lat")$long_name == "Lat", "Testing Lat dimension exists in solar netCDF file")

      # test the netCDF files have the expected dimensions.
      expect_length(ncout$dim$Lat$vals, 691)
      expect_length(ncout$dim$Long$vals, 886)
      expect_length(ncout.solar$dim$Lat$vals, 679)
      expect_length(ncout.solar$dim$Long$vals, 839)

      # Update netcDF grids and expect no errors
      endDate = startDate
      startDate = Sys.Date()-11
      expect_no_error(
        makeNetCDF_file(ncdfFilename=ncdfFilename,
                        ncdfSolarFilename=ncdfSolarFilename,
                        updateFrom=startDate, updateTo=endDate),
        message='Testing updating of netCDF grids by two days prior'
      )

      # Delete temp files and folder
      unlink(ncdfFilename)
      unlink(ncdfSolarFilename)
      unlink(fdir, recursive = TRUE)
    }
)
