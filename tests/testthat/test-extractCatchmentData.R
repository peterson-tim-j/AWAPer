#Test the extraction of data
test_that("netCDF grid can be created",
    {

      Sys.setenv(R_TESTS="")

      expect_no_error(
        {
          # Set dates for building netCDFs and extracting data from yesterday to one week ago.
          startDate = as.Date("2010-08-01","%Y-%m-%d")
          endDate = as.Date("2010-09-30","%Y-%m-%d")

          # define temp direcory for netCDF files
          fdir = tempdir()
          setwd(fdir)

          # Set names for netCDF files (in the system temp. directory).
          ncdfFilename = tempfile(fileext = '.nc')
          ncdfSolarFilename = tempfile(fileext = '.nc')

          # Build netCDF grids for all data but only over the defined time period.
          file.names = makeNetCDF_file(ncdfFilename=ncdfFilename,
                                       ncdfSolarFilename=ncdfSolarFilename,
                                       updateFrom=startDate, updateTo=endDate)
        },
        message='Testing creaion of two month netCDF grids.'
      )

      expect_no_error(
        {
          # Load example catchment boundaries.
          data("catchments")

          # Extract catchment average monthly data P for Bet Bet Creek.
          climateData.P= extractCatchmentData(ncdfFilename=ncdfFilename,
                                            ncdfSolarFilename=ncdfSolarFilename,
                                            extractFrom=startDate, extractTo=endDate,
                                            locations=catchments,
                                            getTmin = F, getTmax = F, getVprp = F, getSolarrad = F, getET = F,
                                            temporal.timestep = 'monthly', temporal.function.name = 'sum',
                                            spatial.function.name='var');
        },
        message='Testing extraction of P monthly data.'
      )

      expect_no_error(
        {
          # Load the ET constants
          data(constants,package='Evapotranspiration')

          # Extract catchment average data for Bet Bet Creek with
          # the Mortons CRAE estimate of potential ET.
          climateData.P_PET= extractCatchmentData(ncdfFilename=ncdfFilename,
                                                               ncdfSolarFilename=ncdfSolarFilename,
                                                               extractFrom=startDate, extractTo=endDate,
                                                               locations=catchments,
                                                               temporal.timestep = 'monthly', temporal.function.name = 'sum',
                                                               spatial.function.name='var',
                                                               ET.function='ET.MortonCRAE',
                                                               ET.timestep='monthly', ET.constants=constants);
        },
        message='Testing extraction of P monthly and PET data.'
      )

      # Check outputs are data frames
      expect_type(climateData.P, 'list')
      expect_type(climateData.P_PET, 'list')

      # Check dimensions of outputs
      expect_true(ncol(climateData.P$catchmentTemporal.sum) == 5, 'Test expected number of columns in precip. data')
      expect_true(ncol(climateData.P_PET$catchmentTemporal.sum) == 11, 'Test expected number of columns in precip. and PET area weighted data')

      expect_true(nrow(climateData.P$catchmentTemporal.sum) == 4, 'Test expected number of rows in precip. point data')
      expect_true(nrow(climateData.P_PET$catchmentTemporal.sum) == 4, 'Test expected number of rows in precip. and PET area weighted data')

      # check data is finite
      expect_true(all(is.finite(climateData.P$catchmentTemporal.sum[,5])), 'Test precip results are finite')
      expect_true(all(is.finite(climateData.P_PET$catchmentTemporal.sum[,11])), 'Test PET results are finite')
    }
)
