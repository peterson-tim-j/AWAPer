#' Get default URLs for loading data.
#'
#' \code{getURLs} get URLS to AWAP and Australian 9s DEM.
#'
#' This function returns a list of default URLs used to download the AWAP and DEM data.
#'
#' @return
#' A list variable of URLs as characters.
#'
#' @examples
#'
#' URLs = getURLs()
#'
#' @export
getURLs <- function() {

  URLs = list(DEM='https://datagovau.s3.amazonaws.com/bioregionalassessments/BA_ALL/ALL/DATA/Geography/Physiography/DEMGA_9second_v3/ebcf6ca2-513a-4ec7-9323-73508c5d7b93.zip',
              precip = 'http://www.bom.gov.au/web03/ncc/www/awap/rainfall/totals/daily/grid/0.05/history/nat/',
              Tmin  =  'http://www.bom.gov.au/web03/ncc/www/awap/temperature/minave/daily/grid/0.05/history/nat/',
              Tmax =   'http://www.bom.gov.au/web03/ncc/www/awap/temperature/maxave/daily/grid/0.05/history/nat/',
              vprp =   'http://www.bom.gov.au/web03/ncc/www/awap/vprp/vprph15/daily/grid/0.05/history/nat/',
              solarrad  = 'http://www.bom.gov.au/web03/ncc/www/awap/solar/solarave/daily/grid/0.05/history/nat/')

  return(URLs)

}

