#' \code{readin.ASCII.file} reads in a grid data.
#'
#' This is an internal function. It reads in an ARCMAP ASCII grid file.
#'
#' @return
#' A matrix of the grid data.
#'
#' @keywords internal
#'
#' @export
readin.ASCII.file <- function(file.name, nRows, noData) {
  OS <- Sys.info()
  OS <- OS[1]
  if (OS=='Windows') {
    #system(paste('7z e -aoa -bso0 ',des.file.name))
    raw<-textConnection(readLines(a<-file(file.name)))
  } else {
    raw<-textConnection(readLines(a<-gzfile(file.name)));
  }

  AWAPgrid<- as.matrix(t(utils::read.table(raw, skip=6, nrows=nRows, na.strings=noData)))
  AWAPgrid <- AWAPgrid[,ncol(AWAPgrid):1]
  close(a)
  close(raw)

  return(AWAPgrid)
}
