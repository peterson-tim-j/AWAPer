#' \code{download.ASCII.file} downloads grid data.
#'
#' This is an internal function. It downloads an ARCMAP ASCII grid file from a URL.
#'
#' @return
#' A list variable of the file name and succes/failure flag.
#'
#' @keywords internal
#'
#' @export
download.ASCII.file <- function (url.string, data.type.label,  workingFolder, datestring) {

  if (!is.character(url.string))
    stop(paste('The input URL for',data.type.label,'must be a URL string.'))

  if (!startsWith(url.string,'http://'))
    stop(paste('The input URL string for',data.type.label,'must start "with http://" '))

  didFail = 1
  url = paste(url.string,datestring, datestring,'.grid.Z',sep='')

  OS <- Sys.info()
  OS <- OS[1]
  if (OS=='Windows') {
    des.file.name = file.path(workingFolder,paste(data.type.label,datestring,'.grid.Z',sep=''))
    didFail = tryCatch({utils::download.file(url,des.file.name, quiet = T, mode = "wb")},error = function(cond) {return(TRUE)})
    if (didFail==0) {
      exitMessage = system(paste0('7z e -aoa -bso0 "',des.file.name, '"'),intern = T)
      file.remove(des.file.name)
      if (!is.null(attr(exitMessage,'status'))) {
        message('------------------------------------------------------------------------------------')
        message('The program "7z" is either not installed or cannot be found. If not installed then')
        message('install it from https://www.7-zip.org/download.html .')
        message('Once installed, do the following step:')
        message('  1. Click "Search Windows", search "Edit environmental variables for your account" and click on it.')
        message('  2. In the "User variables" window, select the "Path", and click "Edit..."')
        message('  3. In the "Edit environmental variable" window, click "New".')
        message('  4. Paste the path to the 7zip application folder, and click OK.')
        message('  5. Restart Windows.')
        message('  6. Open the "Command Prompt" and enter the command "7z".')
        message('     If setup correctly, this should output details such as the version, descriptions of commands, etc.')
        message('------------------------------------------------------------------------------------')
        stop()
      }

      des.file.name = gsub('.Z', '', des.file.name)
    }
  } else {

    des.file.name = file.path(workingFolder,paste(data.type.label,datestring,'.grid.Z',sep=''))
    didFail = tryCatch({utils::download.file(url,des.file.name, quiet = T)},error = function(cond) {return(TRUE)})
    if (didFail==0) {
      system(paste('znew -f ',des.file.name));
      des.file.name = gsub('.Z', '.gz', des.file.name)
    }
  }

  return(list(file.name=des.file.name,didFail=didFail))
}
