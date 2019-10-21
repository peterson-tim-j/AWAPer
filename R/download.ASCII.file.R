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

  didFail = 1
  url = paste(url.string,datestring, datestring,'.grid.Z',sep='')

  OS <- Sys.info()
  OS <- OS[1]
  if (OS=='Windows') {
    des.file.name = file.path(workingFolder,paste(data.type.label,datestring,'.grid.Z',sep=''))
    didFail = tryCatch({utils::download.file(url,des.file.name, quiet = T, mode = "wb")},error = function(cond) {return(TRUE)})
    if (didFail==0) {
      didFail = tryCatch(
        {system(paste0('7z e -aoa -bso0 "',des.file.name, '"'))},
        error = function(cond) {
          message(paste(
            'The program "7z" is either not installed or cannot be found. If not installed then',
            'install it from https://www.7-zip.org/download.html .',
            'Once installed, do the following step:'))
          message('  1. Click "Search Windows", search "Edit environmental variables for your account" and click on it.')
          message('  2. In the "User variables" window, select the "Path", and click "Edit..."')
          message('  3. In the "Edit environmental variable" window, click "New", paste the path to the 7zip application folder, and click OK.')
          message('  4. Restart Windows.')
          message('  5. Open the "Command Prompt", and enter the command "7z". If setup correctly, this should output details such as the version, descriptions of commands, etc.')
          return(TRUE)
        }
      )
      des.file.name = gsub('.Z', '', des.file.name)
    } else {
      message(paste('WARNING: Downloading the following grid failed:',des.file.name,'. Please check the URL, internet connection.'))




    }
  } else {

    des.file.name = file.path(workingFolder,paste(data.type.label,datestring,'.grid.Z',sep=''))
    didFail = tryCatch({utils::download.file(url,des.file.name, quiet = T)},error = function(cond) {return(TRUE)})
    if (didFail==0) {
      system(paste('znew -f ',des.file.name));
      des.file.name = gsub('.Z', '.gz', des.file.name)
    } else {
      message(paste('WARNING: Downloading the following grid failed:',des.file.name,'. Please check the URL and the internet connection.'))
    }
  }

  return(list(file.name=des.file.name,didFail=didFail))
}
