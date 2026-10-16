# internal function managing download of files. Different approach is used depending on the OS.
# In Windows wininet is used by default but does not allow download longer than 30 seconds.
# In order to remove this timeout problem the package check if libcurl is available. If available
# curl package is used. Otherwise wininet is used and can potentially create timeout.
# in Mac and Linux libcurl is used.
#' @import curl
#' @noMd
#' @noRd

bgee_download_file <- function(url, destfile, quiet = TRUE, mode = NULL) {
  if(.Platform$OS.type == "windows") {
    if(capabilities("libcurl")) {
      if(is.null(mode)) {
        success <- curl_download(url = url, destfile = destfile, quiet = quiet)
      } else{
        success <- curl_download(url = url, destfile = destfile, quiet = quiet, mode = mode)
      }
      if(gsub("\\", "/", success, fixed = TRUE) == destfile | success == destfile) {
        return(0)
      } else {
        return(success)
      }
    } else {
      if(is.null(mode)) {
        warning("Windows OS without compatibility with libcurl. Could result to problems to download files.")
        warning("Please read the installation section of the vignette to solve any potential error.")
        return(download.file(url = url, destfile = destfile, quiet = quiet))
      } else{
        return(download.file(url = url, destfile = destfile, quiet = quiet, mode = mode))
      }
    }
  } else{
    if(is.null(mode)) {
      return(download.file(url = url, destfile = destfile, quiet = quiet))
    } else{
      return(download.file(url = url, destfile = destfile, quiet = quiet, mode = mode))
    }
  }
}

count_time <- function(seconds) {
  if (seconds < 60) {
    return(paste(seconds, "second(s)"))
  } else {
    minutes <- floor(seconds / 60)
    remaining_seconds <- seconds %% 60
    
    if (remaining_seconds == 0) {
      return(paste(minutes, "minute(s)"))
    } else {
      return(paste(minutes, "minute(s) and", remaining_seconds, "second(s)"))
    }
  }
}

## The end of a compressed response of the Bgee API, from Bgee 16: a gzip member holding the 5 blank
## lines ending every response, written after the member of the rows. Decompressed whole, a compressed
## response therefore ends with the same blank lines as an uncompressed one.
COMPRESSED_RESPONSE_END <- as.raw(c(0x1f, 0x8b, 0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x04, 0x03,
                                    0xe3, 0xe2, 0x02, 0x02, 0x00, 0x9a, 0x8b, 0x82, 0xf4, 0x05,
                                    0x00, 0x00, 0x00))

## Whether a response of the Bgee API was downloaded completely: every response ends with 5 blank
## lines, in a gzip member of their own when it is compressed. Only the last bytes of the file are
## read, a response can hold tens of millions of rows. A compressed response that does not end
## with that member is decompressed whole: it is truncated, or compressed otherwise.
#' @noRd
.isCompleteResponse <- function(path) {
  size <- file.size(path)
  if (is.na(size) || size == 0) {
    return(FALSE)
  }
  con <- file(path, "rb")
  on.exit(close(con))
  if (!identical(readBin(con, "raw", 2), COMPRESSED_RESPONSE_END[1:2])) {
    seek(con, max(0, size - 64))
    return(.endsWithBlankLines(readBin(con, "raw", 64)))
  }
  if (size >= length(COMPRESSED_RESPONSE_END)) {
    seek(con, size - length(COMPRESSED_RESPONSE_END))
    if (identical(readBin(con, "raw", length(COMPRESSED_RESPONSE_END)), COMPRESSED_RESPONSE_END)) {
      return(TRUE)
    }
  }
  gz <- gzfile(path, "rb")
  on.exit(close(gz), add = TRUE)
  last <- raw(0)
  ## A truncated compressed file makes the reading warn
  suppressWarnings(repeat {
    chunk <- readBin(gz, "raw", 1e7)
    if (length(chunk) == 0) break
    last <- tail(c(last, chunk), 64)
  })
  .endsWithBlankLines(last)
}

## Whether bytes end with the line ending the last row, then 5 empty lines, the lines ending with
## LF or CRLF
#' @noRd
.endsWithBlankLines <- function(bytes) {
  bytes <- bytes[bytes != as.raw(13)]
  length(bytes) >= 6 && all(tail(bytes, 6) == as.raw(10))
}
