context("Completeness of the responses of the Bgee API")

rows <- "GENE_ID\tCONDITION_ID\nG1\tUBERON:0000955\nG2\tCL:0000540-UBERON:0000955\n"
writeResponse <- function(content) {
  path <- tempfile()
  writeBin(charToRaw(content), path)
  path
}
## Compressed as the Bgee API does: the member of the rows, then the member of the blank lines
writeCompressedResponse <- function(content) {
  path <- tempfile()
  con <- gzfile(path, "wb")
  writeBin(charToRaw(content), con)
  close(con)
  con <- file(path, "ab")
  writeBin(COMPRESSED_RESPONSE_END, con)
  close(con)
  path
}

test_that("A complete response ends with 5 blank lines, compressed or not", {
  expect_true(.isCompleteResponse(writeResponse(paste0(rows, "\n\n\n\n\n"))))
  expect_true(.isCompleteResponse(writeResponse(gsub("\n", "\r\n", paste0(rows, "\n\n\n\n\n"),
                                                     fixed = TRUE))))
  compressed <- writeCompressedResponse(rows)
  expect_true(.isCompleteResponse(compressed))
  ## Decompressed whole, it holds the rows then the blank lines, as an uncompressed response
  expect_equal(readLines(gzfile(compressed)), c(strsplit(rows, "\n")[[1]], rep("", 5)))
  expect_equal(nrow(fread(compressed, sep = "\t", header = TRUE, colClasses = "character")), 2)
  ## Compressed otherwise, the blank lines in the member of the rows
  otherwise <- tempfile()
  con <- gzfile(otherwise, "wb")
  writeBin(charToRaw(paste0(rows, "\n\n\n\n\n")), con)
  close(con)
  expect_true(.isCompleteResponse(otherwise))
})

test_that("An incomplete response is detected, compressed or not", {
  expect_false(.isCompleteResponse(writeResponse(rows)))
  expect_false(.isCompleteResponse(writeResponse(paste0(rows, "\n\n\n\n"))))
  compressed <- writeCompressedResponse(paste0(rows, strrep("G3\tUBERON:0000955\n", 1000)))
  bytes <- readBin(compressed, "raw", file.size(compressed))
  truncated <- tempfile()
  writeBin(bytes[seq_len(length(bytes) - 40)], truncated)
  expect_false(.isCompleteResponse(truncated))
  empty <- tempfile()
  file.create(empty)
  expect_false(.isCompleteResponse(empty))
  expect_false(.isCompleteResponse(tempfile()))
})
