if (.Platform$OS.type == "windows") {
  invisible(Sys.setlocale("LC_CTYPE", "Chinese_China.utf8"))
}
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[1] else "."
files <- list.files(root, pattern = "\\.R$", recursive = TRUE, full.names = TRUE)
bad <- character()
for (f in files) tryCatch(parse(f, encoding = "UTF-8"), error = function(e) {
  bad <<- c(bad, f); message(f, ": ", conditionMessage(e))
})
if (length(bad)) stop(length(bad), " R files failed parsing")
cat(length(files), "R files parsed; no analyses were executed.\n")
