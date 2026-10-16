# Run from repository root. Installs missing direct dependencies only.
pkgs <- readLines("environment/direct-packages.txt", warn = FALSE)
pkgs <- pkgs[nzchar(pkgs) & !startsWith(pkgs, "#")]
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
remaining <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(remaining)) stop("Unavailable dependencies: ", paste(remaining, collapse = ", "))
message("All direct dependencies are available. Compare versions with r-packages.csv.")
