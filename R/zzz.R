.onAttach <- function(libname, pkgname) {
  packageStartupMessage("\nTo cite in publications and working papers, please use:\n")
  packageStartupMessage(" Mehlhaff, Isaac D. A Group-Based Approach to Measuring Polarization. Working Paper. May 2022.\n")
}

release_questions <- function() {
  c("Have you incremented the version number?",
    "Have you formatted cran-comments.md appropriately?")
}
