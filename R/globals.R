if (getRversion() >= "2.15.1") {
  utils::globalVariables(
    names = c(
      "Beta",
      "Proportion",
      "Source",
      "Type",
      "density",
      "Group",
      "cov",
      "nsd",
      "psd",
      "samples"
    ),
    package = "cosimmr",
    add = FALSE
  )
}