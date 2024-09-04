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
      "samples",
      "n",
      "Mean_LB"
    ),
    package = "cosimmr",
    add = FALSE
  )
}