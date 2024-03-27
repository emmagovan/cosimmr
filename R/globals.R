if (getRversion() >= "2.15.1") {
  utils::globalVariables(
    names = c(
      "Beta",
      "Proportion",
      "Source",
      "Type",
      "density"
    ),
    package = "cosimmr",
    add = FALSE
  )
}