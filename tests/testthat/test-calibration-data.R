testthat::test_that(
  "calibration.data returns the expected components",
  {
    datos <- robustICP::calibration.data(
      Npatients = 2,
      sigma = 1,
      tm = c(0, 0.1),
      curves.per.patient = 3,
      seed = 1234
    )
    
    testthat::expect_named(
      datos,
      c("Yij.todos", "tm", "sigma", "Nt", "seed")
    )
    
    testthat::expect_length(datos$Yij.todos, 2)
    testthat::expect_equal(datos$Nt, 2)
    testthat::expect_equal(
      dim(datos$Yij.todos[[1]]),
      c(3, 2)
    )
  }
)
