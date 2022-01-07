if (requireNamespace("RUnit", quietly=TRUE)) {
   library(RUnit)
   library(nlsic)

   testSuite <- defineTestSuite(
      name = "nlsic unit tests",
      dirs = system.file("unitTests", package = "nlsic"),
      testFuncRegexp = "^[Tt]est.+"
   )
   Sys.setenv("R_TESTS"="")
   tests <- runTestSuite(testSuite)

   printTextProtocol(tests)
   err=getErrors(tests)
   if (err$nFail > 0) stop("RUnit: ", err$nFail, " test failure(s)")
   if (err$nErr > 0) stop("RUnit: ", err$nErr, " error(s) in RUnit tests")
}
