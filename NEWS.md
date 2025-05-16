## nlsic 1.1.0 (2025-05-14)

 * lsie_ln(): added mnorm and x0 parameters
 * lsi_ln(): rewritten from scratch to use lsie_ln() recursively
 * added pnull() and ls_0() functions
 * ls_ln() now gets mnorm and x0 parameters
 * complited unit tests
 * warning() is called in problematic situations. So user can intercept
   and debug them

## nlsic 1.0.4 (2023-06-26)

 * nlsic(): fixed error reporting when step has NA
 * nlsic(): fixed error reporting when inequalities are not satisfied

## nlsic 1.0.3 (2022-04-29)

 * nlsic(): fixed error when par is a matrix and u=NULL
 * nlsic(): added `paro` field in returned result with original structure of `par`

## nlsic 1.0.2 (2022-04-12)

 * now, lsi() returns NA vector if rhs is all NA.

## nlsic 1.0.1 (2020-01-10)

 * Added a `NEWS.md` file to track changes to the package.
 * minor fixes in documentation
