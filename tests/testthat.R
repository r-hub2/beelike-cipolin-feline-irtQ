# Limit BLAS/OpenMP thread usage before loading testthat/irtQ, so that the
# environment variables take effect before the underlying BLAS/OpenMP
# libraries are loaded into the session.
Sys.setenv(OMP_NUM_THREADS = "1")
Sys.setenv(OMP_THREAD_LIMIT = "1")
Sys.setenv(OPENBLAS_NUM_THREADS = "1")
Sys.setenv(MKL_NUM_THREADS = "1")
Sys.setenv(BLIS_NUM_THREADS = "1")
Sys.setenv(VECLIB_MAXIMUM_THREADS = "1")

library(testthat)
library(irtQ)

test_check("irtQ")
