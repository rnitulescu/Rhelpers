options(mc.cores=parallel::detectCores())
rstan_options(auto_write=TRUE)
Sys.setenv(LOCAL_CPPFLAGS='-march=native')

