module phys_par_core

  real(8),parameter:: tm0=5700.0
  real(8),parameter:: rho0 = 12451.0, delrho_c = 580.0, gamma = 1.5
  real(8),parameter:: cp_c = 750.0, dtldp = 9.0e-9, dtldx = -2.1e4, K0 = 1.403e12
  real(8),parameter:: delta_s = 1.27e+2, x0 = 0.056
  real(8),parameter:: l_rho = 8.039e6
  real(8),parameter:: pi=4.0*atan(1.0)
  real(8),parameter:: kc0=112.0, alpha_c=1.0e-5
  real(8),parameter:: g0=9.8, gc = 6.674e-11
  real(8),parameter:: beta=0.83
  real(8),parameter:: a_rho = 0.484, a_k = 2.39

end module phys_par_core
