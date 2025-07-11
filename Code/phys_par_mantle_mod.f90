module phys_par_mantle
 
  real(8),parameter:: pi=4.0*atan(1.0)
  real(8),parameter:: Mmantle=4.06e24
  real(8),parameter:: rho_um=3300.0, g0=9.8,  & 
                          alpha=3.0e-5, cp_m=1250.0
  real(8),parameter:: k_um=4.2, k_lm=10.0
  real(8),parameter:: nu0=2.0e11, fvisc=3.0, Racr=6.6e2
  real(8),parameter:: Eeta=300.0e3, r_gas=8.314
  real(8),parameter:: qrad0=20.0e12, tau_rad=2.94e9*3.15e7
  real(8),parameter:: a=6371.0e3,b=3486.0e3,Dmantle=a-b
  real(8),parameter:: Vman = (4*pi/3)*(a**3 - b**3)
  real(8),parameter:: rm=4925.0e3, Dtra=660e3
  real(8),parameter:: f_man=((a-Dtra)**3-b**3)/(a**3-b**3) ! Volume ratio of the lower mantle
  real(8),parameter:: a_e=4.0*pi*a**2, a_c=4.0*pi*b**2
  real(8),parameter:: beta=1.0/3.0
  real(8),parameter:: kappa = 1e-6
  real(8),parameter:: f_ocean = 0.7

  real(8),parameter:: tsol0=1423.0,fv0=1.0e-4,gmelt_ad=1.0e-3,gsol=3.9e-3
  real(8),parameter:: gm_ad=1.0e-8
  real(8),parameter:: rho_melt=2700.0
  real(8),parameter:: Lmelt=320.0e3
  real(8),parameter:: Dcarbon=0.0001

end module phys_par_mantle
