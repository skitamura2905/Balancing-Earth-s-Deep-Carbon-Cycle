program earth_evolution

  use phys_par_mantle
  use phys_par_core,only:l_rho,a_rho,gamma
  implicit none
  include 'mpif.h'
  ! mpirun -np 2 ./earth_evo  

  integer,parameter:: ntime=100000
  
  integer :: ierr, rank, nprocs
  integer:: it, isub, imonte, iter
  real(8):: time, delta_t, dtempdt
  real(8):: tm, tcmb, tum, tlm, nu_m, nu_um, nu_lm, delta_um, delta_tlm, uplate
  real(8):: tsurf, tequ
  real(8):: eta_um, eta_lm
  real(8):: qconv, qrad, qcmb, qmelt, qplume, fmelt, zmelt, e_erupt, e_erupt_plume, fplume
  real(8):: xmelt_c
  real(8):: c
  real(8):: cw_s
  real(8):: cc_um,cc_lm,cc_s,cc_total,ingas_c,outgas_c,cc_a,cc_o,cc_oc
  real(8):: meta_c
  real(8):: rate_mantle_c
  real(8):: splume, Sp, S_low_to_up
  real(8):: tsol, pre_i, phi_plume, melt_plume, pre_f
  real(8):: outgas_plume_c, xplume_c, flux_um_c
  real(8):: inter_flux_c
  real(8):: slab_c
  real(8):: part_c
  real(8):: melt_trap
  real(8):: Seff, Psat, fws
  real(8):: f_cw, f_sw, beta_cw, gamma_sw, fmeta
  real(8):: ftrap
  real(8):: random
  real(8),parameter:: albedo=0.3, sigma=5.67e-8
  real(8),parameter:: kact=0.09, krun=0.045, Wref=12e9*44/3.15e7 ! Foley (2015) !7e5*a_e/(1e9*g0)/3.15e7 !3e11/3.15e7
  real(8),parameter:: Esw=80e3, uplate_ref=5e-2/3.15e7, bas_max=4e18*44
  real(8),parameter:: k_carbon = 1e7
  integer:: imin, imax, iter_total

  character(5),external:: cnum
  character(len=20) :: imonte_str, iter_str, rank_str
  character(100):: name1, name2, name3, name4, name5

  ! MPI initialization
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)

  if (rank == 0) then
     call system('mkdir -p output')
     call system('mkdir -p output/parameter')
     call system('mkdir -p output/original')
     call system('mkdir -p output/carbon')
     call system('mkdir -p output/thermal')
  end if

  call MPI_BARRIER(MPI_COMM_WORLD, ierr)

  eta_um = exp(-(a-rm)*alpha*g0/cp_m)
  eta_lm = exp(-(b-rm)*alpha*g0/cp_m)

    name5 = 'output/iteration.dat'
    open(15,file=name5,status='unknown')
    
    iter=1
    imin = 1 + rank * (10000 / nprocs)  ! First iteration for this process
    imax = (rank+1) * (10000 / nprocs)  ! Last iteration for this process
    imonte=imin
    do while (imonte<imax+1)
      write(imonte_str, '(I0)') imonte
      write(rank_str, '(I0)') rank
      name1 = 'output/thermal/thermal_' // trim(imonte_str) // '.dat'
      name2 = 'output/carbon/carbon_' // trim(imonte_str) // '.dat'
      name3 = 'output/parameter/parameter_' // trim(imonte_str) // '.dat'
      open(10,file=name1)   
      open(12,file=name2)
      open(13,file=name3)

    ! Initial setup
    cw_s = 1.4e21
    cc_total = 8.6e20*44/12 ! Taken from Tajika and Matsui (1992)
    cc_oc = 0

    time = 0.0
    Seff = (1 + 0.4*(1 - time/(4.6e9*3.15e7)))**(-1)
    tequ = (Seff*1360*(1-albedo)/(4*sigma))**0.25 - (1360*(1-albedo)/(4*sigma))**0.25
    tsurf = 285 + 4.6*((cc_a/3.3e-4)**0.346 - 1) + 2*tequ
    tm = 3000.0
    tum = eta_um*(tm - tsurf)+tsurf ; tlm=eta_lm*tm
    tcmb = 7000.0*(1.0-(b/l_rho)**2-a_rho*(b/l_rho)**4)**(gamma)
    c = 0.0
  
    ! Initial partitioning between surface and whole mantle
    call random_number(random)
    rate_mantle_c = random

    ! Initial partitioning between upper and lower mantle
    call random_number(random)
    part_c = random

    ! Subducted carbon partitioning between upper and lower mantle
    call random_number(random)
    slab_c = random

    ! Decarbonation ratio of subducting plate at subduction zone
    call random_number(random)
    fmeta = random

    ! Degassing efficiency from mid-ocean ridge
    call random_number(random)
    e_erupt = random

    ! Conversion efficiency from CMB heat flux to plume strength
    call random_number(random)
    fplume = 0.5*random

    ! Degassing efficiency from hot spot
    call random_number(random)
    e_erupt_plume = random

    ! Exponents of continental and seafloor silicate weathering
    beta_cw = 0.5
    gamma_sw = 0.23

    ! Transporting ratio of undegassed and trapped melts in subducting lithosphere
    ftrap = 1.0

    ! Initial distribution of carbon
    cc_um = rate_mantle_c*(cc_total)*part_c/(Vman*(1-f_man))
    cc_lm = rate_mantle_c*(cc_total)*(1-part_c)/(Vman*f_man)
    cc_s = max(0.0,cc_total - Vman*((1-f_man)*cc_um + f_man*cc_lm))

    cc_o = (-(k_carbon*a_e/g0 + 44*cw_s/18 - cc_s) + sqrt((k_carbon*a_e/g0 + 44*cw_s/18 - cc_s)**2 + 4*44*cc_s*cw_s/18))/2
    cc_a = (cc_s - cc_o)*g0*1e-5/a_e

    delta_t = 4.6e9*3.15e7/real(ntime)

    write(13,'(9(1pe13.5))') rate_mantle_c, part_c, slab_c, fmeta, e_erupt, e_erupt_plume, fplume
      write(iter_str, '(I0)') iter
      name4 = 'output/original/parameter_' // trim(rank_str) // '.dat'
      open(14,file=name4)
      write(14,'(9(1pe11.3))') rate_mantle_c, part_c, slab_c, fmeta, e_erupt, e_erupt_plume, fplume

    do it = 1,ntime
        ! Host star evolution and Surface temperature
        Seff = (1 + 0.4*(1 - time/(4.6e9*3.15e7)))**(-1)
        tequ = (Seff*1360*(1-albedo)/(4*sigma))**0.25 - (1360*(1-albedo)/(4*sigma))**0.25
        tsurf = 285 + 4.6*((cc_a/3.3e-4)**0.346 - 1) + 2*tequ

        ! Viscosity
        nu_m = nu0*exp(Eeta/(r_gas*tm))
        nu_um = 0.1*nu_m
        nu_lm = 0.1*fvisc*nu_m

        ! Heat flow - Plate tectonics regime
        qconv = a_e*k_um*(alpha*g0/(racr*kappa))**(beta)* &
                   ((eta_um*(tm-tsurf))**(beta+1))*(nu_um)**(-beta)

        ! Thickness of oceanic lithosphere
        delta_um = a_e*k_um*eta_um*(tm-tsurf)/qconv
        uplate = 5.38*kappa*Dmantle/delta_um**2 ! Speed of subduction

        ! Heat flow across the CMB
        qcmb = a_c*k_lm*(alpha*g0/(racr*kappa))**(beta)* &
               ((tcmb-tlm)**(beta+1))*(nu_lm)**(-beta)

        ! Radioactive heating
        qrad = qrad0*exp((4.6e9*3.15e7-time)/tau_rad)

        ! Plume heat flux
        qplume = qcmb*fplume
        delta_tlm = tm*tm*r_gas/Eeta

        ! Mass flux of plume
        splume = qplume/(delta_tlm*cp_m*rho_um) ! Volume flux
        pre_i = (tum + delta_tlm - tsol0)/(1.20e-7 - gm_ad)
        pre_f = 0

        ! Melt production in plume
        phi_plume = 0.5*(pre_i-pre_f)*1.5e-10
        if (tum+delta_tlm<tsol .OR. phi_plume<0.0) phi_plume=0.0
        if (phi_plume>1.0) phi_plume=1.0
        melt_plume = phi_plume*splume

        ! Carbon content in melt 
        xplume_c = cc_lm*((1.0-(1.0-phi_plume)**(1.0/Dcarbon))/phi_plume)
        if (phi_plume==0) xplume_c=0

        ! Carbon content in produced melts in plume
        outgas_plume_c = melt_plume*xplume_c
        flux_um_c = splume*cc_lm - outgas_plume_c

        ! Heat transport in melt migration
        call compute_qmelt(qmelt,tum,tsurf,delta_um,cc_um,fmelt,xmelt_c,e_erupt)

        ! Water and carbon fluxes
        call compute_volatile_c(outgas_c,ingas_c,meta_c,Sp,qconv,tsurf,delta_um,fmelt,xmelt_c,cc_um,cc_oc,cc_s,fmeta,e_erupt)

        ! Continental weathering
        Psat = exp(-18*2469/r_gas*(1/tsurf - 1/273.15))/exp(-18*2469/r_gas*(1/288.0 - 1/273.15))
        fws = a_e*0.3*0.08*2500*10e-3*44e-3/32e-3/3.15e7
        f_cw = fws*(1 - exp(-Wref/fws*(cc_a/3.3e-4)**beta_cw*Psat**0.3*exp(42e3/r_gas*(1.0/288 - 1/tsurf))))

        ! Seafloor weathering
        f_sw = (uplate/uplate_ref)*1.75e9*44*(cc_a/3.3e-4)**gamma_sw/3.15e7 ! Foley (2015)

        ! Temperature dependent seafloor weathering (Turn on this part) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !f_sw = (uplate/uplate_ref)*1.75e9*44*(cc_a/3.3e-4)**gamma_sw*exp(Esw/r_gas*(1.0/288 - 1/tsurf))/3.15e7 ! Foley (2015)
        !if (cc_bas + (f_sw - subd_bas)*delta_t>=bas_max) then
        !  f_sw = (bas_max - cc_bas)/delta_t + subd_bas
        !elseif (cc_bas==bas_max .AND. f_sw>subd_bas) then
        !  f_sw = subd_bas
        !endif
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (cc_a<=0) f_sw = 0
        
        ! Volume flux difference between slab and plume
        Sp = Sp*delta_um
        S_low_to_up = Sp - splume

        if (S_low_to_up>=0) then
          inter_flux_c = S_low_to_up*cc_lm
        else
          inter_flux_c = S_low_to_up*cc_um
        endif

        ! Mantle temperature change
        dtempdt = (qcmb - qconv - qmelt + qrad - qplume)/(Mmantle*cp_m)
        
        ! Mantle temperature - Averaged temperature
        tm = tm+dtempdt*delta_t

        ! Upper and Lower mantle temperature
        tum = eta_um*(tm - tsurf) + tsurf; tlm=eta_lm*tm

        ! Undegassed and trapped melts in oceanic lithosphere
        melt_trap = (1-e_erupt)*outgas_c/e_erupt + (1-e_erupt_plume)*outgas_plume_c

        ! Carbon content
        cc_um = max(0.0, cc_um + delta_t*(slab_c*ingas_c + inter_flux_c + flux_um_c - outgas_c/e_erupt + (1-ftrap)*melt_trap)/(Vman*(1-f_man)))
        cc_lm = max(0.0, cc_lm + delta_t*((1-slab_c)*ingas_c - outgas_plume_c - flux_um_c - inter_flux_c + ftrap*melt_trap)/(Vman*f_man))
        cc_oc = max(0.0, cc_oc + (0.5*f_cw + f_sw - ingas_c - meta_c)*delta_t)
        cc_s = max(0.0,cc_total - Vman*((1-f_man)*cc_um + f_man*cc_lm) - cc_oc)

        cc_o = (-(k_carbon*a_e/g0 + 44*cw_s/18 - cc_s) + sqrt((k_carbon*a_e/g0 + 44*cw_s/18 - cc_s)**2 + 4*44*cc_s*cw_s/18))/2
        cc_a = (cc_s - cc_o)*g0*1e-5/a_e

        ! Temperature at CMB
        call core_evolution(qcmb,tcmb,delta_t,c)

      !    Output to file

      if (mod(it,10)==0) then
        
      write(10,'(17(1pe13.5))') time/3.15e16, tm, tum, tlm, tcmb, delta_tlm, tsurf, c/1.0e3, & 
        qconv/1.0e12, qcmb/1.0e12, qmelt/1.0e12, qrad/1.0e12, qplume/1.0e12, & 
        uplate, Sp, splume

      write(12,'(20(1pe13.5))') time/3.15e16, &
        cc_um*Vman*(1-f_man), cc_lm*Vman*f_man, cc_s, cc_a, cc_o, cc_oc, &
        slab_c*ingas_c*3.15e7, (1-slab_c)*ingas_c*3.15e7, meta_c*3.15e7, outgas_c*3.15e7, e_erupt_plume*outgas_plume_c*3.15e7, &
        inter_flux_c*3.15e7, flux_um_c*3.15e7, f_cw*3.15e7, f_sw*3.15e7, &
        (1-e_erupt)*outgas_c/e_erupt*3.15e7, (1-e_erupt_plume)*outgas_plume_c*3.15e7, 0.25*melt_trap*3.15e7

      endif

      time = time + delta_t
    end do

    ! Constraints from Isson et al. (2020)
    if (outgas_c*3.15e7 < 5.0 * 44 * 1e9 .AND. outgas_c*3.15e7 > 1.0 * 44 * 1e9 .AND. &
        meta_c*3.15e7 < 11.5 * 44 * 1e9 .AND. meta_c*3.15e7 > 2.3 * 44 * 1e9 .AND. &
        e_erupt_plume*outgas_plume_c*3.15e7 < 3.0 * 44 * 1e9 .AND. e_erupt_plume*outgas_plume_c*3.15e7 > 0.12 * 44 * 1e9 .AND. &
        tsurf < 288*1.01 .AND. tsurf > 288*0.99) then

      close(10)
      close(12)
      close(13)
      imonte=imonte+1
    else
      close(10,status='delete')
      close(12,status='delete')
      close(13,status='delete')
    endif
    iter=iter+1
  enddo

  call MPI_Reduce(iter, iter_total, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

  if (rank == 0) then
    write(15, '(I13)') iter_total
    close(15)
  endif
  close(14)

    ! MPI finish
    call MPI_Finalize(ierr)

end program earth_evolution

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine compute_qmelt(qmelt,tum,tsurf,delta_um,cc,fmelt,xmelt_c,e_erupt)

! Heat transport by the melt migration
! See Fraeman and Korenaga (2010) and Driscoll and Bercovici (2014)

  use phys_par_mantle
  implicit none
  real(8),intent(in):: tum, tsurf, delta_um, cc, e_erupt
  real(8),intent(out):: qmelt, xmelt_c

  real(8):: tmelt, delta_tmelt, fmelt, vdot_up, fm_vol, zm
  real(8):: pre_i,uplate,phi_melt,pre_f
  real(8):: Lplate

  ! Plate velocity
  uplate = 5.38*kappa*Dmantle/delta_um**2

  ! Length of subduction
  Lplate = 1.5*2.0*pi*a 

  ! Melting related pressure and depth
  pre_i = (tum - tsol0)/(1.20e-7 - gm_ad)
  pre_f = 0
  zm = pre_i/(rho_um*g0)
  if (zm>660e3) pre_i = rho_um*g0*660e3
  if (zm>660e3) zm = 660e3

  ! Melt fraction
  phi_melt = 0.5*pre_i*1.5e-10
  if (tum<tsol0) phi_melt=0.0

  tmelt = 0.5*(tum + tsol0)
  delta_tmelt = tmelt - tsurf - 0.5*zm*gmelt_ad

  if (phi_melt>0.0) then
    fmelt = 2*zm*Lplate*uplate*phi_melt
  else
    fmelt = 0.0
  end if

  qmelt = e_erupt*fmelt*rho_melt*(Lmelt + cp_m*delta_tmelt)

! Carbon content in melt 
  xmelt_c = cc*((1.0 - (1.0 - phi_melt)**(1.0/Dcarbon))/phi_melt)

  if (fmelt==0.0) xmelt_c=0.0

end subroutine compute_qmelt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine compute_volatile_c(outgas,ingas,metamorph,Sp,qconv,tsurf,delta_um,fmelt,xmelt,cc_m,cc_oc,cc_s,fmeta,e_erupt)

  use phys_par_mantle
  implicit none
  real(8),intent(out):: outgas, ingas, metamorph, Sp
  real(8),intent(in):: qconv, tsurf, delta_um, fmelt, xmelt, cc_m, cc_oc, cc_s, fmeta, e_erupt

  real(8):: Lplate, uplate
  real(8):: eff_length = 1.5

!  Subduction dynamics
  uplate = 5.38*kappa*Dmantle/delta_um**2 ! Speed of subduction
  Lplate = 2.0*pi*a ! Length of subduction
  Sp = eff_length*Lplate*uplate ! Subduction rate

!  Outgas by the volanic eruptions
  outgas = e_erupt*fmelt*xmelt

!  Ingas by the plate subduction
  ingas = cc_oc/(f_ocean*a_e)*Sp*(1 - fmeta)

!  Metamorphic effect 
  metamorph = max(0.0,fmeta*ingas/(1 - fmeta))
  if (fmeta==1.0) metamorph = cc_oc/(f_ocean*a_e)*Sp

end subroutine compute_volatile_c

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine core_evolution(qcmb,tcmb,dt,c)

  ! See Takehiro and Sasaki (2018) and Nakagawa et al. (2023)
  ! Based on Labrosse (2015)

  use phys_par_core
  use phys_par_mantle,only: b
  implicit none

  integer:: ic_onset
  real(8):: qcmb, tcmb, dt, c
  real(8):: t_core, qs, rho_ic
  real(8):: pc, pl, pg
  real(8):: qc, ql, qg, qth_conv
  real(8):: chem_pot, ent_flux_s, ent_flux_c
  real(8):: dcdt,dtcdt,dtldric
  real(8):: dXdt
  real(8):: dtadr, flux_th, flux_comp
  real(8):: b1, b2, b4, c1, c2, c4
  real(8),parameter:: zero=0.0
  real(8),external:: fc, fx

  b1 = b/l_rho; b2 = b1*b1; b4 = b2*b2
  c1 = c/l_rho; c2 = c1*c1; c4 = c2*c2

  rho_ic=(1.0-c2-a_rho*c4)
  ic_onset = 0
  if (c>0.0) ic_onset = 1
  ent_flux_c=0.0

  if (ic_onset==0) then
    t_core=tcmb/((1.0-b2-a_rho*b4)**(gamma))
    pc = -4.0/3.0*pi*rho0*cp_c*l_rho**3*fc(b1,gamma)
    dtcdt = qcmb/pc
    t_core = t_core + qcmb/pc*dt
    tcmb = t_core*(1.0 - b2 - a_rho*b4)**(gamma)
    dcdt=0.0
    if (t_core<tm0) then
      c = l_rho*sqrt((tm0 - t_core)/(K0*dtldp))
    else
      c = 0.0
    end if
  else if (ic_onset==1) then
    dtldric=-2.0*c/l_rho**2*K0*dtldp+ &
             3.0*dtldx*(x0/fc(b1,zero))*c**2/l_rho**3
    pc = -4.0/3.0*pi*rho0*cp_c*l_rho**3*rho_ic**(-gamma)* &
        (dtldric + 2.0*gamma*t_core*c/l_rho**2* &
           ((1.0 + 2.0*a_rho*(c/l_rho)**2)/rho_ic))*(fc(b1,gamma)-fc(c1,gamma))
    pl = 4.0*pi*c**2*rho0*(1 - c2 - a_rho*c4)*t_core*delta_s
    pg = 8.0*pi**2*x0*gc*rho0**2*beta*c**2*l_rho**2/fc(b1,zero)* &
          (fx(b1,c) - fx(c1,c))

    dcdt = qcmb/(pc+pl+pg)
    dtcdt=dtldric*dcdt
    dXdt = (4.0*pi*c**2*rho0*x0*dcdt)/ &
           (4.0/3.0*pi*rho0*l_rho**3*(fc(b1,zero)-fc(c1,zero)))
    ent_flux_c=4.0*pi*c**2*rho0*x0*dcdt+ &
       dXdt*4.0/3.0*pi*rho0*l_rho**3*(fc(b1,zero)-fc(c1,zero))
    c = c + dcdt*dt
    t_core = tm0-K0*dtldp*(c/l_rho)**2+dtldx*(x0/fc(b1,zero))*(c/l_rho)**3
    tcmb=t_core*((1.0-b2-a_rho*b4)/(1.0-c2-a_rho*c4))**(gamma)
  end if

  dtadr = -2.0*b/l_rho**2*gamma*t_core*rho_ic**(-gamma)* &
                       (1.0+2.0*a_rho*(b/l_rho)**2)* &
                       (1.0-(b/l_rho)**2-a_rho*(b/l_rho)**4)**(gamma-1.0)
  qs = -4.0*pi*b**2*kc0*(1.0-a_k*(b/l_rho)**2)*dtadr

  qc = -4.0/3.0*pi*rho0*cp_c*l_rho**3*dtcdt*(fc(b1,gamma)-fc(c1,gamma))
  ql = 4.0*pi*c**2*rho0*(1-c2-a_rho*c4)*t_core*delta_s*dcdt
  qg =  8.0*pi**2*x0*gc*rho0**2*beta*c**2*l_rho**2/fc(b/l_rho,zero)* &
                   (fx(b1,c)-fx(c1,c))*dcdt

  qth_conv = qc + ql + qg - qs

end subroutine core_evolution

real(8) function fc(x,g)

  use phys_par_core
  real(8),intent(in):: x,g

  fc=x**3*(1.0-0.6*(1.0+g)*x**2-3.0/14.0*(1.0+g)*(2.0*a_rho-g)*x**4)

end function fc

real(8) function fx(x,c)

  use phys_par_core
  real(8),intent(in):: x,c

  fx=x**3*(-1.0/3.0*(c/l_rho)**2+0.2*(1.0+(c/l_rho)**2)*x**2-13.0/70.0*x**4)

end function fx
