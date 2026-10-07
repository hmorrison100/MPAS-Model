! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================

#ifdef MPAS_OPENACC
#define MPAS_ACC_TIMER_START(X) call mpas_timer_start(X)
#define MPAS_ACC_TIMER_STOP(X) call mpas_timer_stop(X)
#else
#define MPAS_ACC_TIMER_START(X)
#define MPAS_ACC_TIMER_STOP(X)
#endif

 module mpas_atmphys_interface
 use mpas_kind_types
 use mpas_pool_routines

 use mpas_atmphys_constants
 use mpas_atmphys_vars
 use mpas_timer, only: mpas_timer_start, mpas_timer_stop
 use mpas_atmphys_utilities, only: physics_error_fatal
 use mpas_atmphys_p3_pblmix, only: allocate_p3pbl,deallocate_p3pbl,p3pbl_from_MPAS

 implicit none
 private
 public:: allocate_forall_physics,   &
          deallocate_forall_physics, &
          pre_microphysics, &
          post_microphysics, &
          MPAS_to_physics,           &
          microphysics_from_MPAS,    &
          microphysics_to_MPAS,      &
          p3_allocate_work

!Interface for conversion between variables used in the MPAS dynamical core and variables needed in the
!physics parameterizations.
!Laura D. Fowler (send comments to laura@ucar.edu).
!2013-05-01.
!
! subroutines in mpas_atmphys_interface_nhyd:
! -------------------------------------------
! allocate_forall_physics  : allocation of all meteorological variables needed in the physics.
! deallocate_forall_physics: deallocation of all meteorological variables needed in the physics.
! MPAS_to_physics          : conversion of input "MPAS" variables to "physics" variables.
! microphysics_from_MPAS   : initialize local arrays needed in cloud microphysics schemes.
! microphysics_to_MPAS     : copy local arrays needed in cloud microphysics schemesto MPAS arrays.
!
! add-ons and modifications to sourcecode:
! ----------------------------------------
! * in subroutine MPAS_to_physics, moved the calculation of the local arrays qv_p,qc_p, and qr_p above the 
!   calculation of th_p so that we do not need to use the pointer qv.
! * in subroutine microphysics_from_MPAS, moved the calculation of the local arrays qv_p,qc_p, and qr_p above 
!   the calculation of th_p so that we do not need to use the pointer qv.
! * in subroutine microphysics_to_MPAS, since microphysics schemes update the temperature and water vapor
!   mixing ratio, we first update the total pressure and exner functions. Then, we update the modified
!   potential temperature and calculate the diabatic tendency due to cloud microphysics processes.
!   Laura D. Fowler (laura@ucar.edu) / 2013-11-07.
! * replaced the variable g (that originally pointed to gravity) with gravity, for simplicity.
!   Laura D. Fowler (laura@ucar.edu) / 2014-03-21.
! * throughout the sourcecode, replaced all "var_struct" defined arrays by local pointers.
!   Laura D. Fowler (laura@ucar.edu) / 2014-04-22.
! * modified sourcecode to use pools.
!   Laura D. Fowler (laura@ucar.edu) / 2014-05-15.
! * added calculation of the surface pressure tendency. Moved the calculation of znu_p below the calculation 
!   of the surface pressure to avoid dividing by zero if the surface pressure is not output in the init file.
!   Laura D. Fowler (birch.mmm.ucar.ecu) / 2014-06-23. 
! * renamed module mpas_atmphys_interface_nhyd to mpas_atmphys_interface.
!   Laura D. Fowler (laura@ucar.edu) / 2014-09-19.
! * in subroutine microphysics_to_MPAS, reverted the calculation of cloud microphysics tendency rt_diabatic_tend, 
!   and update of the state variables to git hash identifier 7a66f273e182f4. This change reflects the fact that
!   we want to compute rt_diabatic_tend at constant volume.
!   Laura D. Fowler (laura@ucar.edu) / 2014-01-015. 
! * added the initialization of local variables needed in the parameterization of the MYNN surface layer scheme
!   and planetary boundary layer scheme from WRF 3.6.1.
!   Laura D. Fowler (laura@ucar.edu) / 2016-04-11.
! * corrected the calculation of the surface pressure, mainly extrapolation of the air density to the surface.
!   Laura D. Fowler (laura@ucar.edu) / 2016-04-25.
! * changed microphysics update such that the microphsyics tendency is ignored above level mp_top_level
!   Bill Skamarock / 2025/10/17 

!P3 predicts its own ice categories (p3_* scalars) and does not allocate the bulk qi, qs, qg and ni scalars.
!Where code shared with the other schemes reads qi/qs/qg/ni for mp_p3, the pointers are redirected to these
!fields instead of scalars(-1,:,:): an all-zero field for qs and qg, and the sums over the P3 ice categories
!of the total ice mass (p3_qit*) and number (p3_qni*) for qi and ni. They are allocated once in
!p3_allocate_work (called from physics_init, outside OpenMP regions) and are never written back to scalars.
 real(kind=RKIND),dimension(:,:),allocatable,target:: p3_zero_qsqg,p3_qi_sum,p3_ni_sum

 contains

!=================================================================================================================
 subroutine p3_allocate_work(state)
!=================================================================================================================
!allocates the mp_p3 work fields p3_zero_qsqg, p3_qi_sum and p3_ni_sum (same shape as scalars(n,:,:)).
 type(mpas_pool_type),intent(in):: state
 real(kind=RKIND),dimension(:,:,:),pointer:: scalars

 call mpas_pool_get_array(state,'scalars',scalars,1)
 if(.not. allocated(p3_zero_qsqg)) allocate(p3_zero_qsqg(size(scalars,2),size(scalars,3)))
 if(.not. allocated(p3_qi_sum)   ) allocate(p3_qi_sum(size(scalars,2),size(scalars,3))   )
 if(.not. allocated(p3_ni_sum)   ) allocate(p3_ni_sum(size(scalars,2),size(scalars,3))   )
 p3_zero_qsqg(:,:) = 0._RKIND
 p3_qi_sum(:,:)    = 0._RKIND
 p3_ni_sum(:,:)    = 0._RKIND

 end subroutine p3_allocate_work

!=================================================================================================================
 subroutine p3_sum_ice(state,scalars,its,ite)
!=================================================================================================================
!sets p3_qi_sum and p3_ni_sum to the sums over the P3 ice categories of p3_qit and p3_qni.
 type(mpas_pool_type),intent(in):: state
 integer,intent(in):: its,ite
 real(kind=RKIND),intent(in),dimension(:,:,:):: scalars

 integer,pointer:: index_p3_qit,index_p3_qni,index_p3_qir
 integer:: c,i,k,nCat

 if(.not. allocated(p3_qi_sum) .or. .not. allocated(p3_ni_sum)) &
    call physics_error_fatal('p3_sum_ice: P3 work fields not allocated (p3_allocate_work not called)')

 call mpas_pool_get_dimension(state,'index_p3_qit',index_p3_qit)
 call mpas_pool_get_dimension(state,'index_p3_qni',index_p3_qni)
 call mpas_pool_get_dimension(state,'index_p3_qir',index_p3_qir)
 nCat = index_p3_qir - index_p3_qni

 do i = its,ite
 do k = kts,kte
    p3_qi_sum(k,i) = 0._RKIND
    p3_ni_sum(k,i) = 0._RKIND
    do c = 1, nCat
       p3_qi_sum(k,i) = p3_qi_sum(k,i) + max(0._RKIND,scalars(index_p3_qit+c-1,k,i))
       p3_ni_sum(k,i) = p3_ni_sum(k,i) + max(0._RKIND,scalars(index_p3_qni+c-1,k,i))
    enddo
 enddo
 enddo

 end subroutine p3_sum_ice


!=================================================================================================================
 subroutine allocate_forall_physics(configs)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs

!local pointers:
 character(len=StrKIND),pointer:: microp_scheme,pbl_scheme
 integer, pointer :: p3_nCat
 logical, pointer :: p3_aerosols

!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_config(configs,'config_microp_scheme',microp_scheme)
 call mpas_pool_get_config(configs,'config_pbl_scheme'   ,pbl_scheme   )

 if(.not.allocated(psfc_p) ) allocate(psfc_p(ims:ime,jms:jme)         )
 if(.not.allocated(ptop_p) ) allocate(ptop_p(ims:ime,jms:jme)         )

 if(.not.allocated(u_p)    ) allocate(u_p(ims:ime,kms:kme,jms:jme)    )
 if(.not.allocated(v_p)    ) allocate(v_p(ims:ime,kms:kme,jms:jme)    )
 if(.not.allocated(fzm_p)  ) allocate(fzm_p(ims:ime,kms:kme,jms:jme)  )
 if(.not.allocated(fzp_p)  ) allocate(fzp_p(ims:ime,kms:kme,jms:jme)  )
 if(.not.allocated(zz_p)   ) allocate(zz_p(ims:ime,kms:kme,jms:jme)   )
 if(.not.allocated(pres_p) ) allocate(pres_p(ims:ime,kms:kme,jms:jme) )
 if(.not.allocated(pi_p)   ) allocate(pi_p(ims:ime,kms:kme,jms:jme)   )
 if(.not.allocated(z_p)    ) allocate(z_p(ims:ime,kms:kme,jms:jme)    )
 if(.not.allocated(zmid_p) ) allocate(zmid_p(ims:ime,kms:kme,jms:jme) )
 if(.not.allocated(dz_p)   ) allocate(dz_p(ims:ime,kms:kme,jms:jme)   )
 if(.not.allocated(t_p)    ) allocate(t_p(ims:ime,kms:kme,jms:jme)    )
 if(.not.allocated(th_p)   ) allocate(th_p(ims:ime,kms:kme,jms:jme)   )
 if(.not.allocated(al_p)   ) allocate(al_p(ims:ime,kms:kme,jms:jme)   )
 if(.not.allocated(rho_p)  ) allocate(rho_p(ims:ime,kms:kme,jms:jme)  )
 if(.not.allocated(rh_p)   ) allocate(rh_p(ims:ime,kms:kme,jms:jme)   )
 if(.not.allocated(znu_p)  ) allocate(znu_p(ims:ime,kms:kme,jms:jme)  )

 if(.not.allocated(w_p)    ) allocate(w_p(ims:ime,kms:kme,jms:jme)    )
 if(.not.allocated(pres2_p)) allocate(pres2_p(ims:ime,kms:kme,jms:jme))
 if(.not.allocated(t2_p)   ) allocate(t2_p(ims:ime,kms:kme,jms:jme)   )
 
 if(.not.allocated(qv_p)   ) allocate(qv_p(ims:ime,kms:kme,jms:jme)   )
 if(.not.allocated(qc_p)   ) allocate(qc_p(ims:ime,kms:kme,jms:jme)   )
 if(.not.allocated(qr_p)   ) allocate(qr_p(ims:ime,kms:kme,jms:jme)   )
 if(.not.allocated(qi_p)   ) allocate(qi_p(ims:ime,kms:kme,jms:jme)   )
 if(.not.allocated(qs_p)   ) allocate(qs_p(ims:ime,kms:kme,jms:jme)   )
 if(.not.allocated(qg_p)   ) allocate(qg_p(ims:ime,kms:kme,jms:jme)   )

 microp_select: select case(trim(microp_scheme))
    case("mp_thompson_aerosols")
       if(.not.allocated(nifa_p)) allocate(nifa_p(ims:ime,kms:kme,jms:jme))
       if(.not.allocated(nwfa_p)) allocate(nwfa_p(ims:ime,kms:kme,jms:jme))

     case("mp_p3")
        call mpas_pool_get_config(configs,'config_p3_nCat',p3_nCat)
        ! P3 always predicts nc. The aerosol numbers nwfa/nifa (shared with Thompson aerosol-aware design)
        ! exist only with predicted aerosols (config_p3_aerosols, package mp_p3_aers_in):
        call mpas_pool_get_config(configs,'config_p3_aerosols',p3_aerosols)
        if(.not.allocated(nc_p))   allocate(nc_p(ims:ime,kms:kme,jms:jme))
        if(p3_aerosols) then
           if(.not.allocated(nifa_p)) allocate(nifa_p(ims:ime,kms:kme,jms:jme))
           if(.not.allocated(nwfa_p)) allocate(nwfa_p(ims:ime,kms:kme,jms:jme))
        endif
        ! P3 category prognostics live only in state%scalars (MPAS-compatible). No 4-D work arrays.
        ! qzi / qli arrays are always allocated for robustness; populate/ignore based on config flags

    case default
 end select microp_select

 pbl_select: select case(trim(pbl_scheme))
    case("bl_mynn")
       if(.not.allocated(nc_p)) allocate(nc_p(ims:ime,kms:kme,jms:jme))
       if(.not.allocated(ni_p)) allocate(ni_p(ims:ime,kms:kme,jms:jme))
       if(.not.allocated(nifa_p)) allocate(nifa_p(ims:ime,kms:kme,jms:jme))
       if(.not.allocated(nwfa_p)) allocate(nwfa_p(ims:ime,kms:kme,jms:jme))

    case default
 end select pbl_select

!... arrays used for calculating the hydrostatic pressure and exner function:
 if(.not.allocated(psfc_hyd_p)  ) allocate(psfc_hyd_p(ims:ime,jms:jme)          )
 if(.not.allocated(psfc_hydd_p) ) allocate(psfc_hydd_p(ims:ime,jms:jme)         )
 if(.not.allocated(pres_hyd_p)  ) allocate(pres_hyd_p(ims:ime,kms:kme,jms:jme)  )
 if(.not.allocated(pres_hydd_p) ) allocate(pres_hydd_p(ims:ime,kms:kme,jms:jme) )
 if(.not.allocated(pres2_hyd_p) ) allocate(pres2_hyd_p(ims:ime,kms:kme,jms:jme) )
 if(.not.allocated(pres2_hydd_p)) allocate(pres2_hydd_p(ims:ime,kms:kme,jms:jme))
 if(.not.allocated(znu_hyd_p)   ) allocate(znu_hyd_p(ims:ime,kms:kme,jms:jme)   )
 
 !P3 prognostic scalars mixed by the YSU or MYNN PBL schemes (mp_p3 only; filled in MPAS_to_physics):
 call allocate_p3pbl(configs)

 end subroutine allocate_forall_physics

!=================================================================================================================
 subroutine deallocate_forall_physics(configs)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs

!local pointers:
 character(len=StrKIND),pointer:: microp_scheme,pbl_scheme

!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_config(configs,'config_microp_scheme',microp_scheme)
 call mpas_pool_get_config(configs,'config_pbl_scheme'   ,pbl_scheme   )

 if(allocated(psfc_p)  ) deallocate(psfc_p  )
 if(allocated(ptop_p)  ) deallocate(ptop_p  )

 if(allocated(u_p)     ) deallocate(u_p     )
 if(allocated(v_p)     ) deallocate(v_p     )
 if(allocated(fzm_p)   ) deallocate(fzm_p   )
 if(allocated(fzp_p)   ) deallocate(fzp_p   )
 if(allocated(zz_p)    ) deallocate(zz_p    )
 if(allocated(pres_p)  ) deallocate(pres_p  )
 if(allocated(pi_p)    ) deallocate(pi_p    )
 if(allocated(z_p)     ) deallocate(z_p     )
 if(allocated(zmid_p)  ) deallocate(zmid_p  )
 if(allocated(dz_p)    ) deallocate(dz_p    )
 if(allocated(t_p)     ) deallocate(t_p     )
 if(allocated(th_p)    ) deallocate(th_p    )
 if(allocated(al_p)    ) deallocate(al_p    )
 if(allocated(rho_p)   ) deallocate(rho_p   ) 
 if(allocated(rh_p)    ) deallocate(rh_p    ) 
 if(allocated(znu_p)   ) deallocate(znu_p   )

 if(allocated(w_p)     ) deallocate(w_p     )
 if(allocated(pres2_p) ) deallocate(pres2_p )
 if(allocated(t2_p)    ) deallocate(t2_p    )

 if(allocated(qv_p)    ) deallocate(qv_p    )
 if(allocated(qc_p)    ) deallocate(qc_p    )
 if(allocated(qr_p)    ) deallocate(qr_p    )
 if(allocated(qi_p)    ) deallocate(qi_p    )
 if(allocated(qs_p)    ) deallocate(qs_p    )
 if(allocated(qg_p)    ) deallocate(qg_p    )

 microp_select: select case(trim(microp_scheme))
    case("mp_thompson_aerosols")
       if(allocated(nifa_p)) deallocate(nifa_p)
       if(allocated(nwfa_p)) deallocate(nwfa_p)

     case("mp_p3")
        if(allocated(nc_p))   deallocate(nc_p)
        if(allocated(nifa_p)) deallocate(nifa_p)
        if(allocated(nwfa_p)) deallocate(nwfa_p)
    case default
 end select microp_select

 pbl_select: select case(trim(pbl_scheme))
    case("bl_mynn")
       if(allocated(nc_p)) deallocate(nc_p)
       if(allocated(ni_p)) deallocate(ni_p)
       if(allocated(nifa_p)) deallocate(nifa_p)
       if(allocated(nwfa_p)) deallocate(nwfa_p)

    case default
 end select pbl_select

 if(allocated(psfc_hyd_p)  ) deallocate(psfc_hyd_p  )
 if(allocated(psfc_hydd_p) ) deallocate(psfc_hydd_p )
 if(allocated(pres_hyd_p)  ) deallocate(pres_hyd_p  )
 if(allocated(pres_hydd_p) ) deallocate(pres_hydd_p )
 if(allocated(pres2_hyd_p) ) deallocate(pres2_hyd_p )
 if(allocated(pres2_hydd_p)) deallocate(pres2_hydd_p)
 if(allocated(znu_hyd_p)   ) deallocate(znu_hyd_p   )
 
 call deallocate_p3pbl

 end subroutine deallocate_forall_physics

!=================================================================================================================
 subroutine MPAS_to_physics(configs,mesh,state,time_lev,diag,diag_physics,its,ite)
!=================================================================================================================

!input variables:
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: state
 type(mpas_pool_type),intent(in):: diag

 integer,intent(in):: its,ite
 integer,intent(in):: time_lev

!inout variables:
 type(mpas_pool_type),intent(inout):: diag_physics

!local pointers:
 character(len=StrKIND),pointer:: microp_scheme,pbl_scheme

 integer,pointer:: index_qv,index_qc,index_qr,index_qi,index_qs,index_qg
 integer,pointer:: index_nc,index_ni,index_nifa,index_nwfa

 real(kind=RKIND),dimension(:),pointer    :: latCell,lonCell
 real(kind=RKIND),dimension(:),pointer    :: fzm,fzp,rdzw
 real(kind=RKIND),dimension(:),pointer    :: surface_pressure,plrad
 real(kind=RKIND),dimension(:,:),pointer  :: zgrid
 real(kind=RKIND),dimension(:,:),pointer  :: zz,exner,pressure_b,rtheta_p,rtheta_b
 real(kind=RKIND),dimension(:,:),pointer  :: rho_zz,theta_m,pressure_p,u,v,w
 real(kind=RKIND),dimension(:,:),pointer  :: qv,qc,qr,qi,qs,qg
 real(kind=RKIND),dimension(:,:),pointer  :: nc,ni,nifa,nwfa
 real(kind=RKIND),dimension(:,:,:),pointer:: scalars

!local variables:
 integer:: i,k,j

 real(kind=RKIND):: z0,z1,z2,w1,w2
 real(kind=RKIND):: rho_a,rho1,rho2,tem1,tem2

!-----------------------------------------------------------------------------------------------------------------
!call mpas_log_write('')
!call mpas_log_write('--- enter subroutine MPAS_to_phys:')
!call mpas_log_write('ims=$i ime=$i',intArgs=(/ims,ime/))
!call mpas_log_write('jms=$i jme=$i',intArgs=(/jms,jme/))
!call mpas_log_write('kms=$i kme=$i',intArgs=(/kms,kme/))
!call mpas_log_write('')
!call mpas_log_write('its=$i ite=$i',intArgs=(/its,ite/))
!call mpas_log_write('jts=$i jte=$i',intArgs=(/jts,jte/))
!call mpas_log_write('kts=$i kte=$i',intArgs=(/kts,kte/))

!initialization:
 call mpas_pool_get_config(configs,'config_microp_scheme',microp_scheme)
 call mpas_pool_get_config(configs,'config_pbl_scheme'   ,pbl_scheme   )

 call mpas_pool_get_array(mesh,'latCell',latCell)
 call mpas_pool_get_array(mesh,'lonCell',lonCell)
 call mpas_pool_get_array(mesh,'fzm'    ,fzm    )
 call mpas_pool_get_array(mesh,'fzp'    ,fzp    )
 call mpas_pool_get_array(mesh,'rdzw'   ,rdzw   )
 call mpas_pool_get_array(mesh,'zgrid'  ,zgrid  )
 call mpas_pool_get_array(mesh,'zz'     ,zz     )

 call mpas_pool_get_array(diag,'surface_pressure'      ,surface_pressure)
 call mpas_pool_get_array(diag,'exner'                 ,exner           )
 call mpas_pool_get_array(diag,'pressure_base'         ,pressure_b      )
 call mpas_pool_get_array(diag,'pressure_p'            ,pressure_p      )
 call mpas_pool_get_array(diag,'rtheta_base'           ,rtheta_b        )
 call mpas_pool_get_array(diag,'rtheta_p'              ,rtheta_p        )
 call mpas_pool_get_array(diag,'uReconstructZonal'     ,u               )
 call mpas_pool_get_array(diag,'uReconstructMeridional',v               )

 call mpas_pool_get_array(state,'rho_zz' ,rho_zz ,time_lev)
 call mpas_pool_get_array(state,'theta_m',theta_m,time_lev)
 call mpas_pool_get_array(state,'w'      ,w      ,time_lev)

 call mpas_pool_get_dimension(state,'index_qv',index_qv)
 call mpas_pool_get_dimension(state,'index_qc',index_qc)
 call mpas_pool_get_dimension(state,'index_qr',index_qr)
 call mpas_pool_get_dimension(state,'index_qi',index_qi)
 call mpas_pool_get_dimension(state,'index_qs',index_qs)
 call mpas_pool_get_dimension(state,'index_qg',index_qg)

 call mpas_pool_get_array(state,'scalars',scalars,time_lev)
 qv => scalars(index_qv,:,:)
 qc => scalars(index_qc,:,:)
 qr => scalars(index_qr,:,:)
 qi => scalars(index_qi,:,:)
 qs => scalars(index_qs,:,:)
 qg => scalars(index_qg,:,:)
 if(trim(microp_scheme) == 'mp_p3') then
    !qi, qs and qg are not allocated with mp_p3 (or are zero, with config_p3_ic_bulk_ice): qi_p is the sum
    !over the P3 ice categories of the total ice mass, and qs_p and qg_p are zero. The P3 scalars mixed by the PBL are copied to p3mix_p.
    call p3_sum_ice(state,scalars,its,ite)
    qi => p3_qi_sum
    qs => p3_zero_qsqg
    qg => p3_zero_qsqg
    call p3pbl_from_MPAS(configs,state,scalars,its,ite)
 endif

 call mpas_pool_get_array(diag_physics,'plrad',plrad)

 do j = jts, jte
 do k = kts, kte
 do i = its, ite

    !water vapor and moist arrays:
    qv_p(i,k,j) = max(0.,qv(k,i))
    qc_p(i,k,j) = max(0.,qc(k,i))
    qr_p(i,k,j) = max(0.,qr(k,i))
    qi_p(i,k,j) = max(0.,qi(k,i))
    qs_p(i,k,j) = max(0.,qs(k,i))
    qg_p(i,k,j) = max(0.,qg(k,i))

    !arrays located at theta points:
    u_p(i,k,j) = u(k,i)
    v_p(i,k,j) = v(k,i)

    zz_p(i,k,j)  = zz(k,i)
    rho_p(i,k,j) = zz(k,i) * rho_zz(k,i)
    rho_p(i,k,j) = rho_p(i,k,j)*(1._RKIND + qv_p(i,k,j))
    th_p(i,k,j)  = theta_m(k,i) / (1._RKIND + R_v/R_d * qv_p(i,k,j))

    t_p(i,k,j)    = th_p(i,k,j) * exner(k,i)
    pi_p(i,k,j)   = exner(k,i)
    pres_p(i,k,j) = pressure_p(k,i) + pressure_b(k,i)
    zmid_p(i,k,j) = 0.5*(zgrid(k+1,i)+zgrid(k,i))
    dz_p(i,k,j)   = zgrid(k+1,i)-zgrid(k,i)

 enddo
 enddo
 enddo

 microp_select: select case(trim(microp_scheme))
    !mp_p3: nwfa/nifa are also needed by the RRTMG shortwave aerosol optics (gt_aod). They exist only
    !with predicted aerosols (config_p3_aerosols; f_nifa and f_nwfa are then true).
    case("mp_thompson_aerosols","mp_p3")
       if(f_nifa .and. f_nwfa) then
          nullify(nifa)
          nullify(nwfa)
          call mpas_pool_get_dimension(state,'index_nifa',index_nifa)
          call mpas_pool_get_dimension(state,'index_nwfa',index_nwfa)
          nifa => scalars(index_nifa,:,:)
          nwfa => scalars(index_nwfa,:,:)
          do j = jts,jte
             do k = kts,kte
                do i = its,ite
                   nifa_p(i,k,j) = max(0.,nifa(k,i))
                   nwfa_p(i,k,j) = max(0.,nwfa(k,i))
                enddo
             enddo
          enddo
       endif

    case default
 end select microp_select

 pbl_select: select case(trim(pbl_scheme))
    case("bl_mynn")
       do j = jts,jte
          do k = kts,kte
             do i = its,ite
                nc_p(i,k,j)   = 0._RKIND
                ni_p(i,k,j)   = 0._RKIND
                nifa_p(i,k,j) = 0._RKIND
                nwfa_p(i,k,j) = 0._RKIND
             enddo
          enddo
       enddo
       !initializes ni_p when running the options "mp_thompson" or "mp_thompson_aerosols":
       if(f_ni) then
          nullify(ni)
          call mpas_pool_get_dimension(state,'index_ni',index_ni)
          ni => scalars(index_ni,:,:)
          do j = jts,jte
             do k = kts,kte
                do i = its,ite
                   ni_p(i,k,j) = max(0.,ni(k,i))
                enddo
             enddo
          enddo
       endif
       !initializes nc_p, nifa_p, and nwfa_p when running the option "mp_thompson_aerosols". mp_p3 always
       !has nc, and has nifa and nwfa only with predicted aerosols (config_p3_aerosols):
       if(f_nc) then
          nullify(nc)
          call mpas_pool_get_dimension(state,'index_nc',index_nc)
          nc   => scalars(index_nc,:,:)
          do j = jts,jte
             do k = kts,kte
                do i = its,ite
                   nc_p(i,k,j)   = max(0.,nc(k,i))
                enddo
             enddo
          enddo
       endif
       if(f_nifa .and. f_nwfa) then
          nullify(nifa)
          nullify(nwfa)
          call mpas_pool_get_dimension(state,'index_nifa',index_nifa)
          call mpas_pool_get_dimension(state,'index_nwfa',index_nwfa)
          nifa => scalars(index_nifa,:,:)
          nwfa => scalars(index_nwfa,:,:)
          do j = jts,jte
             do k = kts,kte
                do i = its,ite
                   nifa_p(i,k,j) = max(0.,nifa(k,i))
                   nwfa_p(i,k,j) = max(0.,nwfa(k,i))
                enddo
             enddo
          enddo
       endif

    case default
 end select pbl_select

!calculation of the surface pressure using hydrostatic assumption down to the surface::
 do j = jts,jte
 do i = its,ite
    tem1 = zgrid(2,i)-zgrid(1,i)
    tem2 = zgrid(3,i)-zgrid(2,i)
    rho1 = rho_zz(1,i) * zz(1,i) * (1. + qv_p(i,1,j))
    rho2 = rho_zz(2,i) * zz(2,i) * (1. + qv_p(i,2,j))
!   surface_pressure(i) = 0.5*gravity*(zgrid(2,i)-zgrid(1,i)) &
!                   * (rho1 + 0.5*(rho2-rho1)*tem1/(tem1+tem2))
    surface_pressure(i) = 0.5*gravity*(zgrid(2,i)-zgrid(1,i)) &
                    * (rho1 - 0.5*(rho2-rho1)*tem1/(tem1+tem2))
    surface_pressure(i) = surface_pressure(i) + pressure_p(1,i) + pressure_b(1,i)
 enddo
 
 do k = kts,kte
 do i = its,ite
    znu_p(i,k,j)  = pres_p(i,k,j) / surface_pressure(i)
 enddo
 enddo
 enddo

!arrays located at w points:
 do j = jts, jte
 do k = kts,kte+1
 do i = its,ite
    w_p(i,k,j) = w(k,i)
    z_p(i,k,j) = zgrid(k,i)
 enddo
 enddo
 enddo

!check that the pressure in the layer above the surface is greater than that in the layer
!above it:
 do j = jts,jte
 do i = its,ite
    if(pres_p(i,1,j) .le. pres_p(i,2,j)) then
       call mpas_log_write('')
       call mpas_log_write('--- subroutine MPAS_to_phys - pressure(1) < pressure(2):')
       call mpas_log_write('i      =$i', intArgs=(/i/))
       call mpas_log_write('latCell=$r', realArgs=(/latCell(i)/degrad/))
       call mpas_log_write('lonCell=$r', realArgs=(/lonCell(i)/degrad/))
       do k = kts,kte
          call mpas_log_write('$i $i $i $r $r $r $r $r $r $r $r', intArgs=(/j,i,k/),&
             realArgs=(/dz_p(i,k,j),pressure_b(k,i),pressure_p(k,i),pres_p(i,k,j), &
                        rho_p(i,k,j),th_p(i,k,j),t_p(i,k,j),qv_p(i,k,j)/))
       enddo
!      call mpas_log_write('pressure increasing with height', messageType=MPAS_LOG_CRIT)
    endif
 enddo
 enddo

!interpolation of pressure and temperature from theta points to w points:
 do j = jts,jte
 do k = kts+1,kte
 do i = its,ite
    tem1 = 1./(zgrid(k+1,i)-zgrid(k-1,i))
    fzm_p(i,k,j) = (zgrid(k,i)-zgrid(k-1,i)) * tem1
    fzp_p(i,k,j) = (zgrid(k+1,i)-zgrid(k,i)) * tem1
    t2_p(i,k,j)    = fzm_p(i,k,j)*t_p(i,k,j) + fzp_p(i,k,j)*t_p(i,k-1,j)
    pres2_p(i,k,j) = fzm_p(i,k,j)*pres_p(i,k,j) + fzp_p(i,k,j)*pres_p(i,k-1,j)
 enddo
 enddo
 enddo

!interpolation of pressure and temperature from theta points to the top-of-the-model: follows
!the calculation of the top-of-the-model pressure and temperature in WRF (subroutine phy_prep
!in ./dyn_em/module_big_step_utilities.F):
 k = kte+1
 do j = jts,jte
 do i = its,ite
    z0 = zgrid(k,i)
    z1 = 0.5*(zgrid(k,i)+zgrid(k-1,i)) 
    z2 = 0.5*(zgrid(k-1,i)+zgrid(k-2,i))
    w1 = (z0-z2)/(z1-z2)
    w2 = 1.-w1
    t2_p(i,k,j) = w1*t_p(i,k-1,j) + w2*t_p(i,k-2,j)
    !use log of pressure to avoid occurrences of negative top-of-the-model pressure.
    pres2_p(i,k,j) = exp(w1*log(pres_p(i,k-1,j))+w2*log(pres_p(i,k-2,j)))
 enddo
 enddo

!ldf (2012-06-22): recalculates the pressure at the surface as an extrapolation of the
!pressures in the 2 layers above the surface, as was originally done:
 k = kts
 do j = jts,jte
 do i = its,ite
    z0 = zgrid(k,i)
    z1 = 0.5*(zgrid(k,i)+zgrid(k+1,i)) 
    z2 = 0.5*(zgrid(k+1,i)+zgrid(k+2,i))
    w1 = (z0-z2)/(z1-z2)
    w2 = 1.-w1
    t2_p(i,k,j)    = w1*t_p(i,k,j)+w2*t_p(i,k+1,j)
    pres2_p(i,k,j) = w1*pres_p(i,k,j)+w2*pres_p(i,k+1,j)
!   psfc_p(i,j) = pres2_p(i,k,j)
    psfc_p(i,j) = surface_pressure(i)
 enddo
 enddo

!calculation of the hydrostatic pressure:
 do j = jts,jte
 do i = its,ite
    !pressure at w-points:
    k = kte+1
    pres2_hyd_p(i,k,j)  = pres2_p(i,k,j)
    pres2_hydd_p(i,k,j) = pres2_p(i,k,j)
    do k = kte,1,-1
       rho_a = rho_p(i,k,j) / (1.+qv_p(i,k,j))
       pres2_hyd_p(i,k,j)  = pres2_hyd_p(i,k+1,j)  + gravity*rho_p(i,k,j)*dz_p(i,k,j)
       pres2_hydd_p(i,k,j) = pres2_hydd_p(i,k+1,j) + gravity*rho_a*dz_p(i,k,j)
    enddo
    !pressure at theta-points:
    do k = kte,1,-1
       pres_hyd_p(i,k,j)  = 0.5*(pres2_hyd_p(i,k+1,j)+pres2_hyd_p(i,k,j))
       pres_hydd_p(i,k,j) = 0.5*(pres2_hydd_p(i,k+1,j)+pres2_hydd_p(i,k,j))
    enddo
    !surface pressure:
    psfc_hyd_p(i,j) = pres2_hyd_p(i,1,j)
    psfc_hydd_p(i,j) = pres2_hydd_p(i,1,j)
    !znu:
    do k = kte,1,-1
       znu_hyd_p(i,k,j) = pres_hyd_p(i,k,j) / psfc_hyd_p(i,j) 
    enddo
 enddo
 enddo

!save the model-top pressure:
 do j = jts,jte
 do i = its,ite
    plrad(i) = pres2_p(i,kte+1,j) 
 enddo
 enddo

 end subroutine MPAS_to_physics

!=================================================================================================================
 subroutine pre_microphysics(configs,state,diag,time_lev)
!=================================================================================================================

!input variables:
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: state
 type(mpas_pool_type),intent(in):: diag

 integer:: time_lev

#ifdef MPAS_OPENACC
!local pointers:
 real(kind=RKIND),dimension(:,:),pointer  :: exner,pressure_b,w
 real(kind=RKIND),dimension(:,:),pointer  :: rho_zz,theta_m,pressure_p
 real(kind=RKIND),dimension(:,:,:),pointer:: scalars


 MPAS_ACC_TIMER_START('pre_microphysics [ACC_data_xfer]')
 call mpas_pool_get_array(diag,'exner'        ,exner     )
 call mpas_pool_get_array(diag,'pressure_base',pressure_b)
 call mpas_pool_get_array(diag,'pressure_p'   ,pressure_p)

 call mpas_pool_get_array(state,'rho_zz' ,rho_zz ,time_lev)
 call mpas_pool_get_array(state,'theta_m',theta_m,time_lev)
 call mpas_pool_get_array(state,'w'      ,w      ,time_lev)
 !$acc update host(exner, pressure_b, pressure_p, rho_zz, theta_m, w)

 call mpas_pool_get_array(state,'scalars',scalars,time_lev)
 !$acc update host(scalars)

 MPAS_ACC_TIMER_STOP('pre_microphysics [ACC_data_xfer]')
#endif

end subroutine pre_microphysics

!=================================================================================================================
 subroutine microphysics_from_MPAS(configs,mesh,state,time_lev,diag,diag_physics,tend_physics,its,ite)
!=================================================================================================================

!input variables:
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: state
 type(mpas_pool_type),intent(in):: diag
 type(mpas_pool_type),intent(in):: diag_physics

 integer,intent(in):: its,ite
 integer:: time_lev

!inout variables:
 type(mpas_pool_type),intent(inout):: tend_physics

!local pointers:
 character(len=StrKIND),pointer:: mp_scheme
 integer,pointer:: index_qv,index_qc,index_qr,index_qi,index_qs,index_qg
 integer,pointer:: index_nc,index_ni,index_nr,index_nifa,index_nwfa
 real(kind=RKIND),dimension(:),pointer    :: nifa2d,nwfa2d,nt_c,mu_c
 real(kind=RKIND),dimension(:,:),pointer  :: zgrid,w
 real(kind=RKIND),dimension(:,:),pointer  :: zz,exner,pressure_b
 real(kind=RKIND),dimension(:,:),pointer  :: rho_zz,theta_m,pressure_p
 real(kind=RKIND),dimension(:,:),pointer  :: qv,qc,qr,qi,qs,qg
 real(kind=RKIND),dimension(:,:),pointer  :: nc,ni,nr,nifa,nwfa
 real(kind=RKIND),dimension(:,:),pointer  :: rainprod,evapprod
 real(kind=RKIND),dimension(:,:),pointer  :: re_cloud,re_ice,re_snow
 real(kind=RKIND),dimension(:,:),pointer  :: rthmpten,rqvmpten,rqcmpten,rqrmpten,rqimpten,rqsmpten,rqgmpten
 real(kind=RKIND),dimension(:,:),pointer  :: rncmpten,rnimpten,rnrmpten,rnifampten,rnwfampten
 real(kind=RKIND),dimension(:,:,:),pointer:: scalars

!local variables:
 integer:: i,k,j,c
 integer :: nScalarTot
 integer :: nCat_layout

 logical, save :: p3dbg_banner_printed = .false.
 logical, save :: p3dbg_to_indices_checked = .false.
 character(len=1024) :: dbgmsg


 real(kind=RKIND) :: p3_kappa, p3_inv1mk
 real(kind=RKIND) :: moistfac, qvpos
 real(kind=RKIND) :: rho_d_tmp
 ! --- P3 config ---
 integer, pointer :: p3_nCat
 logical, pointer :: p3_trplMomI, p3_liqfrac
 logical, pointer :: p3_debug_on
 ! --- P3 tracer base indices in scalars(:,:,:) ---
 integer, pointer :: index_p3_qit, index_p3_qni, index_p3_qir, index_p3_qib
 integer, pointer :: index_p3_qzi, index_p3_qli
 !-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_config(configs,'config_microp_scheme',mp_scheme)
 call mpas_pool_get_config(configs,'config_p3_debug_on',p3_debug_on)



   


 p3_kappa  = R_d/cp
 p3_inv1mk = 1._RKIND / (1._RKIND - p3_kappa)
 call mpas_pool_get_array(mesh,'zgrid',zgrid)
 call mpas_pool_get_array(mesh,'zz'   ,zz   )

 call mpas_pool_get_array(diag,'exner'        ,exner     )
 call mpas_pool_get_array(diag,'pressure_base',pressure_b)
 call mpas_pool_get_array(diag,'pressure_p'   ,pressure_p)

 call mpas_pool_get_array(state,'rho_zz' ,rho_zz ,time_lev)
 call mpas_pool_get_array(state,'theta_m',theta_m,time_lev)
 call mpas_pool_get_array(state,'w'      ,w      ,time_lev)

 call mpas_pool_get_dimension(state,'index_qv',index_qv)
 call mpas_pool_get_dimension(state,'index_qc',index_qc)
 call mpas_pool_get_dimension(state,'index_qr',index_qr)
 call mpas_pool_get_array(state,'scalars',scalars,time_lev)
 qv => scalars(index_qv,:,:)
 qc => scalars(index_qc,:,:)
 qr => scalars(index_qr,:,:)

!initialize variables needed in the cloud microphysics schemes:
 do j = jts, jte
 do k = kts, kte
 do i = its, ite
    qv_p(i,k,j) = qv(k,i)
    qc_p(i,k,j) = qc(k,i)
    qr_p(i,k,j) = qr(k,i)

    rho_p(i,k,j)  = zz(k,i) * rho_zz(k,i)
    th_p(i,k,j)   = theta_m(k,i) / (1._RKIND + R_v/R_d * max(0._RKIND,qv_p(i,k,j)))

    pi_p(i,k,j)   = exner(k,i)
    pres_p(i,k,j) = pressure_b(k,i) + pressure_p(k,i)

    z_p(i,k,j)    = zgrid(k,i)
    dz_p(i,k,j)   = zgrid(k+1,i) - zgrid(k,i)
    w_p(i,k,j)    = w(k,i)
 enddo
 enddo
 enddo

!initialize cloud water species and aerosols as function of cloud microphysics scheme:
 mp_select: select case(trim(mp_scheme))
    case("mp_thompson","mp_thompson_aerosols","mp_wsm6")
       call mpas_pool_get_dimension(state,'index_qi',index_qi)
       call mpas_pool_get_dimension(state,'index_qs',index_qs)
       call mpas_pool_get_dimension(state,'index_qg',index_qg)
       qi => scalars(index_qi,:,:)
       qs => scalars(index_qs,:,:)
       qg => scalars(index_qg,:,:)

       call mpas_pool_get_array(diag_physics,'rainprod',rainprod)
       call mpas_pool_get_array(diag_physics,'evapprod',evapprod)
       call mpas_pool_get_array(diag_physics,'re_cloud',re_cloud)
       call mpas_pool_get_array(diag_physics,'re_ice'  ,re_ice  )
       call mpas_pool_get_array(diag_physics,'re_snow' ,re_snow )

       do j = jts, jte
       do k = kts, kte
       do i = its, ite
          qi_p(i,k,j) = qi(k,i)
          qs_p(i,k,j) = qs(k,i)
          qg_p(i,k,j) = qg(k,i)

          rainprod_p(i,k,j) = rainprod(k,i)
          evapprod_p(i,k,j) = evapprod(k,i)
          recloud_p(i,k,j)  = re_cloud(k,i)
          reice_p(i,k,j)    = re_ice(k,i)
          resnow_p(i,k,j)   = re_snow(k,i)
       enddo
       enddo
       enddo
       mp2_select: select case(trim(mp_scheme))
          case("mp_thompson","mp_thompson_aerosols")
             call mpas_pool_get_dimension(state,'index_ni',index_ni)
             call mpas_pool_get_dimension(state,'index_nr',index_nr)
             ni   => scalars(index_ni,:,:)
             nr   => scalars(index_nr,:,:)

             call mpas_pool_get_array(diag_physics,'nt_c',nt_c)
             call mpas_pool_get_array(diag_physics,'mu_c',mu_c)
             do j = jts,jte
             do i = its,ite
                muc_p(i,j) = mu_c(i)
                ntc_p(i,j) = nt_c(i)
             enddo
             do k = kts, kte
             do i = its, ite
                ni_p(i,k,j) = ni(k,i)
                nr_p(i,k,j) = nr(k,i)
             enddo
             enddo
             enddo

          mp3_select: select case(trim(mp_scheme))
             case("mp_thompson_aerosols")
                call mpas_pool_get_dimension(state,'index_nc'  ,index_nc  )
                call mpas_pool_get_dimension(state,'index_nifa',index_nifa)
                call mpas_pool_get_dimension(state,'index_nwfa',index_nwfa)
                nc   => scalars(index_nc,:,:)
                nifa => scalars(index_nifa,:,:)
                nwfa => scalars(index_nwfa,:,:)

                call mpas_pool_get_array(diag_physics,'nifa2d',nifa2d)
                call mpas_pool_get_array(diag_physics,'nwfa2d',nwfa2d)
                do j = jts,jte
                do i = its,ite
                   nifa2d_p(i,j) = nifa2d(i)
                   nwfa2d_p(i,j) = nwfa2d(i)
                enddo
                do k = kts, kte
                do i = its, ite
                   nc_p(i,k,j) = nc(k,i)
                   nifa_p(i,k,j) = nifa(k,i)
                   nwfa_p(i,k,j) = nwfa(k,i)
                enddo
                enddo
                enddo

             case default
          end select mp3_select

          case default
       end select mp2_select
    case("mp_p3")
       call mpas_pool_get_config(configs,'config_p3_nCat'    ,p3_nCat    )
       call mpas_pool_get_config(configs,'config_p3_trplMomI',p3_trplMomI)
       call mpas_pool_get_config(configs,'config_p3_liqfrac' ,p3_liqfrac )

        if ((associated(p3_debug_on) .and. p3_debug_on) .and. .not. p3dbg_banner_printed) then
           p3dbg_banner_printed = .true.
        end if! MPAS bulk-compatible packing of bulk hydrometeors + number concentrations
! so downstream MPAS diagnostics and radiation coupling see the same fields as Thompson.
! qi, qs, qg and ni are not allocated with mp_p3 (qi, qs, qg are zero with config_p3_ic_bulk_ice): qi and ni are the sums over the P3 ice categories
! of p3_qit and p3_qni, and qs and qg are zero (also used by the tendency diagnostics below).
call p3_sum_ice(state,scalars,its,ite)
qi => p3_qi_sum
qs => p3_zero_qsqg
qg => p3_zero_qsqg
ni => p3_ni_sum

call mpas_pool_get_dimension(state,'index_nr',index_nr)
nr => scalars(index_nr,:,:)

do j = jts, jte
do k = kts, kte
do i = its, ite
   qi_p(i,k,j) = max(0._RKIND, qi(k,i))
   qs_p(i,k,j) = max(0._RKIND, qs(k,i))
   qg_p(i,k,j) = max(0._RKIND, qg(k,i))
   ni_p(i,k,j) = max(0._RKIND, ni(k,i))
   nr_p(i,k,j) = max(0._RKIND, nr(k,i))
end do
end do
end do
       if (.not. associated(scalars)) then
          stop 'scalars not associated in microphysics_from_MPAS P3 branch'
       end if
       if (.not. allocated(nc_p)) then
          stop 'nc_p not allocated in microphysics_from_MPAS'
       end if
       ! nwfa/nifa (and nwfa2d/nifa2d) exist only with predicted aerosols (config_p3_aerosols, package
       ! mp_p3_aers_in); f_nifa and f_nwfa are then true:
       if (f_nifa .and. f_nwfa .and. (.not. allocated(nifa_p) .or. .not. allocated(nwfa_p))) then
          stop 'P3 aerosol number arrays not allocated in microphysics_from_MPAS'
       end if
! Aerosol-aware inputs (reuse MPAS aerosol-aware aerosol-aware tracers/diagnostics)
call mpas_pool_get_dimension(state,'index_nc'  ,index_nc  )
nc   => scalars(index_nc,:,:)
if (f_nifa .and. f_nwfa) then
   call mpas_pool_get_dimension(state,'index_nifa',index_nifa)
   call mpas_pool_get_dimension(state,'index_nwfa',index_nwfa)
   nifa => scalars(index_nifa,:,:)
   nwfa => scalars(index_nwfa,:,:)
end if

call mpas_pool_get_array(diag_physics,'nifa2d'  ,nifa2d  )
call mpas_pool_get_array(diag_physics,'nwfa2d'  ,nwfa2d  )
call mpas_pool_get_array(diag_physics,'rainprod',rainprod)
call mpas_pool_get_array(diag_physics,'evapprod',evapprod)
call mpas_pool_get_array(diag_physics,'re_cloud',re_cloud)
call mpas_pool_get_array(diag_physics,'re_ice'  ,re_ice  )
call mpas_pool_get_array(diag_physics,'re_snow' ,re_snow )
       if (.not. associated(rainprod) .or. &
           .not. associated(evapprod) .or. .not. associated(re_cloud) .or. .not. associated(re_ice) .or. &
           .not. associated(re_snow)) then
          stop 'P3: diag_physics aerosol arrays not associated'
       end if
       if (f_nifa .and. f_nwfa .and. (.not. associated(nifa2d) .or. .not. associated(nwfa2d))) then
          stop 'P3: diag_physics aerosol arrays not associated'
       end if
       ! P3 aerosol / optics diagnostic work arrays must be allocated during
       ! allocate_microphysics. Do not allocate them here: on restart runs this
       ! routine can be reached immediately after restart-time output, when PIO
       ! and restart buffers are still near the memory high-water mark.
       if ((f_nifa .and. f_nwfa .and. (.not. allocated(nifa2d_p) .or. .not. allocated(nwfa2d_p))) .or. &
           .not. allocated(rainprod_p) .or. &
           .not. allocated(evapprod_p) .or. .not. allocated(recloud_p) .or. .not. allocated(reice_p) .or. &
           .not. allocated(resnow_p)) then
          call mpas_log_write('P3_MEMFIX_FATAL: P3 work arrays were not allocated in allocate_microphysics')
          stop 'P3 work arrays must be allocated during allocate_microphysics'
       end if
do j = jts,jte
if (f_nifa .and. f_nwfa) then
   do i = its,ite
      nifa2d_p(i,j) = nifa2d(i)
      nwfa2d_p(i,j) = nwfa2d(i)
   enddo
   do k = kts,kte
   do i = its,ite
      nifa_p(i,k,j)    = nifa(k,i)
      nwfa_p(i,k,j)    = nwfa(k,i)
   enddo
   enddo
end if
do k = kts,kte
do i = its,ite
   nc_p(i,k,j)      = nc(k,i)
   ! P3 overwrites these with post-limiter internal process rates in the driver.
   rainprod_p(i,k,j)= 0._RKIND
   evapprod_p(i,k,j)= 0._RKIND
   recloud_p(i,k,j) = re_cloud(k,i)
   reice_p(i,k,j)   = re_ice(k,i)
   resnow_p(i,k,j)  = re_snow(k,i)
enddo
enddo
enddo

        ! One-time sanity check of P3 scalar layout (avoid per-timestep spam)
        if (.not. p3dbg_to_indices_checked) then
           call mpas_pool_get_dimension(state,'index_p3_qit',index_p3_qit)
           call mpas_pool_get_dimension(state,'index_p3_qni',index_p3_qni)
           call mpas_pool_get_dimension(state,'index_p3_qir',index_p3_qir)
           call mpas_pool_get_dimension(state,'index_p3_qib',index_p3_qib)

           nScalarTot  = size(scalars,1)
           nCat_layout = index_p3_qir - index_p3_qni   ! p3_qit* are in the moist group

           if (nCat_layout < 1 .or. nCat_layout > 3) then
              write(dbgmsg,'(A,I0)') 'P3DBG_FATAL: invalid derived nCat_layout=', nCat_layout
              call mpas_log_write(trim(dbgmsg))
              stop 'P3 scalar layout/config mismatch'
           end if

           if (nCat_layout /= p3_nCat) then
              write(dbgmsg,'(A,I0,A,I0)') 'P3DBG_FATAL: nCat_layout=', nCat_layout, ' != config_p3_nCat=', p3_nCat
              call mpas_log_write(trim(dbgmsg))
              stop 'P3 scalar layout/config mismatch'
           end if

           if ( index_p3_qit + p3_nCat - 1 > nScalarTot .or. &
                index_p3_qni + p3_nCat - 1 > nScalarTot .or. &
                index_p3_qir + p3_nCat - 1 > nScalarTot .or. &
                index_p3_qib + p3_nCat - 1 > nScalarTot ) then
              call mpas_log_write('P3DBG_FATAL: P3 scalar index out of range in microphysics_from_MPAS')
              stop 'P3 scalar index out of range in microphysics_from_MPAS'
           end if

           if (p3_trplMomI) then
              call mpas_pool_get_dimension(state,'index_p3_qzi',index_p3_qzi)
              if ( index_p3_qzi + p3_nCat - 1 > nScalarTot ) then
                 call mpas_log_write('P3DBG_FATAL: P3 qzi index out of range in microphysics_from_MPAS')
                 stop 'P3 qzi index out of range in microphysics_from_MPAS'
              end if
           end if
           if (p3_liqfrac) then
              call mpas_pool_get_dimension(state,'index_p3_qli',index_p3_qli)
              if ( index_p3_qli + p3_nCat - 1 > nScalarTot ) then
                 call mpas_log_write('P3DBG_FATAL: P3 qli index out of range in microphysics_from_MPAS')
                 stop 'P3 qli index out of range in microphysics_from_MPAS'
              end if
           end if

           if (associated(p3_debug_on) .and. p3_debug_on) then
              write(dbgmsg,'(A,I0,A,I0)') 'P3DBG_TO: nCat_layout=', nCat_layout, ' config_p3_nCat=', p3_nCat
              call mpas_log_write(trim(dbgmsg))
           end if

           p3dbg_to_indices_checked = .true.
        end if
       ! P3 categories are consumed directly from state%scalars in the driver (column-local packing).
    case default
 end select mp_select

!begin calculation of cloud microphysics tendencies:
! For P3 these arrays hold the pre-microphysics state and are converted to
! realized, full-step net tendencies in microphysics_to_MPAS.
 mp_tend_select: select case(trim(mp_scheme))
    case("mp_thompson","mp_thompson_aerosols","mp_p3","mp_wsm6")
       call mpas_pool_get_array(tend_physics,'rthmpten',rthmpten)
       call mpas_pool_get_array(tend_physics,'rqvmpten',rqvmpten)
       call mpas_pool_get_array(tend_physics,'rqcmpten',rqcmpten)
       call mpas_pool_get_array(tend_physics,'rqrmpten',rqrmpten)
       call mpas_pool_get_array(tend_physics,'rqimpten',rqimpten)
       call mpas_pool_get_array(tend_physics,'rqsmpten',rqsmpten)
       call mpas_pool_get_array(tend_physics,'rqgmpten',rqgmpten)

       do k = kts,kte
       do i = its,ite
          rthmpten(k,i) = theta_m(k,i)/(1._RKIND+R_v/R_d*max(0._RKIND,qv(k,i)))
          rqvmpten(k,i) = qv(k,i)
          rqcmpten(k,i) = qc(k,i)
          rqrmpten(k,i) = qr(k,i)
          rqimpten(k,i) = qi(k,i)
          rqsmpten(k,i) = qs(k,i)
          rqgmpten(k,i) = qg(k,i)
       enddo
       enddo

       mp2_tend_select: select case(trim(mp_scheme))
          case("mp_thompson","mp_thompson_aerosols","mp_p3")
             call mpas_pool_get_array(tend_physics,'rnimpten',rnimpten)
             call mpas_pool_get_array(tend_physics,'rnrmpten',rnrmpten)

             do k = kts,kte
             do i = its,ite
                rnimpten(k,i) = ni(k,i)
                rnrmpten(k,i) = nr(k,i)
             enddo
             enddo

             mp3_tend_select: select case(trim(mp_scheme))
                case("mp_thompson_aerosols","mp_p3")
                   !mp_p3: nifa and nwfa exist only with predicted aerosols (config_p3_aerosols):
                   call mpas_pool_get_array(tend_physics,'rncmpten',rncmpten)

                   do k = kts,kte
                   do i = its,ite
                      rncmpten(k,i) = nc(k,i)
                   enddo
                   enddo
                   if(f_nifa .and. f_nwfa) then
                   call mpas_pool_get_array(tend_physics,'rnifampten',rnifampten)
                   call mpas_pool_get_array(tend_physics,'rnwfampten',rnwfampten)
                   do k = kts,kte
                   do i = its,ite
                      rnifampten(k,i) = nifa(k,i)
                      rnwfampten(k,i) = nwfa(k,i)
                   enddo
                   enddo
                   endif

                case default
             end select mp3_tend_select

          case default
       end select mp2_tend_select

    case default
 end select mp_tend_select

 end subroutine microphysics_from_MPAS

!=================================================================================================================
 subroutine microphysics_to_MPAS(configs,mesh,state,time_lev,diag,diag_physics,tend_physics,tend,its,ite)
!=================================================================================================================

!input variables:
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: mesh

 integer,intent(in):: time_lev
 integer,intent(in):: its,ite

!inout variables:
 type(mpas_pool_type),intent(inout):: state
 type(mpas_pool_type),intent(inout):: diag
 type(mpas_pool_type),intent(inout):: tend
 type(mpas_pool_type),intent(inout):: diag_physics
 type(mpas_pool_type),intent(inout):: tend_physics


!local pointers:
 character(len=StrKIND),pointer:: mp_scheme
 integer,pointer:: index_qv,index_qc,index_qr,index_qi,index_qs,index_qg
 integer,pointer:: index_nc,index_ni,index_nr,index_nifa,index_nwfa
 real(kind=RKIND),dimension(:),pointer    :: surface_pressure,tend_sfc_pressure
 real(kind=RKIND),dimension(:),pointer    :: nifa2d,nwfa2d
 real(kind=RKIND),dimension(:,:),pointer  :: zgrid
 real(kind=RKIND),dimension(:,:),pointer  :: zz,exner,exner_b,pressure_b,rtheta_p,rtheta_b
 real(kind=RKIND),dimension(:,:),pointer  :: rho_zz,theta_m,pressure_p
 real(kind=RKIND),dimension(:,:),pointer  :: rt_diabatic_tend
 real(kind=RKIND),dimension(:,:),pointer  :: dtheta_dt_mp
 real(kind=RKIND),dimension(:,:),pointer  :: qv,qc,qr,qi,qs,qg
 real(kind=RKIND),dimension(:,:),pointer  :: nc,ni,nr,nifa,nwfa
 real(kind=RKIND),dimension(:,:),pointer  :: rainprod,evapprod,refl10cm
 real(kind=RKIND),dimension(:,:),pointer  :: re_cloud,re_ice,re_snow
 real(kind=RKIND),dimension(:,:),pointer  :: rthmpten,rqvmpten,rqcmpten,rqrmpten,rqimpten,rqsmpten,rqgmpten
 real(kind=RKIND),dimension(:,:),pointer  :: rncmpten,rnimpten,rnrmpten,rnifampten,rnwfampten
 real(kind=RKIND),dimension(:,:,:),pointer:: scalars
 integer,pointer:: mp_top_level

!local variables:

 ! P3 log throttles (avoid multi-GB logs)
  logical, save :: p3dbg_to_indices_checked = .false.

! P3 debug / safety logging locals
 character(len=1024) :: dbgmsg
 integer :: nScalarTot_loc, nScalar1_loc, nScalar2_loc, nScalar3_loc
 logical :: p3_index_ok
 integer:: icount
 integer:: i,k,j,c
 real(kind=RKIND):: rho1,rho2,tem1,tem2

 ! --- P3 config ---
 integer, pointer :: p3_nCat
 logical, pointer :: p3_trplMomI, p3_liqfrac
 logical, pointer :: p3_debug_on
 ! --- P3 tracer base indices in scalars(:,:,:) ---
 integer, pointer :: index_p3_qit, index_p3_qni, index_p3_qir, index_p3_qib
 integer, pointer :: index_p3_qzi, index_p3_qli

!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_config(configs,'config_microp_scheme',mp_scheme)
 call mpas_pool_get_config(configs,'config_p3_debug_on',p3_debug_on)


 call mpas_pool_get_array(mesh,'zz'   ,zz   )
 call mpas_pool_get_array(mesh,'zgrid',zgrid)

 call mpas_pool_get_array(tend_physics,'mp_top_level',mp_top_level)

 call mpas_pool_get_array(diag,'exner'           ,exner           )
 call mpas_pool_get_array(diag,'exner_base'      ,exner_b         )
 call mpas_pool_get_array(diag,'pressure_base'   ,pressure_b      )
 call mpas_pool_get_array(diag,'pressure_p'      ,pressure_p      )
 call mpas_pool_get_array(diag,'rtheta_base'     ,rtheta_b        )
 call mpas_pool_get_array(diag,'rtheta_p'        ,rtheta_p        )
 call mpas_pool_get_array(diag,'surface_pressure',surface_pressure)
 call mpas_pool_get_array(diag,'dtheta_dt_mp'    ,dtheta_dt_mp    )

 call mpas_pool_get_array(tend,'tend_sfc_pressure',tend_sfc_pressure)

 call mpas_pool_get_array(state,'rho_zz' ,rho_zz ,time_lev)
 call mpas_pool_get_array(state,'theta_m',theta_m,time_lev)

 call mpas_pool_get_dimension(state,'index_qv'  ,index_qv  )
 call mpas_pool_get_dimension(state,'index_qc'  ,index_qc  )
 call mpas_pool_get_dimension(state,'index_qr'  ,index_qr  )
 call mpas_pool_get_array(state,'scalars',scalars,time_lev)
 qv   => scalars(index_qv,:,:)
 qc   => scalars(index_qc,:,:)
 qr   => scalars(index_qr,:,:)

! P3DBG_TO: scalar dimension and basic index mapping (compute always; print once when debug enabled)
 nScalar1_loc = size(scalars,1)
 nScalar2_loc = size(scalars,2)
 nScalar3_loc = size(scalars,3)
 if (trim(mp_scheme) == 'mp_p3' .and. associated(p3_debug_on) .and. p3_debug_on .and. &
     (.not. p3dbg_to_indices_checked)) then
    write(dbgmsg,'(A,3I8)') 'P3DBG_TO: scalars dims(nScalar,nVert,nCells)=', &
                            nScalar1_loc, nScalar2_loc, nScalar3_loc
    call mpas_log_write(trim(dbgmsg))
    write(dbgmsg,'(A,3I8)') 'P3DBG_TO: qv,qc,qr indices=', &
                            index_qv, index_qc, index_qr
    call mpas_log_write(trim(dbgmsg))
 end if

 call mpas_pool_get_array(tend,'rt_diabatic_tend',rt_diabatic_tend)

!update variables needed in the dynamical core:
 do j = jts,jte
 do k = kts,mp_top_level
 do i = its,ite

    !initializes tendency of coupled potential temperature potential temperature, and
    !potential temperature heating rate from microphysics:
    rt_diabatic_tend(k,i) = theta_m(k,i)
    dtheta_dt_mp(k,i)     = theta_m(k,i)/(1._RKIND+rvord*qv(k,i))

    !updates water vapor, cloud liquid water, rain mixing ratios, modified potential temperature,
    !tendency of coupled potential temperature, and potential temperature heating rate from microphysics:
    qv(k,i) = qv_p(i,k,j)
    qc(k,i) = qc_p(i,k,j)
    qr(k,i) = qr_p(i,k,j)

    theta_m(k,i) = th_p(i,k,j) * (1._RKIND+rvord*qv_p(i,k,j))
    rt_diabatic_tend(k,i) = (theta_m(k,i) - rt_diabatic_tend(k,i))/dt_dyn
    dtheta_dt_mp(k,i)     = (theta_m(k,i)/(1._RKIND+rvord*qv(k,i))-dtheta_dt_mp(k,i))/(dt_dyn)

    !density-weighted perturbation potential temperature:
    rtheta_p(k,i) = rho_zz(k,i) * theta_m(k,i) - rtheta_b(k,i)

    !exner function:
    exner(k,i) = (zz(k,i)*(R_d/P0)*(rtheta_p(k,i)+rtheta_b(k,i)))**rcv

    !pertubation pressure:
    pressure_p(k,i) = zz(k,i)*R_d*(exner(k,i)*rtheta_p(k,i) &
                    + (exner(k,i)-exner_b(k,i))*rtheta_b(k,i))

 enddo
 enddo
 enddo

 ! no microphysics updates above mp_top_level
 ! why are we not applying this to theta and qv?
 if(kte .gt. mp_top_level) then
    do i = its,ite
    do k = mp_top_level+1,kte
       qc(k,i) = 0.
       qr(k,i) = 0.
       rt_diabatic_tend(k,i) = 0.
    enddo
    enddo
 endif
 
!update surface pressure and calculates the surface pressure tendency:
 do j = jts,jte
 do i = its,ite
    tem1 = zgrid(2,i)-zgrid(1,i)
    tem2 = zgrid(3,i)-zgrid(2,i)
    rho1 = rho_zz(1,i) * zz(1,i) * (1. + qv_p(i,1,j))
    rho2 = rho_zz(2,i) * zz(2,i) * (1. + qv_p(i,2,j))

    tend_sfc_pressure(i) = surface_pressure(i)
!   surface_pressure(i)  = 0.5*gravity*(zgrid(2,i)-zgrid(1,i)) &
!                        * (rho1 + 0.5*(rho2-rho1)*tem1/(tem1+tem2))
    surface_pressure(i)  = 0.5*gravity*(zgrid(2,i)-zgrid(1,i)) &
                         * (rho1 - 0.5*(rho2-rho1)*tem1/(tem1+tem2))
    surface_pressure(i)  = surface_pressure(i) + pressure_p(1,i) + pressure_b(1,i)
    tend_sfc_pressure(i)   = (surface_pressure(i)-tend_sfc_pressure(i)) / dt_dyn
 enddo
 enddo

!update cloud water species and aerosols as functions of cloud microphysics schemes:
 mp_select: select case(trim(mp_scheme))
    case("mp_thompson","mp_thompson_aerosols","mp_wsm6")
       call mpas_pool_get_dimension(state,'index_qi',index_qi)
       call mpas_pool_get_dimension(state,'index_qs',index_qs)
       call mpas_pool_get_dimension(state,'index_qg',index_qg)
       qi => scalars(index_qi,:,:)
       qs => scalars(index_qs,:,:)
       qg => scalars(index_qg,:,:)

       call mpas_pool_get_array(diag_physics,'rainprod',rainprod)
       call mpas_pool_get_array(diag_physics,'evapprod',evapprod)
       call mpas_pool_get_array(diag_physics,'re_cloud',re_cloud)
       call mpas_pool_get_array(diag_physics,'re_ice'  ,re_ice  )
       call mpas_pool_get_array(diag_physics,'re_snow' ,re_snow )
       call mpas_pool_get_array(diag_physics,'refl10cm',refl10cm)

       do j = jts,jte
       do k = kts,mp_top_level
       do i = its,ite
          qi(k,i) = qi_p(i,k,j)
          qs(k,i) = qs_p(i,k,j)
          qg(k,i) = qg_p(i,k,j)

          rainprod(k,i) = rainprod_p(i,k,j)
          evapprod(k,i) = evapprod_p(i,k,j)
          re_cloud(k,i) = recloud_p(i,k,j)
          re_ice(k,i)   = reice_p(i,k,j)
          re_snow(k,i)  = resnow_p(i,k,j)
       enddo
       enddo
       enddo

       ! no microphysics updates above mp_top_level
       if(kte .gt. mp_top_level) then
          do i = its,ite
          do k = mp_top_level+1,kte
             qi(k,i) = 0.
             qs(k,i) = 0.
             qg(k,i) = 0.
             rainprod(k,i) = 0.
             evapprod(k,i) = 0.
             re_cloud(k,i) = re_qc_bg
             re_ice(k,i)   = re_qi_bg
             re_snow(k,i)  = re_qs_bg
          enddo
          enddo
       end if

       mp2_select: select case(trim(mp_scheme))
          case("mp_thompson","mp_thompson_aerosols")
             call mpas_pool_get_dimension(state,'index_ni',index_ni)
             call mpas_pool_get_dimension(state,'index_nr',index_nr)
             ni   => scalars(index_ni,:,:)
             nr   => scalars(index_nr,:,:)

             do j = jts,jte
             do k = kts,kte
             do i = its,ite
                ni(k,i) = ni_p(i,k,j)
                nr(k,i) = nr_p(i,k,j)
             enddo
             enddo
             enddo

             mp3_select: select case(trim(mp_scheme))
                case("mp_thompson_aerosols")
                   call mpas_pool_get_dimension(state,'index_nc'  ,index_nc  )
                   call mpas_pool_get_dimension(state,'index_nifa',index_nifa)
                   call mpas_pool_get_dimension(state,'index_nwfa',index_nwfa)
                   nc   => scalars(index_nc,:,:)
                   nifa => scalars(index_nifa,:,:)
                   nwfa => scalars(index_nwfa,:,:)

                   call mpas_pool_get_array(diag_physics,'nifa2d',nifa2d)
                   call mpas_pool_get_array(diag_physics,'nwfa2d',nwfa2d)
                   do j = jts,jte
                   do i = its,ite
                      nifa2d(i) = nifa2d_p(i,j)
                      nwfa2d(i) = nwfa2d_p(i,j)
                   enddo
                   do k = kts, kte
                   do i = its, ite
                      nc(k,i)   = nc_p(i,k,j)
                      nifa(k,i) = nifa_p(i,k,j)
                      nwfa(k,i) = nwfa_p(i,k,j)
                   enddo
                   enddo
                   enddo

                case default
             end select mp3_select

          case default
       end select mp2_select

    case("mp_p3")

   call mpas_pool_get_config(configs,'config_p3_nCat'    ,p3_nCat    )
   call mpas_pool_get_config(configs,'config_p3_trplMomI',p3_trplMomI)
   call mpas_pool_get_config(configs,'config_p3_liqfrac' ,p3_liqfrac )

! Bulk hydrometeors + number concentrations (match MPAS bulk scalar layout)
! qi, qs, qg and ni are not allocated with mp_p3 (the P3 categories were written to scalars in the driver):
! qi and ni point to the updated sums over the P3 ice categories and qs and qg to zero, for the tendency
! diagnostics below. Nothing is written to these fields.
call p3_sum_ice(state,scalars,its,ite)
qi => p3_qi_sum
qs => p3_zero_qsqg
qg => p3_zero_qsqg
ni => p3_ni_sum

call mpas_pool_get_dimension(state,'index_nr',index_nr)
nr => scalars(index_nr,:,:)

do j = jts,jte
do k = kts,kte
do i = its,ite
   nr(k,i) = max(0._RKIND, nr_p(i,k,j))
end do
end do
end do



! Aerosol-aware outputs/diagnostics (reuse MPAS aerosol-aware aerosol-aware tracers/diagnostics).
! nwfa/nifa (and nwfa2d/nifa2d) exist only with predicted aerosols (f_nifa and f_nwfa true):
call mpas_pool_get_dimension(state,'index_nc'  ,index_nc  )
nc   => scalars(index_nc,:,:)
if (f_nifa .and. f_nwfa) then
   call mpas_pool_get_dimension(state,'index_nifa',index_nifa)
   call mpas_pool_get_dimension(state,'index_nwfa',index_nwfa)
   nifa => scalars(index_nifa,:,:)
   nwfa => scalars(index_nwfa,:,:)
end if

call mpas_pool_get_array(diag_physics,'nifa2d'  ,nifa2d  )
call mpas_pool_get_array(diag_physics,'nwfa2d'  ,nwfa2d  )
call mpas_pool_get_array(diag_physics,'rainprod',rainprod)
call mpas_pool_get_array(diag_physics,'evapprod',evapprod)
call mpas_pool_get_array(diag_physics,'re_cloud',re_cloud)
call mpas_pool_get_array(diag_physics,'re_ice'  ,re_ice  )
call mpas_pool_get_array(diag_physics,'re_snow' ,re_snow )

do j = jts,jte
if (f_nifa .and. f_nwfa) then
   do i = its,ite
      nifa2d(i) = nifa2d_p(i,j)
      nwfa2d(i) = nwfa2d_p(i,j)
   enddo
   do k = kts,kte
   do i = its,ite
      nifa(k,i)    = nifa_p(i,k,j)
      nwfa(k,i)    = nwfa_p(i,k,j)
   enddo
   enddo
end if
do k = kts,kte
do i = its,ite
   nc(k,i)      = nc_p(i,k,j)
   rainprod(k,i)= rainprod_p(i,k,j)
   evapprod(k,i)= evapprod_p(i,k,j)

   ! Keep radii physically non-negative to avoid downstream radiation errors
   re_cloud(k,i)= max(0._RKIND, recloud_p(i,k,j))
   re_ice(k,i)  = max(0._RKIND, reice_p(i,k,j))
   re_snow(k,i) = max(0._RKIND, resnow_p(i,k,j))
enddo
enddo
enddo

   call mpas_pool_get_dimension(state,'index_p3_qit',index_p3_qit)
   call mpas_pool_get_dimension(state,'index_p3_qni',index_p3_qni)
   call mpas_pool_get_dimension(state,'index_p3_qir',index_p3_qir)
   call mpas_pool_get_dimension(state,'index_p3_qib',index_p3_qib)

   if (p3_trplMomI) then
      call mpas_pool_get_dimension(state,'index_p3_qzi',index_p3_qzi)
   end if
   if (p3_liqfrac) then
      call mpas_pool_get_dimension(state,'index_p3_qli',index_p3_qli)
   end if

   ! P3DBG_TO: map and range-check P3 scalar indices against total scalar count (print/check once)
   if (.not. p3dbg_to_indices_checked) then
      nScalarTot_loc = nScalar1_loc
      if (associated(p3_debug_on) .and. p3_debug_on) then
         write(dbgmsg,'(A,4I8)') 'P3DBG_TO: index_p3_qit,qni,qir,qib=', &
                                index_p3_qit, index_p3_qni, index_p3_qir, index_p3_qib
         call mpas_log_write(trim(dbgmsg))
         if (p3_trplMomI) then
            write(dbgmsg,'(A,I8)') 'P3DBG_TO: index_p3_qzi=', index_p3_qzi
            call mpas_log_write(trim(dbgmsg))
         end if
         if (p3_liqfrac) then
            write(dbgmsg,'(A,I8)') 'P3DBG_TO: index_p3_qli=', index_p3_qli
            call mpas_log_write(trim(dbgmsg))
         end if
      end if

      p3_index_ok = .true.
      if (index_p3_qit < 1 .or. index_p3_qit > nScalarTot_loc) p3_index_ok = .false.
      if (index_p3_qni < 1 .or. index_p3_qni > nScalarTot_loc) p3_index_ok = .false.
      if (index_p3_qir < 1 .or. index_p3_qir > nScalarTot_loc) p3_index_ok = .false.
      if (index_p3_qib < 1 .or. index_p3_qib > nScalarTot_loc) p3_index_ok = .false.
      if (p3_trplMomI) then
         if (index_p3_qzi < 1 .or. index_p3_qzi > nScalarTot_loc) p3_index_ok = .false.
      end if
      if (p3_liqfrac) then
         if (index_p3_qli < 1 .or. index_p3_qli > nScalarTot_loc) p3_index_ok = .false.
      end if

      ! Also validate that each P3 category block (base .. base+p3_nCat-1) fits in scalars(:,:,:)
      if (index_p3_qit + p3_nCat - 1 > nScalarTot_loc) p3_index_ok = .false.
      if (index_p3_qni + p3_nCat - 1 > nScalarTot_loc) p3_index_ok = .false.
      if (index_p3_qir + p3_nCat - 1 > nScalarTot_loc) p3_index_ok = .false.
      if (index_p3_qib + p3_nCat - 1 > nScalarTot_loc) p3_index_ok = .false.
      if (p3_trplMomI) then
         if (index_p3_qzi + p3_nCat - 1 > nScalarTot_loc) p3_index_ok = .false.
      end if
      if (p3_liqfrac) then
         if (index_p3_qli + p3_nCat - 1 > nScalarTot_loc) p3_index_ok = .false.
      end if

      if (p3_index_ok) then
         if (associated(p3_debug_on) .and. p3_debug_on) then
            write(dbgmsg,'(A,I8)') 'P3DBG_TO: P3 scalar indices in-range; nScalarTot=', nScalarTot_loc
            call mpas_log_write(trim(dbgmsg))
         end if
      else
         write(dbgmsg,'(A,I8)') 'P3DBG_TO: ERROR P3 scalar index out-of-range; nScalarTot=', nScalarTot_loc
         call mpas_log_write(trim(dbgmsg))
         stop 1
      end if
      p3dbg_to_indices_checked = .true.
   end if
   ! P3 categories are written back directly to state%scalars in the driver (MPAS-compatible; no 4-D intermediates).
    case default
 end select mp_select

!end calculation of cloud microphysics tendencies:
! P3 now follows the same before/after differencing contract as Thompson.
 mp_tend_select: select case(trim(mp_scheme))
    case("mp_thompson","mp_thompson_aerosols","mp_p3","mp_wsm6")
       call mpas_pool_get_array(tend_physics,'rthmpten',rthmpten)
       call mpas_pool_get_array(tend_physics,'rqvmpten',rqvmpten)
       call mpas_pool_get_array(tend_physics,'rqcmpten',rqcmpten)
       call mpas_pool_get_array(tend_physics,'rqrmpten',rqrmpten)
       call mpas_pool_get_array(tend_physics,'rqimpten',rqimpten)
       call mpas_pool_get_array(tend_physics,'rqsmpten',rqsmpten)
       call mpas_pool_get_array(tend_physics,'rqgmpten',rqgmpten)

       do i = its,ite
       do k = kts,mp_top_level
          rthmpten(k,i) = (theta_m(k,i)/(1._RKIND+R_v/R_d*max(0._RKIND,qv(k,i)))-rthmpten(k,i))/dt_dyn
          rqvmpten(k,i) = (qv(k,i)-rqvmpten(k,i))/dt_dyn
          rqcmpten(k,i) = (qc(k,i)-rqcmpten(k,i))/dt_dyn
          rqrmpten(k,i) = (qr(k,i)-rqrmpten(k,i))/dt_dyn
          rqimpten(k,i) = (qi(k,i)-rqimpten(k,i))/dt_dyn
          rqsmpten(k,i) = (qs(k,i)-rqsmpten(k,i))/dt_dyn
          rqgmpten(k,i) = (qg(k,i)-rqgmpten(k,i))/dt_dyn
       enddo
       enddo

       ! no microphysics tendency above level mp_top_level
       if(kte .gt. mp_top_level) then
          do i = its, ite
          do k = mp_top_level+1, kte
             rthmpten(k,i) = 0.
             rqvmpten(k,i) = 0.
             rqcmpten(k,i) = 0.
             rqrmpten(k,i) = 0.
             rqimpten(k,i) = 0.
             rqsmpten(k,i) = 0.
             rqgmpten(k,i) = 0.
          enddo
          enddo
       endif

       mp2_tend_select: select case(trim(mp_scheme))
          case("mp_thompson","mp_thompson_aerosols","mp_p3")
             call mpas_pool_get_array(tend_physics,'rnimpten',rnimpten)
             call mpas_pool_get_array(tend_physics,'rnrmpten',rnrmpten)

             do k = kts,kte
             do i = its,ite
                rnimpten(k,i) = (ni(k,i)-rnimpten(k,i))/dt_dyn
                rnrmpten(k,i) = (nr(k,i)-rnrmpten(k,i))/dt_dyn
             enddo
             enddo

             mp3_tend_select: select case(trim(mp_scheme))
                case("mp_thompson_aerosols","mp_p3")
                   !mp_p3: nifa and nwfa exist only with predicted aerosols (config_p3_aerosols):
                   call mpas_pool_get_array(tend_physics,'rncmpten',rncmpten)

                   do k = kts,kte
                   do i = its,ite
                      rncmpten(k,i)   = (nc(k,i)-rncmpten(k,i))/dt_dyn
                   enddo
                   enddo
                   if(f_nifa .and. f_nwfa) then
                   call mpas_pool_get_array(tend_physics,'rnifampten',rnifampten)
                   call mpas_pool_get_array(tend_physics,'rnwfampten',rnwfampten)
                   do k = kts,kte
                   do i = its,ite
                      rnifampten(k,i) = (nifa(k,i)-rnifampten(k,i))/dt_dyn
                      rnwfampten(k,i) = (nwfa(k,i)-rnwfampten(k,i))/dt_dyn
                   enddo
                   enddo
                   endif

                case default
             end select mp3_tend_select

          case default
       end select mp2_tend_select

    case default
 end select mp_tend_select

 end subroutine microphysics_to_MPAS

!=================================================================================================================
 subroutine post_microphysics(configs,state,diag,tend,time_lev)
!=================================================================================================================

!input variables:
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: state
 type(mpas_pool_type),intent(in):: diag
 type(mpas_pool_type),intent(inout):: tend


 integer:: time_lev

#ifdef MPAS_OPENACC
!local pointers:
 real(kind=RKIND),dimension(:,:),pointer  :: exner,exner_b,pressure_b,rtheta_p,rtheta_b
 real(kind=RKIND),dimension(:,:),pointer  :: rho_zz,theta_m,pressure_p
 real(kind=RKIND),dimension(:,:,:),pointer:: scalars
 real(kind=RKIND),dimension(:,:),pointer  :: rt_diabatic_tend

 call mpas_pool_get_array(diag,'exner'        ,exner     )
 call mpas_pool_get_array(diag,'exner_base'   ,exner_b   )
 call mpas_pool_get_array(diag,'pressure_base',pressure_b)
 call mpas_pool_get_array(diag,'pressure_p'   ,pressure_p)
 call mpas_pool_get_array(diag,'rtheta_base'  ,rtheta_b  )
 call mpas_pool_get_array(diag,'rtheta_p'     ,rtheta_p  )

 call mpas_pool_get_array(state,'rho_zz' ,rho_zz ,time_lev)
 call mpas_pool_get_array(state,'theta_m',theta_m,time_lev)

 call mpas_pool_get_array(state,'scalars',scalars,time_lev)

 call mpas_pool_get_array(tend,'rt_diabatic_tend',rt_diabatic_tend)


 MPAS_ACC_TIMER_START('post_microphysics [ACC_data_xfer]')
 !$acc update device(exner, exner_b, pressure_b, pressure_p, rtheta_b)
 !$acc update device(rtheta_p, rho_zz, theta_m, scalars)
 !$acc update device(rt_diabatic_tend)
 MPAS_ACC_TIMER_STOP('post_microphysics [ACC_data_xfer]')
#endif

end subroutine post_microphysics

!=================================================================================================================
 end module mpas_atmphys_interface
!=================================================================================================================
