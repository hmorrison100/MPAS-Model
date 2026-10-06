! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module mpas_atmphys_p3_cumix
 use mpas_kind_types
 use mpas_pool_routines
 use mpas_constants,only: pii
 use microphy_p3,only: icecat_destination,G_of_mu,p3_deltaD_init

 implicit none
 private
 public:: p3cu_add_tend

!mpas_atmphys_p3_cumix couples the P3 cloud microphysics with the convection schemes. The convection schemes
!(cu_ntiedtke, cu_tiedtke, cu_grell_freitas, cu_kain_fritsch) return the tendencies of the mass mixing ratios
!of the condensate detrained into the grid-scale environment (rqccuten, rqicuten and, for cu_kain_fritsch,
!rqrcuten and rqscuten), but no number tendencies. With P3:
! * cloud water: rqccuten is added to qc (in physics_get_tend_work). If config_p3_predictNc, the corresponding
!   number tendency is added to nc assuming detrained droplets of radius config_p3_cu_rc.
! * rain       : rqrcuten (cu_kain_fritsch only) is added to qr (in physics_get_tend_work). The corresponding
!   number tendency is added to nr assuming detrained drops of mean-volume radius config_p3_cu_rr.
! * cloud ice  : P3 has no bulk qi and qs. The ice tendency qiten = rqicuten (plus rqscuten for cu_kain_fritsch;
!   for P3, f_qs is false so that cu_kain_fritsch already includes its snow in rqicuten) is added to the P3 ice
!   scalars:
!   - qiten > 0 (detrained ice): the detrained ice is added to a single ice category, with the number tendency
!     added to p3_qni assuming detrained ice particles of mean-volume radius config_p3_cu_ri. The detrained ice
!     is unrimed (no p3_qir, p3_qib tendency) and dry (no p3_qli tendency). For config_p3_nCat > 1, the
!     destination category is chosen as for new ice in P3 (icecat_destination): the category whose mean
!     diameter of the existing ice is closest to that of the detrained ice, or a new (empty) category if the
!     difference exceeds the P3 threshold deltaD_init for all non-empty categories. The mean diameter of the
!     existing ice is the P3 volume-equivalent diameter (6 qit/(pi rhoi qni))^(1/3), with the bulk ice density
!     rhoi from the last P3 call (p3_rhoi). For triple-moment ice, the change of the 6th moment Z is that of
!     new ice in P3 (update_zi_proc2) with mu_i = 0 (as mu_r_constant for all new ice in P3):
!     dZ = G(mu) dM3**2/dN, which, for detrained particles of a single size Dv = 2*config_p3_cu_ri, gives
!     dZ = G(mu) 6/(pi rho_i) Dv**3 dQ. The advected scalar p3_qzi is not Z but the transformed variable
!     Y = sqrt(Z N) (P3 converts Z = Y**2/N at the start of p3_main and back at the end); its tendency is
!     (Y_new - Y)/dt, with Y_new = sqrt((Y**2/N + dZ)*(N + dN)) and dQ, dN, dZ the changes over one time step.
!   - qiten < 0 (only cu_kain_fritsch, through the compensating subsidence of the existing ice; the other
!     schemes only detrain ice): the ice is removed from all categories in proportion to their mass, and
!     all P3 scalars of a category (p3_qni, p3_qir, p3_qib and p3_qzi, p3_qli when active) are reduced by the
!     same fraction as p3_qit, so that the mean size, density, rime fraction, shape (mu_i) and liquid fraction
!     of each category are unchanged (p3_qzi = sqrt(Z N) scales by the same fraction as Z and N). The
!     removal over one time step is limited to the existing ice.
!Except for this limit, the mass tendencies are those of the convection schemes, so that the total detrained
!mass and latent heating are the same as with the other microphysics schemes.

 real(kind=RKIND),parameter:: rho_w = 1000._RKIND  !density of detrained droplets and rain (kg m-3), as rhow in P3.
 real(kind=RKIND),parameter:: rho_i =  900._RKIND  !density of detrained ice (kg m-3), as new ice in P3.
 real(kind=RKIND),parameter:: mu_i_cu =  0._RKIND  !mu_i of detrained ice (p3_qzi tendency), as new ice in P3.
 real(kind=RKIND),parameter:: qsmall  = 1.e-14_RKIND !smallest ice mass mixing ratio (kg kg-1), as in P3.
 real(kind=RKIND),parameter:: nsmall  = 1.e-16_RKIND !smallest ice number mixing ratio (kg-1), as in P3.
 integer,parameter:: nCat_max = 3                   !maximum number of P3 ice categories in MPAS.


 contains


!=================================================================================================================
 subroutine p3cu_add_tend(configs,state,diag_physics,tend_physics,mass,scalars,tend_scalars,nCellsSolve,nVertLevels)
!=================================================================================================================
!adds the tendencies of nc and of the P3 ice scalars due to the convective detrainment of cloud water and ice
!to the scalar tendencies.

!input arguments:
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: state
 type(mpas_pool_type),intent(in):: diag_physics
 type(mpas_pool_type),intent(in):: tend_physics
 integer,intent(in):: nCellsSolve,nVertLevels
 real(kind=RKIND),intent(in),dimension(:,:):: mass
 real(kind=RKIND),intent(in),dimension(:,:,:):: scalars  !time level 1, as the tendencies are added once per step.

!inout arguments:
 real(kind=RKIND),intent(inout),dimension(:,:,:):: tend_scalars

!local variables and pointers:
 character(len=StrKIND),pointer:: convection_scheme
 logical,pointer:: p3_predictNc,p3_trplMomI,p3_liqfrac
 integer,pointer:: p3_nCat
 integer,pointer:: index_nc,index_nr,index_p3_qit,index_p3_qni,index_p3_qir,index_p3_qib,index_p3_qzi, &
                   index_p3_qli
 real(kind=RKIND),pointer:: p3_cu_rc,p3_cu_rr,p3_cu_ri,dt
 real(kind=RKIND),dimension(:,:),pointer:: rqccuten,rqrcuten,rqicuten,rqscuten
 real(kind=RKIND),dimension(:,:,:),pointer:: p3_rhoi

 integer:: c,i,k,n,nscal,icat
 integer,dimension(6*nCat_max):: idx
 real(kind=RKIND):: rmass_c,rmass_r,rmass_i,qiten
 real(kind=RKIND):: d_cu,deltaD_init,zfac,qsum,rate,rho_c
 real(kind=RKIND):: dq,n_old,y_old,z_new,n_new
 real(kind=RKIND),dimension(nCat_max):: qi_cat,di_cat
 logical:: add_qr,add_qs

!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_config(configs,'config_convection_scheme',convection_scheme)
 call mpas_pool_get_config(configs,'config_p3_nCat'          ,p3_nCat          )
 call mpas_pool_get_config(configs,'config_p3_trplMomI'      ,p3_trplMomI      )
 call mpas_pool_get_config(configs,'config_p3_liqfrac'       ,p3_liqfrac       )
 call mpas_pool_get_config(configs,'config_p3_predictNc'     ,p3_predictNc     )
 call mpas_pool_get_config(configs,'config_p3_cu_rc'         ,p3_cu_rc         )
 call mpas_pool_get_config(configs,'config_p3_cu_rr'         ,p3_cu_rr         )
 call mpas_pool_get_config(configs,'config_p3_cu_ri'         ,p3_cu_ri         )
 call mpas_pool_get_config(configs,'config_dt'               ,dt               )

 call mpas_pool_get_dimension(state,'index_nc'    ,index_nc    )
 call mpas_pool_get_dimension(state,'index_nr'    ,index_nr    )
 call mpas_pool_get_dimension(state,'index_p3_qit',index_p3_qit)
 call mpas_pool_get_dimension(state,'index_p3_qni',index_p3_qni)
 call mpas_pool_get_dimension(state,'index_p3_qir',index_p3_qir)
 call mpas_pool_get_dimension(state,'index_p3_qib',index_p3_qib)
 call mpas_pool_get_dimension(state,'index_p3_qzi',index_p3_qzi)
 call mpas_pool_get_dimension(state,'index_p3_qli',index_p3_qli)

 nullify(rqccuten)
 nullify(rqrcuten)
 nullify(rqicuten)
 nullify(rqscuten)
 call mpas_pool_get_array(tend_physics,'rqccuten',rqccuten)
 call mpas_pool_get_array(tend_physics,'rqrcuten',rqrcuten)
 call mpas_pool_get_array(tend_physics,'rqicuten',rqicuten)
 call mpas_pool_get_array(tend_physics,'rqscuten',rqscuten)
 add_qr = trim(convection_scheme) == 'cu_kain_fritsch' .and. associated(rqrcuten) .and. index_nr > 0
 add_qs = trim(convection_scheme) == 'cu_kain_fritsch' .and. associated(rqscuten)

!inverse mass of a detrained droplet, rain drop and ice particle (kg-1):
 rmass_c = 1._RKIND / (4._RKIND/3._RKIND*pii*rho_w*p3_cu_rc**3)
 rmass_r = 1._RKIND / (4._RKIND/3._RKIND*pii*rho_w*p3_cu_rr**3)
 rmass_i = 1._RKIND / (4._RKIND/3._RKIND*pii*rho_i*p3_cu_ri**3)

!cloud water number:
 if(p3_predictNc .and. index_nc > 0 .and. associated(rqccuten)) then
    do i = 1, nCellsSolve
    do k = 1, nVertLevels
       tend_scalars(index_nc,k,i) = tend_scalars(index_nc,k,i) + max(0._RKIND,rqccuten(k,i))*rmass_c*mass(k,i)
    enddo
    enddo
 endif

!rain number (cu_kain_fritsch):
 if(add_qr) then
    do i = 1, nCellsSolve
    do k = 1, nVertLevels
       tend_scalars(index_nr,k,i) = tend_scalars(index_nr,k,i) + max(0._RKIND,rqrcuten(k,i))*rmass_r*mass(k,i)
    enddo
    enddo
 endif

 if(.not. associated(rqicuten)) return

!cloud ice:
!indices of the P3 scalars of category c: idx(c+(n-1)*p3_nCat), n = 1 (p3_qit), 2 (p3_qni), 3 (p3_qir),
!4 (p3_qib), then p3_qzi and p3_qli when active (nscal scalars per category):
 nscal = 4
 if(p3_trplMomI) nscal = nscal + 1
 if(p3_liqfrac ) nscal = nscal + 1
 do c = 1, p3_nCat
    idx(c          ) = index_p3_qit + c - 1
    idx(c+  p3_nCat) = index_p3_qni + c - 1
    idx(c+2*p3_nCat) = index_p3_qir + c - 1
    idx(c+3*p3_nCat) = index_p3_qib + c - 1
    n = 4
    if(p3_trplMomI) then
       n = n + 1
       idx(c+(n-1)*p3_nCat) = index_p3_qzi + c - 1
    endif
    if(p3_liqfrac) then
       n = n + 1
       idx(c+(n-1)*p3_nCat) = index_p3_qli + c - 1
    endif
 enddo

!destination category of the detrained ice (config_p3_nCat > 1): diameter of the detrained ice and P3 threshold
!size difference for a new category; bulk density of the existing ice from the last P3 call:
 nullify(p3_rhoi)
 if(p3_nCat > 1) then
    d_cu = 2._RKIND*p3_cu_ri
    deltaD_init = p3_deltaD_init(p3_nCat)
    call mpas_pool_get_array(diag_physics,'p3_rhoi',p3_rhoi)
 endif

!triple-moment ice: change of the 6th moment Z per unit detrained ice mass (m6 kg-1 per kg kg-1):
 zfac = 0._RKIND
 if(p3_trplMomI) zfac = G_of_mu(mu_i_cu)*6._RKIND/(pii*rho_i)*(2._RKIND*p3_cu_ri)**3

 do i = 1, nCellsSolve
 do k = 1, nVertLevels
    qiten = rqicuten(k,i)
    if(add_qs) qiten = qiten + rqscuten(k,i)

    if(qiten > 0._RKIND) then
       !detrained ice, added to the destination category icat:
       icat = 1
       if(p3_nCat > 1) then
          do c = 1, p3_nCat
             qi_cat(c) = max(0._RKIND,scalars(idx(c),k,i))
             di_cat(c) = 0._RKIND
             if(qi_cat(c) >= qsmall) then
                !bulk density from the last P3 call; rho_i if the category was empty then (e.g., ice advected or
                !detrained since):
                rho_c = p3_rhoi(c,k,i)
                if(rho_c <= 0._RKIND) rho_c = rho_i
                di_cat(c) = (6._RKIND*qi_cat(c)/(pii*rho_c*max(nsmall,scalars(idx(c+p3_nCat),k,i))))**(1._RKIND/3._RKIND)
             endif
          enddo
          call icecat_destination(qi_cat(1:p3_nCat),di_cat(1:p3_nCat),d_cu,deltaD_init,icat)
       endif
       tend_scalars(idx(icat),k,i) = tend_scalars(idx(icat),k,i) + qiten*mass(k,i)
       tend_scalars(idx(icat+p3_nCat),k,i) = tend_scalars(idx(icat+p3_nCat),k,i) + qiten*rmass_i*mass(k,i)
       if(p3_trplMomI) then
          !p3_qzi is Y = sqrt(Z N): Z of the existing ice from Y and N, plus that of the detrained ice:
          dq    = qiten*dt
          n_old = max(0._RKIND,scalars(idx(icat+p3_nCat),k,i))
          y_old = max(0._RKIND,scalars(index_p3_qzi+icat-1,k,i))
          z_new = zfac*dq
          if(n_old >= nsmall) z_new = z_new + (y_old/n_old)*y_old
          n_new = n_old + dq*rmass_i
          tend_scalars(index_p3_qzi+icat-1,k,i) = tend_scalars(index_p3_qzi+icat-1,k,i) &
             + (sqrt(z_new)*sqrt(n_new) - y_old)/dt*mass(k,i)
       endif

    elseif(qiten < 0._RKIND) then
       !removal of existing ice (cu_kain_fritsch), proportional to the mass of each category and limited to the
       !existing ice over one time step; all scalars of a category are scaled by the same fraction:
       qsum = 0._RKIND
       do c = 1, p3_nCat
          qsum = qsum + max(0._RKIND,scalars(idx(c),k,i))
       enddo
       if(qsum >= qsmall) then
          rate = max(qiten,-qsum/dt)/qsum   !fractional removal rate (s-1), between -1/dt and 0.
          do c = 1, p3_nCat
             if(scalars(idx(c),k,i) > 0._RKIND) then
                do n = 1, nscal
                   tend_scalars(idx(c+(n-1)*p3_nCat),k,i) = tend_scalars(idx(c+(n-1)*p3_nCat),k,i) &
                      + rate*max(0._RKIND,scalars(idx(c+(n-1)*p3_nCat),k,i))*mass(k,i)
                enddo
             endif
          enddo
       endif
    endif
 enddo
 enddo

 end subroutine p3cu_add_tend

!=================================================================================================================
 end module mpas_atmphys_p3_cumix
!=================================================================================================================
