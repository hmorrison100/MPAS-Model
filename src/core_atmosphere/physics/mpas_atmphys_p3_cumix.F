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
 use mpas_atmphys_utilities,only: physics_error_fatal

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
! * cloud ice  : P3 has no bulk qi and qs. The detrained ice (rqicuten, plus rqscuten for cu_kain_fritsch) is
!   added to p3_qit with the number tendency added to p3_qni, assuming detrained ice particles of mean-volume
!   radius config_p3_cu_ri. The detrained ice is unrimed (no p3_qir, p3_qib tendency) and dry (no p3_qli
!   tendency). For P3, f_qs is false so that cu_kain_fritsch already includes its snow in rqicuten.
!The mass tendencies are those of the convection schemes, so that the total detrained mass and latent heating
!are the same as with the other microphysics schemes.
!Only config_p3_nCat = 1 and double-moment ice are supported (the destination category for config_p3_nCat > 1
!and the p3_qzi tendency for triple-moment ice are not implemented yet).

 real(kind=RKIND),parameter:: rho_w = 1000._RKIND  !density of detrained droplets and rain (kg m-3), as rhow in P3.
 real(kind=RKIND),parameter:: rho_i =  900._RKIND  !density of detrained ice (kg m-3), as new ice in P3.


 contains


!=================================================================================================================
 subroutine p3cu_add_tend(configs,state,tend_physics,mass,tend_scalars,nCellsSolve,nVertLevels)
!=================================================================================================================
!adds the tendencies of nc and of the P3 ice scalars due to the convective detrainment of cloud water and ice
!to the scalar tendencies.

!input arguments:
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: state
 type(mpas_pool_type),intent(in):: tend_physics
 integer,intent(in):: nCellsSolve,nVertLevels
 real(kind=RKIND),intent(in),dimension(:,:):: mass

!inout arguments:
 real(kind=RKIND),intent(inout),dimension(:,:,:):: tend_scalars

!local variables and pointers:
 character(len=StrKIND),pointer:: convection_scheme
 logical,pointer:: p3_cu_detrain,p3_predictNc,p3_trplMomI
 integer,pointer:: p3_nCat
 integer,pointer:: index_nc,index_nr,index_p3_qit,index_p3_qni
 real(kind=RKIND),pointer:: p3_cu_rc,p3_cu_rr,p3_cu_ri
 real(kind=RKIND),dimension(:,:),pointer:: rqccuten,rqrcuten,rqicuten,rqscuten

 integer:: i,k
 real(kind=RKIND):: rmass_c,rmass_r,rmass_i,qiten
 logical:: add_qr,add_qs

!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_config(configs,'config_p3_cu_detrain',p3_cu_detrain)
 if(.not. p3_cu_detrain) return

 call mpas_pool_get_config(configs,'config_convection_scheme',convection_scheme)
 call mpas_pool_get_config(configs,'config_p3_nCat'          ,p3_nCat          )
 call mpas_pool_get_config(configs,'config_p3_trplMomI'      ,p3_trplMomI      )
 call mpas_pool_get_config(configs,'config_p3_predictNc'     ,p3_predictNc     )
 call mpas_pool_get_config(configs,'config_p3_cu_rc'         ,p3_cu_rc         )
 call mpas_pool_get_config(configs,'config_p3_cu_rr'         ,p3_cu_rr         )
 call mpas_pool_get_config(configs,'config_p3_cu_ri'         ,p3_cu_ri         )

 if(p3_nCat > 1 .or. p3_trplMomI) call physics_error_fatal('p3cu_add_tend: convective detrainment into P3 ' // &
    'is only implemented for config_p3_nCat = 1 and double-moment ice; set config_p3_cu_detrain = false')

 call mpas_pool_get_dimension(state,'index_nc'    ,index_nc    )
 call mpas_pool_get_dimension(state,'index_nr'    ,index_nr    )
 call mpas_pool_get_dimension(state,'index_p3_qit',index_p3_qit)
 call mpas_pool_get_dimension(state,'index_p3_qni',index_p3_qni)

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

!cloud ice mass and number (ice category 1):
 if(associated(rqicuten)) then
    do i = 1, nCellsSolve
    do k = 1, nVertLevels
       qiten = rqicuten(k,i)
       if(add_qs) qiten = qiten + rqscuten(k,i)
       tend_scalars(index_p3_qit,k,i) = tend_scalars(index_p3_qit,k,i) + qiten*mass(k,i)
       tend_scalars(index_p3_qni,k,i) = tend_scalars(index_p3_qni,k,i) + max(0._RKIND,qiten)*rmass_i*mass(k,i)
    enddo
    enddo
 endif

 end subroutine p3cu_add_tend

!=================================================================================================================
 end module mpas_atmphys_p3_cumix
!=================================================================================================================
