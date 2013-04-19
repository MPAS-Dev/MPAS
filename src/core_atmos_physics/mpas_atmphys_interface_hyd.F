!=============================================================================================
 module mpas_atmphys_interface_hyd
 use mpas_configure
 use mpas_grid_types

 use mpas_atmphys_constants
 use mpas_atmphys_vars

 implicit none
 private
 public:: allocate_forall_physics,   &
          deallocate_forall_physics, &
          MPAS_to_physics,           &
          microphysics_from_MPAS,    &
          microphysics_to_MPAS

!... local variables:
 integer,private:: i,k,j
 real(kind=RKIND),dimension(:,:),pointer,private:: u,v,w
 real(kind=RKIND),dimension(:,:),pointer,private:: geopotential,pressure,theta,qv

 contains

!=============================================================================================
 subroutine allocate_forall_physics
!=============================================================================================

 if(.not.allocated(psfc_p))  allocate(psfc_p(ims:ime,jms:jme)           )
 if(.not.allocated(ptop_p))  allocate(ptop_p(ims:ime,jms:jme)           )

 if(.not.allocated(u_p)    )  allocate(u_p(ims:ime,kms:kme,jms:jme)     )
 if(.not.allocated(v_p)    )  allocate(v_p(ims:ime,kms:kme,jms:jme)     )
 if(.not.allocated(zz_p)   )  allocate(zz_p(ims:ime,kms:kme,jms:jme)    )
 if(.not.allocated(pres_p) )  allocate(pres_p(ims:ime,kms:kme,jms:jme)  )
 if(.not.allocated(pi_p)   )  allocate(pi_p(ims:ime,kms:kme,jms:jme)    )
 if(.not.allocated(z_p)    )  allocate(z_p(ims:ime,kms:kme,jms:jme)     )
 if(.not.allocated(dz_p)   )  allocate(dz_p(ims:ime,kms:kme,jms:jme)    )
 if(.not.allocated(t_p)    )  allocate(t_p(ims:ime,kms:kme,jms:jme)     )
 if(.not.allocated(th_p)   )  allocate(th_p(ims:ime,kms:kme,jms:jme)    )
 if(.not.allocated(al_p)   )  allocate(al_p(ims:ime,kms:kme,jms:jme)    )
 if(.not.allocated(rho_p)  )  allocate(rho_p(ims:ime,kms:kme,jms:jme)   )
 if(.not.allocated(rh_p)   )  allocate(rh_p(ims:ime,kms:kme,jms:jme)    )

 if(.not.allocated(w_p)    )  allocate(w_p(ims:ime,kms:kme,jms:jme)     )
 if(.not.allocated(pres2_p))  allocate(pres2_p(ims:ime,kms:kme,jms:jme) )
 if(.not.allocated(t2_p)   )  allocate(t2_p(ims:ime,kms:kme,jms:jme)    )

 if(.not.allocated(pres_hyd_p) ) allocate(pres_hyd_p(ims:ime,kms:kme,jms:jme)  )
 if(.not.allocated(pres2_hyd_p)) allocate(pres2_hyd_p(ims:ime,kms:kme,jms:jme) )
 
 do j = jms,jme
 do i = ims,ime
    psfc_p(i,j)   = 0.
    ptop_p(i,j)   = 0.
 enddo
 enddo

 do j = jms,jme
 do k = kms,kme
 do i = ims,ime
    u_p(i,k,j)     = 0.
    v_p(i,k,j)     = 0.
    w_p(i,k,j)     = 0.
    pres_p(i,k,j)  = 0.
    pi_p(i,k,j)    = 0.
    z_p(i,k,j)     = 0.
    dz_p(i,k,j)    = 0.
    t_p(i,k,j)     = 0.
    th_p(i,k,j)    = 0.
    al_p(i,k,j)    = 0.
    rho_p(i,k,j)   = 0.
    rh_p(i,k,j)    = 0.

    w_p(i,k,j)     = 0.
    pres2_p(i,k,j) = 0.
    t2_p(i,k,j)    = 0.

    pres_hyd_p(i,k,j)  = 0.
    pres2_hyd_p(i,k,j) = 0.
 enddo
 enddo
 enddo

!allocate moist species (to be revisited!):
 if(.not.allocated(qv_p)   ) allocate(qv_p(ims:ime,kms:kme,jms:jme)     )
 if(.not.allocated(qc_p)   ) allocate(qc_p(ims:ime,kms:kme,jms:jme)     )
 if(.not.allocated(qr_p)   ) allocate(qr_p(ims:ime,kms:kme,jms:jme)     )
 if(.not.allocated(qi_p)   ) allocate(qi_p(ims:ime,kms:kme,jms:jme)     )
 if(.not.allocated(qs_p)   ) allocate(qs_p(ims:ime,kms:kme,jms:jme)     )
 if(.not.allocated(qg_p)   ) allocate(qg_p(ims:ime,kms:kme,jms:jme)     )

 end subroutine allocate_forall_physics

!=============================================================================================
 subroutine deallocate_forall_physics
!=============================================================================================

!de-allocation of all physics arrays:
 if(allocated(psfc_p)  ) deallocate(psfc_p  )
 if(allocated(ptop_p)  ) deallocate(ptop_p  )

 if(allocated(u_p)     ) deallocate(u_p     )
 if(allocated(v_p)     ) deallocate(v_p     )
 if(allocated(zz_p)    ) deallocate(zz_p    )
 if(allocated(pres_p)  ) deallocate(pres_p  )
 if(allocated(pi_p)    ) deallocate(pi_p    )
 if(allocated(z_p)     ) deallocate(z_p     )
 if(allocated(dz_p)    ) deallocate(dz_p    )
 if(allocated(t_p)     ) deallocate(t_p     )
 if(allocated(th_p)    ) deallocate(th_p    )
 if(allocated(al_p)    ) deallocate(al_p    )
 if(allocated(rho_p)   ) deallocate(rho_p   ) 
 if(allocated(rh_p)    ) deallocate(rh_p    ) 

 if(allocated(w_p)     ) deallocate(w_p     )
 if(allocated(pres2_p) ) deallocate(pres2_p )
 if(allocated(t2_p)    ) deallocate(t2_p    )

 if(allocated(pres_hyd_p) ) deallocate(pres_hyd_p  )
 if(allocated(pres2_hyd_p)) deallocate(pres2_hyd_p )
 
 if(allocated(qv_p)    ) deallocate(qv_p    )
 if(allocated(qc_p)    ) deallocate(qc_p    )
 if(allocated(qr_p)    ) deallocate(qr_p    )
 if(allocated(qi_p)    ) deallocate(qi_p    )
 if(allocated(qs_p)    ) deallocate(qs_p    )
 if(allocated(qg_p)    ) deallocate(qg_p    )

 end subroutine deallocate_forall_physics

!=============================================================================================
 subroutine MPAS_to_physics(state,diag)
!=============================================================================================

!input variables:
 type(state_type),intent(in):: state
 type(diag_type),intent(in):: diag
 
!local variables:
 real(kind=RKIND):: tm

!---------------------------------------------------------------------------------------------

!initialization:
 w => state % w % array
 u => diag % uReconstructZonal % array
 v => diag % uReconstructMeridional % array

 geopotential => state % geopotential % array
 pressure     => state % pressure % array
 theta        => state % theta % array
 qv           => state % scalars % array(state%index_qv,:,:)

 do j = jts, jte
 do k = kts, kte
 do i = its, ite
    
    u_p(i,k,j) = u(k,i)
    v_p(i,k,j) = v(k,i)
    w_p(i,k,j) = w(k,i)
 
    dz_p(i,k,j)   = (geopotential(k+1,i) - geopotential(k,i)) / g
    pres_p(i,k,j) = (pressure(k+1,i) + pressure(k,i)) / 2
    th_p(i,k,j)   = theta(k,i)
    qv_p(i,k,j)   = qv(k,i)

    pi_p(i,k,j) = (pres_p(i,k,j)/p0)**(R_d/cp)
    t_p(i,k,j)  = th_p(i,k,j)*pi_p(i,k,j)

    tm = (1.+1.61*qv_p(i,k,j))*th_p(i,k,j) 
    al_p(i,k,j)  = R_d/P0*tm*(pres_p(i,k,j)/p0)**cvpm
    rho_p(i,k,j) = 1./al_p(i,k,j)
    
 enddo
 enddo
 enddo

!formats: 
 201 format(3i8,10(1x,e15.8))
 202 format(2i6,10(1x,e15.8))

 end subroutine MPAS_to_physics

!=============================================================================================
 subroutine microphysics_from_MPAS(state)
!=============================================================================================

!input variables:
 type(state_type),intent(in):: state

!---------------------------------------------------------------------------------------------

!initialization:
 write(0,*)
 write(0,*) '--- enter subroutine microphysics_from_MPAS:'

 geopotential => state % geopotential % array
 pressure     => state % pressure % array
 theta        => state % theta % array

!copy sounding variables from the geodesic grid to the wrf-physics grid:
 do j = jts, jte
 do k = kts, kte
 do i = its, ite
    dz_p(i,k,j)   = (geopotential(k+1,i) - geopotential(k,i)) / g
    pres_p(i,k,j) = (pressure(k+1,i) +  pressure(k,i) ) / 2
    th_p(i,k,j)   = theta(k,i)

    pi_p(i,k,j) = (pres_p(i,k,j)/p0)**(R_d/cp)
    t_p(i,k,j)  = th_p(i,k,j)*pi_p(i,k,j)

!   qv_p(i,k,j) = max(0.,state % scalars % array(state%index_qv,k,i))
!   qc_p(i,k,j) = max(0.,state % scalars % array(state%index_qc,k,i))
!   qr_p(i,k,j) = max(0.,state % scalars % array(state%index_qr,k,i))
    qv_p(i,k,j) = state % scalars % array(state%index_qv,k,i)
    qc_p(i,k,j) = state % scalars % array(state%index_qc,k,i)
    qr_p(i,k,j) = state % scalars % array(state%index_qr,k,i)
 enddo
 enddo
 enddo

!additional initialization as function of cloud microphysics scheme:
 microp_select_init: select case(microp_scheme)

    case ("thompson")

       do j = jts, jte
       do k = kts, kte
       do i = its, ite
          !mass mixing ratios:
!         qi_p(i,k,j) = max(0.,state % scalars % array(state%index_qi,k,i))
!         qs_p(i,k,j) = max(0.,state % scalars % array(state%index_qs,k,i))
!         qg_p(i,k,j) = max(0.,state % scalars % array(state%index_qg,k,i))
          qi_p(i,k,j) = state % scalars % array(state%index_qi,k,i)
          qs_p(i,k,j) = state % scalars % array(state%index_qs,k,i)
          qg_p(i,k,j) = state % scalars % array(state%index_qg,k,i)

          !number concentrations:
!         qnr_p(i,k,j) = max(0., state % scalars % array(state%index_qnr,k,i))
!         qni_p(i,k,j) = max(0., state % scalars % array(state%index_qni,k,i))
          qnr_p(i,k,j) = state % scalars % array(state%index_qnr,k,i)
          qni_p(i,k,j) = state % scalars % array(state%index_qni,k,i)
       enddo
       enddo
       enddo

    case default

 end select microp_select_init

 write(0,*) '--- end subroutine microphysics_from_MPAS:'

!formats:
 201 format(3i6,10(1x,e15.8))

 end subroutine microphysics_from_MPAS

!=============================================================================================
 subroutine microphysics_to_MPAS(state,tend)
!=============================================================================================

!output variables:
 type(state_type),intent(inout):: state
 type(tend_type) ,intent(inout):: tend

!---------------------------------------------------------------------------------------------

 write(0,*) '--- enter subroutine microphysics_to_MPAS:'

!initialization:
 tend % h_diabatic % array(:,:) = state % theta % array(:,:)

!copy all arrays back to the MPAS grid:
 do j = jts, jte
 do k = kts, kte
 do i = its, ite
    t_p(i,k,j) = th_p(i,k,j) * pi_p(i,k,j)
    state % theta   % array(k,i) = th_p(i,k,j)
    state % scalars % array(state%index_qv,k,i)  = qv_p(i,k,j)
    state % scalars % array(state%index_qc,k,i)  = qc_p(i,k,j)
    state % scalars % array(state%index_qr,k,i)  = qr_p(i,k,j)

    tend % h_diabatic % array(k,i) = (state % theta % array(k,i) &
                              -  tend % h_diabatic % array(k,i)) / dt_dyn
 enddo
 enddo
 enddo

!copy all arrays back to the MPAS grid as functions of cloud microphysics scheme:
 microp_select_init: select case(microp_scheme)

    case ("thompson")

       do j = jts, jte
       do k = kts, kte
       do i = its, ite
          !mass mixing ratios:
          state % scalars % array(state%index_qi,k,i)  = qi_p(i,k,j)
          state % scalars % array(state%index_qs,k,i)  = qs_p(i,k,j)
          state % scalars % array(state%index_qg,k,i)  = qg_p(i,k,j)

          !number concentrations:
          state % scalars % array(state%index_qnr,k,i) = qnr_p(i,k,j)
          state % scalars % array(state%index_qni,k,i) = qni_p(i,k,j)
       enddo
       enddo
       enddo

    case default

 end select microp_select_init

 write(0,*) '--- end subroutine microphysics_to_MPAS:'

 end subroutine microphysics_to_MPAS

!=============================================================================================
 end module mpas_atmphys_interface_hyd
!=============================================================================================


