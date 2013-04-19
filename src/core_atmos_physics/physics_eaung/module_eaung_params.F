!==============================================================================
!WRF:MEDIATION LAYER:EAU_NG
 MODULE module_eaung_params
 USE kinds

 IMPLICIT NONE
 PUBLIC
 SAVE

!module_eaung_params is a modified version of the module eau_params.F developed
!for the CSU GCM. It contains all the constants needed for the eau_cup and eau-
!liq parameterizations of convection and cloud microphysics.

!Laura D. Fowler (CIRA/Colorado State University).
!fowler@cira.colostate.edu.

!REFERENCES:
!Fowler, L.D, D.A. Randall, and S.A. Rutledge, 1996: Liquid and Ice Cloud Micro
!physics in the CSU General Circulation Model: Model description and simulated
!cloud microphysical processes.
!J. Climate, 9, 489-529.

!Randall, D.A., and L.D. Fowler, 1999: EAUliq: The next generation. Dept. of
!Atmospheric Science Paper 673, Dept. of Atmospheric Science. Colorado State
!University, Fort Collins, Colorado, 65 pp.

!==============================================================================

!EAU_CUP PARAMETERS:
INTEGER (KIND=INT_KIND), PARAMETER:: &
   ncup   = 1,                      &!
   loadw  = 0                        !set the water loading in the computation
                                     !of the cloud work function and virtual
                                     !temperature.
!REAL (KIND=DBL_KIND), PARAMETER:: &
DOUBLE PRECISION, PARAMETER:: &
   alpham = 1.e+08_dbl_kind,        &!
   ckemin = 5.0_dbl_kind,           &!minimum value of cke.
   taudis = 600.0_dbl_kind,         &!cke dissipation time scale.
   amiu   = 0.0_dbl_kind             !cloud-top entrainment.

!EAU_LIQ PARAMETERS:
INTEGER (KIND=INT_KIND), PARAMETER:: &
!  nc1elq = 4                        !number of microphysics time-steps per
!                                    !dynamical time-step.
   nc1elq = 1                        !number of microphysics time-steps per
                                     !dynamical time-step.
!REAL (KIND=DBL_KIND), PARAMETER:: &
DOUBLE PRECISION, PARAMETER:: &
   a0elq    = -.267_dbl_kind,       &!
   a1elq    = 5.15e+03_dbl_kind,    &!
   a2elq    = -1.0225e+06_dbl_kind, &!
   a3elq    = 7.55e+07_dbl_kind,    &!
   alphaelq = 0.001_dbl_kind,       &!
   aprime   = 3.e+03_dbl_kind,      &!
   asecond  = 1.139_dbl_kind,       &!
   belq     = 0.11_dbl_kind,        &!
   betaelq  = 0.001_dbl_kind,       &!
   diffelq  = 2.26e-05_dbl_kind,    &!
   erccoef  = 1._dbl_kind,          &!
   esccoef  = 1._dbl_kind,          &!
   esicoef  = 0.1_dbl_kind,         &!
   gam3     = 2._dbl_kind,          &!
   gams1    = 2.21891_dbl_kind,     &!
   gams2    = 1.38784_dbl_kind,     &!
   gams3    = 6.90080_dbl_kind,     &!
   kap      = .2861328125_dbl_kind, &!
   muelq    = 1.718e-05_dbl_kind,   &!
   nzeror   = 8.e+06_dbl_kind,      &!
   nzeros   = 2.e+07_dbl_kind,      &!
   pielq    = 3.14159265,           &!
   pzero    = 1.e+05_dbl_kind,      &!
!  qci0     = 0.01e-03_dbl_kind,    &!
   qci0     = 0.1e-03_dbl_kind,     &!
!  qcw0     = 0.25e-03_dbl_kind,    &!
   qcw0     = 0.7e-03_dbl_kind,     &!
   rhor     = 1.e03_dbl_kind,       &!
   rhos     = 1.e02_dbl_kind,       &!
   taul     = 5._dbl_kind,          &!
   tauf     = 5._dbl_kind,          &!
   therco   = 2.43e-02_dbl_kind      !

!SHARED EAU_LIQ AND EAU_CUP PARAMETERS:
LOGICAL (KIND=LOG_KIND), PARAMETER:: &
   w_diag   = .true.

!REAL (KIND=DBL_KIND), PARAMETER:: &
DOUBLE PRECISION, PARAMETER:: &
   eauc0    = 0.0_dbl_kind,         &!
   eauc1    = 1.0_dbl_kind,         &!
   ethresh  = 1.e-12_dbl_kind        !threshold below which there is no cloud
                                     !microphysics processes.

!REAL (KIND=DBL_KIND), PARAMETER:: &
DOUBLE PRECISION, PARAMETER:: &
   delta    = 0.622_dbl_kind,       &!
   gravelq  = 9.81_dbl_kind,        &!
   cpelq    = 1005.000_dbl_kind,    &!
   r_d      = 287.000_dbl_kind,     &!
   rgas     = 4.61e+02_dbl_kind,    &!
   rcp      = r_d/cpelq,            &!
   lcond    = 2.5e+06_dbl_kind,     &!
   lfus     = 0.3336e+06_dbl_kind,  &!
   lsub     = lcond+lfus,           &!
   t00      = 273.15_dbl_kind,      &!
!  tbgmin   = 253.15_dbl_kind,      &!
   tbgmin   = 233.15_dbl_kind,      &!
   tbgmax   = 273.15_dbl_kind        !

!EAU_NG PARAMETERS:
!REAL (KIND=DBL_KIND), PARAMETER:: &
DOUBLE PRECISION, PARAMETER:: &
   cClrCld  = eauc1,                &!exchange terms between cld and clr.
   cCldClr  = eauc1                  !exchange terms between cld and clr.

!==============================================================================

 END MODULE module_eaung_params

!==============================================================================
