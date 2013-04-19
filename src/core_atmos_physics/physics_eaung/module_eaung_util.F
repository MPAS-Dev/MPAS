!==============================================================================
 MODULE module_eaung_util

 USE kinds
 USE module_eaung_params

 IMPLICIT NONE
 PRIVATE
 PUBLIC :: eau_sat

!==============================================================================

 CONTAINS

!==============================================================================
 SUBROUTINE eau_sat &
                (   len ,  pl ,    tl ,  qwvl ,  qsats , &
                 dqsats , ess , dtess ,  sats ,   wght   &
                )
!------------------------------------------------------------------------------
!eau_sat computes the saturation mixing ratio, vapor pressure, and saturation,
!and their derivatives with respect to temperature over water, ice and mixed-
!phase clouds. the parameterization is that used in the Community Climate Com-
!munity Climate Model at NCAR.
!Laura D. Fowler /slikrock (07-01-01).

!send comments to laura@atmos.colostate.edu.

!subroutines called:
!none.

!ARGUMENT LIST VARIABLES:
!INPUT ARGUMENTS:
!----------------
 integer (kind=int_kind), intent(in):: &
    len               !length of vector.

 real(kind=dbl_kind), intent(in), dimension(len):: &
    pl,              &!pressure                                           (Pa).
    tl                !temperature                                         (K).

 real(kind=dbl_kind), intent(in), dimension(len), optional :: &
    qwvl              !water vapor mixing ratio                        (kg/kg).

!OUTPUT ARGUMENTS:
!-----------------
 real (kind=dbl_kind), intent(out), dimension(len), optional :: &
    qsats,           &!saturation mixing ratio                         (kg/kg).
    dqsats,          &!derivative of qsat with respect to temperature(kg/kg/K).
    ess,             &!saturation vapor pressure                          (Pa).
    dtess,           &!derivative of es with respect to temperature     (Pa/K).
    sats,            &!saturation                                          (-).
    wght              !partitioning between cloud water and cloud ice
                      !for mixed-phase clouds                              (-).

!LOCAL VARIABLES:
INTEGER (KIND=INT_KIND):: i

REAL (KIND=DBL_KIND):: &
   twmin=173.16,     &!lowest allowed temperature boundary for water       (K).
   twmax=373.16,     &!highest allowed temperature boundary for water      (K).
   timin=173.16,     &!lowest allowed temperature boundary for ice         (K).
   timax=273.16,     &!highest allowed temperature boundary for ice        (K).
   tnull=273.16       !freezing temperature                                (K).

REAL (KIND=DBL_KIND) :: tstl , t0tl

REAL (KIND=DBL_KIND), DIMENSION(len):: &
      esw ,      esi ,      esm , &
      tl0 ,    wghtm , qsat_tem , &
   es_tem , dtes_tem , wght_tem

!CCM PARAMETERIZATION OF SATURATION VAPOR PRESSURE OVER WATER AND ICE:
REAL (KIND=DBL_KIND), PARAMETER:: &
   ps = 1013.246,    &!reference pressure                                (hPa).
   ts = 373.16,      &!reference temperature                               (K).
   t0 = 273.16        !freezing temperature                                (K)

REAL (KIND=DBL_KIND):: &
       e1 ,   e2 ,     f ,    f1 , &
       f2 ,   f3 ,    f4 ,    f5 , &
   lphase , term , term1 , term2 , &
   term3

REAL (KIND=DBL_KIND):: tem

!------------------------------------------------------------------------------
!write(6,*)
!write(6,*) '--- enter subroutine EAU_SAT: len=', len
 do i = 1, len
    esi(i) = eauc0
    esm(i) = eauc0
    esw(i) = eauc0
    es_tem(i) = eauc0
 enddo

!SATURATION OVER WATER:
 do i = 1, len
    tl0(i)    = max(twmin,tl(i))
    tl0(i)    = min(twmax,tl0(i))
    if(tl0(i).ge.tbgmin) then
       tstl   = ts / tl0(i)
       e1     = 11.344_dbl_kind*(1.0_dbl_kind - tl0(i)/ts)
       e2     = -3.49149_dbl_kind*(tstl - 1.0_dbl_kind)
       f1     = -7.90298_dbl_kind*(tstl - 1.0_dbl_kind)
       f2     = 5.02808_dbl_kind*log10(tstl)
       f3     = -1.3816_dbl_kind*(10.0_dbl_kind**e1-1.0_dbl_kind) * 0.0000001_dbl_kind
       f4     = 8.1328_dbl_kind*(10.0_dbl_kind**e2-1.0_dbl_kind) * 0.001_dbl_kind
       f5     = log10(ps)
       f      = f1 + f2 + f3 + f4 + f5

       esw(i) = (10.0_dbl_kind**f)*1.e+02_dbl_kind
       tem    = esw(i)
       esw(i) = min(esw(i),pl(i)*0.9_dbl_kind)
       wght_tem(i)  = eauc1
       es_tem(i)    = esw(i)
    endif
 enddo

!SATURATION OVER ICE:
 do i = 1, len
    if(tl0(i).lt.timax) then
       tl0(i)    = max(tl0(i),timin)
       t0tl      = t0 / tl0(i)
       term1     = 2.01889049_dbl_kind/(t0tl)
       term2     = 3.56654_dbl_kind*log(t0tl)
       term3     = 20.947031_dbl_kind*(t0tl)
       term      = term1 + term2 + term3

       esi(i)    = 575.185606e10_dbl_kind*exp(-(term1 + term2 + term3))
       esi(i)    = min(esi(i),pl(i)*0.9_dbl_kind)
       wght_tem(i)  = eauc0
       es_tem(i)    = esi(i)
    endif
 enddo

!interpolated saturation variables:

 do i = 1, len
    if(tl0(i).lt.tbgmax .and. tl0(i).ge.tbgmin) then
       wghtm(i)  = (tl0(i)-tbgmin)/(tbgmax-tbgmin)
       esm(i)    = wghtm(i)*esw(i) + (1.-wghtm(i))*esi(i)
       esm(i)    = min(esm(i),pl(i)*0.9)
       es_tem(i)    = esm(i)
       wght_tem(i)  = wghtm(i)
    endif
 enddo

 if(present(ess)) ess = es_tem

 if(present(dtess).or.present(dqsats)) then
    do i = 1, len
       lphase    = lcond*wght_tem(i)+lsub*(1.-wght_tem(i))
       dtes_tem(i) = lphase*es_tem(i)/(rgas*tl0(i)*tl0(i))
    enddo
    if(present(dtess)) dtess = dtes_tem
 endif

 if(present(qsats) .or. present(dqsats) .or. present(sats))  then
    qsat_tem = 0.622_dbl_kind * es_tem / (pl - es_tem)
    if(present(qsats) ) qsats = qsat_tem
    if(present(dqsats)) dqsats = qsat_tem*dtes_tem*pl / (es_tem* &
                 (pl-es_tem))
    if(present(sats)) then
       if(present(qwvl)) then
          sats = qwvl / qsat_tem
       else
          sats = eauc0
       endif
    endif
 endif

 if(present(wght)) wght = wght_tem

 END SUBROUTINE eau_sat

!==============================================================================

 END MODULE module_eaung_util

!==============================================================================
