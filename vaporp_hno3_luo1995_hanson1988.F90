! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!!  Calculates the vapor pressure for Nitric acid for STS.
!!  <pvap_liq> and <pvap_ice> are vapor pressures in units of [dyne/cm^2]
!!  This subroutine is used when NH4 is not considered in the model
!!  If NH4 is present in the particle, should use Lin&Tabazadeh 2001 instead
!!
!! @author Yunqian
!! @version Oct-2011

subroutine vaporp_hno3_luo1995_hanson1988(carma, cstate, iz, rc, pvap_liq, pvap_ice)
!     types
  use carma_precision_mod
  use carma_enums_mod
  use carma_constants_mod
  use carma_types_mod
  use carmastate_mod
  use carma_mod
  use sulfate_utils
  
  implicit none

  type(carma_type), intent(inout)      :: carma     !! the carma object
  type(carmastate_type), intent(inout) :: cstate    !! the carma state object
  integer, intent(in)                  :: iz        !! z index
  real(kind=f), intent(out)            :: pvap_liq  !! vapor pressure wrt liquid [dyne/cm2]
  real(kind=f), intent(out)            :: pvap_ice  !! vapor pressure wrt ice [dyne[cm2]
  integer, intent(inout)               :: rc        !! return code, negative indicates failure

  ! Local declarations
  real(kind=f)            :: A_hno3,B_hno3
  real(kind=f)            :: wt1,wt2
  real(kind=f)            :: mt, bt
  real(kind=f)		  :: rvap
  real(kind=f) 		  :: gc_cgs
  real(kind=f), parameter :: a0    = 22.74_f
  real(kind=f), parameter :: b0    = -7689.8_f
  real(kind=f), parameter :: a1    = 29.0_f
  real(kind=f), parameter :: b1    = -2896.1_f
  real(kind=f), parameter :: a2    = -0.457_f
  real(kind=f), parameter :: b2    = 2859.8_f
  real(kind=f), parameter :: a3    = -5.03_f
  real(kind=f), parameter :: b3    = -274.2_f
  real(kind=f), parameter :: a4    = 1.447_f
  real(kind=f), parameter :: b4    = -389.5_f
  real(kind=f), parameter :: a5    = -29.72_f
  real(kind=f), parameter :: b5    = 7281.0_f
  real(kind=f), parameter :: a6    = -13.9_f
  real(kind=f), parameter :: b6    = 6475.0_f
  real(kind=f), parameter :: a7    = 6.1_f
  real(kind=f), parameter :: b7    = 801.0_f

  ! wtpct is NH4+H2SO4, but since NH4 is not present, wtpct is H2SO4
  call wtpct_luo_total(carma, cstate, iz, rc, wtpct(iz),wtpct_hno3(iz))
  if (rc < 0) return
     
  wt1 = max(wtpct_hno3(iz)/100._f,0.001_f)
  wt2 = wtpct(iz)/100._f

  A_hno3 = a0 + a1*wt1 + a2*wt2 + a3*sqrt(wt1) + a4*sqrt(wt2) &
              + a5*wt1*wt1 + a6*wt1*wt2 + a7*wt2*wt2
  B_hno3 = b0 + b1*wt1 + b2*wt2 + b3*sqrt(wt1) + b4*sqrt(wt2) &
              + b5*wt1*wt1 + b6*wt1*wt2 + b7*wt2*wt2
	      
  mt = -2.7836_f - 0.00088_f * t(iz)
  bt = 38.9855_f - 11397.0_f/t(iz) +0.009179_f * t(iz)

  !Vapor pressure   
  pvap_liq = exp(A_hno3 + B_hno3/t(iz) + log(wt1*(wt1+0.09_f*wt2))) ! mbar
    
  !calculate the ambience water vapor pressure
  gc_cgs = gc(iz,1) / zmet(iz)
  rvap = RGAS / gwtmol(1)
       
  pvap_ice = 10._f**(mt * log10(gc_cgs * rvap * t(iz)*0.00075_f) + bt)  ! torr
       
  ! Convert ==> dynes/cm2
  pvap_liq = pvap_liq * 1000._f
  pvap_ice = pvap_ice * 1333.2237_f   !torr to dynes/cm2

  return
end
