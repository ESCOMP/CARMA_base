! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!! This routine calculates aerosol hygroscopicity
!!
!! @author  Pengfei Yu, Mike Mills
!! @version Oct 2020
subroutine hygroscopicity(carma, cstate, rc)
  ! types
  use carma_precision_mod
  use carma_enums_mod
  use carma_constants_mod
  use carma_types_mod
  use carmastate_mod
  use carma_mod
  use carmaelement_mod

  implicit none

  type(carma_type), intent(in)         :: carma   !! the carma object
  type(carmastate_type), intent(inout) :: cstate  !! the carma state object
  integer, intent(inout) :: rc      !! return code, negative indicates failure

  ! Local declarations
  integer         :: igroup         !! group index
  integer         :: ibin           !! bin index
  integer         :: iepart         !! element in group containing the particle concentration
  integer         :: icore, i, z
  real(kind=f)    :: coremass(NZ), shellmass(NZ)
  
  1 format('hygroscopicity::ibin=',i4,',core mass=',e10.3,',shell mass=',e10.3,',hygroscopicity=',f9.4)
              
  hygro(:NZ,:NBIN,:NGROUP) = 0.5_f ! default
  
  ! loop through all bins, groups, and elements to calculate hygro for each group:  
  do igroup = 1,NGROUP
    ! Only calculate hygro for groups that use it
    if (irhswell(igroup) == I_KAPPA) then 
      iepart = ienconc(igroup)     ! element of particle number concentration             
      do ibin = 1, NBIN
        do z = 1, NZ
          if (pc(z, ibin, iepart).gt.0._f) then
            hygro(z,ibin,igroup) = 0._f

            ! Weight hygro by mass of each core
            coremass = 0._f      
            do i = 1, ncore(igroup)
              icore = icorelem(i, igroup)
              coremass(z) = coremass(z) + pc(z, ibin, icore)

              hygro(z,ibin,igroup) = hygro(z,ibin,igroup) + pc(z,ibin,icore) * kappaelem(icore)
            end do ! i = 1, ncore(igroup)

            ! Add shell mass to hygro weighting      
            shellmass(z) = max((pc(z, ibin, iepart) * rmass(ibin, igroup)) - coremass(z), 0._f)
            hygro(z,ibin,igroup) = hygro(z,ibin,igroup) + shellmass(z) * kappaelem(iepart)

            !Divide by total mass of all particles in the bin to normalize:
            hygro(z,ibin,igroup) = hygro(z,ibin,igroup) / pc(z, ibin, iepart) / rmass(ibin, igroup)
          end if
!          write(LUNOPRT,1) ibin,coremass(z),shellmass(z),hygro(z,ibin,igroup)
        end do ! z = 1, NZ
      end do ! ibin = 1, NBIN
    end if ! irhswell(igroup) == I_KAPPA
  end do ! igroup = 1,NGROUP

  rc = RC_OK
  
  return 
end subroutine !hygroscopicity
