! Include shortname defintions, so that the pre-existing F77 code does not have to be modified to
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
  real(kind=f)    :: coremass, shellmass
  real(kind=f), parameter :: thresh  = 1e-14_f     ! relative threshold for core mass and kappa checking
  real(kind=f), parameter :: thresh1 = 1._f+thresh
  real(kind=f), parameter :: thresh0 = 0._f-thresh

  1 format('hygroscopicity::ibin=',i4,',core mass=',e10.3,',shell mass=',e10.3,',hygroscopicity=',f9.4)

  kappahygro(:NZ,:NBIN,:NGROUP) = -huge(1._f) ! default

  ! loop through all bins, groups, and elements to calculate hygro for each group:
  do igroup = 1,NGROUP
    ! Only calculate hygro for groups that use it
    if (irhswell(igroup) == I_PETTERS) then
      iepart = ienconc(igroup)     ! element of particle number concentration
      do ibin = 1, NBIN
        do z = 1, NZ
          kappahygro(z,ibin,igroup) = 0._f
          if (pc(z, ibin, iepart).gt.0._f) then
            ! Weight hygro by mass of each core
            coremass = 0._f
            do i = 1, ncore(igroup)
              icore = icorelem(i, igroup)
              coremass = coremass + pc(z, ibin, icore)
              kappahygro(z,ibin,igroup) = kappahygro(z,ibin,igroup) + pc(z,ibin,icore) * kappaelem(icore)
            end do ! i = 1, ncore(igroup)

            ! Add shell mass to hygro weighting
            shellmass = max((pc(z, ibin, iepart) * rmass(ibin, igroup)) - coremass, 0._f)
            kappahygro(z,ibin,igroup) = kappahygro(z,ibin,igroup) + shellmass * kappaelem(iepart)
            !Divide by total mass of all particles in the bin to normalize:
            kappahygro(z,ibin,igroup) = kappahygro(z,ibin,igroup) / pc(z, ibin, iepart) / rmass(ibin, igroup)
          end if
          if (kappahygro(z,ibin,igroup).gt.thresh1.or.kappahygro(z,ibin,igroup).lt.thresh0) then
            write(LUNOPRT,*) "hygro77: z,ibin,kappahygro,pc,rmass,shellmass,coremass", &
              z,ibin,kappahygro(z,ibin,igroup),pc(z, ibin, iepart),rmass(ibin, igroup),shellmass,coremass
            rc=RC_ERROR
            return
          end if
          kappahygro(z,ibin,igroup) = min(kappahygro(z,ibin,igroup),1._f)
          kappahygro(z,ibin,igroup) = max(kappahygro(z,ibin,igroup),0._f)
        end do ! z = 1, NZ
      end do ! ibin = 1, NBIN
    end if ! irhswell(igroup) == I_PETTERS
  end do ! igroup = 1,NGROUP

  rc = RC_OK

  return
end subroutine !hygroscopicity
