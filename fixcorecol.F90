! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!! This routine checks the column for errors where the sum of the core masses
!! is larger than the concentration element. This implies a negative mass
!! of the concentration element, which is not physical.
!!
!! This routine attempts to conserve mass by using mass from the positive
!! values to offset the negative values.
!!
!! @author Charles Bardeen
!! version Feb-2023
subroutine fixcorecol(carma, cstate, rc)

  ! types
  use carma_precision_mod
  use carma_enums_mod
  use carma_constants_mod
  use carma_types_mod
  use carmastate_mod
  use carma_mod

  implicit none

  type(carma_type), intent(in)         :: carma   !! the carma object
  type(carmastate_type), intent(inout) :: cstate  !! the carma state object
  integer, intent(inout)               :: rc      !! return code, negative indicates failure

  ! Declare local variables
  integer        :: ibin
  integer        :: igroup
  integer        :: iepart
  integer        :: iz
  integer        :: icore
  integer        :: i
  real(kind=f)   :: total_core(NZ)
  real(kind=f)   :: concgas_md(NZ)
  real(kind=f)   :: total_mass
  real(kind=f)   :: missing_mass
  real(kind=f)   :: factor

  rc = RC_OK

  ! Advection can cause error in tracer/tracer relationship that can cause
  ! negative values for the concentration element. Use a mass conserving
  ! fixer to make sure there are no negative values.

  ! Find out the total amount of the concentration mass and the amount needed to
  ! fill in any negative values.
  !
  do igroup =  1, NGROUP

    ! We only need to do this where there are core masses.
    if (ncore(igroup) .gt. 0) then

      iepart = ienconc(igroup)

      do ibin = 1, NBIN

        ! Find the total of the core masses
        total_core(:) = 0._f

        do i = 1, ncore(igroup)
          icore = icorelem(i,igroup)
          total_core(:) = total_core(:) + pc(:, ibin, icore)
        end do

        ! Are there any places where the sum of the core masses is
        ! greater than the concentration element?
        if (any(total_core(:) > pc(:, ibin, iepart)*rmass(ibin, igroup))) then

          ! Determine the total mass and the missing mass.
          total_mass = 0._f
          missing_mass = 0._f

          do iz = 1, NZ
            concgas_md(iz) = pc(iz, ibin, iepart)*rmass(ibin, igroup) - total_core(iz)
            if (concgas_md(iz) > 0._f) then
              total_mass = total_mass + concgas_md(iz) * dz(iz)
            else if (concgas_md(iz) < 0._f) then
              missing_mass = missing_mass - concgas_md(iz) * dz(iz)
            end if
          end do

          ! Is there enough mass to fill in the holes?
          if (total_mass >= missing_mass) then

            ! Do we need to add any fudge factor to this for
            ! roundoff, ...?
            factor = (total_mass - missing_mass) / total_mass

            do iz = 1, NZ

              ! Scale the positive concentration mass by enough to fill
              ! the gaps.
              if (concgas_md(iz) > 0._f) then

                pc(iz, ibin, iepart) = (total_core(iz) + factor * concgas_md(iz)) / rmass(ibin, igroup)

              ! Fill the negative values with zero.
              else if (concgas_md(iz) < 0._f) then

                pc(iz, ibin, iepart) = total_core(iz) / rmass(ibin, igroup)
              end if
            end do
          else

            ! Since there isn't enough mass in the column to fix the negative values, just
            ! zero out the column.
            pc(:, ibin, iepart) = total_core(:) / rmass(ibin, igroup)
          end if
        end if
      end do
    end if
  end do

  return
end subroutine fixcorecol
