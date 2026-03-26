! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!! This routine calculates the total amount of condensate associated with each gas.
!!
!! @author Chuck Bardeen
!! @version Nov-2009
subroutine totalcondensate(carma, cstate, iz, total_ice, total_liquid, rc)

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
  integer, intent(in)                  :: iz      !! z index
  real(kind=f), intent(out)            :: total_ice(NGAS)      !! total ice at the start of substep
  real(kind=f), intent(out)            :: total_liquid(NGAS)   !! total liquid at the start of substep
  integer, intent(inout)               :: rc      !! return code, negative indicates failure

  ! Local declarations
  integer                              :: igroup  ! group index
  integer                              :: icore   ! core index
  integer                              :: igas    ! gas index
  integer                              :: ibin    ! bin index
  integer                              :: ielem   ! element index
  integer                              :: i
  real(kind=f)                         :: coremass
  real(kind=f)                         :: volatilemass


  ! Initialize local variables for keeping track of gas changes due
  ! to nucleation and growth in each particle group.
  total_ice(:)    = 0._f
  total_liquid(:) = 0._f

  ! Iterate over each particle type and total up that ones that interact
  ! with the gases.
  !
  ! This code assumes that all changes in condensate are associated with
  ! growth in a particular gas. This doesn't handle all possible changes
  ! associated with nucleation, if the group do not also participate in
  ! growth.
  do igroup = 1,NGROUP

    ielem = ienconc(igroup)     ! element of particle number concentration

    igas = igrowgas(ielem)      ! condensing gas

    if ((itype(ielem) == I_VOLATILE) .and. (igas /= 0)) then

      do ibin = 1, NBIN
   
        ! If this group has core masses, then determine the involatile component.
        coremass = 0._f
      
        do i = 1, ncore(igroup)
          icore = icorelem(i, igroup)
          coremass = coremass + pc(iz, ibin, icore)
        end do
        
        volatilemass = (pc(iz, ibin, ielem) * rmass(ibin, igroup)) - coremass
        
        ! There seem to be times when the coremass becomes larger than the total
        ! mass. This shouldn't happen, but check for it here.
        !
        ! NOTE: This can be caused by advection in the parent model or sedimentation
        ! in this model.
        if (volatilemass > 0._f) then
          if (is_grp_ice(igroup)) then
            total_ice(igas) = total_ice(igas) + volatilemass
          else
            total_liquid(igas) = total_liquid(igas) + volatilemass
          end if
        else
          if(abs(volatilemass) .gt. FEW_PC)then
            if(do_print) write(*,*) 'totalcondensate.F90:ERROR - VOLATILE mass is negative.'
            if(do_print) write(*,*) 'ibin,',ibin,'iz',iz
            if(do_print) write(*,*) 'H2SO4',pc(iz, ibin, ielem),'rmass',rmass(ibin, igroup)
            if(do_print) write(*,*) 'H2SO4 mass',pc(iz, ibin, ielem) * rmass(ibin, igroup), 'coremass',coremass
            if(do_print) write(*,*) 'FEW_PC',FEW_PC,'SMALL_PC',SMALL_PC
            if(do_print) write(*,*) "ncore(igroup)",ncore(igroup)
            do i = 1, ncore(igroup)
              if(do_print) write(*,*) "icorelem(i, igroup)",icorelem(i, igroup)
              if(do_print) write(*,*) "pc(iz, ibin, icore)", pc(iz, ibin, icore)
            enddo
            rc = RC_WARNING_RETRY
            return
          else
            total_liquid(igas) = total_liquid(igas)+0._f
            total_ice(igas) = total_ice(igas)+0._f
          end if
        end if
      end do
    end if

    do i = 1, ncore(igroup)
      icore = icorelem(i,igroup)
      igas = igrowgas(icore)

      if((itype(icore) == I_GROWCORE) .and. (igas /= 0)) then
        do ibin = 1, NBIN
          if (pc(iz,ibin,icore) > 0._f) then
            if(is_grp_ice(igroup)) then
              total_ice(igas) = total_ice(igas) + pc(iz,ibin,icore)
            else
              total_liquid(igas) = total_liquid(igas) + pc(iz,ibin,icore)
            end if
            ! else
            ! rc = -1
            ! if(do_print) write(*,*) 'totalcondensate.F90:ERROR - GROWCORE mass is negative.'
          end if
        end do
      end if
    end do

  end do
      
  return
end
