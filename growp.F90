! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!! This routine calculates particle source terms due to growth <growpe>
!! for one particle size bin at one spatial grid point per call.
!!
!! @author Andy Ackerman
!! @version Dec-1995
subroutine growp(carma, cstate, iz, ibin, ielem, rc)

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
  integer, intent(in)                  :: ibin    !! bin index
  integer, intent(in)                  :: ielem   !! element index
  integer, intent(inout)               :: rc      !! return code, negative indicates failure

  ! Local declarations
  integer                              :: igroup  ! group index
  integer                              :: iepart
  integer                              :: igas,icore,isub,growindex
  real(kind=f)                         :: add,dmdtetotal,growcoremass,dmdtratio

  ! Define group & particle # concentration indices for current element
  igroup = igelem(ielem)      ! target particle group 
  iepart = ienconc(igroup)	  ! target particle number concentration element

  ! Calculate production terms due to condensational growth <growpe>
  ! only if group to which element belongs grows.
  growindex = igrowgas(iepart)

  if (growindex .eq. 0)then
    growindex = igrowgas(iepart+ncore(igroup))
  end if

  if( growindex .ne. 0 .and. ibin .ne. 1 )then

    ! Bypass calculation if few droplets are present 
    if( pconmax(iz,igroup) .gt. FEW_PC )then
      growpe(ibin,ielem) = pc(iz,ibin-1,ielem) * growlg(ibin-1,igroup) 

      ! zhu 2015 09 09
      growcoremass = 0.
      dmdtetotal = 0.
      do isub = iepart,iepart+nelemg(igroup)-1
        if (itype(isub) .eq. I_GROWCORE)then
          growcoremass = growcoremass + pc(iz,ibin-1,isub)
          if(dmdte(ibin-1,isub-iepart,igroup) .gt. 0._f)&
            dmdtetotal = dmdtetotal + dmdte(ibin-1,isub-iepart,igroup)
        end if
      end do

      if (itype(ielem) .eq. I_GROWCORE) then
        igas = igrowgas(ielem)
        icore = ielem-iepart
        add = 0.0_f

        ! <add> is g/cm3 of growcore gas which is removed from vapor phase
        ! (and added to particle in ibin) due to growth so that the fraction
        ! will stay constant as the particle moves across bins

        if (dmdtmain(ibin-1,igroup) .ge. 0.0_f)then
          if (dmdte(ibin-1,icore,igroup) .gt. 0.0_f)then
            ! Maingas and growcore gas should both grow on particle
            ! (for now will grow at constant fraction)
            dmdtratio = dmdte(ibin-1,icore,igroup)/(dmdtetotal+dmdtmain(ibin-1,igroup))
            add = pc(iz,ibin-1,iepart) * growlg(ibin-1,igroup) * &
            diffmass(ibin,igroup,ibin-1,igroup) * dmdtratio
          else
            ! Maingas should grow, but growcore gas should evaporate
            add = 0.0_f
          endif
        else
          ! Grow the growcore, but keep the primary cloud volatile
          ! constant
          ! (might need to change this part for multiple growcores, ie.
          !  add condition if current growcore should not grow)
          if(dmdtetotal .ne. 0)then
            dmdtratio = dmdte(ibin-1,icore,igroup)/dmdtetotal
          else
            dmdtratio = 0
          end if
          add = dmdtratio*pc(iz,ibin-1,iepart)*growlg(ibin-1,igroup)* &
          diffmass(ibin,igroup,ibin-1,igroup)
        endif
       
        growpe(ibin,ielem) = growpe(ibin,ielem) + add
        gprod_grow(igroup,igas) = gprod_grow(igroup,igas) - add

      endif !GROWCORE
    endif
  endif

  ! Return to caller with growth production terms evaluated.
  return
end
