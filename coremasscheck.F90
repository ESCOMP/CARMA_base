! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!!  This routine checks if the coremass exceeds the total
!!
!! @author Yunqian Zhu
!! @version Apr-2021
subroutine coremasscheck(carma, cstate, iz, fixcoremass,logmsg,abort, packagename, rc)

  ! types
  use carma_precision_mod
  use carma_enums_mod
  use carma_constants_mod
  use carma_types_mod
  use carmastate_mod
  use carma_mod

  implicit none

  type(carma_type), intent(inout)      :: carma       !! the carma object
  type(carmastate_type), intent(inout) :: cstate      !! the carma state object
  integer, intent(in)                  :: iz
  logical, intent(in)                  :: fixcoremass !! if we fix the coremass if exceeds total
  logical, intent(in)                  :: logmsg      !! write message to log file if core mass check fails
  logical, intent(in)                  :: abort       !! set return code to ERROR and return if  if core mass check fails
  character(len=*), intent(in)         :: packagename !! string to add to error message
  integer, intent(inout)               :: rc          !! return code, negative

  ! Local declarations
  integer                        :: igroup,ibin
  integer                        :: iepart,i,icore
  real(kind=f)                   :: coremass
  real(kind=f)                   :: factor
  real(kind=f),parameter         :: roundoff = 1.e-14_f

  ! check the coremass exceeding the total mass
  do igroup = 1,NGROUP
     iepart = ienconc(igroup)
     do ibin = 1,NBIN
        if (pc(iz, ibin, iepart) > 0._f) then
           coremass = 0._f
           do i = 1, ncore(igroup)
              icore = icorelem(i,igroup)
              coremass = coremass + pc(iz, ibin, icore)
           end do ! i = 1, ncore(igroup)
           if (coremass > pc(iz, ibin, iepart) * rmass(ibin, igroup)) then
             if(coremass < pc(iz, ibin, iepart) * rmass(ibin, igroup)*(1._f+roundoff))then
                pc(iz, ibin, iepart) = coremass/rmass(ibin,igroup) + ALMOST_ZERO
             else
               if(fixcoremass)then
                 pc(iz, ibin, iepart) = coremass/rmass(ibin,igroup) + ALMOST_ZERO
               endif
               if (logmsg) then
                 write(LUNOPRT,*) "Error - coremass exceeds total: ",packagename
                 write(LUNOPRT,*) "coremass",coremass,"total",pc(iz, ibin, iepart) * rmass(ibin,igroup)
               end if
               if (abort) then
                 rc = RC_ERROR
                 return
               end if
             end if
           end if
        end if
     end do
  end do

end subroutine coremasscheck
