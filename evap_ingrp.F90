! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!! This routine calculates particle source terms <evappe> of droplets
!! evaporating within a particle group.
!!
!! Distinct evaporation of cores has not been treated.
!!
!! @author Andy Ackerman
!! @version Aug-2001
subroutine evap_ingrp(carma,cstate,iz,ibin,ig,ip,rc)

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
  integer, intent(in)                  :: ig      !! group index
  integer, intent(in)                  :: ip      !! the number density element
  integer, intent(inout)               :: rc      !! return code, negative indicates failure

  ! Local declarations
  integer                              :: ie, ib, ibloop ,icore,nothercore
  integer                              :: isub
  real(kind=f)                         :: ndelemmass, dmdtratio, growcoremass,dmdtetotal,coremassmass
  real(kind=f)                         :: togas(nelemg(ig)), remove(nelemg(ig)), rest
  real(kind=f)                         :: restparticle,sumtogas
  integer                              :: growindex(nelemg(ig)),index_num,iindex,smallevapindex,smallndindex

  ! For a single group, the core mass fraction is 0.
  cmf(ibin,ig) = 0.0_f
  
  ! The smallest bin cannot be a source to smaller bins in same group
  if( ibin .eq. 1 )then
    return
  endif

  ! Evaluate evaporation source term <evappe> for all elements in group
  if (itype(ip+nelemg(ig)-1) .ne. I_GROWCORE) then
    do isub = 1, nelemg(ig)
      ie = ip + isub - 1
      evappe(ibin-1,ie) = evappe(ibin-1,ie) + pc(iz,ibin,ie)*evaplg(ibin,ig)
    enddo
  else
    ! the particle contains GROWCORE.

    ! first, calculate how many mass of particles have the type of GROWCORE and
    ! how much growcore want to evaporate.
    growcoremass = 0.       !GROWCORE particle concentration g/unit
    coremassmass = 0.       !COREMASS particle concentration g/unit
    dmdtetotal = 0.
    nothercore = 0
    do isub = ip, ip+nelemg(ig)-1
      if(itype(isub) .ne. I_GROWCORE) then
        nothercore = nothercore + 1
      end if
      if (itype(isub) .eq. I_GROWCORE)then
        growcoremass = growcoremass + pc(iz,ibin,isub)
        if(dmdte(ibin-1,isub-ip,ig) .lt. 0._f)dmdtetotal = dmdtetotal + dmdte(ibin-1,isub-ip,ig)
        !write(*,*) "evep_ingrp: dmdte",dmdte(ibin-1,isub-ip,ig) 
      end if
      if (itype(isub) .eq. I_COREMASS)then
        coremassmass = coremassmass + pc(iz,ibin,isub)
      end if
    end do

    ! calculate the mass of number density element in ONE particle.
    ndelemmass = rmass(ibin,ig) - (growcoremass+coremassmass)/pc(iz,ibin,ip)
    if (ndelemmass .lt. 0._f) then
      if (abs(ndelemmass) .gt. rmass(ibin,ig)*1.e-8_f .and. pc(iz,ibin,ip) .gt. FEW_PC)then
        write(*,*) "ndelemmass less than ZERO in group ",ig, " , adjust ",ndelemmass," to ZERO."
        write(*,*) "rmass(ibin,ig)",rmass(ibin,ig),"ibin",ibin,"ip",ip,"pc(iz,ibin,ip)",pc(iz,ibin,ip),"FEW_PC",FEW_PC
        write(*,*) "total mass in the bin",pc(iz,ibin,ip)*rmass(ibin,ig),"HNO3 in the bin",pc(iz,ibin,ip+1)
        write(*,*) "growcoremass",growcoremass
        rc = RC_ERROR
        return
      else
        ndelemmass = 0._f
      endif
    endif

    ! index_num: How many particles in a group can grow and evaporate.
    ! growindex: what is the element index that can grow and evaporate.
    growindex = 0
    index_num = 0

    do isub = 1, nelemg(ig)
      ie = ip + isub - 1
      if(itype(ie) .eq. I_VOLATILE)then
        index_num = index_num + 1
        growindex(index_num) = ie
      else
        if(itype(ie) .eq. I_GROWCORE)then
          index_num = index_num + 1
          growindex(index_num) = ie
        end if
      end if
    end do
    
    !  <remove> is the total g/cm3 of growcore gas that was lost from the
    !  particle in ibin due to evaporation.
    !  <togas> is the g/cm3 of growcore gas that will be liberated to
    !  the vapor phase (from the ibin particle).
    !  <remove> - <togas> is the g/cm3 of growcore gas which will be
    !  transfered with the particle to ibin-1
    smallevapindex = 0          !when togas > remove, add 1
    smallndindex = 0
    togas = 0._f
    remove = 0._f

    do iindex = 1, index_num
      if (itype(growindex(iindex)) .eq. I_VOLATILE)then
        remove(iindex) = pc(iz,ibin,growindex(iindex)) * evaplg(ibin,ig)*ndelemmass
        
        if (abs(dmdtmain(ibin-1,ig)) .lt. abs(dmdtetotal)*1.e-8_f)then
          dmdtmain(ibin,ig) = 0.0_f
        endif
    
        if (dmdtmain(ibin-1,ig) .le. 0.0_f) then
          dmdtratio = dmdtmain(ibin-1,ig)/(dmdtetotal+dmdtmain(ibin-1,ig))
          togas(iindex) = dmdtratio*pc(iz,ibin,ip)*evaplg(ibin,ig)*diffmass(ibin,ig,ibin-1,ig)
        else
          togas(iindex) = 0.0_f
        endif

        if(togas(iindex) .gt. remove(iindex))then
          smallndindex = smallndindex + 1
        end if

      else if (itype(growindex(iindex)) .eq. I_GROWCORE) then
        remove(iindex) = pc(iz,ibin,growindex(iindex)) * evaplg(ibin,ig)

        icore = growindex(iindex)-ip
        if (abs(dmdtmain(ibin-1,ig)) .lt. abs(dmdte(ibin-1,icore,ig))*1.e-8_f)then
          dmdtmain(ibin-1,ig) = 0.0_f
        endif
        !write(*,*) "ie",growindex(iindex),"icore-nothercore+1",growindex(iindex)-ip-nothercore+1

        if (dmdtmain(ibin-1,ig) .le. 0.0_f) then
          if (dmdte(ibin-1,icore,ig) .lt. 0.0_f) then
            ! Both maingas and growcore gas evaporate from particle
            ! (keep at constant fraction for now)
            dmdtratio = dmdte(ibin-1,icore,ig)/(dmdtetotal+dmdtmain(ibin-1,ig))
            togas(iindex) = dmdtratio*pc(iz,ibin,ip)*evaplg(ibin,ig)*diffmass(ibin,ig,ibin-1,ig)
          else
            ! Maingas evaporates, but not growcore
            togas(iindex) = 0.0_f
            dmdtratio = 0.
          endif
        else
          if (dmdte(ibin-1,icore,ig) .lt. 0.0_f ) then
            ! Growcore evaporates, but not primary cloud volatile
            ! (may need to change this for multiple growcores)
            if(dmdtetotal .ne. 0)then
              dmdtratio = dmdte(ibin-1,icore,ig)/dmdtetotal
            else
              dmdtratio = 0._f
            end if
            togas(iindex) = dmdtratio*pc(iz,ibin,ip)*diffmass(ibin,ig,ibin-1,ig)*evaplg(ibin,ig)
          else
            togas(iindex) = 0.0_f
          endif
        endif

        if(togas(iindex) .gt. remove(iindex))then
          smallevapindex = smallevapindex + 1
        end if

      endif
    enddo

!    if(ig .eq. 5)write(*,*) "togas",togas,"remove",remove

    ! if remove greater than the mass size of lower bin,
    ! HNO3 will exceed the total mass
    rest = 0._f

    sumtogas = 0
    do iindex = 1, index_num
      if(togas(iindex) .gt. remove(iindex))then
        togas(iindex) = remove(iindex)
        sumtogas = sumtogas +togas(iindex)
      else
        sumtogas = sumtogas +togas(iindex)
      endif
    enddo

    !write(*,*) "ig",ig,"ibin",ibin,"smallevapindex",smallevapindex

    ! Evaporating number density element is different from evaporating the GROWCORE
    ! becasue we cannot evaporate everything for number density element
    ! If the number density element is smallevapindex, we evaporate 90% of the mass
    ! rest split between ibin and ibin-1
    if(dmdtmain(ibin-1,ig) .lt. 0._f .and. smallndindex .gt. 0)then
      ! zhu 2015 09 18 if H2SO4 is a small amount and want to evaporate to lower bin
!      write(*,*) "here 1..."

      rest = rmass(ibin,ig)*pc(iz,ibin,ip)*evaplg(ibin,ig)-sumtogas+ndelemmass*pc(iz,ibin,ip)*evaplg(ibin,ig)*0.1_f
      
      if(rest .gt. rmass(ibin,ig)*pc(iz,ibin,ip)*evaplg(ibin,ig) .or. rest .lt. rmass(ibin-1,ig)*pc(iz,ibin,ip)*evaplg(ibin,ig))then
        if(rest .lt. rmass(ibin,ig)*pc(iz,ibin,ip)*evaplg(ibin,ig)*(1._f+1.e-8_f))then
          rest = rmass(ibin,ig)*pc(iz,ibin,ip)*evaplg(ibin,ig)
        else if(rest .gt. rmass(ibin-1,ig)*pc(iz,ibin,ip)*evaplg(ibin,ig)*(1._f-1.e-8_f))then
          rest = rmass(ibin-1,ig)*pc(iz,ibin,ip)*evaplg(ibin,ig)
        else
          write(*,*) "rest is wrong: rest",rest,"upper bin",rmass(ibin,ig)*pc(iz,ibin,ip)*evaplg(ibin,ig),"lower bin",rmass(ibin-1,ig)*pc(iz,ibin,ip)*evaplg(ibin,ig)
          rc = RC_ERROR
          return
        end if
      end if

      ! then splite the particle between bin and bin-1
      evappe(ibin-1,ip) = evappe(ibin-1,ip)+(rmass(ibin,ig)*pc(iz,ibin,ip)*evaplg(ibin,ig)-rest)/(rmass(ibin,ig)-rmass(ibin-1,ig))
      evappe(ibin,ip) = evappe(ibin,ip)+(rest-rmass(ibin-1,ig)*pc(iz,ibin,ip)*evaplg(ibin,ig))/(rmass(ibin,ig)-rmass(ibin-1,ig))
!      write(*,*) "ip",ip,"pc mass ",pc(iz,ibin,ip)*evaplg(ibin,ig)*rmass(ibin,ig)-sumtogas+ndelemmass*pc(iz,ibin,ip)*evaplg(ibin,ig)*0.1_f
!      write(*,*) "add mass",(rmass(ibin,ig)*pc(iz,ibin,ip)*evaplg(ibin,ig)-rest)/(rmass(ibin,ig)-rmass(ibin-1,ig))*rmass(ibin-1,ig)+(rest-rmass(ibin-1,ig)*pc(iz,ibin,ip)*evaplg(ibin,ig))/(rmass(ibin,ig)-rmass(ibin-1,ig))*rmass(ibin,ig)
!      write(*,*) "pc num ",pc(iz,ibin,ip)*evaplg(ibin,ig)
!      write(*,*) "add num",(rmass(ibin,ig)*pc(iz,ibin,ip)*evaplg(ibin,ig)-rest)/(rmass(ibin,ig)-rmass(ibin-1,ig))+(rest-rmass(ibin-1,ig)*pc(iz,ibin,ip)*evaplg(ibin,ig))/(rmass(ibin,ig)-rmass(ibin-1,ig))
!      write(*,*) "evappe(ibin-1,ip)",evappe(ibin-1,ip)
!      write(*,*) "evappe(ibin,ip)",evappe(ibin,ip)

      do ie = ip, ip+nelemg(ig)-1
        if(growindex(iindex) .eq. ie)then
          do iindex = 1, index_num
            if(itype(ie) .eq. I_GROWCORE)then
              evappe(ibin-1,ie) = evappe(ibin-1,ie)+(rmass(ibin,ig)*pc(iz,ibin,ip)*evaplg(ibin,ig) &
                                -rest)/(rmass(ibin,ig)-rmass(ibin-1,ig))*rmass(ibin-1,ig)*(remove(iindex)-togas(iindex))/rest
              evappe(ibin,ie) = evappe(ibin,ie)+(rest-rmass(ibin-1,ig)*pc(iz,ibin,ip) &
                                *evaplg(ibin,ig))/(rmass(ibin,ig)-rmass(ibin-1,ig))*rmass(ibin,ig)*(remove(iindex)-togas(iindex))/rest
!              write(*,*) "ie",ie
!              write(*,*) "pc(iz,ibin,ie)*evaplg(ibin,ig)",pc(iz,ibin,ie)*evaplg(ibin,ig)-togas(iindex)
!              write(*,*) "add...",(rmass(ibin,ig)*pc(iz,ibin,ip)*evaplg(ibin,ig) &
!                -rest)/(rmass(ibin,ig)-rmass(ibin-1,ig))*rmass(ibin-1,ig)*(remove(iindex)-togas(iindex))/rest+(rest-rmass(ibin-1,ig)*pc(iz,ibin,ip) &
!                *evaplg(ibin,ig))/(rmass(ibin,ig)-rmass(ibin-1,ig))*rmass(ibin,ig)*(remove(iindex)-togas(iindex))/rest
!              write(*,*) "evappe(ibin-1,ie)",evappe(ibin-1,ie)
!              write(*,*) "evappe(ibin,ie)",evappe(ibin,ie)
            end if
          end do
        end if
        if(itype(ie) .eq. I_COREMASS)then
          evappe(ibin-1,ie) = evappe(ibin-1,ie)+(rmass(ibin,ig)*pc(iz,ibin,ip)*evaplg(ibin,ig)-rest)/(rmass(ibin,ig)-rmass(ibin-1,ig))*rmass(ibin-1,ig)*(pc(iz,ibin,ie)*evaplg(ibin,ig))/rest
          evappe(ibin,ie) = evappe(ibin,ie)+(rest-rmass(ibin-1,ig)*pc(iz,ibin,ip)*evaplg(ibin,ig))/(rmass(ibin,ig)-rmass(ibin-1,ig))*rmass(ibin,ig)*(pc(iz,ibin,ie)*evaplg(ibin,ig))/rest
!          write(*,*) "ie",ie
!          write(*,*) "pc(iz,ibin,ie)*evaplg(ibin,ig)",pc(iz,ibin,ie)*evaplg(ibin,ig)
!          write(*,*) "add...",(rmass(ibin,ig)*pc(iz,ibin,ip)*evaplg(ibin,ig)-rest)/(rmass(ibin,ig)-rmass(ibin-1,ig))*rmass(ibin-1,ig)*(pc(iz,ibin,ie)*evaplg(ibin,ig))/rest+(rest-rmass(ibin-1,ig)*pc(iz,ibin,ip)*evaplg(ibin,ig))/(rmass(ibin,ig)-rmass(ibin-1,ig))*rmass(ibin,ig)*(pc(iz,ibin,ie)*evaplg(ibin,ig))/rest
!          write(*,*) "evappe(ibin-1,ie)",evappe(ibin-1,ie)
!          write(*,*) "evappe(ibin,ie)",evappe(ibin,ie)
        end if
      enddo

    else
!      write(*,*) "here 2..."
      rest = rmass(ibin,ig)*pc(iz,ibin,ip)*evaplg(ibin,ig)-sumtogas

      if(rest .gt. rmass(ibin,ig)*pc(iz,ibin,ip)*evaplg(ibin,ig) .or. rest .lt. rmass(ibin-1,ig)*pc(iz,ibin,ip)*evaplg(ibin,ig))then
        if(rest .gt. rmass(ibin,ig)*pc(iz,ibin,ip)*evaplg(ibin,ig) .and. rest .lt. rmass(ibin,ig)*pc(iz,ibin,ip)*evaplg(ibin,ig)*(1._f+1.e-8_f))then
          rest = rmass(ibin,ig)*pc(iz,ibin,ip)*evaplg(ibin,ig)
        else if(rest .lt. rmass(ibin-1,ig)*pc(iz,ibin,ip)*evaplg(ibin,ig) .and. rest .gt. rmass(ibin-1,ig)*pc(iz,ibin,ip)*evaplg(ibin,ig)*(1._f-1.e-8_f))then
          rest = rmass(ibin-1,ig)*pc(iz,ibin,ip)*evaplg(ibin,ig)
        else
          write(*,*) "rest is wrong: rest",rest,"upper bin",rmass(ibin,ig)*pc(iz,ibin,ip)*evaplg(ibin,ig),"lower bin",rmass(ibin-1,ig)*pc(iz,ibin,ip)*evaplg(ibin,ig)
          rc = RC_ERROR
          return
        end if
      end if

      evappe(ibin-1,ip) = evappe(ibin-1,ip)+(rmass(ibin,ig)*pc(iz,ibin,ip) &
                        *evaplg(ibin,ig)-rest)/(rmass(ibin,ig)-rmass(ibin-1,ig))
      evappe(ibin,ip) = evappe(ibin,ip)+(rest-rmass(ibin-1,ig)*pc(iz,ibin,ip) &
                        *evaplg(ibin,ig))/(rmass(ibin,ig)-rmass(ibin-1,ig))
!    if(ibin .eq. 13 .and. ig .eq. 5)then
!      write(*,*) "ip",ip,"pc mass ",pc(iz,ibin,ip)*evaplg(ibin,ig)*rmass(ibin,ig)-sumtogas
!      write(*,*) "add mass",(rmass(ibin,ig)*pc(iz,ibin,ip)*evaplg(ibin,ig)-rest)/(rmass(ibin,ig)-rmass(ibin-1,ig))*rmass(ibin-1,ig)+(rest-rmass(ibin-1,ig)*pc(iz,ibin,ip)*evaplg(ibin,ig))/(rmass(ibin,ig)-rmass(ibin-1,ig))*rmass(ibin,ig)
!      write(*,*) "pc num ",pc(iz,ibin,ip)*evaplg(ibin,ig)
!      write(*,*) "add num",(rmass(ibin,ig)*pc(iz,ibin,ip)*evaplg(ibin,ig)-rest)/(rmass(ibin,ig)-rmass(ibin-1,ig))+(rest-rmass(ibin-1,ig)*pc(iz,ibin,ip)*evaplg(ibin,ig))/(rmass(ibin,ig)-rmass(ibin-1,ig)), "to ibin-1",(rmass(ibin,ig)*pc(iz,ibin,ip) &
!      *evaplg(ibin,ig)-rest)/(rmass(ibin,ig)-rmass(ibin-1,ig)),"to ibin",(rest-rmass(ibin-1,ig)*pc(iz,ibin,ip) &
!      *evaplg(ibin,ig))/(rmass(ibin,ig)-rmass(ibin-1,ig))
!      write(*,*) "rest, rmass-1,rmass",rest,rmass(ibin-1,ig)*pc(iz,ibin,ip)*evaplg(ibin,ig),rmass(ibin,ig)*pc(iz,ibin,ip)*evaplg(ibin,ig)
!      write(*,*) "evappe(ibin-1,ip)",evappe(ibin-1,ip)
!      write(*,*) "evappe(ibin,ip)",evappe(ibin,ip)
!     end if

      do ie = ip, ip+nelemg(ig)-1
        if(itype(ie) .eq. I_GROWCORE)then
          do iindex = 1, index_num
            if(growindex(iindex) .eq. ie)then
              evappe(ibin-1,ie) = evappe(ibin-1,ie)+(rmass(ibin,ig)*pc(iz,ibin,ip) &
                            *evaplg(ibin,ig)-rest)/(rmass(ibin,ig)-rmass(ibin-1,ig))*rmass(ibin-1,ig)*(remove(iindex)-togas(iindex))/rest
              evappe(ibin,ie) = evappe(ibin,ie)+(rest-rmass(ibin-1,ig)*pc(iz,ibin,ip) &
                            *evaplg(ibin,ig))/(rmass(ibin,ig)-rmass(ibin-1,ig))*rmass(ibin,ig)*(remove(iindex)-togas(iindex))/rest
!              write(*,*) "ie",ie,"ibin",ibin
!              write(*,*) "pc(iz,ibin,ie)*evaplg(ibin,ig)",pc(iz,ibin,ie)*evaplg(ibin,ig)-togas(iindex)
!              write(*,*) "add...",(rmass(ibin,ig)*pc(iz,ibin,ip) &
!                *evaplg(ibin,ig)-rest)/(rmass(ibin,ig)-rmass(ibin-1,ig))*rmass(ibin-1,ig)*(remove(iindex)-togas(iindex))/rest+(rest-rmass(ibin-1,ig)*pc(iz,ibin,ip) &
!                *evaplg(ibin,ig))/(rmass(ibin,ig)-rmass(ibin-1,ig))*rmass(ibin,ig)*(remove(iindex)-togas(iindex))/rest
!              write(*,*) "evappe(ibin-1,ie)",evappe(ibin-1,ie)
!              write(*,*) "evappe(ibin,ie)",evappe(ibin,ie)
            endif
          enddo
        end if
        if(itype(ie) .eq. I_COREMASS)then
          evappe(ibin-1,ie) = evappe(ibin-1,ie)+(rmass(ibin,ig)*pc(iz,ibin,ip)*evaplg(ibin,ig)-rest)/(rmass(ibin,ig)-rmass(ibin-1,ig))*rmass(ibin-1,ig)*(pc(iz,ibin,ie)*evaplg(ibin,ig))/rest
          evappe(ibin,ie) = evappe(ibin,ie)+(rest-rmass(ibin-1,ig)*pc(iz,ibin,ip)*evaplg(ibin,ig))/(rmass(ibin,ig)-rmass(ibin-1,ig))*rmass(ibin,ig)*(pc(iz,ibin,ie)*evaplg(ibin,ig))/rest
!          write(*,*) "ie",ie
!          write(*,*) "pc(iz,ibin,ie)*evaplg(ibin,ig)",pc(iz,ibin,ie)*evaplg(ibin,ig)
!          write(*,*) "add...",(rmass(ibin,ig)*pc(iz,ibin,ip)*evaplg(ibin,ig)-rest)/(rmass(ibin,ig)-rmass(ibin-1,ig))*rmass(ibin-1,ig)*(pc(iz,ibin,ie)*evaplg(ibin,ig))/rest+(rest-rmass(ibin-1,ig)*pc(iz,ibin,ip)*evaplg(ibin,ig))/(rmass(ibin,ig)-rmass(ibin-1,ig))*rmass(ibin,ig)*(pc(iz,ibin,ie)*evaplg(ibin,ig))/rest
!          write(*,*) "evappe(ibin-1,ie)",evappe(ibin-1,ie)
!          write(*,*) "evappe(ibin,ie)",evappe(ibin,ie)
        end if
      enddo
    endif
  endif

  do ie = ip, ip+nelemg(ig)-1
!    if(ig .eq. 5)write(*,*) "ibin",ibin,"ie",ie,"evappe(ibin,ie)",evappe(ibin,ie),"evappe(ibin-1,ie)",evappe(ibin-1,ie)
    if(evappe(ibin,ie) .lt. 0._f .or. evappe(ibin-1,ie) .lt. 0._f)then
      write(*,*) "negative evappe: ibin ",ibin,"ie ",ie,"evappe ",evappe(ibin,ie)
      rc = RC_ERROR
      return
    end if
  end do

  return
end
