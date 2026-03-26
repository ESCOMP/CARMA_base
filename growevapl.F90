! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!! This routine evaluate particle loss rates due to condensational
!! growth and evaporation for all condensing gases.
!!
!! The loss rates for each group are <growlg> and <evaplg>.
!!
!! Units are [s^-1].
!!
!! @author Andy Ackerman
!! @version Dec-1995
subroutine growevapl(carma, cstate, iz, rc)

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
  integer, intent(inout)               :: rc      !! return code, negative indicates failure

  ! Local declarations
  integer                        :: igroup
  integer                        :: iepart
  integer                        :: igas,icoregas,icore
  integer                        :: ibin
  integer                        :: isol
  integer                        :: nother
  integer                        :: ieoth_rel
  integer                        :: ieoth_abs
  integer                        :: jother
  integer                        :: iesub,ielemg
  real(kind=f)                   :: argsol
  real(kind=f)                   :: othermtot
  real(kind=f)                   :: condm
  real(kind=f)                   :: akas
  real(kind=f)                   :: expon
  real(kind=f)                   :: g0
  real(kind=f)                   :: g1
  real(kind=f)                   :: g2
  real(kind=f)                   :: ss
  real(kind=f)                   :: pvap
  real(kind=f)                   :: dpc
  real(kind=f)                   :: dpc1
  real(kind=f)                   :: dpcm1
  real(kind=f)                   :: rat1
  real(kind=f)                   :: rat2
  real(kind=f)                   :: rat3
  real(kind=f)                   :: rat4
  real(kind=f)                   :: ratt1
  real(kind=f)                   :: ratt2
  real(kind=f)                   :: ratt3
  real(kind=f)                   :: den1
  real(kind=f)                   :: test1
  real(kind=f)                   :: test2
  real(kind=f)                   :: x
  integer                        :: ieother(NELEM)
  real(kind=f)                   :: otherm(NELEM)
  real(kind=f)                   :: dela(NBIN)
  real(kind=f)                   :: delma(NBIN)
  real(kind=f)                   :: aju(NBIN)
  real(kind=f)                   :: ar(NBIN)
  real(kind=f)                   :: al(NBIN)
  real(kind=f)                   :: a6(NBIN)
  real(kind=f)                   :: dmdt(NBIN)
  real(kind=f)                   :: dmdt_grow(NBIN),dmdt_evap(NBIN)  

  dmdt(:) = 0.0_f               ! total growth rate
  dmdte(:,:,:) = 0.0_f       ! growcore growth rate
  dmdtmain(:,:) = 0.0_f    ! number density element growth rate
  dmdt_grow(:) = 0.0_f          ! particle growth rate
  dmdt_evap(:) = 0.0_f          ! particle evaporation rate

  do igroup = 1,NGROUP

    ! element of particle number concentration 
    iepart = ienconc(igroup)
    
    ! condensing gas
    igas = igrowgas(iepart)

    ! if there are no condensing gas for particle number concentration, 
    ! check if there is condensing gas due to nucleation
    if( igas .eq. 0 )then
      icoregas = igrowgas(iepart+ncore(igroup))
    end if

    if (igas .ne. 0 .or. icoregas .ne. 0) then
      ! Only valid for condensing liquid water and sulfric acid currently.
      if ((igas /= igash2o) .and. (igas .ne. igash2so4) .and. (icoregas /= igash2o) .and. (icoregas .ne. igashno3) .and. (icoregas .ne. igasnh3)) then
        if (do_print) write(LUNOPRT,*) 'growevapl::ERROR - Invalid gas (', igas, ').'
        rc = -1
        return
      endif           

      !write(*,*) "checkpoint1",igroup, pconmax(iz,igroup),FEW_PC

      ! Treat condensation of gas <igas> to/from particle group <igroup>.
      !
      ! Bypass calculation if few particles are present 
      if( pconmax(iz,igroup) .gt. FEW_PC )then
      
        ! initialize the growth rate for different purpose
        dmdt(:) = 0.0_f               ! total growth rate
        dmdte(:,:,igroup) = 0.0_f       ! growcore growth rate
        dmdtmain(:,igroup) = 0.0_f    ! number density element growth rate
        dmdt_grow(:) = 0.0_f          ! particle growth rate
        dmdt_evap(:) = 0.0_f          ! particle evaporation rate

        do ibin = 1,NBIN-1

          ! initialize the growth rate for different purpose
          dmdt(ibin) = 0.0_f               ! total growth rate
          dmdte(ibin,:,igroup) = 0.0_f       ! growcore growth rate
          dmdtmain(ibin,igroup) = 0.0_f    ! number density element growth rate
          dmdt_grow(ibin) = 0.0_f          ! particle growth rate
          dmdt_evap(ibin) = 0.0_f          ! particle evaporation rate

          ! Determine the growth rate (dmdt). This calculation may take into account
          ! radiative effects on the particle which can affect the growth rates.
          if(igas .eq. 0)then
            dmdt(ibin) = 0._f
            dmdtmain(ibin,igroup) = 0._f
          else
            call pheat(carma, cstate, iz, igroup, iepart, ibin, igas, dmdt(ibin), rc)
            dmdtmain(ibin,igroup) = dmdt(ibin)    ! growth rate for number density element
          end if

          ! zhu 2011 10 19
          ! Include primary's growth rate in appropriate dmdt term
          if (dmdt(ibin) .gt. 0.0_f) then
            dmdt_grow(ibin) = dmdtmain(ibin,igroup)
          else
            dmdt_evap(ibin) = dmdtmain(ibin,igroup)
          endif

          ! calculate growcore element growth rate(dmdt).
          ! Growth rate for particle in <ibin> is sum of growth rates for each gas.
          ! Calculate <dmdt_gro> and <dmdte_gro> for use later to get mass fraction rate of change.

          if(nelemg(igroup) .gt. 1)then
          do ielemg = 2,nelemg(igroup)
            iesub = ielemg - 1 + iepart
            icore = ielemg - 1
!           ! yu test, calculate NH4 first, then NO3
!            iesub = nelemg(igroup) + 2 - ielemg - 1 + iepart
!            icore = nelemg(igroup) + 2 - ielemg - 1
!           ! end yu test

            if (itype(iesub) .eq. I_GROWCORE) then
              icoregas = igrowgas(iesub)    !growcore gas

              if(icoregas .eq. 0)then
                if (do_print) write(LUNOPRT,*) 'growevapl::ERROR - Invalid gas for GROWCORE.'
                rc = RC_ERROR
                return
              end if

              !write(*,*) "iepart",iepart,"icoregas",icoregas,"icore",icore
              ! in cold environment, uses Luo and Tabazadeh (2001) for UTLS
              call pheat(carma, cstate, iz, igroup, iepart, ibin, icoregas, dmdte(ibin,icore,igroup), rc)
              !write(*,*) "igroup",igroup,"iepart",iepart,"icoregas",icoregas,"dmdte(ibin,icore,igroup)",dmdte(ibin,icore,igroup)
             
              ! NH3 condensation / evaporation
              if (icoregas .eq. igasnh3) then
                 !
                 dmdte(ibin,icore,igroup) = 0._f
                 !dmdte(ibin,icore,igroup) = max(0._f,dmdtmain(ibin,igroup)*17._f/98._f) + dmdtnh3b(ibin)
                 call pheat(carma, cstate, iz, igroup, iepart, ibin, icoregas, dmdte(ibin,icore,igroup), rc)
                 !write(*,*) "icore",icore,"ibin",ibin,"igroup",igroup
                 !write(*,*) "icoregas",icoregas,"dmdte(ibin,icore,igroup)",dmdte(ibin,icore,igroup)

                 !
              endif
              
              ! add constrains following Yunqian's suggestion
              !if (dmdte(ibin,icore,igroup) .lt. 0._f .and. pc(iz,ibin,iesub)/pc(iz,ibin,iepart)/rmass(ibin,igroup) .le. 1.e-2_f) then
              if (dmdte(ibin,icore,igroup) .lt. 0._f .and. sum(pc(iz,:,iesub))/sum(pc(iz,:,iepart)*rmass(:,igroup)) .le. 1.e-2_f) then
                 dmdte(ibin,icore,igroup) = 0._f
              endif
        !write(*,*) "igroup",igroup,"dmdtmain",dmdtmain(ibin,igroup)
        !if(iesub .eq. 7) write(*,*) "igroup",igroup,"dmdte",dmdte(ibin,icore,igroup)


              dmdt(ibin) = dmdtmain(ibin,igroup) + dmdte(ibin,icore,igroup)      ! growth rate of whole particle

              ! Total dmdt's for growth and evap
              if (dmdte(ibin,icore,igroup) .gt. 0.0_f) then
                dmdt_grow(ibin) = dmdt_grow(ibin) + dmdte(ibin,icore,igroup)
              else
                dmdt_evap(ibin) = dmdt_evap(ibin) + dmdte(ibin,icore,igroup)
              endif
            endif   ! (itype(iesub) .eq. I_GROWCORE)
          enddo    ! sub elements in group
          end if


        enddo     ! ibin = 1,NBIN-1


        ! Now calculate condensation/evaporation production and loss rates.
        ! Use Piecewise Polynomial Method [Colella and Woodward, J. Comp. Phys.,
        ! 54, 174-201, 1984]
        !
        ! First, use cubic fits to estimate concentration values at bin
        ! boundaries
        do ibin = 2,NBIN-1

          dpc = pc(iz,ibin,iepart) / dm(ibin,igroup)
          dpc1 = pc(iz,ibin+1,iepart) / dm(ibin+1,igroup)
          dpcm1 = pc(iz,ibin-1,iepart) / dm(ibin-1,igroup)
          ratt1 = pratt(1,ibin,igroup)
          ratt2 = pratt(2,ibin,igroup)
          ratt3 = pratt(3,ibin,igroup)
          dela(ibin) = ratt1 * ( ratt2*(dpc1-dpc) + ratt3*(dpc-dpcm1) )
          delma(ibin) = 0._f

          if( (dpc1-dpc)*(dpc-dpcm1) .gt. 0._f ) &
            delma(ibin) = min( abs(dela(ibin)), 2._f*abs(dpc-dpc1), &
                 2._f*abs(dpc-dpcm1) ) * sign(1._f, dela(ibin))

        enddo     ! ibin = 2,NBIN-2

        do ibin = 2,NBIN-2

          dpc = pc(iz,ibin,iepart) / dm(ibin,igroup)
          dpc1 = pc(iz,ibin+1,iepart) / dm(ibin+1,igroup)
          dpcm1 = pc(iz,ibin-1,iepart) / dm(ibin-1,igroup)
          rat1 = prat(1,ibin,igroup)
          rat2 = prat(2,ibin,igroup)
          rat3 = prat(3,ibin,igroup)
          rat4 = prat(4,ibin,igroup)
          den1 = pden1(ibin,igroup)

          ! <aju(ibin)> is the estimate for concentration (dn/dm) at bin
          ! boundary <ibin>+1/2.
          aju(ibin) = dpc + rat1*(dpc1-dpc) + 1._f/den1 * &
                  ( rat2*(rat3-rat4)*(dpc1-dpc) - &
                  dm(ibin,igroup)*rat3*delma(ibin+1) + &
                  dm(ibin+1,igroup)*rat4*delma(ibin) )
        enddo     ! ibin = 2,NBIN-2

        ! Now construct polynomial functions in each bin
        do ibin = 3,NBIN-2
          al(ibin) = aju(ibin-1)
          ar(ibin) = aju(ibin)
        enddo

        ! Use linear functions in first two and last two bins
        if( NBIN .gt. 1 )then
          ibin = NBIN
          
          ar(2) = aju(2)
          al(2) = pc(iz,1,iepart)/dm(1,igroup) + &
                      palr(1,igroup) * &
                      (pc(iz,2,iepart)/dm(2,igroup)- &
                      pc(iz,1,iepart)/dm(1,igroup))
          ar(1) = al(2)
          al(1) = pc(iz,1,iepart)/dm(1,igroup) + &
                      palr(2,igroup) * &
                      (pc(iz,2,iepart)/dm(2,igroup)- &
                      pc(iz,1,iepart)/dm(1,igroup))
          
          al(ibin-1) = aju(ibin-2)
          ar(ibin-1) = pc(iz,ibin-1,iepart)/dm(ibin-1,igroup) + &
                      palr(3,igroup) * &
                      (pc(iz,ibin,iepart)/dm(ibin,igroup)- &
                      pc(iz,ibin-1,iepart)/dm(ibin-1,igroup))
          al(ibin) = ar(ibin-1)
          ar(ibin) = pc(iz,ibin-1,iepart)/dm(ibin-1,igroup) + &
                      palr(4,igroup) * &
                      (pc(iz,ibin,iepart)/dm(ibin,igroup)- &
                      pc(iz,ibin-1,iepart)/dm(ibin-1,igroup))
        endif

        ! Next, ensure that polynomial functions do not deviate beyond the
        ! range [<al(ibin)>,<ar(ibin)>]
        do ibin = 1,NBIN

          dpc = pc(iz,ibin,iepart) / dm(ibin,igroup)

          if( (ar(ibin)-dpc)*(dpc-al(ibin)) .le. 0._f )then
            al(ibin) = dpc
            ar(ibin) = dpc
          endif

          test1 = (ar(ibin)-al(ibin))*(dpc - 0.5_f*(al(ibin)+ar(ibin)))
          test2 = 1._f/6._f*(ar(ibin)-al(ibin))**2

          if( test1 .gt. test2 )then
             al(ibin) = 3._f*dpc - 2._f*ar(ibin)
          elseif( test1 .lt. -test2 )then
             ar(ibin) = 3._f*dpc - 2._f*al(ibin)
          endif
        enddo

        !  Lastly, calculate fluxes across each bin boundary.
        !
        !  Use upwind advection when courant number > 1.
        do ibin = 1,NBIN
          dpc = pc(iz,ibin,iepart) / dm(ibin,igroup)
          dela(ibin) = ar(ibin) - al(ibin)
          a6(ibin) = 6._f * ( dpc - 0.5_f*(ar(ibin)+al(ibin)) )
        enddo

        do ibin = 1,NBIN-1

          if( dmdt_grow(ibin) .gt. 0._f .and. &
              pc(iz,ibin,iepart) .gt. SMALL_PC )then

            x = dmdt_grow(ibin)*dtime/dm(ibin,igroup)

            if( x .lt. 1._f )then
              growlg(ibin,igroup) = dmdt_grow(ibin)/pc(iz,ibin,iepart) &
                       * ( ar(ibin) - 0.5*dela(ibin)*x + &
                       (x/2._f - x**2/3._f)*a6(ibin) )
            else
              growlg(ibin,igroup) = dmdt_grow(ibin) / dm(ibin,igroup)
            endif

          elseif( dmdt_evap(ibin) .lt. 0._f .and. &
              pc(iz,ibin+1,iepart) .gt. SMALL_PC )then

            x = -dmdt_evap(ibin)*dtime/dm(ibin+1,igroup)

            if( x .lt. 1._f )then
              evaplg(ibin+1,igroup) = -dmdt_evap(ibin)/ &
                      pc(iz,ibin+1,iepart) &
                      * ( al(ibin+1) + 0.5_f*dela(ibin+1)*x + &
                      (x/2._f - (x**2)/3._f)*a6(ibin+1) )
            else
              evaplg(ibin+1,igroup) = -dmdt_evap(ibin) / dm(ibin+1,igroup)
            endif

            ! Boundary conditions: for evaporation out of first bin (with cores), 
            ! use evaporation rate from second bin.
!            if( ibin .eq. 1 .and. ncore(igroup) .gt. 0 )then
            if( ibin .eq. 1)then
              evaplg(1,igroup) = -dmdt_evap(1) / dm(1,igroup)
            endif
          endif

          !if( ibin .eq. NBIN-1 .and. dmdt_evap(ibin) .lt. 0._f ) then
          !  evaplg(ibin+1,igroup) = -dmdt_evap(ibin+1) / dm(ibin+1,igroup)
          !end if

          !zhu 2010 10 02 limit the growth rate to "reasonable" values
          !growlg(ibin,igroup) = min(growlg(ibin,igroup),1.e10_f)
          !evaplg(ibin+1,igroup) = min(evaplg(ibin+1,igroup),1.e10_f)
          !evaplg(1,igroup) = min(evaplg(1,igroup),1.e10_f)

        enddo    ! ibin = 1,NBIN-1
      endif     ! (pconmax .gt. FEW_PC)
    endif      ! (igas = igrowgas(ielem)) .ne. 0 
  enddo       ! igroup = 1,NGROUP

  !write(*,*) "dmdt",dmdt,"growlg",growlg,"evaplg",evaplg
  
  ! Return to caller with particle loss rates for growth and evaporation
  ! evaluated.
   return
end
