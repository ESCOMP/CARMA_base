! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!!  Calculates the wtpct for HNO3 and NH4+H2SO4.
!!  wtpct_h2so4 is NH4+H2SO4 when NH4 is present based on Lin&Tabazadeh 2001
!!  wtpct_h2so4 is H2SO4 when NH4 is not present
!!
!! @author Yunqian
!! @version Aug-2012

subroutine wtpct_luo_total(carma, cstate, iz, rc, wt_h2so4, wt_hno3)
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
  real(kind=f), intent(out)            :: wt_h2so4 
  real(kind=f), intent(out)            :: wt_hno3 
  integer, intent(inout)               :: rc        !! return code, negative indicates failure

  ! Local declarations
  integer                        :: bin_pick
  integer                        :: ielem, hno3group, iepart, igroup, pursulgroup, pursulelem
  real(kind=f)                   :: ratio   !the ratio of HNO3/H2SO4 inside particles
  real(kind=f)                   :: gc_cgs  ! water vapor mass concentration [g/cm3]
  real(kind=f)                   :: pc_hno3, pc_h2so4, pc_oc, pc_bc, pc_dst, pc_slt, pc_nh4, pc_soa
  real(kind=f)                   :: ph2o, ph2olog
  real(kind=f)                   :: a1,b1,c1,d1, AA,BB,CC
  real(kind=f)                   :: Delta
  real(kind=f)                   :: x1,x2,x3
  real(kind=f)                   :: Y1,Y2
  real(kind=f)                   :: K
  real(kind=f)                   :: TT,Theta1,temp
  real(kind=f)                   :: rat1

 
  wt_h2so4 = 0.0_f
  wt_hno3 = 0.0_f
  ratio = 0.0_f  
  rat1 = 1.0_f

  ! select the group with hno3 as condensing gas 
  !(only works if there's one group condenses HNO3, 
  ! need to modify the code if multiple groups condense HNO3)
  do ielem = 1,NELEM
    if( igrowgas(ielem) .eq. igashno3)then
      hno3group = igelem(ielem)
    end if
  end do

  ! sum and calculate the total mass for different particle composition
  pc_hno3 = 0._f
  pc_soa = 0._f
  pc_oc  = 0._f
  pc_bc  = 0._f
  pc_dst = 0._f
  pc_slt = 0._f
  pc_nh4 = 0._f
  pc_h2so4 = 0._f

  iepart = ienconc(hno3group) !the concentration element of the group that contains hno3

  do ielem = iepart, iepart+ncore(hno3group)
    if(icomp(ielem) .eq. I_H2SO4)then
      pc_h2so4 = pc_h2so4+ sum(cstate%f_pc(iz,:,ielem)*rmass(:,hno3group))
    elseif(icomp(ielem) .eq. I_HNO3)then
      pc_hno3 = pc_hno3+ sum(cstate%f_pc(iz,:,ielem))
    elseif(icomp(ielem) .eq. I_SOA1 .or. icomp(ielem) .eq. I_SOA2 .or. icomp(ielem) .eq. I_SOA3 &
        .or. icomp(ielem) .eq. I_SOA4 .or. icomp(ielem) .eq. I_SOA5 .or. icomp(ielem) .eq. I_SOA)then
      pc_soa = pc_soa + sum(cstate%f_pc(iz,:,ielem))
    elseif(icomp(ielem) .eq. I_OC)then
      pc_oc  = pc_oc + sum(cstate%f_pc(iz,:,ielem))
    elseif (icomp(ielem) .eq. I_BC)then
      pc_bc  = pc_bc + sum(cstate%f_pc(iz,:,ielem))
    elseif (icomp(ielem) .eq. I_DUST)then
      pc_dst = pc_dst + sum(cstate%f_pc(iz,:,ielem))
    elseif (icomp(ielem) .eq. I_SALT)then
      pc_slt = pc_slt + sum(cstate%f_pc(iz,:,ielem))
    elseif (icomp(ielem) .eq. I_NH3)then
      pc_nh4 = pc_nh4 + sum(cstate%f_pc(iz,:,ielem))
    else
      if (do_print) write(LUNOPRT,*) "unidenfied element for wtpct calculation, composition ",icomp(ielem)
      rc = RC_ERROR
      return
    end if 
  end do

  ! add pure sulfate group and mixed sulfate group total H2SO4 together
  do igroup = 1,NGROUP
    if (nelemg(igroup) .eq. 1 .and. icomp(ienconc(igroup)) .eq. I_H2SO4)then
       pursulgroup = igroup
       pursulelem = ienconc(igroup)
    end if
  end do
  
  pc_h2so4 = sum(cstate%f_pc(iz,:,hno3group)*rmass(:,hno3group)) + &
       			  sum(cstate%f_pc(iz,:,pursulelem)*rmass(:,pursulgroup)) - &
       			  pc_oc-pc_bc-pc_dst-pc_slt-pc_nh4-pc_soa-pc_hno3

  !write(*,*) "pc_h2so4",pc_h2so4,"pc_hno3",pc_hno3

  ! Don't allow saturation vapor pressure to underflow at very low temperatures
  temp=max(t(iz),140._f)
 
  ! Convert water vapor concentration to g/cm3:
  gc_cgs = gc(iz, igash2o) / zmet(iz)
       
  ! Using Luo's equation to calculate the weight percent of H2SO4 and HNO3
  ! If NH4 is zero, it is identical with Luo 1995
  ! If NH4 is larger than zero, we consider NH4 in the solution as Lin&Tabazadeh2001
  if (pc_h2so4 .gt. 0.0_f) then
    !ratio= pc_hno3/(pc_h2so4+pc_nh4)
    ratio= pc_hno3/pc_h2so4
    !rat1= (pc_h2so4+pc_nh4)/pc_h2so4
    !rat1 = max(1._f,rat1) 
    !write(*,*) "ratio",ratio

    if (ratio .lt. 0.0_f) then
      write(*,*) 'iz',iz,'rmass', rmass(:,2),'ratio',ratio,'pc h2so4', pc_h2so4,'pc hno3',pc_hno3
      write(*,*) 'temperature',t(iz),' pressure[dyne/cm^2] = ', p(iz)
      write(*,*) 'invalid wtpct in wtpct_luo.f'
	  rc = RC_ERROR
      return
	endif
	
	if (ratio .lt. 0.0001_f .and. t(iz) .gt. 235._f) then
	  gc_cgs = gc(iz,1)/ zmet(iz)
	  wt_h2so4 = wtpct_tabaz(carma,t(iz),gc_cgs,pvapl(iz,1),rc)
	  wt_hno3 = wt_h2so4 * ratio
	  return
	endif
	
  elseif(pc_h2so4 .eq. 0.0_f)then
    wt_hno3 = 50._f
	wt_h2so4 = 0.1_f
	return
  else
	wt_hno3 = 0.1_f
	wt_h2so4 = 0.1_f
	return
  endif

  ! Luo1995 first calculate the gas phase vapor pressure
  ! which is equivalent with H2O vapor pressure over the particle
  gc_cgs = gc(iz,1)/ zmet(iz)
  ! h2o pressure in mbar     
  ph2o = max(gc_cgs*RGAS*t(iz)/gwtmol(1)/1000._f,1.e-30_f)
  ph2olog = log(ph2o)

  !%shengjin's formulas
  ! solve for H2SO4 w2 wtpct
  a1 = (ratio+1.4408_f)**2._f*(4._f*ratio+8.19_f)
  b1 = 1._f/t(iz)*(ratio+1.4408_f)*(2309._f*ratio+1876.7_f)- &
           (ratio+1.4408_f)*(ratio*7.451_f+12._f)
  c1 = 4.5261_f * ratio + 5.3465_f - 1._f / t(iz) * (1033._f * ratio + 928.9_f)
  d1 = ph2olog+1._f/t(iz)*5814._f-23.306_f

  AA = b1**2._f - 3._f*a1*c1
  BB = b1*c1 - 9._f*a1*d1
  CC = c1**2._f - 3._f*b1*d1

  Delta = BB**2._f-4._f*AA*CC

  if (AA .eq. 0._f .and. BB .eq. 0._f)then
    x1= - b1/3._f/a1
    x2= x1
    x3= x2
  elseif(Delta .gt. 0._f)then
    Y1= AA*b1+3._f*a1*(-BB+Delta**(1._f/2._f))/2._f
    Y2= AA*b1+3._f*a1*(-BB-Delta**(1._f/2._f))/2._f
	if (Y1 .gt. 0._f .and. Y2 .gt. 0._f) then
      x1= (-b1- Y1**(1._f/3._f)- Y2**(1._f/3._f))/(3._f*a1)
    elseif(Y1 .gt. 0._f) then
	  x1= (-b1- Y1**(1._f/3._f)+ (-Y2)**(1._f/3._f))/(3._f*a1)
	elseif (Y2 .gt. 0._f) then
	  x1= (-b1+ (-Y1)**(1._f/3._f)- Y2**(1._f/3._f))/(3._f*a1)
	else
	  x1= (-b1+ (-Y1)**(1._f/3._f)+ (-Y2)**(1._f/3._f))/(3._f*a1)
	endif
    x2= -999._f;
    x3= -999._f;
  elseif(Delta .eq. 0._f .and. AA .ne. 0._f)then
    K = BB/AA
    x1= -b1/a1+K
    x2= -K/2._f
    x3= x2
  elseif(Delta .lt. 0._f .and. AA .gt. 0._f)then
    TT = (2._f*AA*b1-3._f*a1*BB)/(2._f*sqrt(AA**3._f))
          theta1= acos(TT);
    x1= (-b1+sqrt(AA)*(cos(theta1/3._f)+ &
              sqrt(3._f)*sin(theta1/3._f)))/(3._f*a1)
    x2= (-b1-2._f*sqrt(AA)*cos(theta1/3._f))/(3._f*a1)
    x3= (-b1+sqrt(AA)*(cos(theta1/3._f)- &
              sqrt(3._f)*sin(theta1/3._f)))/(3._f*a1)
  else
    write(*,*)'Warining condition: A=',AA,'Delta=',Delta
	write(*,*) 'iz',iz,'rmass',rmass(:,1),'pc h2so4', pc_h2so4,'pc hno3',pc_hno3
    write(*,*) 'temperature',t(iz)
	write(*,*) '   pressure[dyne/cm^2]  = ', p(iz)
	write(*,*) 'gc_cgs',gc_cgs,'ph2o',ph2o
	write(*,*) 'ph2olog',ph2olog
	write(*,*) 'subroutine parameter:','Y1',Y1,'Y2',Y2
	write(*,*) 'a1',a1,'b1',b1,'c1',c1,'d1',d1,'ratio',ratio
	write(*,*) 'AA',AA,'BB',BB,'CC',CC,'Delta',Delta
	write(*,*) 'sqrt Delta',sqrt(Delta)
	write(*,*)  'Delta**(1./2.)',Delta**(1._f/2._f),'AA^2',AA**2.
	rc = RC_ERROR
	return
  endif
      
  if (x1 .gt. 0._f)then
    if (x2 .gt. 0._f .and. x3 .gt. 0._f)then
	  wt_h2so4 = max(x1*100._f,x2*100._f,x3*100._f)
	elseif (x2 .gt. 0._f .and. x3 .le. 0._f)then
	     
    else
      wt_h2so4 = x1*100._f
    endif
  elseif (x2 .gt. 0._f)then
    if (x3 .gt. 0._f)then
	  wt_h2so4 = x2*100._f
    else
      wt_h2so4 = x2*100._f
    endif
  elseif (x3 .gt. 0._f)then
    wt_h2so4 = x3*100._f
  else
    wt_h2so4 = 1/(ratio+1)*wtpct_tabaz(carma,t(iz),gc_cgs,pvapl(iz,1),rc)
  endif

  wt_hno3 = wt_h2so4 * ratio
	  
  if((wt_hno3+wt_h2so4) .ge. 100.0_f) then
    wt_h2so4 = 1/(ratio+1)*wtpct_tabaz(carma,t(iz),gc_cgs,pvapl(iz,1),rc)
	wt_hno3 = wt_h2so4*ratio
  endif
	  
  if((wt_hno3+wt_h2so4) .gt. 100._f) then
	write(*,*) "wtpct not correct: H2SO4",wt_h2so4,"HNO3",wt_hno3
	write(*,*) wtpct_tabaz(carma,t(iz),gc_cgs,pvapl(iz,1),rc)
	rc = RC_ERROR
    return
  endif

  return
end
