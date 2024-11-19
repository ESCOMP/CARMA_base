! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!!  Calculates particle production rates due to nucleation <rhompe>:
!!  binary homogeneous nucleation of sulfuric acid and water only
!!  Numerical method follows one of the following:
!!   1. Zhao & Turco, JAS, V.26, No.5, 1995.
!!   2. Vehkamaki, H., M. Kulmala, I. Napari, K.E.J. Lehtinen,
!!       C. Timmreck, M. Noppel and A. Laaksonen, 2002,
!!       An improved parameterization for sulfuric acid-water nucleation
!!       rates for tropospheric and stratospheric conditions,
!!       J. Geophys. Res., 107, 4622, doi:10.1029/2002jd002184
!!
!!
!!  @author Mike Mills, Chuck Bardeen
!!  @version May-2022
subroutine sulfnucrate(carma, cstate, iz, igroup, h2o, h2so4, beta1, beta2, radius_cluster, nucbin, nucrate, rc)
  use carma_precision_mod
  use carma_enums_mod
  use carma_constants_mod
  use carma_types_mod
  use carmastate_mod
  use carma_mod
  use sulfate_utils

  implicit none

  type(carma_type), intent(in)         :: carma          !! the carma object
  type(carmastate_type), intent(inout) :: cstate         !! the carma state object
  integer, intent(in)                  :: iz             !! level index
  integer, intent(in)                  :: igroup         !! group index
  real(kind=f), intent(out)            :: h2o            !! H2O concentrations in molec/cm3
  real(kind=f), intent(out)            :: h2so4          !! H2SO4 concentrations in molec/cm3
  real(kind=f), intent(out)            :: beta1
  real(kind=f), intent(out)            :: beta2
  real(kind=f), intent(out)            :: radius_cluster !! critical radius (cm)
  integer, intent(out)                 :: nucbin         !! bin in which nucleation occurs
  real(kind=f), intent(out)            :: nucrate        !! nucleation rate #/x/y/z/s
  integer, intent(inout)               :: rc             !! return code, negative indicates failure

  !  Local declarations
  integer           :: i, ibin, ie
  real(kind=f)      :: nucrate_cgs      ! binary nucleation rate, j (# cm-3 s-1)
  real(kind=f)      :: cnum_h2so4       ! number of h2so4 molecules in the critical nucleus
  real(kind=f)      :: cnum_tot         ! total number of molecules in the critical nucleus
  real(kind=f)      :: rb               ! [erg/mol]
  real(kind=f)      :: h2so4_cgs        ! H2SO4 densities in g/cm3
  real(kind=f)      :: h2o_cgs          ! H2O densities in g/cm3
  real(kind=f)      :: rh               ! relative humidity (0-1)
  real(kind=f)      :: mass_cluster_dry ! dry mass of the cluster ()
  real(kind=f)      :: h2so4_bb         ! bounded value of H2SO4 concentration in molec/cm3
  real(kind=f)      :: temp_bb          ! bounded value of temperature in Kelvins
  real(kind=f)      :: rh_bb            ! bounded value of relative humidity
  real(kind=f)      :: ftry

 5 format(/,'microfast::WARNING - nucleation rate exceeds 5.e1: ie=', i2,', iz=',i4,',lat=', &
              f7.2,',lon=',f7.2, ', rhompe=', e10.3)

  ! default values for outputs
  nucbin  = 1
  nucrate = 0.0_f
  nucrate_cgs = 0.0_f
  radius_cluster = 0.0_f
  mass_cluster_dry = 0.0_f
  ftry = NOTSET

  !--------------------------------------------------------------

  ! beta1 and beta2 are calculated and used in sulfnucrate, and output for use in sulfhetnucrate
  !   kT/(2*Pi*M) = [erg/mol/K]*[K]/[g/mol] = [erg/g] = [cm2/s2]
  !   RB[erg/mol] = RGAS[erg/mol/K] * T[K] / (2Pi)
  rb = RGAS * t(iz) / 2._f / PI

  ! Beta[cm/s] = sqrt(RB[erg/mol] / WTMOL[g/mol])
  beta1 = sqrt(rb / gwtmol(igash2so4)) ! H2SO4
  beta2 = sqrt(rb / gwtmol(igash2o))   ! H2O

  !--------------------------------------------------------------

  ! Compute H2SO4 densities in g/cm3
  h2so4_cgs = gc(iz, igash2so4) / zmet(iz)

  ! Compute H2O densities in g/cm3
  h2o_cgs   = gc(iz, igash2o)   / zmet(iz)

  ! Compute H2SO4 concentrations in molec/cm3
  h2so4 = h2so4_cgs * AVG / gwtmol(igash2so4)

  ! Compute H2O concentrations in molec/cm3
  h2o   = h2o_cgs   * AVG / gwtmol(igash2o)

  ! Compute relative humidity of water wrt liquid water
  rh = (supsatl(iz, igash2o) + 1._f) !* 100._f

  ! Select nucleation method
  select case (sulfnuclmethod)
  case('ZhaoTurco')
    call binary_nuc_zhao1995( carma, cstate, t(iz), wtpct(iz), rh, h2so4, h2so4_cgs, h2o, h2o_cgs, beta1, &
                  nucrate_cgs, mass_cluster_dry, radius_cluster, ftry, rc )
  case('Vehkamaki')
    if (h2so4 >= 1.0e4_f) then
      temp_bb = max( 230.0_f, min( 305.0_f, t(iz) ) )
      rh_bb = max( 1.0e-4_f, min( 1.0_f, rh ) )
      h2so4_bb = max( 1.0e4_f, min( 1.0e11_f, h2so4 ) )

      call binary_nuc_vehk2002( carma, temp_bb, rh_bb, h2so4_bb, nucrate_cgs, &
                  mass_cluster_dry, radius_cluster )
    end if
  case default
    write(LUNOPRT,*)'sulfnucrate: '//trim(sulfnuclmethod)//' nucleation method no recognized'
    rc = RC_ERROR
    return
  end select

  !   Calc bin # of crit nucleus
  if (mass_cluster_dry.lt.rmassup(1,igroup)) then
    nucbin = 1
  else
    nucbin = 2 + int(log(mass_cluster_dry / rmassup(1,igroup)) / log(rmrat(igroup)))
  endif

  ! If none of the bins are large enough for the critical radius, then
  ! no nucleation will occur.
  if (nucbin <= NBIN) then
    ! Scale to #z/s
    nucrate = nucrate_cgs * zmet(iz)
  endif

  return
end subroutine sulfnucrate


!----------------------------------------------------------------------
!-----------------------------------------------------------------------

subroutine binary_nuc_vehk2002( carma, temp, rh, h2so4,   &
      nucrate_cgs, mass_cluster_dry, radius_cluster )
!
! calculates binary nucleation rate and critical cluster size
! using the parameterization in
!   Vehkamäki, H., M. Kulmala, I. Napari, K.E.J. Lehtinen,
!        C. Timmreck, M. Noppel and A. Laaksonen, 2002,
!        An improved parameterization for sulfuric acid-water nucleation
!        rates for tropospheric and stratospheric conditions,
!        J. Geophys. Res., 107, 4622, doi:10.1029/2002jd002184
!
  use carma_precision_mod
  use carma_enums_mod
  use carma_constants_mod
  use carma_types_mod
  use carmastate_mod
  use carma_mod
  use sulfate_utils

  implicit none

! subr arguments (in)
  type(carma_type), intent(in)         :: carma          !! the carma object
  real(kind=f), intent(in) :: temp              ! temperature (k)
  real(kind=f), intent(in) :: rh                ! relative humidity (0-1)
  real(kind=f), intent(in) :: h2so4             ! concentration of h2so4 (molecules cm-3)

! subr arguments (out)
  real(kind=f), intent(out) :: nucrate_cgs      ! binary nucleation rate, j (# cm-3 s-1)
  real(kind=f), intent(out) :: mass_cluster_dry ! the mass of cluster (g)
  real(kind=f), intent(out) :: radius_cluster   ! the radius of cluster (cm)

! local variables
  real(kind=f) :: crit_x
  real(kind=f) :: acoe, bcoe, ccoe, dcoe, ecoe, fcoe, gcoe, hcoe, icoe, jcoe
  real(kind=f) :: tmpa, tmpb
  real(kind=f) :: cnum_h2so4       ! number of h2so4 molecules in the critical nucleus
  real(kind=f) :: cnum_tot         ! total number of molecules in the critical nucleus

! executable

! calc sulfuric acid mole fraction in critical cluster
  crit_x = 0.740997_f - 0.00266379_f * temp   &
         - 0.00349998_f * log (h2so4)   &
         + 0.0000504022_f * temp * log (h2so4)   &
         + 0.00201048_f * log (rh)   &
         - 0.000183289_f * temp * log (rh)   &
         + 0.00157407_f * (log (rh)) ** 2.0_f   &
         - 0.0000179059_f * temp * (log (rh)) ** 2.0_f   &
         + 0.000184403_f * (log (rh)) ** 3.0_f   &
         - 1.50345e-6_f * temp * (log (rh)) ** 3.0_f

! calc nucleation rate
  acoe    = 0.14309_f+2.21956_f*temp   &
          - 0.0273911_f * temp**2.0_f   &
          + 0.0000722811_f * temp**3.0_f + 5.91822_f/crit_x

  bcoe    = 0.117489_f + 0.462532_f *temp   &
          - 0.0118059_f * temp**2.0_f   &
          + 0.0000404196_f * temp**3.0_f + 15.7963_f/crit_x

  ccoe    = -0.215554_f-0.0810269_f * temp   &
          + 0.00143581_f * temp**2.0_f   &
          - 4.7758e-6_f * temp**3.0_f   &
          - 2.91297_f/crit_x

  dcoe    = -3.58856_f+0.049508_f * temp   &
          - 0.00021382_f * temp**2.0_f   &
          + 3.10801e-7_f * temp**3.0_f   &
          - 0.0293333_f/crit_x

  ecoe    = 1.14598_f - 0.600796_f * temp   &
          + 0.00864245_f * temp**2.0_f   &
          - 0.0000228947_f * temp**3.0_f   &
          - 8.44985_f/crit_x

  fcoe    = 2.15855_f + 0.0808121_f * temp   &
          -0.000407382_f * temp**2.0_f   &
          -4.01957e-7_f * temp**3.0_f   &
          + 0.721326_f/crit_x

  gcoe    = 1.6241_f - 0.0160106_f * temp   &
          + 0.0000377124_f * temp**2.0_f   &
          + 3.21794e-8_f * temp**3.0_f   &
          - 0.0113255_f/crit_x

  hcoe    = 9.71682_f - 0.115048_f * temp   &
          + 0.000157098_f * temp**2.0_f   &
          + 4.00914e-7_f * temp**3.0_f   &
          + 0.71186_f/crit_x

  icoe    = -1.05611_f + 0.00903378_f * temp   &
          - 0.0000198417_f * temp**2.0_f   &
          + 2.46048e-8_f  * temp**3.0_f   &
          - 0.0579087_f/crit_x

  jcoe    = -0.148712_f + 0.00283508_f * temp   &
          - 9.24619e-6_f  * temp**2.0_f   &
          + 5.00427e-9_f * temp**3.0_f   &
          - 0.0127081_f/crit_x

  tmpa     =     (   &
            acoe   &
          + bcoe * log (rh)   &
          + ccoe * ( log (rh))**2.0_f   &
          + dcoe * ( log (rh))**3.0_f   &
          + ecoe * log (h2so4)   &
          + fcoe * (log (rh)) * (log (h2so4))   &
          + gcoe * ((log (rh) ) **2.0_f)   &
                 * (log (h2so4))   &
          + hcoe * (log (h2so4)) **2.0_f   &
          + icoe * log (rh)   &
                 * ((log (h2so4)) **2.0_f)   &
          + jcoe * (log (h2so4)) **3.0_f   &
          )

  tmpa = min( tmpa, log(1.0e38_f) )
  nucrate_cgs = exp ( tmpa )

! calc number of molecules in critical cluster
  acoe    = -0.00295413_f - 0.0976834_f*temp   &
          + 0.00102485_f * temp**2.0_f   &
          - 2.18646e-6_f * temp**3.0_f - 0.101717_f/crit_x

!  write(LUNOPRT,*)'291 acoe=',acoe

  bcoe    = -0.00205064_f - 0.00758504_f*temp   &
          + 0.000192654_f * temp**2.0_f   &
          - 6.7043e-7_f * temp**3.0_f - 0.255774_f/crit_x

  ccoe    = +0.00322308_f + 0.000852637_f * temp   &
          - 0.0000154757_f * temp**2.0_f   &
          + 5.66661e-8_f * temp**3.0_f   &
          + 0.0338444_f/crit_x

  dcoe    = +0.0474323_f - 0.000625104_f * temp   &
          + 2.65066e-6_f * temp**2.0_f   &
          - 3.67471e-9_f * temp**3.0_f   &
          - 0.000267251_f/crit_x

  ecoe    = -0.0125211_f + 0.00580655_f * temp   &
          - 0.000101674_f * temp**2.0_f   &
          + 2.88195e-7_f * temp**3.0_f   &
          + 0.0942243_f/crit_x

  fcoe    = -0.038546_f - 0.000672316_f * temp   &
          + 2.60288e-6_f * temp**2.0_f   &
          + 1.19416e-8_f * temp**3.0_f   &
          - 0.00851515_f/crit_x

  gcoe    = -0.0183749_f + 0.000172072_f * temp   &
          - 3.71766e-7_f * temp**2.0_f   &
          - 5.14875e-10_f * temp**3.0_f   &
          + 0.00026866_f/crit_x

  hcoe    = -0.0619974_f + 0.000906958_f * temp   &
          - 9.11728e-7_f * temp**2.0_f   &
          - 5.36796e-9_f * temp**3.0_f   &
          - 0.00774234_f/crit_x

  icoe    = +0.0121827_f - 0.00010665_f * temp   &
          + 2.5346e-7_f * temp**2.0_f   &
          - 3.63519e-10_f * temp**3.0_f   &
          + 0.000610065_f/crit_x

  jcoe    = +0.000320184_f - 0.0000174762_f * temp   &
          + 6.06504e-8_f * temp**2.0_f   &
          - 1.4177e-11_f * temp**3.0_f   &
          + 0.000135751_f/crit_x

  cnum_tot = acoe + bcoe * log (rh)

  cnum_tot = exp (   &
            acoe   &
          + bcoe * log (rh)   &
          + ccoe * ( log (rh))**2.0_f   &
          + dcoe * ( log (rh))**3.0_f   &
          + ecoe * log (h2so4)   &
          + fcoe * (log (rh)) * (log (h2so4))   &
          + gcoe * ((log (rh) ) **2.0_f)   &
                 * (log (h2so4))   &
          + hcoe * (log (h2so4)) **2.0_f   &
          + icoe * log (rh)   &
                 * ((log (h2so4)) **2.0_f)   &
          + jcoe * (log (h2so4)) **3.0_f   &
          )

  cnum_h2so4 = cnum_tot * crit_x

!   calc radius (nm) of critical cluster
  radius_cluster = exp( -1.6524245_f + 0.42316402_f*crit_x   &
                        + 0.3346648_f*log(cnum_tot) )

  radius_cluster = radius_cluster * 1e-7_f ! nm -> cm

  mass_cluster_dry = cnum_h2so4 * gwtmol(igash2so4) / AVG ! cluster dry mass in g

  return
end subroutine binary_nuc_vehk2002

!----------------------------------------------------------------------
!-----------------------------------------------------------------------
subroutine binary_nuc_zhao1995( carma, cstate, temp, weight_percent, rh, h2so4, h2so4_cgs, h2o, h2o_cgs, beta1, &
          nucrate_cgs, mass_cluster_dry, radius_cluster, ftry, rc )
!!  Calculates particle production rates due to nucleation <rhompe>:
!!  binary homogeneous nucleation of sulfuric acid and water only
!!  Numerical method follows Zhao & Turco, JAS, V.26, No.5, 1995.
  use carma_precision_mod
  use carma_enums_mod
  use carma_constants_mod
  use carma_types_mod
  use carmastate_mod
  use carma_mod
  use sulfate_utils

  implicit none
! subr arguments (in)
  type(carma_type), intent(in)         :: carma          !! the carma object
  type(carmastate_type), intent(inout) :: cstate         !! the carma state object
  real(kind=f), intent(in) :: temp              ! temperature (k)
  real(kind=f), intent(in) :: weight_percent             ! weight percent H2SO4 (0-100)
  real(kind=f), intent(in) :: rh                ! relative humidity of water wrt liquid water (0-1)
  real(kind=f), intent(in) :: h2so4             ! concentration of H2SO4 (molecules cm-3)
  real(kind=f), intent(in) :: h2o               ! concentration of H2SO4 (molecules cm-3)
  real(kind=f), intent(in) :: h2so4_cgs         ! H2SO4 densities in g/cm3
  real(kind=f), intent(in) :: h2o_cgs           ! H2O densities in g/cm3
  real(kind=f), intent(in) :: beta1

! subr arguments (out)
  real(kind=f), intent(out) :: nucrate_cgs      ! binary nucleation rate, j (# cm-3 s-1)
  real(kind=f), intent(out) :: radius_cluster   ! the radius of cluster (cm)
  real(kind=f), intent(out) :: mass_cluster_dry ! dry mass of cluster (g)
  real(kind=f), intent(out) :: ftry
  integer, intent(inout)    :: rc              !! return code, negative indicates failure

  !  Local declarations
  integer           :: i, ibin, ie
  real(kind=f)      :: dens(46)
  real(kind=f)      :: pa(46)
  real(kind=f)      :: pb(46)
  real(kind=f)      :: c1(46)
  real(kind=f)      :: c2(46)
  real(kind=f)      :: fct(46)
  real(kind=f)      :: wtmolr         ! molecular weight ration of H2SO4/H2O
  real(kind=f)      :: h2oln          ! H2O ambient vapor pressures [dynes/cm2]
  real(kind=f)      :: h2so4ln        ! H2SO4 ambient vapor pressures [dynes/cm2]
  real(kind=f)      :: SA             ! total surface area of pre-existing wet particles
  real(kind=f)      :: SAbin          ! bin surface area of pre-existing wet particles
  real(kind=f)      :: cw
  real(kind=f)      :: dw
  real(kind=f)      :: wvp            ! water eq.vp over solution
  real(kind=f)      :: wvpln
  real(kind=f)      :: t0_kulm
  real(kind=f)      :: seqln
  real(kind=f)      :: t_crit_kulm
  real(kind=f)      :: factor_kulm
  real(kind=f)      :: dw1, dw2
  real(kind=f)      :: dens1
  real(kind=f)      :: dens11
  real(kind=f)      :: dens12
  real(kind=f)      :: xfrac
  real(kind=f)      :: wstar
  real(kind=f)      :: dstar
  real(kind=f)      :: rhln
  real(kind=f)      :: raln
  real(kind=f)      :: wfstar
  real(kind=f)      :: sigma
  real(kind=f)      :: ystar
  real(kind=f)      :: r2
  real(kind=f)      :: gstar
  real(kind=f)      :: rpr
  real(kind=f)      :: rpre
  real(kind=f)      :: fracmol
  real(kind=f)      :: zphi
  real(kind=f)      :: zeld
  real(kind=f)      :: cfac
  real(kind=f)      :: ahom
  real(kind=f)      :: exhom
  real(kind=f)      :: frac_h2so4
  real(kind=f)      :: rhomlim
  real(kind=f)      :: dnpot(46), dnwf(46)
  real(kind=f)      :: rho_H2SO4_wet

  radius_cluster = -1._f

  !  Parameterized fit developed by Mike Mills in 1994 to the partial molal
  !  Gibbs energies (F2|o-F2) vs. weight percent H2SO4 table in Giauque et al.,
  !  J. Am. Chem. Soc, 82, 62-70, 1960.  The parameterization gives excellent
  !  agreement.  Ayers (GRL, 7, 433-436, 1980) refers to F2|o-F2 as mu - mu_0
  !  (chemical potential).  This parameterization may be replaced by a lookup
  !  table, as was done ultimately in the Garcia-Solomon sulfate code.
  do i = 1, 46
    dnpot(i) = 4.184_f * (23624.8_f - 1.14208e8_f / ((dnwtp(i) - 105.318_f)**2 + 4798.69_f))
    dnwf(i) = dnwtp(i) / 100._f
  end do

  ! Molecular weight ratio of H2SO4 / H2O:
  wtmolr = gwtmol(igash2so4) / gwtmol(igash2o)

  ! Compute ln of H2O and H2SO4 ambient vapor pressures [dynes/cm2]
  h2oln   = log(h2o_cgs   * (RGAS / gwtmol(igash2o))   * temp)
  h2so4ln = log(h2so4_cgs * (RGAS / gwtmol(igash2so4)) * temp)

  ! loop through wt pcts and calculate vp/composition for each
  do i = 1, 46
    dens(i) = dnc0(i) + dnc1(i) * temp

    ! Calc. water eq.vp over solution using (Lin & Tabazadeh eqn 5, JGR, 2001)
    cw = 22.7490_f + 0.0424817_f * dnwtp(i) - 0.0567432_f * dnwtp(i)**0.5_f - 0.000621533_f * dnwtp(i)**2
    dw = -5850.24_f + 21.9744_f * dnwtp(i) - 44.5210_f * dnwtp(i)**0.5_f - 0.384362_f * dnwtp(i)**2

    ! pH20 | eq[mb]
    wvp   = exp(cw + dw / temp)

    ! Ln(pH2O | eq [dynes/cm2])
    wvpln = log(wvp * 1013250._f / 1013.25_f)

    ! Save the water eq.vp over solution at each wt pct into this array:
    !
    ! Ln(pH2O/pH2O|eq) with both terms in dynes/cm2
    pb(i) = h2oln - wvpln

    ! Calc. sulfuric acid eq.vp over solution using (Ayers et. al., GRL, V.7, No.6, June 1980)
    !
    ! T0 set in the low end of the Ayers measurement range (338-445K)
    t0_kulm = 340._f
    seqln   = -10156._f / t0_kulm + 16.259_f

    ! Now calc. Kulmala correction (J. CHEM. PHYS. V.93, No.1, 1 July 1990)
    !
    ! Critical temperature = 1.5 * Boiling point
    t_crit_kulm = 905._f
    factor_kulm = -1._f / temp + 1._f / t0_kulm + 0.38_f / (t_crit_kulm - t0_kulm) * &
      (1.0_f + log(t0_kulm / temp) - t0_kulm / temp)

    ! For pure sulfuric acid
    seqln = seqln + 10156._f * factor_kulm

    ! Now adjust vp based on weight % composition using parameterization of Giauque 1960
    !
    ! Adjust for weight percent composition
    seqln = seqln - dnpot(i) / (8.3143_f * temp)

    ! Convert atmospheres => dynes/cm2
    seqln = seqln + log(1013250._f)

    ! Save the sulfuric acid eq.vp over solution at each wt pct into this array:
    !
    ! Ln(pH2SO4/pH2SO4|eq) with both terms in dynes/cm2
    pa(i) = h2so4ln - seqln

    ! Create 2-component solutions of varying composition c1 and c2
    c1(i) = pa(i) - pb(i) * wtmolr
    c2(i) = pa(i) * dnwf(i) + pb(i) * (1._f - dnwf(i)) * wtmolr
  end do  ! end of loop through weight percents

  ! Now loop through until we find the c1+c2 combination with minimum Gibbs free energy
  dw2     = dnwtp(46) - dnwtp(45)
  dens1   = (dens(46) - dens(45)) / dw2
  fct(46) = c1(46) + c2(46) * 100._f * dens1 / dens(46)
  dens12 = dens1

  do i = 45, 2, -1
    dw1    = dw2
    dens11 = dens12
    dw2    = dnwtp(i) - dnwtp(i-1)
    dens12 = (dens(i) - dens(i-1)) / dw2
    dens1  = (dens11 * dw2 + dens12 * dw1) / (dw1 + dw2)

    fct(i) = c1(i) + c2(i) * 100._f * dens1 / dens(i)

    ! Find saddle where fct(i)<0<fct(i+1)
    if (fct(i) * fct(i+1) <= 0._f) exit
  end do

  if (i == 1) then
    dens1  = (dens(2) - dens(1)) / (dnwtp(2) - dnwtp(1))
    fct(1) = c1(1) + c2(1) * 100._f * dens1 / dens(1)
  end if

  ! Possibility 1: loop finds no saddle, so no nucleation occurs:
  if (fct(i) * fct(i+1) > 0._f) then
    nucrate_cgs = 0.0_f
    radius_cluster = 0.0_f
    mass_cluster_dry = 0.0_f
    ftry = 0.0_f
    return

  ! Possibility 2: loop crossed the saddle; interpolate to find exact value:
  else if (fct(i) * fct(i+1) < 0._f) then
    xfrac = fct(i+1) / (fct(i+1) - fct(i))
    wstar = dnwtp(i+1) * (1.0_f - xfrac) + dnwtp(i) * xfrac ! critical weight percent
    dstar = dens(i+1)  * (1.0_f - xfrac) + dens(i)  * xfrac
    rhln  = pb(i+1) * (1.0_f - xfrac) + pb(i) * xfrac
    raln  = pa(i+1) * (1.0_f - xfrac) + pa(i) * xfrac

  ! Possibility 3: loop found the saddle point exactly
  else
    dstar = dens(i)

    ! critical weight percent
    wstar = dnwtp(i)
    rhln  = pb(i)
    raln  = pa(i)
  end if

  ! Critical weight fraction
  wfstar = wstar / 100._f

  if ((wfstar < 0._f) .or. (wfstar > 1._f)) then
    write(LUNOPRT,*)'sulfnuc: wstar out of bounds!'
    rc = RC_ERROR
    return
  end if

  ! Critical surface tension  [erg/cm2]
  sigma = sulfate_surf_tens(carma, wstar, temp, rc)

  ! Critical Y (eqn 13 in Zhao & Turco 1993) [erg/cm3]
  ystar = dstar * RGAS * temp * (wfstar / gwtmol(igash2so4) &
      * raln + (1._f - wfstar) / gwtmol(igash2o) * rhln)
  if (ystar < 1.e-20_f) then
    nucrate_cgs = 0.0_f
    radius_cluster = 0.0_f
    mass_cluster_dry = 0.0_f
    ftry = 0.0_f
    return
  end if

  ! Critical cluster radius [cm]
  radius_cluster = 2._f * sigma / ystar
  radius_cluster = max(radius_cluster, 0.0_f)
  r2    = radius_cluster * radius_cluster

  ! Critical Gibbs free energy [erg]
  gstar = (4._f * PI / 3._f) * r2 * sigma

  ! RPR[molecules/s] = 4Pi * R2[cm2] * H2O[molecules/cm3] * Beta[cm/s]
  rpr = 4._f * PI * r2 * h2o * beta1

  ! RPRE[/cm3/s] = RPR[/s] * H2SO4[/cm3]; first part of Zhao & Turco eqn 16
  rpre = rpr * h2so4

  ! Zeldovitch non-equilibrium correction factor [unitless]
  ! Jaecker-Voirol & Mirabel, 1988 (not considered in Zhao & Turco)
  fracmol = 1._f /(1._f + wtmolr * (1._f - wfstar) / wfstar)
  zphi    = atan(fracmol)
  zeld    = 0.25_f / (sin(zphi))**2

  ! Empirical correction factor:
  cfac = 0.0_f

  ! Gstar exponential term in Zhao & Turco eqn 16 [unitless]
  ftry = (-gstar / BK / temp)
  ahom = ftry + cfac
  if (ahom .lt. -500._f) then
    exhom=0.0_f
  else
    exhom = exp(min(ahom, 28.0_f))
  endif

  !   Calculate mass of critical nucleus
  rho_H2SO4_wet = sulfate_density(carma, weight_percent, temp, rc)
  mass_cluster_dry = (4._f * PI / 3._f) * rho_H2SO4_wet * r2 * radius_cluster

  ! Calculate dry mass of critical nucleus
  mass_cluster_dry = mass_cluster_dry * wfstar

  ! Calculate the nucleation rate [#/cm3/s], Zhao & Turco eqn 16.
  nucrate_cgs = rpre * zeld * exhom

  return
end subroutine binary_nuc_zhao1995
