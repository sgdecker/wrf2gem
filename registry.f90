module registry
  use wrf2gem_subs, only: read_values, mem_error, freq_divisible,  &
                          open_wrf_file, close_wrf_file
  use gempak,       only: write_gempak
  use diagnostics,  only: vint, hydro_interp, mix_ratio_to_spec_hum,  &
                          temp_from_theta_p, surface_pres, density, sea_pres, &
                          lifted_index, mean_rh, unstag_hght, cape,  &
                          sfcape => sbcape, mlcape, lcl_lfc, storm_motion,  &
                          storm_rel_hel, calcdbz, calcgust
  use diagnostics,  only: Grav
  implicit none
  save
  private
  public :: init_reg, output_var

  ! The currently programmed WRF variables are:
  ! 1  - HGT: Terrain height (m)
  ! 2  - T2: 2-m temperature (K)
  ! 3  - Q2: 2-m water vapor mixing ratio (kg/kg)
  ! 4  - ZNU: Eta coordinate on mass (half) levels
  ! 5  - PB: Base state pressure (Pa)
  ! 6  - P: Perturbation pressure (Pa)
  ! 7  - T: Perterbation potential temperature (K)
  ! 8  - U: U component of wind (m/s)
  ! 9  - V: V component of wind (m/s)
  ! 10 - W: W component of wind (m/s)
  ! 11 - U10: U at 10 m (m/s)
  ! 12 - V10: V at 10 m (m/s)
  ! 13 - QVAPOR: Water vapor mixing ratio (kg/kg)
  ! 14 - QCLOUD: Cloud water mixing ratio (kg/kg)
  ! 15 - QICE: Ice mixing ratio (kg/kg)
  ! 16 - RAINC: Accumulated cumulus precip (mm)
  ! 17 - RAINNC: Accumulated grid-scale precip (mm)
  ! 18 - PH: Perturbation geopotential (m^2/s^2)
  ! 19 - PHB: Base-state geopotential (m^2/s^2)
  ! 20 - ZNW: Eta coordinate on w (full) levels
  ! 21 - LU_INDEX: Land use category
  ! 22 - MU: Perturbation dry air mass in column (Pa)
  ! 23 - MUB: Base-state dry air mass in column (Pa)
  ! 24 - TH2: 2-m potential temperature (K)
  ! 25 - QRAIN: Rain water mixing ratio (kg/kg)
  ! 26 - QSNOW: Snow mixing ratio (kg/kg)
  ! 27 - QGRAUP: Graupel mixing ratio (kg/kg)
  ! 28 - PBLH: PBL height (m)
  ! 29 - SST: Sea surface temperature (K)
  ! 30 - TSK: Surface skin temperature (K)
  ! 31 - SWDOWN: Downward short wave flux at ground surface (W/m^2)
  ! 32 - GLW: Downward long wave flux at ground surface (W/m^2)
  ! 33 - GRDFLX: Ground heat flux (W/m^2)
  ! 34 - HFX: Upward heat flux at the surface (W/m^2)
  ! 35 - QFX: Upward moisture flux at the surface (W/m^2)
  ! 36 - LH: Latent heat flux at the surface (W/m^2)
  ! 37 - PBLH: PBL Height (m)

  ! The currently programmed diagnostics are:
  ! 101 - Surface pressure (hPa)
  ! 102 - Pressure at 2 m (hPa)
  ! 103 - Sea level pressure (hPa)
  ! 104 - Convective precipitation accumulated over one hour (mm)
  ! 105 - Total precipitation accumulated over one hour (mm)
  ! 106 - Total precipitation accumulated over three hours (mm)
  ! 107 - Total precipitation accumulated over six hours (mm)
  ! 108 - Total accumulated precipitation (mm)
  ! 109 - Pressure (hPa)
  ! 110 - Potential temperature (K)
  ! 111 - Unstaggered U component of wind (m/s)
  ! 112 - Unstaggered V component of wind (m/s)
  ! 113 - Cloud water + Ice mixing ratio (kg/kg)
  ! 114 - Unstaggered W component of wind (m/s)
  ! 115 - 2-m temperature including at initial time (K)
  ! 116 - 2-m water vapor mixing ratio including at initial time (kg/kg)
  ! 117 - U at 10 m including at initial time (m/s)
  ! 118 - V at 10 m including at initial time (m/s)
  ! 119 - Lifted Index using lowest 4 levels (K)
  ! 120 - Mean relative humidity (%) in 850-500-hPa layer
  ! 121 - Geopotential height (m) on w (full) levels
  ! 122 - 2-m specific humidity (kg/kg)
  ! 123 - Specific humidity (kg/kg)
  ! 124 - Temperature (K)
  ! 125 - Geopotential height (m) on half levels
  ! 126 - Dry air mass in column (Pa)
  ! 127 - Most unstable CAPE w/ virt. temp. correction (J/kg)
  ! 128 - Lifted parcel level of the most unstable parcel (m)
  ! 129 - Surface-based CAPE w/ virt. temp. correction (J/kg)
  ! 130 - Surface-based CIN w/ virt. temp. correction (J/kg)
  ! 131 - Mixed-layer CAPE w/ virt. temp. correction (J/kg)
  ! 132 - Mixed-layer CIN w/ virt. temp. correction (J/kg)
  ! 133 - Mixed-layer LCL (m)
  ! 134 - Mixed-layer LFC (m)
  ! 135 - U component of Bunkers storm motion (m/s)
  ! 136 - V component of Bunkers storm motion (m/s)
  ! 137 - 0-1-km storm relative helicity (m2/s2)
  ! 138 - Column integrated precipitable water (cm)
  ! 139 - Reflectivity (dBZ)
  ! 140 - Composite reflectivity (dBZ)
  ! 141 - Surface Wind Gusts (m/s)
  ! 142 - Height above ground on half levels (m)
  ! 143 - Reflectivity assuming no graupel (dBZ)

  ! The currently programmed diagnostics interpolated to pressure are:
  ! 501 - Temperature (K)
  ! 502 - U component of wind (m/s)
  ! 503 - V component of wind (m/s)
  ! 504 - W component of wind (m/s)
  ! 505 - Geopotential height (m)
  ! 506 - Specific humidity (kg/kg)

  ! Encode GEMPAK vertical coordinate IDs
  integer, parameter :: HAGL = 1279738184, None = 0, Sgma = 4, Pres = 1

  ! Set up enumeration for the available variables
  integer, parameter :: LastWRF = 37, FirstDiag = 101, LastDiag = 143,  &
                        FirstPresDiag = 501, LastPresDiag = 506
  integer, parameter :: TerHght   = 1,  SfcTemp   = 2,  SfcMixR   = 3,   &
                        EtaMass   = 4,  PBase     = 5,  PPert     = 6,   &
                        ThetaPert = 7,  UStag     = 8,  VStag     = 9,   &
                        Wstag     = 10, SfcU      = 11, SfcV      = 12,  &
                        MixR      = 13, MixRCloud = 14, MixRIce   = 15,  &
                        ConvPre   = 16, StabPre   = 17, GeopPert  = 18,  &
                        GeopBase  = 19, EtaW      = 20, LandUse   = 21,  &
                        MuPert    = 22, MuBase    = 23, SfcTheta  = 24,  &
                        MixRRain  = 25, MixRSnow  = 26, MixRGraup = 27,  &
                        PBLHgt    = 28, SeaTemp   = 29, SkinTemp  = 30,  &
                        SWGrd     = 31, LWGrd     = 32, GrdFx     = 33,  &
                        HeatFx    = 34, MoisFx    = 35, LaHeat    = 36,  &
                        PBLHght   = 37

  integer, parameter :: SfcP      = 101, SfcP2m    = 102, SeaLevP   = 103,  &
                        ConvPre1h = 104, TotPre1h  = 105, TotPre3h  = 106,  &
                        TotPre6h  = 107, TotPre    = 108, Pressure  = 109,  &
                        Theta     = 110, UUnstag   = 111, VUnstag   = 112,  &
                        CWtr      = 113, WUnstag   = 114, SfcTempF0 = 115,  &
                        SfcMixRF0 = 116, SfcUF0    = 117, SfcVF0    = 118,  &
                        LiftIdx   = 119, MeanRH    = 120, GeopHghtS = 121,  &
                        SfcSpHum  = 122, SpHum     = 123, TmpK      = 124,  &
                        GeopHght  = 125, DryAir    = 126, MUCAPE    = 127,  &
                        MULPL     = 128, SBCAPE    = 129, SBCIN     = 130,  &
                        MixCAPE   = 131, MixCIN    = 132, MixLCL    = 133,  &
                        MixLFC    = 134, UStorm    = 135, VStorm    = 136,  &
                        StormRHel = 137, PreWater  = 138, Refl      = 139,  &
                        CompRefl  = 140, WindGust  = 141, ZAbvGrnd  = 142
  integer, parameter :: ReflNoG   = 143

  integer, parameter :: TmpKPres     = 501, UUnstagPres  = 502,  &
                        VUnstagPres  = 503, WUnstagPres  = 504,  &
                        GeopHghtPres = 505, SpHumPres    = 506

  ! This structure contains information about the various WRF variables.
  type wrf_var
     integer :: dimType,  &       ! 0-z,t; 1-x,y,t; 2-x,y,z,t; 3-xs,y,z,t;
                                  ! 4-x,ys,z,t; 5-x,y,zs,t
                lev1, lev2,  &    ! GEMPAK level parameters. Relevant only for
                                  ! dimType 1
                vertCordID        ! GEMPAK vertical coordinate ID
     character(len=8) :: wrfName  ! WRF variable name
     character(len=4) :: gemName  ! GEMPAK variable name
     logical          :: const    ! T - Fixed with time; F - Varies with time
  end type wrf_var

  ! This is an analogous structure for diagnosed variables.
  type diag_var
     integer           :: dimType, lev1, lev2, vertCordID, diagID
     character(len=12) :: gemName
  end type diag_var

  ! And finally a structure for diagnosed variables interpolated to pressure.
  type pres_diag_var
     integer           :: interpType,  &  ! 1-linear interp; 2-log interp
                          diagID
     character(len=12) :: gemName
  end type pres_diag_var

  type(wrf_var),  dimension(LastWRF)                         :: wrfVar
  type(diag_var), dimension(FirstDiag:LastDiag)              :: diagVar
  type(pres_diag_var), dimension(FirstPresDiag:LastPresDiag) :: presDiagVar
  integer :: timesPer12hr, pBot, pTop, dp, np
  
contains  ! ===================================================================

  ! init_reg initializes the registry, which contains the vital parameters for
  ! each variable or diagnostic, with the number of outputs per 12 hours (freq)
  ! and pressure interpolation information (pb, pt, deltap) provided as input.
  subroutine init_reg(freq, pb, pt, deltap)
    integer, intent(in) :: freq, pb, pt, deltap

    wrfVar(TerHght)   = wrf_var(1, 0, -1, None, "HGT", "HGHT", .true.)
    wrfVar(SfcTemp)   = wrf_var(1, 2, -1, HAGL, "T2", "TMPK", .false.)
    wrfVar(SfcMixR)   = wrf_var(1, 2, -1, HAGL, "Q2", "MIXR", .false.)
    wrfVar(EtaMass)   = wrf_var(0, -1, -1, -1, "ZNU", "SGMA", .true.)
    wrfVar(PBase)     = wrf_var(2, -1, -1, Sgma, "PB", "PBAR", .true.)
    wrfVar(PPert)     = wrf_var(2, -1, -1, Sgma, "P", "PP", .false.)
    wrfVar(ThetaPert) = wrf_var(2, -1, -1, Sgma, "T", "TMPO", .false.)
    wrfVar(UStag)     = wrf_var(3, -1, -1, Sgma, "U", "USTG", .false.)
    wrfVar(VStag)     = wrf_var(4, -1, -1, Sgma, "V", "VSTG", .false.)
    wrfVar(WStag)     = wrf_var(5, -1, -1, Sgma, "W", "WSTG", .false.)
    wrfVar(SfcU)      = wrf_var(1, 10, -1, HAGL, "U10", "UREL", .false.)
    wrfVar(SfcV)      = wrf_var(1, 10, -1, HAGL, "V10", "VREL", .false.)
    wrfVar(MixR)      = wrf_var(2, -1, -1, Sgma, "QVAPOR", "MIXR", .false.)
    wrfVar(MixRCloud) = wrf_var(2, -1, -1, Sgma, "QCLOUD", "WC", .false.)
    wrfVar(MixRIce)   = wrf_var(2, -1, -1, Sgma, "QICE", "WI", .false.)
    wrfVar(ConvPre)   = wrf_var(1, 0, -1, None, "RAINC", "CTOT", .false.)
    wrfVar(StabPre)   = wrf_var(1, 0, -1, None, "RAINNC", "ETOT", .false.)
    wrfVar(GeopPert)  = wrf_var(5, -1, -1, Sgma, "PH", "HP", .false.)
    wrfVar(GeopBase)  = wrf_var(5, -1, -1, Sgma, "PHB", "HBAR", .true.)
    wrfVar(EtaW)      = wrf_var(0, -1, -1, -1, "ZNW", "SGMA", .true.)
    wrfVar(LandUse)   = wrf_var(1, 0, -1, None, "LU_INDEX", "LUC", .true.)
    wrfVar(MuPert)    = wrf_var(1, 0, -1, None, "MU", "MUP", .false.)
    wrfVar(MuBase)    = wrf_var(1, 0, -1, None, "MUB", "MUB", .true.)
    wrfVar(SfcTheta)  = wrf_var(1, 2, -1, HAGL, "TH2", "THTA", .false.)
    wrfVar(MixRRain)  = wrf_var(2, -1, -1, Sgma, "QRAIN", "WR", .false.)
    wrfVar(MixRSnow)  = wrf_var(2, -1, -1, Sgma, "QSNOW", "WS", .false.)
    wrfVar(MixRGraup) = wrf_var(2, -1, -1, Sgma, "QGRAUP", "WG", .false.)
    wrfVar(PBLHgt)    = wrf_var(1, 0, -1, None, "PBLH", "PBLH", .false.)
    wrfVar(SeaTemp)   = wrf_var(1, 0, -1, None, "SST", "SST", .false.)
    wrfVar(SkinTemp)  = wrf_var(1, 0, -1, None, "TSK", "TSK", .false.)
    wrfVar(SWGrd)     = wrf_var(1, 0, -1, None, "SWDOWN", "SWD", .false.)
    wrfVar(LWGrd)     = wrf_var(1, 0, -1, None, "GLW", "LWD", .false.)
    wrfVar(GrdFx)     = wrf_var(1, 0, -1, None, "GRDFLX", "GRDF", .false.)
    wrfVar(HeatFx)    = wrf_var(1, 0, -1, None, "HFX", "HFX", .false.)
    wrfVar(MoisFx)    = wrf_var(1, 0, -1, None, "QFX", "QFX", .false.)
    wrfVar(LaHeat)    = wrf_var(1, 0, -1, None, "LH", "LH", .false.)
    wrfVar(PBLHght)   = wrf_var(1, 0, -1, None, "PBLH", "ZPBL", .false.)
    
    diagVar(SfcP)      = diag_var(1, 0, -1, None, SfcP, "PRES")
    diagVar(SfcP2m)    = diag_var(1, 2, -1, HAGL, SfcP2m, "PRES")
    diagVar(SeaLevP)   = diag_var(1, 0, -1, None, SeaLevP, "PMSL")
    diagVar(ConvPre1h) = diag_var(1, 0, -1, None, ConvPre1h, "C01M")
    diagVar(TotPre1h)  = diag_var(1, 0, -1, None, TotPre1h, "P01M")
    diagVar(TotPre3h)  = diag_var(1, 0, -1, None, TotPre3h, "P03M")
    diagVar(TotPre6h)  = diag_var(1, 0, -1, None, TotPre6h, "P06M")
    diagVar(TotPre)    = diag_var(1, 0, -1, None, TotPre, "PTOT")
    diagVar(Pressure)  = diag_var(2, -1, -1, Sgma, Pressure, "PRES")
    diagVar(Theta)     = diag_var(2, -1, -1, Sgma, Theta, "THTA")
    diagVar(UUnstag)   = diag_var(2, -1, -1, Sgma, UUnstag, "UREL")
    diagVar(VUnstag)   = diag_var(2, -1, -1, Sgma, VUnstag, "VREL")
    diagVar(CWtr)      = diag_var(2, -1, -1, Sgma, CWtr, "CWTR")
    diagVar(WUnstag)   = diag_var(2, -1, -1, Sgma, WUnstag, "W")
    diagVar(SfcTempF0) = diag_var(1, 2, -1, HAGL, SfcTempF0, "TMPK")
    diagVar(SfcMixRF0) = diag_var(1, 2, -1, HAGL, SfcMixRF0, "MIXR")
    diagVar(SfcUF0)    = diag_var(1, 10, -1, HAGL, SfcUF0, "UREL")
    diagVar(SfcVF0)    = diag_var(1, 10, -1, HAGL, SfcVF0, "VREL")
    diagVar(LiftIdx)   = diag_var(1, 0, -1, None, LiftIdx, "LFT4")
    diagVar(MeanRH)    = diag_var(1, 850, 500, Pres, MeanRH, "RELH")
    diagVar(GeopHghtS) = diag_var(5, -1, -1, Sgma, GeopHghtS, "HGHT")
    diagVar(SfcSpHum)  = diag_var(1, 2, -1, HAGL, SfcSpHum, "SPFH")
    diagVar(SpHum)     = diag_var(2, -1, -1, Sgma, SpHum, "SPFH")
    diagVar(TmpK)      = diag_var(2, -1, -1, Sgma, TmpK, "TMPK")
    diagVar(GeopHght)  = diag_var(2, -1, -1, Sgma, GeopHght, "HGHT")
    diagVar(DryAir)    = diag_var(1, 0, -1, None, DryAir, "MU")
    diagVar(MUCAPE)    = diag_var(1, 0, -1, None, MUCAPE, "MUCAPE")
    diagVar(MULPL)     = diag_var(1, 0, -1, None, MULPL, "LPL")
    diagVar(SBCAPE)    = diag_var(1, 0, -1, None, SBCAPE, "SBCAPE")
    diagVar(SBCIN)     = diag_var(1, 0, -1, None, SBCIN, "SBCIN")
    diagVar(MixCAPE)   = diag_var(1, 0, -1, None, MixCAPE, "MLCAPE")
    diagVar(MixCIN)    = diag_var(1, 0, -1, None, MixCIN, "MLCIN")
    diagVar(MixLCL)    = diag_var(1, 0, -1, None, MixLCL, "LCL")
    diagVar(MixLFC)    = diag_var(1, 0, -1, None, MixLFC, "LFC")
    diagVar(UStorm)    = diag_var(1, 0, -1, None, UStorm, "USTM")
    diagVar(VStorm)    = diag_var(1, 0, -1, None, VStorm, "VSTM")
    diagVar(StormRHel) = diag_var(1, 0, 1000, HAGL, StormRHel, "SRH")
    diagVar(PreWater)  = diag_var(1, 0, -1, None, PreWater, "PW")
    diagVar(Refl)      = diag_var(2, -1, -1, Sgma, Refl, "REFL")
    diagVar(CompRefl)  = diag_var(1, 0, -1, None, CompRefl, "CREF")
    diagVar(WindGust)  = diag_var(1, 0, -1, None, WindGust, "GUST")
    diagVar(ZAbvGrnd)  = diag_var(2, -1, -1, Sgma, ZAbvGrnd, "HAGL")
    diagVar(ReflNoG)   = diag_var(2, -1, -1, Sgma, ReflNoG, "REFL")

    presDiagVar(TmpKPres)     = pres_diag_var(1, TmpKPres, "TMPK")
    presDiagVar(UUnstagPres)  = pres_diag_var(1, UUnstagPres, "UREL")
    presDiagVar(VUnstagPres)  = pres_diag_var(1, VUnstagPres, "VREL")
    presDiagVar(WUnstagPres)  = pres_diag_var(1, WUnstagPres, "W")
    presDiagVar(GeopHghtPres) = pres_diag_var(2, GeopHghtPres, "HGHT")! 1 or 2?
    presDiagVar(SpHumPres)    = pres_diag_var(1, SpHumPres, "SPFH")

    timesPer12hr = freq

    pBot = pb
    pTop = pt
    dp = deltap
    np = (pBot - pTop) / dp + 1
  end subroutine init_reg

  ! ===========================================================================

  ! output_var writes a particular WRF variable or diagnostic to the GEMPAK
  ! file associated with gemid.  ncid contains the file handle for the WRF
  ! history file from which the WRF data is being read.  The WRF model
  ! dimensions are specified by nx, ny, nz.  var is the ID for the requested
  ! diagnostic.  ind is the current index of histFiles, so that histFiles(ind)
  ! is the pathname for the GEMPAK history file currently being read.  times 
  ! is an array of output times in GEMPAK form, which corresponds to all of the
  ! times in the history file.
  subroutine output_var(ncid, gemid, nx, ny, nz, var, ind, histFiles, times)
    integer, intent(in) :: ncid, gemid, nx, ny, nz, var, ind
    character(len=80), dimension(:), intent(in) :: histFiles
    character(len=15), dimension(:), intent(in) :: times

    character(len=*), parameter :: fmt1 = "(2x,a,i3,a)"

    write (*,fmt1) "Writing field ", var, " to GEMPAK..."

    ! Call the appropriate subroutine to produce the required output for GEMPAK
    select case (var)
    case (1:LastWRF)
       if (wrfVar(var)%dimType == 1 .or. wrfVar(var)%dimType == 2 .or.  &
            wrfVar(var)%dimType == 5) then
          call simple_out(ncid, gemid, nx, ny, nz, wrfVar(var), times)
       else
          print *, "It does not make sense to output this variable to GEMPAK. &
               &Skipping!"
       end if
    case (FirstDiag:LastDiag)
       if (diagVar(var)%dimType == 1 .or. diagVar(var)%dimType == 2 .or.  &
            diagVar(var)%dimType == 5) then
          call diag_out(ncid, gemid, nx, ny, nz, ind, histFiles,  &
               diagVar(var), times)
       else
          print *, "This diagnostic is not appropriate for output to GEMPAK."
          print *, "It's either staggered in x or y, or not a function of x &
               &and y."
          print *, "Skipping!"
       end if
    case (FirstPresDiag:LastPresDiag)
       call pres_diag_out(ncid, gemid, nx, ny, nz, presDiagVar(var), times)
    case default
       write (*,fmt1) "Variable ", var, " is unknown. Skipping!"
    end select
  end subroutine output_var

  ! ===========================================================================

  ! simple_out reads the variable var in from the WRF history file associated
  ! with ncid and writes it out to the GEMPAK file associated with gemid.
  ! Other inputs to the subroutine are the WRF model dimensions (nx,ny,nz) and
  ! the array of times in GEMPAK form to be output, which corresponds to all of
  ! the times in the history file.
  subroutine simple_out(ncid, gemid, nx, ny, nz, var, times)
    integer,                         intent(in) :: ncid, gemid, nx, ny, nz
    type(wrf_var),                   intent(in) :: var
    character(len=15), dimension(:), intent(in) :: times
     
    real, dimension(:,:,:,:), allocatable :: out4d
    real, dimension(:,:,:),   allocatable :: out3d
    real, dimension(:,:),     allocatable :: elvl
    integer                               :: kMax, eta, nt, status, t, k, lev
    character(len=20)                     :: gdat1
    character(len=12)                     :: parm

    ! Initialization
    select case (var%dimType)
    case (1)  ! x,y,t
       kMax = 0
    case (2)  ! x,y,z,t
       kMax = nz
       eta = EtaMass
    case (5)  ! x,y,zs,t
       kMax = nz + 1
       eta = EtaW
    case default
       ! GEMPAK only understands variables that depend on (unstaggered) x and y
       print *, "Invalid dimType in simple_out!"
       stop
    end select
    
    nt = size(times)
    parm = var%gemName
    
    ! Read variable, and, if necessary, the vertical levels it's at
    if (var%dimType == 1) then
       allocate (out3d(nx,ny,nt), stat=status)
       call mem_error(status, 1, "simple_out")
       call read_values(ncid, var%wrfName, out3d)
    else
       allocate (out4d(nx,ny,kMax,nt), elvl(kMax,nt), stat=status)
       call mem_error(status, 1, "simple_out")
       call read_values(ncid, var%wrfName, out4d)
       call read_values(ncid, wrfVar(eta)%wrfName, elvl)
    end if

    ! Write the variable to GEMPAK
    if (var%const) nt = 1  ! Output constant variables only at initial time
    do t = 1, nt
       gdat1 = times(t)
       if (var%dimType == 1) then
          call write_gempak(gemid, out3d(:,:,t), nx, ny, gdat1, var%lev1,  &
               var%lev2, var%vertCordID, parm)
       else
          do k = 1, kMax
             lev = nint(elvl(k,1) * 10000)
             call write_gempak(gemid, out4d(:,:,k,t), nx, ny, gdat1, lev,  &
                  var%lev2, var%vertCordID, parm)
          end do
       end if
    end do
        
    ! Clean up
    if (allocated(out3d)) then
       deallocate (out3d, stat=status)
    else
       deallocate (out4d, elvl, stat=status)
    end if
    call mem_error(status, 2, "simple_out")
  end subroutine simple_out
  
  ! ===========================================================================

  ! diag_out is a driver routine that produces the diagnostic var based on data
  ! from the WRF history file associated with ncid.  The diagnostic is written
  ! out to the GEMPAK file associated with gemid.  Other inputs to the
  ! subroutine are the WRF model dimensions (nx,ny,nz), the current index (ind)
  ! of histFiles, so that histFiles(ind) is the pathname for the GEMPAK history
  ! file currently being read, and the array of times in GEMPAK form to be
  ! output, which corresponds to all of the times in the history file.
  subroutine diag_out(ncid, gemid, nx, ny, nz, ind, histFiles, var, times)
    integer,                         intent(in) :: ncid, gemid, nx, ny, nz, ind
    character(len=80), dimension(:), intent(in) :: histFiles
    type(diag_var),                  intent(in) :: var
    character(len=15), dimension(:), intent(in) :: times
    
    real, dimension(:,:,:,:), allocatable :: out4d
    real, dimension(:,:,:),   allocatable :: out3d
    real, dimension(:,:),     allocatable :: elvl
    integer           :: kMax, eta, nt, status, prev, diff, t, k, lev
    character(len=20) :: gdat1
    logical           :: ok
    
    ! Initialization
    select case (var%dimType)
    case (1)  ! x,y,t
       kMax = 0
    case (2)  ! x,y,z,t
       kMax = nz
       eta = EtaMass
    case (5)  ! x,y,zs,t
       kMax = nz + 1
       eta = EtaW
    case default
       print *, "Invalid dimType in diag_out!"
       stop
    end select

    nt = size(times)
    ok = .true.
  
    if (var%dimType == 1) then
       allocate (out3d(nx,ny,nt), stat=status)
       call mem_error(status, 1, "diag_out")
    else
       allocate (out4d(nx,ny,kMax,nt), elvl(kMax,nt), stat=status)
       call mem_error(status, 1, "diag_out")
       call read_values(ncid, wrfVar(eta)%wrfName, elvl)
    end if

    ! Compute appropriate diagnostic
    select case (var%diagID)
    case (SfcP)
       out3d = diag_SfcP(ncid, nx, ny, nz, nt)
    case (SfcP2m)
       out3d = diag_SfcP2m(ncid, nx, ny, nz, nt)
    case (SeaLevP)
       out3d = diag_SeaLevP(ncid, nx, ny, nz, nt)
    case (ConvPre1h)
       if (.not. freq_divisible(timesPer12hr, 12, 1)) then
          print *, "Calculating 1-hr conv. precip makes no sense. Skipping..."
          ok = .false.
       else
          if (ind > 1) call open_wrf_file(histFiles(ind-1), prev)
          out3d = diag_ConvPre1h(ncid, nx, ny, nt, ind, prev)
          if (ind > 1) call close_wrf_file(prev)
       end if
    case (TotPre1h)
       if (.not. freq_divisible(timesPer12hr, 12, 1)) then
          print *, "Calculating 1-hr total precip makes no sense. Skipping..."
          ok = .false.
       else
          if (ind > 1) call open_wrf_file(histFiles(ind-1), prev)
          out3d = diag_TotPreHr(ncid, nx, ny, nt, 1, ind, prev)
          if (ind > 1) call close_wrf_file(prev)
       end if
    case (TotPre3h)
       if (.not. freq_divisible(timesPer12hr, 12, 3)) then
          print *, "Calculating 3-hr total precip makes no sense. Skipping..."
          ok = .false.
       else
          diff = timesPer12hr / 4
          if (ind > diff) call open_wrf_file(histFiles(ind-diff), prev)
          out3d = diag_TotPreHr(ncid, nx, ny, nt, 3, ind, prev)
          if (ind > diff) call close_wrf_file(prev)
       end if
    case (TotPre6h)
       if (.not. freq_divisible(timesPer12hr, 12, 6)) then
          print *, "Calculating 6-hr total precip makes no sense. Skipping..."
          ok = .false.
       else
          diff = timesPer12hr / 2
          if (ind > diff) call open_wrf_file(histFiles(ind-diff), prev)
          out3d = diag_TotPreHr(ncid, nx, ny, nt, 6, ind, prev)
          if (ind > diff) call close_wrf_file(prev)
       end if
    case (TotPre)
       out3d = diag_TotPre(ncid, nx, ny, nt)
    case (Pressure)
       out4d = diag_Pressure(ncid, nx, ny, nz, nt)
    case (Theta)
       out4d = diag_Theta(ncid, nx, ny, nz, nt)
    case (UUnstag)
       out4d = diag_UUnstag(ncid, nx, ny, nz, nt)
    case (VUnstag)
       out4d = diag_VUnstag(ncid, nx, ny, nz, nt)
    case (CWtr)
       out4d = diag_CWtr(ncid, nx, ny, nz, nt)
    case (WUnstag)
       out4d = diag_WUnstag(ncid, nx, ny, nz, nt)
    case (SfcTempF0)
       out3d = diag_SfcTempF0(ncid, nx, ny, nz, nt)
    case (SfcMixRF0)
       out3d = diag_SfcMixRF0(ncid, nx, ny, nz, nt)
    case (SfcUF0)
       out3d = diag_SfcUF0(ncid, nx, ny, nz, nt)
    case (SfcVF0)
       out3d = diag_SfcVF0(ncid, nx, ny, nz, nt)
    case (LiftIdx)
       out3d = diag_LiftIdx(ncid, nx, ny, nz, nt)
    case (MeanRH)
       out3d = diag_MeanRH(ncid, nx, ny, nz, nt)
    case (GeopHghtS)
       out4d = diag_GeopHghtS(ncid, nx, ny, nz, nt)
    case (SfcSpHum)
       out3d = diag_SfcSpHum(ncid, nx, ny, nz, nt)
    case (SpHum)
       out4d = diag_SpHum(ncid, nx, ny, nz, nt)
    case (TmpK)
       out4d = diag_TmpK(ncid, nx, ny, nz, nt)
    case (GeopHght)
       out4d = diag_GeopHght(ncid, nx, ny, nz, nt)
    case (DryAir)
       out3d = diag_DryAir(ncid, nx, ny, nt)
    case (MUCAPE)
       out3d = diag_MUCAPE(ncid, nx, ny, nz, nt)
    case (MULPL)
       out3d = diag_MULPL(ncid, nx, ny, nz, nt)
    case (SBCAPE)
       out3d = diag_SBCAPE(ncid, nx, ny, nz, nt)
    case (SBCIN)
       out3d = diag_SBCIN(ncid, nx, ny, nz, nt)
    case (MixCAPE)
       out3d = diag_MixCAPE(ncid, nx, ny, nz, nt)
    case (MixCIN)
       out3d = diag_MixCIN(ncid, nx, ny, nz, nt)
    case (MixLCL)
       out3d = diag_MixLCL(ncid, nx, ny, nz, nt)
    case (MixLFC)
       out3d = diag_MixLFC(ncid, nx, ny, nz, nt)
    case (UStorm)
       out3d = diag_UStorm(ncid, nx, ny, nz, nt)
    case (VStorm)
       out3d = diag_VStorm(ncid, nx, ny, nz, nt)
    case (StormRHel)
       out3d = diag_StormRHel(ncid, nx, ny, nz, nt)
    case (PreWater)
       out3d = diag_PreWater(ncid, nx, ny, nz, nt)
    case (Refl)
       out4d = diag_Refl(ncid, nx, ny, nz, nt)
    case (CompRefl)
       out3d = diag_CompRefl(ncid, nx, ny, nz, nt)
    case (WindGust)
       out3d = diag_WindGust(ncid, nx, ny, nz, nt)
    case (ZAbvGrnd)
       out4d = diag_ZAbvGrnd(ncid, nx, ny, nz, nt)
    case (ReflNoG)
       out4d = diag_Refl(ncid, nx, ny, nz, nt, .false.)
    case default
       write (*,"(a,i3,a)") "WARNING: Field ", var, " unknown in diag_out!"
    end select

    ! Write the diagnostic out to GEMPAK
    if (ok) then
       do t = 1, nt
          gdat1 = times(t)
          if (var%dimType == 1) then
             call write_gempak(gemid, out3d(:,:,t), nx, ny, gdat1, var%lev1,  &
                  var%lev2, var%vertCordID, var%gemName)
          else
             do k = 1, kMax
                lev = nint(elvl(k,1) * 10000)
                call write_gempak(gemid, out4d(:,:,k,t), nx, ny, gdat1, lev,  &
                     var%lev2, var%vertCordID, var%gemName)
             end do
          end if
       end do
    end if

    ! Clean up
    if (allocated(out3d)) then
       deallocate (out3d, stat=status)
    else
       deallocate (out4d, elvl, stat=status)
    end if
    call mem_error(status, 2, "diag_out")
  end subroutine diag_out

  ! ===========================================================================

  ! pres_diag_out is a driver routine that produces the diagnostic var
  ! interpolated to pressure surfaces based on data from the WRF history file
  ! associated with ncid.  The diagnostic is written out to the GEMPAK file
  ! associated with gemid.  Other inputs to the subroutine are the WRF model
  ! dimensions (nx,ny,nz) and the array of times in GEMPAK form to be output,
  ! which corresponds to all of the times in the history file.  np is obtained
  ! through use association.
  subroutine pres_diag_out(ncid, gemid, nx, ny, nz, var, times)
    integer,                         intent(in) :: ncid, gemid, nx, ny, nz
    type(pres_diag_var),             intent(in) :: var
    character(len=15), dimension(:), intent(in) :: times

    real, dimension(:,:,:,:), allocatable :: tk
    real, dimension(nx,ny,nz+1,size(times)) :: outNoInterp, p
    real, dimension(nx,ny,np,size(times))   :: outInterp
    integer                                 :: nt, t, j, i, status, k, lev
    character(len=20)                       :: gdat1
    
    nt = size(times)
    
    ! Compute appropriate diagnostic
    select case (var%diagID)
    case (TmpKPres)
       outNoInterp(:,:,1,:) = diag_SfcTempF0(ncid, nx, ny, nz, nt)
       outNoInterp(:,:,2:nz+1,:) = diag_TmpK(ncid, nx, ny, nz, nt)
    case (UUnstagPres)
       outNoInterp(:,:,1,:) = diag_SfcUF0(ncid, nx, ny, nz, nt)
       outNoInterp(:,:,2:nz+1,:) = diag_UUnstag(ncid, nx, ny, nz, nt)
    case (VUnstagPres)
       outNoInterp(:,:,1,:) = diag_SfcVF0(ncid, nx, ny, nz, nt)
       outNoInterp(:,:,2:nz+1,:) = diag_VUnstag(ncid, nx, ny, nz, nt)
    case (WUnstagPres)
       outNoInterp(:,:,1,:) = 0
       outNoInterp(:,:,2:nz+1,:) = diag_WUnstag(ncid, nx, ny, nz, nt)
    case (GeopHghtPres)
       call read_values(ncid, wrfVar(TerHght)%wrfName, outNoInterp(:,:,1,:))
       outNoInterp(:,:,2:nz+1,:) = diag_GeopHght(ncid, nx, ny, nz, nt)
    case (SpHumPres)
       outNoInterp(:,:,1,:) = diag_SfcSpHum(ncid, nx, ny, nz, nt)
       outNoInterp(:,:,2:nz+1,:) = diag_SpHum(ncid, nx, ny, nz, nt)
    case default
       write (*,"(a,i3,a)") "WARNING: Field ", var,  &
            " unknown in pres_diag_out!"
    end select
    
    ! Need pressure for interpolation
    p(:,:,1,:) = diag_SfcP(ncid, nx, ny, nz, nt)
    p(:,:,2:nz+1,:) = diag_Pressure(ncid, nx, ny, nz, nt)

    do t = 1, nt
       do j = 1, ny
          do i = 1, nx
             call vint(outNoInterp(i,j,:,t), p(i,j,:,t), pBot, -dp,  &
                  var%interpType, outInterp(i,j,:,t))
          end do
       end do
    end do

    ! Interpolate underground for heights
    if (var%diagID == GeopHghtPres) then
       allocate (tk(nx,ny,nz,nt), stat=status)
       call mem_error(status, 1, "pres_diag_out")
       tk = diag_TmpK(ncid, nx, ny, nz, nt)

       call hydro_interp(outNoInterp(:,:,1,:), p(:,:,1,:), tk(:,:,1,:), pBot, &
            -dp, outInterp)
       
       deallocate (tk, stat=status)
       call mem_error(status, 2, "pres_diag_out")
    end if
    
    ! Write the diagnostic out to GEMPAK
    do t = 1, nt
       gdat1 = times(t)
       do k = 1, np
          lev = pBot - (k-1) * dp
          call write_gempak(gemid, outInterp(:,:,k,t), nx, ny, gdat1, lev,   &
               -1, Pres, var%gemName)
       end do
    end do
  end subroutine pres_diag_out

  ! ===========================================================================

  ! The next group of functions produce the various diagnostics.  They all take
  ! as input the WRF history file handle (ncid), the WRF domain size 
  ! (nx, ny, and sometimes nz), and the number of times contained in the
  ! history file (nt).

  ! diag_SfcP diagnoses the surface pressure (in hPa).
  function diag_SfcP(ncid, nx, ny, nz, nt) result (ps)
    integer, intent(in)       :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nt) :: ps

    character(len=9), parameter :: funName = "diag_SfcP"

    real, dimension(:,:,:,:), allocatable :: p, w, t, ht
    real, dimension(nx,ny,2,nt)           :: ht12
    real, dimension(nx,ny,nt)             :: p1, q1, sfcq, t1
    integer                               :: status
    
    ! Read in pressure, saving only lowest level
    allocate (p(nx,ny,nz,nt), stat=status)
    call mem_error(status, 1, funName//" p")
    p = diag_Pressure(ncid, nx, ny, nz, nt, si=.true.)  ! Do calculations in Pa
    p1 = p(:,:,1,:)
    deallocate (p, stat=status)
    call mem_error(status, 2, funName//" p")

    ! Handle specific humidity
    allocate (w(nx,ny,nz,nt), stat=status)
    call mem_error(status, 1, funName//" w")
    call read_values(ncid, wrfVar(MixR)%wrfName, w)
    q1 = mix_ratio_to_spec_hum(w(:,:,1,:))
    deallocate (w, stat=status)
    call mem_error(status, 2, funName//" w")
    sfcq = mix_ratio_to_spec_hum(diag_SfcMixRF0(ncid, nx, ny, nz, nt))

    ! Handle temperature
    allocate (t(nx,ny,nz,nt), stat=status)
    call mem_error(status, 1, funName//" t")
    t = diag_Theta(ncid, nx, ny, nz, nt)
    t1 = temp_from_theta_p(t(:,:,1,:), p1)
    deallocate (t, stat=status)
    call mem_error(status, 2, funName//" t")

    ! Handle geopotential height
    allocate (ht(nx,ny,nz+1,nt), stat=status)
    call mem_error(status, 1, funName//" ht")
    ht = diag_GeopHghtS(ncid, nx, ny, nz, nt)
    ht12 = ht(:,:,1:2,:)
    deallocate (ht, stat=status)
    call mem_error(status, 2, funName//" ht")

    ps = .01 * surface_pres(ht12, p1, q1, t1, sfcq)  ! Convert to hPa
  end function diag_SfcP

  ! ===========================================================================

  ! diag_SfcP2m diagnoses the pressure (in hPa) 2 meters above the surface.
  ! The density is approximated by using surface pressure and 2-m temperature.
  function diag_SfcP2m(ncid, nx, ny, nz, nt) result (ps2m)
    integer, intent(in)       :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nt) :: ps, t2, q, rhog, ps2m
    
    ps = diag_SfcP(ncid, nx, ny, nz, nt) * 100  ! In Pa
    t2 = diag_SfcTempF0(ncid, nx, ny, nz, nt)
    q = diag_SfcMixRF0(ncid, nx, ny, nz, nt)
    rhog = Grav * density(ps, t2, q)

    ! Use hydrostatic equation
    ps2m = .01 * (ps - 2*rhog)  ! In hPa
  end function diag_SfcP2m

  ! ===========================================================================

  ! diag_SeaLevP diagnoses the sea level pressure (in hPa).
  function diag_SeaLevP(ncid, nx, ny, nz, nt) result (slp)
    integer, intent(in)       :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nt) :: slp
    
    real, dimension(nx,ny,nz+1,nt) :: ht
    real, dimension(nx,ny,nz,nt)   :: th, p, w
    
    ! Read in WRF data
    ht = diag_GeopHghtS(ncid, nx, ny, nz, nt)
    th = diag_Theta(ncid, nx, ny, nz, nt)
    p = diag_Pressure(ncid, nx, ny, nz, nt, si=.true.)  ! Do calculations in Pa
    call read_values(ncid, wrfVar(MixR)%wrfName, w)  

    ! Compute sea level pressure
    slp = .01 * sea_pres(ht, th, p, mix_ratio_to_spec_hum(w))  ! Convert to hPa
  end function diag_SeaLevP

  ! ===========================================================================
  
  ! diag_ConvPre1h diagnoses 1-hr accumulated convective precipitation (mm).
  ! histInd is the index of the WRF history file being read, and m1h is the
  ! netCDF file handle for the previous WRF history file.
  function diag_ConvPre1h(ncid, nx, ny, nt, histInd, m1h) result (convPrecip1h)
    integer, intent(in)       :: ncid, nx, ny, nt, histInd, m1h
    real, dimension(nx,ny,nt) :: convPrecip1h

    real, dimension(nx,ny,nt) :: accConvPrecip, accConvPrecipm1h
    
    call read_values(ncid, wrfVar(ConvPre)%wrfName, accConvPrecip)

    ! If block handles case where history file contains only one output time
    if (nt > 1) then
       convPrecip1h(:,:,1) = 0
       convPrecip1h(:,:,2:nt) =  &
            accConvPrecip(:,:,2:nt) - accConvPrecip(:,:,1:nt-1)
    else if (histInd > 1) then
       call read_values(m1h, wrfVar(ConvPre)%wrfName, accConvPrecipm1h)
       convPrecip1h = accConvPrecip - accConvPrecipm1h
    else
       convPrecip1h = 0
    end if
  end function diag_ConvPre1h

  ! ===========================================================================

  ! diag_TotPreHr diagnoses x-hr accumulated precipitation (mm), where x is
  ! provided by the input variable hr.  histInd is the index of the WRF history
  ! file being read, and prev is the netCDF file handle for the WRF history
  ! file that contains the precipitation information from x-hr previous.
  function diag_TotPreHr(ncid, nx, ny, nt, hr, histInd, prev) result (hrPrecip)
    integer, intent(in)       :: ncid, nx, ny, nt, hr, histInd, prev
    real, dimension(nx,ny,nt) :: hrPrecip

    real, dimension(nx,ny,nt) :: accPrecip
    integer                   :: diff

    ! Calculate the time index difference between "hr"-hr periods
    diff = hr * timesPer12hr / 12

    accPrecip = diag_TotPre(ncid, nx, ny, nt)
    
    ! If-block handles case where history file contains only one output time
    if (nt > 1) then
       hrPrecip(:,:,1:diff) = 0
       hrPrecip(:,:,diff+1:nt) =  &
            accPrecip(:,:,diff+1:nt) - accPrecip(:,:,1:nt-diff)
    else if (histInd > diff) then
       hrPrecip = accPrecip - diag_TotPre(prev, nx, ny, nt)
    else
       hrPrecip = 0
    end if
  end function diag_TotPreHr

  ! ===========================================================================

  ! diag_TotPre diagnoses total accumulated precipitation (mm).
  function diag_TotPre(ncid, nx, ny, nt) result (totPrecip)
    integer, intent(in)       :: ncid, nx, ny, nt
    real, dimension(nx,ny,nt) :: totPrecip

    real, dimension(nx,ny,nt) :: accConvPrecip
    real, dimension(nx,ny,nt) :: accStabPrecip
    
    call read_values(ncid, wrfVar(ConvPre)%wrfName, accConvPrecip)
    call read_values(ncid, wrfVar(StabPre)%wrfName, accStabPrecip)
    totPrecip = accConvPrecip + accStabPrecip
  end function diag_TotPre

  ! ===========================================================================

  ! diag_Pressure diagnoses pressure. By default, the output is in hPa, but
  ! if si is true, the output is in Pa.
  function diag_Pressure(ncid, nx, ny, nz, nt, si) result (p)
    integer, intent(in)           :: ncid, nx, ny, nz, nt
    logical, intent(in), optional :: si
    real, dimension(nx,ny,nz,nt)  :: p

    real, dimension(nx,ny,nz,nt) :: pb
    logical                      :: hPa
    
    ! Set appropriate output unit flag
    if (present(si)) then
       hPa = .not. si
    else
       hPa = .true.
    end if

    ! Read in WRF data
    call read_values(ncid, wrfVar(PBase)%wrfName, pb)
    call read_values(ncid, wrfVar(PPert)%wrfName, p)

    ! Calculate pressure
    p = p + pb
    if (hPa) p = .01 * p
  end function diag_Pressure
  
  ! ===========================================================================

  ! diag_Theta diagnoses potential temperature (K).
  function diag_Theta(ncid, nx, ny, nz, nt) result (theta)
    integer, intent(in)          :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nz,nt) :: theta!, q
    
    call read_values(ncid, wrfVar(ThetaPert)%wrfName, theta)
    theta = theta + 300  ! WRF stores theta as theta-300
  end function diag_Theta

  ! ===========================================================================

  ! diag_UUnstag diagnoses the wind component (m/s) in the x-coordinate
  ! direction at the mass points (grid box centers).
  function diag_UUnstag(ncid, nx, ny, nz, nt) result (u)
    integer, intent(in)          :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nz,nt) :: u

    real, dimension(nx+1,ny,nz,nt) :: us

    ! Read staggered wind
    call read_values(ncid, wrfVar(UStag)%wrfName, us)

    ! Unstagger it
    u(1:nx,:,:,:) = .5 * (us(1:nx,:,:,:) + us(2:nx+1,:,:,:))
  end function diag_UUnstag
  
  ! ===========================================================================

  ! diag_VUnstag diagnoses the wind component (m/s) in the y-coordinate
  ! direction at the mass points (grid box centers).
  function diag_VUnstag(ncid, nx, ny, nz, nt) result (v)
    integer, intent(in)          :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nz,nt) :: v

    real, dimension(nx,ny+1,nz,nt) :: vs

    ! Read staggered wind
    call read_values(ncid, wrfVar(VStag)%wrfName, vs)

    ! Unstagger it
    v(:,1:ny,:,:) = .5 * (vs(:,1:ny,:,:) + vs(:,2:ny+1,:,:))
  end function diag_VUnstag
  
  ! ===========================================================================

  ! diag_CWtr diagnoses the cloud water mixing ratio (kg/kg), including ice.
  function diag_CWtr(ncid, nx, ny, nz, nt) result (cw)
    integer, intent(in)          :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nz,nt) :: cw

    real, dimension(nx,ny,nz,nt) :: cloud, ice
    
    call read_values(ncid, wrfVar(MixRCloud)%wrfName, cloud)
    call read_values(ncid, wrfVar(MixRIce)%wrfName, ice)
    cw = cloud + ice
  end function diag_CWtr

  ! ===========================================================================
  
  ! diag_WUnstag diagnoses the wind component (m/s) in the vertical direction
  ! at the mass points (grid box centers).
  function diag_WUnstag(ncid, nx, ny, nz, nt) result (w)
    integer, intent(in)          :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nz,nt) :: w

    real, dimension(nx,ny,nz+1,nt) :: ws

    ! Read staggered wind
    call read_values(ncid, wrfVar(WStag)%wrfName, ws)

    ! Unstagger it
    w(:,:,1:nz,:) = .5 * (ws(:,:,1:nz,:) + ws(:,:,2:nz+1,:))
  end function diag_WUnstag

  ! ===========================================================================

  ! diag_GeopHghtS diagnoses geopotential height (m) on w (full) levels.
  function diag_GeopHghtS(ncid, nx, ny, nz, nt) result (height)
    integer, intent(in)            :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nz+1,nt) :: height
    
    real, dimension(nx,ny,nz+1,nt) :: gpb, gpp
    
    call read_values(ncid, wrfVar(GeopBase)%wrfName, gpb)
    call read_values(ncid, wrfVar(GeopPert)%wrfName, gpp)
    height = (gpb + gpp) / Grav
  end function diag_GeopHghtS

  ! ===========================================================================

  ! diag_SfcTempF0 diagnoses 2-m temperature (K), including the initial time.
  function diag_SfcTempF0(ncid, nx, ny, nz, nt) result (temp)
    integer, intent(in)       :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nt) :: temp

    real, dimension(nx,ny,nz,nt) :: potTemp, p
    real, dimension(nx,ny)       :: tk1

    ! For initial time, calculate temperature in lowest grid box.
    potTemp = diag_Theta(ncid, nx, ny, nz, nt)
    p = diag_Pressure(ncid, nx, ny, nz, nt, si = .true.)  ! p in Pa
    tk1 = temp_from_theta_p(potTemp(:,:,1,1), p(:,:,1,1))

    ! Read in WRF 2-m temperature, but WRF doesn't include this variable at
    ! initial time, so assign the temperature in the lowest grid box to it.
    call read_values(ncid, wrfVar(SfcTemp)%wrfName, temp)
    temp(:,:,1) = tk1
  end function diag_SfcTempF0

  ! ===========================================================================

  ! diag_SfcMixRF0 diagnoses 2-m mixing ratio (kg/kg), including the initial
  ! time.
  function diag_SfcMixRF0(ncid, nx, ny, nz, nt) result (w2m)
    integer, intent(in)       :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nt) :: w2m

    real, dimension(nx,ny,nz,nt) :: w

    ! Read in WRF 2-m mixing ratio, but WRF doesn't include this variable at
    ! initial time, so assign the mixing ratio in the lowest grid box to it.
    call read_values(ncid, wrfVar(MixR)%wrfName, w)
    call read_values(ncid, wrfVar(SfcMixR)%wrfName, w2m)
    w2m(:,:,1) = w(:,:,1,1)
  end function diag_SfcMixRF0

  ! ===========================================================================

  ! diag_SfcUF0 diagnoses the 2-m wind component (m/s) in the x-coordinate
  ! direction, including the initial time.
  function diag_SfcUF0(ncid, nx, ny, nz, nt) result (u2m)
    integer, intent(in)       :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nt) :: u2m
    
    real, dimension(nx,ny,nz,nt) :: u
    
    ! Read in WRF 2-m u-component wind, but WRF doesn't include this variable
    ! at initial time, so assign the u-component wind in the lowest grid box to
    ! it.
    u = diag_UUnstag(ncid, nx, ny, nz, nt)
    call read_values(ncid, wrfVar(SfcU)%wrfName, u2m)
    u2m(:,:,1) = u(:,:,1,1)
  end function diag_SfcUF0

  ! ===========================================================================

  ! diag_SfcVF0 diagnoses the 2-m wind component (m/s) in the y-coordinate
  ! direction, including the initial time.
  function diag_SfcVF0(ncid, nx, ny, nz, nt) result (v2m)
    integer, intent(in)       :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nt) :: v2m

    real, dimension(nx,ny,nz,nt) :: v
    
    ! Read in WRF 2-m v-component wind, but WRF doesn't include this variable
    ! at initial time, so assign the v-component wind in the lowest grid box to
    ! it.
    v = diag_VUnstag(ncid, nx, ny, nz, nt)
    call read_values(ncid, wrfVar(SfcV)%wrfName, v2m)
    v2m(:,:,1) = v(:,:,1,1)
  end function diag_SfcVF0

  ! ===========================================================================

  ! diag_LiftIdx diagnoses the Lifted Index (K).
  function diag_LiftIdx(ncid, nx, ny, nz, nt) result (li)
    integer, intent(in)       :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nt) :: li

    real, dimension(nx,ny,nz,nt) :: th, p, w
    integer :: i, j, t

    ! Read in the necessary data from WRF output.
    th = diag_Theta(ncid, nx, ny, nz, nt)
    p = diag_Pressure(ncid, nx, ny, nz, nt)  ! p in hPa
    call read_values(ncid, wrfVar(MixR)%wrfName, w)
    
    ! Calculate lifted index, one profile at a time.
    do t = 1, nt
       do j = 1, ny
          do i = 1, nx
             li(i,j,t) = lifted_index(th(i,j,:,t), p(i,j,:,t), w(i,j,:,t))
          end do
       end do
    end do
  end function diag_LiftIdx

  ! ===========================================================================
  
  ! diag_MeanRH diagnoses the mean relative humidity (%) in the 850-500 mb
  ! layer.
  function diag_MeanRH(ncid, nx, ny, nz, nt) result (mrh)
    integer, intent(in)       :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nt) :: mrh

    real, dimension(nx,ny,nz,nt) :: th, p, w
    real, dimension(nx,ny,nt)    :: ps
    integer :: i, j, t
    
    ! Read in the necessary data from WRF output.
    th = diag_Theta(ncid, nx, ny, nz, nt)
    p = diag_Pressure(ncid, nx, ny, nz, nt)  ! p in hPa
    call read_values(ncid, wrfVar(MixR)%wrfName, w)
    ps = diag_SfcP(ncid, nx, ny, nz, nt)  ! ps in hPa
    
    ! Calculate mean RH, one profile at a time.
    do t = 1, nt
       do j = 1, ny
          do i = 1, nx
             mrh(i,j,t) = mean_rh(th(i,j,:,t), p(i,j,:,t), w(i,j,:,t),  &
                  ps(i,j,t), 850., 500.)
          end do
       end do
    end do
  end function diag_MeanRH

  ! ===========================================================================

  ! diag_SfcSpHum diagnoses 2-m specific humidity (kg/kg).
  function diag_SfcSpHum(ncid, nx, ny, nz, nt) result (q2m)
    integer, intent(in)       :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nt) :: q2m

    q2m = mix_ratio_to_spec_hum(diag_SfcMixRF0(ncid, nx, ny, nz, nt))
  end function diag_SfcSpHum

  ! ===========================================================================
  
  ! diag_SpHum diagnoses specific humidity (kg/kg).
  function diag_SpHum(ncid, nx, ny, nz, nt) result (q)
    integer, intent(in)          :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nz,nt) :: q

    real, dimension(nx,ny,nz,nt) :: w

    call read_values(ncid, wrfVar(MixR)%wrfName, w)
    q = mix_ratio_to_spec_hum(w)
  end function diag_SpHum

  ! ===========================================================================

  ! diag_TmpK diagnoses temperature (K).
  function diag_TmpK(ncid, nx, ny, nz, nt) result (tk)
    integer, intent(in)          :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nz,nt) :: tk

    tk = temp_from_theta_p(diag_Theta(ncid, nx, ny, nz, nt),  &
         diag_Pressure(ncid, nx, ny, nz, nt, si=.true.))
  end function diag_TmpK

  ! ===========================================================================

  ! diag_GeopHght diagnoses geopotential height (m) on half (theta) levels.
  function diag_GeopHght(ncid, nx, ny, nz, nt) result (hght)
    integer, intent(in)          :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nz,nt) :: hght

    real, dimension(:,:), allocatable :: etaTemp
    real, dimension(nx,ny,nz+1,nt)    :: hghtStag
    real, dimension(nx,ny,nt)         :: mu
    real, dimension(nt)               :: pt  ! really a constant
    real, dimension(nz+1)             :: znw
    real, dimension(nz)               :: znu
    integer :: status, i, j, t

    hghtStag = diag_GeopHghtS(ncid, nx, ny, nz, nt)
    mu = diag_DryAir(ncid, nx, ny, nt)

    allocate(etaTemp(nz+1,nt), stat=status)
    call mem_error(status, 1, "diag_GeopHght")
    call read_values(ncid, wrfVar(EtaW)%wrfName, etaTemp)
    znw = etaTemp(:,1)
    deallocate(etaTemp, stat=status)
    call mem_error(status, 2, "diag_GeopHght")

    allocate(etaTemp(nz,nt), stat=status)
    call mem_error(status, 1, "diag_GeopHght")
    call read_values(ncid, wrfVar(EtaMass)%wrfName, etaTemp)
    znu = etaTemp(:,1)
    deallocate(etaTemp, stat=status)
    call mem_error(status, 2, "diag_GeopHght")

    call read_values(ncid, "P_TOP", pt)

    ! Unstagger height by column.
    do t = 1, nt
       do j = 1, ny
          do i = 1, nx
             hght(i,j,:,t) =  &
                  unstag_hght(hghtStag(i,j,:,t), znw, znu, mu(i,j,t), pt(1))
          end do
       end do
    end do
  end function diag_GeopHght

  ! ===========================================================================

  ! diag_DryAir diagnoses the dry air mass in the column (Pa).
  function diag_DryAir(ncid, nx, ny, nt) result (mu)
    integer, intent(in)       :: ncid, nx, ny, nt
    real, dimension(nx,ny,nt) :: mu
    
    real, dimension(nx,ny,nt) :: mub
    
    call read_values(ncid, wrfVar(MuBase)%wrfName, mub)
    call read_values(ncid, wrfVar(MuPert)%wrfName, mu)
    mu = mub + mu
  end function diag_dryAir

  ! ===========================================================================

  ! diag_MUCAPE diagnoses the convective available potential energy following
  ! Doswell and Rasmussen (1994).
  function diag_MUCAPE(ncid, nx, ny, nz, nt) result (c)
    integer, intent(in)       :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nt) :: c

    real, dimension(nx,ny,nz,nt) :: th, p, w, z
    integer :: i, j, t

    th = diag_Theta(ncid, nx, ny, nz, nt)
    p = diag_Pressure(ncid, nx, ny, nz, nt)  ! p in hPa
    call read_values(ncid, wrfVar(MixR)%wrfName, w)
    z = diag_GeopHght(ncid, nx, ny, nz, nt)

    do t = 1, nt
       do j = 1, ny
          do i = 1, nx
             c(i,j,t) = cape(th(i,j,:,t), p(i,j,:,t), w(i,j,:,t), z(i,j,:,t))
          end do
       end do
    end do
  end function diag_MUCAPE

  ! ===========================================================================

  ! diag_MULPL diagnoses the lifted parcel level of the most unstable parcel.
  function diag_MULPL(ncid, nx, ny, nz, nt) result (lpl)
    integer, intent(in)       :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nt) :: lpl

    real, dimension(nx,ny,nz,nt) :: th, p, w, z
    real, dimension(nx,ny,nt)    :: ter  ! Really constant in time
    integer :: i, j, t

    th = diag_Theta(ncid, nx, ny, nz, nt)
    p = diag_Pressure(ncid, nx, ny, nz, nt)  ! p in hPa
    call read_values(ncid, wrfVar(MixR)%wrfName, w)
    z = diag_GeopHght(ncid, nx, ny, nz, nt)
    call read_values(ncid, wrfVar(TerHght)%wrfName, ter)

    do t = 1, nt
       do j = 1, ny
          do i = 1, nx
             lpl(i,j,t) = cape(th(i,j,:,t), p(i,j,:,t), w(i,j,:,t),  &
                  z(i,j,:,t), ter(i,j,1), lplOutI = .true.)
          end do
       end do
    end do
  end function diag_MULPL

  ! ===========================================================================
  
  ! diag_SBCAPE diagnoses the convective available potential energy of a 
  ! surface parcel following Doswell and Rasmussen (1994).
  function diag_SBCAPE(ncid, nx, ny, nz, nt) result (c)
    integer, intent(in)       :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nt) :: c
    
    real, dimension(nx,ny,nz,nt) :: th, p, w, z
    integer :: i, j, t

    th = diag_Theta(ncid, nx, ny, nz, nt)
    p = diag_Pressure(ncid, nx, ny, nz, nt)  ! p in hPa
    call read_values(ncid, wrfVar(MixR)%wrfName, w)
    z = diag_GeopHght(ncid, nx, ny, nz, nt)

    do t = 1, nt
       do j = 1, ny
          do i = 1, nx
             c(i,j,t) = sfcape(th(i,j,:,t), p(i,j,:,t), w(i,j,:,t), z(i,j,:,t))
          end do
       end do
    end do
  end function diag_SBCAPE
  
  ! ===========================================================================

  ! diag_SBCIN diagnoses the convective inhibition of a surface parcel.
  function diag_SBCIN(ncid, nx, ny, nz, nt) result (c)
    integer, intent(in)       :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nt) :: c
    
    real, dimension(nx,ny,nz,nt) :: th, p, w, z
    integer :: i, j, t

    th = diag_Theta(ncid, nx, ny, nz, nt)
    p = diag_Pressure(ncid, nx, ny, nz, nt)  ! p in hPa
    call read_values(ncid, wrfVar(MixR)%wrfName, w)
    z = diag_GeopHght(ncid, nx, ny, nz, nt)
    
    do t = 1, nt
       do j = 1, ny
          do i = 1, nx
             c(i,j,t) = sfcape(th(i,j,:,t), p(i,j,:,t), w(i,j,:,t),  &
                  z(i,j,:,t), .true.)
          end do
       end do
    end do
  end function diag_SBCIN

  ! ===========================================================================
  
  ! diag_MixCAPE diagnoses the convective available potential energy of a 
  ! near-surface mixed layer parcels.
  function diag_MixCAPE(ncid, nx, ny, nz, nt) result (c)
    integer, intent(in)       :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nt) :: c
    
    real, dimension(nx,ny,nz,nt) :: th, p, w, z
    integer :: i, j, t

    th = diag_Theta(ncid, nx, ny, nz, nt)
    p = diag_Pressure(ncid, nx, ny, nz, nt)  ! p in hPa
    call read_values(ncid, wrfVar(MixR)%wrfName, w)
    z = diag_GeopHght(ncid, nx, ny, nz, nt)

    do t = 1, nt
       do j = 1, ny
          do i = 1, nx
             c(i,j,t) = mlcape(th(i,j,:,t), p(i,j,:,t), w(i,j,:,t), z(i,j,:,t))
          end do
       end do
    end do
  end function diag_MixCAPE
  
  ! ===========================================================================

  ! diag_MixCIN diagnoses the convective inhibition of near-surface mixed
  ! layer parcels.
  function diag_MixCIN(ncid, nx, ny, nz, nt) result (c)
    integer, intent(in)       :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nt) :: c
    
    real, dimension(nx,ny,nz,nt) :: th, p, w, z
    integer :: i, j, t

    th = diag_Theta(ncid, nx, ny, nz, nt)
    p = diag_Pressure(ncid, nx, ny, nz, nt)  ! p in hPa
    call read_values(ncid, wrfVar(MixR)%wrfName, w)
    z = diag_GeopHght(ncid, nx, ny, nz, nt)
    
    do t = 1, nt
       do j = 1, ny
          do i = 1, nx
             c(i,j,t) = mlcape(th(i,j,:,t), p(i,j,:,t), w(i,j,:,t),  &
                  z(i,j,:,t), .true.)
          end do
       end do
    end do
  end function diag_MixCIN

  ! ===========================================================================

  ! diag_MixLCL diagnoses the lifing condensation level of near-surface mixed
  ! layer parcels.
  function diag_MixLCL(ncid, nx, ny, nz, nt) result (lcl)
    integer, intent(in)       :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nt) :: lcl
    
    real, dimension(nx,ny,nz,nt) :: th, p, w, z
    real, dimension(nx,ny,nt)    :: ter  ! Really constant in time
    integer :: i, j, t
    
    th = diag_Theta(ncid, nx, ny, nz, nt)
    p = diag_Pressure(ncid, nx, ny, nz, nt)  ! p in hPa
    call read_values(ncid, wrfVar(MixR)%wrfName, w)
    z = diag_GeopHght(ncid, nx, ny, nz, nt)
    call read_values(ncid, wrfVar(TerHght)%wrfName, ter)

    do t = 1, nt
       do j = 1, ny
          do i = 1, nx
             lcl(i,j,t) = lcl_lfc(th(i,j,:,t), p(i,j,:,t), w(i,j,:,t),  &
                  z(i,j,:,t), ter(i,j,1), .true.)
          end do
       end do
    end do
  end function diag_MixLCL

  ! ===========================================================================

  ! diag_MixLFC diagnoses the level of free convection of near-surface mixed
  ! layer parcels.
  function diag_MixLFC(ncid, nx, ny, nz, nt) result (lfc)
    integer, intent(in)       :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nt) :: lfc
    
    real, dimension(nx,ny,nz,nt) :: th, p, w, z
    real, dimension(nx,ny,nt)    :: ter  ! Really constant in time
    integer :: i, j, t

    th = diag_Theta(ncid, nx, ny, nz, nt)
    p = diag_Pressure(ncid, nx, ny, nz, nt)  ! p in hPa
    call read_values(ncid, wrfVar(MixR)%wrfName, w)
    z = diag_GeopHght(ncid, nx, ny, nz, nt)
    call read_values(ncid, wrfVar(TerHght)%wrfName, ter)
  
    do t = 1, nt
       do j = 1, ny
          do i = 1, nx
             lfc(i,j,t) = lcl_lfc(th(i,j,:,t), p(i,j,:,t), w(i,j,:,t),  &
                  z(i,j,:,t), ter(i,j,1), .false.)
          end do
       end do
    end do
  end function diag_MixLFC

  ! ===========================================================================

  ! diag_UStorm diagnoses the u component (m) of the Bunkers storm motion
  ! vector.
  function diag_UStorm(ncid, nx, ny, nz, nt) result (us)
    integer, intent(in)       :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nt) :: us

    real, dimension(nx,ny,nz,nt) :: u, v, z
    real, dimension(nx,ny,nt) :: u10, v10, ter

    u = diag_UUnstag(ncid, nx, ny, nz, nt)
    v = diag_VUnstag(ncid, nx, ny, nz, nt)
    z = diag_GeopHght(ncid, nx, ny, nz, nt)
    u10 = diag_SfcUF0(ncid, nx, ny, nz, nt)
    v10 = diag_SfcVF0(ncid, nx, ny, nz, nt)
    call read_values(ncid, wrfVar(TerHght)%wrfName, ter)
    
    us = storm_motion(u, v, z, u10, v10, ter(:,:,1), .true.)
  end function diag_UStorm

  ! ===========================================================================
  
  ! diag_VStorm diagnoses the v component (m) of the Bunkers storm motion
  ! vector.
  function diag_VStorm(ncid, nx, ny, nz, nt) result (vs)
    integer, intent(in)       :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nt) :: vs

    real, dimension(nx,ny,nz,nt) :: u, v, z
    real, dimension(nx,ny,nt)    :: u10, v10, ter

    u = diag_UUnstag(ncid, nx, ny, nz, nt)
    v = diag_VUnstag(ncid, nx, ny, nz, nt)
    z = diag_GeopHght(ncid, nx, ny, nz, nt)
    u10 = diag_SfcUF0(ncid, nx, ny, nz, nt)
    v10 = diag_SfcVF0(ncid, nx, ny, nz, nt)
    call read_values(ncid, wrfVar(TerHght)%wrfName, ter)
    
    vs = storm_motion(u, v, z, u10, v10, ter(:,:,1), .false.)
  end function diag_VStorm

  ! ===========================================================================

  ! diag_StormRHel diagnoses the 0-1-km storm relative helicity (m2/s2).
  function diag_StormRHel(ncid, nx, ny, nz, nt) result (srh)
    integer, intent(in)       :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nt) :: srh

    real, dimension(5), parameter :: Heights = (/ 0, 250, 500, 750, 1000 /)

    real, dimension(nx,ny,nz,nt) :: u, v, z
    real, dimension(nx,ny,nt)    :: u10, v10, ter

    u = diag_UUnstag(ncid, nx, ny, nz, nt)
    v = diag_VUnstag(ncid, nx, ny, nz, nt)
    z = diag_GeopHght(ncid, nx, ny, nz, nt)
    u10 = diag_SfcUF0(ncid, nx, ny, nz, nt)
    v10 = diag_SfcVF0(ncid, nx, ny, nz, nt)
    call read_values(ncid, wrfVar(TerHght)%wrfName, ter)

    srh = storm_rel_hel(u, v, z, u10, v10, ter(:,:,1), Heights)
  end function diag_StormRHel

  ! ===========================================================================

  ! diag_PreWater diagnoses the column integrated precipitable water (cm).
  ! This function is largely taken from that in RIP4.
  function diag_PreWater(ncid, nx, ny, nz, nt) result(pw)
    integer, intent(in)          :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nt)    :: pw
    
    ! This parameter implicitly includes a conversion from m to cm
    real, parameter              :: GravRhoWInv = 1. / (Grav * 10.)

    real, dimension(nx,ny,nz,nt) :: p, w

    p = diag_Pressure(ncid, nx, ny, nz, nt, si=.true.) ! units of Pa
    call read_values(ncid, wrfVar(MixR)%wrfName, w)

    pw = GravRhoWInv * sum( .5 * (w(:,:,1:nz-1,:) + w(:,:,2:nz,:))  &
         * (p(:,:,1:nz-1,:) - p(:,:,2:nz,:)), dim=3, mask=  &
         p(:,:,2:nz,:) > 30000. )
  end function diag_PreWater
  
  ! ===========================================================================
  
  ! diag_Refl diagnoses reflectivity (dBZ).
  function diag_Refl(ncid, nx, ny, nz, nt, graupel) result(dbz)
    integer, intent(in)           :: ncid, nx, ny, nz, nt
    logical, intent(in), optional :: graupel
    real, dimension(nx,ny,nz,nt)  :: dbz
    
    real, dimension(:,:,:,:), allocatable :: p, qvp
    real, dimension(nx,ny,nz,nt)          :: tmk, qra, qsn, qgr, rho
    integer                               :: status
    logical                               :: graup

    if (present(graupel)) then
       graup = graupel
    else
       graup = .true.
    end if

    allocate(p(nx,ny,nz,nt), qvp(nx,ny,nz,nt), stat=status)
    call mem_error(status, 1, "diag_Refl")
    
    p = diag_Pressure(ncid, nx, ny, nz, nt, si=.true.) ! units of Pa
    tmk = diag_TmpK(ncid, nx, ny, nz, nt)
    call read_values(ncid, wrfVar(MixR)%wrfName, qvp)
    
    ! Make sure hydrometeor mixing ratios are nonnegative
    call read_values(ncid, wrfVar(MixRRain)%wrfName, qra)
    where (qra < 0.) qra = 0.    
    call read_values(ncid, wrfVar(MixRSnow)%wrfName, qsn)
    where (qsn < 0.) qsn = 0.
    if (graup) then
       call read_values(ncid, wrfVar(MixRGraup)%wrfName, qgr)
       where (qgr < 0.) qgr = 0.
    else
       qgr = 0.
    end if

    rho = density(p, tmk, qvp)

    deallocate(p, qvp, stat=status)
    call mem_error(status, 2, "diag_Refl")

    dbz = calcdbz(rho, tmk, qra, qsn, qgr)
  end function diag_Refl
  
  ! ===========================================================================
  
  ! diag_CompRefl diagnoses composite reflectivity (dBZ), but requires that the
  ! WRF output contain the graupel mixing ratio.
  function diag_CompRefl(ncid, nx, ny, nz, nt) result(cref)
    integer, intent(in)       :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nt) :: cref

    cref = maxval(diag_Refl(ncid, nx, ny, nz, nt), dim=3)
  end function diag_CompRefl

  ! ===========================================================================

  ! diag_WindGust diagnoses the surface wind gusts (m/s).
  function diag_WindGust(ncid, nx, ny, nz, nt) result(gust)
    integer, intent(in)       :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nt) :: gust

    real, dimension(nx,ny,nz,nt) :: u, v, zagl
    real, dimension(nx,ny,nt)    :: u10, v10, zpbl

    u = diag_UUnstag(ncid, nx, ny, nz, nt)
    v = diag_VUnstag(ncid, nx, ny, nz, nt)
    zagl = diag_ZAbvGrnd(ncid, nx, ny, nz, nt)
    u10 = diag_SfcUF0(ncid, nx, ny, nz, nt)
    v10 = diag_SfcVF0(ncid, nx, ny, nz, nt)
    call read_values(ncid, wrfVar(PBLHght)%wrfName, zpbl)

    gust = calcgust(u, v, zagl, u10, v10, zpbl)
  end function diag_WindGust

  ! ===========================================================================

  ! diag_ZAbvGrnd diagnoses the height of the model level above the ground (m).
  function diag_ZAbvGrnd(ncid, nx, ny, nz, nt) result(zag)
    integer, intent(in)          :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nz,nt) :: zag

    real, dimension(nx,ny,nz,nt) :: z
    real, dimension(nx,ny,nt)    :: ter

    z = diag_GeopHght(ncid, nx, ny, nz, nt)
    call read_values(ncid, wrfVar(TerHght)%wrfName, ter)

    zag = z - spread(ter, 3, nz)
  end function diag_ZAbvGrnd
end module registry
