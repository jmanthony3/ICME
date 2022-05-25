c
c=====================================================================72
c
c------- crystal Type
c
      integer   kFCC, kBCC, kHCP
      parameter (kFCC=1, kBCC=2, kHCP=3)
c
c------- number of Crystal Slip Systems
c
      integer	MAX_SLIP
      parameter (MAX_SLIP = 48)

      integer   kSlipFCC, kSlipBCC, kSlipHCP
      parameter (kSlipFCC=12, kSlipBCC=12, kSlipHCP=12)
c
c------- maximum number of vertex stresses of SCYS
c
      integer	MAX_VTX
      parameter (MAX_VTX  = 241)
c
c------- code for distribution of crystal lattice orientations
c
      integer   kODF_same_all,   kODF_diff_elems, 
     &          kODF_diff_inpts, kODF_from_file, 
     &          kODF_from_abq
      parameter (kODF_same_all=1,   kODF_diff_elems=2, 
     &           kODF_diff_inpts=3, kODF_from_file=4,
     &           kODF_from_abq=5)
c
c------- size/subscripts of state variable array
c
      integer   NSTAV
      integer   kEVOL, kPRES, kDETVe, kWPNORM, kSSACT
      parameter (NSTAV=5)
      parameter (kEVOL=1, kPRES=2, kDETVe=3, kWPNORM=4, kSSACT=5)
c
c------- size/subscripts for array of useful scalar variables
c
      integer   NEQVA
      integer   kEQP_n, kEQP, kMISES_n, kMISES, kSHRATE_n, kSHRATE,
     &          kGAMTOT_n, kGAMTOT
      parameter (NEQVA=8)
c      parameter (kEQP_n=1, kEQP=4, kMISES_n=2, kMISES=5, 
c     &           kSHRATE_n=3, kSHRATE=6, kGAMTOT_n=7, kGAMTOT=8)
      parameter (kEQP_n=1, kMISES_n=2, kSHRATE_n=3, kGAMTOT_n=4, 
     &           kEQP  =5, kMISES  =6, kSHRATE  =7, kGAMTOT  =8)
c
c------- parameters for crystal elasticity
c
      integer   kELAS_ISO, kELAS_ANI
      parameter (kELAS_ISO=1, kELAS_ANI=2)
c
c------- parameters for crystal kinetic equation
c
      integer   kGAMDOT, kdGAMdTAU, kdGAMdKAPP
      parameter (kGAMDOT=1, kdGAMdTAU=2, kdGAMdKAPP=3)
c
c------- parameters for integration of hardness
c
      integer   kHARD_EXPL, kHARD_MIDP, kHARD_ANAL
      parameter (kHARD_EXPL=1, kHARD_MIDP=2, kHARD_ANAL=3)
c
c------- other dimensions
c
      integer    NVECS, NVECA, NKAPP, NPROPS, DIMS, DIMS2
      parameter  (NVECS=5, NVECA=3, NKAPP=MAX_SLIP,
     &            NPROPS=15, DIMS=3, DIMS2=6)
c
c------- number stages and points in load history (mps)
c
      integer   NSTAGE, NLOAD
      parameter (NSTAGE=100, NLOAD=10000)
c
c------- number of stress/strain components (for ABAQUS)
c
      integer    NTENS, NDI, NSHR
      parameter  (NTENS=6, NDI=3, NSHR=3)
c      parameter (NTENS=4, NDI=3, NSHR=1)
c
c------- I/O units (mps)
c
c---------- general
      integer   INPUT, OUTPT
      parameter (INPUT=10, OUTPT=12)

      integer XTAL_I, XTAL_O, XTAL_E, XTAL_TXT_IN, XTAL_TXT_OUT
      parameter (XTAL_I=80)
      parameter (XTAL_O=81)
      parameter (XTAL_E=82)
      parameter (XTAL_TXT_IN=83)
      parameter (XTAL_TXT_OUT=84)

c      integer SXFILE_I           ! define in SUB CrystalMaterialData
c      parameter (SXFILE_I=91)
c
c---------- for a crystal
      integer   XTAL_STRESS_O, XTAL_STRAIN_O, XTAL_EFFSS_O, 
     &          XTAL_TRUESS_O, XTAL_ITER_O
      parameter (XTAL_STRESS_O=85, XTAL_STRAIN_O=86, XTAL_EFFSS_O=87, 
     &           XTAL_TRUESS_O=88, XTAL_ITER_O=89)
c
c---------- for an aggregate
      integer   AGG_EFFSS_O
      parameter (AGG_EFFSS_O=90)
c
c------- parameters to check singularities
c
      real*8    TINY, BIG
      parameter (TINY=1.0d-16, BIG=1.0d30)
c
c------- Crystal material model number
c
      integer XTAL_MODEL
      parameter(XTAL_MODEL = 110)
c
c------- general array dimensions (some from maxdim.default)
c
      integer  NUMEL_T, NUMVP_T, NUMQPT_T
      parameter (NUMEL_T  = 1)
      parameter (NUMVP_T  = 1)
      parameter (NUMQPT_T = 1)

      integer MAX_GRN, MAX_ORIEN
      parameter(MAX_GRN=1000, MAX_ORIEN=10000)
c
c------- parameters for printing from material routine
c
      integer kPRINT_MODEL, kPRINT_ELEM, kPRINT_QPT, kPRINT_GRN
      parameter (kPRINT_MODEL=1, kPRINT_ELEM=1, kPRINT_QPT=1,
     &           kPRINT_GRN=1)
c
c------- parameter for secant line search in global iterations
c
      integer MAX_SEARCH_ITERS
      parameter (MAX_SEARCH_ITERS=50)
c
c------- error codes for numerical integration of XTAL model
c
      integer XTAL_CONVERGED, XTAL_SING_JACOBIAN, XTAL_LS_FAILED, 
     &        XTAL_MAX_ITERS_HIT
      parameter (XTAL_CONVERGED=0, XTAL_SING_JACOBIAN=1, 
     &           XTAL_LS_FAILED=2, XTAL_MAX_ITERS_HIT=3)
c
c------- error codes for global convergence of velocity iterations
c
       integer kGLOBAL_CONVERGED, kGLOBAL_CONTINUE, kGLOBAL_FAILED
       parameter (kGLOBAL_CONVERGED=0, kGLOBAL_CONTINUE=1, 
     &            kGLOBAL_FAILED=2)
c
c------- codes for increasing/decreasing time step in global iters
c
       integer kDECREASE_STEP, kSAME_STEP, kINCREASE_STEP
       parameter (kDECREASE_STEP=-1, kSAME_STEP=0, kINCREASE_STEP=1)
c
c------- counter for crystals no converging
c
       integer NXTALS_NC_LIMIT
       parameter (NXTALS_NC_LIMIT=50)
c
c=====================================================================72
c
