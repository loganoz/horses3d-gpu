!
!////////////////////////////////////////////////////////////////////////
!
!      ProblemFile.f90
!      Created: June 26, 2015 at 8:47 AM 
!      By: David Kopriva  
!
!      The Problem File contains user defined procedures
!      that are used to "personalize" i.e. define a specific
!      problem to be solved. These procedures include initial conditions,
!      exact solutions (e.g. for tests), etc. and allow modifications 
!      without having to modify the main code.
!
!      The procedures, *even if empty* that must be defined are
!
!      UserDefinedSetUp
!      UserDefinedInitialCondition(mesh)
!      UserDefinedPeriodicOperation(mesh)
!      UserDefinedFinalize(mesh)
!      UserDefinedTermination
!
!//////////////////////////////////////////////////////////////////////// 
!
!  ********************************************************************
!  This module contains the definitions for the channel flow simulation
!  ********************************************************************
   module OscarChannel
#ifdef _HAS_MPI_
      use mpi
#endif
      use SMConstants
#ifdef NAVIERSTOKES
      use FluidData_NS, only: RefValues_t
#endif
      implicit none
!
!     Parameters
!     ----------
      real(kind=RP), parameter :: alpha    =  0.9_RP
!
!     Module variables 
!     ----------------
      real(kind=RP) :: flowDir(NDIM)            ! Flow direction
      ! real(kind=RP) :: f_x(2), f_y(2), f_z(3)   ! Correction factors in x and y direction
      real(kind=RP) :: f_x(2), f_y(2), f_z(2)   ! Correction factors in x and y direction
      real(kind=RP) :: dpdx, dpdy, dpdz         ! Pressure gradient
      
      real(kind=RP), allocatable :: dp_part(:)  ! dpdx, dpdy or dpdz in the different partitions (MPI)
      integer                    :: nProcs      ! Number of MPI partitions
      
      contains
#ifdef NAVIERSTOKES
         subroutine InitializeOscarChannel(refValues_)
            implicit none
            type(RefValues_t),      intent(in)  :: refValues_
            integer :: ierr
!
!           Initialize the flow direction
!           -----------------------------
            flowDir(1) = cos(refValues_ % AoATheta*pi/180._RP) * cos(refValues_ % AoAPhi*pi/180._RP)
            flowDir(2) = sin(refValues_ % AoAtheta*pi/180._RP) * cos(refValues_ % AoAphi*pi/180._RP)
            flowDir(3) = sin(refValues_ % AoAphi*pi/180._RP)
!
!           Initialize dpdx and dpdy
!           ------------------------
            dpdx  = 3.08e-3_RP
            ! dpdy  = 3.08e-3_RP
            ! dpdz  = 0._RP 
            dpdy  = 0._RP 
            ! dpdz  = 3.08e-3_RP
            dpdz  = 3.08e-5_RP
!
!           Initialize MPI
!           --------------
#ifdef _HAS_MPI_
            call mpi_comm_size(MPI_COMM_WORLD, nProcs, ierr)
            allocate( dp_part(nProcs) )
#else
            nProcs = 1
#endif
         end subroutine InitializeOscarChannel
#endif
   end module OscarChannel
!
!//////////////////////////////////////////////////////////////////////// 
!
         SUBROUTINE UserDefinedStartup
!
!        --------------------------------
!        Called before any other routines
!        --------------------------------
!
            IMPLICIT NONE  
         END SUBROUTINE UserDefinedStartup
!
!//////////////////////////////////////////////////////////////////////// 
! 
         SUBROUTINE UserDefinedFinalSetup(mesh &
#if defined(NAVIERSTOKES)
                                        , thermodynamics_ &
                                        , dimensionless_  &
                                        , refValues_ & 
#endif
#if defined(CAHNHILLIARD)
                                        , multiphase_ &
#endif
                                        )
!
!           ----------------------------------------------------------------------
!           Called after the mesh is read in to allow mesh related initializations
!           or memory allocations.
!           ----------------------------------------------------------------------
!
            USE HexMeshClass
            use PhysicsStorage
            use FluidData
            use OscarChannel
            IMPLICIT NONE
            CLASS(HexMesh)                      :: mesh
#if defined(NAVIERSTOKES)
            type(Thermodynamics_t), intent(in)  :: thermodynamics_
            type(Dimensionless_t),  intent(in)  :: dimensionless_
            type(RefValues_t),      intent(in)  :: refValues_
#endif
#if defined(CAHNHILLIARD)
            type(Multiphase_t),     intent(in)  :: multiphase_
#endif
#if defined(NAVIERSTOKES)
!
!           Initialize Channel variables
!           ----------------------------
            call InitializeOscarChannel(refValues_)
#endif
         END SUBROUTINE UserDefinedFinalSetup
!
!//////////////////////////////////////////////////////////////////////// 
! 
         subroutine UserDefinedInitialCondition(mesh &
#if defined(NAVIERSTOKES)
                                        , thermodynamics_ &
                                        , dimensionless_  &
                                        , refValues_ & 
#endif
#if defined(CAHNHILLIARD)
                                        , multiphase_ &
#endif
                                        )
!
!           ------------------------------------------------
!           called to set the initial condition for the flow
!              - by default it sets an uniform initial
!                 condition.
!           ------------------------------------------------
!
            use smconstants
            use physicsstorage
            use hexmeshclass
            use fluiddata
#if defined(NAVIERSTOKES)
            use OscarChannel, only: flowDir
#endif
            implicit none
            class(hexmesh)                      :: mesh
#if defined(NAVIERSTOKES)
            type(Thermodynamics_t), intent(in)  :: thermodynamics_
            type(Dimensionless_t),  intent(in)  :: dimensionless_
            type(RefValues_t),      intent(in)  :: refValues_
#endif
#if defined(CAHNHILLIARD)
            type(Multiphase_t),     intent(in)  :: multiphase_
#endif
!
!           ---------------
!           Local variables
!           ---------------
!
#if defined(NAVIERSTOKES)
            real(kind=RP) :: x(3), ydivL, L
            INTEGER       :: i, j, k, eID
            real(kind=RP) :: rho , u , v , w , p
            real(kind=RP) :: rho_0, p_0
            integer       :: Nx, Ny, Nz
            
            rho_0 = 1.0_RP 
            p_0   = 1.0_RP
            L     = 2.0_RP
            
            associate( gamma   => thermodynamics_ % gamma , &
                       gammaM2 => dimensionless_ % gammaM2  ) 
            do eID = 1, size(mesh % elements)
               Nx = mesh % elements(eID) % Nxyz(1)
               Ny = mesh % elements(eID) % Nxyz(2)
               Nz = mesh % elements(eID) % Nxyz(3)

               do k = 0, Nz   ;  do j = 0, Ny   ;  do i = 0, Nx 
                  
                  x = mesh % elements(eID) % geom % x(:,i,j,k)
                  ydivL = x(2) / L


                  rho = rho_0
                  u   =  ( 6.0_RP * (ydivL - ydivL **2) ) * flowDir(1) 
                  ! v   =  ( 3.0_RP / 2.0_RP * (1.0_RP-(x(3))**2) ) * flowDir(2)
                  ! w   =  0.0_RP
                  v   =  0.0_RP
                  w   =  ( 6.0_RP * (ydivL - ydivL **2) ) * flowDir(3)
                  ! w   =  ( 3.0_RP / 2.0_RP * (1.0_RP-(x(2))**2) ) * flowDir(3)
                  p   =  p_0 / ( gammaM2 ) 

                  !Add random noise
                  !rho = rho + RandomNormal(0._RP,0.5_RP)
                  !u   = u   + RandomNormal(0._RP,0.5_RP)
                  !v   = v   + RandomNormal(0._RP,0.5_RP)
                  !w   = w   + RandomNormal(0._RP,0.5_RP)
                  
                  mesh % elements(eID) % storage % Q(1,i,j,k) = rho
                  mesh % elements(eID) % storage % Q(2,i,j,k) = rho*u
                  mesh % elements(eID) % storage % Q(3,i,j,k) = rho*v
                  mesh % elements(eID) % storage % Q(4,i,j,k) = rho*w
                  mesh % elements(eID) % storage % Q(5,i,j,k) = p / (gamma - 1.0_RP) + 0.5_RP * rho * (u*u + v*v + w*w)
                  
               end do         ;  end do         ;  end do
               
            end do
            end associate
#endif            
            contains
               function RandomNormal(mean,sigma)
                  ! RandomNormal function taken from a StackOverflow post

                  ! mean  : mean of distribution
                  ! sigma : number of standard deviations

                  use SMConstants

                  implicit none

                  real(kind=RP) :: RandomNormal
                  real(kind=RP) :: mean
                  real(kind=RP) :: sigma 

                  real          :: rand_num
                  real(kind=RP) :: tmp
                  real(kind=RP) :: fac
                  real(kind=RP) :: gsave
                  real(kind=RP) :: rsq
                  real(kind=RP) :: r1
                  real(kind=RP) :: r2
                  INTEGER       :: flag

                  SAVE :: flag
                  SAVE :: gsave
                  DATA flag /0/

                  
                  IF (flag .EQ. 0) THEN
                        rsq = 2.0_RP

                        DO WHILE(rsq .GE. 1.0_RP .OR. rsq .EQ. 0.0_RP) 
                              CALL RANDOM_NUMBER(rand_num)
                              r1  = 2.0_RP * rand_num - 1.0_RP
                              CALL RANDOM_NUMBER(rand_num)
                              r2  = 2.0_RP * rand_num - 1.0_RP
                              rsq = r1 * r1 + r2 * r2
                        ENDDO

                        fac   = SQRT( -2.0_RP * LOG( rsq ) / rsq )
                        gsave = r1 * fac
                        tmp   = r2 * fac
                        flag  = 1
                  ELSE
                        tmp   = gsave
                        flag  = 0
                  ENDIF

                  RandomNormal = tmp * sigma + mean
               end function RandomNormal    
         end subroutine UserDefinedInitialCondition
#if defined(NAVIERSTOKES)
         subroutine UserDefinedState1(x, t, nHat, Q, thermodynamics_, dimensionless_, refValues_)
!
!           -------------------------------------------------
!           Used to define an user defined boundary condition
!           -------------------------------------------------
!
            use SMConstants
            use PhysicsStorage
            use FluidData
            implicit none
            real(kind=RP), intent(in)     :: x(NDIM)
            real(kind=RP), intent(in)     :: t
            real(kind=RP), intent(in)     :: nHat(NDIM)
            real(kind=RP), intent(inout)  :: Q(NCONS)
            type(Thermodynamics_t),    intent(in)  :: thermodynamics_
            type(Dimensionless_t),     intent(in)  :: dimensionless_
            type(RefValues_t),         intent(in)  :: refValues_
         end subroutine UserDefinedState1

         subroutine UserDefinedNeumann(x, t, nHat, U_x, U_y, U_z)
!
!           --------------------------------------------------------
!           Used to define a Neumann user defined boundary condition
!           --------------------------------------------------------
!
            use SMConstants
            use PhysicsStorage
            use FluidData
            implicit none
            real(kind=RP), intent(in)     :: x(NDIM)
            real(kind=RP), intent(in)     :: t
            real(kind=RP), intent(in)     :: nHat(NDIM)
            real(kind=RP), intent(inout)  :: U_x(NGRAD)
            real(kind=RP), intent(inout)  :: U_y(NGRAD)
            real(kind=RP), intent(inout)  :: U_z(NGRAD)
         end subroutine UserDefinedNeumann
#endif
!
!//////////////////////////////////////////////////////////////////////// 
! 
         SUBROUTINE UserDefinedPeriodicOperation(mesh, time, dt, Monitors)
!
!           -------------------------------------------------------------------------
!           Called before every time-step to allow periodic operations to be performed
!           Here:
!              * Reading mean momentum/velocity in the volume
!              * Obtaining fTurbulentChannel according to it
!           -------------------------------------------------------------------------
!
            use SMConstants
            USE HexMeshClass
#if defined(NAVIERSTOKES)
            use MonitorsClass
            use OscarChannel
#endif
            IMPLICIT NONE
            !-arguments---------------------------------------------------
            CLASS(HexMesh)               :: mesh
            real(kind=RP)                :: time
            real(kind=RP)                :: dt
#if defined(NAVIERSTOKES)
            type(Monitor_t), intent(in) :: monitors
#else
            logical, intent(in) :: monitors
#endif
#if defined(NAVIERSTOKES)
            !-local-variables---------------------------------------------
            real(kind=RP) :: meanMomentum(NDIM)   ! Mean momentum
#ifdef _HAS_MPI_
            integer       :: ierr, dp_ind
            !-------------------------------------------------------------
!
!           Get the same global dpdx,dpdy,dpdz in all partitions
!           ----------------------------------------------------
            if (nProcs > 1) then   
!
!              dpdx
!              ----
               call MPI_Allgather( dpdx, 1, MPI_DOUBLE, dp_part, 1, MPI_DOUBLE, MPI_COMM_WORLD, ierr)
               dp_ind  = maxloc(abs(dp_part),1)
               dpdx    = dp_part(dp_ind)
!
!              dpdy
!              ----
               call MPI_Allgather( dpdy, 1, MPI_DOUBLE, dp_part, 1, MPI_DOUBLE, MPI_COMM_WORLD, ierr)
               dp_ind  = maxloc(abs(dp_part),1)
               dpdy    = dp_part(dp_ind)
!
!              dpdz
!              ----
               call MPI_Allgather( dpdz, 1, MPI_DOUBLE, dp_part, 1, MPI_DOUBLE, MPI_COMM_WORLD, ierr)
               dp_ind  = maxloc(abs(dp_part),1)
               dpdz    = dp_part(dp_ind)
            end if
#endif      
!
!           Compute the mean velocity/momentum
!           ----------------------------------
            meanMomentum = monitors % volumeMonitors(1) % getLast()
            
!
!           Compute the correction parameter
!           --------------------------------

            f_x(1) = ( flowDir(1) - meanMomentum(1) ) * alpha / dt
            f_x(2) = (1._RP - alpha) * dpdx
            
            f_y(1) = ( flowDir(2) - meanMomentum(2) ) * alpha / dt
            f_y(2) = (1._RP - alpha) * dpdy

            f_z(1) = ( flowDir(3) - meanMomentum(3) ) * alpha / dt
            f_z(2) = (1._RP - alpha) * dpdz
#endif
         END SUBROUTINE UserDefinedPeriodicOperation
!
!//////////////////////////////////////////////////////////////////////// 
! 
#if defined(NAVIERSTOKES)
         subroutine UserDefinedSourceTermNS(x, Q, time, S, thermodynamics_, dimensionless_, refValues_) ! , Q, dt
!
!           --------------------------------------------
!           Called to apply source terms to the equation
!           --------------------------------------------
!
            use SMConstants
            USE HexMeshClass
            use PhysicsStorage
            use FluidData
            use OscarChannel
            IMPLICIT NONE
            !-arguments--------------------------------------------------
            real(kind=RP),             intent(in)  :: x(NDIM)
            real(kind=RP),             intent(in)  :: Q(NCONS)
            real(kind=RP),             intent(in)  :: time
            real(kind=RP),             intent(out) :: S(NCONS)
            type(Thermodynamics_t),    intent(in)  :: thermodynamics_
            type(Dimensionless_t),     intent(in)  :: dimensionless_
            type(RefValues_t),         intent(in)  :: refValues_
            
         end subroutine UserDefinedSourceTermNS
#endif
!
!//////////////////////////////////////////////////////////////////////// 
! 
         SUBROUTINE UserDefinedFinalize(mesh, time, iter, maxResidual &
#if defined(NAVIERSTOKES)
                                                    , thermodynamics_ &
                                                    , dimensionless_  &
                                                    , refValues_ & 
#endif   
#if defined(CAHNHILLIARD)
                                                    , multiphase_ &
#endif
                                                    , monitors, &
                                                      elapsedTime, &
                                                      CPUTime   )
!
!           --------------------------------------------------------
!           Called after the solution computed to allow, for example
!           error tests to be performed
!           --------------------------------------------------------
!
            use SMConstants
            use FTAssertions
            USE HexMeshClass
            use PhysicsStorage
            use FluidData
            use MonitorsClass
            IMPLICIT NONE
            CLASS(HexMesh)                        :: mesh
            real(kind=RP)                         :: time
            integer                               :: iter
            real(kind=RP)                         :: maxResidual
#if defined(NAVIERSTOKES)
            type(Thermodynamics_t), intent(in)    :: thermodynamics_
            type(Dimensionless_t),  intent(in)    :: dimensionless_
            type(RefValues_t),      intent(in)    :: refValues_
#endif
#if defined(CAHNHILLIARD)
            type(Multiphase_t),     intent(in)    :: multiphase_
#endif
            type(Monitor_t),        intent(in)    :: monitors
            real(kind=RP),             intent(in) :: elapsedTime
            real(kind=RP),             intent(in) :: CPUTime
!
!           ---------------
!           Local variables
!           ---------------
!
            CHARACTER(LEN=29)                  :: testName           = "Re_tau 5200 standard wall function channel with BR1"
            REAL(KIND=RP)                      :: maxError
            REAL(KIND=RP), ALLOCATABLE         :: QExpected(:,:,:,:)
            INTEGER                            :: eID
            INTEGER                            :: i, j, k, N
            TYPE(FTAssertionsManager), POINTER :: sharedManager
            LOGICAL                            :: success
!
!           -----------------------------------------------------------------------------------------
!           Expected solutions. 
!           Front periodic coupled Back
!           Back periodic coupled Front            
!           Left NoSlipWall
!           Right NoSlipWall
!           Top periodic coupled Bottom
!           Bottom periodic coupled Top
!           use channel true
!           -----------------------------------------------------------------------------------------
!
!
!           ------------------------------------------------
!           Expected Solutions: Turbulent channel flow
!           Number of iterations are for CFL of 0.4, for
!           the roe solver and mach = 0.1
!           ------------------------------------------------
!
#if defined(NAVIERSTOKES)
            INTEGER                            :: iterations(3:7) = [100, 0, 0, 0, 0]
            real(kind=RP), parameter :: residuals(5) = [  0.0000006546358336_RP, &
                                                          0.0624406738429053_RP, &
                                                          0.0000012718216866_RP, &
                                                          0.0000000004842832_RP, &
                                                          0.0009699391660335_RP] 


            CALL initializeSharedAssertionsManager
            sharedManager => sharedAssertionsManager()

            CALL FTAssertEqual(expectedValue = residuals(1)+1.0_RP, &
                               actualValue   = monitors % residuals % values(1,1)+1.0_RP, &
                               tol           = 1.d-10, &
                               msg           = "Continuity residual")

            CALL FTAssertEqual(expectedValue = residuals(2)+1.0_RP, &
                               actualValue   = monitors % residuals % values(2,1)+1.0_RP, &
                               tol           = 1.d-10, &
                               msg           = "X-Momentum residual")

            CALL FTAssertEqual(expectedValue = residuals(3)+1.0_RP, &
                               actualValue   = monitors % residuals % values(3,1)+1.0_RP, &
                               tol           = 1.d-10, &
                               msg           = "Y-Momentum residual")

            CALL FTAssertEqual(expectedValue = residuals(4)+1.0_RP, &
                               actualValue   = monitors % residuals % values(4,1)+1.0_RP, &
                               tol           = 1.d-10, &
                               msg           = "Z-Momentum residual")

            CALL FTAssertEqual(expectedValue = residuals(5)+1.0_RP, &
                               actualValue   = monitors % residuals % values(5,1)+1.0_RP, &
                               tol           = 1.d-10, &
                               msg           = "Energy residual")


            CALL sharedManager % summarizeAssertions(title = testName,iUnit = 6)
   
            IF ( sharedManager % numberOfAssertionFailures() == 0 )     THEN
               WRITE(6,*) testName, " ... Passed"
               WRITE(6,*) "This test case has no expected solution yet, only checks the residual after 100 iterations."
            ELSE
               WRITE(6,*) testName, " ... Failed"
                error stop 99
            END IF 
            WRITE(6,*)
            
            CALL finalizeSharedAssertionsManager
            CALL detachSharedAssertionsManager
#endif

         END SUBROUTINE UserDefinedFinalize
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE UserDefinedTermination
!
!        -----------------------------------------------
!        Called at the the end of the main driver after 
!        everything else is done.
!        -----------------------------------------------
!
         IMPLICIT NONE  
      END SUBROUTINE UserDefinedTermination
      
