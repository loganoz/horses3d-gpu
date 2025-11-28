#include "Includes.h"
module SpatialDiscretization
      use SMConstants
      use HyperbolicDiscretizations
      use EllipticDiscretizations
      use DGIntegrals
      use MeshTypes
      use LESModels
      use HexMeshClass
      use ElementClass
      use PhysicsStorage
      use Physics
      use MPI_Face_Class
      use MPI_Process_Info
      use DGSEMClass
      use FluidData
      use VariableConversion, only: mGradientVariables, GetmOneFluidViscosity,&
                                    GetmTwoFluidsViscosity, chGradientVariables,&
                                    GetCHViscosity
      use BoundaryConditions
      use ProblemFileFunctions, only: UserDefinedSourceTermNS_f
      use ParticlesClass
      use IBMClass
#ifdef _HAS_MPI_
      use mpi
#endif

      private
      public   ComputeTimeDerivative, ComputeTimeDerivativeIsolated, viscousDiscretizationKey
      public   Initialize_SpaceAndTimeMethods, Finalize_SpaceAndTimeMethods

      abstract interface
         SUBROUTINE computeElementInterfaceFluxF(f)
            use FaceClass
            IMPLICIT NONE
            TYPE(Face)   , INTENT(inout) :: f   
         end subroutine computeElementInterfaceFluxF

         SUBROUTINE computeMPIFaceFluxF(f)
            use FaceClass
            IMPLICIT NONE
            TYPE(Face)   , INTENT(inout) :: f   
         end subroutine computeMPIFaceFluxF

         SUBROUTINE computeBoundaryFluxF(f, time)
            use SMConstants
            use FaceClass,  only: Face
            IMPLICIT NONE
            type(Face),    intent(inout) :: f
            REAL(KIND=RP)                :: time
         end subroutine computeBoundaryFluxF
      end interface
      
      character(len=LINE_LENGTH), parameter  :: viscousDiscretizationKey = "viscous discretization"
      character(len=LINE_LENGTH), parameter     :: CHDiscretizationKey = "cahn-hilliard discretization"
      character(len=LINE_LENGTH), parameter  :: FLUID1_COMPRESSIBILITY_KEY = "fluid 1 sound speed square (m/s)"


      real(kind=RP), protected :: IMEX_S0 = 0.0_RP 
      real(kind=RP), protected :: IMEX_K0 = 1.0_RP
      logical                  :: use_non_constant_speed_of_sound = .false.
      integer                  :: speed_of_sound_model = 0
!
!     ========      
      CONTAINS 
!     ========      
!
!////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine Initialize_SpaceAndTimeMethods(controlVariables, mesh)
         use FTValueDictionaryClass
         use Utilities, only: toLower
         use mainKeywordsModule
         use Headers
         use MPI_Process_Info
         implicit none
         class(FTValueDictionary),  intent(in)  :: controlVariables
         class(HexMesh)                         :: mesh
!
!        ---------------
!        Local variables
!        ---------------
!
         character(len=LINE_LENGTH) :: inviscidDiscretizationName
         character(len=LINE_LENGTH) :: viscousDiscretizationName
         character(len=LINE_LENGTH) :: CHDiscretizationName

         if (.not. mesh % child) then ! If this is a child mesh, all these constructs were already initialized for the parent mesh
         
            if ( MPI_Process % isRoot ) then
               write(STD_OUT,'(/)')
               call Section_Header("Spatial discretization scheme")
               write(STD_OUT,'(/)')
            end if
   !
   !        Initialize inviscid discretization
   !        ----------------------------------
            inviscidDiscretizationName = controlVariables % stringValueForKey(inviscidDiscretizationKey,requestedLength = LINE_LENGTH)

            call toLower(inviscidDiscretizationName)
         
            select case ( trim(inviscidDiscretizationName) )
            case ( "standard" )
               if (.not. allocated(HyperbolicDiscretization)) allocate( StandardDG_t  :: HyperbolicDiscretization )

            case ( "split-form")
               print*, "There are no split-forms available for the Multiphase Solver"
               errorMessage(STD_OUT)
               error stop
            case default
               write(STD_OUT,'(A,A,A)') 'Requested inviscid discretization "',trim(inviscidDiscretizationName),'" is not implemented.'
               write(STD_OUT,'(A)') "Implemented discretizations are:"
               write(STD_OUT,'(A)') "  * Standard"
               errorMessage(STD_OUT)
               error stop 

            end select
               
            call HyperbolicDiscretization % Initialize(controlVariables)
   !
   !        Initialize viscous discretization
   !        ---------------------------------         
               if ( .not. controlVariables % ContainsKey(viscousDiscretizationKey) ) then
                  print*, "Input file is missing entry for keyword: viscous discretization"
                  errorMessage(STD_OUT)
                  error stop
               end if

               viscousDiscretizationName = controlVariables % stringValueForKey(viscousDiscretizationKey, requestedLength = LINE_LENGTH)
               call toLower(viscousDiscretizationName)
               
               select case ( trim(viscousDiscretizationName) )
               case("br1")
                  allocate(BassiRebay1_t     :: ViscousDiscretization)

               case("br2")
                  allocate(BassiRebay2_t     :: ViscousDiscretization)

               case("ip")
                  allocate(InteriorPenalty_t :: ViscousDiscretization)

               case default
                  write(STD_OUT,'(A,A,A)') 'Requested viscous discretization "',trim(viscousDiscretizationName),'" is not implemented.'
                  write(STD_OUT,'(A)') "Implemented discretizations are:"
                  write(STD_OUT,'(A)') "  * BR1"
                  write(STD_OUT,'(A)') "  * BR2"
                  write(STD_OUT,'(A)') "  * IP"
                  errorMessage(STD_OUT)
                  error stop 

               end select

               call ViscousDiscretization % Construct(controlVariables, ELLIPTIC_MU)
               call ViscousDiscretization % Describe
               call ViscousDiscretization % CreateDeviceData
!
!           Compute wall distances
!           ----------------------
            call mesh % ComputeWallDistances


!           Initialize models
!           -----------------
            call InitializeLESModel(LESModel, controlVariables)
!
!           Initialize Cahn--Hilliard discretization
!           ----------------------------------------         
            if ( .not. controlVariables % ContainsKey(CHDiscretizationKey) ) then
               print*, "Input file is missing entry for keyword: Cahn-Hilliard discretization"
               errorMessage(STD_OUT)
               error stop
            end if
   
            CHDiscretizationName = controlVariables % stringValueForKey(CHDiscretizationKey, requestedLength = LINE_LENGTH)
            call toLower(CHDiscretizationName)
            
            select case ( trim(CHDiscretizationName) )
            case("br1")
               allocate(BassiRebay1_t     :: CHDiscretization)
   
            case("br2")
               allocate(BassiRebay2_t     :: CHDiscretization)
   
            case("ip")
               allocate(InteriorPenalty_t :: CHDiscretization)
   
            case default
               write(STD_OUT,'(A,A,A)') 'Requested viscous discretization "',trim(CHDiscretizationName),'" is not implemented.'
               write(STD_OUT,'(A)') "Implemented discretizations are:"
               write(STD_OUT,'(A)') "  * BR1"
               write(STD_OUT,'(A)') "  * BR2"
               write(STD_OUT,'(A)') "  * IP"
               errorMessage(STD_OUT)
               error stop 
   
            end select
            
            use_non_constant_speed_of_sound = controlVariables % ContainsKey(FLUID1_COMPRESSIBILITY_KEY)

            if ( controlVariables % ContainsKey("speed of sound profile") .and. (trim(controlVariables % stringValueForKey('speed of sound profile', requestedLength = LINE_LENGTH)) == 'linear quadratic')) then
                speed_of_sound_model = 1
            else
                speed_of_sound_model = 0
            end if 

            call CHDiscretization % Construct(controlVariables, ELLIPTIC_CH)
            call CHDiscretization % Describe
            call CHDiscretization % CreateDeviceData

            if ( .not. MPI_Process % isRoot ) return
            
            if(use_non_constant_speed_of_sound) then
               write(STD_OUT,'(A)') "  Implementing artificial compressibility with a non-constant speed of sound in each phase"
               if (speed_of_sound_model.eq.1) then
                write(STD_OUT,'(A)') "         Speed of sound profile: linear quadratic along interface"
               else 
                write(STD_OUT,'(A)') "         Speed of sound profile: linear along interface"
               end if 
               
            else
               write(STD_OUT,'(A)') "  Implementing artificial compressibility with a constant ACM factor in each phase"
            endif
         
         end if

      end subroutine Initialize_SpaceAndTimeMethods
!
!////////////////////////////////////////////////////////////////////////
!
      subroutine Finalize_SpaceAndTimeMethods
         implicit none
         IF ( ALLOCATED(HyperbolicDiscretization) ) DEALLOCATE( HyperbolicDiscretization )
      end subroutine Finalize_SpaceAndTimeMethods
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE ComputeTimeDerivative( mesh, particles, time, mode, HO_Elements, Level)
         use openacc
         IMPLICIT NONE 
!
!        ---------
!        Arguments
!        ---------
!
         type(HexMesh), target           :: mesh
         type(Particles_t)               :: particles
         REAL(KIND=RP)                   :: time
         integer, intent(in)             :: mode
         logical, intent(in), optional   :: HO_Elements
         integer, intent(in), optional   :: Level
!
!        ---------------
!        Local variables
!        ---------------
!
         INTEGER                 :: k, eID, fID, i, j, locLevel, lID, side
         real(kind=RP)           :: sqrtRho, invMa2, jacobian
         class(Element), pointer :: e
         logical                 :: set_mu   
         real(kind=RP)           :: cs_1, cs_2, Qclip, use_mask, model_mask, invMa2_lin, invMa2_quad, invMa2_const
         real(kind=RP)           :: mu_smag, delta

         if (present(Level)) then
            locLevel = Level
         else
            locLevel = 1
         end if
!$acc wait
!$omp parallel shared(mesh, time, locLevel, mode)       
!
!///////////////////////////////////////////////////
!        1st step: Get chemical potential
!///////////////////////////////////////////////////
!
!        ------------------------------------------
!        Update concentration with the state vector
!        ------------------------------------------
!
         select case (mode)
         case (CTD_IGNORE_MODE,CTD_IMEX_EXPLICIT)
!$omp do schedule(runtime) private(eID, i, j, k)
!$acc parallel loop gang vector_length(128) present(mesh) copyin(locLevel) private(eID) async(1)
            do lID = 1, mesh % MLRK % MLIter(locLevel,8)
               eID = mesh % MLRK % MLIter_eIDN(lID)
               !$acc loop vector collapse(3)
               do k = 0, mesh % elements(eID) % Nxyz(3) ; do j = 0, mesh % elements(eID) % Nxyz(2) ; do i = 0, mesh % elements(eID) % Nxyz(1)
                  mesh % elements(eID) % storage % c(1,i,j,k) = mesh % elements(eID) % storage % Q(IMC,i,j,k)
                  mesh % elements(eID) % storage % S_NS(1,i,j,k) = 0.0_RP                                             ! Assign initial value to source term
                  mesh % elements(eID) % storage % S_NS(2,i,j,k) = 0.0_RP      
                  mesh % elements(eID) % storage % S_NS(3,i,j,k) = 0.0_RP 
                  mesh % elements(eID) % storage % S_NS(4,i,j,k) = 0.0_RP   
                  mesh % elements(eID) % storage % S_NS(5,i,j,k) = 0.0_RP        
                  
               end do               ; end do                ; end do
            end do
!$acc end parallel loop 
!$omp end do
         end select
!
!        -------------------------------
!        Set memory to Cahn-Hilliard (C)
!        -------------------------------
!
!$omp single
         !call mesh % SetStorageToEqn(C_BC)
         select case (mode)
         case (CTD_IGNORE_MODE,CTD_IMEX_EXPLICIT)
            call SetBoundaryConditionsEqn(C_BC)

         case (CTD_IMEX_IMPLICIT,CTD_LAPLACIAN)
            call SetBoundaryConditionsEqn(MU_BC)

         end select
!$omp end single
!
!        --------------------------------------------
!        Prolong Cahn-Hilliard concentration to faces
!        --------------------------------------------
!
         call HexMesh_ProlongSolToFaces(mesh, NCOMP, Level=locLevel)
!
!        ----------------
!        Update MPI Faces
!        ----------------
!
#ifdef _HAS_MPI_
!$omp single
         call HexMesh_UpdateMPIFacesSolution(mesh, NCOMP)
         call HexMesh_GatherMPIFacesSolution(mesh, NCOMP)   
!$omp end single
#endif

!$omp do schedule(runtime) private(fID, i, j)
         !$acc parallel loop gang vector_length(32) present(mesh) copyin(locLevel) private(fID) async(1)
         do lID = 1, mesh % MLRK % MLIter(locLevel,2)
            fID = mesh % MLRK % MLIter_fID(lID)
            !$acc loop vector collapse(2)
            do j = 0, mesh % faces(fID) % Nf(2)  ; do i = 0, mesh % faces(fID) % Nf(1)
               mesh % faces(fID) % storage(1) % c(1,i,j) = mesh % faces(fID) % storage(1) % Q(IMC,i,j)
               mesh % faces(fID) % storage(2) % c(1,i,j) = mesh % faces(fID) % storage(2) % Q(IMC,i,j)
            end do               ; end do 
         end do
         !$acc end parallel loop
!$omp end do
! !        ----------------
! !        Update MPI Faces
! !        ----------------
! !
! #ifdef _HAS_MPI_
! !$omp single
         ! call HexMesh_UpdateMPIFacesSolution(mesh, NCOMP)
! !$omp end single
! #endif
!
!        ------------------------------------------------------------
!        Get concentration (lifted) gradients (also prolong to faces)
!        ------------------------------------------------------------
!
!$acc wait
         set_mu = .false.
         call HexMesh_ComputeLocalGradientCH(mesh, set_mu, Level=locLevel)

         ! This is chGradientVariables in master but it is not used here - dummy input
         call CHDiscretization % ComputeGradient(NCOMP, NCOMP, mesh, time, mGradientVariables, Level = locLevel) 
!
!        --------------------
!        Update MPI Gradients
!        --------------------
!
!$acc wait
#ifdef _HAS_MPI_
!$omp single
         call HexMesh_UpdateMPIFacesGradients(mesh, NCOMP)
         call HexMesh_GatherMPIFacesGradients(mesh, NCOMP)
!$omp end single
#endif
!
!        ----------------------
!        Get chemical potential
!        ----------------------
!
!        Get the concentration Laplacian (into QDot => cDot)
         call ComputeLaplacian(mesh, time, Level=locLevel)

         select case (mode)
         case (CTD_IGNORE_MODE, CTD_IMEX_EXPLICIT)
!$omp do schedule(runtime) private(eID, i, j, k)
            !$acc parallel loop gang vector_length(128) present(mesh, multiphase) copyin(locLevel) private(eID) async(1)
            do lID = 1, mesh % MLRK % MLIter(locLevel,1)
               eID = mesh % MLRK % MLIter_eID(lID)
!
               !$acc loop vector collapse(3)
               do k = 0, mesh % elements(eID) % Nxyz(3) ; do j = 0, mesh % elements(eID) % Nxyz(2) ; do i = 0, mesh % elements(eID) % Nxyz(1)
!              + Linear part
                  !The 1st index is only one because the equation is set to CH
                  !mesh % elements(eID) % storage % mu = - POW2(multiphase % eps)* mesh % elements(eID) % storage % QDot
                  mesh % elements(eID) % storage % mu(1,i,j,k) = - 1.5_RP * multiphase % eps * multiphase % sigma * mesh % elements(eID) % storage % QDot(1,i,j,k)
!              + NonLinear part
                  !call AddQuarticDWPDerivative(mesh % elements(eID) % storage % c, mesh % elements(eID) % storage % mu)
                  call Multiphase_AddChemFEDerivative(mesh % elements(eID) % storage % c(1,i,j,k), mesh % elements(eID) % storage % mu(1,i,j,k))
               end do               ; end do                ; end do
!
            end do
            !$acc end parallel loop 
!$omp end do    
   
         case (CTD_IMEX_IMPLICIT)
!$omp do schedule(runtime) private(i, j, k)
            !$acc parallel loop gang vector_length(128) present(mesh, multiphase) async(1)
            do eID = 1, size(mesh % elements)
               !$acc loop vector collapse(3)
               do k = 0, mesh % elements(eID) % Nxyz(3) ; do j = 0, mesh % elements(eID) % Nxyz(2) ; do i = 0, mesh % elements(eID) % Nxyz(1)
!              + Linear part
                  !mesh % elements(eID) % storage % mu = - IMEX_K0 * POW2(multiphase % eps) * mesh % elements(eID) % storage % QDot &
                  mesh % elements(eID) % storage % mu(1,i,j,k) = - 1.5_RP * IMEX_K0 * multiphase % eps * multiphase % sigma * mesh % elements(eID) % storage % QDot(1,i,j,k) &
                                                     + IMEX_S0 * mesh % elements(eID) % storage % c(1,i,j,k)
!              + Multiply by mobility
                  mesh % elements(eID) % storage % mu(1,i,j,k) = multiphase % M0 * mesh % elements(eID) % storage % mu(1,i,j,k)

               end do               ; end do                ; end do
            end do
            !$acc end parallel loop 
!$omp end do         
         end select
!
!        -----------------------------------
!        Prolong chemical potential to faces
!        -----------------------------------
!
         select case(mode)
         case(CTD_LAPLACIAN)
         case default
!!$omp single
         !   call mesh % SetStorageToEqn(MU_BC)
!!$omp end single
            ! copy mu to Q(1)
!$omp do schedule(runtime) private(eID, i, j, k)
            !$acc parallel loop gang vector_length(128) present(mesh, mesh % elements) copyin(locLevel) private(eID) async(1)
            do lID = 1, mesh % MLRK % MLIter(locLevel,8)
               eID = mesh % MLRK % MLIter_eIDN(lID)
               !$acc loop vector collapse(3)
               do k = 0, mesh % elements(eID) % Nxyz(3) ; do j = 0, mesh % elements(eID) % Nxyz(2) ; do i = 0, mesh % elements(eID) % Nxyz(1)
                  mesh % elements(eID) % storage % Q(1,i,j,k) = mesh % elements(eID) % storage % mu(1,i,j,k) 
               end do               ; end do                ; end do
            end do
            !$acc end parallel loop 
!$omp end do         

            call HexMesh_ProlongSolToFaces(mesh, NCOMP, Level=locLevel)
!
!         ----------------
!           Update MPI Faces
!           ----------------
!
!$acc wait
#ifdef _HAS_MPI_
!$omp single
         call HexMesh_UpdateMPIFacesSolution(mesh, NCOMP)
         call HexMesh_GatherMPIFacesSolution(mesh, NCOMP)   
!$omp end single
#endif

            !!! copy c to Q(1) &&&  copy Q faces to mu faces !!!!
!$omp do schedule(runtime) private(eID, i, j, k)
            !$acc parallel loop gang vector_length(128) present(mesh) copyin(locLevel) private(eID) async(1)
            do lID = 1, mesh % MLRK % MLIter(locLevel,8)
               eID = mesh % MLRK % MLIter_eIDN(lID)
               !$acc loop vector collapse(3)
               do k = 0, mesh % elements(eID) % Nxyz(3) ; do j = 0, mesh % elements(eID) % Nxyz(2) ; do i = 0, mesh % elements(eID) % Nxyz(1)
                  mesh % elements(eID) % storage % Q(1,i,j,k) = mesh % elements(eID) % storage % c(1,i,j,k) 
               end do               ; end do                ; end do
            end do
            !$acc end parallel loop 
!$omp end do  

!$omp do schedule(runtime) private(fID, i, j)
         !$acc parallel loop gang present(mesh) copyin(locLevel) private(fID) async(1)
         do lID = 1, mesh % MLRK % MLIter(locLevel,2)
            fID = mesh % MLRK % MLIter_fID(lID)
            !$acc loop vector collapse(2)
            do j = 0, mesh % faces(fID) % Nf(2)  ; do i = 0, mesh % faces(fID) % Nf(1)
               mesh % faces(fID) % storage(1) % mu(1,i,j) = mesh % faces(fID) % storage(1) % Q(1,i,j)
               mesh % faces(fID) % storage(2) % mu(1,i,j) = mesh % faces(fID) % storage(2) % Q(1,i,j)
               
               mesh % faces(fID) % storage(1) % Q(1,i,j) = mesh % faces(fID) % storage(1) % c(1,i,j)
               mesh % faces(fID) % storage(2) % Q(1,i,j) = mesh % faces(fID) % storage(2) % c(1,i,j)
            end do               ; end do 
         end do
         !$acc end parallel loop
!$omp end do
         end select
!
!/////////////////////////////////////////////////////////////////////////////////
!        2nd step: If IMEX_IMPLCIIT, get the chemical potential laplacian and exit
!/////////////////////////////////////////////////////////////////////////////////
!
         select case (mode)
         case (CTD_IMEX_IMPLICIT)
!
!           ------------------------------------------------------------
!           Get concentration (lifted) gradients (also prolong to faces)
!           ------------------------------------------------------------
!
            !set_mu is always false for CH
            set_mu = .false.
            call HexMesh_ComputeLocalGradientCH(mesh, set_mu)
            ! This is chGradientVariables in master but it is not used here - dummy input
            call CHDiscretization % ComputeGradient(NCOMP, NCOMP, mesh, time, mGradientVariables)
!
!           --------------------
!           Update MPI Gradients
!           --------------------
!
#ifdef _HAS_MPI_
!$omp single
            call HexMesh_UpdateMPIFacesGradients(mesh, NCOMP)
!$omp end single
#endif
!
!           ----------------------
!           Get chemical potential
!           ----------------------
!
!           Get the concentration Laplacian (into QDot => cDot)

            call ComputeLaplacian(mesh, time)
!
!           ------------------------------------------
!           *** WARNING! The storage leaves set to CH!
!           ------------------------------------------
!
!$omp single
            call mesh % SetStorageToEqn(C_BC)
            call SetBoundaryConditionsEqn(C_BC)
!$omp end single
         end select
!
!///////////////////////////////////////////////
!        3rd step: Navier-Stokes time derivative
!///////////////////////////////////////////////
!
         select case (mode)
         case (CTD_IGNORE_MODE, CTD_IMEX_EXPLICIT)
!$omp single         
         !call mesh % SetStorageToEqn(NS_BC)
         !$acc wait
         call SetBoundaryConditionsEqn(NS_BC)
!$omp end single
!
!        -------------------------
!        Prolong solution to faces        
!        -------------------------
!

         call HexMesh_ProlongSolToFaces(mesh, NCONS, Level=locLevel)
!
!        ----------------
!        Update MPI Faces
!        ----------------
!
!$acc wait
#ifdef _HAS_MPI_
!$omp single
         call HexMesh_UpdateMPIFacesSolution(mesh, NCONS)
         call HexMesh_GatherMPIFacesSolution(mesh, NCONS)
!$omp end single
#endif
!
!        -------------------------------------
!        Get the density and invMa2 in faces and elements
!        -------------------------------------
!
        ! Compute masks as 0.0_RP or 1.0_RP
         use_mask = merge(1.0_RP, 0.0_RP, use_non_constant_speed_of_sound)
         model_mask = merge(1.0_RP, 0.0_RP, speed_of_sound_model.eq.1)
         
         cs_1 = sqrt(dimensionless % invMa2(1) / dimensionless % rho(1))
         cs_2 = sqrt(dimensionless % invMa2(2) / dimensionless % rho(2))

!$omp do schedule(runtime) private(eID, Qclip, invMa2_lin, invMa2_quad, invMa2_const, i, j, k)
         !$acc parallel loop gang vector_length(128) present(mesh, dimensionless, mesh % elements, mesh % MLRK) copyin(locLevel, use_mask, model_mask, cs_1, cs_2) private(eID) async(1)
         do lID = 1, mesh % MLRK % MLIter(locLevel,1)
            eID = mesh % MLRK % MLIter_eID(lID)
            !$acc loop vector collapse(3) private(Qclip, invMa2_lin, invMa2_quad, invMa2_const)
            do k = 0, mesh % elements(eID) % Nxyz(3) ; do j = 0, mesh % elements(eID) % Nxyz(2) ; do i = 0, mesh % elements(eID) % Nxyz(1)
            
                ! Clip Q between 0 and 1
                Qclip = min(max(mesh % elements(eID) % storage % Q(IMC,i,j,k), 0.0_RP), 1.0_RP)
                ! Rho
                mesh % elements(eID) % storage % rho(i,j,k) = dimensionless % rho(2) + (dimensionless % rho(1)-dimensionless % rho(2))*Qclip
                ! Linear model
                invMa2_lin = ( cs_1 * Qclip + cs_2 * (1.0_RP - Qclip) )**2
                invMa2_lin = invMa2_lin * mesh % elements(eID) % storage % rho(i,j,k)
                ! Linear-quadratic model
                invMa2_quad = dimensionless % invMa2(2) + (dimensionless % invMa2(1) - dimensionless % invMa2(2)) * Qclip
                ! Constant model
                invMa2_const = dimensionless % invMa2(1)
                ! Blend them arithmetically
                mesh % elements(eID) % storage % invMa2 (i,j,k) = &
                    use_mask * ( (1.0_RP - model_mask) * invMa2_lin + model_mask * invMa2_quad ) + &
                    (1.0_RP - use_mask) * invMa2_const
                    
            end do ; end do ; end do
         end do
   !$acc end parallel loop
!$omp end do nowait

!$omp do schedule(runtime) private(fID, Qclip, invMa2_lin, invMa2_quad, invMa2_const, i, j, side)
         !$acc parallel loop gang vector_length(64) present(mesh, mesh % faces, dimensionless) copyin(locLevel, use_mask, model_mask, cs_1, cs_2) private(fID) async(1)
         do lID = 1, mesh % MLRK % MLIter(locLevel,2)
            fID = mesh % MLRK % MLIter_fID(lID)
            !$acc loop vector collapse(2) private(Qclip, invMa2_lin, invMa2_quad, invMa2_const, side)
            do j = 0, mesh % faces(fID) % Nf(2)  ; do i = 0, mesh % faces(fID) % Nf(1) 
            
                !$acc loop seq 
                do side = 1, 2
                    ! Clip Q between 0 and 1
                    Qclip = min(max(mesh % faces(fID) % storage(side) % Q(IMC,i,j), 0.0_RP), 1.0_RP)
                    ! Rho
                    mesh % faces(fID) % storage(side) % rho(i,j) = dimensionless % rho(2) + (dimensionless % rho(1)-dimensionless % rho(2))*Qclip
                    ! Linear model
                    invMa2_lin = ( cs_1 * Qclip + cs_2 * (1.0_RP - Qclip) )**2
                    invMa2_lin = invMa2_lin * mesh % faces(fID) % storage(side) % rho(i,j)
                    ! Linear-quadratic model
                    invMa2_quad = dimensionless % invMa2(2) + &
                                  (dimensionless % invMa2(1) - dimensionless % invMa2(2)) * Qclip
                    ! Constant model
                    invMa2_const = dimensionless % invMa2(1)
                    ! Blend them arithmetically
                    mesh % faces(fID) % storage(side) % invMa2 (i,j) = &
                        use_mask * ( (1.0_RP - model_mask) * invMa2_lin + model_mask * invMa2_quad ) + &
                        (1.0_RP - use_mask) * invMa2_const
                end do
            end do ; end do 
         end do 
         !$acc end parallel loop
!$omp end do
!$acc wait
!
!        ----------------------------------------
!        Compute local entropy variables gradient
!        ----------------------------------------
!
         !set_mu is always true for MU
         set_mu = .true.
         call HexMesh_ComputeLocalGradientMU(mesh, set_mu, Level=locLevel)
!
!        -------------------------------------
!        Add the Non-Conservative term to QDot
!        -------------------------------------
!
!$omp do schedule(runtime) private(i,j,k,sqrtRho,invMa2,eID,jacobian)
!$acc parallel loop gang vector_length(128) present(mesh, mesh % elements) copyin(locLevel) private(eID) async(1)
         do lID = 1, mesh % MLRK % MLIter(locLevel,1) ! 
            eID = mesh % MLRK % MLIter_eID(lID)
            !$acc loop vector collapse(3) private(jacobian,sqrtRho,invMa2)
            do k = 0, mesh % elements(eID) % Nxyz(3) ; do j = 0, mesh % elements(eID) % Nxyz(2) ; do i = 0, mesh % elements(eID) % Nxyz(1)
               jacobian = mesh % elements(eID) % geom % jacobian(i,j,k)
               sqrtRho = sqrt(mesh % elements(eID) % storage % rho(i,j,k))
               invMa2 = mesh % elements(eID) % storage % invMa2(i,j,k)
               mesh % elements(eID) % storage % QDot(IMC,i,j,k)      = 0.0_RP
               mesh % elements(eID) % storage % QDot(IMSQRHOU,i,j,k) = (-0.5_RP*sqrtRho*( mesh % elements(eID) % storage % Q(IMSQRHOU,i,j,k)*mesh % elements(eID) % storage % U_x(IGU,i,j,k) & 
                                                                                        + mesh % elements(eID) % storage % Q(IMSQRHOV,i,j,k)*mesh % elements(eID) % storage % U_y(IGU,i,j,k) &   
                                                                                        + mesh % elements(eID) % storage % Q(IMSQRHOW,i,j,k)*mesh % elements(eID) % storage % U_z(IGU,i,j,k) ) &
                                                                                        - mesh % elements(eID) % storage % Q(IMC,i,j,k)*mesh % elements(eID) % storage % U_x(IGMU,i,j,k) )* jacobian

               mesh % elements(eID) % storage % QDot(IMSQRHOV,i,j,k) = (-0.5_RP*sqrtRho*( mesh % elements(eID) % storage % Q(IMSQRHOU,i,j,k)*mesh % elements(eID) % storage % U_x(IGV,i,j,k) & 
                                                                                        + mesh % elements(eID) % storage % Q(IMSQRHOV,i,j,k)*mesh % elements(eID) % storage % U_y(IGV,i,j,k) &   
                                                                                        + mesh % elements(eID) % storage % Q(IMSQRHOW,i,j,k)*mesh % elements(eID) % storage % U_z(IGV,i,j,k) ) &
                                                                                        - mesh % elements(eID) % storage % Q(IMC,i,j,k)*mesh % elements(eID) % storage % U_y(IGMU,i,j,k)) * jacobian

               mesh % elements(eID) % storage % QDot(IMSQRHOW,i,j,k) = (-0.5_RP*sqrtRho*( mesh % elements(eID) % storage % Q(IMSQRHOU,i,j,k)*mesh % elements(eID) % storage % U_x(IGW,i,j,k) & 
                                                                                        + mesh % elements(eID) % storage % Q(IMSQRHOV,i,j,k)*mesh % elements(eID) % storage % U_y(IGW,i,j,k) &   
                                                                                        + mesh % elements(eID) % storage % Q(IMSQRHOW,i,j,k)*mesh % elements(eID) % storage % U_z(IGW,i,j,k) ) &
                                                                                        - mesh % elements(eID) % storage % Q(IMC,i,j,k)*mesh % elements(eID) % storage % U_z(IGMU,i,j,k)) * jacobian
               
               mesh % elements(eID) % storage % QDot(IMP,i,j,k) = -invMa2*(  mesh % elements(eID) % storage % U_x(IGU,i,j,k) + mesh % elements(eID) % storage % U_y(IGV,i,j,k) &
                                                                                       + mesh % elements(eID) % storage % U_z(IGW,i,j,k)) * jacobian

            end do                ; end do                ; end do
         end do
!$acc end parallel loop
!$omp end do
!$acc wait
         call ViscousDiscretization % LiftGradients( NCONS, NCONS, mesh , time , mGradientVariables)
!$acc wait
#ifdef _HAS_MPI_
!$omp single
         call HexMesh_UpdateMPIFacesGradients(mesh,NCONS)
!$omp end single
#endif  
!
!        -----------------------
!        Compute time derivative
!        -----------------------
!
         select case (mode)
         case(CTD_IMEX_EXPLICIT)
            call multiphase % SetStarMobility(0.0_RP)
         case(CTD_IGNORE_MODE)
            call multiphase % SetStarMobility(multiphase % M0)
         end select
!$acc wait
            call ComputeNSTimeDerivative(mesh, time, Level=locLevel)

            call multiphase % SetStarMobility(multiphase % M0)
!$acc wait
         end select
!
!        -------------------------------------------------------------------------------
!        If IMEX_Explicit, compute cDot with the explicit part of the chemical potential
!        -------------------------------------------------------------------------------
!
         select case (mode)
         case(CTD_IMEX_EXPLICIT)
!$omp do schedule(runtime) private(i, j, k)
            !$acc parallel loop gang vector_length(128) present(mesh, multiphase) async(1)
            do eID = 1, size(mesh % elements)
               !$acc loop vector collapse(3)
               do k = 0, mesh % elements(eID) % Nxyz(3) ; do j = 0, mesh % elements(eID) % Nxyz(2) ; do i = 0, mesh % elements(eID) % Nxyz(1)
!                 + Linear part
                  mesh % elements(eID) % storage % mu(1,i,j,k) = - IMEX_S0 * mesh % elements(eID) % storage % c(1,i,j,k) &
                                                        - 1.5_RP*(1.0_RP - IMEX_K0)*multiphase % sigma*multiphase % eps*mesh % elements(eID) % storage % cDot(1,i,j,k)
                  !mesh % elements(eID) % storage % mu = - IMEX_S0 * mesh % elements(eID) % storage % c &
                  !                                      - (1.0_RP - IMEX_K0)* POW2(multiphase % eps)*mesh % elements(eID) % storage % cDot
!
!                 + NonLinear part
                  !call AddQuarticDWPDerivative(mesh % elements(eID) % storage % c, mesh % elements(eID) % storage % mu)
                  call Multiphase_AddChemFEDerivative(mesh % elements(eID) % storage % c(1,i,j,k), mesh % elements(eID) % storage % mu(1,i,j,k))
               end do               ; end do                ; end do
            end do
            !$acc end parallel loop
!$omp end do         
!
!           -----------------------------------
!           Prolong chemical potential to faces
!           -----------------------------------
!
!$omp single
            call mesh % SetStorageToEqn(MU_BC)
            call SetBoundaryConditionsEqn(MU_BC)
!$omp end single
            call HexMesh_ProlongSolToFaces(mesh, NCOMP)
!
!           ------------------------------------------------------------
!           Get concentration (lifted) gradients (also prolong to faces)
!           ------------------------------------------------------------
!
            !set_mu is always false for CH
            set_mu = .false.
            call HexMesh_ComputeLocalGradientCH(mesh, set_mu)
            ! This is chGradientVariables in master but it is not used here - dummy input
            call CHDiscretization % ComputeGradient(NCOMP, NCOMP, mesh, time, mGradientVariables)
!           --------------------------------
!           Get chemical potential laplacian
!           --------------------------------
!
!           Get the concentration Laplacian (into QDot => cDot)
            call ComputeLaplacian(mesh, time)
!
!$omp single
            call mesh % SetStorageToEqn(NS_BC)
            call SetBoundaryConditionsEqn(NS_BC)
!$omp end single
!
!           -----------------------------------------
!           Add the Chemical potential to the NS QDot
!           -----------------------------------------
!
!$omp do schedule(runtime)private(i, j, k)
            !$acc parallel loop gang vector_length(128) present(mesh, mesh % elements) async(1)
            do eID = 1, size(mesh % elements)
               !$acc loop vector collapse(3)
               do k = 0, mesh % elements(eID) % Nxyz(3) ; do j = 0, mesh % elements(eID) % Nxyz(2) ; do i = 0, mesh % elements(eID) % Nxyz(1)
                  mesh % elements(eID) % storage % QDot(IMC,i,j,k) = mesh % elements(eID) % storage % QDot(IMC,i,j,k) &
                                                                                   + multiphase % M0*mesh % elements(eID) % storage % cDot(1,i,j,k)
               end do               ; end do                ; end do
            end do
            !$acc end parallel loop
!$omp end do
         end select
!$omp end parallel
!

!        TODO: This is not useful - leaving here for debug
         !$acc wait

      END SUBROUTINE ComputeTimeDerivative
!
!////////////////////////////////////////////////////////////////////////////////////
!
!           Navier--Stokes procedures
!           -------------------------
!
!////////////////////////////////////////////////////////////////////////////////////
!
      subroutine ComputeNSTimeDerivative( mesh , t, level )
         use SpongeClass, only: sponge, addSourceSponge
         use ActuatorLine, only: farm, ForcesFarm
         use AcousticSourceClass, only: AcousticSource, addSourceAcoustic
         implicit none
         type(HexMesh)              :: mesh
         real(kind=RP)              :: t
         integer, intent(in), optional   :: Level
         procedure(UserDefinedSourceTermNS_f) :: UserDefinedSourceTermNS
!
!        ---------------
!        Local variables
!        ---------------
         integer     :: eID , i, j, k, ierr, locLevel,lID
         integer     :: fID, side, iFace, iEl, eq
         real(kind=RP) :: sqrtRho, invSqrtRho, delta, factor, minvalC, maxvalC
         real(kind=RP) :: Source(NCONS), prod
         
         if (present(Level)) then
            locLevel = Level
         else
            locLevel = 1
         end if
         
         if ( LESModel % active) then
!$omp do schedule(runtime) private(i,j,k,delta,eID,prod)
            !$acc parallel loop gang present(mesh, LESModel) vector_length(128) copyin(locLevel) private(delta,eID,prod) async(1)
            do lID = 1, mesh % MLRK % MLIter(locLevel,1)
               eID = mesh % MLRK % MLIter_eID(lID)
			   prod  = (mesh % elements(eID) % Nxyz(1) + 1.0_RP)*(mesh % elements(eID) % Nxyz(2) + 1.0_RP)*(mesh % elements(eID) % Nxyz(3) + 1.0_RP)
               delta = (mesh % elements(eID) % geom % Volume / prod) ** (1.0_RP / 3.0_RP)
               
               !$acc loop vector collapse(3) 
               do k = 0, mesh % elements(eID) % Nxyz(3) ; do j = 0, mesh % elements(eID) % Nxyz(2) ; do i = 0, mesh % elements(eID) % Nxyz(1)
                  call LESModel_Selector(LESModel, delta, mesh % elements(eID) % geom % dWall(i,j,k), &
                                                          mesh % elements(eID) % storage % Q(:,i,j,k),   &
                                                          mesh % elements(eID) % storage % U_x(:,i,j,k), &
                                                          mesh % elements(eID) % storage % U_y(:,i,j,k), &
                                                          mesh % elements(eID) % storage % U_z(:,i,j,k), &
                                                          mesh % elements(eID) % storage % mu_turb_NS(i,j,k) )

                  mesh % elements(eID) % storage % mu_NS(1,i,j,k) = mesh % elements(eID) % storage % mu_turb_NS(i,j,k)
               end do                ; end do                ; end do
            end do
            !$acc end parallel loop
!$omp end do 
      end if
!
!        Compute viscosity at interior and boundary faces
!        ------------------------------------------------
         call compute_viscosity_at_faces(mesh % MLRK % MLIter(locLevel,3), 2, mesh % MLRK % MLIter_fID_Interior(1:mesh % MLRK % MLIter(locLevel,3)), mesh)
         call compute_viscosity_at_faces(mesh % MLRK % MLIter(locLevel,4), 1, mesh % MLRK % MLIter_fID_Boundary(1:mesh % MLRK % MLIter(locLevel,4)), mesh)
!$acc wait                
!
!        ****************
!        Volume integrals
!        ****************
!
         call TimeDerivative_VolumetricContribution(mesh, Level=locLevel)
!
!        ******************************************
!        Compute Riemann solver of non-shared faces
!        ******************************************
!
!$omp do schedule(runtime) private(fID, side)
!$acc parallel loop gang collapse(2) vector_length(64) present(mesh, multiphase) copyin(locLevel) private(fID) async(1)
         do iFace = 1, mesh % MLRK % MLIter(locLevel,3) ; do side = 1,2
            fID = mesh % MLRK % MLIter_fID_Interior(iFace)
            call computeElementInterfaceFlux_MUviscous(mesh % faces(fID), side)
         end do ; end do 
!$acc end parallel loop
!$omp end do
!$acc wait
         call computeElementInterfaceFlux_MU(mesh, Level=locLevel)
!$acc wait
         call computeBoundaryFlux_MU(mesh,t)
!$acc wait
! !
! !        *************************************************************************************
! !        Element without shared faces: Surface integrals, scaling of elements with Jacobian, 
! !                                      sqrt(rho), and add source terms
! !        *************************************************************************************
! !
!$omp do schedule(runtime) private(i,j,k,sqrtRho,invSqrtRho, eID)
         !$acc parallel loop gang num_gangs(mesh % MLRK % MLIter(locLevel,5)) vector_length(128) present(mesh, dimensionless) copyin(locLevel) private(eID) async(1)
         do iEl = 1, mesh % MLRK % MLIter(locLevel,5)
            eID = mesh % MLRK % MLIter_eID_Seq(iEl)
            
            call TimeDerivative_FacesContribution(mesh % elements(eID), mesh)
            
            !$acc loop vector collapse(3) private(sqrtRho,invSqrtRho)
            do k = 0, mesh % elements(eID) % Nxyz(3) ; do j = 0, mesh % elements(eID) % Nxyz(2) ; do i = 0, mesh % elements(eID) % Nxyz(1) 
               sqrtRho = sqrt(mesh % elements(eID) % storage % rho(i,j,k))
               invSqrtRho = 1.0_RP / sqrtRho
!
!           ++ Scale with sqrt(Rho) and add gravity
               mesh % elements(eID) % storage % QDot(IMSQRHOU,i,j,k) = mesh % elements(eID) % storage % QDot(IMSQRHOU,i,j,k) * invSqrtRho
               mesh % elements(eID) % storage % QDot(IMSQRHOU,i,j,k) = mesh % elements(eID) % storage % QDot(IMSQRHOU,i,j,k) & 
                                                                              + sqrtRho * dimensionless % invFr2 * dimensionless % gravity_dir(1)
               mesh % elements(eID) % storage % QDot(IMSQRHOV,i,j,k) = mesh % elements(eID) % storage % QDot(IMSQRHOV,i,j,k) * invSqrtRho
               mesh % elements(eID) % storage % QDot(IMSQRHOV,i,j,k) = mesh % elements(eID) % storage % QDot(IMSQRHOV,i,j,k) & 
                                                                              + sqrtRho * dimensionless % invFr2 * dimensionless % gravity_dir(2)
               mesh % elements(eID) % storage % QDot(IMSQRHOW,i,j,k) = mesh % elements(eID) % storage % QDot(IMSQRHOW,i,j,k) * invSqrtRho
               mesh % elements(eID) % storage % QDot(IMSQRHOW,i,j,k) = mesh % elements(eID) % storage % QDot(IMSQRHOW,i,j,k) & 
                                                                              + sqrtRho * dimensionless % invFr2 * dimensionless % gravity_dir(3)

            end do         ; end do          ; end do 
         end do
         !$acc end parallel loop
!$omp end do
! !
! !        ******************************************
! !        Do the same for elements with shared faces
! !        ******************************************
! !
!$acc wait
#ifdef _HAS_MPI_
         if ( MPI_Process % doMPIAction ) then
!$omp single
             call HexMesh_GatherMPIFacesGradients(mesh, NCONS)
!$omp end single
! !
! !           Compute viscosity at MPI faces
! !           ------------------------------
             call compute_viscosity_at_faces(mesh % MLRK % MLIter(locLevel,7), 2, mesh % MLRK % MLIter_fID_MPI(1:mesh % MLRK % MLIter(locLevel,7)), mesh)
! !           **************************************
! !           Compute Riemann solver of shared faces
! !           **************************************
! !
!$omp do schedule(runtime) private(fID)
!$acc parallel loop gang vector_length(64) present(mesh, multiphase, ViscousDiscretization) copyin(locLevel) private(fID) async(1)
            do iFace = 1, mesh % MLRK % MLIter(locLevel,7)
               fID = mesh % MLRK % MLIter_fID_MPI(iFace)
               CALL computeMPIFaceFlux_MU( mesh % faces(fID) )
            end do
!$acc end parallel loop
!$omp end do 
!$acc wait
! !
! !           ***********************************************************
! !           Surface integrals and scaling of elements with shared faces
! !           ***********************************************************
! ! 
!$omp do schedule(runtime) private(i,j,k, eID, sqrtRho, invSqrtRho)
!$acc parallel loop gang vector_length(128) present(mesh,dimensionless) copyin(locLevel) private(eID) async(1)
            do iEl = 1, mesh % MLRK % MLIter(locLevel,6)
               eID = mesh % MLRK % MLIter_eID_MPI(iEl)
               
               call TimeDerivative_FacesContribution(mesh % elements(eID), mesh)
               
            !$acc loop vector collapse(3) private(sqrtRho,invSqrtRho)
                do k = 0, mesh % elements(eID) % Nxyz(3) ; do j = 0, mesh % elements(eID) % Nxyz(2) ; do i = 0, mesh % elements(eID) % Nxyz(1) 
                   sqrtRho = sqrt(mesh % elements(eID) % storage % rho(i,j,k))
                   invSqrtRho = 1.0_RP / sqrtRho
    !
    !           ++ Scale with sqrt(Rho) and add gravity
                   mesh % elements(eID) % storage % QDot(IMSQRHOU,i,j,k) = mesh % elements(eID) % storage % QDot(IMSQRHOU,i,j,k) * invSqrtRho
                   mesh % elements(eID) % storage % QDot(IMSQRHOU,i,j,k) = mesh % elements(eID) % storage % QDot(IMSQRHOU,i,j,k) & 
                                                                                  + sqrtRho * dimensionless % invFr2 * dimensionless % gravity_dir(1)
                   mesh % elements(eID) % storage % QDot(IMSQRHOV,i,j,k) = mesh % elements(eID) % storage % QDot(IMSQRHOV,i,j,k) * invSqrtRho
                   mesh % elements(eID) % storage % QDot(IMSQRHOV,i,j,k) = mesh % elements(eID) % storage % QDot(IMSQRHOV,i,j,k) & 
                                                                                  + sqrtRho * dimensionless % invFr2 * dimensionless % gravity_dir(2)
                   mesh % elements(eID) % storage % QDot(IMSQRHOW,i,j,k) = mesh % elements(eID) % storage % QDot(IMSQRHOW,i,j,k) * invSqrtRho
                   mesh % elements(eID) % storage % QDot(IMSQRHOW,i,j,k) = mesh % elements(eID) % storage % QDot(IMSQRHOW,i,j,k) & 
                                                                                  + sqrtRho * dimensionless % invFr2 * dimensionless % gravity_dir(3)
               end do         ; end do          ; end do 
            end do
!$omp end do
!$acc wait
!
!           Add an MPI Barrier
!           ------------------
!$omp single
            call mpi_barrier(MPI_COMM_WORLD, ierr)
!$omp end single
         end if
#endif
!           ***************
!           Add source term
!           ***************


!           User-defined source term for CPU 
#ifndef _OPENACC
!$omp do schedule(runtime) private(i,j,k, invSqrtRho, eID)
         do lID = 1, mesh % MLRK % MLIter(locLevel,1)
            eID = mesh % MLRK % MLIter_eID(lID)
            
            do k = 0, mesh % elements(eID) % Nxyz(3)   ; do j = 0, mesh % elements(eID) % Nxyz(2) ; do i = 0, mesh % elements(eID) % Nxyz(1) 
                invSqrtRho = 1.0_RP / sqrt(mesh % elements(eID) % storage % rho(i,j,k))
                call UserDefinedSourceTermNS(mesh % elements(eID) % geom % x(:,i,j,k), mesh % elements(eID) % storage % Q(:,i,j,k), t, mesh % elements(eID) % storage % S_NS(:,i,j,k), thermodynamics, dimensionless, refValues, multiphase)
            end do   ;  end do   ;  end do   
         end do
!$omp end do
#endif

!The scale with sqrtRho is done in the subroutines, not done againg here
         !$acc wait
         call addSourceSponge(sponge,mesh)
         !$acc wait
         call CompilerDefinedSourceTerm(mesh, t)
         !$acc wait
         call ForcesFarm(farm, mesh, t, Level=locLevel) 
         !$acc wait
         call addSourceAcoustic(AcousticSource, mesh, t)
         !$acc wait 
!
!        ****************************
!        Now add all the source terms
!        ****************************
!$omp do schedule(runtime) private(i,j,k,eID,eq)
!$acc parallel loop gang vector_length(128) present(mesh) copyin(locLevel) private(eID) async(1)
         do lID = 1, mesh % MLRK % MLIter(locLevel,1)
            eID = mesh % MLRK % MLIter_eID(lID)
            !$acc loop vector collapse(4)
            do k = 0, mesh % elements(eID) % Nxyz(3)   ; do j = 0, mesh % elements(eID) % Nxyz(2) ; do i = 0, mesh % elements(eID) % Nxyz(1) ; do eq = 1, NCONS
               mesh % elements(eID) % storage % QDot(eq,i,j,k) = mesh % elements(eID) % storage % QDot(eq,i,j,k) + mesh % elements(eID) % storage % S_NS(eq,i,j,k)
            end do   ;  end do   ;  end do   ;  end do
         end do
!$acc end parallel loop
!$omp end do
!
!        *********************
!        Add IBM source term
!        *********************
! no wall function for MULTIPHASE
!$acc wait
         if( mesh% IBM% active ) then
            if( .not. mesh% IBM% semiImplicit ) then 
!$omp do schedule(runtime) private(i,j,k,Source,eID,eq)
                  ! Check if update(t) is required
                  !$acc parallel loop gang vector_length(128) present(mesh) copyin(t, locLevel) private(Source,eID) async(1)
                  do lID = 1, mesh % MLRK % MLIter(locLevel,1)
                     eID = mesh % MLRK % MLIter_eID(lID)
                     !$acc loop vector collapse(3) private(Source)
                     do k = 0, mesh % elements(eID) % Nxyz(3)   ; do j = 0, mesh % elements(eID) % Nxyz(2) ; do i = 0, mesh % elements(eID) % Nxyz(1)
                        if( mesh % elements(eID) % isInsideBody(i,j,k) ) then
                           !$acc loop seq
                           do eq = 1, NCONS
                               Source(eq) = 0.0_RP
                           end do
                           ! only without moving for now in MULTIPHASE
                           if( .not. mesh % IBM % stl(mesh % elements(eID) % STL(i,j,k)) % move ) then 
                              call IBM_SourceTerm(mesh % IBM, eID = eID, Q = mesh % elements(eID) % storage % Q(:,i,j,k), Source = Source, wallfunction = .false. )
                           end if 
                           !$acc loop seq
                           do eq = 1, NCONS
                                mesh % elements(eID) % storage % QDot(eq,i,j,k) = mesh % elements(eID) % storage % QDot(eq,i,j,k) + Source(eq)
                           end do
                        end if
                     end do                  ; end do                ; end do
                  end do
                  !$acc end parallel loop
!$omp end do       
            end if 
         end if

!$acc wait
!
      end subroutine ComputeNSTimeDerivative
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine TimeDerivative_VolumetricContribution(mesh, Level)
         use HexMeshClass
         use ElementClass
         use DGIntegrals
         implicit none
         type(HexMesh), intent(inout)             :: mesh
         integer      , intent(in)   , optional   :: Level
!
!        ---------------
!        Local variables
!        ---------------
!
         integer       :: i, j, k, eq, eID, lID, locLevel
         real(kind=RP) :: mu, kappa, beta
         real(kind=RP) :: inviscidFlux(1:NCONS, 1:NDIM)
         real(kind=RP) :: viscousFlux(1:NCONS, 1:NDIM)
         real(kind=RP) :: jGradXi(1:3), jGradEta(1:3), jGradZeta(1:3)
         
         if (present(Level)) then
            locLevel = Level
         else
            locLevel = 1
         end if
!
!        *************************************
!        Compute interior contravariant fluxes
!        *************************************
!

!        Compute inviscid - viscous contravariant flux
!        ---------------------------------------------
!$omp do schedule(runtime) private(eID, i, j, k, eq, beta, kappa, mu, inviscidFlux, viscousFlux, jGradXi, jGradEta, jGradZeta)
         !$acc parallel loop gang vector_length(128) present(mesh, multiphase) copyin(locLevel) private(eID) async(1)
         do lID = 1, mesh % MLRK % MLIter(locLevel,1)
            eID = mesh % MLRK % MLIter_eID(lID)

            !$acc loop vector collapse(3) private(beta, kappa, mu, inviscidFlux, viscousFlux, jGradXi, jGradEta, jGradZeta)
            do k = 0, mesh % elements(eID) % Nxyz(3) ; do j = 0, mesh % elements(eID) % Nxyz(2) ; do i = 0, mesh % elements(eID) % Nxyz(1)
                  
               call mEulerFlux(mesh % elements(eID) % storage % Q(:,i,j,k), inviscidFlux, mesh % elements(eID) % storage % rho(i,j,k))

               call GetmTwoFluidsViscosity(mesh % elements(eID) % storage % Q(IMC,i,j,k), mu)
               beta  = multiphase % M0_star
               kappa = 0.0_RP
               mu    = mu + mesh % elements(eID) % storage % mu_NS (1,i,j,k)  ! Add subgrid LES viscosity

               call mViscousFlux( NCONS, NGRAD, mesh % elements(eID) % storage % Q(:,i,j,k) , mesh % elements(eID) % storage % U_x(:,i,j,k) , & 
                                       mesh % elements(eID) % storage % U_y(:,i,j,k) , mesh % elements(eID) % storage % U_z(:,i,j,k), mu, beta, kappa, viscousFlux)
!
               jGradXi(1) = mesh % elements(eID) % geom % jGradXi(1,i,j,k)
               jGradXi(2) = mesh % elements(eID) % geom % jGradXi(2,i,j,k)
               jGradXi(3) = mesh % elements(eID) % geom % jGradXi(3,i,j,k)
               
               jGradEta(1) = mesh % elements(eID) % geom % jGradEta(1,i,j,k)
               jGradEta(2) = mesh % elements(eID) % geom % jGradEta(2,i,j,k)
               jGradEta(3) = mesh % elements(eID) % geom % jGradEta(3,i,j,k)
               
               jGradZeta(1) = mesh % elements(eID) % geom % jGradZeta(1,i,j,k)
               jGradZeta(2) = mesh % elements(eID) % geom % jGradZeta(2,i,j,k)
               jGradZeta(3) = mesh % elements(eID) % geom % jGradZeta(3,i,j,k)
               
               !$acc loop seq
               do eq =1, NCONS

               inviscidFlux(eq,1) = inviscidFlux(eq,1) - viscousFlux(eq,1)
               inviscidFlux(eq,2) = inviscidFlux(eq,2) - viscousFlux(eq,2)
               inviscidFlux(eq,3) = inviscidFlux(eq,3) - viscousFlux(eq,3)
                  
               mesh % elements(eID) % storage % contravariantFlux(eq,i,j,k,IX)  = &
                                                           inviscidFlux(eq,IX) * jGradXi(IX)  &
                                                         + inviscidFlux(eq,IY) * jGradXi(IY)  &
                                                         + inviscidFlux(eq,IZ) * jGradXi(IZ)

               mesh % elements(eID) % storage % contravariantFlux(eq,i,j,k,IY)  = &
                                                           inviscidFlux(eq,IX) * jGradEta(IX)  &
                                                         + inviscidFlux(eq,IY) * jGradEta(IY)  &
                                                         + inviscidFlux(eq,IZ) * jGradEta(IZ)
                  
               mesh % elements(eID) % storage % contravariantFlux(eq,i,j,k,IZ)  = &
                                                           inviscidFlux(eq,IX) * jGradZeta(IX)  &
                                                         + inviscidFlux(eq,IY) * jGradZeta(IY)  &
                                                         + inviscidFlux(eq,IZ) * jGradZeta(IZ)
               end do
            end do               ; end do                ; end do
!
!           Perform volume integrals
!           ------------------------
            call ScalarWeakIntegrals_StdVolumeGreen( mesh % elements(eID) % Nxyz, NCONS, mesh % elements(eID) % storage % contravariantFlux, &
                                                     mesh % elements(eID) % storage % QDot)

         end do
         !$acc end parallel loop 
!$omp end do
!
      end subroutine TimeDerivative_VolumetricContribution
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine TimeDerivative_FacesContribution(e, mesh)
         !$acc routine vector
         use HexMeshClass
         use DGIntegrals, only: ScalarWeakIntegrals_StdFace
         implicit none
         type(Element)           :: e
         type(HexMesh)           :: mesh

         integer                 :: i,j,k,eq

         call ScalarWeakIntegrals_StdFace( NCONS, e % Nxyz, &
                      mesh % faces(e % faceIDs(EFRONT))  % storage(e % faceSide(EFRONT))  % fStar, &
                      mesh % faces(e % faceIDs(EBACK))   % storage(e % faceSide(EBACK))   % fStar, &
                      mesh % faces(e % faceIDs(EBOTTOM)) % storage(e % faceSide(EBOTTOM)) % fStar, &
                      mesh % faces(e % faceIDs(ERIGHT))  % storage(e % faceSide(ERIGHT))  % fStar, &
                      mesh % faces(e % faceIDs(ETOP))    % storage(e % faceSide(ETOP))    % fStar, &
                      mesh % faces(e % faceIDs(ELEFT))   % storage(e % faceSide(ELEFT))   % fStar, &
                      -1, e % storage % QDot )

         !$acc loop vector collapse(3)
         do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
            !$acc loop seq
            do eq = 1,NCONS
				e % storage % QDot(eq,i,j,k) = e % storage % QDot(eq,i,j,k)  / e % geom % jacobian(i,j,k)
            enddo
         end do         ; end do          ; end do
         
      end subroutine TimeDerivative_FacesContribution
!
!///////////////////////////////////////////////////////////////////////////////////////////// 
! 
!        Riemann solver drivers 
!        ---------------------- 
! 
!///////////////////////////////////////////////////////////////////////////////////////////// 
! 
      subroutine computeElementInterfaceFlux_MUviscous(fc, side)
         !$acc routine vector
         use FaceClass
         use RiemannSolvers_MU
         implicit none
         type(Face)   , intent(inout) :: fc
         integer,       intent(in)    :: side

         integer       :: i, j
         real(kind=RP) :: mu
         real(kind=RP) :: U_xS(1:NCONS), U_yS(1:NCONS), U_zS(1:NCONS)
         
         !$acc loop vector collapse(2) private(U_xS, U_yS, U_zS, mu)
         do j = 0, fc % Nf(2) ; do i = 0, fc % Nf(1)

            call GetmTwoFluidsViscosity(fc % storage(side) % Q(IMC,i,j), mu)
			mu = mu + fc % storage(side) % mu_NS(1,i,j)   ! Add subgrid viscosity

            U_xS = fc % storage(side) % U_x(:,i,j)*[1.0_RP,mu,mu,mu,1.0_RP]
            U_yS = fc % storage(side) % U_y(:,i,j)*[1.0_RP,mu,mu,mu,1.0_RP]
            U_zS = fc % storage(side) % U_z(:,i,j)*[1.0_RP,mu,mu,mu,1.0_RP] 

            call mViscousFlux( NCONS, NGRAD, fc % storage(side) % Q(:,i,j) , &
                              U_xS, U_yS, U_zS, 1.0_RP, multiphase % M0_star, &
                              0.0_RP, fc % storage(side) % unStar(:,:,i,j))   
                              
         end do ; end do
      end subroutine computeElementInterfaceFlux_MUviscous

      SUBROUTINE computeElementInterfaceFlux_MU(mesh, Level)
         use HexMeshClass
         use RiemannSolvers_MU
         use EllipticBR1
         use FaceClass
         implicit none
         type(HexMesh), intent(inout)  :: mesh
         integer, intent(in), optional :: Level
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                       :: i, j, eq, iFace, fID, locLevel
         
         if (present(Level)) then
            locLevel = Level
         else
            locLevel = 1
         end if

!$omp do schedule(runtime) private(fID, i, j, eq)
         !$acc parallel loop gang vector_length(32) present(mesh, ViscousDiscretization, multiphase) copyin(locLevel) private(fID)
         do iFace = 1, mesh % MLRK % MLIter(locLevel,3)
            fID = mesh % MLRK % MLIter_fID_Interior(iFace)
         
            call BR1_RiemannSolver_acc(mesh % faces(fID), NCONS, NGRAD, [multiphase % M0_star, 0.0_RP, 0.0_RP, 0.0_RP, 0.0_RP], &
                                         ViscousDiscretization % sigma, mesh % faces(fID) % storage(2) % FStar)

            call RiemannSolver_Selector_MU(mesh % faces(fID) % Nf(1), &                         
                                           mesh % faces(fID) % Nf(2), &
                                           mesh % faces(fID) % storage(1) % Q, &
                                           mesh % faces(fID) % storage(2) % Q, &
                                           mesh % faces(fID) % storage(1) % rho, &
                                           mesh % faces(fID) % storage(2) % rho, &
                                           mesh % faces(fID) % storage(1) % mu,&
                                           mesh % faces(fID) % storage(2) % mu,&
                                           mesh % faces(fID) % geom % normal, &
                                           mesh % faces(fID) % geom % t1, &
                                           mesh % faces(fID) % geom % t2, &
                                           mesh % faces(fID) % storage(1) % Q_aux,&
                                           mesh % faces(fID) % storage(2) % Q_aux,&
                                           mesh % faces(fID) % storage(1) % invMa2, &
                                           mesh % faces(fID) % storage(2) % invMa2)

!           ------------------------
!           Multiply by the Jacobian
!           ------------------------
            !$acc loop vector collapse(3)
            do j = 0, mesh % faces(fID) % Nf(2) ; do i = 0, mesh % faces(fID) % Nf(1) ; do eq = 1, NCONS
               mesh % faces(fID) % storage(1) % FStar(eq,i,j) = ( mesh % faces(fID) % storage(1) % Q_aux(eq,i,j) - mesh % faces(fID) % storage(2) % FStar(eq,i,j)) * mesh % faces(fID) % geom % jacobian(i,j)
               mesh % faces(fID) % storage(2) % Q_aux(eq,i,j) = ( mesh % faces(fID) % storage(2) % Q_aux(eq,i,j) - mesh % faces(fID) % storage(2) % FStar(eq,i,j)) * mesh % faces(fID) % geom % jacobian(i,j)
            end do ; end do ;  end do
!
!           ---------------------------
!           Return the flux to elements
!           ---------------------------
            call Face_ProjectFluxToElements(mesh % faces(fID), NCONS, mesh % faces(fID) % storage(1) % FStar, 1)   ! For storage(1) FStar is fluxL
            call Face_ProjectFluxToElements(mesh % faces(fID), NCONS, mesh % faces(fID) % storage(2) % Q_aux, 2)   ! For storage(2) Q_aux is fluxR

         end do
         !$acc end parallel loop
!$omp end do nowait

      END SUBROUTINE computeElementInterfaceFlux_MU

      SUBROUTINE computeMPIFaceFlux_MU(f)
         !$acc routine vector
         use FaceClass
         use RiemannSolvers_MU
         use EllipticBR1
         TYPE(Face)   , INTENT(inout) :: f   
!
!        ---------------
!        Local variables
!        ---------------
!
         integer       :: i, j, Sidearray, maxId, eq, k
         real(kind=RP) :: muL, muR, flux(2)
         real(kind=RP) :: UxL(1:NGRAD), UyL(1:NGRAD), UzL(1:NGRAD)
         real(kind=RP) :: UxR(1:NGRAD), UyR(1:NGRAD), UzR(1:NGRAD)
         real(kind=RP) :: scaleMu(1:NGRAD)

         associate(s1 => f%storage(1), s2 => f%storage(2))
       
         !$acc loop vector collapse(2) private(flux, UxL, UyL, UzL, UxR, UyR, UzR, muL, muR, scaleMu, k)
         DO j = 0, f % Nf(2) ; DO i = 0, f % Nf(1)

              ! compute viscosity on each side and add subgrid term
              call GetmTwoFluidsViscosity(s1%Q(IMC,i,j), muL)
              call GetmTwoFluidsViscosity(s2%Q(IMC,i,j), muR)
              
              muL = muL +f % storage(1) % mu_NS(1,i,j)   ! Add subgrid viscosity
              muR = muR +f % storage(2) % mu_NS(1,i,j)   ! Add subgrid viscosity     

              ! Left side
              !-----------------------------------
              ! prepare scaling factors: 1 on first/last components, mu on interior gradient components                         
              scaleMu = [1.0_RP, muL, muL, muL, 1.0_RP]

              ! fill left-side premultiplied gradient arrays
			  !$acc loop seq
              do k = 1, NGRAD
                 UxL(k) = s1%U_x(k,i,j) * scaleMu(k)
                 UyL(k) = s1%U_y(k,i,j) * scaleMu(k)
                 UzL(k) = s1%U_z(k,i,j) * scaleMu(k)
              end do

              ! Right side
              !-----------------------------------
              scaleMu = [1.0_RP, muR, muR, muR, 1.0_RP]

              ! fill right-side premultiplied gradient arrays
			  !$acc loop seq
              do k = 1, NGRAD                        
                 UxR(k) = s2%U_x(k,i,j) * scaleMu(k)
                 UyR(k) = s2%U_y(k,i,j) * scaleMu(k)
                 UzR(k) = s2%U_z(k,i,j) * scaleMu(k)
              end do

              ! compute viscous flux (left and right)
              call mViscousFlux( NCONS, NGRAD, s1%Q(:,i,j) , &
                                 UxL, UyL, UzL, 1.0_RP, multiphase % M0_star, &
                                                         
                                 0.0_RP, s1%unStar(:,:,i,j))

              call mViscousFlux( NCONS, NGRAD, s2%Q(:,i,j) , &
                                 UxR, UyR, UzR, 1.0_RP, multiphase % M0_star, &
                                 0.0_RP, s2%unStar(:,:,i,j))

       END DO ; END DO

       ! Riemann solver 
       call BR1_RiemannSolver_acc(f, NCONS, NGRAD, [multiphase % M0_star, 0.0_RP, 0.0_RP, 0.0_RP, 0.0_RP], &
                                  ViscousDiscretization % sigma, s2%Fstar)

       ! inviscid fluxes (unchanged). Ensure f%storage(*) arrays used here are present on device.
       call RiemannSolver_Selector_MU(f % Nf(1), &                         
                                     f % Nf(2), &
                                     s1%Q, &
                                     s2%Q, &
                                     s1%rho, &
                                     s2%rho, &
                                     s1%mu, &
                                     s2%mu, &
                                     f % geom % normal, &
                                     f % geom % t1, &
                                     f % geom % t2, &
                                     s1%Q_aux, &
                                     s2%Q_aux, &
                                     s1%invMa2, &
                                     s2%invMa2)
                                     
!      ------------------------
!      Multiply by the Jacobian -- Q_aux is the inviscid flux and Viscous Flux is f % storage(2) % FStar
!      ------------------------

       !$acc loop vector collapse(3)
       do j = 0, f % Nf(2) ; do i = 0, f % Nf(1) ; do eq = 1, NCONS
          s1%Q_aux(eq,i,j) = ( s1%Q_aux(eq,i,j) - s2%Fstar(eq,i,j)) * f % geom % jacobian(i,j)
          s2%Q_aux(eq,i,j) = ( s2%Q_aux(eq,i,j) - s2%Fstar(eq,i,j)) * f % geom % jacobian(i,j)
       end do ; end do ; end do

       ! Find the largest element ID and the side array index 
       maxId = -1
       !$acc loop seq
       do i = 1, size(f % elementIDs)
          if (f%elementIDs(i) > maxId) maxId = f%elementIDs(i)
       end do

       !$acc loop seq
       do i = 1, size(f % elementIDs)
          if (f % elementIDs(i) == maxId) then
             Sidearray = i
             exit
          end if
       end do

       call Face_ProjectFluxToElements(f, NCONS, f % storage(Sidearray) % Q_aux, Sidearray)

       end associate
    END SUBROUTINE computeMPIFaceFlux_MU

      SUBROUTINE computeBoundaryFlux_MU(mesh, time, Level)
      USE ElementClass
      use FaceClass
      USE RiemannSolvers_MU
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      type(HexMesh), intent(inout)    :: mesh
      REAL(KIND=RP)                   :: time
      integer, intent(in), optional   :: Level
!
!     ---------------
!     Local variables
!     ---------------
!
      INTEGER                         :: i, j, eq, locLevel
      INTEGER                         :: nZones, zoneID, zonefID, fID
      real(kind=RP)                   :: mu
      
      if (present(Level)) then
         locLevel = Level
      else
         locLevel = 1
      end if

      nZones = size(mesh % zones)
       do zoneID=1, nZones
         !$acc wait
         call BCs(zoneID) % bc % FlowState(mesh, mesh % zones(zoneID))  
         !$acc wait
!$omp do schedule(runtime) private(fID, i, j, eq, mu )
         !$acc parallel loop gang vector_length(32) present(mesh, dimensionless, multiphase) copyin(locLevel,zoneID) private(fID) async(1) 
         do zonefID = 1, mesh % zones(zoneID) % no_of_faces
             fID =  mesh % zones(zoneID) % faces(zonefID)
    
            !$acc loop vector collapse(2) private(i,j)
            do j = 0, mesh % faces(fID) % Nf(2) ;  do i = 0, mesh % faces(fID) % Nf(1)
               mesh % faces(fID) % storage(2) % rho(i,j) = dimensionless % rho(2) + (dimensionless % rho(1)-dimensionless % rho(2))*mesh % faces(fID) % storage(2) % Q(IMC,i,j)
               mesh % faces(fID) % storage(2) % rho(i,j) = min(max(mesh % faces(fID) % storage(2) % rho(i,j), dimensionless % rho_min),dimensionless % rho_max)
            enddo ; enddo

            !$acc loop vector collapse(2) private(i,j,eq,mu)
            do j = 0, mesh % faces(fID) % Nf(2) ;  do i = 0, mesh % faces(fID) % Nf(1)
               call GetmTwoFluidsViscosity(mesh % faces(fID) % storage(1) % Q(IMC,i,j), mu)
               mu = mu + mesh % faces(fID) % storage(1) % mu_NS(1,i,j)   ! Add subgrid viscosity
               call mViscousFlux(NCONS, NCONS, mesh % faces(fID) % storage(1) % Q(:,i,j), &
                                 mesh % faces(fID) % storage(1) % U_x(:,i,j), &
                                 mesh % faces(fID) % storage(1) % U_y(:,i,j), &
                                 mesh % faces(fID) % storage(1) % U_z(:,i,j), &
                                 mu, multiphase % M0_star, 0.0_RP, &
                                 mesh % faces(fID) % storage(1) % unStar(:,:,i,j))
                                 
                !$acc loop seq
                do eq = 1, NCONS
                    mesh % faces(fID) % storage(2) % FStar(eq,i,j) = mesh % faces(fID) % storage(1) % unStar(eq,IX,i,j)* mesh % faces(fID) % geom % normal(IX,i,j) &
                                                                 + mesh % faces(fID) % storage(1) % unStar(eq,IY,i,j)* mesh % faces(fID) % geom % normal(IY,i,j) &
                                                                 + mesh % faces(fID) % storage(1) % unStar(eq,IZ,i,j)* mesh % faces(fID) % geom % normal(IZ,i,j)
                end do 

            enddo ; enddo
            
         end do
         !$acc end parallel loop 
!$omp end do 

         !$acc wait
         CALL BCs(zoneID) % bc % FlowNeumann(mesh, mesh % zones(zoneID))    
         !$acc wait
         
!$omp do schedule(runtime) private(fID, i, j, eq )
         !$acc parallel loop gang vector_length(32) present(mesh) copyin(locLevel,zoneID) private(fID, i, j, eq) async(1)
         do zonefID = 1, mesh % zones(zoneID) % no_of_faces
            fID =  mesh % zones(zoneID) % faces(zonefID)

            call RiemannSolver_Selector_MU(mesh % faces(fID) % Nf(1), &                         
                                          mesh % faces(fID) % Nf(2), &
                                          mesh % faces(fID) % storage(1) % Q, &
                                          mesh % faces(fID) % storage(2) % Q, &
                                          mesh % faces(fID) % storage(1) % rho, &
                                          mesh % faces(fID) % storage(2) % rho, &
                                          mesh % faces(fID) % storage(1) % mu,&
                                          mesh % faces(fID) % storage(2) % mu,&
                                          mesh % faces(fID) % geom % normal, &
                                          mesh % faces(fID) % geom % t1, &
                                          mesh % faces(fID) % geom % t2, &
                                          mesh % faces(fID) % storage(1) % Q_aux,&
                                          mesh % faces(fID) % storage(2) % Q_aux,&
                                          mesh % faces(fID) % storage(1) % invMa2, &
                                          mesh % faces(fID) % storage(2) % invMa2)
!           ------------------------
!           Multiply by the Jacobian
!           ------------------------
            !$acc loop vector collapse(3)
            do j = 0, mesh % faces(fID) % Nf(2) ; do i = 0, mesh % faces(fID) % Nf(1) ; do eq = 1, NCONS               
                                                                    
                  mesh % faces(fID) % storage(1) % FStar(eq,i,j) = (mesh % faces(fID) % storage(1) % Q_aux(eq,i,j)  - &
                                                                    mesh % faces(fID) % storage(2) % FStar(eq,i,j)) * &
                                                                    mesh % faces(fID) % geom % jacobian(i,j)
            end do ; end do ; enddo
            !
            !           ---------------------------
            !           Return the flux to elements
            !           ---------------------------
            !
            call Face_ProjectFluxToElements(mesh % faces(fID), NCONS, mesh % faces(fID) % storage(1) % FStar, 1)
         enddo
         !$acc end parallel loop 
!$omp end do 
!$acc wait
      end do 

      END SUBROUTINE computeBoundaryFlux_MU
!
!////////////////////////////////////////////////////////////////////////////////////////
!
!        Laplacian procedures
!        --------------------
!
!////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine ComputeLaplacian( mesh , t, Level)
         implicit none
         type(HexMesh)                   :: mesh
         real(kind=RP)                   :: t
         integer, intent(in), optional   :: Level
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: eID , i, j, k, ierr, fID, locLevel
         integer     :: iFace, iEl, lID
         
         
         if (present(Level)) then
            locLevel = Level
         else
            locLevel = 1 
         end if
!
!        ****************
!        Volume integrals
!        ****************
!
!$acc wait
         call Laplacian_VolumetricContribution(mesh, Level=locLevel)
!
!        ******************************************
!        Compute Riemann solver of non-shared faces
!        ******************************************
!
!$acc wait
         call Laplacian_computeElementInterfaceFlux(mesh, Level=locLevel)
!$acc wait
         call Laplacian_computeBoundaryFlux(mesh, t)
!
!        ***************************************************************
!        Surface integrals and scaling of elements with non-shared faces
!        ***************************************************************
! 
!$acc wait
!$omp do schedule(runtime) private(i,j,k,eID)
!$acc parallel loop gang num_gangs(mesh % MLRK % MLIter(locLevel,5)) vector_length(128) present(mesh, mesh % MLRK) copyin(locLevel) private(eID) async(1)
         do iEl = 1, mesh % MLRK % MLIter(locLevel,5)
            eID = mesh % MLRK % MLIter_eID_Seq(iEl)
            call Laplacian_FacesContribution(mesh, eID)
         end do
!$acc end parallel loop 
!$omp end do
!
!        ***********************************************************
!        Surface integrals and scaling of elements with shared faces
!        ***********************************************************
! 
!$acc wait
#ifdef _HAS_MPI_
         if ( MPI_Process % doMPIAction ) then
!$omp single
            call HexMesh_GatherMPIFacesGradients(mesh, NCOMP)
!$omp end single
!$acc wait
!
!           **************************************
!           Compute Riemann solver of shared faces
!           **************************************
!
!$omp do schedule(runtime) private(fID)
!$acc parallel loop gang vector_length(32) present(mesh, CHDiscretization, mesh % MLRK, mesh % faces) copyin(locLevel) private(fID) async(1)       
            do iFace = 1, mesh % MLRK % MLIter(locLevel,7)
               fID = mesh % MLRK % MLIter_fID_MPI(iFace)
               CALL Laplacian_computeMPIFaceFlux( mesh % faces(fID) )
            end do
!$acc end parallel loop             
!$omp end do 
!$acc wait
!
!           ***********************************************************
!           Surface integrals and scaling of elements with shared faces
!           ***********************************************************
!$omp do schedule(runtime) private(eID)  
!$acc parallel loop gang vector_length(128) present(mesh, mesh % MLRK) copyin(locLevel) private(eID) async(1) 
            do iEl = 1, mesh % MLRK % MLIter(locLevel,6)
               eID = mesh % MLRK % MLIter_eID_MPI(iEl)
               call Laplacian_FacesContribution(mesh, eID)  
            end do
!$acc end parallel loop     
!$omp end do
!$acc wait
!
!           Add a MPI Barrier
!           -----------------
!$omp single
            call mpi_barrier(MPI_COMM_WORLD, ierr)
!$omp end single
         end if
#endif

      end subroutine ComputeLaplacian
!
!////////////////////////////////////////////////////////////////////////////////////////
!
!     Ger: This is not called somewhere
!      subroutine ComputeLaplacianNeumannBCs( mesh , t)
!         implicit none
!         type(HexMesh)              :: mesh
!         real(kind=RP)              :: t
!
!        ---------------
!        Local variables
!        ---------------
!
!         integer     :: eID , i, j, k, ierr, fID
!
!        **************************
!        Reset QDot and face fluxes
!        **************************
!
!         do eID = 1, mesh % no_of_elements
!            mesh % elements(eID) % storage % QDot = 0.0_RP
!         end do
!   
!         do fID = 1, size(mesh % faces)
!            mesh % faces(fID) % storage(1) % genericInterfaceFluxMemory = 0.0_RP
!            mesh % faces(fID) % storage(2) % genericInterfaceFluxMemory = 0.0_RP
!         end do
!
!        ******************************************
!        Compute Riemann solver of non-shared faces
!        ******************************************
!
!!$omp do schedule(runtime) 
!         do fID = 1, size(mesh % faces) 
!            associate( f => mesh % faces(fID)) 
!            select case (f % faceType) 
!            case (HMESH_BOUNDARY) 
!               CALL Laplacian_computeBoundaryFlux(f, t) 
!            end select 
!            end associate 
!         end do 
!!$omp end do 
!
!        ***************************************************************
!        Surface integrals and scaling of elements with non-shared faces
!        ***************************************************************
! 
!!$omp do schedule(runtime) private(i, j, k)
!         do eID = 1, size(mesh % elements) 
!            associate(e => mesh % elements(eID)) 
!            if ( e % hasSharedFaces ) cycle
!            call Laplacian_FacesContribution(e, t, mesh) 
! 
!            do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1) 
!               e % storage % QDot(:,i,j,k) = e % storage % QDot(:,i,j,k) / e % geom % jacobian(i,j,k) 
!            end do         ; end do          ; end do 
!            end associate 
!         end do
!!$omp end do
!
!      end subroutine ComputeLaplacianNeumannBCs
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine Laplacian_VolumetricContribution(mesh, Level)
         use HexMeshClass
         use ElementClass
         use DGIntegrals
         implicit none
         type(HexMesh), intent (inout)             :: mesh
         integer,       intent(in)     , optional  :: Level
!
!        ---------------
!        Local variables
!        ---------------
!
         integer       :: eID, i, j, k, lID, locLevel
         real(kind=RP) :: mu, beta, kappa
         real(kind=RP) :: cartesianFlux(1:NDIM)

         if (present(Level)) then
            locLevel = Level
         else
            locLevel = 1
         end if
!
!        *************************************
!        Compute interior contravariant fluxes
!        *************************************
!
!        Compute contravariant flux
!        --------------------------
!$omp do schedule(runtime) private(eID, kappa, beta,mu,i,j,k, cartesianFlux)
         !$acc parallel loop gang vector_length(128) present(mesh, mesh % elements, mesh % MLRK) copyin(locLevel) private(eID) async(1)
         do lID = 1, mesh % MLRK % MLIter(locLevel,1)
            eID = mesh % MLRK % MLIter_eID(lID)

            !$acc loop vector collapse(3) private(cartesianFlux, kappa, beta,mu)
            do k = 0, mesh % elements(eID) % Nxyz(3) ; do j = 0, mesh % elements(eID) % Nxyz(2) ; do i = 0, mesh % elements(eID) % Nxyz(1)
               
               call GetCHViscosity( mesh % elements(eID) % storage % Q(IMC,i,j,k), mu)      
               kappa = 0.0_RP
               beta  = multiphase % M0_star

               call CHDivergenceFlux( NCOMP, NCOMP, mesh % elements(eID) % storage % Q(1:IMC,i,j,k) , mesh % elements(eID) % storage % U_x(1:IMC,i,j,k) , & 
                                      mesh % elements(eID) % storage % U_y(1:IMC,i,j,k) , mesh % elements(eID) % storage % U_z(1:IMC,i,j,k), mu, beta, kappa, cartesianFlux)
            
               mesh % elements(eID) % storage % contravariantFlux(IMC,i,j,k,IX)  = &
                                                         - cartesianFlux(IX) * mesh % elements(eID) % geom % jGradXi(IX,i,j,k)  &
                                                         - cartesianFlux(IY) * mesh % elements(eID) % geom % jGradXi(IY,i,j,k)  &
                                                         - cartesianFlux(IZ) * mesh % elements(eID) % geom % jGradXi(IZ,i,j,k)

               mesh % elements(eID) % storage % contravariantFlux(IMC,i,j,k,IY)  = &
                                                         - cartesianFlux(IX) * mesh % elements(eID) % geom % jGradEta(IX,i,j,k)  &
                                                         - cartesianFlux(IY) * mesh % elements(eID) % geom % jGradEta(IY,i,j,k)  &
                                                         - cartesianFlux(IZ) * mesh % elements(eID) % geom % jGradEta(IZ,i,j,k)
                  
               mesh % elements(eID) % storage % contravariantFlux(IMC,i,j,k,IZ)  = &
                                                         - cartesianFlux(IX) * mesh % elements(eID) % geom % jGradZeta(IX,i,j,k)  &
                                                         - cartesianFlux(IY) * mesh % elements(eID) % geom % jGradZeta(IY,i,j,k)  &
                                                         - cartesianFlux(IZ) * mesh % elements(eID) % geom % jGradZeta(IZ,i,j,k)
               
               mesh % elements(eID) % storage % Qdot(NCOMP,i,j,k)  = 0.0_RP

            end do               ; end do            ; end do
!
!           ************************
!           Perform volume integrals
!           ************************
            call ScalarWeakIntegrals_StdVolumeGreen( mesh % elements(eID) % Nxyz, NCOMP, mesh % elements(eID) % storage % contravariantFlux, &
                                                     mesh % elements(eID) % storage % QDot)
         end do
         !$acc end parallel loop 
!$omp end do
!$acc wait

      end subroutine Laplacian_VolumetricContribution
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine Laplacian_FacesContribution(mesh, eID)
         !$acc routine vector
         use HexMeshClass
         use PhysicsStorage
         use DGIntegrals, only: ScalarWeakIntegrals_StdFace
         implicit none
         type(HexMesh)           :: mesh
         integer, intent(in)     :: eID

         integer                 :: i,j,k

         call ScalarWeakIntegrals_StdFace( NCOMP, mesh % elements(eID) % Nxyz, &
                      mesh % faces(mesh % elements(eID) % faceIDs(EFRONT))  % storage(mesh % elements(eID) % faceSide(EFRONT))  % fStar, &
                      mesh % faces(mesh % elements(eID) % faceIDs(EBACK))   % storage(mesh % elements(eID) % faceSide(EBACK))   % fStar, &
                      mesh % faces(mesh % elements(eID) % faceIDs(EBOTTOM)) % storage(mesh % elements(eID) % faceSide(EBOTTOM)) % fStar, &
                      mesh % faces(mesh % elements(eID) % faceIDs(ERIGHT))  % storage(mesh % elements(eID) % faceSide(ERIGHT))  % fStar, &
                      mesh % faces(mesh % elements(eID) % faceIDs(ETOP))    % storage(mesh % elements(eID) % faceSide(ETOP))    % fStar, &
                      mesh % faces(mesh % elements(eID) % faceIDs(ELEFT))   % storage(mesh % elements(eID) % faceSide(ELEFT))   % fStar, &
                      1, mesh % elements(eID) % storage % QDot )

         !$acc loop vector collapse(3)
         do k = 0, mesh % elements(eID) % Nxyz(3) ; do j = 0, mesh % elements(eID) % Nxyz(2) ; do i = 0, mesh % elements(eID) % Nxyz(1)
            mesh % elements(eID) % storage % QDot(IMC,i,j,k) = mesh % elements(eID) % storage % QDot(IMC,i,j,k)  / mesh % elements(eID) % geom % jacobian(i,j,k)
         end do         ; end do          ; end do
         
      end subroutine Laplacian_FacesContribution
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine compute_viscosity_at_faces(no_of_faces, no_of_sides, face_ids, mesh)
         implicit none
         integer, intent(in)           :: no_of_faces
         integer, intent(in)           :: no_of_sides
         integer, intent(in)           :: face_ids(no_of_faces)
         class(HexMesh), intent(inout) :: mesh
!
!        ---------------
!        Local variables
!        ---------------
!
         integer       :: iFace, i, j, side
         real(kind=RP) :: delta, mu_smag, factor, minvalC, maxvalC, minvalC2, maxvalC2
         real(kind=RP) :: dWall_ij, prod

         if ( LESModel % Active ) then
!$omp do schedule(runtime) private(i,j,delta,mu_smag,side, dWall_ij, prod)
            !$acc parallel loop gang vector_length(32) present(mesh, LESModel, mesh % faces) private(delta,prod) copyin(no_of_faces,no_of_sides,face_ids) async(1)
            do iFace = 1, no_of_faces
               prod  = (mesh % faces(face_ids(iFace)) % Nf(1) + 1.0_RP)*(mesh % faces(face_ids(iFace)) % Nf(2) + 1.0_RP)
			   delta = sqrt(mesh % faces(face_ids(iFace)) % geom % surface / prod)
              
               !$acc loop vector collapse(2) private(mu_smag, dWall_ij)
               do j = 0, mesh % faces(face_ids(iFace)) % Nf(2) ; do i = 0, mesh % faces(face_ids(iFace)) % Nf(1)
                  dWall_ij = mesh % faces(face_ids(iFace)) % geom % dWall(i,j)
				  
				  !$acc loop seq
                  do side = 1, no_of_sides
                     call LESModel_Selector(LESModel, delta, dWall_ij, &
                                                             mesh % faces(face_ids(iFace)) % storage(side) % Q(:,i,j),   &
                                                             mesh % faces(face_ids(iFace)) % storage(side) % U_x(:,i,j), &
                                                             mesh % faces(face_ids(iFace)) % storage(side) % U_y(:,i,j), &
                                                             mesh % faces(face_ids(iFace)) % storage(side) % U_z(:,i,j), &
                                                                                    mu_smag)
                     mesh % faces(face_ids(iFace)) % storage(side) % mu_NS(1,i,j) = mu_smag 
                  end do
               end do              ; end do
            end do
            !$acc end parallel loop
!$omp end do
!$acc wait
         end if

      end subroutine compute_viscosity_at_faces
!
!///////////////////////////////////////////////////////////////////////////////////////////// 
! 
!        Riemann solver drivers 
!        ---------------------- 
! 
!///////////////////////////////////////////////////////////////////////////////////////////// 
! 
      subroutine Laplacian_computeElementInterfaceFlux(mesh, Level)
         use FaceClass
         use Physics
         use PhysicsStorage
         use EllipticBR1
         IMPLICIT NONE
         type(HexMesh), intent (inout)           :: mesh
         integer, intent(in), optional :: Level
!
!        ---------------
!        Local variables
!        ---------------
!
         integer        :: i, j, iFace, fID, locLevel
         real(kind=RP)  :: mu
         
         if (present(Level)) then
            locLevel = Level
         else
            locLevel = 1
         end if
         
!$omp do schedule(runtime) private(fID, i, j, mu)
!$acc parallel loop gang vector_length(128) present(mesh, CHDiscretization, mesh % faces) copyin(locLevel) private(fID) async(1)
         do iFace = 1, mesh % MLRK % MLIter(locLevel,3)
            fID = mesh % MLRK % MLIter_fID_Interior(iFace)
            
            !$acc loop vector collapse(2) private(mu)
            do j = 0, mesh % faces(fID) % Nf(2)  ;  do i = 0, mesh % faces(fID) % Nf(1)

               call GetCHViscosity(0.0_RP, mu)

               call CHDivergenceFlux( NCONS, NCONS, mesh % faces(fID) % storage(1) % Q(1:IMC,i,j) , mesh % faces(fID) % storage(1) % U_x(1:IMC,i,j) , & 
                                      mesh % faces(fID) % storage(1) % U_y(1:IMC,i,j) , mesh % faces(fID) % storage(1) % U_z(1:IMC,i,j), mu, 0.0_RP, 0.0_RP, mesh % faces(fID) % storage(1) % unStar(:,:,i,j))

               call CHDivergenceFlux( NCONS, NCONS, mesh % faces(fID) % storage(2) % Q(1:IMC,i,j) , mesh % faces(fID) % storage(2) % U_x(1:IMC,i,j) , & 
                                      mesh % faces(fID) % storage(2) % U_y(1:IMC,i,j) , mesh % faces(fID) % storage(2) % U_z(1:IMC,i,j), mu, 0.0_RP, 0.0_RP, mesh % faces(fID) % storage(2) % unStar(:,:,i,j))
            
            end do   ;  end do 

            call BR1_RiemannSolver_acc(mesh % faces(fID), NCOMP, NCOMP, [1.0_RP], CHDiscretization % sigma, mesh % faces(fID) % storage(1) % FStar)

!           ------------------------
!           Multiply by the Jacobian
!           ------------------------
            !$acc loop vector collapse(2)
            do j = 0, mesh % faces(fID) % Nf(2) ; do i = 0, mesh % faces(fID) % Nf(1) 
               mesh % faces(fID) % storage(1) % FStar(1,i,j) = mesh % faces(fID) % storage(1) % FStar(1,i,j) * mesh % faces(fID) % geom % jacobian(i,j)
            end do ; end do
!
!           ---------------------------
!           Return the flux to elements
!           ---------------------------
!
            call Face_ProjectFluxToElements(mesh % faces(fID), NCOMP, mesh % faces(fID) % storage(1) % FStar, 1)
            call Face_ProjectFluxToElements(mesh % faces(fID), NCOMP, mesh % faces(fID) % storage(1) % FStar, 2)

         end do
!$acc end parallel loop
!$omp end do

      end subroutine Laplacian_computeElementInterfaceFlux

      subroutine Laplacian_computeMPIFaceFlux(f)
         !$acc routine vector
         use FaceClass
         use Physics
         use PhysicsStorage
         use EllipticBR1
         IMPLICIT NONE
         TYPE(Face)   , INTENT(inout) :: f   
         integer       :: i, j, m
         integer       :: maxId, Sidearray
         real(kind=RP) :: mu
         real(kind=RP) :: sigma0
         
         !$acc loop vector collapse(2) private(mu)
         do j = 0, f % Nf(2)  ;  do i = 0, f % Nf(1)
           call GetCHViscosity(0.0_RP, mu)
                                                                                  
           call CHDivergenceFlux( NCONS, NCONS, f % storage(1) % Q(1:IMC,i,j) , f % storage(1) % U_x(1:IMC,i,j) , & 
                                  f % storage(1) % U_y(1:IMC,i,j) , f % storage(1) % U_z(1:IMC,i,j), mu, 0.0_RP, 0.0_RP, f % storage(1) % unStar(:,:,i,j))

           call CHDivergenceFlux( NCONS, NCONS, f % storage(2) % Q(1:IMC,i,j) , f % storage(2) % U_x(1:IMC,i,j) , & 
                                  f % storage(2) % U_y(1:IMC,i,j) , f % storage(2) % U_z(1:IMC,i,j), mu, 0.0_RP, 0.0_RP, f % storage(2) % unStar(:,:,i,j))
         end do   ;  end do 
         
#ifdef _OPENACC
         call BR1_RiemannSolver_acc(f, NCOMP, NCOMP, 1.0_RP, CHDiscretization % sigma, f % storage(1) % FStar)
#else
         call BR1_RiemannSolver_acc(f, NCOMP, NCOMP, [1.0_RP], CHDiscretization % sigma, f % storage(1) % FStar)
#endif

!        ------------------------
!        Multiply by the Jacobian
!        ------------------------
         !$acc loop vector collapse(2)
         do j = 0, f % Nf(2) ; do i = 0, f % Nf(1) 
           f % storage(1) % FStar(1,i,j) = f % storage(1) % FStar(1,i,j) * f % geom % jacobian(i,j)
         end do ; end do
!
!        ---------------------------
!        Return the flux to elements
!        ---------------------------
!
         maxId=-1
         !$acc loop seq
         do i = 1, size(f % elementIDs)
            if (f%elementIDs(i) > maxId) then
                maxId = f%elementIDs(i)
            end if
         end do

         !$acc loop seq
         do i=1, size(f % elementIDs)
             if(f % elementIDs(i)==maxId)THEN
               Sidearray=i
               exit
             endif
         end do
         
         call Face_ProjectFluxToElements(f, NCOMP, f % storage(1) % FStar, Sidearray)

      end subroutine Laplacian_ComputeMPIFaceFlux

      subroutine Laplacian_computeBoundaryFlux(mesh, time)
         USE ElementClass
         use FaceClass
         USE EllipticDiscretizations
         USE Physics
         use PhysicsStorage
         IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
         type(HexMesh), intent(inout)  :: mesh
         REAL(KIND=RP)                 :: time
!
!     ---------------
!     Local variables
!     ---------------
!
         INTEGER                         :: i, j
         INTEGER, DIMENSION(2)           :: N
         real(kind=RP)                   :: mu
         integer                         :: nZones, zoneID, zonefID, fID
!
!     -------------------
!     Get external states
!     -------------------
!

         nZones = size(mesh % zones)
         do zoneID=1, nZones
!$omp do schedule(runtime) private(fID, i, j, mu)         
            !$acc parallel loop gang vector_length(128) present(mesh, mesh % zones) copyin(zoneID) private(fID) async(1)
            do zonefID = 1, mesh % zones(zoneID) % no_of_faces
               fID =  mesh % zones(zoneID) % faces(zonefID)
               
               !$acc loop vector collapse(2) private(mu)
               do j = 0, mesh % faces(fID) % Nf(2) ; do i = 0, mesh % faces(fID) % Nf(1)
               mesh % faces(fID) % storage(2) % Q(IMC,i,j) = mesh % faces(fID) % storage(1) % Q(IMC,i,j)
               call GetCHViscosity(0.0_RP, mu)
               call CHDivergenceFlux(NCONS, NCONS, mesh % faces(fID) % storage(1) %   Q(1:IMC,i,j), &
                                                   mesh % faces(fID) % storage(1) % U_x(1:IMC,i,j), &
                                                   mesh % faces(fID) % storage(1) % U_y(1:IMC,i,j), &
                                                   mesh % faces(fID) % storage(1) % U_z(1:IMC,i,j), &
                                                   mu, 0.0_RP, 0.0_RP, &
                                                   mesh % faces(fID) % storage(1) % unStar(:,:,i,j))
               enddo ; enddo


               !$acc loop vector collapse(2)
               do j = 0, mesh % faces(fID) % Nf(2) ; do i = 0, mesh % faces(fID) % Nf(1)
                  mesh % faces(fID) % storage(2) % FStar(NCOMP,i,j) = mesh % faces(fID) % storage(1) % unStar(NCOMP,IX,i,j)* mesh % faces(fID) % geom % normal(IX,i,j) &
                                                                    + mesh % faces(fID) % storage(1) % unStar(NCOMP,IY,i,j)* mesh % faces(fID) % geom % normal(IY,i,j) &
                                                                    + mesh % faces(fID) % storage(1) % unStar(NCOMP,IZ,i,j)* mesh % faces(fID) % geom % normal(IZ,i,j)
               enddo ; enddo

            enddo
            !$acc end parallel loop 
!$omp end do 

            CALL BCs(zoneID) % bc % NeumannForEqn(mesh, mesh % zones(zoneID))   
            
!$omp do schedule(runtime) private(fID, i, j)  
            !$acc parallel loop gang vector_length(32) present(mesh, mesh % zones) copyin(zoneID) private(fID) async(1)
            do zonefID = 1, mesh % zones(zoneID) % no_of_faces
               fID =  mesh % zones(zoneID) % faces(zonefID)

!              ------------------------------------------------
!              Multiply by the Jacobian
!              ------------------------------------------------
               !$acc loop vector collapse(2)
               do j = 0, mesh % faces(fID) % Nf(2) ; do i = 0, mesh % faces(fID) % Nf(1)
                  mesh % faces(fID) % storage(1) % FStar(1,i,j) = (mesh % faces(fID) % storage(2) % FStar(1,i,j)) * &
                                                                   mesh % faces(fID) % geom % jacobian(i,j)
               end do ;  end do
               
!              ---------------------------
!              Return the flux to elements
!              ---------------------------
               call Face_ProjectFluxToElements(mesh % faces(fID), NCOMP, mesh % faces(fID) % storage(1) % FStar, 1)
            enddo
            !$acc end parallel loop 
!$omp end do 
!$acc wait
         enddo

      end subroutine Laplacian_computeBoundaryFlux

      SUBROUTINE ComputeTimeDerivativeIsolated( mesh, particles, time, mode, HO_Elements, Level)
         use EllipticDiscretizationClass
         IMPLICIT NONE 
!
!        ---------
!        Arguments
!        ---------
!
         TYPE(HexMesh), target           :: mesh
         type(Particles_t)               :: particles
         REAL(KIND=RP)                   :: time
         integer,             intent(in) :: mode
         logical, intent(in), optional   :: HO_Elements
         integer, intent(in), optional   :: Level

      end subroutine ComputeTimeDerivativeIsolated
      
      subroutine CompilerDefinedSourceTerm(mesh, time)
            use SMConstants
            USE HexMeshClass
            use PhysicsStorage
            use FluidData
            IMPLICIT NONE 
!
!           -------------------------------------------------------------
!           Aux. GPU source terms -- Alternative to UserDefinedSourceTerm
!           -------------------------------------------------------------
!
            type(HexMesh)               :: mesh
            real(kind=RP), intent(in)   :: time



      end subroutine CompilerDefinedSourceTerm

end module SpatialDiscretization
