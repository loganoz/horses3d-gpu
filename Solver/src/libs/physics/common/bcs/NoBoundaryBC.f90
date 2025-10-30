#include "Includes.h"
module NoBoundaryBCClass
   use SMConstants
   use PhysicsStorage
   use VariableConversion
   use FileReaders,            only: controlFileName
   use FileReadingUtilities,   only: GetKeyword, GetValueAsString, PreprocessInputLine, CheckIfBoundaryNameIsContained
   use FTValueDictionaryClass, only: FTValueDictionary
   use GenericBoundaryConditionClass
   use FluidData
   use FileReadingUtilities, only: getRealArrayFromString
   use Utilities, only: toLower, almostEqual
   use HexMeshClass
   use ZoneClass
   implicit none
!
!  *****************************
!  Default everything to private
!  *****************************
!
   private
!
!  ****************
!  Public variables
!  ****************
!
!
!  ******************
!  Public definitions
!  ******************
!
   public NoBoundaryBC_t
!
!  ****************************
!  Static variables definitions
!  ****************************
!
!
!  ****************
!  Class definition
!  ****************
!
   type, extends(GenericBC_t) ::  NoBoundaryBC_t

      contains
         procedure         :: Destruct          => NoBoundaryBC_Destruct
         procedure         :: Describe          => NoBoundaryBC_Describe
#ifdef FLOW
         procedure         :: FlowState         => NoBoundaryBC_FlowState
         procedure         :: CreateDeviceData  => NoBoundaryBC_CreateDeviceData
         procedure         :: ExitDeviceData    => NoBoundaryBC_ExitDeviceData
         procedure         :: FlowGradVars      => NoBoundaryBC_FlowGradVars
         procedure         :: FlowNeumann       => NoBoundaryBC_FlowNeumann
#endif
#ifdef CAHNHILLIARD
         procedure         :: PhaseFieldState   => NoBoundaryBC_PhaseFieldState
         procedure         :: PhaseFieldGradVars=> NoBoundaryBC_PhaseFieldGradVars
         procedure         :: PhaseFieldNeumann => NoBoundaryBC_PhaseFieldNeumann
         procedure         :: ChemPotState      => NoBoundaryBC_ChemPotState
         procedure         :: ChemPotGradVars   => NoBoundaryBC_ChemPotGradVars
         procedure         :: ChemPotNeumann    => NoBoundaryBC_ChemPotNeumann
#endif
   end type NoBoundaryBC_t
!
!  *******************************************************************
!  Traditionally, constructors are exported with the name of the class
!  *******************************************************************
!
   interface NoBoundaryBC_t
      module procedure ConstructNoBoundaryBC
   end interface NoBoundaryBC_t
!
!  *******************
!  Function prototypes
!  *******************
!
   abstract interface
   end interface
!
!  ========
   contains
!  ========
!
!/////////////////////////////////////////////////////////
!
!        Class constructor
!        -----------------
!
!/////////////////////////////////////////////////////////
!
      function ConstructNoBoundaryBC(bname)
!
!        ********************************************************************
!        · Definition of the noSlipWall boundary condition in the control file:
!              #define boundary bname
!                 type             = noSlipWall
!                 velocity         = #value        (only in incompressible NS)
!                 Mach number      = #value        (only in compressible NS)
!                 AoAPhi           = #value
!                 AoATheta         = #value
!                 density          = #value        (only in monophase)
!                 pressure         = #value        (only in compressible NS)
!                 multiphase type  = mixed/layered
!                 phase 1 layer x  > #value
!                 phase 1 layer y  > #value
!                 phase 1 layer z  > #value
!                 phase 1 velocity = #value
!                 phase 2 velocity = #value
!              #end
!        ********************************************************************
!
         implicit none
         type(NoBoundaryBC_t)             :: ConstructNoBoundaryBC
         character(len=*), intent(in) :: bname
!
!        ---------------
!        Local variables
!        ---------------
!
         integer        :: fid, io
         character(len=LINE_LENGTH) :: boundaryHeader
         character(len=LINE_LENGTH) :: currentLine
         character(len=LINE_LENGTH) :: keyword, keyval
         logical                    :: inside
         type(FTValueDIctionary)    :: bcdict

         open(newunit = fid, file = trim(controlFileName), status = "old", action = "read")

         call bcdict % InitWithSize(16)

         ConstructNoBoundaryBC % BCType = "noboundary"
         ConstructNoBoundaryBC % bname  = bname
         call toLower(ConstructNoBoundaryBC % bname)

         write(boundaryHeader,'(A,A)') "#define boundary ",trim(bname)
         call toLower(boundaryHeader)
!
!        Navigate until the "#define boundary bname" sentinel is found
!        -------------------------------------------------------------
         inside = .false.
         do 
            read(fid, '(A)', iostat=io) currentLine

            IF(io .ne. 0 ) EXIT

            call PreprocessInputLine(currentLine)
            call toLower(currentLine)

            if ( index(trim(currentLine),"#define boundary") .ne. 0 ) then
               inside = CheckIfBoundaryNameIsContained(trim(currentLine), trim(ConstructNoBoundaryBC % bname)) 
            end if
!
!           Get all keywords inside the zone
!           --------------------------------
            if ( inside ) then
               if ( trim(currentLine) .eq. "#end" ) exit

               keyword  = ADJUSTL(GetKeyword(currentLine))
               keyval   = ADJUSTL(GetValueAsString(currentLine))
               call ToLower(keyword)
      
               call bcdict % AddValueForKey(keyval, trim(keyword))

            end if

         end do

         close(fid)
         call bcdict % Destruct

      end function ConstructNoBoundaryBC

      subroutine NoBoundaryBC_Describe(self)
!
!        ***************************************************
!              Describe the inflow boundary condition
         implicit none
         class(NoBoundaryBC_t),  intent(in)  :: self
         write(STD_OUT,'(30X,A,A28,A)') "->", " Boundary condition type: ", "NoBoundary"
         
      end subroutine NoBoundaryBC_Describe

!
!/////////////////////////////////////////////////////////
!
!        Class destructors
!        -----------------
!
!/////////////////////////////////////////////////////////
!
      subroutine NoBoundaryBC_Destruct(self)
         implicit none
         class(NoBoundaryBC_t)    :: self

      end subroutine NoBoundaryBC_Destruct
!
!////////////////////////////////////////////////////////////////////////////
!
!        Subroutines for all flow equations
!        -----------------------------------
!
!////////////////////////////////////////////////////////////////////////////
!
#if defined(FLOW)

      subroutine NoBoundaryBC_CreateDeviceData(self)
         implicit none 
         class(NoBoundaryBC_t), intent(in)    :: self

         !$acc enter data copyin(self)

      end subroutine NoBoundaryBC_CreateDeviceData

      subroutine NoBoundaryBC_ExitDeviceData(self)
         implicit none 
         class(NoBoundaryBC_t), intent(in)    :: self

         !$acc exit data delete(self)

      end subroutine NoBoundaryBC_ExitDeviceData
#endif
!
!////////////////////////////////////////////////////////////////////////////
!
!        Subroutines for compressible Navier--Stokes equations
!        -----------------------------------------------------
!
!////////////////////////////////////////////////////////////////////////////
!
#ifdef NAVIERSTOKES


      subroutine NoBoundaryBC_FlowState(self, mesh, zone)
!
!        *************************************************************
!           Compute the state variables for a general wall
!
!           · Density is computed from the interior state
!           · Wall velocity is set to v_interior - 2v_normal:
!                 -> Maintains tangential speed.
!                 -> Removes normal speed (weakly).
!           · Internal energy is either the interior state for 
!              adiabatic walls, or the imposed for isothermal.
!              eBC = eInt + kWallType (eIso - eInt)
!              where kWallType = 0 for adiabatic and 1 for isothermal.        
!        *************************************************************
!
         use HexMeshClass
         implicit none
         class(NoBoundaryBC_t), intent(in)    :: self
         type(HexMesh), intent(inout)           :: mesh
         type(Zone_t), intent(in)               :: zone

!         integer,                 intent(in)    :: zoneID                              
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: qNorm, pressure_aux
         real(kind=RP) :: Q(NCONS)
         integer       :: i,j,zonefID,fID
         
         !$acc parallel loop gang present(mesh, self, zone) async(1)
         do zonefID = 1, zone % no_of_faces
            fID = zone % faces(zonefID)
            !$acc loop vector collapse(2) private(Q)            
            do j = 0, mesh % faces(fID) % Nf(2)  ; do i = 0, mesh % faces(fID) % Nf(1)
               
               Q = mesh % faces(fID) % storage(1) % Q(:,i,j)

               qNorm = mesh % faces(fID) % geom % normal(IX,i,j) * Q(IRHOU) + &
                       mesh % faces(fID) % geom % normal(IY,i,j) * Q(IRHOV) + &
                       mesh % faces(fID) % geom % normal(IZ,i,j) * Q(IRHOW) 
         
               Q(IRHOU:IRHOW) = Q(IRHOU:IRHOW) - 2.0_RP * qNorm * mesh % faces(fID) % geom % normal(:,i,j)
               
               mesh % faces(fID) % storage(2) % Q(1:5,i,j) = Q(1:5)
               mesh % faces(fID) % storage(2) % Q(2,i,j) = mesh % faces(fID) % storage(1) % Q(2,i,j)
               
            enddo ; enddo
         enddo
         !$acc end parallel loop

      end subroutine NoBoundaryBC_FlowState

      subroutine NoBoundaryBC_FlowGradVars(self, mesh, zone)
!
!        *****************************************************************
!           Only set the temperature, velocity is Neumann, use interior!
!        *****************************************************************
!
         implicit none
         class(NoBoundaryBC_t), intent(in)    :: self
         type(HexMesh), intent(inout)           :: mesh
         type(Zone_t), intent(in)               :: zone

!         integer,                 intent(in)    :: zoneID 
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)  :: rhou_n
         real(kind=RP)  :: Q_aux(NCONS),Q(NCONS)
         real(kind=RP)  :: u_int(NGRAD), u_star(NGRAD)
         integer        :: i,j,zonefID,fID
         !$acc parallel loop gang present(mesh, self, zone) private(fID) async(1)
         do zonefID = 1, zone % no_of_faces
            fID = zone % faces(zonefID)
            !$acc loop vector collapse(2) private(Q, Q_aux, u_star, u_int)            
            do j = 0, mesh % faces(fID) % Nf(2)  ; do i = 0, mesh % faces(fID) % Nf(1)

               mesh % faces(fID) % storage(1) % unStar(:,1,i,j) = 0.0_RP
               mesh % faces(fID) % storage(1) % unStar(:,2,i,j) = 0.0_RP    
               mesh % faces(fID) % storage(1) % unStar(:,3,i,j) = 0.0_RP

            enddo ; enddo
         enddo
         !$acc end parallel loop           
      end subroutine NoBoundaryBC_FlowGradVars

      subroutine NoBoundaryBC_FlowNeumann(self, mesh, zone)
!
!        ***********************************************************
!        ***********************************************************
!
         implicit none
         class(NoBoundaryBC_t), intent(in)    :: self
         type(HexMesh), intent(inout)           :: mesh
         type(Zone_t), intent(in)               :: zone
!
!        ---------------
!        Local Variables
!        ---------------
!   
         integer        :: i,j,zonefID,fID
         real(kind=RP)  :: viscWork, heatFlux
         real(kind=RP)  :: flux(NCONS),Q(NCONS)

         !$acc parallel loop gang present(mesh, self, zone) private(fID) async(1)
         do zonefID = 1, zone % no_of_faces
            fID = zone % faces(zonefID)
            !$acc loop vector collapse(2) private(Q, flux, viscWork, heatFlux)     
            do j = 0, mesh % faces(fID) % Nf(2)  ; do i = 0, mesh % faces(fID) % Nf(1)

               mesh % faces(fID) % storage(2) % FStar(1:5,i,j) = 0.0_RP
            enddo ; enddo
         enddo
         !$acc end parallel loop  

      end subroutine NoBoundaryBC_FlowNeumann
#endif
!
!////////////////////////////////////////////////////////////////////////////
!
!        Subroutines for incompressible Navier--Stokes equations
!        -------------------------------------------------------
!
!////////////////////////////////////////////////////////////////////////////
!
#ifdef INCNS
      subroutine NoBoundaryBC_FlowState(self, x, t, nHat, Q)
!
!        *************************************************************
!           Compute the state variables for a general wall
!
!           · Density is computed from the interior state
!           · Wall velocity is set to 2v_wall - v_interior
!           · Pressure is computed from the interior state
!        *************************************************************
!

         implicit none
         class(NoBoundaryBC_t),  intent(in)    :: self
         real(kind=RP),       intent(in)    :: x(NDIM)
         real(kind=RP),       intent(in)    :: t
         real(kind=RP),       intent(in)    :: nHat(NDIM)
         real(kind=RP),       intent(inout) :: Q(NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)  :: vn
!
!        -----------------------------------------------
!        Generate the external flow along the face, that
!        represents a solid wall.
!        -----------------------------------------------
!
         vn = sum(Q(INSRHOU:INSRHOW)*nHat)
                                                                                                    

         Q(INSRHO)          = Q(INSRHO)
         Q(INSRHOU:INSRHOW) = Q(INSRHOU:INSRHOW) - 2.0_RP * vn * nHat
         Q(INSP)            = Q(INSP)

      end subroutine NoBoundaryBC_FlowState

      subroutine NoBoundaryBC_FlowGradVars(self, x, t, nHat, Q, U, GetGradients)
!
!        **************************************************************
!           Use the interior velocity: Neumann BC!
!        **************************************************************
!
         implicit none
         class(NoBoundaryBC_t),  intent(in)  :: self
         real(kind=RP),          intent(in)    :: x(NDIM)
         real(kind=RP),          intent(in)    :: t
         real(kind=RP),          intent(in)    :: nHat(NDIM)
         real(kind=RP),          intent(in)    :: Q(NCONS)
         real(kind=RP),          intent(inout) :: U(NGRAD)
         procedure(GetGradientValues_f)        :: GetGradients

      end subroutine NoBoundaryBC_FlowGradVars

      subroutine NoBoundaryBC_FlowNeumann(self, x, t, nHat, Q, U_x, U_y, U_z, flux)
!
!        ***************************************************************
!           Set homogeneous Neumann BCs everywhere
!        ***************************************************************
!        
         implicit none
         class(NoBoundaryBC_t),  intent(in) :: self
         real(kind=RP),       intent(in)      :: x(NDIM)
         real(kind=RP),       intent(in)      :: t
         real(kind=RP),       intent(in)      :: nHat(NDIM)
         real(kind=RP),       intent(in)      :: Q(NCONS)
         real(kind=RP),       intent(in)      :: U_x(NCONS)
         real(kind=RP),       intent(in)      :: U_y(NCONS)
         real(kind=RP),       intent(in)      :: U_z(NCONS)
         real(kind=RP),       intent(inout)   :: flux(NCONS)


      end subroutine NoBoundaryBC_FlowNeumann
#endif
!
!////////////////////////////////////////////////////////////////////////////
!
!        Subroutines for multiphase solver
!        ---------------------------------
!
!////////////////////////////////////////////////////////////////////////////
!
#ifdef MULTIPHASE

      subroutine NoBoundaryBC_FlowState(self, mesh, zone)
!
!        *************************************************************
!           Compute the state variables for a general wall
!
!           · Density is computed from the interior state
!           · Wall velocity is set to 2v_wall - v_interior
!           · Pressure is computed from the interior state
!        *************************************************************
!

         use HexMeshClass
         implicit none
         class(NoBoundaryBC_t), intent(in)    :: self
         type(HexMesh), intent(inout)              :: mesh
         type(Zone_t), intent(in)               :: zone 
!
!        ---------------
!        Local variables
!        ---------------
!
        real(RP) :: q1, q2, q3, q4, q5
        real(RP) :: vn
        integer :: fID
        integer :: i, j, zonefID
        ! local pointers/aliases (Fortran pointer syntax omitted — we copy values/elements to locals)
        ! Note: adjust vector_length/gang/workers to your hardware.

        !$acc parallel loop gang collapse(1) present(mesh, self, zone) async(1)
        do zonefID = 1, zone%no_of_faces
           fID = zone%faces(zonefID)

           ! optionally, copy pointers/offsets to local variables if your compiler supports
           ! For portability, we'll use direct indexing but minimize repeated long expressions.

           !$acc loop vector collapse(2) private(i,j,q1,q2,q3,q4,q5,vn)
           do j = 0, mesh%faces(fID)%Nf(2)
              do i = 0, mesh%faces(fID)%Nf(1)

                 ! load the 5-component state into locals
                 q1 = mesh%faces(fID)%storage(1)%Q(1,i,j)
                 q2 = mesh%faces(fID)%storage(1)%Q(2,i,j)
                 q3 = mesh%faces(fID)%storage(1)%Q(3,i,j)
                 q4 = mesh%faces(fID)%storage(1)%Q(4,i,j)
                 q5 = mesh%faces(fID)%storage(1)%Q(5,i,j)

                 ! compute vn using local scalars and direct normals (load normals once)
                 vn = mesh%faces(fID)%geom%normal(1,i,j) * q2 + &
                      mesh%faces(fID)%geom%normal(2,i,j) * q3 + &
                      mesh%faces(fID)%geom%normal(3,i,j) * q4

                 ! apply wall reflection update to local q's
                 ! modify third and fourth as in original code:
                 !q2 = q2 - 2.0_RP * vn * mesh%faces(fID)%geom%normal(1,i,j)
                 q3 = q3 - 2.0_RP * vn * mesh%faces(fID)%geom%normal(2,i,j)
                 q4 = q4 - 2.0_RP * vn * mesh%faces(fID)%geom%normal(3,i,j)

                 ! single write-back of the slice
                 mesh%faces(fID)%storage(2)%Q(1,i,j) = q1
                 mesh%faces(fID)%storage(2)%Q(2,i,j) = q2
                 mesh%faces(fID)%storage(2)%Q(3,i,j) = q3
                 mesh%faces(fID)%storage(2)%Q(4,i,j) = q4
                 mesh%faces(fID)%storage(2)%Q(5,i,j) = q5

                 ! copy mu (single scalar)
                 mesh%faces(fID)%storage(2)%mu(1,i,j) = mesh%faces(fID)%storage(1)%mu(1,i,j)

              end do
           end do
        end do
        !$acc end parallel loop

      end subroutine NoBoundaryBC_FlowState

      subroutine NoBoundaryBC_FlowGradVars(self, mesh, zone)
         !TODO slightly unsure about this one 
!
!        **************************************************************
!           Use the interior velocity: Neumann BC!
!        **************************************************************
!
         implicit none
         class(NoBoundaryBC_t), intent(in)    :: self
         type(HexMesh), intent(inout)           :: mesh
         type(Zone_t), intent(in)               :: zone

!
!        ---------------
!        Local variables
!        ---------------
!        
         real(kind=RP)  :: Q(NCONS), rho
         integer        :: i,j,zonefID,fID, m
         real(kind=RP)  :: u_int(NGRAD), u_star(NGRAD)

         !$acc parallel loop gang present(mesh, self, zone) private(fID) async(1)
         do zonefID = 1, zone % no_of_faces
            fID = zone % faces(zonefID)
            !$acc loop vector collapse(3) private(Q, u_int)            
            do j = 0, mesh % faces(fID) % Nf(2)  ; do i = 0, mesh % faces(fID) % Nf(1) ; do m=1, NCONS

               mesh % faces(fID) % storage(1) % unStar(m,1,i,j) = 0.0_RP
               mesh % faces(fID) % storage(1) % unStar(m,2,i,j) = 0.0_RP  
               mesh % faces(fID) % storage(1) % unStar(m,3,i,j) = 0.0_RP
               
            enddo ; enddo ; enddo
         enddo
         !$acc end parallel loop
         
      end subroutine NoBoundaryBC_FlowGradVars

      subroutine NoBoundaryBC_FlowNeumann(self, mesh, zone)
         implicit none
         class(NoBoundaryBC_t), intent(in)    :: self
         type(HexMesh), intent(inout)              :: mesh
         type(Zone_t), intent(in)               :: zone

         integer        :: i,j,zonefID,fID, m
         real(kind=RP)  :: flux(NCONS),Q(NCONS)


         !$acc parallel loop gang present(mesh, self, zone) private(fID) async(1)
         do zonefID = 1, zone % no_of_faces
            fID = zone % faces(zonefID)
            !$acc loop vector collapse(3) independent private(Q,flux)  
            do j = 0, mesh % faces(fID) % Nf(2) ; do i = 0, mesh % faces(fID) % Nf(1) ; do m=1, NCONS
               mesh % faces(fID) % storage(2) % FStar(m,i,j) = 0.0_RP
            enddo 
            enddo
          enddo
         enddo
         !$acc end parallel loop

      end subroutine NoBoundaryBC_FlowNeumann
#endif
!
!////////////////////////////////////////////////////////////////////////////
!
!        Subroutines for Cahn--Hilliard
!        ------------------------------
!
!////////////////////////////////////////////////////////////////////////////
!
#if defined(CAHNHILLIARD)


      subroutine NoBoundaryBC_PhaseFieldState(self, mesh, zone)
         use HexMeshClass
         implicit none
         class(NoBoundaryBC_t), intent(in)    :: self
         type(HexMesh), intent(inout)           :: mesh
         type(Zone_t), intent(in)               :: zone
!
!        ---------------
!        Local variables
!        ---------------
!
        integer :: Nf1, Nf2
        integer :: i, j, zonefID, fID

        ! ensure mesh,zone are already on device
        !$acc parallel loop gang present(mesh, self, zone) async(1)
        do zonefID = 1, zone%no_of_faces
           fID = zone%faces(zonefID)

           ! copy extents to local scalars (cheap)
           Nf1 = mesh%faces(fID)%Nf(1)
           Nf2 = mesh%faces(fID)%Nf(2)

           !$acc loop vector collapse(2) private(i,j)
           do j = 0, Nf2
              do i = 0, Nf1

                 ! direct contiguous slice copy — avoids extra private array and temporaries
                 mesh%faces(fID)%storage(2)%Q(1,i,j) = mesh%faces(fID)%storage(1)%Q(1,i,j)
                 mesh%faces(fID)%storage(2)%Q(2,i,j) = mesh%faces(fID)%storage(1)%Q(2,i,j)
                 mesh%faces(fID)%storage(2)%Q(3,i,j) = mesh%faces(fID)%storage(1)%Q(3,i,j)
                 mesh%faces(fID)%storage(2)%Q(4,i,j) = mesh%faces(fID)%storage(1)%Q(4,i,j)
                 mesh%faces(fID)%storage(2)%Q(5,i,j) = mesh%faces(fID)%storage(1)%Q(5,i,j)

                 ! if mu needs copying:
                 mesh%faces(fID)%storage(2)%mu(1,i,j) = mesh%faces(fID)%storage(1)%mu(1,i,j)

              end do
           end do
        end do
        !$acc end parallel loop
      end subroutine NoBoundaryBC_PhaseFieldState

      subroutine NoBoundaryBC_PhaseFieldGradVars(self, mesh, zone)
         implicit none
         class(NoBoundaryBC_t), intent(in)    :: self
         type(HexMesh), intent(inout)           :: mesh
         type(Zone_t), intent(in)               :: zone

!
!        ---------------
!        Local variables
!        ---------------
!        
         integer        :: i,j,zonefID,fID,m

         !$acc parallel loop gang present(mesh, self, zone) async(1)
         do zonefID = 1, zone % no_of_faces
            fID = zone % faces(zonefID)
            !$acc loop vector collapse(3)            
            do j = 0, mesh % faces(fID) % Nf(2)  ; do i = 0, mesh % faces(fID) % Nf(1) ; do m=1, NCONS
               
               mesh % faces(fID) % storage(1) % unStar(m,1,i,j) = 0.0_RP
               mesh % faces(fID) % storage(1) % unStar(m,2,i,j) = 0.0_RP
               mesh % faces(fID) % storage(1) % unStar(m,3,i,j) = 0.0_RP
               
            enddo ; enddo ; enddo
         enddo
         !$acc end parallel loop
         
      end subroutine NoBoundaryBC_PhaseFieldGradVars

      subroutine NoBoundaryBC_PhaseFieldNeumann(self, mesh, zone)
         implicit none
         class(NoBoundaryBC_t), intent(in)    :: self
         type(HexMesh), intent(inout)           :: mesh
         type(Zone_t), intent(in)               :: zone
!
!        ---------------
!        Local variables
!        ---------------
!  
         real(kind=RP), parameter   :: MIN_ = 1.0e-1_RP
         real(kind=RP)              :: prod, Q(NCOMP)
         integer                    :: i,j,zonefID,fID, m
         
         !$acc parallel loop gang present(mesh, self, zone) private(fID) async(1)
         do zonefID = 1, zone % no_of_faces
            fID = zone % faces(zonefID)
            !$acc loop vector collapse(3) independent private(Q, prod)
            do j = 0, mesh % faces(fID) % Nf(2) ; do i = 0, mesh % faces(fID) % Nf(1) ; do m=1, NCONS
             
               mesh % faces(fID) % storage(2) % FStar(m,i,j) = 0.0_RP
            enddo 
            enddo
          enddo
         enddo
         !$acc end parallel loop

      end subroutine NoBoundaryBC_PhaseFieldNeumann

      subroutine NoBoundaryBC_ChemPotState(self, mesh, zone)
         use HexMeshClass
         implicit none
         class(NoBoundaryBC_t), intent(in)    :: self
         type(HexMesh), intent(inout)           :: mesh
         type(Zone_t), intent(in)               :: zone
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: Q(NCONS)
         integer       :: i,j,zonefID,fID,m
         
         !$acc parallel loop gang present(mesh, self, zone) async(1)
         do zonefID = 1, zone % no_of_faces
            fID = zone % faces(zonefID)
            !$acc loop vector collapse(3) private(Q)            
            do j = 0, mesh % faces(fID) % Nf(2)  ; do i = 0, mesh % faces(fID) % Nf(1) ; do m=1,NCONS

               mesh % faces(fID) % storage(2) % Q(m,i,j) = mesh % faces(fID) % storage(1) % Q(m,i,j)
            
            enddo ; enddo ; enddo
         enddo
         !$acc end parallel loop
      end subroutine NoBoundaryBC_ChemPotState

      subroutine NoBoundaryBC_ChemPotGradVars(self, mesh, zone)
         implicit none
         class(NoBoundaryBC_t), intent(in)    :: self
         type(HexMesh), intent(inout)           :: mesh
         type(Zone_t), intent(in)               :: zone
!
!        ---------------
!        Local variables
!        ---------------
!        
         integer        :: i,j,zonefID,fID, m

         !$acc parallel loop gang present(mesh, self, zone) private(fID) async(1)
         do zonefID = 1, zone % no_of_faces
            fID = zone % faces(zonefID)
            !$acc loop vector collapse(3)            
            do j = 0, mesh % faces(fID) % Nf(2)  ; do i = 0, mesh % faces(fID) % Nf(1) ; do m=1,NCONS
               
               mesh % faces(fID) % storage(1) % unStar(m,1,i,j) = 0.0
               mesh % faces(fID) % storage(1) % unStar(m,2,i,j) = 0.0    
               mesh % faces(fID) % storage(1) % unStar(m,3,i,j) = 0.0
               
            enddo ; enddo; enddo
         enddo
         !$acc end parallel loop
         
      end subroutine NoBoundaryBC_ChemPotGradVars

      subroutine NoBoundaryBC_ChemPotNeumann(self, mesh, zone)
         implicit none
         class(NoBoundaryBC_t), intent(in)    :: self
         type(HexMesh), intent(inout)              :: mesh
         type(Zone_t), intent(in)               :: zone

         integer        :: i,j,zonefID,fID, m
         real(kind=RP)  :: flux(NCONS),Q(NCONS)

         !$acc parallel loop gang present(mesh, self, zone) private(fID) async(1)
         do zonefID = 1, zone % no_of_faces
            fID = zone % faces(zonefID)
            !$acc loop vector collapse(3) independent private(Q,flux)  
            do j = 0, mesh % faces(fID) % Nf(2) ; do i = 0, mesh % faces(fID) % Nf(1) ; do m=1,NCONS
               mesh % faces(fID) % storage(2) % FStar(m,i,j) = 0.0_RP
            enddo 
            enddo
          enddo
         enddo
         !$acc end parallel loop

      end subroutine NoBoundaryBC_ChemPotNeumann

#endif
!
end module NoBoundaryBCClass
