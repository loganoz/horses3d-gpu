!
!//////////////////////////////////////////////////////
!
! This module stores the connection of each face of the wall that will be used in the 
! Wall Function with the first normal neighbour
! element and other needed variables
!
!//////////////////////////////////////////////////////
!
#include "Includes.h"
#if defined(NAVIERSTOKES)
Module WallFunctionConnectivity  !

    use SMConstants
    use HexMeshClass
    Implicit None
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
    public useWallFunc, wallFunBCs
!
!  ******************
!  Public definitions
!  ******************
!
    public Initialize_WallConnection, WallUpdateMeanV, WallStartMeanV, WallGetFaceConnectedQ, meanVelocity
!
    logical                                                          :: useWallFunc
    character(len=BC_STRING_LENGTH), dimension(:), allocatable       :: wallFunBCs
    integer, dimension(:), allocatable                               :: wallFaceOppositeID, wallSideID, wallFaceID
    real(kind=RP), dimension(:,:,:,:), allocatable                   :: meanVelocity
    real(kind=RP)                                                    :: timeCont
   
    !$acc declare create(meanVelocity, timeCont, wallFaceOppositeID, wallSideID, wallFaceID)

    contains 
!   
!------------------------------------------------------------------------------------------------------------------------
!
    Subroutine Initialize_WallConnection(controlVariables, mesh)
!     *******************************************************************
!        This subroutine initializes and store the arrays index that are
!        used to represent the connection of the face to the neighbour element.
!        Also calls the definitions of the wall function to update the 
!        parameter of the models based on the controlVariables.
!     *******************************************************************
!
        use FTValueDictionaryClass
        use Utilities,            only: toLower
        use FileReadingUtilities, only: getCharArrayFromString
        use ElementConnectivityDefinitions, only: normalAxis, FACES_PER_ELEMENT
        use MPI_Process_Info
        use WallFunctionDefinitions, only: Initialize_Wall_Function, wallFuncIndex, STD_WALL, ABL_WALL, u_tau0, useAverageV
        use Headers
#ifdef _HAS_MPI_
       use mpi
#endif
        implicit none
        class(FTValueDictionary),  intent(in)  :: controlVariables
        class(HexMesh), intent(inout)          :: mesh

!
!       ---------------
!       Local variables
!       ---------------
!
        character(len=LINE_LENGTH)             :: wallBC_str
        integer                                :: numberFacesWall, numberBC
        integer, dimension(:), allocatable     :: zonesWall
        integer                                :: i, j, nz, k, fID, efID
        integer                                :: linkedfID, linkedFaceSide, side, currentElementID, linkedElementglobID, normalIndex, oppositeIndex
        integer                                :: allFaces, ierr, element_index
        integer                                :: ghostCountLocal, ghostCountGlobal

        call Initialize_Wall_Function(controlVariables, useWallFunc)
        if (.not. useWallFunc) then
            return
        end if

        ghostCountLocal = 0
        
        ! get BC where the Wall Function will be applied
        wallBC_str = controlVariables % stringValueForKey("wall function boundaries", LINE_LENGTH)
        call toLower(wallBC_str)
        call getCharArrayFromString(wallBC_str, BC_STRING_LENGTH, wallFunBCs)

        ! get zones and number of faces for wall function
        numberBC = size(wallFunBCs)
        numberFacesWall = 0
        allocate(zonesWall(numberBC))

        do i = 1, numberBC
            do nz = 1, size(mesh % zones)
                if (trim(mesh % zones(nz) % Name) .eq. trim(wallFunBCs(i))) then
                    zonesWall(i) = nz
                    mesh % zones(nz) % useWallFunction = .true.
                    numberFacesWall = numberFacesWall + mesh % zones(nz) % no_of_faces
                    exit
                end if
            end do
        end do

        if (MPI_Process % doMPIAction) then
#ifdef _HAS_MPI_
            call mpi_allreduce(numberFacesWall, allFaces, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
        else
            allFaces = numberFacesWall
        end if
        if (allFaces .eq. 0) then
            useWallFunc = .false.
            if (MPI_Process % isRoot) write(STD_OUT,'(A)') "No wall BC found, the wall function will be deactivated"
            return
        end if

        allocate( wallFaceOppositeID(numberFacesWall), wallSideID(numberFacesWall), wallFaceID(numberFacesWall) )

        !get for each face of the wall, the linked element, normalDirection and index
        k = 0
     
        do j = 1, numberBC
            nz = zonesWall(j)
            do i = 1, mesh % zones(nz) % no_of_faces
                k = k + 1
                fID = mesh % zones(nz) % faces(i)
                
                currentElementID = mesh % faces(fID) % ElementIDs(1)
                associate ( e => mesh % elements(currentElementID) )
                    elem_loop:do efID = 1, FACES_PER_ELEMENT
                        if ( trim(wallFunBCs(j)) .eq. trim(e % boundaryName(efID)) ) then
                            normalIndex = normalAxis(efID)
                            exit elem_loop
                        end if
                    end do elem_loop
          
                    oppositeIndex = -1 * normalIndex
                    ! use the maxloc line if the compiler doesn't support findloc
                    oppositeIndex = maxloc(merge(1.0, 0.0, normalAxis == oppositeIndex),dim=1)
                    
                    linkedElementglobID = e % Connection(oppositeIndex) % globID
                    linkedfID = e % faceIDs(oppositeIndex)
                    
                    associate ( f => mesh % faces(linkedfID))
                        do side = 1, 2
                            element_index = f % elementIDs(side)
                            if (element_index .ne. 0) then
                                ! if not equal to zero means that the element exists at the same MPI partition
                                if ( linkedElementglobID .eq. mesh%elements(f % elementIDs(side))%globID) then
                                    linkedFaceSide = side
                                    exit
                                endif
                            else
                                ! if element_index is zero, then the element is in another MPI partition
                                ! so we accept it as the linked element (Ghost handling)
                                linkedFaceSide = side
                                ghostCountLocal = ghostCountLocal + 1
                            end if     
                        end do
                    end associate
                end associate

                wallFaceID(k) = fID
                wallFaceOppositeID(k) = linkedfID
                wallSideID(k) = linkedFaceSide

            end do
        end do

!
!       Initialize u_tau storage for the wall faces
!       -------------------------------------------
        do i = 1, numberFacesWall
            fID = wallFaceID(i)
            associate( f => mesh%faces(fID) )
                f % storage(1) % u_tau_NS = u_tau0
            end associate
        end do 

!
!       Allocate average velocity if requested, only for same p in all faces
!       -------------------------------------------
        if (useAverageV) then
            fID = wallFaceID(1)
            associate( f => mesh%faces(fID) )
                allocate( meanVelocity(NDIM,numberFacesWall,0:f % Nf(1),0:f % nf(2)) )
            end associate
            ! meanVelocity = 0.0_RP
            call WallStartMeanV(mesh)
            timeCont = 0.0_RP
            !$acc update device(meanVelocity, timeCont)
        end if

        !$acc update device(wallFaceID, wallFaceOppositeID, wallSideID)

!
!       Ghost faces count across all MPI processes
!       ------------------------------------------
#ifdef _HAS_MPI_
        if (MPI_Process % doMPIAction) then
            call mpi_allreduce(ghostCountLocal, ghostCountGlobal, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD, ierr)
        else
            ghostCountGlobal = ghostCountLocal
        end if
#else
        ghostCountGlobal = ghostCountLocal
#endif

       
!
!       Describe the Wall function
!       --------------------------
        if ( .not. MPI_Process % isRoot ) return
        write(STD_OUT,'(/)')
        call Subsection_Header("Wall function")

        write(STD_OUT,'(30X,A,A28,I0)') "->", "Number of faces: ", allFaces
        select case (wallFuncIndex)
            case (STD_WALL)
                write(STD_OUT,'(30X,A,A28,A)') "->", "Wall Function Law: ", "Reichardt"
            case (ABL_WALL)
                write(STD_OUT,'(30X,A,A28,A10)') "->", "Wall Function Law: ", "ABL"
        end select
        if (useAverageV) write(STD_OUT,'(30X,A,A28)') "->", "Wall Function using time averaged values"
        if (MPI_Process % isRoot .and. ghostCountGlobal > 0) then
            write(STD_OUT,'(/,30X,A)') "WARNING: Wall Function connectivity crosses MPI partitions"
            write(STD_OUT,'(32X,A,I8)') "Number of ghost-handled faces: ", ghostCountGlobal
            write(STD_OUT,'(32X,A)')   "Results may differ slightly from serial execution or from"
            write(STD_OUT,'(32X,A)')   "different partitionings."
            write(STD_OUT,'(32X,A)')   "Neighbour-face Q taken via ghost storage from previous step."
            write(STD_OUT,'(/)') 
        end if



    End Subroutine Initialize_WallConnection

    Subroutine WallGetFaceConnectedQ(mesh,f,Q,x,faceIndex,i,j)
        !$acc routine seq
!
!        *******************************************************************
!        This subroutine get the flow solution of the neighbour element
!        of the face.
!        *******************************************************************
        use PhysicsStorage
        use FaceClass
        type(HexMesh), intent(in)                                        :: mesh
        class(Face), intent(in)                                          :: f
        real(kind=RP), dimension(NCONS), intent(out)                     :: Q
        real(kind=RP), dimension(NDIM), intent(out)                      :: x
        integer, intent(in)                                              :: i, j
        integer, intent(in)                                              :: faceIndex
!
!       Local variables
        integer                                                          :: fID, side
        
        fID = wallFaceOppositeID(faceIndex)
        side = wallSideID(faceIndex)

        ! Accessing storage from the neighbor FACE instead of element
        Q = mesh % faces(fID) % storage(side) % Q(:,i,j)
        x = mesh % faces(fID) % geom % x(:,i,j)

    End Subroutine WallGetFaceConnectedQ

    Subroutine WallUpdateMeanV(mesh, dt)
        use PhysicsStorage
        use WallFunctionDefinitions, only: useAverageV
        implicit none
        type(HexMesh), intent(in)                       :: mesh
        real(kind=RP),  intent(in)                      :: dt
!
!       ---------------
!       Local variables
!       ---------------
!
        integer                                         :: fIndex, fID, i, j
        real(kind=RP), dimension(NCONS)                 :: Q
        real(kind=RP), dimension(NDIM)                  :: x
        real(kind=RP)                                   :: invRho
        real(kind=RP), dimension(NDIM)                  :: localV

        !
        ! create separate to set initial conditions
        if (.not. useAverageV) return

        !$acc parallel loop gang async(1) present(mesh) firstprivate(dt) copyin(timeCont) copyout(timeCont)
        do fIndex = 1, size(wallFaceID)
            fID = wallFaceID(fIndex)
            !$acc loop vector collapse(2) private(Q, x, localV)
            do j = 0, mesh % faces(fID) % Nf(2)
                do i = 0, mesh % faces(fID) % Nf(1)
                    call WallGetFaceConnectedQ(mesh, mesh%faces(fID), Q, x, fIndex, i, j)
                    invRho = 1.0_RP / Q(IRHO)
                    localV(:) = Q(IRHOU:IRHOW) * invRho
     
                    meanVelocity(:,fIndex,i,j) = ( (meanVelocity(:,fIndex,i,j) * timeCont) + localV(:) * dt ) / (timeCont+dt)
                end do
            end do
        end do 
        !$acc end parallel loop 

        timeCont = timeCont + dt

    End Subroutine WallUpdateMeanV

   
    Subroutine WallStartMeanV(mesh)
        use PhysicsStorage
        use WallFunctionDefinitions, only: useAverageV
        implicit none
        type(HexMesh), intent(in)                       :: mesh
!
!       ---------------
!       Local variables
!       ---------------
!
        integer                                         :: fIndex, fID, i, j
        real(kind=RP), dimension(NCONS)                 :: Q
        real(kind=RP), dimension(NDIM)                  :: x        
        real(kind=RP)                                   :: invRho
        real(kind=RP), dimension(NDIM)                  :: localV

        !
        ! create separate to set initial conditions
        if (.not. useAverageV) return

        do fIndex = 1, size(wallFaceID)
            fID = wallFaceID(fIndex)
            associate( f => mesh%faces(fID) )
                do j = 0, f % Nf(2)
                    do i = 0, f % Nf(1)
                        call WallGetFaceConnectedQ(mesh, f, Q, x, fIndex, i, j)
                        invRho = 1.0_RP / Q(IRHO)
                        localV(:) = Q(IRHOU:IRHOW) * invRho
       
                         meanVelocity(:,fIndex,i,j) = localV(:)
                    end do
                end do
            end associate
        end do 

    End Subroutine WallStartMeanV

End Module WallFunctionConnectivity
#endif