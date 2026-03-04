#include "Includes.h"
module MeshPartitioning
   use SMConstants
   use MeshTypes
   use HexMeshClass
   use PartitionedMeshClass
   use FileReadingUtilities            , only: RemovePath, getFileName

   private
   public   PerformMeshPartitioning

   contains
      subroutine PerformMeshPartitioning(mesh, no_of_allElements, no_of_domains, partitions, useWeights, Nx, Ny, Nz)
         implicit none
         type(HexMesh), intent(in)  :: mesh
         integer,       intent(in)  :: no_of_allElements
         integer,       intent(in)  :: no_of_domains
         type(PartitionedMesh_t)    :: partitions(no_of_domains)
         logical,       intent(in)  :: useWeights
         integer,       intent(in)  :: Nx(no_of_allElements), Ny(no_of_allElements), Nz(no_of_allElements)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer               :: fID, domain
         integer               :: elementsDomain(no_of_allElements)
!
!        Initialize partitions
!        ---------------------
         do domain = 1, no_of_domains
            partitions(domain) = PartitionedMesh_t(domain)
         end do
!
!        Get each domain elements and nodes
!        ----------------------------------
         call GetElementsDomain(mesh, no_of_allElements, no_of_domains, elementsDomain, partitions, useWeights, Nx, Ny, Nz)
!
!        Get the partition boundary faces
!        --------------------------------
         call GetPartitionBoundaryFaces(mesh, no_of_domains, elementsDomain, partitions)

!
!        Export partitions file
!        ----------------------
         call WritePartitionsFile(mesh, elementsDomain)
         call WriteMeshPartitioningFile(mesh, no_of_domains, partitions)

      end subroutine PerformMeshPartitioning

      subroutine GetElementsDomain(mesh, no_of_allElements, no_of_domains, elementsDomain, partitions, useWeights, Nx, Ny, Nz)
         use IntegerDataLinkedList
         use MPI_Process_Info
         implicit none
         type(HexMesh), intent(in)              :: mesh
         integer,       intent(in)              :: no_of_allElements
         integer,       intent(in)              :: no_of_domains
         integer,       intent(out)             :: elementsDomain(mesh % no_of_elements)
         type(PartitionedMesh_t), intent(inout) :: partitions(no_of_domains)      
         logical,       intent(in)  :: useWeights
         integer,       intent(in)  :: Nx(no_of_allElements), Ny(no_of_allElements), Nz(no_of_allElements)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer, allocatable :: nodesDomain(:)


         allocate(nodesDomain(size(mesh % nodes)))
!
!        **********************************************
!        Set which elements belong to each domain using
!        **********************************************
!
         select case (MPI_Partitioning)
!     
!           Space-filling curve partitioning
!           --------------------------------
            case (SFC_PARTITIONING)
               call GetSFCElementsPartition(no_of_domains, no_of_allElements, elementsDomain, useWeights, Nx, Ny, Nz)
!     
!           METIS partitioning
!           ------------------
            case (METIS_PARTITIONING)
               call GetMETISElementsPartition(mesh, no_of_domains, elementsDomain, nodesDomain, useWeights)
         end select
!
!        ****************************************
!        Get which nodes belong to each partition
!        ****************************************
!
         call GetNodesPartition(mesh, no_of_domains, elementsDomain, partitions)   

         deallocate(nodesDomain)   
         
      end subroutine GetElementsDomain

!
!////////////////////////////////////////////////////////////////////////
!
      subroutine GetNodesPartition(mesh, no_of_domains, elementsDomain, partitions)
         use Utilities, only: Qsort
         implicit none
     
         type(HexMesh), intent(in)              :: mesh
         integer,       intent(in)              :: no_of_domains
         integer,       intent(in)              :: elementsDomain(mesh%no_of_elements)
         type(PartitionedMesh_t), intent(inout) :: partitions(no_of_domains)
     
         integer              :: nvertex, i, j, k, ipoint, jpoint
         integer              :: idomain, npoints, ielem
         logical              :: isnewpoint, meshIsHOPR
         integer, allocatable :: points(:), HOPRpoints(:)
         logical, allocatable :: mask(:)  ! Mask array for node uniqueness
         integer              :: max_nodes
     
         nvertex = 8
         meshIsHOPR = allocated(mesh%HOPRnodeIDs)
         max_nodes = 0
         do ielem = 1, mesh%no_of_elements
            max_nodes = max(maxval(mesh%elements(ielem)%nodeIDs(:)), max_nodes)
         end do

         allocate(mask(max_nodes))  ! One-time allocation
     
         do idomain = 1, no_of_domains
             partitions(idomain)%no_of_elements = count(elementsDomain == idomain)
             allocate(partitions(idomain)%elementIDs(partitions(idomain)%no_of_elements))
             allocate(points(nvertex * partitions(idomain)%no_of_elements))
             points = 0
     
             if (meshIsHOPR) then
                 allocate(HOPRpoints(nvertex * partitions(idomain)%no_of_elements))
                 HOPRpoints = 0
             end if
     
             mask = .false.  ! Reset mask for new domain
             k = 0
             npoints = 0
             do ielem = 1, mesh%no_of_elements
                 if (elementsDomain(ielem) == idomain) then
                     k = k + 1
                     partitions(idomain)%elementIDs(k) = ielem
                     do j = 1, nvertex
                         jpoint = mesh%elements(ielem)%nodeIDs(j)
                         if (.not. mask(jpoint)) then
                             npoints = npoints + 1
                             points(npoints) = jpoint
                             mask(jpoint) = .true.
                             if (meshIsHOPR) HOPRpoints(npoints) = mesh%HOPRnodeIDs(jpoint)
                         end if
                     end do
                 end if
             end do
     
             allocate(partitions(idomain)%nodeIDs(npoints))
             partitions(idomain)%nodeIDs = points(1:npoints)
     
             if (meshIsHOPR) then
                 allocate(partitions(idomain)%HOPRnodeIDs(npoints))
                 partitions(idomain)%HOPRnodeIDs = HOPRpoints(1:npoints)
             end if
     
             if (.not. meshIsHOPR) call Qsort(partitions(idomain)%nodeIDs)
     
             partitions(idomain)%no_of_nodes = npoints
     
             deallocate(points)
             if (allocated(HOPRpoints)) deallocate(HOPRpoints)
         end do
     
         deallocate(mask)
     end subroutine GetNodesPartition
!
!////////////////////////////////////////////////////////////////////////
!
      subroutine GetPartitionBoundaryFaces(mesh, no_of_domains, elementsDomain, partitions)
 
         implicit none
         type(HexMesh), intent(in)  :: mesh
         integer,       intent(in)  :: no_of_domains
         integer,       intent(in)  :: elementsDomain(mesh % no_of_elements)
         type(PartitionedMesh_t)    :: partitions(no_of_domains)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: fID, eL, eR, dL, dR, domain
         integer  :: bfaceID(no_of_domains)
!
!        *******************************************
!        Get the number of boundary faces per domain 
!        *******************************************
!
         do fID = 1, size(mesh % faces)
            associate(f => mesh % faces(fID))
!
!           Cycle non-interior faces
!           ------------------------
            if ( f % faceType .ne. HMESH_INTERIOR ) cycle
!
!           Create references to left and right elements
!           --------------------------------------------
            associate( eL => mesh % elements(f % elementIDs(1)), &
                       eR => mesh % elements(f % elementIDs(2)))
!
!           Get each elements domain
!           ------------------------
            dL = elementsDomain(eL % eID)
            dR = elementsDomain(eR % eID)
!
!           Cycle if both elements belong to the same domain
!           ------------------------------------------------
            if ( dL .eq. dR ) cycle 
!
!           Otherwise, the face is a domain boundary face for domains dL and dR
!           -------------------------------------------------------------------
            partitions(dL) % no_of_mpifaces = partitions(dL) % no_of_mpifaces + 1
            partitions(dR) % no_of_mpifaces = partitions(dR) % no_of_mpifaces + 1

            end associate
            end associate
         end do
!
!        **************************************
!        Allocate boundary faces-related memory
!        **************************************
!
         do domain = 1, no_of_domains
            associate(nFaces => partitions(domain) % no_of_mpifaces)
            allocate(partitions(domain) % mpiface_elements(nFaces))
            allocate(partitions(domain) % element_mpifaceSide(nFaces))
            allocate(partitions(domain) % element_mpifaceSideOther(nFaces))
            allocate(partitions(domain) % mpiface_rotation(nFaces))
            allocate(partitions(domain) % mpiface_elementSide(nFaces))
            allocate(partitions(domain) % mpiface_sharedDomain(nFaces))
            end associate
         end do
!
!        ***************************
!        Get each boundary face data
!        ***************************
!
         bfaceID = 0
         do fID = 1, size(mesh % faces)
            associate(f => mesh % faces(fID))
!
!           Cycle non-interior faces
!           ------------------------
            if ( f % faceType .ne. HMESH_INTERIOR ) cycle
!
!           Create references to left and right elements
!           --------------------------------------------
            associate( eL => mesh % elements(f % elementIDs(1)), &
                       eR => mesh % elements(f % elementIDs(2)))
!
!           Get each elements domain
!           ------------------------
            dL = elementsDomain(eL % eID)
            dR = elementsDomain(eR % eID)
!
!           Cycle if both elements belong to the same domain
!           ------------------------------------------------
            if ( dL .eq. dR ) cycle 
!
!           Otherwise, the face is a domain boundary face for domains dL and dR
!           -------------------------------------------------------------------
            bfaceID(dL) = bfaceID(dL) + 1 
            bfaceID(dR) = bfaceID(dR) + 1 
!
!           Get the elements
!           ----------------
            partitions(dL) % mpiface_elements(bfaceID(dL)) = eL % eID
            partitions(dR) % mpiface_elements(bfaceID(dR)) = eR % eID
!
!           Get the face sides in the elements (current partition)
!           ------------------------------------------------------
            partitions(dL) % element_mpifaceSide(bfaceID(dL)) = f % elementSide(1)
            partitions(dR) % element_mpifaceSide(bfaceID(dR)) = f % elementSide(2)
!
!           Get the face sides in the elements (neighbor partition)
!           ------------------------------------------------------
            partitions(dL) % element_mpifaceSideOther(bfaceID(dL)) = f % elementSide(2)
            partitions(dR) % element_mpifaceSideOther(bfaceID(dR)) = f % elementSide(1)
!
!           Get the face rotation
!           ---------------------
            partitions(dL) % mpiface_rotation(bfaceID(dL)) = f % rotation 
            partitions(dR) % mpiface_rotation(bfaceID(dR)) = f % rotation 
!
!           Get the element face side
!           -------------------------
            partitions(dL) %  mpiface_elementSide(bfaceID(dL)) = 1
            partitions(dR) %  mpiface_elementSide(bfaceID(dR)) = 2
!
!           Get the shared domain
!           ---------------------
            partitions(dL) % mpiface_sharedDomain(bfaceID(dL)) = dR
            partitions(dR) % mpiface_sharedDomain(bfaceID(dR)) = dL

            end associate
            end associate
         end do

      end subroutine GetPartitionBoundaryFaces
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     --------------------------------
!     Space-filling curve partitioning
!     --------------------------------
      subroutine GetSFCElementsPartition(no_of_domains, no_of_allElements, elementsDomain, useWeights, Nx, Ny, Nz)
         implicit none
         !-arguments--------------------------------------------------
         integer, intent(in)    :: no_of_domains
         integer, intent(in)    :: no_of_allElements
         integer, intent(inout) :: elementsDomain(no_of_allElements)
         logical, intent(in)    :: useWeights
         integer, intent(in)    :: Nx(no_of_allElements), Ny(no_of_allElements), Nz(no_of_allElements)
         !-local-variables--------------------------------------------
         integer :: elems_per_domain(no_of_domains)
         integer :: biggerdomains
         integer :: first, last, domain
         integer :: ielem, ndof, max_dof, dof_in_domain
         integer :: dof_per_domain(no_of_domains), start_index(no_of_domains+1)
         logical                :: neddWeights
         integer, allocatable, target   :: weights(:)
         !------------------------------------------------------------

         if (useWeights) then
             allocate(weights(no_of_allElements))
             do ielem=1,no_of_allElements
                 weights(ielem) = (Nx(ielem) + 1) * (Ny(ielem) + 1) * (Nz(ielem) + 1)
             end do
             if (maxval(weights) .eq. minval(weights)) then
                 neddWeights = .false.
                 deallocate(weights)
             else
                 neddWeights = .true.
                 ndof = sum(weights)
             endif
         end if 
         
         elems_per_domain = no_of_allElements / no_of_domains
         biggerdomains = mod(no_of_allElements,no_of_domains)
         elems_per_domain(1:biggerdomains) = elems_per_domain(1:biggerdomains) + 1
         
         first = 1
         do domain = 1, no_of_domains
            last = first + elems_per_domain(domain) - 1
            elementsDomain(first:last) = domain
            first = last + 1
         end do

         if (neddWeights) then

             max_dof = ndof / no_of_domains
             start_index = 1
             do domain = 1, no_of_domains
                 start_index(domain+1) = start_index(domain) + elems_per_domain(domain)
             end do

             do domain = 1, no_of_domains-1
                 if (start_index(domain) .ge. start_index(domain+1)) start_index(domain+1) = start_index(domain) + 1
                 dof_in_domain = sum(weights(start_index(domain):start_index(domain+1)))
                 do ielem=1,no_of_allElements
                     if (dof_in_domain .lt. max_dof) then
                         start_index(domain+1) = start_index(domain+1) + 1
                         dof_in_domain = sum(weights(start_index(domain):start_index(domain+1)))
                         if (abs(dof_in_domain-max_dof) .le. abs(dof_in_domain-max_dof+weights(start_index(domain+1)+1))) exit
                     else
                         start_index(domain+1) = start_index(domain+1) - 1
                         dof_in_domain = sum(weights(start_index(domain):start_index(domain+1)))
                         if (abs(dof_in_domain-max_dof) .le. abs(dof_in_domain-max_dof-weights(start_index(domain+1)-1))) exit
                     end if
                 end do
             end do

             dof_per_domain = 0
             do domain = 1, no_of_domains
                 dof_per_domain(domain) = sum(weights(start_index(domain):start_index(domain+1)-1))
             end do

             do domain = 1, no_of_domains
                elementsDomain(start_index(domain):start_index(domain+1)-1) = domain
             end do

         end if
         
      end subroutine GetSFCElementsPartition
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     ---------------------------------
!     Write the partitions' information
!     ---------------------------------
      subroutine WritePartitionsFile(mesh,elementsDomain)
         implicit none
         !-arguments--------------------------------------------------
         type(HexMesh), intent(in)  :: mesh
         integer      , intent(in)  :: elementsDomain(mesh % no_of_elements)
         !-local-variables--------------------------------------------
         character(LINE_LENGTH)     :: pmeshName
         integer                    :: fID, eID
         !------------------------------------------------------------
         
         pmeshName = "./MESH/" // trim(removePath(getFileName(mesh % meshFileName))) // ".pmesh"
         
         open(newunit = fID, file=trim(pmeshName),action='write')
         
         write(fID,*) mesh % no_of_elements
         do eID = 1, mesh % no_of_elements
            write(fID,*) elementsDomain(eID)
         end do
            
         close(fID)
         
      end subroutine WritePartitionsFile

      ! This routine writes one file per partition (in the MESH folder) containing all the
      ! mesh partitioning information. Currently, this is done in ASCII format.
      ! This could be improved by using the hdf5 library and using a single binary file.
      subroutine WriteMeshPartitioningFile(mesh, no_of_domains, partitions)
         implicit none
         !-arguments--------------------------------------------------
         type(HexMesh), intent(in)  :: mesh
         integer,       intent(in)  :: no_of_domains
         type(PartitionedMesh_t), intent(in) :: partitions(no_of_domains)
         !-local-variables--------------------------------------------
         character(LINE_LENGTH)     :: pmeshName
         character(LINE_LENGTH)     :: pID_ch
         integer                    :: fileID, eID, pID, nID, fID
         logical                    :: meshIsHOPR
         !------------------------------------------------------------
         
         meshIsHOPR = allocated(mesh % HOPRnodeIDs)

         do pID = 1, no_of_domains
            write(pID_ch, '(I0)') pID

            pmeshName = "./MESH/" // trim(removePath(getFileName(mesh % meshFileName))) // "_" // trim(pID_ch) // ".partitioning"
            open(newunit = fileID, file=trim(pmeshName),action='write')
            
            write(fileID,*) no_of_domains, pID

            write(fileID,*) partitions % no_of_nodes, partitions % no_of_elements, partitions % no_of_mpifaces
            ! Write nodeIDs
            do nID = 1, no_of_nodes
               write(fileID,*) partitions(pID) % nodeIDs(nID)
            end do
            ! Write elementIDs
            do eID = 1, no_of_elements
               write(fileID,*) partitions(pID) % elementIDs(eID)
            end do
            ! Write mpiface_elements
            do fID = 1, no_of_mpifaces
               write(fileID,*) partitions(pID) % mpiface_elements(fID)
            end do
            ! Write element_mpifaceSide
            do fID = 1, no_of_mpifaces
               write(fileID,*) partitions(pID) % element_mpifaceSide(fID)
            end do

            ! Write element_mpifaceSideOther
            do fID = 1, no_of_mpifaces
               write(fileID,*) partitions(pID) % element_mpifaceSideOther(fID)
            end do

            ! Write mpiface_rotation
            do fID = 1, no_of_mpifaces
               write(fileID,*) partitions(pID) % mpiface_rotation(fID)
            end do

            ! Write mpiface_elementSide
            do fID = 1, no_of_mpifaces
               write(fileID,*) partitions(pID) % mpiface_elementSide(fID)
            end do

            ! Write mpiface_sharedDomain
            do fID = 1, no_of_mpifaces
               write(fileID,*) partitions(pID) % mpiface_sharedDomain(fID)
            end do
            
            write(fileID,*) "meshIsHOPR", meshIsHOPR

            if (meshIsHOPR) then
               ! Write HOPRnodeIDs
               do nID = 1, no_of_nodes
                  write(fileID,*) partitions(pID) % HOPRnodeIDs(nID)
               end do
            end if

            close(fileID)
         end do
            
      end subroutine WriteMeshPartitioningFile
         
end module MeshPartitioning
