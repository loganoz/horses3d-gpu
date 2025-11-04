#include "Includes.h"
Module ChannelForcing  !
#if defined(NAVIERSTOKES)
    use SMConstants
    use MPI_Process_Info
    use PhysicsStorage
#ifdef _HAS_MPI_
    use mpi
#endif
    Implicit None

private
public initializeChannel, updateChannel, channelSource

!
!     Module variables 
!     ----------------
!
      real(kind=RP) :: flowDir(NDIM)            ! Flow direction
      real(kind=RP) :: f_x(2), f_y(2), f_z(2)   ! Correction factors in stream-wise and span-wise direction
      real(kind=RP) :: dpdx, dpdy, dpdz         ! Pressure gradient
      
      real(kind=RP), allocatable :: dp_part(:)  ! dpdx, dpdy or dpdz in the different partitions (MPI)

      logical                    :: channelIsActive

!     Paramaters
!     ----------
      real(kind=RP), parameter :: alpha    =  0.9_RP


!!$acc declare create(f_x, f_y, f_z, dpdx, dpdy, dpdz, dp_part, flowDir, alpha)
!$acc declare create(f_x, f_y, f_z, dpdx, dpdy, dpdz)

!  ========
contains
!  ========
!
!///////////////////////////////////////////////////////////////////////////////////////
!
     Subroutine initializeChannel(controlVariables)
       use FTValueDictionaryClass
       use FluidData, only: refValues
       use Headers
       implicit none
       type(FTValueDictionary), intent(in)                     :: controlVariables

       integer :: ierr

       channelIsActive = controlVariables % logicalValueForKey("use channel")

       if (.not. channelIsActive) return
!
!      Initialize the flow direction
!      -----------------------------
       flowDir(1) = cos(refValues % AoATheta*pi/180._RP) * cos(refValues % AoAPhi*pi/180._RP)
       flowDir(2) = sin(refValues % AoAtheta*pi/180._RP) * cos(refValues % AoAphi*pi/180._RP)
       flowDir(3) = sin(refValues % AoAphi*pi/180._RP)
!
!      Initialize dpdx and dpdy
!      ------------------------
       dpdx  = 3.08e-3_RP
       ! dpdy  = 3.08e-3_RP
       ! dpdz  = 0._RP 
       dpdy  = 0._RP 
       ! dpdz  = 3.08e-3_RP
       dpdz  = 3.08e-5_RP
!      ----------------------
       f_x(1) = 0.0_RP
       f_x(2) = (1.0_RP - alpha) * dpdx
       f_y    = 0.0_RP
       f_z    = 0.0_RP
!
!      Initialize MPI
!      --------------
       if (MPI_Process % doMPIAction) then
           allocate( dp_part(MPI_Process % nProcs) )
       end if

       channelIsActive = .true.

!!$acc update device(alpha)
!!$acc update device(flowDir)
!$acc update device(dpdx)
!$acc update device(dpdy)
!$acc update device(dpdz)
!$acc update device(f_x)
!$acc update device(f_y)
!$acc update device(f_z)
       if ( .not. MPI_Process % isRoot ) return
       call Subsection_Header("Channel")
       write(STD_OUT,'(30X,A,A28)') "->", "Channel activated"

     End Subroutine initializeChannel
!
!//////////////////////////////////////////////////////////////////////// 
! 
     Subroutine updateChannel(mesh, dt, Monitors)
!          * Reading mean momentum/velocity in the volume
!          * Obtaining fTurbulentChannel according to it
!       -------------------------------------------------------------------------
!
        use HexMeshClass
        use MonitorsClass
        implicit none

        type(HexMesh)               :: mesh
        real(kind=RP)                :: dt
        type(Monitor_t), intent(in) :: monitors

        real(kind=RP) :: meanMomentum(NDIM)   ! Mean momentum
        integer       :: ierr, dp_ind
        !-------------------------------------------------------------

        if (.not. channelIsActive) return

!$acc update self(dpdx)
!$acc update self(dpdy)
!$acc update self(dpdz)
!
!       Get the same global dpdx,dpdy,dpdz in all partitions
!       ----------------------------------------------------
        if (MPI_Process % doMPIAction) then
#ifdef _HAS_MPI_
!          dpdx
!          ----
           call MPI_Allgather( dpdx, 1, MPI_DOUBLE, dp_part, 1, MPI_DOUBLE, MPI_COMM_WORLD, ierr)
           dp_ind  = maxloc(abs(dp_part),1)
           dpdx    = dp_part(dp_ind)
!
!          dpdy
!          ----
           call MPI_Allgather( dpdy, 1, MPI_DOUBLE, dp_part, 1, MPI_DOUBLE, MPI_COMM_WORLD, ierr)
           dp_ind  = maxloc(abs(dp_part),1)
           dpdy    = dp_part(dp_ind)
!
!          dpdz
!          ----
           call MPI_Allgather( dpdz, 1, MPI_DOUBLE, dp_part, 1, MPI_DOUBLE, MPI_COMM_WORLD, ierr)
           dp_ind  = maxloc(abs(dp_part),1)
           dpdz    = dp_part(dp_ind)
#endif
!$acc update device(dpdx)
!$acc update device(dpdy)
!$acc update device(dpdz)
!
        end if
!
!       Compute the mean velocity/momentum
!       ----------------------------------
        meanMomentum = monitors % volumeMonitors(1) % getLast()
!
!       Compute the correction parameter
!       --------------------------------

        f_x(1) = ( flowDir(1) - meanMomentum(1) ) * alpha / dt
        f_x(2) = (1._RP - alpha) * dpdx
        
        f_y(1) = ( flowDir(2) - meanMomentum(2) ) * alpha / dt
        f_y(2) = (1._RP - alpha) * dpdy

        f_z(1) = ( flowDir(3) - meanMomentum(3) ) * alpha / dt
        f_z(2) = (1._RP - alpha) * dpdz
!
!$acc update device(f_x)
!$acc update device(f_y)
!$acc update device(f_z)
     End Subroutine updateChannel
!
!//////////////////////////////////////////////////////////////////////// 
! 
     Subroutine channelSource(mesh)
         USE HexMeshClass
         Implicit None
         type(HexMesh), intent(in)         :: mesh

         real(kind=RP) :: dp_vals(2), rho
         integer       :: dp_ind, i, j, k, eID
         real(kind=RP) :: Su,Sv,Sw,Se

        if (.not. channelIsActive) return

!$acc parallel loop gang vector_length(128) present(mesh) async(1) private(dp_vals)
         do eID = 1, mesh % no_of_elements
            !$acc loop vector collapse(3) reduction(max:dpdx,dpdy,dpdz)
            do k = 0, mesh % elements(eID) % Nxyz(3) ; do j = 0, mesh % elements(eID) % Nxyz(2) ; do i = 0, mesh % elements(eID) % Nxyz(1)

                rho = mesh % elements(eID) % storage % Q(IRHO,i,j,k)
                Su = f_x(1) * rho + f_x(2)
                Sv = 0.0_RP
                Sw = f_z(1) * rho + f_z(2)
                Se = ( mesh % elements(eID) % storage % Q(IRHOU,i,j,k) / rho ) * Su + (mesh % elements(eID) % storage % Q(IRHOV,i,j,k) / rho ) * Sv + &
                     ( mesh % elements(eID) % storage % Q(IRHOW,i,j,k) / rho ) * Sw

                mesh % elements(eID) % storage % S_NS(IRHOU,i,j,k) = Su
                mesh % elements(eID) % storage % S_NS(IRHOV,i,j,k) = Sv
                mesh % elements(eID) % storage % S_NS(IRHOW,i,j,k) = Sw
                mesh % elements(eID) % storage % S_NS(IRHOE,i,j,k) = Se
                
                dp_vals = [dpdx,Su]
                dp_ind  = maxloc(abs(dp_vals),1)
                dpdx    = dp_vals(dp_ind)
                
                dp_vals = [dpdy,Sv]
                dp_ind  = maxloc(abs(dp_vals),1)
                dpdy    = dp_vals(dp_ind)
                
                dp_vals = [dpdz,Sw]
                dp_ind  = maxloc(abs(dp_vals),1)
                dpdz    = dp_vals(dp_ind)

            end do                  ; end do                ; end do
         end do
!$acc end parallel loop

     End Subroutine channelSource

#endif
End Module ChannelForcing
