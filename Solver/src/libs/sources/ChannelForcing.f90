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

!     Parameters
!     ----------
      real(kind=RP), parameter :: alpha    =  0.9_RP


!$acc declare create(f_x, f_y, f_z)

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
         type(HexMesh), intent(inout)         :: mesh

         real(kind=RP) :: rho, rho_max, rho_min, S_max, S_min
         real(kind=RP) :: local_rho_max, local_rho_min
         integer       :: i, j, k, eID, dp_ind
         real(kind=RP) :: Su,Sv,Sw,Se
         real(kind=RP) :: dp_vals(3)

        if (.not. channelIsActive) return
!
        rho_min = 2.0_RP
        rho_max = 0.0_RP

!$acc parallel loop gang present(mesh) num_gangs(9700) reduction(max:rho_max) reduction(min:rho_min)
         do eID = 1, mesh % no_of_elements
            local_rho_min = 2.0_RP
            local_rho_max = 0.0_RP
            !$acc loop vector collapse(3) reduction(max:local_rho_max) reduction(min:local_rho_min)
            do k = 0, mesh % elements(eID) % Nxyz(3) ; do j = 0, mesh % elements(eID) % Nxyz(2) ; do i = 0, mesh % elements(eID) % Nxyz(1)

                rho = mesh % elements(eID) % storage % Q(IRHO,i,j,k)
                Su = f_x(1) * rho + f_x(2)
                Sv = 0.0_RP
                Sw = f_z(1) * rho + f_z(2)
                Se = ( mesh % elements(eID) % storage % Q(IRHOU,i,j,k) * Su + mesh % elements(eID) % storage % Q(IRHOV,i,j,k)* Sv + &
                       mesh % elements(eID) % storage % Q(IRHOW,i,j,k) * Sw ) / rho 

                mesh % elements(eID) % storage % S_NS(IRHOU,i,j,k) = Su
                mesh % elements(eID) % storage % S_NS(IRHOV,i,j,k) = Sv
                mesh % elements(eID) % storage % S_NS(IRHOW,i,j,k) = Sw
                mesh % elements(eID) % storage % S_NS(IRHOE,i,j,k) = Se

                local_rho_min = min(local_rho_min, rho)
                local_rho_max = max(local_rho_max, rho)
                
            end do                  ; end do                ; end do
            rho_min = min(rho_min, local_rho_min)
            rho_max = max(rho_max, local_rho_min)
         end do
!$acc end parallel loop
!
! get the maximum absolute value using the extreme values of rho (only variable changing in the loop)
        S_min   = f_x(1) * rho_min + f_x(2)
        S_max   = f_x(1) * rho_max + f_x(2)
        dp_vals = [dpdx,S_max,S_min]
        dp_ind  = maxloc(abs(dp_vals),1)
        dpdx    = dp_vals(dp_ind)

        ! S_min   = f_y(1) * rho_min + f_y(2)
        ! S_max   = f_y(1) * rho_max + f_y(2)
        S_min   = 0.0_RP
        S_max   = 0.0_RP
        dp_vals = [dpdy,S_max,S_min]
        dp_ind  = maxloc(abs(dp_vals),1)
        dpdy    = dp_vals(dp_ind)

        S_min   = f_z(1) * rho_min + f_z(2)
        S_max   = f_z(1) * rho_max + f_z(2)
        dp_vals = [dpdz,S_max,S_min]
        dp_ind  = maxloc(abs(dp_vals),1)
        dpdz    = dp_vals(dp_ind)
!
     End Subroutine channelSource

#endif
End Module ChannelForcing
