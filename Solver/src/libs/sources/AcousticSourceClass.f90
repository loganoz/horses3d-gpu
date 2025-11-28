!//////////////////////////////////////////////////////
!
! This class represents the numerical source term that represents a Acoustic Source
! The acoustic source is based on the multipole expansion 

!Limitations:

#include "Includes.h"
module AcousticSourceClass 
    use SMConstants
    use HexMeshClass
    use SolutionFile
    use PhysicsStorage, only: NCONS
    Implicit None

    private
    public AcousticSource, addSourceAcoustic, ConstructAcousticSource, DestructAcousticSource 
    integer                    :: Nfreq                               ! Number of frequency YET EVERYTHING IS HARDCODED, SO ADJUST freq_Hz accordingly 
    integer, parameter         :: numberDim = NDIM                    ! 2D Source 
    integer, parameter         :: ratio = 24                          ! wL / sigma ratio    
    real(kind=RP), parameter   :: first_central_frequency  = 10.0_RP  ! First central frequency for the 1/3 diagram 
    real(kind=RP), parameter   :: R = sqrt(250**2+87.6**2)            ! Distance at which the WT SPL is received 
    real(kind=RP), parameter   :: fMAX = 150.0_RP                     ! Maximum frequency (Hz) to consider (ideally the maximum that the mesh is able to capture)
    real(kind=RP), dimension(34) :: SPL_WT = [ &
                3.983823E+001_RP, 4.175864E+001_RP, 4.397531E+001_RP, 4.629396E+001_RP, &
                4.893568E+001_RP, 5.178563E+001_RP, 5.466095E+001_RP, 5.719623E+001_RP, &
                5.914524E+001_RP, 5.969899E+001_RP, 5.892264E+001_RP, 5.752926E+001_RP, &
                5.550930E+001_RP, 5.305830E+001_RP, 5.036188E+001_RP, 4.780741E+001_RP, &
                4.562016E+001_RP, 4.397583E+001_RP, 4.254709E+001_RP, 4.117360E+001_RP, &
                3.981533E+001_RP, 3.830182E+001_RP, 3.647226E+001_RP, 3.479272E+001_RP, &
                3.305135E+001_RP, 3.113701E+001_RP, 2.900127E+001_RP, 2.684142E+001_RP, &
                2.442864E+001_RP, 2.174479E+001_RP, 1.905975E+001_RP, 1.619207E+001_RP, &
                1.278492E+001_RP, 9.467272E+000_RP ]                                          ! Data from WT 

    !definition of Monopole Source
    type AcousticSource_t
        real(kind=RP), allocatable                               :: Amplitude(:)              ! amplitude of the source term
        real(kind=RP), allocatable                               :: freq_Hz(:)                ! Frequency of the source (Hz)
        real(kind=RP), allocatable                               :: freq(:)                   ! Frequency of the source (-)  
        real(kind=RP), allocatable                               :: phase(:)                  ! Phase of the wave (rad)
        real(kind=RP), allocatable                               :: omega(:)                  ! Angular frequency (-) 
        real(kind=RP), allocatable                               :: sigma(:)                  ! Gaussian standard deviation 
        real(kind=RP)                                            :: x0(NDIM)                  ! Source position 
        real(kind=RP)                                            :: ratio                     ! wL / sigma ratio 
        real(kind=RP)                                            :: c0                        ! Sound speed (-)
        logical                                                  :: isActive = .false.     

        contains

        procedure :: CreateDeviceData  => sourceCreateDeviceData
        procedure :: ExitDeviceData    => sourceExitDeviceData

    end type AcousticSource_t

    type(AcousticSource_t)                                       :: AcousticSource 

    contains 
        function hankel_h0_1(x) result(H0)
                real(kind=RP), intent(in) :: x
                complex(kind=RP) :: H0
        
                H0 = cmplx(bessel_j0(x), bessel_y0(x), kind=RP)
        end function hankel_h0_1

!/////////////////////////////////////////////////////////////////////////
!           ACOUSTIC SOURCE PROCEDURES --------------------------
!/////////////////////////////////////////////////////////////////////////
    subroutine sourceCreateDeviceData(self)
#ifdef _OPENACC
        use cudafor
        use openacc
#endif
        implicit none
        !-----------------------------------------------------------
        class(AcousticSource_t)       :: self
    
#ifdef _OPENACC
        !$acc enter data copyin(self)
        !$acc enter data copyin(self%x0)
        !$acc enter data copyin(self%c0)
        !$acc enter data copyin(self%Amplitude)
        !$acc enter data copyin(self%sigma)
        !$acc enter data copyin(self%omega)
        !$acc enter data copyin(self%phase)
#endif
        !
    end subroutine sourceCreateDeviceData
        !
    subroutine sourceExitDeviceData(self)
#ifdef _OPENACC
        use cudafor
        use openacc
#endif
        implicit none
        !-----------------------------------------------------------
        class(AcousticSource_t)               :: self
        !-----------------------------------------------------------
        
#ifdef _OPENACC
        !$acc exit data delete(self)
        !$acc exit data delete(self%x0)
        !$acc exit data delete(self%c0)
        !$acc exit data delete(self%Amplitude)
        !$acc exit data delete(self%sigma)
        !$acc exit data delete(self%omega)
        !$acc exit data delete(self%phase)
#endif
        !
    end subroutine sourceExitDeviceData

    Subroutine ConstructAcousticSource(self, controlVariables)
        use FluidData
        use FileReadingUtilities, only: getRealArrayFromString
        use FTValueDictionaryClass
        use MPI_Process_Info
        use Headers
        use mainKeywordsModule
        use FileReadingUtilities, only: getFileName
        
        Implicit None
        class(AcousticSource_t), intent(inout)                :: self
        type(FTValueDictionary), intent(in)                   :: controlVariables
        
        ! Local variables 
        real(kind=RP)                                         :: xwrap(Nfreq), dummy(Nfreq)
        real(kind=RP)                                         :: p_ref, rho_ref, U_ref
        real(kind=RP)                                         :: k, prms_band
        integer                                               :: i  
#ifdef FLOW
        ! Acoustic reference pressure value 
        p_ref = 20E-06_RP ! Acoustic pressure reference value 
        rho_ref = refValues%rho ! Reference density in Horses 
        U_ref = refValues%V ! Reference velocity in Horses 

        ! To begin, everything hardcoded 
        ! Just the active flag is read from the control file
        self % isActive = .false.
        if (.not. controlVariables % logicalValueForKey("use monopole")) return ! What happens if there is no use monopole in the control file 
        self % isActive = .true. 


#ifdef NAVIERSTOKES
        self % c0 = 1.0 / dimensionless % Mach     ! Dimensionless speed of sound 
#endif
#ifdef MULTIPHASE 
        self % c0 = 1 / sqrt(dimensionless%Ma2(1)) ! Dimensionless speed of sound of fluid one 
#endif 

        safedeallocate(self % Amplitude)
        safedeallocate(self % freq_Hz)
        safedeallocate(self % freq)
        safedeallocate(self % phase)
        safedeallocate(self % omega)
        safedeallocate(self % sigma)

        ! Number of frequencies depending on fMAX 
        Nfreq = 1+int(3.0_RP*log(fMAX*U_ref/first_central_frequency)/log(2.0_RP))

        allocate(self%Amplitude(Nfreq))
        allocate(self%freq_Hz(Nfreq))
        allocate(self%freq(Nfreq))
        allocate(self%phase(Nfreq))
        allocate(self%omega(Nfreq))
        allocate(self%sigma(Nfreq))

        ! Source position 
        self % x0(1) = -15.0_RP 
		self % x0(2) = 120.0_RP 
		self % x0(3) = 0.0_RP 

        ! Amplitude for each frequency 
        do i=1,Nfreq
                self % freq_Hz(i) = first_central_frequency*2.0_RP**(real(i-1,RP)/3.0_RP)
                self % freq(i) = self % freq_Hz(i) / U_ref 
                self % omega(i) = 2 * PI * self % freq(i)
                self % sigma(i) = self % c0 / (self % freq(i) * ratio)                         
                k = self % omega(i) / self % c0 
                prms_band = 10**(SPL_WT(i)/20.0_RP) * p_ref / (rho_ref * U_ref * U_ref)
                self % Amplitude(i) = prms_band*4*self%c0*sqrt(2.0_RP)/k/abs(hankel_h0_1(k*R)) ! Monopole Theoretical
        enddo 


        if (Nfreq .eq. 1) then
            self % phase = 0.0_RP
        else
        ! phase using parabolic distribution to avoid over increases
            dummy = 1.0
            self % phase = [(i,i=1,Nfreq)]
            self % phase = 1 - self % phase * self % phase
            self % phase = -PI/(real(Nfreq,RP))*self % phase
            xwrap = mod(self % phase, 2.0_RP*PI)
            self % phase = xwrap + merge(-2.0_RP*PI,0.0_RP,abs(xwrap)>PI)*sign(dummy, xwrap)
        end if
#endif 
        ! We send the variables to the GPU 
#ifdef _OPENACC
        call self % CreateDeviceData()
#endif 
        ! TODO, define some user parametrization from control file 

        if ( .not. MPI_Process % isRoot ) return
        call Subsection_Header("Monopole Acoustic Source")
        write(STD_OUT,'(30X,A,A28,I0)')    "->", "Source Dim: ",              numberDim
        write(STD_OUT,'(30X,A,A28,*(F10.2))') "->", "Amplitude: ",            self % Amplitude
        write(STD_OUT,'(30X,A,A28,*(F10.2))') "->", "Phase (rad): ",          self % phase
        write(STD_OUT,'(30X,A,A28,*(F10.2))') "->", "Frequency (Hz): ",       self % freq_Hz
        write(STD_OUT,'(30X,A,A28,*(F10.2))') "->", "Frequency (-): ",        self % freq
        write(STD_OUT,'(30X,A,A28,*(F10.2))') "->", "Gaussian width (-): ",   self % sigma
        write(STD_OUT,'(30X,A,A28,*(F10.2))') "->", "Ratio: ",                self % ratio

    end subroutine ConstructAcousticSource

    Subroutine DestructAcousticSource(self)
        Implicit None
        class(AcousticSource_t), intent(inout)             :: self

!       Check if is activated
!       ------------------------
        if (.not. self % isActive) return

!       Deallocate all the arrays in CPU and GPU
!       ------------------------
#ifdef _OPENACC
        call self % ExitDeviceData()
#endif

        safedeallocate(self % Amplitude)
        safedeallocate(self % freq_Hz)
        safedeallocate(self % freq)
        safedeallocate(self % phase)
        safedeallocate(self % omega)
        safedeallocate(self % sigma)

    End Subroutine DestructAcousticSource
!   

    Subroutine addSourceAcoustic(self,mesh, time)
        use FluidData 
        Implicit None
        class(AcousticSource_t), intent(in)              :: self
        type(HexMesh), intent(inout)                     :: mesh
        real(kind=RP), intent(in)                        :: time 
    
        !local variables
        integer                                          :: i, j, k, eID, Nxyz(NDIM), dim, f_index
        real(kind=RP)                                    :: Sp, r2
        real(kind=RP)                                    :: s2, factor, GM 

!       Check if is activated
!       ------------------------
        if (.not. self % isActive) return

        !$acc parallel loop gang vector_length(128) present(self,mesh, mesh % elements) private(eID,i,j,k,dim,f_index, r2,s2,factor,GM,Sp) async(1)
        do eID = 1, mesh % no_of_elements 
            !$acc loop vector collapse(3) 
            do k = 0, mesh % elements(eID) % Nxyz(3) 
                do j = 0, mesh % elements(eID) % Nxyz(2) 
                    do i = 0, mesh % elements(eID) % Nxyz(1)
                        ! Distance to the source 
                        r2 = 0.0_RP 
                        !$acc loop seq
                        do dim = 1, numberDim
                                r2 = r2 + POW2(mesh % elements(eID) % geom % x(dim,i,j,k)-self%x0(dim))
                        end do 
                        ! Acoustic source corresponding to each frequency 
                        Sp = 0.0_RP 
                        !$acc loop seq
                        do f_index = 1, Nfreq
                                s2 = POW2(self%sigma(f_index)) 
                                factor = 1.0_RP / sqrt((2.0_RP*PI*s2)**numberDim) ! Normalize gaussian 
                                GM = exp(-r2/(2.0_RP*s2)) ! Gaussian distribution of x0 mu and sigma deviation
                                Sp = Sp + self%Amplitude(f_index) * GM * cos(self%omega(f_index)*time+self%phase(f_index)) * factor
                        end do 
#if defined(NAVIERSTOKES) 
                ! In NS source terms in density and energy equations 
                mesh % elements(eID) % storage % S_NS(1,i,j,k) =  mesh % elements(eID) % storage % S_NS(1,i,j,k) + Sp / POW2(self % c0)
                mesh % elements(eID) % storage % S_NS(5,i,j,k) =  mesh % elements(eID) % storage % S_NS(5,i,j,k) + Sp / (thermodynamics % gamma - 1.0_RP)
#endif
#if defined(ACOUSTIC)
                mesh % elements(eID) % storage % S_NS(5,i,j,k) = mesh % elements(eID) % storage % S_NS(5,i,j,k) + Sp
#endif
                    
#if defined(MULTIPHASE)
                mesh % elements(eID) % storage % S_NS(5,i,j,k) = mesh % elements(eID) % storage % S_NS(5,i,j,k) + Sp 
                ! Here the 1 /c0^2 is included intrinsically due to the equations 
#endif 
        end do         ; end do          ; end do       ; end do 
        !$acc end parallel loop  
    end subroutine addSourceAcoustic
end module AcousticSourceClass