#include "Includes.h"
module MonitorsClass
   use SMConstants
   use NodalStorageClass
   use HexMeshClass
   use MonitorDefinitions
   use ResidualsMonitorClass
   use VolumeMonitorClass
   use LoadBalancingMonitorClass
   use FileReadingUtilities      , only: getFileName
#ifdef FLOW
   use ProbeClass
#endif
#if defined(NAVIERSTOKES) || defined(INCNS)
   use StatisticsMonitor
   use SurfaceMonitorClass
#endif
   implicit none
!

   private
   public      Monitor_t
!
!  *****************************
!  Main monitor class definition
!  *****************************
!  
   type Monitor_t
      character(len=LINE_LENGTH)                 :: solution_file
      character(len=LINE_LENGTH)                 :: probes_solution_file = ""
      integer                                    :: no_of_probes
      integer                                    :: no_of_surfaceMonitors
      integer                                    :: no_of_volumeMonitors
      integer                                    :: no_of_loadBalancingMonitors
      integer                                    :: no_of_fileProbes = 0
      character(len=LINE_LENGTH)                 :: probesFileName = ""
      character(len=STR_LEN_MONITORS), allocatable :: probesVariables(:)
      real(kind=RP)                              :: probeFileSaveTimestep = 0.0_RP
      character(len=8)                           :: probeFileOutputFormat = "ASCII"
      real(kind=RP)                              :: fp_lastSavedTime = -huge(0.0_RP)
      integer                                    :: bufferLine
      integer                      , allocatable :: iter(:)
      integer                                    :: dt_restriction
      logical                                    :: write_dt_restriction
      real(kind=RP)                , allocatable :: t(:)
      real(kind=RP)                , allocatable :: SolverSimuTime(:)
      real(kind=RP)                , allocatable :: TotalSimuTime(:)
      type(Residuals_t)                          :: residuals
      class(VolumeMonitor_t)       , allocatable :: volumeMonitors(:)
      class(LoadBalancingMonitor_t), allocatable :: loadBalancingMonitors(:)
#ifdef FLOW
      class(Probe_t)               , allocatable :: probes(:)
#endif
#if defined(NAVIERSTOKES) || defined(INCNS)
      class(SurfaceMonitor_t)      , allocatable :: surfaceMonitors(:)
      type(StatisticsMonitor_t)                  :: stats
#endif
      contains
         procedure   :: Construct       => Monitors_Construct
         procedure   :: WriteLabel      => Monitor_WriteLabel
         procedure   :: WriteUnderlines => Monitor_WriteUnderlines
         procedure   :: WriteValues     => Monitor_WriteValues
         procedure   :: UpdateValues    => Monitor_UpdateValues
         procedure   :: WriteToFile     => Monitor_WriteToFile
         procedure   :: WriteProbesFileSummary => Monitor_WriteProbesFileSummary
         procedure   :: destruct        => Monitor_Destruct
         procedure   :: copy            => Monitor_Assign
         generic     :: assignment(=)   => copy
   end type Monitor_t
!
!  ========
   contains
!  ========
!
!///////////////////////////////////////////////////////////////////////////////////////
!
      subroutine Monitors_Construct( Monitors, mesh, controlVariables )
         use FTValueDictionaryClass
         use mainKeywordsModule
         use MPI_Process_Info
         implicit none
         class(Monitor_t)                     :: Monitors
         class(HexMesh), intent(in)           :: mesh
         class(FTValueDictionary), intent(in) :: controlVariables
         
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                         :: fID , io
         integer                         :: i
         character(len=STR_LEN_MONITORS) :: line
         character(len=STR_LEN_MONITORS) :: solution_file
         logical, save                   :: FirstCall = .TRUE.
         logical                         :: saveGradients
         character(len=LINE_LENGTH)      :: probesFileName
         character(len=LINE_LENGTH)      :: probesVariablesLine
         character(len=STR_LEN_MONITORS), allocatable :: probesVariables(:)
         integer                         :: no_of_fileProbes
         integer                         :: no_of_probesVariables
         real(kind=RP)                   :: probeFileSaveTimestep
         character(len=8)                :: probeFileOutputFormat
         character(len=LINE_LENGTH)      :: probes_solution_file
         character(len=LINE_LENGTH)      :: probes_dir
         integer                         :: last_slash
!
!        Setup the buffer
!        ----------------
         if (controlVariables % containsKey("monitors flush interval") ) then
            BUFFER_SIZE = controlVariables % integerValueForKey("monitors flush interval")
         end if
         
         allocate ( Monitors % TotalSimuTime(BUFFER_SIZE), &
                    Monitors % SolverSimuTime(BUFFER_SIZE), &
                    Monitors % t(BUFFER_SIZE), &
                    Monitors % iter(BUFFER_SIZE) )
!
!        Get the solution file name
!        --------------------------
         solution_file = controlVariables % stringValueForKey( solutionFileNameKey, requestedLength = STR_LEN_MONITORS )
!
!        Remove the *.hsol termination
!        -----------------------------
         solution_file = trim(getFileName(solution_file))
         Monitors % solution_file = trim(solution_file)
!
!        Build the probes output subdirectory: <dir>/probes/<base>
!        ----------------------------------------------------------
         last_slash = index(trim(solution_file), '/', back=.true.)
         if (last_slash .gt. 0) then
            probes_dir           = trim(solution_file(1:last_slash)) // "probes"
            probes_solution_file = trim(solution_file(1:last_slash)) // "probes/" // trim(solution_file(last_slash+1:))
         else
            probes_dir           = "probes"
            probes_solution_file = "probes/" // trim(solution_file)
         end if
         if (MPI_Process % isRoot) then
            call execute_command_line('mkdir -p "' // trim(probes_dir) // '"', wait=.true.)
         end if
         Monitors % probes_solution_file = trim(probes_solution_file)
!
!        Search in case file for probes, surface monitors, and volume monitors
!        ---------------------------------------------------------------------
         no_of_fileProbes = 0
         if (mesh % child) then ! Return doing nothing if this is a child mesh
            Monitors % no_of_probes = 0
            Monitors % no_of_surfaceMonitors = 0
            Monitors % no_of_volumeMonitors = 0
            Monitors % no_of_loadBalancingMonitors = 0
         else
            call getNoOfMonitors( Monitors % no_of_probes, Monitors % no_of_surfaceMonitors, Monitors % no_of_volumeMonitors, Monitors % no_of_loadBalancingMonitors )
!
!           Check for an additional probes definition file (one probe per
!           line: "x y z variable1 [variable2 ...]"), allowing several
!           variables to be sampled and saved for the same probe location
!           ---------------------------------------------------------------
#ifdef FLOW
            call readProbesFileBlock( probesFileName, probesVariablesLine, no_of_probesVariables, probeFileSaveTimestep, probeFileOutputFormat )

            if ( len_trim(probesFileName) .gt. 0 ) then
               call countProbesInFile( trim(probesFileName), no_of_fileProbes )

               if ( len_trim(probesVariablesLine) .gt. 0 ) then
                  call splitIntoTokens( trim(probesVariablesLine), probesVariables, no_of_probesVariables )
               else
                  write(STD_OUT,*) "Error: 'variables' must be specified inside the '#define probe file' block."
                  stop
               end if
            end if
#endif
         end if
!
!        Initialize the Monitors class in the GPU
!        ----------------------------------------
         !!$acc enter data copyin(Monitors)
         !$acc update device(Monitors)

!        Pro tip: This is necessary to avoid the compiler doing behind the scenes copies of the class.
!        Pro tip(cont): Its not really necessary for the class itself, but for the arrays inside the class.
!
!        Initialize
!        ----------
         call Monitors % residuals % Initialization( solution_file , FirstCall )

         allocate ( Monitors % volumeMonitors ( Monitors % no_of_volumeMonitors )  )
         !$acc update device(Monitors)
         do i = 1 , Monitors % no_of_volumeMonitors
            call Monitors % volumeMonitors(i) % Initialization ( mesh , i, solution_file , FirstCall  )
         end do

         allocate ( Monitors % loadBalancingMonitors ( Monitors % no_of_loadBalancingMonitors )  )
         !$acc update device(Monitors)
         do i = 1 , Monitors % no_of_loadBalancingMonitors
            call Monitors % loadBalancingMonitors(i) % Initialization ( mesh , i, solution_file , FirstCall )
         end do

#ifdef FLOW
         allocate ( Monitors % probes ( Monitors % no_of_probes + no_of_fileProbes )  )
         !$acc update device(Monitors)
         do i = 1 , Monitors % no_of_probes
            call Monitors % probes(i) % Initialization ( mesh , i, probes_solution_file , FirstCall )
         end do

         if ( no_of_fileProbes .gt. 0 ) then
            call InitializeProbesFromFile( trim(probesFileName), Monitors % probes, Monitors % no_of_probes, &
                                            mesh, probes_solution_file, FirstCall, probesVariables, probeFileSaveTimestep, &
                                            probeFileOutputFormat )
            Monitors % probesFileName         = trim(probesFileName)
            Monitors % probeFileSaveTimestep  = probeFileSaveTimestep
            Monitors % probeFileOutputFormat  = probeFileOutputFormat
            Monitors % fp_lastSavedTime       = -huge(0.0_RP)
            allocate( Monitors % probesVariables(size(probesVariables)) )
            Monitors % probesVariables = probesVariables
#ifdef HAS_HDF5
            if ( trim(probeFileOutputFormat) .eq. "HDF5" ) then
               call Monitor_InitFileProbesHDF5( Monitors, no_of_fileProbes )
            end if
#endif
         end if

         Monitors % no_of_fileProbes = no_of_fileProbes
         Monitors % no_of_probes = Monitors % no_of_probes + no_of_fileProbes
#endif

#if defined(NAVIERSTOKES) || defined(INCNS)
         saveGradients    = controlVariables % logicalValueForKey(saveGradientsToSolutionKey)
         call Monitors % stats     % Construct(mesh, saveGradients)


         allocate ( Monitors % surfaceMonitors ( Monitors % no_of_surfaceMonitors )  )
         !$acc update device(Monitors)
         do i = 1 , Monitors % no_of_surfaceMonitors
            call Monitors % surfaceMonitors(i) % Initialization ( mesh , i, solution_file , FirstCall )
         end do
#endif

         Monitors % write_dt_restriction = controlVariables % logicalValueForKey( "write dt restriction" )
         
         Monitors % bufferLine = 0
         
         FirstCall = .FALSE.
!
!        Include the latest changes in the GPU
!        ----------------------------------------
         !$acc update device(Monitors)

      end subroutine Monitors_Construct

      subroutine Monitor_WriteLabel ( self )
!
!        ***************************************************
!           This subroutine prints the labels for the time
!         integrator Display procedure.
!        ***************************************************
!
         use MPI_Process_Info
         implicit none
         class(Monitor_t)              :: self
         integer                       :: i 
      
         if ( .not. MPI_Process % isRoot ) return
!
!        Write "Iteration" and "Time"
!        ----------------------------
         write ( STD_OUT , ' ( A10    ) ' , advance = "no" ) "Iteration"
         write ( STD_OUT , ' ( 3X,A10 ) ' , advance = "no" ) "Time"
!
!        Write residuals labels
!        ----------------------
         call self % residuals % WriteLabel
!
!        Write volume monitors labels
!        -----------------------------
         do i = 1 , self % no_of_volumeMonitors
            call self % volumeMonitors(i) % WriteLabel
         end do
!
!        Write load balancing monitor labels
!        ------------------------------------
         do i = 1 , self % no_of_loadBalancingMonitors
            call self % loadBalancingMonitors(i) % WriteLabel
         end do

#ifdef FLOW
!
!        Write probes labels (file-based probes excluded for readability)
!        ---------------------------------------------------------------
         do i = 1 , self % no_of_probes - self % no_of_fileProbes
            call self % probes(i) % WriteLabel
         end do
#endif

#if defined(NAVIERSTOKES) || defined(INCNS)
!
!        Write surface monitors labels
!        -----------------------------
         do i = 1 , self % no_of_surfaceMonitors
            call self % surfaceMonitors(i) % WriteLabel
         end do

         call self % stats % WriteLabel
#endif
!
!        Write label for dt restriction
!        ------------------------------
         if (self % write_dt_restriction) write ( STD_OUT , ' ( 3X,A10 ) ' , advance = "no" ) "dt restr."

         write(STD_OUT , *) 

      end subroutine Monitor_WriteLabel

      subroutine Monitor_WriteUnderlines( self ) 
!
!        ********************************************************
!              This subroutine displays the underlines for the
!           time integrator Display procedure.
!        ********************************************************
!
         use PhysicsStorage
         use MPI_Process_Info
         implicit none
         class(Monitor_t)                         :: self
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                                  :: i, j
         character(len=MONITOR_LENGTH), parameter :: dashes = "----------"

         if ( .not. MPI_Process % isRoot ) return
!
!        Print dashes for "Iteration" and "Time"
!        ---------------------------------------
         write ( STD_OUT , ' ( A10    ) ' , advance = "no" ) trim ( dashes ) 
         write ( STD_OUT , ' ( 3X,A10 ) ' , advance = "no" ) trim ( dashes ) 
!
!        Print dashes for residuals
!        --------------------------
         do i = 1 , NCONS
            write(STD_OUT , '(3X,A10)' , advance = "no" ) trim(dashes)
         end do
!
!        Print dashes for volume monitors
!        --------------------------------
         do i = 1 , self % no_of_volumeMonitors  ; do j=1, size ( self % volumeMonitors(i) % values, 1 )
            write(STD_OUT , '(3X,A10)' , advance = "no" ) dashes(1 : min(10 , len_trim( self % volumeMonitors(i) % monitorName ) + 2 ) )
         end do                                  ; end do
!
!        Print dashes for load balancing monitor
!        --------------------------------------
         do i = 1 , self % no_of_loadBalancingMonitors ; do j=1, size ( self % loadBalancingMonitors(i) % values, 1 )
            write(STD_OUT , '(3X,A10)' , advance = "no" ) dashes(1 : min(10 , len_trim( self % loadBalancingMonitors(i) % monitorName ) + 2 ) )
         end do                                  ; end do

#ifdef FLOW
!
!        Print dashes for probes (file-based probes excluded for readability)
!        ----------------------------------------------------------------
         do i = 1 , self % no_of_probes - self % no_of_fileProbes
            if ( self % probes(i) % active ) then
               write(STD_OUT , '(3X,A10)' , advance = "no" ) dashes(1 : min(10 , len_trim( self % probes(i) % monitorName ) + 2 ) )
            end if
         end do
#endif

#if defined(NAVIERSTOKES) || defined(INCNS)
!
!        Print dashes for surface monitors
!        ---------------------------------
         do i = 1 , self % no_of_surfaceMonitors
            write(STD_OUT , '(3X,A10)' , advance = "no" ) dashes(1 : min(10 , len_trim( self % surfaceMonitors(i) % monitorName ) + 2 ) )
         end do

         if ( self % stats % state .ne. 0 ) write(STD_OUT,'(3X,A10)',advance="no") trim(dashes)
#endif

         
!
!        Print dashes for dt restriction
!        -------------------------------
         if (self % write_dt_restriction) write ( STD_OUT , ' ( 3X,A10 ) ' , advance = "no" ) trim ( dashes ) 
         
         write(STD_OUT , *) 

      end subroutine Monitor_WriteUnderlines

      subroutine Monitor_WriteValues ( self )
!
!        *******************************************************
!              This subroutine prints the values for the time
!           integrator Display procedure.
!        *******************************************************
!
         use MPI_Process_Info
         implicit none
         class(Monitor_t)           :: self
         integer                    :: i

         if ( .not. MPI_Process % isRoot ) return
!
!        Print iteration and time
!        ------------------------
         write ( STD_OUT , ' ( I10            ) ' , advance = "no" ) self % iter    ( self % bufferLine ) 
         write ( STD_OUT , ' ( 1X,A,1X,ES10.3 ) ' , advance = "no" ) "|" , self % t ( self % bufferLine ) 
!
!        Print residuals
!        ---------------
         call self % residuals % WriteValues( self % bufferLine )
!
!        Print volume monitors
!        ---------------------
         do i = 1 , self % no_of_volumeMonitors
            call self % volumeMonitors(i) % WriteValues ( self % bufferLine )
         end do
!
!        Print load balancing monitors
!        -----------------------------
         do i = 1 , self % no_of_loadBalancingMonitors
            call self % loadBalancingMonitors(i) % WriteValues ( self % bufferLine )
         end do

#ifdef FLOW
!
!        Print probes (file-based probes excluded for readability)
!        ----------------------------------------------------------
         do i = 1 , self % no_of_probes - self % no_of_fileProbes
            call self % probes(i) % WriteValues ( self % bufferLine )
         end do
#endif

#if defined(NAVIERSTOKES) || defined(INCNS)
!
!        Print surface monitors
!        ----------------------
         do i = 1 , self % no_of_surfaceMonitors
            call self % surfaceMonitors(i) % WriteValues ( self % bufferLine )
         end do

         call self % stats % WriteValue
#endif
!
!        Print dt restriction
!        --------------------
         if (self % write_dt_restriction) then
            select case (self % dt_restriction)
               case (DT_FIXED) ; write ( STD_OUT , ' ( 1X,A,1X,A10) ' , advance = "no" ) "|" , 'Fixed'
               case (DT_DIFF)  ; write ( STD_OUT , ' ( 1X,A,1X,A10) ' , advance = "no" ) "|" , 'Diffusive'
               case (DT_CONV)  ; write ( STD_OUT , ' ( 1X,A,1X,A10) ' , advance = "no" ) "|" , 'Convective'
            end select
         end if

         write(STD_OUT , *) 

      end subroutine Monitor_WriteValues

      subroutine Monitor_UpdateValues ( self, mesh, t , iter, maxResiduals, Autosave, dt )
!
!        ***************************************************************
!              This subroutine updates the values for the residuals,
!           for the probes, surface and volume monitors.
!        ***************************************************************
!        
         use PhysicsStorage
         use StopwatchClass
         implicit none
         class(Monitor_t)    :: self
         class(HexMesh)      :: mesh
         real(kind=RP)       :: t
         integer             :: iter
         real(kind=RP)       :: maxResiduals(NCONS), dt
         logical             :: Autosave
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                       :: i 
!
!        Move to next buffer line
!        ------------------------
         self % bufferLine = self % bufferLine + 1
!
!        Save time, iteration and CPU-time
!        -----------------------
         self % t       ( self % bufferLine )  = t
         self % iter    ( self % bufferLine )  = iter
         self % SolverSimuTime ( self % bufferLine )  = Stopwatch % ElapsedTime("Solver")
         self % TotalSimuTime ( self % bufferLine )  = Stopwatch % ElapsedTime("TotalTime")
!
!        Compute current residuals
!        -------------------------
         call self % residuals % Update( mesh, maxResiduals, self % bufferLine )
!
!        Update volume monitors
!        ----------------------
         do i = 1 , self % no_of_volumeMonitors
            call self % volumeMonitors(i) % Update( mesh , self % bufferLine )
         end do
!
!        Update load balancing monitors
!        ------------------------------
         do i = 1 , self % no_of_loadBalancingMonitors
            call self % loadBalancingMonitors(i) % Update( mesh , self % bufferLine )
         end do

#ifdef FLOW
!
!        Update standard probes (GPU path with per-probe MPI)
!        -----------------------------------------------------
         do i = 1 , self % no_of_probes - self % no_of_fileProbes
            call self % probes(i) % Update( mesh , self % bufferLine )
         end do
!
!        Update file probes: OpenMP compute + single MPI_Allreduce
!        ----------------------------------------------------------
         if ( self % no_of_fileProbes .gt. 0 ) then
            call Monitor_UpdateFileProbes( self, mesh, self % bufferLine )
         end if
#endif

#if defined(NAVIERSTOKES) || defined(INCNS)
!
!        Update surface monitors
!        -----------------------
         do i = 1 , self % no_of_surfaceMonitors
            call self % surfaceMonitors(i) % Update( mesh , self % bufferLine, iter, autosave, dt )
         end do
!
!        Update statistics
!        -----------------
         call self % stats % Update(mesh, iter, t, trim(self % solution_file) )
#endif

!
!        Update dt restriction
!        ---------------------
         if (self % write_dt_restriction) self % dt_restriction = mesh % dt_restriction 
         
      end subroutine Monitor_UpdateValues

      subroutine Monitor_WriteToFile ( self , mesh, force) 
!
!        ******************************************************************
!              This routine has a double behaviour:
!           force = .true.  -> Writes to file and resets buffers
!           force = .false. -> Just writes to file if the buffer is full
!        ******************************************************************
!
         use MPI_Process_Info
         implicit none
         class(Monitor_t)        :: self
         class(HexMesh)          :: mesh
         logical, optional       :: force
!        ------------------------------------------------
         integer                 :: i 
         logical                 :: forceVal

         if ( present ( force ) ) then
            forceVal = force

         else
            forceVal = .false.

         end if

         if ( forceVal ) then 
!
!           In this case the monitors are exported to their files and the buffer is reset
!           -----------------------------------------------------------------------------
            call self % residuals % WriteToFile ( self % iter , self % t, self % TotalSimuTime, self % SolverSimuTime , self % bufferLine )
   
            do i = 1 , self % no_of_volumeMonitors
               call self % volumeMonitors(i) % WriteToFile ( self % iter , self % t , self % bufferLine )
            end do

            do i = 1 , self % no_of_loadBalancingMonitors
               call self % loadBalancingMonitors(i) % WriteToFile ( self % iter , self % t , self % bufferLine )
            end do

#ifdef FLOW
            do i = 1 , self % no_of_probes - self % no_of_fileProbes
               call self % probes(i) % WriteToFile ( self % iter , self % t , self % bufferLine )
            end do
            if ( self % no_of_fileProbes .gt. 0 ) then
#ifdef HAS_HDF5
               if ( trim(self % probeFileOutputFormat) .eq. "HDF5" ) then
                  call Monitor_WriteFileProbesHDF5( self, self % iter, self % t, self % bufferLine )
               else
#endif
                  do i = self % no_of_probes - self % no_of_fileProbes + 1, self % no_of_probes
                     call self % probes(i) % WriteToFile ( self % iter , self % t , self % bufferLine )
                  end do
#ifdef HAS_HDF5
               end if
#endif
            end if
#endif

#if defined(NAVIERSTOKES) || defined(INCNS)
            do i = 1 , self % no_of_surfaceMonitors
               call self % surfaceMonitors(i) % WriteToFile ( self % iter , self % t , self % bufferLine )
            end do
!
!              Write statistics
!              ----------------
            if ( self % bufferLine .eq. 0 ) then
               i = 1
            else
               i = self % bufferLine
            end if
            call self % stats % WriteFile(mesh, self % iter(i), self % t(i), self % solution_file)
#endif
!
!           Reset buffer
!           ------------
            self % bufferLine = 0

         else
!
!           The monitors are exported just if the buffer is full
!           ----------------------------------------------------
            if ( self % bufferLine .eq. BUFFER_SIZE ) then

               call self % residuals % WriteToFile ( self % iter , self % t, self % TotalSimuTime, self % SolverSimuTime, BUFFER_SIZE )

               do i = 1 , self % no_of_volumeMonitors
                  call self % volumeMonitors(i) % WriteToFile ( self % iter , self % t , self % bufferLine )
               end do

               do i = 1 , self % no_of_loadBalancingMonitors
                  call self % loadBalancingMonitors(i) % WriteToFile ( self % iter , self % t , self % bufferLine )
               end do

#ifdef FLOW
               do i = 1 , self % no_of_probes - self % no_of_fileProbes
                  call self % probes(i) % WriteToFile ( self % iter , self % t , self % bufferLine )
               end do
               if ( self % no_of_fileProbes .gt. 0 ) then
#ifdef HAS_HDF5
                  if ( trim(self % probeFileOutputFormat) .eq. "HDF5" ) then
                     call Monitor_WriteFileProbesHDF5( self, self % iter, self % t, self % bufferLine )
                  else
#endif
                     do i = self % no_of_probes - self % no_of_fileProbes + 1, self % no_of_probes
                        call self % probes(i) % WriteToFile ( self % iter , self % t , self % bufferLine )
                     end do
#ifdef HAS_HDF5
                  end if
#endif
               end if
#endif

#if defined(NAVIERSTOKES) || defined(INCNS)
               do i = 1 , self % no_of_surfaceMonitors
                  call self % surfaceMonitors(i) % WriteToFile ( self % iter , self % t , self % bufferLine )
               end do
#endif
!
!              Reset buffer
!              ------------
               self % bufferLine = 0

            end if
         end if

      end subroutine Monitor_WriteToFile

      subroutine Monitor_WriteProbesFileSummary ( self )
!
!        ********************************************************************
!              Prints a startup-log summary of the bulk probes-file monitor,
!           shown instead of including the file-based probes in the
!           per-iteration screen log.
!        ********************************************************************
!
         use Headers
         use MPI_Process_Info
         implicit none
         class(Monitor_t)        :: self
!
!        ---------------
!        Local variables
!        ---------------
!
         integer :: j

         if ( .not. MPI_Process % isRoot ) return
         if ( self % no_of_fileProbes .le. 0 ) return

         write(STD_OUT,'(/)')
         call Section_Header("Probes file monitor")
         write(STD_OUT,'(/)')

         write(STD_OUT,'(30X,A,A28,A)')   "->" , "File: " , trim(self % probesFileName)
         write(STD_OUT,'(30X,A,A28,I10)') "->" , "Number of probes: " , self % no_of_fileProbes
         write(STD_OUT,'(30X,A,A28)',advance="no") "->" , "Variables: "
         do j = 1 , size(self % probesVariables)
            write(STD_OUT,'(A)',advance="no") trim(self % probesVariables(j)) // " "
         end do
         write(STD_OUT,*)
         if ( self % probeFileSaveTimestep .gt. 0.0_RP ) then
            write(STD_OUT,'(30X,A,A28,ES14.6)') "->" , "Save timestep: " , self % probeFileSaveTimestep
         end if
         write(STD_OUT,'(30X,A,A28,A)') "->" , "Output format: " , trim(self % probeFileOutputFormat)

      end subroutine Monitor_WriteProbesFileSummary

      subroutine Monitor_Destruct (self)
         implicit none
         class(Monitor_t)        :: self
         
         deallocate (self % iter)
         deallocate (self % t)
         deallocate (self % TotalSimuTime)
         deallocate (self % SolverSimuTime)
         
         call self % residuals % destruct
         
         call self % volumeMonitors % destruct
         safedeallocate(self % volumeMonitors)

         call self % loadBalancingMonitors % destruct
         safedeallocate(self % loadBalancingMonitors)
         
#ifdef FLOW
         call self % probes % destruct
         safedeallocate (self % probes)
#endif
         
#if defined(NAVIERSTOKES) || defined(INCNS)
         call self % surfaceMonitors % destruct
         safedeallocate (self % surfaceMonitors)
         
         !call self % stats % destruct
#endif         
      end subroutine
      
      impure elemental subroutine Monitor_Assign ( to, from )
         implicit none
         !-arguments--------------------------------------
         class(Monitor_t), intent(inout)  :: to
         type(Monitor_t) , intent(in)     :: from
         !-local-variables--------------------------------
         !------------------------------------------------
         
         to % solution_file               = from % solution_file
         to % probes_solution_file        = from % probes_solution_file
         to % no_of_probes                = from % no_of_probes
         to % no_of_surfaceMonitors       = from % no_of_surfaceMonitors
         to % no_of_volumeMonitors        = from % no_of_volumeMonitors
         to % no_of_loadBalancingMonitors = from % no_of_loadBalancingMonitors
         to % no_of_fileProbes            = from % no_of_fileProbes
         to % probesFileName              = from % probesFileName
         to % probeFileSaveTimestep       = from % probeFileSaveTimestep
         to % probeFileOutputFormat       = from % probeFileOutputFormat
         to % fp_lastSavedTime            = from % fp_lastSavedTime
         if ( allocated(from % probesVariables) ) then
            safedeallocate ( to % probesVariables )
            allocate ( to % probesVariables ( size(from % probesVariables) ) )
            to % probesVariables = from % probesVariables
         end if
         to % bufferLine                  = from % bufferLine
         
         safedeallocate ( to % iter )
         allocate ( to % iter ( size(from % iter) ) )
         to % iter = from % iter
         
         to % dt_restriction        = from % dt_restriction
         to % write_dt_restriction  = from % write_dt_restriction
         
         safedeallocate (to % t)
         allocate (to % t (size (from % t) ) ) 
         to % t = from % t
         
         safedeallocate ( to % TotalSimuTime )
         allocate ( to % TotalSimuTime ( size(from % TotalSimuTime) ) )
         to % TotalSimuTime = from % TotalSimuTime
         
         safedeallocate ( to % SolverSimuTime )
         allocate ( to % SolverSimuTime ( size(from % SolverSimuTime) ) )
         to % SolverSimuTime = from % SolverSimuTime
         
         to % residuals = from % residuals
         
         safedeallocate ( to % volumeMonitors )
         allocate ( to % volumeMonitors ( size(from % volumeMonitors) ) )
         to % volumeMonitors = from % volumeMonitors

         safedeallocate ( to % loadBalancingMonitors )
         allocate ( to % loadBalancingMonitors ( size(from % loadBalancingMonitors) ) )
         to % loadBalancingMonitors = from % loadBalancingMonitors
      
#ifdef FLOW
         safedeallocate ( to % probes )
         allocate ( to % probes ( size(from % probes) ) )
         to % probes = from % probes
#endif
         
#if defined(NAVIERSTOKES) || defined(INCNS)
         safedeallocate ( to % surfaceMonitors )
         allocate ( to % surfaceMonitors ( size(from % surfaceMonitors) ) )
         to % surfaceMonitors = from % surfaceMonitors
         
         to % stats = from % stats
#endif
         
      end subroutine Monitor_Assign
      
!
!//////////////////////////////////////////////////////////////////////////////
!
!        Auxiliars
!
!//////////////////////////////////////////////////////////////////////////////
!
   subroutine getNoOfMonitors(no_of_probes, no_of_surfaceMonitors, no_of_volumeMonitors, no_of_loadBalancingMonitors)
      use ParamfileRegions
      implicit none
      integer, intent(out)    :: no_of_probes
      integer, intent(out)    :: no_of_surfaceMonitors
      integer, intent(out)    :: no_of_volumeMonitors
      integer, intent(out)    :: no_of_loadBalancingMonitors
!
!     ---------------
!     Local variables
!     ---------------
!
      character(len=LINE_LENGTH) :: case_name, line
      integer                    :: fID
      integer                    :: io
!
!     Initialize
!     ----------
      no_of_probes = 0
      no_of_surfaceMonitors = 0
      no_of_volumeMonitors = 0
      no_of_loadBalancingMonitors = 0
!
!     Get case file name
!     ------------------
      call get_command_argument(1, case_name)

!
!     Open case file
!     --------------
      open ( newunit = fID , file = case_name , status = "old" , action = "read" )

!
!     Read the whole file to find monitors
!     ------------------------------------
readloop:do 
         read ( fID , '(A)' , iostat = io ) line

         if ( io .lt. 0 ) then
!
!           End of file
!           -----------
            line = ""
            exit readloop

         elseif ( io .gt. 0 ) then
!
!           Error
!           -----
            errorMessage(STD_OUT)
            error stop "Stopped."

         else
!
!           Succeeded
!           ---------
            line = getSquashedLine( line )

            if ( index ( line , '#defineprobefile' ) .gt. 0 ) then
!
!              The probe-file block is not an individual probe definition
!              -----------------------------------------------------------

            elseif ( index ( line , '#defineprobe' ) .gt. 0 ) then
               no_of_probes = no_of_probes + 1

            elseif ( index ( line , '#definesurfacemonitor' ) .gt. 0 ) then
               no_of_surfaceMonitors = no_of_surfaceMonitors + 1 

            elseif ( index ( line , '#definevolumemonitor' ) .gt. 0 ) then
               no_of_volumeMonitors = no_of_volumeMonitors + 1 

            elseif ( index ( line , '#defineloadbalancingmonitor' ) .gt. 0 ) then
               no_of_loadBalancingMonitors = no_of_loadBalancingMonitors + 1

            end if
            
         end if

      end do readloop
!
!     Close case file
!     ---------------
      close(fID)                             

end subroutine getNoOfMonitors

!
!///////////////////////////////////////////////////////////////////////////////////
!
!     Probes-from-file auxiliary routines
!
!     The probes file is a plain text file. Blank lines and lines starting
!     with "#" are ignored. Every remaining line holds the coordinates of
!     one probe:
!
!        x  y  z
!
!     The list of variables to sample (shared by every probe in the file)
!     is given through the "probes file variables" control-file keyword.
!
!     A separate output file "<solution_file>.probe_<N>.probe" is created for
!     each probe, with one column per variable, and one row written per
!     saved time step.
!
!///////////////////////////////////////////////////////////////////////////////////
!
   subroutine readProbesFileBlock(fileName, variablesLine, no_of_probesVariables, saveTimestep, outputFormat)
!
!     ******************************************************************
!        Reads the "#define probe file ... #end" block from the case
!     file, if present:
!
!        #define probe file
!           file                = Probe.dat
!           variables           = u
!           probe save timestep = 1.0E-3
!           output format       = HDF5
!        #end
!
!     Note: this is a hand-rolled parser (rather than readValueInRegion)
!     because readValueInRegion lower-cases every line it reads, which
!     would corrupt a case-sensitive file path.
!     ******************************************************************
!
      use ParamfileRegions, only: getSquashedLine
      implicit none
      character(len=LINE_LENGTH), intent(out) :: fileName
      character(len=LINE_LENGTH), intent(out) :: variablesLine
      integer,                    intent(out) :: no_of_probesVariables
      real(kind=RP),              intent(out) :: saveTimestep
      character(len=8),           intent(out) :: outputFormat
!
!     ---------------
!     Local variables
!     ---------------
!
      character(len=LINE_LENGTH) :: paramFile
      character(len=LINE_LENGTH) :: line, squashed, valStr
      integer                    :: fID, io, position

      logical                    :: inside

      fileName              = ""
      variablesLine         = ""
      no_of_probesVariables = 0
      saveTimestep          = 0.0_RP
      outputFormat          = "ASCII"
      inside                = .false.

      call get_command_argument(1, paramFile)

      open ( newunit = fID , file = trim(paramFile) , status = "old" , action = "read" )

      do
         read ( fID , '(A)' , iostat = io ) line
         if ( io .ne. 0 ) exit

         squashed = getSquashedLine(line)

         if ( squashed .eq. getSquashedLine("#define probe file") ) then
            inside = .true.
            cycle
         elseif ( squashed .eq. getSquashedLine("#end") ) then
            inside = .false.
            cycle
         end if

         if ( .not. inside ) cycle
!
!        Strip a trailing comment, keeping the value's original case
!        -------------------------------------------------------------
         position = index(line , '!')
         if ( position .gt. 0 ) line = line(1:position-1)

         position = max( index(line,'='), index(line,':') )
         if ( position .eq. 0 ) cycle

         if ( getSquashedLine(line(1:position-1)) .eq. getSquashedLine("file") ) then
            fileName = adjustl( removeQuotes( line(position+1:) ) )
         elseif ( getSquashedLine(line(1:position-1)) .eq. getSquashedLine("variables") ) then
            variablesLine = adjustl( removeQuotes( line(position+1:) ) )
         elseif ( getSquashedLine(line(1:position-1)) .eq. getSquashedLine("probe save timestep") ) then
            valStr = adjustl( removeQuotes( line(position+1:) ) )
            read( valStr , * ) saveTimestep
         elseif ( getSquashedLine(line(1:position-1)) .eq. getSquashedLine("output format") ) then
            valStr = adjustl( removeQuotes( line(position+1:) ) )
            valStr = adjustl( getSquashedLine(valStr) )
            if ( index(valStr, "HDF5") .gt. 0 ) then
               outputFormat = "HDF5"
            else
               outputFormat = "ASCII"
            end if
         end if
      end do

      close(fID)

   end subroutine readProbesFileBlock

   function removeQuotes(str) result(res)
      implicit none
      character(len=*), intent(in) :: str
      character(len=LINE_LENGTH)   :: res
      character(len=LINE_LENGTH)   :: auxstr
      integer                      :: i, j

      auxstr = trim(adjustl(str))
      res    = ""
      j      = 0

      do i = 1 , len_trim(auxstr)
         if ( auxstr(i:i) .eq. '"' .or. auxstr(i:i) .eq. "'" ) cycle
         j = j + 1
         res(j:j) = auxstr(i:i)
      end do

   end function removeQuotes

   subroutine countProbesInFile(fileName, n)
      implicit none
      character(len=*), intent(in)  :: fileName
      integer,          intent(out) :: n
!
!     ---------------
!     Local variables
!     ---------------
!
      integer                    :: fID, io
      character(len=LINE_LENGTH) :: line

      n = 0
      open ( newunit = fID , file = fileName , status = "old" , action = "read" )

      do
         read ( fID , '(A)' , iostat = io ) line
         if ( io .ne. 0 ) exit

         line = adjustl(line)
         if ( len_trim(line) .eq. 0 ) cycle
         if ( line(1:1) .eq. '#' )    cycle

         n = n + 1
      end do

      close(fID)

   end subroutine countProbesInFile

   subroutine splitIntoTokens(line, tokens, n)
      implicit none
      character(len=*),                              intent(in)  :: line
      character(len=STR_LEN_MONITORS), allocatable,   intent(out) :: tokens(:)
      integer,                                        intent(out) :: n
!
!     ---------------
!     Local variables
!     ---------------
!
      character(len=LINE_LENGTH)      :: auxline
      character(len=STR_LEN_MONITORS) :: buffer(64)
      integer                         :: pos

      auxline = adjustl(line)
      n = 0

      do while ( len_trim(auxline) .gt. 0 )
         pos = index(trim(auxline), " ")
         n = n + 1

         if ( pos .eq. 0 ) then
            buffer(n) = trim(auxline)
            auxline   = ""
         else
            buffer(n) = auxline(1:pos-1)
            auxline   = adjustl(auxline(pos+1:))
         end if
      end do

      allocate( tokens(n) )
      tokens(1:n) = buffer(1:n)

   end subroutine splitIntoTokens

#ifdef FLOW
   subroutine Monitor_UpdateFileProbes(self, mesh, bufferPos)
!
!     Evaluates all file-probes using OpenMP threads (CPU path, no OpenACC),
!     then performs a single MPI_Allreduce(SUM) to share results across ranks.
!     This replaces N_fileProbes individual MPI_Bcast calls with one collective.
!
      use MPI_Process_Info
#ifdef _HAS_MPI_
      use mpi
#endif
      implicit none
      class(Monitor_t), intent(inout) :: self
      class(HexMesh),   intent(in)    :: mesh
      integer,          intent(in)    :: bufferPos
!
!     ---------------
!     Local variables
!     ---------------
!
      integer        :: i, v, j, nfp, nv, fp_offset, ierr
      real(kind=RP), allocatable :: buf(:)

      nfp       = self % no_of_fileProbes
      nv        = size(self % probesVariables)
      fp_offset = self % no_of_probes - nfp

      ! Parallel CPU computation — each probe writes to its own values(:,bufferPos).
      ! Non-owning ranks store 0 so MPI_Allreduce(SUM) gives the correct result.
      !$omp parallel do schedule(dynamic,16) default(shared)
      do i = fp_offset + 1, self % no_of_probes
         call self % probes(i) % ComputeLocal(mesh, bufferPos)
      end do
      !$omp end parallel do

#ifdef _HAS_MPI_
      if ( MPI_Process % doMPIAction ) then
         ! Pack all file-probe values into a flat buffer
         allocate( buf(nfp * nv) )
         j = 0
         do i = fp_offset + 1, self % no_of_probes
            do v = 1, nv
               j = j + 1
               buf(j) = self % probes(i) % values(v, bufferPos)
            end do
         end do

         ! Single collective instead of nfp individual Bcasts
         call MPI_Allreduce(MPI_IN_PLACE, buf, nfp * nv, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

         ! Unpack
         j = 0
         do i = fp_offset + 1, self % no_of_probes
            do v = 1, nv
               j = j + 1
               self % probes(i) % values(v, bufferPos) = buf(j)
            end do
         end do
         deallocate(buf)
      end if
#endif

   end subroutine Monitor_UpdateFileProbes

   subroutine InitializeProbesFromFile(fileName, probes, offset, mesh, solution_file, FirstCall, variables, saveTimestep, outputFormat)
      implicit none
      character(len=*),   intent(in)    :: fileName
      class(Probe_t),     intent(inout) :: probes(:)
      integer,            intent(in)    :: offset
      class(HexMesh),     intent(in)    :: mesh
      character(len=*),   intent(in)    :: solution_file
      logical,            intent(in)    :: FirstCall
      character(len=*),   intent(in)    :: variables(:)
      real(kind=RP),      intent(in)    :: saveTimestep
      character(len=*),   intent(in)    :: outputFormat
!
!     ---------------
!     Local variables
!     ---------------
!
      integer                                      :: fID, io, idx, nTok
      character(len=LINE_LENGTH)                   :: line
      real(kind=RP)                                :: x(NDIM)
      character(len=STR_LEN_MONITORS), allocatable  :: tokens(:)
      character(len=STR_LEN_MONITORS)               :: pname

      idx = offset
      open ( newunit = fID , file = fileName , status = "old" , action = "read" )

      do
         read ( fID , '(A)' , iostat = io ) line
         if ( io .ne. 0 ) exit

         line = adjustl(line)
         if ( len_trim(line) .eq. 0 ) cycle
         if ( line(1:1) .eq. '#' )    cycle

         call splitIntoTokens(line, tokens, nTok)

         if ( nTok .lt. 3 ) then
            deallocate(tokens)
            cycle
         end if

         read(tokens(1),*) x(1)
         read(tokens(2),*) x(2)
         read(tokens(3),*) x(3)

         idx = idx + 1
         write(pname,'(A,I0)') "probe_", idx

         call probes(idx) % Initialization( mesh, idx, solution_file, FirstCall, &
                                             x_in = x, variables_in = variables, name_in = trim(pname), &
                                             isFileProbe_in = .true., outputFormat_in = trim(outputFormat) )
         probes(idx) % saveTimestep = saveTimestep

         deallocate(tokens)
      end do

      close(fID)

   end subroutine InitializeProbesFromFile

! ============================================================
!  HDF5 collective I/O for file-based probes
!  Compiled only when HAS_HDF5 is defined.
! ============================================================

#ifdef HAS_HDF5
   subroutine Monitor_InitFileProbesHDF5(self, no_of_fileProbes)
!
!     Creates the HDF5 output file for bulk file-probes with
!     the following structure:
!
!       /coordinates  (3, nProbes)   — fixed, written once
!       /time         (extendible)   — one value per saved step
!       /iteration    (extendible)   — one value per saved step
!       /<varname>    (nProbes, ext) — one row per saved step
!
      use HDF5
      use MPI_Process_Info
      implicit none
      class(Monitor_t), intent(inout) :: self
      integer,          intent(in)    :: no_of_fileProbes
!
!     ---------------
!     Local variables
!     ---------------
!
      integer(HID_T)   :: file_id, dset_id, dspace_id, dcpl_id
      integer(HSIZE_T) :: dims2(2), maxdims2(2), chunk2(2)
      integer(HSIZE_T) :: dims1(1), maxdims1(1), chunk1(1)
      integer          :: iError, j, v, fp_offset, nv
      character(len=LINE_LENGTH) :: fname
      real(kind=RP), allocatable :: coords(:,:)

      if ( .not. MPI_Process % isRoot ) return

      nv        = size(self % probesVariables)
      fp_offset = self % no_of_probes  ! probes are allocated 1..no_of_probes+no_of_fileProbes
                                        ! but at this point no_of_probes has not yet been bumped

      write(fname,'(A,A)') trim(self % probes_solution_file), ".probes.h5"

      call h5open_f(iError)
      call h5fcreate_f(trim(fname), H5F_ACC_TRUNC_F, file_id, iError)

      ! /coordinates  (3, nProbes) — fixed at creation
      allocate( coords(3, no_of_fileProbes) )
      do j = 1, no_of_fileProbes
         coords(:, j) = self % probes(fp_offset + j) % x
      end do
      dims2 = [ int(3, HSIZE_T), int(no_of_fileProbes, HSIZE_T) ]
      call h5screate_simple_f(2, dims2, dspace_id, iError)
      call h5dcreate_f(file_id, "coordinates", H5T_NATIVE_DOUBLE, dspace_id, dset_id, iError)
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, coords, dims2, iError)
      call h5dclose_f(dset_id, iError)
      call h5sclose_f(dspace_id, iError)
      deallocate(coords)

      ! /time  (extendible 1-D)
      dims1(1)    = 0
      maxdims1(1) = H5S_UNLIMITED_F
      chunk1(1)   = int(max(BUFFER_SIZE, 1), HSIZE_T)
      call h5screate_simple_f(1, dims1, dspace_id, iError, maxdims1)
      call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl_id, iError)
      call h5pset_chunk_f(dcpl_id, 1, chunk1, iError)
      call h5dcreate_f(file_id, "time", H5T_NATIVE_DOUBLE, dspace_id, dset_id, iError, dcpl_id)
      call h5dclose_f(dset_id, iError)
      call h5sclose_f(dspace_id, iError)
      call h5pclose_f(dcpl_id, iError)

      ! /iteration  (extendible 1-D)
      call h5screate_simple_f(1, dims1, dspace_id, iError, maxdims1)
      call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl_id, iError)
      call h5pset_chunk_f(dcpl_id, 1, chunk1, iError)
      call h5dcreate_f(file_id, "iteration", H5T_NATIVE_INTEGER, dspace_id, dset_id, iError, dcpl_id)
      call h5dclose_f(dset_id, iError)
      call h5sclose_f(dspace_id, iError)
      call h5pclose_f(dcpl_id, iError)

      ! /<varname>  (extendible × nProbes) — h5ls shows {nProbes, Inf}
      ! chunk2(1)=1: one time step per chunk avoids pre-allocation waste when
      ! the save-timestep filter makes n_write << BUFFER_SIZE per flush.
      dims2(1)    = 0
      dims2(2)    = int(no_of_fileProbes, HSIZE_T)
      maxdims2(1) = H5S_UNLIMITED_F
      maxdims2(2) = int(no_of_fileProbes, HSIZE_T)
      chunk2(1)   = int(1, HSIZE_T)
      chunk2(2)   = int(no_of_fileProbes, HSIZE_T)

      do v = 1, nv
         call h5screate_simple_f(2, dims2, dspace_id, iError, maxdims2)
         call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl_id, iError)
         call h5pset_chunk_f(dcpl_id, 2, chunk2, iError)
         call h5dcreate_f(file_id, trim(self % probesVariables(v)), H5T_NATIVE_DOUBLE, &
                          dspace_id, dset_id, iError, dcpl_id)
         call h5dclose_f(dset_id, iError)
         call h5sclose_f(dspace_id, iError)
         call h5pclose_f(dcpl_id, iError)
      end do

      call h5fclose_f(file_id, iError)
      call h5close_f(iError)

   end subroutine Monitor_InitFileProbesHDF5

   subroutine Monitor_WriteFileProbesHDF5(self, iter, t, no_of_lines)
!
!     Appends one buffer of file-probe data to the HDF5 output file.
!     Only lines that satisfy the saveTimestep filter are written.
!     Called from Monitor_WriteToFile in place of the per-probe ASCII loop.
!
      use HDF5
      use MPI_Process_Info
      implicit none
      class(Monitor_t), intent(inout) :: self
      integer,          intent(in)    :: iter(:)
      real(kind=RP),    intent(in)    :: t(:)
      integer,          intent(in)    :: no_of_lines
!
!     ---------------
!     Local variables
!     ---------------
!
      integer(HID_T)   :: file_id, dset_id, dspace_id, mspace_id
      integer(HSIZE_T) :: cur1(1), max1(1), off1(1), cnt1(1)
      integer(HSIZE_T) :: new1(1)
      integer(HSIZE_T) :: cur2(2), max2(2), off2(2), cnt2(2)
      integer(HSIZE_T) :: new2(2)
      integer          :: iError, i, j, k, v, nfp, nv, fp_offset, n_write
      character(len=LINE_LENGTH) :: fname
      logical,      allocatable :: wmask(:)
      integer,      allocatable :: ibuf(:)
      real(kind=RP), allocatable :: rbuf(:), vbuf(:)

      if ( .not. MPI_Process % isRoot ) return

      nfp       = self % no_of_fileProbes
      nv        = size(self % probesVariables)
      fp_offset = self % no_of_probes - nfp

      ! Build write mask (apply saveTimestep filter)
      allocate( wmask(no_of_lines) )
      wmask = .false.
      do i = 1, no_of_lines
         if ( self % probeFileSaveTimestep .gt. 0.0_RP ) then
            if ( t(i) .lt. self % fp_lastSavedTime + self % probeFileSaveTimestep ) cycle
         end if
         wmask(i) = .true.
         self % fp_lastSavedTime = t(i)
      end do
      n_write = count(wmask)
      if ( n_write .eq. 0 ) then
         deallocate(wmask)
         return
      end if

      write(fname,'(A,A)') trim(self % probes_solution_file), ".probes.h5"

      call h5open_f(iError)
      call h5fopen_f(trim(fname), H5F_ACC_RDWR_F, file_id, iError)

      ! Query current extent of /time to get append offset
      call h5dopen_f(file_id, "time", dset_id, iError)
      call h5dget_space_f(dset_id, dspace_id, iError)
      call h5sget_simple_extent_dims_f(dspace_id, cur1, max1, iError)
      call h5sclose_f(dspace_id, iError)
      call h5dclose_f(dset_id, iError)
      ! cur1(1) = number of time steps already written

      cnt1(1) = int(n_write, HSIZE_T)
      off1(1) = cur1(1)
      new1(1) = cur1(1) + cnt1(1)

      ! Pack filtered iteration values
      allocate( ibuf(n_write) )
      j = 0
      do i = 1, no_of_lines
         if ( .not. wmask(i) ) cycle
         j = j + 1
         ibuf(j) = iter(i)
      end do

      call h5dopen_f(file_id, "iteration", dset_id, iError)
      call h5dextend_f(dset_id, new1, iError)
      call h5dget_space_f(dset_id, dspace_id, iError)
      call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, off1, cnt1, iError)
      call h5screate_simple_f(1, cnt1, mspace_id, iError)
      call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, ibuf, cnt1, iError, mspace_id, dspace_id)
      call h5sclose_f(mspace_id, iError)
      call h5sclose_f(dspace_id, iError)
      call h5dclose_f(dset_id, iError)
      deallocate(ibuf)

      ! Pack filtered time values
      allocate( rbuf(n_write) )
      j = 0
      do i = 1, no_of_lines
         if ( .not. wmask(i) ) cycle
         j = j + 1
         rbuf(j) = t(i)
      end do

      call h5dopen_f(file_id, "time", dset_id, iError)
      call h5dextend_f(dset_id, new1, iError)
      call h5dget_space_f(dset_id, dspace_id, iError)
      call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, off1, cnt1, iError)
      call h5screate_simple_f(1, cnt1, mspace_id, iError)
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, rbuf, cnt1, iError, mspace_id, dspace_id)
      call h5sclose_f(mspace_id, iError)
      call h5sclose_f(dspace_id, iError)
      call h5dclose_f(dset_id, iError)
      deallocate(rbuf)

      ! Write each variable: one extend + one write per flush (avoids chunk waste).
      ! vbuf2(nfp, n_write) packs all filtered timesteps into a contiguous 2D block.
      off2(1) = cur1(1)             ! append after existing time steps
      off2(2) = int(0,   HSIZE_T)
      cnt2(1) = int(n_write, HSIZE_T)
      cnt2(2) = int(nfp,    HSIZE_T)
      new2(1) = cur1(1) + int(n_write, HSIZE_T)
      new2(2) = int(nfp, HSIZE_T)

      allocate( vbuf(nfp * n_write) )

      do v = 1, nv
         call h5dopen_f(file_id, trim(self % probesVariables(v)), dset_id, iError)

         ! Pack into column-major 2D layout vbuf(n_write, nfp):
         ! time index j varies fastest (Fortran column-major with cnt2=[n_write,nfp])
         j = 0
         do i = 1, no_of_lines
            if ( .not. wmask(i) ) cycle
            j = j + 1
            do k = 1, nfp
               vbuf( j + (k-1)*n_write ) = self % probes(fp_offset + k) % values(v, i)
            end do
         end do

         call h5dextend_f(dset_id, new2, iError)
         call h5dget_space_f(dset_id, dspace_id, iError)
         call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, off2, cnt2, iError)
         call h5screate_simple_f(2, cnt2, mspace_id, iError)
         call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, vbuf, cnt2, iError, mspace_id, dspace_id)
         call h5sclose_f(mspace_id, iError)
         call h5sclose_f(dspace_id, iError)

         call h5dclose_f(dset_id, iError)
      end do

      deallocate(vbuf)
      deallocate(wmask)

      call h5fclose_f(file_id, iError)
      call h5close_f(iError)

   end subroutine Monitor_WriteFileProbesHDF5
#endif   ! HAS_HDF5

#endif   ! FLOW

end module MonitorsClass
!
!///////////////////////////////////////////////////////////////////////////////////
!
