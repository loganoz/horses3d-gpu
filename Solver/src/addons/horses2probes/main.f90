program horses2probes
   use SMConstants
   use Headers
   use MPI_Process_Info
   use SharedSpectralBasis
#ifdef FLOW
   use Solution2ProbesModule, only: RunSolution2ProbesFromControl
#endif
   implicit none

   character(len=LINE_LENGTH) :: controlFile

   call MPI_Process % Init
   call Main_Header("HORSES probes post-processing utility",__DATE__,__TIME__)
   call ConstructSpectralBasis

   call get_command_argument(1, controlFile)
   if (len_trim(controlFile) == 0) then
      if (MPI_Process % isRoot) write(STD_OUT,'(A)') "Usage: horses2probes <control_file>"
      call DestructSpectralBasis
      call MPI_Process % Close
      error stop
   end if

#ifndef FLOW
   if (MPI_Process % isRoot) write(STD_OUT,'(A)') "ERROR: horses2probes compiled without -DFLOW"
   call DestructSpectralBasis
   call MPI_Process % Close
   error stop
#else
   call RunSolution2ProbesFromControl(trim(controlFile))
   call DestructSpectralBasis
   call MPI_Process % Close
#endif

end program horses2probes
