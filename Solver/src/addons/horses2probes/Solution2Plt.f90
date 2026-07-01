module Solution2ProbesModule
#ifdef FLOW
   use SMConstants
   use Storage                    , only: Mesh_t, Element_t
   use PhysicsStorage              , only: IRHO, IRHOU, IRHOV, IRHOW
   use VariableConversion         , only: Pressure
   use NodalStorageClass          , only: NodalStorage
   use SharedSpectralBasis        , only: spA, addNewSpectralBasis
   use Utilities                  , only: SolveThreeEquationLinearSystem
   use MPI_Process_Info
   implicit none

   private
   public :: RunSolution2ProbesFromControl

   integer, parameter :: MAX_VARS = 16

contains

   subroutine RunSolution2ProbesFromControl(controlFile)
      character(len=*), intent(in) :: controlFile

      character(len=LINE_LENGTH) :: meshName, solName, probesName, outName, outFmt
      character(len=LINE_LENGTH) :: varsLine
      character(len=32)          :: varNames(MAX_VARS)
      integer :: nVars

      type(Mesh_t) :: mesh
      real(kind=RP), allocatable :: probes(:,:)    ! (3,nP)
      real(kind=RP), allocatable :: outVals(:,:)   ! (nVars,nP)
      logical,      allocatable  :: found(:)
      integer :: nP, fid

      call ReadControlFile(controlFile, meshName, solName, probesName, outName, outFmt, varsLine)
      call ParseVarList(varsLine, varNames, nVars)

      if (MPI_Process % isRoot) then
         write(STD_OUT,'(A)') "horses2probes: mesh   = "//trim(meshName)
         write(STD_OUT,'(A)') "horses2probes: hsol   = "//trim(solName)
         write(STD_OUT,'(A)') "horses2probes: probes = "//trim(probesName)
      end if

      call mesh % ReadMesh(trim(meshName))
      call mesh % ReadSolution(trim(solName))

      call ReadProbesXYZ(trim(probesName), probes, nP)
      allocate(outVals(nVars, nP))
      allocate(found(nP))

      call InterpolateAll(mesh, probes, varNames, nVars, outVals, found)

      if (MPI_Process % isRoot) then
         open(newunit=fid, file=trim(outName), status="unknown", action="write")
         if (trim(ToLower(outFmt)) == "tec") then
            call WriteTecplot(fid, probes, varNames, nVars, outVals, found)
         else
            call WriteTxt(fid, probes, varNames, nVars, outVals, found)
         end if
         close(fid)
         write(STD_OUT,'(A)') "horses2probes: wrote "//trim(outName)
      end if
   end subroutine RunSolution2ProbesFromControl


   ! ==========================================================
   ! Control file parsing
   ! ==========================================================
   subroutine ReadControlFile(fname, meshName, solName, probesName, outName, outFmt, varsLine)
      character(len=*), intent(in)  :: fname
      character(len=*), intent(out) :: meshName, solName, probesName, outName, outFmt, varsLine

      integer :: fid, ios
      character(len=LINE_LENGTH) :: line, key, val

      meshName   = ""
      solName    = ""
      probesName = ""
      outName    = "probes.tec"
      outFmt     = "tec"
      varsLine   = "rho,u,v,w,p"

      open(newunit=fid, file=trim(fname), status="old", action="read")
      do
         read(fid,'(A)', iostat=ios) line
         if (ios /= 0) exit
         line = trim(StripComment(line))
         if (len_trim(line) == 0) cycle

         call SplitKeyValue(line, key, val)
         key = trim(ToLower(key))
         val = trim(StripQuotes(val))

         select case (key)
         case ("hmesh file")
            meshName = val
         case ("hsol file")
            solName = val
         case ("probes file")
            probesName = val
         case ("output file")
            outName = val
         case ("output format")
            outFmt = val
         case ("output variables")
            varsLine = val
         end select
      end do
      close(fid)

      if (len_trim(meshName) == 0) error stop "horses2probes: missing 'hmesh file'"
      if (len_trim(solName)  == 0) error stop "horses2probes: missing 'hsol file'"
      if (len_trim(probesName) == 0) error stop "horses2probes: missing 'probes file'"
   end subroutine ReadControlFile

   subroutine SplitKeyValue(line, key, val)
      character(len=*), intent(in)  :: line
      character(len=*), intent(out) :: key, val
      integer :: p
      p = index(line, "=")
      if (p <= 0) then
         key = trim(line)
         val = ""
      else
         key = adjustl(line(1:p-1))
         val = adjustl(line(p+1:))
      end if
   end subroutine SplitKeyValue

   pure function StripQuotes(s) result(out)
      character(len=*), intent(in) :: s
      character(len=len(s)) :: out
      integer :: L
      out = trim(s)
      L = len_trim(out)
      if (L >= 2) then
         if ((out(1:1) == '"' .and. out(L:L) == '"') .or. &
             (out(1:1) == "'" .and. out(L:L) == "'")) then
            out = out(2:L-1)
         end if
      end if
      out = trim(out)
   end function StripQuotes

   pure function StripComment(s) result(out)
      character(len=*), intent(in) :: s
      character(len=len(s)) :: out
      integer :: p
      out = s
      p = index(out, "!")
      if (p > 0) out = out(1:p-1)
      p = index(out, "#")
      if (p > 0) out = out(1:p-1)
      out = trim(out)
   end function StripComment

   pure function ToLower(s) result(out)
      character(len=*), intent(in) :: s
      character(len=len(s)) :: out
      integer :: i, c
      out = s
      do i=1,len(s)
         c = iachar(out(i:i))
         if (c >= iachar('A') .and. c <= iachar('Z')) out(i:i) = achar(c + 32)
      end do
   end function ToLower


   ! ==========================================================
   ! Parse variable list (ONLY: rho,u,v,w,p)
   ! ==========================================================
   subroutine ParseVarList(line, names, n)
      character(len=*), intent(in) :: line
      character(len=32), intent(out) :: names(MAX_VARS)
      integer, intent(out) :: n

      character(len=LINE_LENGTH) :: tmp
      integer :: p

      n = 0
      names = ""

      tmp = trim(line)
      do
         if (len_trim(tmp) == 0) exit
         p = index(tmp, ",")
         if (p == 0) then
            call PushVar(trim(tmp), names, n)
            exit
         else
            call PushVar(trim(tmp(1:p-1)), names, n)
            tmp = adjustl(tmp(p+1:))
         end if
      end do

      if (n == 0) then
         call PushVar("rho", names, n)
         call PushVar("u",   names, n)
         call PushVar("v",   names, n)
         call PushVar("w",   names, n)
         call PushVar("p",   names, n)
      end if
   end subroutine ParseVarList

   subroutine PushVar(v, names, n)
      character(len=*), intent(in) :: v
      character(len=32), intent(inout) :: names(MAX_VARS)
      integer, intent(inout) :: n
      character(len=32) :: vv

      vv = trim(ToLower(v))
      if (len_trim(vv) == 0) return
      if (n >= MAX_VARS) error stop "horses2probes: too many variables"

      ! Optional alias:
      if (vv == "pressure") vv = "p"

      select case (vv)
      case ("rho","u","v","w","p")
         n = n + 1
         names(n) = vv
      case default
         error stop "horses2probes: unsupported variable '"//trim(vv)//"' (use rho,u,v,w,p)"
      end select
   end subroutine PushVar


   ! ==========================================================
   ! Probes file reader: one line per probe -> x y z
   ! ==========================================================
   subroutine ReadProbesXYZ(fname, xyz, n)
      character(len=*), intent(in)  :: fname
      real(kind=RP), allocatable, intent(out) :: xyz(:,:)  ! (3,n)
      integer, intent(out) :: n
      integer :: fid, ios
      real(kind=RP) :: x,y,z

      n = 0
      open(newunit=fid, file=trim(fname), status="old", action="read")
      do
         read(fid,*,iostat=ios) x,y,z
         if (ios /= 0) exit
         n = n + 1
      end do
      rewind(fid)

      allocate(xyz(3,n))
      n = 0
      do
         read(fid,*,iostat=ios) x,y,z
         if (ios /= 0) exit
         n = n + 1
         xyz(1,n) = x
         xyz(2,n) = y
         xyz(3,n) = z
      end do
      close(fid)
   end subroutine ReadProbesXYZ


   ! ==========================================================
   ! Interpolation: FindPointInMesh (Newton) + lj() + interpolate Q
   ! ==========================================================
   subroutine InterpolateAll(mesh, probes, varNames, nVars, outVals, found)
      type(Mesh_t), intent(inout) :: mesh
      real(kind=RP), intent(in)    :: probes(:,:)      ! (3,nP)
      character(len=32), intent(in):: varNames(MAX_VARS)
      integer, intent(in)          :: nVars
      real(kind=RP), intent(out)   :: outVals(:,:)     ! (nVars,nP)
      logical, intent(out)         :: found(:)         ! (nP)

      integer :: nP, ip

      if (size(probes,1) /= 3) error stop "horses2probes: probes must have shape (3,nP)"
      nP = size(probes,2)

      if (size(outVals,1) /= nVars .or. size(outVals,2) /= nP) error stop "horses2probes: outVals wrong shape"
      if (size(found) /= nP) error stop "horses2probes: found wrong size"

      do ip = 1, nP
         call InterpolateOne(mesh, probes(:,ip), varNames, nVars, outVals(:,ip), found(ip))
      end do
   end subroutine InterpolateAll

   subroutine InterpolateOne(mesh, x, varNames, nVars, vals, ok)
      type(Mesh_t), intent(inout)  :: mesh
      real(kind=RP), intent(in)    :: x(3)
      character(len=32), intent(in):: varNames(MAX_VARS)
      integer, intent(in)          :: nVars
      real(kind=RP), intent(out)   :: vals(nVars)
      logical, intent(out)         :: ok

      integer :: eID, i, j, k, iv
      real(kind=RP) :: xi(3), w
      real(kind=RP), allocatable :: lxi(:), leta(:), lzeta(:)
      real(kind=RP), allocatable :: Qp(:)

      vals = 0.0_RP
      call FindPointInMesh(mesh, x, eID, xi, ok)
      if (.not. ok) return

      associate(e => mesh % elements(eID))

         call addNewSpectralBasis(spA, e % Nsol, mesh % nodeType)

         allocate(lxi(0:e%Nsol(1)))
         allocate(leta(0:e%Nsol(2)))
         allocate(lzeta(0:e%Nsol(3)))

         lxi   = NodalStorage(e%Nsol(1)) % lj(xi(1))
         leta  = NodalStorage(e%Nsol(2)) % lj(xi(2))
         lzeta = NodalStorage(e%Nsol(3)) % lj(xi(3))

         allocate(Qp(size(e%Q,1)))
         Qp = 0.0_RP

         do k=0,e%Nsol(3)
            do j=0,e%Nsol(2)
               do i=0,e%Nsol(1)
                  w = lxi(i)*leta(j)*lzeta(k)
                  Qp(:) = Qp(:) + e%Q(:,i,j,k) * w
               end do
            end do
         end do

         do iv = 1, nVars
            select case (trim(varNames(iv)))
            case ("rho")
               vals(iv) = Qp(IRHO)
            case ("u")
               vals(iv) = Qp(IRHOU) / Qp(IRHO)
            case ("v")
               vals(iv) = Qp(IRHOV) / Qp(IRHO)
            case ("w")
               vals(iv) = Qp(IRHOW) / Qp(IRHO)
            case ("p")
               vals(iv) = Pressure(Qp)
            end select
         end do

         deallocate(Qp, lxi, leta, lzeta)
      end associate

   end subroutine InterpolateOne


   ! ==========================================================
   ! Point location: loop over elements, Newton-invert geometry
   ! mapping x(xi) = sum x_ijk * l_i(xi)*l_j(eta)*l_k(zeta)
   ! ==========================================================
   subroutine FindPointInMesh(mesh, x, eID, xi, found)
      type(Mesh_t),  intent(in)  :: mesh
      real(kind=RP), intent(in)  :: x(3)
      integer,       intent(out) :: eID
      real(kind=RP), intent(out) :: xi(3)
      logical,       intent(out) :: found

      integer :: e
      logical :: ok

      found = .false.
      eID   = 0
      xi    = 0.0_RP

      do e = 1, mesh % no_of_elements
         call FindPointInElement(mesh % elements(e), mesh % nodeType, x, xi, ok)
         if (ok) then
            found = .true.
            eID   = e
            return
         end if
      end do
   end subroutine FindPointInMesh

   subroutine FindPointInElement(e, nodeType, x, xi, ok)
      type(Element_t), intent(in)  :: e
      integer,          intent(in) :: nodeType
      real(kind=RP),    intent(in) :: x(3)
      real(kind=RP),    intent(out):: xi(3)
      logical,          intent(out):: ok

      integer,       parameter :: MAX_ITER    = 50
      real(kind=RP), parameter :: TOL         = 1.0e-12_RP
      real(kind=RP), parameter :: INSIDE_TOL  = 1.0e-8_RP

      integer :: i, j, k, iter
      real(kind=RP), allocatable :: lxi(:), leta(:), lzeta(:)
      real(kind=RP), allocatable :: dlxi(:), dleta(:), dlzeta(:)
      real(kind=RP) :: F(3), Jac(3,3), dx(3)

      call addNewSpectralBasis(spA, e % Nmesh, nodeType)

      allocate(lxi(0:e%Nmesh(1)),  leta(0:e%Nmesh(2)),  lzeta(0:e%Nmesh(3)))
      allocate(dlxi(0:e%Nmesh(1)), dleta(0:e%Nmesh(2)), dlzeta(0:e%Nmesh(3)))

      xi = 0.0_RP

      do iter = 1, MAX_ITER
         lxi   = NodalStorage(e%Nmesh(1)) % lj(xi(1))
         leta  = NodalStorage(e%Nmesh(2)) % lj(xi(2))
         lzeta = NodalStorage(e%Nmesh(3)) % lj(xi(3))

         F = 0.0_RP
         do k=0,e%Nmesh(3) ; do j=0,e%Nmesh(2) ; do i=0,e%Nmesh(1)
            F = F + e%x(:,i,j,k) * lxi(i)*leta(j)*lzeta(k)
         end do            ; end do            ; end do
         F = F - x

         if (maxval(abs(F)) .lt. TOL) exit
         if (any(abs(xi) .ge. 2.5_RP)) exit

         dlxi   = NodalStorage(e%Nmesh(1)) % dlj(xi(1))
         dleta  = NodalStorage(e%Nmesh(2)) % dlj(xi(2))
         dlzeta = NodalStorage(e%Nmesh(3)) % dlj(xi(3))

         Jac = 0.0_RP
         do k=0,e%Nmesh(3) ; do j=0,e%Nmesh(2) ; do i=0,e%Nmesh(1)
            Jac(:,1) = Jac(:,1) + e%x(:,i,j,k) * dlxi(i)*leta(j)*lzeta(k)
            Jac(:,2) = Jac(:,2) + e%x(:,i,j,k) * lxi(i)*dleta(j)*lzeta(k)
            Jac(:,3) = Jac(:,3) + e%x(:,i,j,k) * lxi(i)*leta(j)*dlzeta(k)
         end do            ; end do            ; end do

         dx = SolveThreeEquationLinearSystem(Jac, -F)
         xi = xi + dx
      end do

      ok = all(abs(xi) .lt. 1.0_RP + INSIDE_TOL)

      deallocate(lxi, leta, lzeta, dlxi, dleta, dlzeta)

   end subroutine FindPointInElement


   ! ==========================================================
   ! Output writers
   ! ==========================================================
   subroutine WriteTecplot(fid, probes, varNames, nVars, vals, found)
      integer, intent(in) :: fid
      real(kind=RP), intent(in) :: probes(:,:)   ! (3,nP)
      character(len=32), intent(in) :: varNames(MAX_VARS)
      integer, intent(in) :: nVars
      real(kind=RP), intent(in) :: vals(:,:)     ! (nVars,nP)
      logical, intent(in) :: found(:)            ! (nP)

      integer :: i, iv, nP
      nP = size(probes,2)

      write(fid,'(A)') 'TITLE="horses2probes"'
      write(fid,'(A)', advance="no") 'VARIABLES="x","y","z"'
      do iv=1,nVars
         write(fid,'(A)', advance="no") ',"'//trim(varNames(iv))//'"'
      end do
      write(fid,'(A)') ',"found"'
      write(fid,'(A,I0,A)') 'ZONE T="PROBES" I=', nP, ', F=POINT'

      do i=1,nP
         write(fid,'(3ES24.16)', advance="no") probes(1,i), probes(2,i), probes(3,i)
         do iv=1,nVars
            write(fid,'(1X,ES24.16)', advance="no") vals(iv,i)
         end do
         write(fid,'(1X,I1)') merge(1,0,found(i))
      end do
   end subroutine WriteTecplot

   subroutine WriteTxt(fid, probes, varNames, nVars, vals, found)
      integer, intent(in) :: fid
      real(kind=RP), intent(in) :: probes(:,:)   ! (3,nP)
      character(len=32), intent(in) :: varNames(MAX_VARS)
      integer, intent(in) :: nVars
      real(kind=RP), intent(in) :: vals(:,:)     ! (nVars,nP)
      logical, intent(in) :: found(:)            ! (nP)

      integer :: i, iv, nP
      nP = size(probes,2)

      write(fid,'(A)', advance="no") "# x y z"
      do iv=1,nVars
         write(fid,'(A)', advance="no") " "//trim(varNames(iv))
      end do
      write(fid,'(A)') " found"

      do i=1,nP
         write(fid,'(3ES24.16)', advance="no") probes(1,i), probes(2,i), probes(3,i)
         do iv=1,nVars
            write(fid,'(1X,ES24.16)', advance="no") vals(iv,i)
         end do
         write(fid,'(1X,I1)') merge(1,0,found(i))
      end do
   end subroutine WriteTxt

#endif
end module Solution2ProbesModule
