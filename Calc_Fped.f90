! This code is adapted from meuw.f from PEDIG
! by Jisca Huisman in 2024
! - translation from Fortran77 to Fortran2003 standard
! - removed dependencies on other source files in the PEDIG program
! - translation of messages from French to English


! ==============================================================================
! ==============================================================================

module utils
  implicit none
  
contains
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! round number x to d significant digits
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  integer function roundit(x,d)
    double precision, intent(IN) :: x
    integer, intent(IN) :: d
    double precision :: z
    integer :: i
    
    z = x
    i = 0
    do
      i = i+1
      z = z/10.0
      if (z < 10**d)  exit
    enddo
  
    roundit = NINT(z) * 10**i     
  end function roundit
  
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! determine the number of rows in a file
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  integer function FileNumRow(FileName)
    character(len=*), intent(IN) :: FileName
    integer :: nrow, i, maxRow, ios
    character(len=42) :: dumC

    maxRow = 1e8  ! fail safe
    nrow = 0
    open(unit=102, file=trim(FileName), status="old", action='read')
      do i=1, maxRow
        read(102,*,IOSTAT=ios) dumC
        if (ios < 0) then
          exit  ! end of file
        else
          nrow = nrow +1  
        end if
      enddo
    close(102)
    if (nrow >= maxRow) print *, 'WARNING: '//trim(FileName)//' EXCEEDING MAXIMUM NUMBER OF ROWS!'
    FileNumRow = nrow

  end function FileNumRow 
  
end module utils

! ==============================================================================


! ==============================================================================
! ==============================================================================
!  meuw.f in PEDIG
! ==============================================================================
! This code is by Didier Boichard
! and is originallly part of PEDIG, described in 
! Pedig: a fortran package for pedigree analysis suited for large populations.
! 7th World Congress on Genetics Applied to Livestock Production, August 19-23, 2002, Montpellier, France
! https://hal.inrae.fr/hal-02833573/document

! downloaded on 2024-06-24 from 
! https://gabi.jouy.hub.inrae.fr/services-ressources/logiciels-et-outils/pedig 
! publicly available without any license

! The algorithm is described in
! Meuwissen & Luo, 1992, Computing inbreeding coefficients in large populations, GSE

module Fped_module
  use utils, only: roundit
  implicit none
  
  contains 

function F_meuwissen(ped, N)  
  implicit none
  integer, intent(IN) :: N  ! number of individuals
  double precision :: F_meuwissen(1:N), F(0:N)
  integer, intent(IN) :: ped(2,N) 
  integer :: i, id, is, j, k, kd, ks  ! individual indices
  integer :: point(N)   ! ??
  integer :: np, npp(200)  ! effective number of ancestors (not sure why the 200)
  integer :: ninbr       ! number of individuals with non-zero inbreeding coefficient
  integer :: print_chunk, t  ! for printing progress to console
  double precision :: fi  ! temp for F(i)
  double precision :: l(N), d(N)  ! ??
  double precision :: r    ! relatedness?
  
  ! Meuwissen & Luo 1992:
  ! A = LDL' where
  ! L: lower triangular matrix containing the fraction of the genes that animals 
  !    derive from their ancestors
  ! D: diagonal matrix containing the within family additive genetic variances of animals
  ! animals are ordered so that parents precede offspring (ancestor j < focal i)
  ! it follows that (Quaas, 1976)
  ! A_{ii} = \sum_{j=i}^i L_{ij}^2 D_{jj}
  
  ! genes shared with ancestor j = half of what dam shares + half of what sire shares
  ! L_{ij} = (L_{s_{i}j} + L_{d{i}j})/2 
  ! within-family V_A:
  ! D_{jj} = 0.5 - (F_{s_j} + F_{d_j})/4
    
  npp = 0
  ninbr = 0 
  F(0) = -1.d0
  l = 0.d0  
  do i=1,N
    if (i < Ped(1,i) .or. i < Ped(2,i)) then
      print *, 'Problems with pedigree coding: parents should be before offspring!'
      print *, i, Ped(:,i)
      stop
    endif
    point(i) = 0
  enddo
  
  ! print progress to console at approx every 5%
  print_chunk = roundit(n/20D0,2)  
  t = 1  
  
  do i=1,n
    if (mod(i,print_chunk)==0) then
      print *, i, '   ', t*5,'%'
      t = t+1
    endif
    is=ped(1,i)
    id=ped(2,i)
    ! calculate within-family V_A
    d(i)=.5d0 - .25d0*(f(is)+f(id))
    ! if either parent is unknown, individual is by definition not inbred
    if (is == 0 .or. id == 0) then
      f(i) = 0.d0
    else
      np = 0
      fi=-1.d0
      l(i)=1.d0  
      j=i
      do while(j /= 0)
        k=j
        r=.5d0 * l(k)
        ks = maxval(ped(:,k))
        kd = minval(ped(:,k))    
        if (ks > 0) then
          l(ks)=l(ks) + r
          do while(point(k) > ks)
            k=point(k)
          end do
          if (ks /= point(k)) then
            point(ks)=point(k)
            point(k)=ks
          end if
          if (kd > 0) then
            l(kd)=l(kd) + r
            do while(point(k) > kd)
              k=point(k)
            end do
            if (kd /= point(k)) then
              point(kd)=point(k)
              point(k)=kd
            end if
          end if
        end if
        
        fi=fi + l(j)*l(j)*d(j) 

        l(j)=0.d0
        k=j
        j=point(j)
        point(k)=0
        np=np+1
      end do   ! end of while loop 
      
      f(i)=fi
      if (fi > 0.000001d0) ninbr=ninbr + 1
      if (np > 200) np=200   ! ??
      if (np == 0) stop 'strange error...'
      npp(np)=npp(np)+1
    end if  ! end of: both parents are known
  enddo
  
  F_meuwissen = F(1:N)  ! exclude F(0) from output
    
!    print *, 'Effective number of ancestors:'
!    do i=1,200
!      if (npp(i) > 0)  write(*, '(i6, i10)') i, npp(i)
!    enddo

  write(*, '("Number of inbred animals: ", i8)')  ninbr
  
end function F_meuwissen

end module Fped_module




! ==============================================================================
! ==============================================================================
! - sort pedigree to ensure parents are always listed before their offspring
! - read pedigree from file
! ==============================================================================

module pedigree_tools
  implicit none
 ! public :: sort_pedigree
 ! private :: QsortC, Partition, get_generations

  integer, parameter :: nchar_ID = 40

contains

  !=============================================================================
  ! Recursive Fortran 95 quicksort routine
  ! sorts real numbers into ascending numerical order
  ! Author: Juli Rew, SCD Consulting (juliana@ucar.edu), 9/03
  ! Based on algorithm from Cormen et al., Introduction to Algorithms,1997
  ! Made F conformant by Walt Brainerd

  ! Adapted by J Huisman (jisca.huisman@gmail.com) to output rank, to
  ! enable sorting of parallel vectors, and changed to decreasing rather 
  ! than increasing order
  recursive subroutine QsortC(A, Rank)
    double precision, intent(in out), dimension(:) :: A
    integer, intent(in out), dimension(:) :: Rank
    integer :: iq

    if(size(A) > 1) then
     call Partition(A, iq, Rank)
     call QsortC(A(:iq-1), Rank(:iq-1))
     call QsortC(A(iq:), Rank(iq:))
    endif
  end subroutine QsortC

  subroutine Partition(A, marker, Rank)
    double precision, intent(in out), dimension(:) :: A
    integer, intent(in out), dimension(:) :: Rank
    integer, intent(out) :: marker
    integer :: i, j, TmpI
    double precision :: temp
    double precision :: x      ! pivot point
    x = A(1)
    i= 0
    j= size(A) + 1
    do
     j = j-1
     do
      if (j < 1) exit
      if (A(j) <= x) exit
      j = j-1
     end do
     i = i+1
     do
      if (i >= size(A)) exit
      if (A(i) >= x) exit
      i = i+1
     end do
     if (i < j) then
      ! exchange A(i) and A(j)
      temp = A(i)
      A(i) = A(j)
      A(j) = temp 
      
      TmpI = Rank(i) 
      Rank(i) = Rank(j)
      Rank(j) = TmpI
     elseif (i == j) then
      marker = i+1
      return
     else
      marker = i
      return
     endif
    end do

  end subroutine Partition
  !=============================================================================


  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! get generation numbers
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
  function get_generations(Ped, N)  result(Gen)
    integer, intent(IN) :: N
    integer, intent(IN), allocatable :: Ped(:,:)  ! 2,N
    integer, allocatable :: Gen(:), Gen_tmp(:)
    integer :: i, g
    integer, parameter :: gen_max = 200   ! maximum pedigree depth
    
    allocate(Gen(0:N))
    Gen = gen_max
    Gen(0) = 0
    
    do i=1,N
      if (all( Ped(:,i)==0 ))   Gen(i) = 0  ! founder
    enddo
    
    do g=0, gen_max
      do i=1,N
        if (Gen(i) < gen_max)  cycle  ! generation number already determined
        if (all( Gen( Ped(:,i) ) <= g ))  Gen(i) = g+1
      enddo
      if (all(Gen < gen_max))  exit  ! all done
    enddo
    
    ! remove index 0
    allocate(Gen_tmp(N))
    Gen_tmp = Gen(1:N)
    call move_alloc(Gen_tmp, Gen)  ! from, to
    
  end function get_generations
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! sort pedigree & ID vector
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
  subroutine sort_pedigree(Ped, IDs, N)
    integer, intent(IN) :: N
    integer, intent(INOUT), allocatable :: Ped(:,:)   ! 2,N
    character(len=nchar_ID), intent(INOUT), allocatable :: IDs(:)
    integer :: Rank(N), i, k, j
    double precision, allocatable :: Gen_dbl(:)
    integer, allocatable :: Gen(:), Ped_new(:,:)   ! 2,N
    character(len=nchar_ID), allocatable :: IDs_new(:)
    
    Gen = get_generations(Ped, N)

    Rank = (/ (i, i=1, N, 1) /)  ! vector to be sorted
    allocate(Gen_dbl(N))
    Gen_dbl = REAL(Gen, 8)        ! sorting algorithm uses double precision input
    call QsortC(Gen_dbl, Rank)
    deallocate(Gen_dbl)
    
    ! sort individual names in new order
    allocate(IDs_new(N))
    do i=1,N
      IDs_new(i) = Ids(Rank(i))
    enddo 
    
    ! renumber parents to match new ID order
    do i=1,N
      do k=1,2
        if (Ped(k,i) == 0) then
          Ped(k,i) = 0
        else
          do j=1,N
            if (IDs(Ped(k,i)) == IDs_new(j)) then
              Ped(k,i) = j
              exit
            endif         
          enddo  ! j
        endif
      enddo
    enddo

    ! reorder rows in pedigree
    allocate(Ped_new(2,N))
    do i=1,N
      Ped_new(:,i) = Ped(:,Rank(i))
    enddo
      
    ! move new pedigree & ids to input/output
    call move_alloc(Ped_new, Ped)  ! from, to
    call move_alloc(Ids_new, IDs)
  
  end subroutine sort_pedigree  
  
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! read pedigree from text file
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
  subroutine read_pedigree(Ped, IDs, FileName, hasHeader) 
    use utils, only: FileNumRow
    integer, allocatable, intent(OUT) :: Ped(:,:)
    character(len=nchar_ID), allocatable, intent(OUT) :: IDs(:)
    character(len=*), intent(IN) :: FileName
    logical, intent(IN), optional :: hasHeader   ! default=TRUE
    logical :: file_exists, do_read_header
    character(len=nchar_ID), allocatable :: NamePed(:,:), IDs_extra(:), IDs_tmp(:)
    character(len=nchar_ID) :: tmpC(3)
    integer :: N, i, j, k, ios, N_extra
    
    inquire(file=trim(FileName), exist = file_exists)
    if (.not. file_exists)  stop 'File '//trim(FileName)//' not found!'
    
    N = FileNumRow(FileName) 
    
    if (present(hasHeader)) then
      do_read_header = hasHeader
    else
      do_read_header = .TRUE.
    endif
    if (do_read_header)  N = N-1
    
    allocate(IDs(N))  
    allocate(NamePed(2,N))
    
    ! read in names from file
    open(unit=103, file=trim(FileName), status="old")
      if (do_read_header)  read(103,*)                             
      do i=1,N
        read(103, *,IOSTAT=ios)  tmpC
!        call IOstat_handler(ios, i, FileName)   in sqa_fileIO.f90
        IDs(i) = tmpC(1)
        NamePed(:,i) = tmpC(2:3) 
      enddo
    close(103)
    
    ! add extra entries for parents that don't have own entry
    N_extra = 0
    allocate(IDs_extra(N))
    do i=1,N
      do k=1,2
        if (NamePed(k,i) == 'NA') cycle
        if (.not. any(IDs == NamePed(k,i)) .and. .not. any(IDs_extra == NamePed(k,i))) then
          N_extra = N_extra + 1
          IDs_extra(N_extra) = NamePed(k,i)
        endif
      enddo
    enddo
    
    print *, 'N_extra: ', N_extra
    
    if (N_extra > 0) then
      allocate(IDs_tmp(N + N_extra))
      ! putting parents without own entry first (by definition founders) makes next step too confusing
      IDs_tmp(1:N) = IDs
      IDs_tmp((N+1):(N+N_extra)) = IDs_extra(1:N_extra)
      call move_alloc(IDs_tmp, IDs)  ! from, to
    endif


    ! transform names to numbers
    ! this is the slow part, but alternatives (e.g. using WHERE) are slower.
    allocate(Ped(2,N+N_extra))
    Ped = 0
    do i=1,N
      do k=1,2
        if (NamePed(k,i) == 'NA') then
          Ped(k,i) = 0
        else if (any(IDs == NamePed(k,i))) then
          do j=1,N
            if (NamePed(k,i) == IDs(j)) then
              Ped(k,i) = j
              exit
            endif         
          enddo  ! j
        else
          print *, 'read_pedigree error: ID not found! ', i, k, NamePed(:,i)
          stop
        endif
      enddo  ! k
    enddo  ! i
    
    deallocate(NamePed)
    deallocate(IDs_extra)
  
  end subroutine read_pedigree
  

end module pedigree_tools




! ==============================================================================
! ==============================================================================
!             main program
! ==============================================================================

program calc_Fped
  use pedigree_tools, only: read_pedigree, sort_pedigree, nchar_ID
  use Fped_module
  implicit none
  
  character(len=2000) :: pedigree_file, Fped_file, sorted_pedigree_file
  integer, allocatable :: pedigree(:,:)
  double precision, allocatable :: Fped(:)
  character(len=nchar_ID), allocatable :: IDs(:), IDz(:)
  integer :: N, i
  logical :: pedigree_is_sorted, return_sorted_pedigree
  
  ! set defaults
  pedigree_file = 'pedigree.txt'   ! TODO: read from command line
  Fped_file = 'Fped.txt'
  return_sorted_pedigree = .FALSE.
  
  ! read arguments from command line
  call read_args()
  
  ! read pedigree file
  print *, 'reading pedigree file...'
  call read_pedigree(pedigree, IDs, pedigree_file)
  N = size(pedigree, dim=2)
  print *, 'read ', N , ' individuals'
  
  ! sort pedigree (parents always before their offspring) if needed
  pedigree_is_sorted = .TRUE.
  do i=1,N
    if (pedigree(1,i) > i .or. pedigree(2,i) > i) then
      pedigree_is_sorted = .FALSE.
      exit
    endif
  enddo
  if (.not. pedigree_is_sorted) then
    print *, 'sorting pedigree ...'
    call sort_pedigree(pedigree, IDs, N)
  endif
  
  if (return_sorted_pedigree) then
    allocate(IDz(0:N))
    IDz(0) = 'NA'
    IDz(1:N) = IDs
    sorted_pedigree_file = pedigree_file(1:(LEN_TRIM(pedigree_file)-4))//'_sorted.txt'    
    open(unit=101, file=trim(sorted_pedigree_file), status='unknown')
      write(101, '(3a20)') 'ID', 'dam', 'sire'
      do i=1,N
        write(101, '(3a20, i6)')  IDs(i), IDz(pedigree(:,i))
      enddo
    close(101)
  
  endif
 
  
  ! calculate inbreeding coefficients
  print *, 'calculating inbreeding coefficients...'
  Fped = F_meuwissen(pedigree, N)  
   
  ! write to file
  open(unit=101, file=trim(Fped_file), status='unknown')
    write(101, '(a40, 2x, a10)') 'ID', 'F.ped'
    do i=1,N
      write(101, '(a40, 2x, f10.8)')  IDs(i), Fped(i)   
    enddo
  close(101)
  
  
contains
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! read command line arguments
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine read_args()
    integer :: nArg, i
    character(len=32) :: arg
    
    nArg = command_argument_count()
    i = 0
    do while (i < nArg)
      i = i+1
      call get_command_argument(i, arg)
      
      select case (arg)
        ! case ('-h', '--help')
          ! call print_help()
          ! stop
        case ('--pedigree')
          i = i+1
          call get_command_argument(i, pedigree_file)
        case ('--out')
          i = i+1
          call get_command_argument(i, Fped_file)
        case ('--return-sorted-pedigree')
          return_sorted_pedigree = .TRUE.
          
        case default
          print *, ''
          print *, 'Unrecognised command-line argument: ', arg 
          stop
      end select
    end do

  end subroutine read_args

end program calc_Fped
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  
! ==============================================================================
! ==============================================================================