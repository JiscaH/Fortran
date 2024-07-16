! This code is adapted from meuw.f from PEDIG
! by Jisca Huisman in 2024
! - translation from Fortran77 to Fortran2003 standard
! - removed dependencies on other source files in the PEDIG program
! - translation of messages from French to English


! ==============================================================================
! ==============================================================================

module utils
  implicit none
  
  logical :: quiet
  
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


  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! print timestamp
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine timestamp(add_blank)
    ! NOTE: print *, & write(*,*) have leading blank, write(*,'(...)') does not
    logical, intent(IN), OPTIONAL :: add_blank
    logical :: do_blank
    integer :: date_time_values(8)
    
    if(present(add_blank))then
      do_blank=add_blank
    else
      do_blank=.FALSE.
    endif
       
    call date_and_time(VALUES=date_time_values)
    write(*,'(i2.2,":",i2.2,":",i2.2, 1X)', advance='no') date_time_values(5:7)
    if (do_blank) write(*, '(1X)', advance='no')
     
  end subroutine timestamp
    
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! print text, preceded by timestamp
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
  subroutine printt(text, int)
    character(len=*), intent(IN), optional :: text
    integer, intent(IN), optional :: int   ! e.g. counter in a loop
    character(len=10) :: intchar
    
    call timestamp()
    if (present(text) .and. present(int)) then
      write(intchar,'(i10)') int
      print *, text//adjustl(trim(intchar))
    else if (present(text)) then
      print *, text
    else if (present(int)) then
      print *, int
    endif
  
  end subroutine printt  
  
  
end module utils

! ==============================================================================


! ==============================================================================
! ==============================================================================
!  meuw.f in PEDIG
! ==============================================================================
! 
! This method to calculate pedigree inbreeding coefficients is described in
! Meuwissen & Luo, 1992, Computing inbreeding coefficients in large populations, GSE
! and a version of this Fortran algorithm was provided as appendix in that paper

! further explanation of the mathematical terms can be found in 
! Henderson, 1976, A Simple Method for Computing the Inverse of a Numerator Relationship 
!    Matrix Used in Prediction of Breeding Values

!
! This current version is an adaptation of meuw.f in PEDIG by Didier Boichard, described in 
! Pedig: a fortran package for pedigree analysis suited for large populations.
! 7th World Congress on Genetics Applied to Livestock Production, August 19-23, 2002, Montpellier, France
! https://hal.inrae.fr/hal-02833573/document

! downloaded on 2024-06-24 from 
! https://gabi.jouy.hub.inrae.fr/services-ressources/logiciels-et-outils/pedig 
! publicly available without any license

! Adapted by Jisca Huisman, 2024.

! NOTE: this algorithm (probably) doesn't work with selfing


module Fped_module
  use utils, only: roundit, quiet
  implicit none
  
  contains 

function F_meuwissen(ped, N)  
  implicit none
  integer, intent(IN) :: N  ! number of individuals
  double precision :: F_meuwissen(1:N)  ! output
  double precision :: F(0:N)  ! working vector for inbreeding coefficient, starting from 0
  integer, intent(IN) :: ped(2,N)   ! pedigree matrix
  integer :: i, id, is, j, k, kd, ks  ! individual indices
  integer :: point(N)   ! link list with next oldest ancestor; 0 if i is oldest ancestor
  integer :: np(N)     ! effective number of ancestors (=number of paths)
  double precision :: fi  ! temp for F(i)
  double precision :: L(N)  ! L(j) = element ij of matrix L, when animal i is evaluated
                            ! = fraction of genes of animal i that derive from j
  double precision :: D(N)  ! within family variance of animal i. 
      ! = similar to the F_A (ancestor inbreeding coefficient) term in the path method
  integer :: print_chunk, t  ! for printing progress to console


  ! animals are ordered so that parents precede offspring (ancestor j < focal i)
  do i=1,N
    if (i == Ped(1,i) .or. i == Ped(2,i)) then
      print *, 'Problems with pedigree: individual is its own parent!'
      print *, i, Ped(:,i)
      stop
    else if (i < Ped(1,i) .or. i < Ped(2,i)) then
      print *, 'Problems with pedigree coding: parents should be before offspring!'
      print *, i, Ped(:,i)
      stop
    endif
  enddo
   
  ! From Meuwissen & Luo, 1992:
  ! additive relationship matrix A = LDL' where
  ! L: lower triangular matrix containing the fraction of the genes that animals 
  !    derive from their ancestors
  ! D: diagonal matrix containing the within family additive genetic variances of animals
  
  ! it follows that (Quaas, 1976)
  ! A_{ii} = \sum_{j=i}^i L_{ij}^2 D_{jj}
  
  ! genes shared with ancestor j = half of what dam shares + half of what sire shares
  ! L_{ij} = (L_{s_{i}j} + L_{d{i}j})/2 
  ! within-family V_A:
  ! D_{jj} = 0.5  - (F_{s_j} + F_{d_j})/4   if both parents known
  !        = 0.75 - F_{k_j}/4               if one parent known
  !        = 1.0                            if no parents known  
  F(0) = -1.d0  ! ensures correct within-family variances for unknown parents
  ! (also more accurate than calculating GRM diagonal F+1 and then subtracting 1, wrt rounding errors)
     
  ! print progress to console at approx every 5%
  print_chunk = roundit(n/20D0,2)  
  t = 1  
  
  ! initiate linked list vector
  point = 0
  ! initiate vector with number of paths per individual
  np = 0
  
  do i=1,n
   ! print progress to console
   if (.not. quiet .and. mod(i,print_chunk)==0) then
     print *, i, '   ', t*5,'%'
     t = t+1
   endif
   
    ! get sire & dam of focal individual i
    is=ped(1,i)
    id=ped(2,i)
    ! reset vector l: this is the i'th row of the L matrix 
    l = 0.d0   
    ! calculate within-family V_A
    D(i) = .5d0 - .25d0*(f(is)+f(id))  
    
    ! if either parent is unknown, individual is by definition not inbred
    if (is == 0 .or. id == 0) then
      f(i) = 0.d0
    
    ! if individual is full sib of preceding individual, they have identical inbreeding coefficients
    else if (ped(1,i) == ped(1,i-1) .and. ped(2,i) == ped(2,i-1)) then
      f(i) = f(i-1)
    
    else
      fi=-1.d0   ! inbreeding coefficient of individual i (working number)
      l(i)=1.d0  ! fraction of genes that i derives from itself (starting point)
      j=i        ! start at diagonal
      do while(j /= 0)   ! loop over all earlier individuals in the pedigree
        k=j              ! k is a temporary variable
        ! use oldest parent of k as 'sire', and youngest as 'dam'
        ks = maxval(ped(:,k))   
        kd = minval(ped(:,k))
        ! increase fraction of genes that i derives from ks and kd by 1/2*[fraction derived from k] 
        if (ks > 0)  l(ks) = l(ks) + 0.5d0 * l(k)
        if (kd > 0)  l(kd) = l(kd) + 0.5d0 * l(k)
        
        ! add contribution L*D*L of animal j to inbreeding coefficient of i
        fi=fi + l(j)*l(j)*d(j) 
        
        ! update link list
        if (ks > 0) then    ! find slot in link list for sire
          do while(point(k) > ks)   ! stop if next ancestor is older than sire, or is sire.
            k=point(k)
          end do
          if (ks /= point(k)) then  ! include sire in link list, if not already in there
            point(ks)=point(k)
            point(k)=ks
          end if
          if (kd > 0) then   ! do same for dam. 
            ! no need to reinitialise variable k, because initially point(k)==ks, and ks > kd
            do while(point(k) > kd)             
              k=point(k)
            end do
            if (kd /= point(k)) then
              point(kd)=point(k)
              point(k)=kd
            end if
          end if
        end if
    
        ! prep for next loop
        l(j)=0.d0    ! clear l(j)
        k=j          ! store 'old' j in temporary variable k
        j=point(j)   ! 'new' j is next oldest animal in list
        point(k)=0   ! clear point('old' j)
        np(i)=np(i)+1      ! counter of number of paths/ancestors
      end do   ! end of while loop 
      
      F(i)=fi   ! store result of iterations over paths for this individual in vector
      if (np(i) == 0) stop 'strange error...'   
    end if  ! end of: both parents are known
  enddo
  
  F_meuwissen = F(1:N)  ! exclude F(0) from output
  
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
    ! sorting algorithm uses double precision input
    allocate(Gen_dbl(N))
    ! add original order as fractional element, to keep original order within each generation
    Gen_dbl = REAL(Gen, 8) + Rank/REAL(N,8)    
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
            if (Ped(k,i) == Rank(j)) then
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
      
    ! move new pedigree & new ids to input/output arrays
    call move_alloc(Ped_new, Ped)  ! from, to
    call move_alloc(Ids_new, IDs)
  
  end subroutine sort_pedigree  
  
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! read pedigree from text file
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
  subroutine read_pedigree(Ped, IDs, FileName, hasHeader) 
    use utils, only: FileNumRow, timestamp, printt, quiet
    integer, allocatable, intent(OUT) :: Ped(:,:)
    character(len=nchar_ID), allocatable, intent(OUT) :: IDs(:)
    character(len=*), intent(IN) :: FileName
    logical, intent(IN), optional :: hasHeader   ! default=TRUE
    logical :: file_exists, do_read_header, found
    character(len=nchar_ID), allocatable :: NamePed(:,:), IDs_tmp(:)
    integer, allocatable :: Ped_tmp(:,:)
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
      
    
    ! read in names from file
    allocate(NamePed(3,N))
    open(unit=103, file=trim(FileName), status="old")
      if (do_read_header)  read(103,*)                             
      do i=1,N
        read(103, *,IOSTAT=ios)  NamePed(:,i)  ! id - parent1 (dam) - parent2 (sire)
        if (ios > 0) then
          print *, "ERROR: Wrong input in file "//trim(FileName)//" on line ", i
          stop
        endif
      enddo
    close(103)

    ! transform names to numbers
    ! this is the slow part, but alternatives (e.g. using WHERE) are slower.
    if (.not. quiet)  call printt('transforming pedigree from character to numeric ...')
    allocate(IDs(-(2*N):N))
    IDs = 'NA'
    IDs(1:N) = NamePed(1,:)
    
    allocate(Ped(2,N))
    Ped = 0
    N_extra = 0   ! number of parents without own entry
    do i=1,N
      do k=1,2
       if (NamePed(k+1,i) == 'NA') then
         Ped(k,i) = 0
       else if (NamePed(k+1,i) == NamePed(1,i)) then
         print *, 'individual is its own parent! ', i
         stop
       else 
          found = .FALSE.
          do j=-N_extra,N
            if (NamePed(k+1,i) == IDs(j)) then
              Ped(k,i) = j
              found = .TRUE.
              exit
            endif         
          enddo  ! j
         if (.not. found) then   ! parent does not have own entry pedigree
           if (.not. quiet)  print *, 'extra : ', N_extra, k, NamePed(:,i) 
           N_extra = N_extra +1
           IDs(-N_extra) = NamePed(k+1,i)
           Ped(k,i) = -N_extra
         endif
       endif
      enddo  ! k
    enddo  ! i
        
    if (N_extra > 0) then
      if (.not. quiet)  print *, 'Adding parents without own entry, N= ', N_extra, ' ...'
      ! renumber pedigree
      do k=1,2
        WHERE(ped(k,:) > 0)  Ped(k,:) = Ped(k,:) + N_extra
        WHERE(ped(k,:) < 0)  Ped(k,:) = Ped(k,:) + N_extra +1  ! hop over 0 
      enddo      
      ! add rows to pedigree
      allocate(ped_tmp(2,N+N_extra))
      Ped_tmp(:,1:N_extra) = 0   ! parents that had no entry, have no parents themselves
      Ped_tmp(:,(N_extra+1):(N+N_extra)) = Ped
      call move_alloc(Ped_tmp, Ped)  ! from, to; deallocates Ped_tmp     
      
      ! add to ID vector
      allocate(IDs_tmp(N+N_extra))
      IDs_tmp(1:N_extra) = IDs(-N_extra : -1)
      IDs_tmp((N_extra+1):(N+N_extra)) = IDs(1:N)
      call move_alloc(IDs_tmp, IDs)  ! from, to; deallocates IDs_tmp
            
      ! update total number of individuals     
      N = N + N_extra
    endif
    
    deallocate(NamePed)
    
  
  end subroutine read_pedigree
  

end module pedigree_tools




! ==============================================================================
! ==============================================================================
!             main program
! ==============================================================================

program calc_Fped
  use pedigree_tools, only: read_pedigree, sort_pedigree, nchar_ID
  use Fped_module
  use utils, only: timestamp, printt, quiet
  implicit none
  
  character(len=2000) :: pedigree_file, Fped_file, sorted_pedigree_file
  integer, allocatable :: pedigree(:,:)
  double precision, allocatable :: Fped(:)
  character(len=nchar_ID), allocatable :: IDs(:), IDz(:)
  integer :: N, i
  logical :: pedigree_is_sorted, return_sorted_pedigree
  
  ! set defaults
  pedigree_file = 'pedigree.txt'  
  Fped_file = 'Fped.txt'
  return_sorted_pedigree = .FALSE.
  quiet = .FALSE.
  
  ! read arguments from command line
  call read_args()
  
  ! read pedigree file
  if (.not. quiet)  call printt('reading pedigree file...')
  call read_pedigree(pedigree, IDs, pedigree_file)
  N = size(pedigree, dim=2)
  call timestamp()
  if (.not. quiet)  print *, 'read ', N , ' individuals'
  
  ! sort pedigree (parents always before their offspring) if needed
  pedigree_is_sorted = .TRUE.
  do i=1,N
    if (pedigree(1,i) > i .or. pedigree(2,i) > i) then
      pedigree_is_sorted = .FALSE.
      exit
    endif
  enddo
  if (.not. pedigree_is_sorted) then
    if (.not. quiet)  call printt('sorting pedigree ...')
    call sort_pedigree(pedigree, IDs, N)
  endif
  
  if (return_sorted_pedigree) then
    allocate(IDz(0:N))
    IDz(0) = 'NA'
    IDz(1:N) = IDs
    sorted_pedigree_file = pedigree_file(1:(LEN_TRIM(pedigree_file)-4))//'_sorted.txt'    
    open(unit=101, file=trim(sorted_pedigree_file), status='unknown')
      write(101, '(3a20)') 'id', 'dam', 'sire'
      do i=1,N
        write(101, '(3a20, i6)')  IDs(i), IDz(pedigree(:,i))
      enddo
    close(101) 
  endif
 
  
  ! calculate inbreeding coefficients
  if (.not. quiet)  call printt('calculating inbreeding coefficients...')
  Fped = F_meuwissen(pedigree, N)  
  call timestamp()
  write(*, '("Number of inbred animals: ", i8)')  count(Fped > 0.000001d0)
  
   
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
        case ('-h', '--help')
          call print_help()
          stop
        case ('--pedigree')
          i = i+1
          call get_command_argument(i, pedigree_file)
        case ('--out')
          i = i+1
          call get_command_argument(i, Fped_file)
        case ('--return-sorted-pedigree')
          return_sorted_pedigree = .TRUE.
        case ('--quiet')
          quiet = .TRUE.
          
        case default
          print *, ''
          print *, 'Unrecognised command-line argument: ', arg 
          stop
      end select
    end do

  end subroutine read_args
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! print help info
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine print_help()
    print *, 'still to be writen'   ! TODO
    stop
  
  
  end subroutine print_help
  

end program calc_Fped
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  
! ==============================================================================
! ==============================================================================