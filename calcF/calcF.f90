! calc genomic inbreeding coefficients
! using equations in Zhang, Goudet & Weir. Rank-invariant estimation of inbreeding coefficients
! Heredity (2022) 128:1–10; https://doi.org/10.1038/s41437-021-00471-4

! calculate various estimators:
! - Funiu  = Fhat3 in plink --ibc
! - Funiw  = F_uni_w = as Fhat3, but weighted (w) ratio of averages over loci rather
!            than unweighted (u) average of ratios
! - FAS = Allele Sharing estimator
! maybe add later:
! - psi    = average coancestry of individual j with other members of
! the study sample   (FAS = Funiw + 2psi = Fhat1 + 4psi)
! - FSTD = Fhat1 in plink --ibc

! ==============================================================================
! ==============================================================================

module Global
  use sqa_fileIO, ONLY: ishort, nchar_ID
  implicit none
  
  integer :: Ng, nSnp, Np, Nc
  integer(kind=ishort), allocatable :: Geno(:,:)    ! genotype minor allele counts (0/1/2)
  integer, allocatable :: pedigree(:,:)
  character(len=nchar_ID), allocatable, dimension(:) :: IDp, IDg, SNP_names, IDc
  double precision, allocatable, dimension(:) :: p, F_uni_u, F_uni_w, R_parents
  logical :: quiet


!=================================================
contains 
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! help
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine print_help()
    print '(a)',  'Calculate inbreeding coefficient from pedigree and/or SNP data'
    print '(a)',  'if both --pedigree & --geno specified, genomic relatedness between parents is calculated too'
    print '(a, /)', 'command-line options:'
    print '(a)',    '  --help               print usage information and exit'
    print '(a)',    '  --pedigree <filename> file with pedigree, columns id-parent1-parent2. Column names ignored.'
    print '(a)',    '  --geno <filename>    file with genotype data, file extension will be added based on',&
                    '                         --informat. Default: Geno'
    print '(a)',    '  --informat <x>       SEQ: no header, 0/1/2/-9, IDs in column 1; .txt (default)', &  
                    '                       PED: no header, 11/12/22/00, IDs in column 2 of 6 ', &
                    '                        non-SNP columns; .ped (+.map)', &
                    '                       RAW: header, 0/1/2/NA, IDs in column 2 of 6 non-SNP columns; .raw', &
                    '                       LMT: no header, 0/1/2 without spacing, IDs in separate file;', &
                    '                        .geno + .id'    
    print '(a)',    '  --maf <x>            filter on minimum allele frequency (default: 0.01)'   
    print '(a)',    '  --freq, --af <filename>      optional input file with allele frequencies; only relevant',&
                    '                        in combination with --min_prob. Either 1 column and no header,',&
                    '                        or multiple columns with a column MAF, AF, or Frequency',&
                    '                        E.g. output from plink --freq (cannot contain NAs!).'                     
    print '(a)',    '  --out <filename>     output file name. Default: inbreeding_coefficients.txt'
    print '(a)',    '  --return-sorted-pedigree  write pedigree to file after sorting parents before offspring'
    print '(a)',    '  --quiet              suppress all messages'
  end subroutine print_help
  
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! get allele frequencies, from separate file or estimated from data
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function getAF(FileName)
    use sqa_fileIO, ONLY: readAF
    
    character(len=*), intent(IN), optional :: FileName
    double precision :: getAF(nSnp)
    double precision, allocatable :: AF_tmp(:)
    logical :: FromFile
       
    if (.not. present(FileName)) then
      FromFile = .FALSE.
    else if (FileName == 'NoFile') then
      FromFile = .FALSE.
    else 
      FromFile = .TRUE.
    endif
    
    if (FromFile) then  
      
      if (.not. quiet) print *, "Reading allele frequencies in "//trim(FileName)//" ..."   
      AF_tmp = readAF(trim(FileName))
      if (SIZE(AF_tmp) /= nSnp) then
        stop "MAF file "//trim(FileName)//" has different number of SNPs than genotype file!"
      else
        getAF = AF_tmp
        deallocate(AF_tmp)
      endif
      
    else
    
      getAF = calcAF()   
    
    endif

  end function getAF
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! calculate allele frequencies (SNPs in D2)
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  pure function calcAF()
    use sqa_fileIO, ONLY: ishort, ilong
    
    double precision :: calcAF(nSnp)
    integer :: l
    
    calcAF = 1D0
    do l=1,nSnp
      if (ALL(Geno(l,:)==-1)) cycle
      calcAF(l) = dble(SUM(int(Geno(l,1:),kind=ilong), MASK=Geno(l,1:)/=-1))/(COUNT(Geno(l,1:)/=-1)*2)
    enddo
  end function calcAF
  
  
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


end module Global



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
  use Global, only: roundit, quiet
  implicit none
  
  contains 

function F_meuwissen(ped)   result(F)
  integer, intent(IN), allocatable :: ped(:,:)   ! 2,N pedigree matrix 
  double precision, allocatable :: F(:), L(:), D(:)
  integer, allocatable :: point(:), np(:)
  integer :: N  ! number of individuals
  integer :: i, id, is, j, k, kd, ks  ! individual indices  
  double precision :: fi  ! temp for F(i)
  integer :: print_chunk, t  ! for printing progress to console

  N = size(Ped,2)
!  allocate(F_meuwissen(N))  ! output
  allocate(F(0:N))   ! working vector
  allocate(point(N))   ! link list with next oldest ancestor; 0 if i is oldest ancestor
  allocate(np(N))     ! effective number of ancestors (=number of paths)
  allocate(L(N))  ! L(j) = element ij of matrix L, when animal i is evaluated
                            ! = fraction of genes of animal i that derive from j
  allocate(D(N))  ! within family variance of animal i. 
      ! = similar to the F_A (ancestor inbreeding coefficient) term in the path method

  ! check that animals are ordered so that parents precede offspring (ancestor j < focal i for all j & i)
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
  
  !F_meuwissen = F(1:N)  ! exclude F(0) from output
  F(0) = 0d0
  
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
  subroutine sort_pedigree(Ped, IDp, N)
    integer, intent(IN) :: N
    integer, intent(INOUT), allocatable :: Ped(:,:)   ! 2,N
    character(len=nchar_ID), intent(INOUT), allocatable :: IDp(:)
    integer :: Rank(N), i, k, j
    double precision, allocatable :: Gen_dbl(:)
    integer, allocatable :: Gen(:), Ped_new(:,:)   ! 2,N
    character(len=nchar_ID), allocatable :: IDp_new(:)
    
    allocate(Gen(N))
    Gen = get_generations(Ped, N)

    Rank = (/ (i, i=1, N, 1) /)  ! vector to be sorted
    ! sorting algorithm uses double precision input
    allocate(Gen_dbl(N))
    ! add original order as fractional element, to keep original order within each generation
    Gen_dbl = REAL(Gen, 8) + Rank/REAL(N,8)    
    call QsortC(Gen_dbl, Rank)
    deallocate(Gen_dbl)
    
    ! sort individual names in new order
    allocate(IDp_new(N))
    do i=1,N
      IDp_new(i) = IDp(Rank(i))
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
    call move_alloc(IDp_new, IDp)
  
  end subroutine sort_pedigree  
  
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! read pedigree from text file
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
  subroutine read_pedigree(Ped, IDp, FileName, hasHeader) 
    use Global, only: timestamp, printt, quiet
    use sqa_fileIO, only: FileNumRow
    integer, allocatable, intent(OUT) :: Ped(:,:)
    character(len=nchar_ID), allocatable, intent(OUT) :: IDp(:)
    character(len=*), intent(IN) :: FileName
    logical, intent(IN), optional :: hasHeader   ! default=TRUE
    logical :: file_exists, do_read_header, found
    character(len=nchar_ID), allocatable :: NamePed(:,:), IDp_tmp(:)
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
    allocate(IDp(-(2*N):N))
    IDp = 'NA'
    IDp(1:N) = NamePed(1,:)
    
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
            if (NamePed(k+1,i) == IDp(j)) then
              Ped(k,i) = j
              found = .TRUE.
              exit
            endif         
          enddo  ! j
         if (.not. found) then   ! parent does not have own entry pedigree
!           if (.not. quiet)  print *, 'extra : ', N_extra, k, NamePed(:,i) 
           N_extra = N_extra +1
           IDp(-N_extra) = NamePed(k+1,i)
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
      allocate(IDp_tmp(N+N_extra))
      IDp_tmp(1:N_extra) = IDp(-N_extra : -1)
      IDp_tmp((N_extra+1):(N+N_extra)) = IDp(1:N)
      call move_alloc(IDp_tmp, IDp)  ! from, to; deallocates IDp_tmp
            
      ! update total number of individuals     
      N = N + N_extra
    endif
    
    deallocate(NamePed)
    
  
  end subroutine read_pedigree
  

end module pedigree_tools


! ==============================================================================
! ==============================================================================

module Fsnp
  implicit none

!=================================================
contains
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! calc F_uni_u (=Fhat3, mean of ratios) and F_uni_w (ratio of means)
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine calc_Funi(maf_min)
    use Global
    double precision, intent(IN) :: maf_min
    integer :: i,l, y
    ! temp arrays for denominator & numerator for F_uni
    double precision, allocatable :: uni_denom(:), uni_num(:,:), F_uni_tmp(:)
    logical, allocatable :: SNP_OK(:)

    allocate(uni_denom(nSnp))
    uni_denom = 0d0
    allocate(uni_num(-1:2, nSnp))
    uni_num = 0d0 
    
    do l=1,nSnp
      if (p(l) < maf_min)  cycle
      uni_denom(l) = 2*p(l) * (1-p(l))
      do y=0,2  ! number of copies of reference allele
        uni_num(y,l) = y**2 - (1 + 2*p(l)) * y + 2*p(l)**2
      enddo
    enddo

    allocate(F_uni_tmp(nSnp))
    F_uni_tmp = 0d0
    allocate(SNP_OK(nSnp))
    do i=1,Ng
      do l=1,nSnp
        F_uni_tmp(l) = uni_num(Geno(l,i), l)
      enddo
      SNP_OK = p > maf_min .and. Geno(:,i)>=0
      F_uni_u(i) = sum(F_uni_tmp / uni_denom, MASK=SNP_OK) / COUNT(SNP_OK)
      F_uni_w(i) = sum(F_uni_tmp, MASK=SNP_OK) / sum(uni_denom, MASK=SNP_OK)
      ! TODO?: add a shrink towards zero option for non-genotyped SNPs (i.o. excluding)
    enddo
    
    deallocate(uni_denom)
    deallocate(uni_num)
    deallocate(F_uni_tmp)

  end subroutine calc_Funi
  
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! calc genomic relatedness of parent pairs (identical to GCTA & PLINK, so _uni_u)
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine calc_R_parents(maf_min)
    use Global
    double precision, intent(IN) :: maf_min
    integer :: i,l,j, xi, xj, o, z
    logical, allocatable :: SNP_OK(:)
    double precision, allocatable :: RelA(:,:,:), R_tmp(:)  ! lookup table; per-SNP contributions
      
    ! equation 3 in Yang et al, Am J Hum Genet 2011
    allocate(RelA(-1:2,-1:2,nSnp))   
    RelA = 0d0   ! if either genotype is missing, contribution of that SNP is 0 (shrink towards zero)
    do l=1,nSnp
      if (p(l) < maf_min)  cycle
      do xi=0,2  ! genotype of 1st indiv
        do xj=0,2  ! genotype of 2nd indiv
          RelA(xj,xi,l) = (xi - 2*p(l)) * (xj - 2*p(l)) / (2*p(l)*(1-p(l)))
        enddo
      enddo
    enddo
 
    ! loop over the pedigree, find SNP data indices from parents & calc R
    R_parents = 0d0
    allocate(R_tmp(nSnp))
    R_tmp = 0d0
    allocate(SNP_OK(nSnp))
    do o=1, Np
      if (any(pedigree(:,o) == 0))  cycle   ! either parent unknown --> by definition unrelated
      i = 0
      j = 0
      do z=1,Ng
        if (IDg(z) == IDp(pedigree(1,o)))  i = z
        if (IDg(z) == IDp(pedigree(2,o)))  j = z
      enddo
      if (i==0 .or. j==0)  cycle   ! either parent not genotyped
      do l=1,nSnp
        R_tmp(l) = RelA(Geno(l,i), Geno(l,j), l)  ! get value from lookup table
      enddo    
!      R_parents(o) = SUM(R_tmp) / nSnp    ! shrink towards zero if genotype unknown
      SNP_OK = p > maf_min .and. Geno(:,i)>=0 .and. Geno(:,j)>=0  ! TODO: transpose geno
      R_parents(o) = SUM(R_tmp, MASK=SNP_OK) / COUNT(SNP_OK)
      
      ! TODO: add weighed version. 
    enddo
  
  end subroutine calc_R_parents


end module Fsnp



! ==============================================================================
! main program
! ==============================================================================
program main
  use sqa_fileIO, ONLY: read_geno, write_geno, nchar_ID
  use Global
  use pedigree_tools, only: read_pedigree, sort_pedigree
  use Fped_module
  use Fsnp, ONLY: calc_Funi, calc_R_parents
  implicit none
  
  character(len=200) :: geno_file, af_file, pedigree_file, sorted_pedigree_file, out_file
  character(len=3) :: geno_format
  double precision :: maf_min 
  double precision, allocatable :: Fped(:)
  character(len=nchar_ID), allocatable :: IDz(:)
  integer :: i, ig, ip, y
  logical :: pedigree_is_sorted, return_sorted_pedigree
  
  call read_args()   ! set defaults & read command line arguments

  ! read pedigree 
  if (pedigree_file /= 'nofile') then
    if (.not. quiet)  call printt('reading pedigree file...')
    call read_pedigree(pedigree, IDp, pedigree_file)
    Np = size(pedigree, dim=2)
    call timestamp()
    if (.not. quiet)  print *, 'read ', Np, ' individuals'
  
   ! sort pedigree (parents always before their offspring) if needed
    pedigree_is_sorted = .TRUE.
    do i=1,Np
      if (pedigree(1,i) > i .or. pedigree(2,i) > i) then
        pedigree_is_sorted = .FALSE.
        exit
      endif
    enddo
    if (.not. pedigree_is_sorted) then
      if (.not. quiet)  call printt('sorting pedigree ...')
      call sort_pedigree(pedigree, IDp, Np)
    endif
  else
    Np = 0
  endif
  
  ! read genomic data 
  if (geno_file /= 'nofile') then
    if (.not. quiet)   call printt('reading genotypes ...')
    call read_geno(Geno=Geno, ID=IDg, SNP_names=SNP_names, FileName=geno_file, &
      FileFormat = geno_format, transp=.TRUE.)  
    Ng = SIZE(Geno, DIM=2) -1 ! dim 0:Ng
    nSnp  = SIZE(Geno, DIM=1)
    call timestamp()
    if (.not. quiet)  print *, 'read ', Ng, ' individuals by ', nSnp, ' SNPs'
  else
    Ng = 0
  endif  
  
  
  
  ! combine pedigree & genomic IDs into single vector 
  ! (first all IDp in order, then any genotyped individuals not in pedigree) 
  call combine_ped_geno()
    
  if (Np>0 .and. return_sorted_pedigree) then
    allocate(IDz(0:Np))
    IDz(0) = 'NA'
    IDz(1:Np) = IDp
    sorted_pedigree_file = pedigree_file(1:(LEN_TRIM(pedigree_file)-4))//'_sorted.txt'    
    open(unit=101, file=trim(sorted_pedigree_file), status='unknown')
      write(101, '(3a20)') 'id', 'dam', 'sire'
      do i=1,Np
        write(101, '(3a20, i6)')  IDz(i), IDz(pedigree(:,i))
      enddo
    close(101) 
  endif


  ! calculate pedigree inbreeding coefficients
  allocate(Fped(0:Np))
  if (Np > 0) then
    if (.not. quiet)  call printt('calculating inbreeding coefficients...')
    Fped = F_meuwissen(pedigree)  
    call timestamp()
    write(*, '("Number of inbred animals in pedigree: ", i8)')  count(Fped > 0.000001d0)
  else
    Fped = 0d0
  endif

  
  ! calculate genomic inbreeding coefficients  
  allocate(F_uni_u(0:Ng))
  allocate(F_uni_w(0:Ng))
  F_uni_u = 0d0
  F_uni_w = 0d0
  if (Ng > 0) then
    p = getAF(af_file) 
    call calc_Funi(maf_min)   ! TODO: return 0:N not 1:N
  endif
  
  write(*, '("MAF >= ", f6.4, "  # SNPs: ", i10)') maf_min, COUNT(p >= maf_min)

  
  ! calculate genomic relatedness between parents
  allocate(R_parents(0:Nc))
  R_parents = 0d0
  call calc_R_parents(maf_min)
  

  ! write output to file
  open(unit=101, file=trim(out_file), status='unknown')
    write(101, '(a2, 40x, 5a12)') 'ID', 'F_ped', 'F_uni_u', 'F_uni_w', 'R_uni_u' !, 'R_uni_w'
    do i=1,Nc
      ip = 0
      ig = 0
      if (any(IDp == IDc(i))) ip = i
      do y=1,Ng
        if (IDg(y) == IDc(i)) ig = y 
      enddo   
      write(101, '(a40, 2x, 5f12.8)')  IDc(i), Fped(ip), F_uni_u(ig), F_uni_w(ig), R_parents(i)
    enddo
  close(101)


  
!=================================================
contains 

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! read command line arguments
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
  subroutine read_args()
    use sqa_fileIO, only: valid_formats    ! genotype formats: SEQ, PED, RAW, LMT
    integer :: nArg, i
    character(len=32) :: arg, argOption
    
    ! set defaults
    pedigree_file = 'nofile'  !'pedigree.txt'     
    geno_file = 'nofile'  ! 'geno.txt'
    geno_format = 'SEQ'
    af_file = 'NoFile'
    out_file = 'inbreeding_coefficients.txt'
    return_sorted_pedigree = .FALSE.
    quiet = .FALSE.
    maf_min = 0.01d0
    
    
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
        case ('--geno')
          i = i+1
          call get_command_argument(i, geno_file)
        case ('--informat', '--inFormat', '--genoformat', '--genoFormat')
          i = i+1
          call get_command_argument(i, geno_format)
          if (.not. any(valid_formats == geno_format)) then
            print *, 'ERROR: informat must be one of: ', valid_formats
            stop
          endif
        case ('--maf')
          i = i+1
          call get_command_argument(i, argOption)
          read(argOption, '(a2)') maf_min
        case ('--freq', '--af')
          i = i+1
          call get_command_argument(i, af_file)     
        case ('--out')
          i = i+1
          call get_command_argument(i, out_file)    
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
    
    if (pedigree_file == 'nofile' .and. geno_file == 'nofile') then
      print *, 'please specify pedigree and/or genotype file'
      stop
    endif
    

  end subroutine read_args  
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! combine pedigree & genomic ID vectors into a single vector
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine combine_ped_geno()
    integer :: i
    character(len=nchar_ID), allocatable :: IDc_tmp(:) 
  
    if (Ng == 0) then
      Nc = Np
      IDc = IDp 
      return
    else if (Np == 0) then
      Nc = Ng
      IDc = IDg
      return
    endif
    
    Nc = Np  ! total number, combined pedigree + genotype data
    allocate(IDc_tmp(Np + Ng))   ! maximum Nc: none in pedigree are genotyped
    IDc_tmp(1:Np) = IDp(1:Np)
    ! add any genotyped-but-not-in-pedigree individuals
    do i=1,Ng
      if (any(IDp == IDg(i)))  cycle
      Nc = Nc +1
      IDc_tmp(Nc) = IDg(i)
    enddo
    
    allocate(IDc(Nc))
    IDc(1:Nc) = IDc_tmp(1:Nc)
    deallocate(IDc_tmp)
  
  end subroutine combine_ped_geno
  
  
end program main
! ==============================================================================