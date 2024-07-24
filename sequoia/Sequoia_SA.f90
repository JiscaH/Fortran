! Author: Jisca Huisman,  jisca.huisman@gmail.com
! Most of this code was written as a post doc in Evolutionary biology 
! at the University of Edinburgh, UK.
!
! This code is available under GNU General Public License v3
!
! The program is described in the paper
! "Pedigree reconstruction from SNP data: 
! Parentage assignment, sibship clustering, and beyond", 
! in Molecular Ecology Resources, 2017
!
! Updates are made available at https://github.com/JiscaH , 
! as well as the R version. 
!
! to compile:
! gfortran -fall-intrinsics -O3 Sequoia_SA.f90 -o sequoia
! to debug:
! gfortran -fall-intrinsics -Wall -pedantic -fbounds-check -g -Og Sequoia_SA.f90 -o sequoia
!
! ####################################################################
! @@@@   MODULES   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! ####################################################################
module Global 
use sqa_fileIO, ONLY: ishort   ! ishort=selected_int_kind(1)
use sqa_general, ONLY: AHWE, OHWE, OcA, AKAP, OKAP, OKOP, AKA2P, OKA2P, &
  timestamp, printt
implicit none

character(len=*), parameter :: version = 's2.11.1 (24 July 2024)'
integer :: nInd, nSnp, nIndLH, maxSibSize, MaxOppHom, MaxMendelE, MaxMismatchDup, &
  Hermaphrodites, nC(2), nYears, maxAgePO, nPairs, Complx, quiet, AgePhase, BYzero, &
  ID_len, mxCP, viginti(20)
integer, parameter :: mxA=2**6, & ! max no. ancestors considered when testing for pedigree loop
!   mxCP = 50, &  ! max no. candidate parents per sex   -- now readspecs()
   MaxMaxAgePO = 101, &  ! maximum of MaxAgePO
   nchar_filename = 2000, &
   nchar_ID = 40, &     ! max. ID nchar 
   XP = 5  ! multiplier nInd --> max no. candidate sib pairs 
logical :: DoMtDif, DoSibs, AnyYearLast
logical, allocatable, dimension(:) :: ToCheck(:), SelfedIndiv(:), skip(:), IsNewSibship(:,:), mtDif(:,:)
integer, allocatable, dimension(:) :: Sex, BY, nFS, Mate, YearLast
integer,allocatable,dimension(:,:) :: Parent, nS, FSID, DumMate, DumClone
integer(kind=ishort), allocatable :: Genos(:,:)
integer, allocatable, dimension(:,:,:) :: SibID, GpID
double precision :: TF, TA, zero
double precision, parameter ::  missing = 999D0, impossible=777D0, &
  NotImplemented = 444D0, MaybeOtherParent = 222D0, AlreadyAss = 888D0
double precision, allocatable, dimension(:) ::  Lind
double precision, allocatable, dimension(:,:) :: CLL
double precision, allocatable, dimension(:,:,:) :: PHS, PFS, IndBY, AgePriorA, LindX
double precision, allocatable, dimension(:,:,:,:) :: DumP, DumBY!, FSLik
double precision, allocatable, dimension(:,:,:,:,:) :: XPr
 character(len=20) :: DumPrefix(3)
 character(len=nchar_ID), allocatable :: Id(:), DummyNamesIO(:,:)


 !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
  contains
pure function MaxLL(V)
double precision, intent(IN) :: V(:)
double precision :: MaxLL
MaxLL = missing
if (ANY(V < 0 .and. V>-HUGE(0.0D0))) then
  MaxLL = MAXVAL(V, mask = (V<0 .and. V>-HUGE(0.0D0)), DIM=1)
else
  MaxLL = MINVAL(V, mask = (V>-HUGE(0.0D0)), DIM=1)  
  ! impossible: can't do; AlreadyAss: already is; missing: not calc'd
  ! V should not ever be -INF
endif
end function MaxLL

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pure function addALR(LLg, ALR)   
double precision, intent(IN) :: LLg, ALR
double precision :: addALR

if (LLg > 0) then
  addALR = LLg
else if (ALR == impossible) then
  addALR = impossible
else
  addALR = LLg + ALR
endif

end function addALR

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function getPar(A, kA)
integer, intent(IN) :: A, kA
integer :: getPar(2)

if (A > 0) then
  getPar = Parent(A,:)
else if (A < 0) then
  getPar = GpID(:, -A, kA)
else
  getPar = 0
endif

end function getPar
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function getAP(AgeD, Rel, k, m, noGo)  
! Rel: 1=PO, 2=FS, 3=HS, 4=GP, 5=FA, 6=HA
! k: sex of 2nd indiv (except for HA)
! m: related via mat/pat

integer, intent(IN) :: AgeD, Rel, k, m
double precision :: getAP
double precision :: AM(5,5), noGo
integer :: D2, D3

getAP = zero
if (AgeD == 999) return
if (Rel < 1 .or. Rel > 6)  call Erstop('getAP: illegal Rel', .TRUE.)
if (Rel == 1 .and. AgeD <=0)  getAP = LOG10(zero)
if (Rel == 4 .and. AgeD <=1)  getAP = LOG10(zero)
if (AgeD < -MaxAgePO)  getAP = LOG10(zero)
if (Rel == 1 .and. AgeD > MaxAgePO)  getAP = LOG10(zero)
if (getAP < -HUGE(0.0D0)) then
  getAP = noGo
  return
endif

if (((m<1 .or.m>4) .and. Rel>2) .or. &
  ((k<1 .or.k>4) .and. (Rel==1 .or. Rel==4 .or. Rel==6)))  then
  call Erstop("getAP: illegal k or m!", .TRUE.)
endif

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! AgePriorA D2 + D3:
!  1    2     3
!  M   MGM   PGM
!  P   MGP   PGP
! FS   MFA   PFA
! MS  MMHA  PMHA
! PS  MPHA  PPHA
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if (Rel < 4) then
  D3 = 1
else
  D3 = 1+m  ! via mother/father
endif
if (Rel==1 .or. Rel==4) then 
  D2 = k  ! sex of parent/GP
else if (Rel==2 .or. Rel==5) then  ! FS/FA
  D2 = 3
else if (Rel==3) then
  D2 = 3+m       
else if (Rel==6) then   ! HA (k HS of parent m)
  D2 = 3+k
else
  D2 = 0
  call ErStop("getAP: illegal Rel", .TRUE.)
endif

AM(:,1:3) = AgePriorA(AgeD,:,:)
AM(:,4) = (AM(:,2) + AM(:,3)) / 2.0    ! m=3: unknown via which parent
AM(:,5) = AM(:,4)    ! m=4: hermaphrodite = unknown sex

if ((Rel==1 .or. Rel==4) .and. k>2) then  ! unknown (grand)parent sex
  getAP = SUM(AM(1:2, D3)) / 2.0
else if ((Rel==3 .and. m>2) .or. (Rel==6 .and. k>2)) then
  getAP = SUM(AM(4:5, D3)) / 2.0
else 
  getAP = AM(D2, D3)
endif

if (getAP/=getAP .or. getAP==zero) then
  getAP = noGo
else
  getAP = LOG10(getAP)
endif

end function getAP

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pure function AgeDiff(i,j)   

integer, intent(IN) :: i,j
integer :: AgeDiff

if (BY(i)>=0 .and. BY(j)>=0) then
  AgeDiff = BY(i) - BY(j)
else
  AgeDiff = 999
endif

end function AgeDiff

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
subroutine Rprint(message, IntData, DblData, DataType)
implicit none

character(len=*), intent(IN) :: message
integer, intent(IN) :: IntData(:)
double precision, intent(IN) :: DblData(:)
character(3), intent(IN) :: DataType
!character(len=200) :: dblfmt
!integer :: nchar, ndata
!integer :: IntDummy(0)

! nchar = LEN(trim(message))

! if (DataType == "DBL") then
  ! ndata = SIZE(DblData)
  ! call dblepr(trim(message), nchar, DblData, ndata)
! else if (DataType == "INT") then
  ! ndata = SIZE(IntData)
  ! call intpr(trim(message), nchar, IntData, ndata)
! else if (DataType == "NON") then
  ! call intpr(trim(message), nchar, IntDummy, 0) 
! else
  ! call ErStop("invalid DataType for Rprint")
! endif

call timestamp(.FALSE.)

if (DataType == "DBL") then
!  dblfmt = "'("//trim(message)//", 50f10.3)'"
!  write(*, '(a200, 50f10.3)') trim(message), DblData
  print *, message, DblData
else if (DataType == "INT") then
  print *, message, IntData
else if (DataType == "NON") then
  print *, adjustl(trim(message))
endif

end subroutine Rprint

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
subroutine rchkusr    ! stand-in for R subroutine
! do nothing
end subroutine rchkusr

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
! print progress in (slow) loop over all individuals
subroutine print_progress(xi,xc)
  integer, intent(IN) :: xi  ! individual number
  integer, intent(INOUT) :: xc  ! chunk number

  call timestamp(.FALSE.)
  print *, xi, '  ', xc*5, '%'
  xc = xc +1
end subroutine print_progress


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! create a sequence from 1 to n, with length equal to n_steps
function mk_seq(n, n_steps) result(seq) 
  integer, intent(IN) :: n, n_steps
  integer, allocatable :: seq(:)
  integer :: i
  real, allocatable :: probs(:)
  
  allocate(probs(n_steps))
  allocate(seq(n_steps))
  
  probs = (/ (i,i=1,n_steps) /) / real(n_steps)
  seq = NINT( probs * n )
  if (seq(1) == 0)  seq(1) = 1
  if (seq(n_steps) > n)  seq(n_steps) = n
  
end function mk_seq
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

end module Global

!===============================================================================

module OHfun
  use Global, ONLY : nSnp, Genos, maxOppHom, OHWE, OKOP
  implicit none
  
  integer :: IsOppHom(-1:2,-1:2), Ecnt(3,3,3)
  double precision, allocatable :: PPO(:,:,:)

  contains
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine init_OH()
  
    integer :: l,i,j
  
    IsOppHom = 0
    IsOppHom(0, 2) = 1
    IsOppHom(2, 0) = 1
    
    ! offspring - dam - sire
    Ecnt(:,1,1) = (/ 0, 1, 2 /)
    Ecnt(:,1,2) = (/ 0, 0, 1 /)
    Ecnt(:,1,3) = (/ 1, 0, 1 /)

    Ecnt(:,2,1) = (/ 0, 0, 1 /)
    Ecnt(:,2,2) = (/ 0, 0, 0 /)
    Ecnt(:,2,3) = (/ 1, 0, 0 /)

    Ecnt(:,3,1) = (/ 1, 0, 1 /)
    Ecnt(:,3,2) = (/ 1, 0, 0 /)
    Ecnt(:,3,3) = (/ 2, 1, 0 /)
    
    allocate(PPO(-1:2,-1:2,nSnp))
    PPO = 1D0
    do l=1,nSnp
      do i=0,2  ! obs offspring 
        do j=0,2    ! obs parent
          PPO(i,j,l) = OKOP(i,j,l) / OHWE(i,l)
        enddo
      enddo
    enddo 
  
  end subroutine init_OH
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  integer function calcOH(i,j)  ! no. opp. hom. loci

    integer, intent(IN) :: i,j
    integer :: l, OH

    OH = 0
    do l=1,nSnp
      OH = OH + IsOppHom(Genos(l,i), Genos(l,j))
      if (OH > maxOppHom) exit                      
    enddo
    calcOH = OH

  end function calcOH

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  pure function QLR_PO(i,j)  !quick check, not conditioning on parents.

  integer, intent(IN) :: i,j
  double precision :: QLR_PO
  integer :: l
  double precision :: PrL(nSnp)

  PrL = 0D0
  do l=1,nSnp
    PrL(l) = LOG10(PPO(Genos(l,i), Genos(l,j), l))  
  enddo
  QLR_PO = SUM(PrL)

  end function QLR_PO
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  integer function CalcTrioErr(A,par)

    integer, intent(IN) :: A, Par(2)
    integer :: l,k, ME

    ME = 0
    do l=1,nSnp
      if (Genos(l,A)==-1 .or. ALL(Genos(l, Par)==-1)) then
        cycle
      else if (ANY(Genos(l, Par)==-1)) then
        do k=1,2
          if (Genos(l, Par(k))==-1) then
            if (((Genos(l,A)==0).and.(Genos(l,Par(3-k))==2)) .or. &
             ((Genos(l,A)==2).and.(Genos(l,Par(3-k))==0))) then
              ME = ME +1
              cycle
            endif
          endif
        enddo
      else
        ME = ME + Ecnt(Genos(l,A)+1, Genos(l, Par(1))+1, Genos(l, Par(2))+1)
      endif
    enddo
    CalcTrioErr = ME

  end function CalcTrioErr
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

end module OHfun

!===============================================================================

module CalcLik
  use Global, ONLY : Genos, Parent, nFS, FSID, maxSibSize, hermaphrodites, DumClone,&
    AKA2P, OKA2P
  implicit none

  contains
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine CalcLind(i)
    use Global
    implicit none

    integer, intent(IN) :: i
    integer :: nFSi, FSi(MaxSibSize), Inbr(2), Anc(2,mxA), PIK, l, x, y, k, z, j
    double precision :: PrL(nSnp), Px(3,2), PrG(3), PrXYZ(3,3,3), PrX(3)
    logical :: Selfed

    if (real(COUNT(Genos(:,i)/=-1)) < nSnp/20.0) then
      return 
    endif

    nFSi = 1 
    FSi = 0                
    if (Parent(i,1)<0 .or. Parent(i,2)<0) then
      nFSi = nFS(FSID(maxSibSize+1,i))
      FSi = FSID(1:maxSibSize, FSID(maxSibSize+1,i))  
    endif

    ! PO- and GP-mating
    Inbr = 0
    Anc = 0
    PIK = 0
    if (Parent(i,1)/=0 .and. Parent(i,2)/=0) then
      call getAncest(i,1,Anc)
      do k=1,2
        if (Anc(3-k,5-k) == Parent(i,3-k)) then
          Inbr(k) = 1
          if (Parent(i,3-k)>0) then
            PIK = Parent(i,3-k)
          endif
        endif
      enddo
    endif

    Selfed = .FALSE.
    if (hermaphrodites/=0) then
      if (Parent(i,1)==Parent(i,2) .and. Parent(i,1)>0) then
        Selfed = .TRUE.
      else if (all(Parent(i,:) < 0)) then
        if (DumClone(-Parent(i,1),1) == -parent(i,2)) then  
          Selfed = .TRUE.
        endif
      endif
    endif

    PrL = 0D0
    do l=1,nSnp
      do k=1,2
        call ParProb(l, Parent(i,k), k, i,-1, Px(:,k))
        if (Inbr(k)==1) then
          call ParProb(l, Anc(k,k+2), k, PIK,0, PrG)
        endif
      enddo
      PrXYZ = 0D0
      do y=1,3
        do z=1,3
          if (Selfed .and. y/=z)  cycle
          if (ANY(Inbr==1)) then
            do k=1,2
              if (Inbr(k)==-1) then
                PrXYZ(:,y,z) =AKA2P(:,y,z) * SUM(AKA2P(y,z,:)*PrG) *&
                 Px(y,k) * Px(z,3-k)
              endif
            enddo
          else if (Selfed) then   
            PrXYZ(:,y,y) = AKA2P(:, y, y) * Px(y,1)
          else
            PrXYZ(:,y,z) = AKA2P(:, y, z) * Px(y,1) * Px(z,2)
          endif 
          if (nFSi > 1) then
            do j=1, nFSi
              if (FSi(j) == i) cycle
              PrXYZ(:,y,z) = PrXYZ(:,y,z) * OKA2P(Genos(l,FSi(j)),y, z)
            enddo
          endif
        enddo
      enddo

      PrX = 0D0
      if(SUM(PrXYZ) > 0D0) then  
        do x=1,3
          PrX(x) = SUM(PrXYZ(x,:,:))/SUM(PrXYZ)
        enddo
        PrX = PrX / SUM(PrX)
      endif
      PrX = OcA(:,Genos(l,i)) * PrX
      PrL(l) = LOG10(SUM(PrX))
      LindX(:,l, i) = PrX
    enddo
    Lind(i) = SUM(PrL)

    if (Lind(i)> 0D0 .or. Lind(i)/=Lind(i)  .or. Lind(i) < -HUGE(1D0) .or. &   
      any(LindX(:,:,i)/=LindX(:,:,i))) then  
      call Rprint("",(/i,Parent(i,:)/), (/0D0/), "INT")
      if (any(LindX(:,:,i)/=LindX(:,:,i)))  call Rprint("LindX NA", (/0/), (/0D0/), "NON")
      call Erstop("Invalid individual LL - try increasing Err", .FALSE.)
    endif

  end subroutine CalcLind

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
  function FSLik(l,i)    ! LL of FS set

    integer, intent(IN) :: l,i
    double precision :: FSLik(3,3), FSLikTMP(3,3)
    integer :: j, x, y

    FSLik = 1D0
    if (nFS(i)==0) then
      FSLik = 1D0
    else
      FSLikTMP = 1D0
      do j=1, nFS(i)
        do y=1,3
          do x=1,3
            FSLikTMP(x,y) = FSLikTMP(x,y) * OKA2P(Genos(l,FSID(j,i)), x, y)
          enddo
        enddo
      enddo
      if (any(FSLikTMP/=FSLikTMP) .or. any(FSLikTMP>1.0D0))  then 
        call Erstop("Invalid FS LL", .TRUE.)
      endif
      FSLik = FSLikTMP
    endif
    
  end function FSLik

end module CalcLik


! ####################################################################
! Recursive Fortran 95 quicksort routine
! sorts real numbers into ascending numerical order
! Author: Juli Rew, SCD Consulting (juliana@ucar.edu), 9/03
! Based on algorithm from Cormen et al., Introduction to Algorithms,1997
! Made F conformant by Walt Brainerd

! Adapted by J Huisman (jisca.huisman@gmail.com) to output rank, to
! enable sorting of parallel vectors, and changed to decreasing rather 
! than increasing order

module qsort_c_module
implicit none
public :: QsortC
private :: Partition

 contains
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

end module qsort_c_module

! #####################################################################
 
subroutine Erstop(message, bug)
use Global
implicit none

character(len=*), intent(IN) :: message
logical, intent(IN) :: bug
! Error = 1
 call DeAllocAll
 !call rexit("  ERROR! ***"//message//"***")
print *, ""
print *, " *** ERROR! ***"
print *, message
print *, ""
if (bug)  print *, "Please report bug"
print *, ""
stop 'error'   ! String of no more that 5 digits or a character constant 

end subroutine Erstop

! ####################################################################

! @@@@   PROGRAM   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

! ####################################################################

program Main
use Global
use sqa_fileIO, ONLY: checkFile, FileNumRow
implicit none

integer :: x, i, CalcLLR, AgeEffect, FindMaybe(2), ResumePed, &
  nAmbMax(2), FindMaybeX, NP, date_time_values(8)
double precision :: TotLL(42), Er, ErV(3)
character(len=2) :: ResumePedC, HermC
character(len=3) :: ErrFlavour, GenoFormat
character(len=nchar_filename) :: PedigreeFileName, PairsFileName, OutFileName, &
  GenoFileName, LifehistFileName, AgePriorFileName, OnlyListFileName, &
  AF_FileName, mt_FileName, SpecsFileName
character(len=nchar_filename) :: FN(7)
character(len=nchar_ID), allocatable :: ID_tmp(:)
logical :: DoDup, DoPar, DoPairs, FileOK, DoReadParents, dupQuiet, &
  MaybePO_onlyOH, NotDup, withAssignmentLog, file_exists


call set_defaults()
call get_SpecsFileName()
inquire(file=trim(SpecsFileName), exist = file_exists)
if (file_exists)  call ReadSpecs()  ! read arguments from SpecsFile
call read_args()  ! read command line arguments

if (GenoFileName == 'NoFile') then
  write(*,*)  "Please specify genotype file, either in SequoiaSpecs.txt or with --geno"
  stop
endif

inquire(file=trim(LifeHistFileName), exist = file_exists)
if (.not. file_exists .and. LifehistFileName == 'LifeHist.txt') LifehistFileName = 'NoFile'  ! default file name
if (LifeHistFileName == 'NoFile') then
  if (AgePriorFileName/="AgePriors.txt") then
    write(*,*) "WARNING: the specified ageprior file will not be used if there is no LifeHistFile"
  endif
  AgePriorFileName = 'NoFile'
  if (quiet < 1 .and. (DoPar .or. DoSibs)) then
    print *, "NOTE: expect lower assignment rate when no LifeHistFile is provided"
  endif
endif

! check if files exist
FN = (/PedigreeFileName, PairsFileName, LifehistFileName, &
  AgePriorFileName, OnlyListFileName, AF_FileName, mt_FileName/)
! GenoFileName: checked & extension added by read_geno in sqa_fileIO.f90
do i=1,7
  if (FN(i) /= 'NoFile')   call checkFile(FN(i))
enddo

!=========================
call date_and_time(VALUES=date_time_values)  ! start time                                                        

if (quiet < 1) then
  write(*,*) ""
  write(*,*) "@@@@@@@@@@@@@@@@@@@@@@@"
  write(*,*) "@@@@    SEQUOIA    @@@@"
  write(*,*) "@@@@@@@@@@@@@@@@@@@@@@@"
  write(*,*) "||  Jisca Huisman  |  https://github.com/JiscaH  |  GNU GPLv3  ||"
  write(*,*) "" 
endif

! initiate all arrays
call Initiate() !GenoFileName, GenoFormat, LifehistFileName, AgePriorFileName, PedigreeFileName, &
! OnlyListFileName, Er, ErrFlavour, AF_FileName, mt_FileName)  

if (DoDup .or. DoPar .or. DoSibs .or. DoPairs .or. any(FindMaybe==1)) then  !  .or. CalcLLR==1
  if (quiet < 1)  call printt('updating all probs ... ')  ! timestamp() + print *,
  call UpdateAllProbs()
else
  MaxOppHom = nSnp   ! only calculate OH for all parents in --pedigreeIN
endif

! overview of settings
select case(Hermaphrodites)
  case (1)
    HermC = "A"
  case (2)
    HermC = "B"
  case default
    HermC = "X"
end select

if (quiet < 1) then
  write(*,*) ""
  write(*,*) "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
  write(*, '(" @ start time    : ", i4,"-",i2.2,"-",i2.2, 2x, i2.2,":",i2.2)') date_time_values(1:3), date_time_values(5:6)
  write(*, '(" @ Geno          : ", i7, " individuals x ", i6, " SNPs")')  nInd, nSnp
  write(*, '(" @ Sex           : Female: ", i6, "  Male:", i6)') count(sex==1), count(sex==2) 
  write(*, '(" @               : Unknown:", i6, "  Hermaphrodite:", i6)') count(sex==3), count(sex==4)
  write(*, '(" @ Birth year    : min* - max: ", 2i5)') BYzero +1, BYzero +nYears 
  write(*,*) "@  (*: earliest birthyear for a grandparent of oldest individual)"
  write(*, '(" @ Max age parent: ", i7)') maxAgePO
  if (any(skip)) write(*,*) "@ N in --only   : ", COUNT(.not. skip)
  write(*,*) "@ Pedigree-IN   : ", trim(PedigreeFileName)
  write(*,'(" @ Genotyping error rate : ", f7.5)') Er  
  write(*,'(" @ Max mismatch  :   DUP:", i4, "  OH:", i4,  "  ME:", i4)') &
   MaxMismatchDup, MaxOppHom, MaxMendelE
  write(*,'(" @ Duplicates    : ", l3)')  DoDup
  write(*,'(" @ Parentage     : ", l3)')  DoPar
  write(*,'(" @ Full Pedigree : ", l3)')  DoSibs
  write(*,'(" @ Complexity    : ", i3)')  Complx
  write(*,'(" @ Age effect    : ", i3)')  AgeEffect
  write(*,'(" @ Hermaphrodite : ",2X,a1)')  HermC
  write(*,'(" @ Parent LLR    : ", l3)')  CalcLLR ==1 
  write(*,'(" @ Maybe P-O     : ", l3)')  FindMaybe(1) == 1
  write(*,'(" @ Maybe Related : ", l3)')  FindMaybe(2) == 1
  write(*,'(" @ Pair LLs      : ", l3)')  DoPairs
  write(*,*) "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
  write(*,*) ""
endif


!====================================
! Duplicate check
!====================================
if (DoDup) then
  if (.not. dupQuiet) then
    print *, ""
    call printt("~~~~~~~  Duplicate Check  ~~~~~~~~")
    print *, ""
  endif
  call duplicates(dupQuiet)
endif

!====================================
! parentage assignment
!====================================
if (DoPar) then
  print *, ""
  call printt("~~~~~~~  Parentage Assignment  ~~~~~~~~")
  if (PedigreeFileName /= "NoFile") then
    ! ReadPedFile() called by initiate()
    if (quiet==-1) call printt("Checking pedigreeIN ...")
    call CheckPedigree(.TRUE., .TRUE.)   ! genotyped parents only ; drop parents of any in --only
    if (quiet==-1) call timestamp()
    if(quiet==-1)  write(*, '("Total LogLik: ", f12.1, "  # parents:", 2i6)') &
          SUM(Lind), count(Parent/=0, DIM=1)
  endif
  print *, ""
  call parents(TotLL, withAssignmentLog)
  if (DoSibs .or. OutFileName == "NoFile") then
    call writeped("Parents.txt", CalcLLR==1)
  else
    call writeped(OutFileName, CalcLLR==1) 
  endif
  call writeAgePrior
  call WriteBYProb
endif
  
!====================================
! full pedigree reconstruction
!====================================
if (DoSibs) then
  print *, ""
  call printt("~~~~~~~  Full Pedigree Reconstruction  ~~~~~~~~")
  print *, ""

  if (ResumePedC /= "XX") then
    call timestamp()
    print *, "Resuming at round: ", ResumePed
  else if (DoReadParents) then
    if (trim(PedigreeFileName) == "NoFile")  PedigreeFileName = "Parents.txt"
    call checkFile(PedigreeFileName)
    call ReadPedFile(PedigreeFileName)
    call UpdateAllProbs()
!    if (quiet<1)  print *, "# parents: ", COUNT(Parent /= 0, DIM=1)
    if(quiet==-1)  call Rprint("Checking parents / pedigreeIN ...", (/0/), (/0.0D0/), "NON")
    if (trim(PedigreeFileName) == "Parents.txt") then   ! or if (DoPar) ?
      call CheckPedigree(.FALSE., .FALSE.)  ! double check parents, using updated ageprior 
    else 
      call CheckPedigree(.FALSE., .TRUE.)  ! drop parents of any in --only
    endif
    call UpdateAllProbs()
    if(quiet==-1) call timestamp(.TRUE.)
    if(quiet==-1)  write(*, '("Total LogLik: ", f12.1, "  # parents:", 2i6)') &
            SUM(Lind), count(Parent/=0, DIM=1)
  else if (quiet<1) then
    call Rprint("NOTE: not starting from assigned parents / pedigreeIN ...", (/0/), (/0.0D0/), "NON")  
  endif
  
  call sibships(AgeEffect, ResumePed)
!  if (CalcLLR == 1)   call CalcParentLLR
  if (OutFileName == "NoFile") OutFileName = "Pedigree_seq.txt"
  call writeped(OutFileName, CalcLLR==1)
  call writeAgePrior
  if(quiet<1) call printt("Write dummy parents ...")
  call WriteDummies  ! slow?
  call WriteBYProb    
endif

!====================================
! Calculate parent LLR's for read-in pedigree (analogous to function calcOHLLR() in R)
!====================================
if (trim(PedigreeFileName) /= "NoFile" .and. .not. &
  (DoPar .or. DoSibs .or. DoPairs .or. any(FindMaybe==1))) then
  if (CalcLLR==1) then
    print *, ""
    call printt("~~~~~~~  Parent LLR only  ~~~~~~~~")
    print *, ""
!    call CalcParentLLR   ! done by writeped()
    if (OutFileName == "NoFile") then
      i = index(PedigreeFileName, ".txt")   ! find location of ".txt"
      OutFileName = PedigreeFileName(1:(i-1))//"_LLR.txt"
    endif
    call writeped(OutFileName, .TRUE.)
  else
    print *, ""
    call printt("~~~~~~~  Parent OH only  ~~~~~~~~")
    print *, ""
    if (OutFileName == "NoFile") then
      i = index(PedigreeFileName, ".txt")   
      OutFileName = PedigreeFileName(1:(i-1))//"_OH.txt"
    endif
    call writeped(OutFileName, .FALSE.)
  endif
endif

!====================================
! Identify likely (remaining) relatives
!====================================
if (ANY(FindMaybe==1)) then
  if (OutFileName == "NoFile" .and. PedigreeFileName /= "NoFile") then
    OutFileName = PedigreeFileName
  endif
  
  do x=1,2
    if (FindMaybe(x)==1) then
      if (quiet < 1) then
        print *, ""
        if (x==1) then
          call printt("~~~~~~~  Checking for likely Parent - Offspring pairs  ~~~~~~~~")
        else
          call printt("~~~~~~~  Checking for likely relatives  ~~~~~~~~")
        endif
        if (PedigreeFileName /= "NoFile") then
          call printt("~~~~~~~  (Conditional on pedigree "// trim(OutFileName) //") ~~~~~~~~")
        endif
        print *, ""
      endif
      
      if (nAmbMax(x) == -99)  nAmbMax(x) = 7*nInd
      call findambig(x, nAmbMax(x), MaybePO_onlyOH) 
    endif
  enddo
endif


!====================================
! calc likelihoods for specified pairs
!====================================
if (DoPairs) then
  NP = FileNumRow(PairsFileName) -1  ! number of pairs; file expected to have header
  print *, ""
  call printt("~~~~~~~  Likelihoods for pairs  ~~~~~~~~")
  call timestamp()
  print *, "~~~~~~~  (n = ", NP, ") ~~~~~~~~"
  print *, ""
  
  call getPairLL(PairsFileName, NP)
endif

!====================================
if (quiet < 1)  print *, ""
call printt("Done.")
call DeAllocAll()



contains
  !====================================
  ! set default values
  !====================================
  subroutine set_defaults()
    quiet = 0  ! 0=not quiet, 1=quiet, -1=verbose
    dupQuiet = .TRUE.

    PedigreeFileName  = "NoFile"
    PairsFileName = "NoFile"
    OutFileName = "NoFile"
    AgePriorFileName = "AgePriors.txt"
    OnlyListFileName = 'NoFile'
    AF_FileName = 'NoFile'
    mt_FileName = 'NoFile'   
    GenoFormat = 'SEQ'                   

    ! MaxSibIter = 42  ! deprecated
    FindMaybe = -1
    nAmbMax = -99
    ResumePedC = "XX"
    ResumePed = -99
    Hermaphrodites = -99

    DoDup = .FALSE.
    DoPar = .FALSE.
    DoSibs = .FALSE.
    DoPairs = .FALSE.
    DoReadParents = .TRUE.
    MaybePO_onlyOH = .FALSE.
    NotDup = .FALSE.
    withAssignmentLog = .FALSE.
    
    Er = 0.001
    ErrFlavour = '2.9'
    ErV = 0D0   ! default: use Er + ErrFlavour
    MaxMismatchDup = -1
    MaxOppHom = -1
    MaxMendelE = -1
    TF = -2.0
    TA = 0.5

    GenoFileName = 'Geno.txt'
    LifehistFileName = 'LifeHist.txt'
    DumPrefix = (/'F', 'M', 'H'/)
    maxSibSize = 100
    Complx = 2
    Hermaphrodites = 0
    AgeEffect = 1
    FindMaybe = 0
    CalcLLR = 1   
    mxCP = 50  
  end subroutine set_defaults
  
  !====================================
  ! get filename of specs file
  !====================================
  subroutine get_SpecsFileName()  
    integer :: nArg, i, x
    character(len=32) :: arg
    logical :: file_exists

    SpecsFileName = 'SequoiaSpecs.txt'

    nArg = command_argument_count() 
    i = 0
    do x = 1, nArg
      i = i+1
      if (i > nArg)  exit
      call get_command_argument(i, arg)
      
      if (arg == '--specs') then
        i = i+1
        call get_command_argument(i, SpecsFileName)
        inquire(file=SpecsFileName, exist = file_exists)
        if (.not. file_exists) then
          print *, '--specs: file '//trim(SpecsFileName)//' not found!'
          stop
        endif  
      endif
    enddo
  end subroutine get_SpecsFileName 
  
  !====================================
  ! read arguments from file
  !====================================
  subroutine ReadSpecs()
    use sqa_fileIO, ONLY: FileNumRow
    character(len=200) :: tag
    character(len=nchar_filename)  :: tagvalue
    integer :: x, ntags

    ntags = FileNumRow(SpecsFileName)

    open (unit=101, file=trim(SpecsFileName), status="old")
    do x=1, ntags
      read(101, *)  tag, tagvalue
      
      select case (tag)
        case ('GenoFile', 'Genofile')
          GenoFileName = tagvalue
        case ('LHfile', 'LifeHist', 'LifeHistFile')
          LifehistFileName = tagvalue
        case ('AgePriorFile', 'AgePrior', 'ageprior')
          AgePriorFileName = tagvalue
        case ('GenotypingErrorRate', 'Err', 'ErrorRate')
          read(tagvalue(1:20), '(f20.0)')  Er   
        case ('MaxMismatchDUP', 'Max_dif_dup')
          read(tagvalue(1:20), '(i20)') MaxMismatchDup
        case ('MaxMismatchOH', 'Max_OH_PO')
          read(tagvalue(1:20), '(i20)') MaxOppHom
        case ('MaxMismatchME', 'Max_ME_PPO')
          read(tagvalue(1:20), '(i20)') MaxMendelE
        case ('MaxMismatch')  ! used in version 1.3; both incorrect; now calculated in R
          read(tagvalue(1:20), '(i20)') MaxMismatchDup    
          MaxOppHom = MaxMismatchDup - FLOOR(-nSNP * Er)  
          MaxMendelE = 3*MaxOppHom  
        case ('Tfilter')
          read(tagvalue(1:20), '(f20.0)') TF
        case ('Tassign')
          read(tagvalue(1:20), '(f20.0)') TA
        case ('MaxSibshipSize', 'Max_n_offspring')
          read(tagvalue(1:20), '(i20)') maxSibSize
        case ('DummyPrefixFemale')
          read(tagvalue(1:20), '(a20)')  DumPrefix(1) 
        case ('DummyPrefixMale')
          read(tagvalue(1:20), '(a20)')  DumPrefix(2)
        case ('DummyPrefixHerm')
          read(tagvalue(1:20), '(a20)')  DumPrefix(3)
        case ('Complexity', 'Complex', 'Complx')
          read(tagvalue(1:2), '(i2)') Complx
        case ('Hermaphrodites', 'Herm')
          read(tagvalue(1:2), '(i2)') Hermaphrodites
        case ('UseAge', 'AgeEffect')
          read(tagvalue(1:2), '(i2)') AgeEffect
        case ('FindMaybeRel', 'GetMaybeRel', 'MaybeRel')
          read(tagvalue(1:2), '(i2)') FindMaybe
        case ('CalcLLR')
          read(tagvalue(1:2), '(i2)') CalcLLR
        case ('ErrFlavour')
          read(tagvalue(8:10), '(a3)')  ErrFlavour  
        case ('MaxCandParents')
          read(tagvalue(1:20), '(i20)') mxCP     
        case default
          ! ignore the rest       
      end select
    enddo
    close(101)

    ! backwards compatability (version 2.1 & 2.2)
    if (Complx > 2) then    
      select case (Complx)
        case (4, 41)
          Hermaphrodites = 1
        case (5, 51)
          Hermaphrodites = 2
        case default
          Hermaphrodites = 0  
      end select
      if (Complx == 4 .or. Complx == 5)  Complx = 2
      if (Complx == 41 .or. Complx == 51)  Complx = 1
    endif

  end subroutine ReadSpecs


  !====================================
  ! read command line arguments
  !====================================
  subroutine read_args()
    use sqa_fileIO, only: valid_formats    ! genotype formats: SEQ, PED, RAW, LMT
    integer :: nArg, i, x, z
    character(len=32) :: arg, argOption

    nArg = command_argument_count()

    if (nArg == 0) then
      print *, "please provide at least 1 argument"
      print *, ""
      call print_help()
      stop
    endif

    i = 0
    do x = 1, nArg
      i = i+1
      if (i > nArg)  exit
      call get_command_argument(i, arg)
      
      select case (arg)
        case ('-v', '--version')
          print '(2a)', 'version ', version
          stop
        
        case ('-h', '--help')
          call print_help()
          stop
        
        case ('--dup')
          DoDup = .TRUE.
          
        case('--nodup')
          NotDup = .TRUE.
          
        case ('--par')
          DoPar = .TRUE.  
          
        case ('--ped')
          DoSibs = .TRUE. 

        case('--nopar')
          DoReadParents = .FALSE.
          
        case ('--maybePO')
          FindMaybe(1) = 1
          call get_command_argument(i+1, argOption)
          if (Len_Trim(argOption) > 0 .and. argOption(1:2)/="--") then   ! optional argument
            i = i+1
            read(argOption(1:10), '(i10)')  nAmbMax(1)   ! default: 7*nInd
            if (nAmbMax(1) == 0)   FindMaybe(1) = 0
          endif
         
        case ('--maybeRel')
          FindMaybe(2) = 1
          call get_command_argument(i+1, argOption)
          if (Len_Trim(argOption) > 0 .and. argOption(1:2)/="--") then
            i = i+1
            read(argOption(1:10), '(i10)')  nAmbMax(2)   
            if (nAmbMax(2) == 0)   FindMaybe(2) = 0
          endif
           
        case ('--geno')  
          i = i+1
          call get_command_argument(i, GenoFileName)
          
       case ('--informat', '--inFormat', '--genoformat', '--genoFormat')
          i = i+1
          call get_command_argument(i, GenoFormat)
          if (.not. any(valid_formats == GenoFormat)) then
            print *, 'ERROR: informat must be one of: ', valid_formats
            stop
          endif
          
        case ('--lifehist')    
          i = i+1
          call get_command_argument(i, LifehistFileName)
          
        case ('--ageprior', '--agepriors')    
          i = i+1
          call get_command_argument(i, AgePriorFileName)
          
        case ('--pedigreeIN')  
          i = i+1
          call get_command_argument(i, PedigreeFileName)
          
        case ('--resume')  
          DoSibs = .TRUE.        
          i = i+1
          call get_command_argument(i, argOption)
          read(argOption(1:2), '(i2)')  ResumePed
          if (ResumePed > 0) then
            write(ResumePedC, '(i2.2)') ResumePed   ! trailing 0
            inquire(file="Pedigree_round"//ResumePedC//".txt", exist = FileOK)
            if (.not. FileOK) then
              write(*,*)  "--resume: file ", "Pedigree_round"//ResumePedC//".txt", " not found"
              stop
            endif
          else 
            inquire(file="Pairs_01.txt", exist = FileOK)
            if (.not. FileOK) then
              write(*,*)  "--resume 0: file Pairs_01.txt not found"
              stop
            endif
         endif
         
        case ('--only')
          i = i+1
          call get_command_argument(i, OnlyListFileName)
          
        case ('--err')
          i = i+1
          call get_command_argument(i, argOption)
          read(argOption, *)  Er   
          if (Er <= 0.0 .or. Er > 0.5)  stop 'please provide a genotyping error rate --err >0 and <=0.5'        

        case ('--errV')
          do z=1,3
            i = i+1
            call get_command_argument(i, argOption)
            if (Len_Trim(argOption) == 0 .or. argOption(1:2)=="--") then
              stop '--errV requires 3 values: hom|other hom, het|hom, and hom|het'
            endif
            read(argOption, *)  ErV(z)
            if (ErV(z) <= 0.0 .or. ErV(z) > 0.5)  stop 'Genotyping error rates must be >0 and <=0.5'
          enddo
          
        case ('--Tfilter')
          i = i+1
          call get_command_argument(i, argOption)
          read(argOption, *)  TF   
          if (TF <= -50.0 .or. TF > 0.0)  stop 'please provide a --Tfilter >-50 and <=0'
          
        case ('--Tassign')
          i = i+1
          call get_command_argument(i, argOption)
          read(argOption, *)  TA   
          if (TA <= 0.0 .or. TA > 10.0)  stop 'please provide a --Tassign >0 and <=10'     
                           
        case ('--maxsibsize')
          i = i+1
          call get_command_argument(i, argOption)
          read(argOption(1:10), '(i10)')  maxSibSize
          
        case ('--maxCP')
          i = i+1
          call get_command_argument(i, argOption)
          read(argOption(1:10), '(i10)')  mxCP   ! default: 50 (ReadSpecs())
          
        case ('--dummyprefix')
          i = i+1
          call get_command_argument(i, argOption)
          read(argOption(1:20), '(i2)')  DumPrefix(1)
          i = i+1
          call get_command_argument(i, argOption)
          if (Len_Trim(argOption) == 0 .or. argOption(1:2)=="--") then
            stop '--dummyprefix requires 2 values: prefix for females and for males'
          else
            read(argOption(1:20), '(i2)')  DumPrefix(2)
          endif
          i = i+1
          call get_command_argument(i, argOption)
          if (Len_Trim(argOption) /= 0 .and. argOption(1:2)/="--") then 
            read(argOption(1:20), '(i2)')  DumPrefix(3)   ! optional: prefix for hermaphrodites        
          else
            i = i-1
          endif
          
        case ('--af', '--maf', '--freq')
          i = i+1
          call get_command_argument(i, AF_FileName)
          
         case ('--mt')
          i = i+1
          call get_command_argument(i, mt_FileName)
          
        case ('-o', '--out')
          i = i+1
          call get_command_argument(i, OutFileName)
          
        case ('--pairs')
          DoPairs = .TRUE.
          i = i+1
          call get_command_argument(i, PairsFileName)
          
        case ('--complex') 
          i = i+1
          call get_command_argument(i, argOption)
          select case (argOption)
            case ('mono')
              Complx = 0
            case ('simp') 
              Complx = 1
            case ('full')
              Complx = 2
            case default
              print '(2a, /)', '--Complex must be "mono", "simp" or "full", got: ', argOption
              call print_help()
              stop  
          end select
          
        case ('--herm')    
          i = i+1
          call get_command_argument(i, argOption)  
          select case (argOption)
            case ('no')
              hermaphrodites = 0
            case ('A') 
              hermaphrodites = 1
            case ('B')
              hermaphrodites = 2
            case default
              print '(2a, /)', '--herm must be "no", "A" or "B", got: ', argOption
              call print_help()
              stop  
          end select
                  
        case ('--age')
          i = i+1
          call get_command_argument(i, argOption)
          select case (argOption)
            case ('no')
              AgeEffect = 0
            case ('yes') 
              AgeEffect = 1
            case ('extra')
              AgeEffect = 2
            case default
              print '(2a, /)', '--age must be "no", "yes" or "extra", got: ', argOption
              call print_help()
              stop  
          end select
          
         case ('--noLLR')
          CalcLLR = 0
          
        case ('--log')
          withAssignmentLog = .TRUE.
          
        case ('--quiet')
          if (quiet /= 0) then
            write(*,*)  "You can't specify both quiet & verbose!"
            stop
          endif
          quiet = 1
          
        case ('--verbose')
          if (quiet /= 0) then
            write(*,*)  "You can't specify both quiet & verbose!"
            stop
          endif
          quiet = -1
        
        case('--specs') 
          i = i+1
          ! handled by get_SpecsFileName()  
          
        case default
            print *, ''
            print *, 'Unrecognised command-line argument: ', arg 
            print *, 'See --help for valid arguments'
            !call print_help()
            stop

      end select
    end do
    
    ! check for illogical combinations
    if (DoPar .and. DoSibs .and. .not. DoReadParents) then
      write(*,*)  "Cannot combine --nopar with --par --ped"
      stop
    endif

    if (PedigreeFileName/='NoFile' .and. .not. DoReadParents) then
      write(*,*)  "Cannot combine --pedigreeIN with --nopar;" 
      write(*,*)  "parents.txt is not used when --pedigreeIN specifies a different file."
      stop
    endif

    if (DoSibs .and. .not. DoPar .and. PedigreeFileName/='NoFile') then
      DoReadParents = .FALSE.   ! else read in twice
    endif

    if (DoPairs .and. (DoDup .or. DoPar .or. DoSibs) .and. .not. ANY(FindMaybe==1)) then
      write(*,*)  "Cannot combine --pairs with --dup/--par/--ped"
      stop
    endif

    if ((DoPar .or. DoSibs) .and. .not. NotDup) then
      DoDup = .TRUE.
    endif
    if (DoDup) then
      dupQuiet = quiet == 1  
    endif
    

    ! ensure consistency between parameters
    if (ALL(FindMaybe == -1)) then
      if (DoSibs) then
        FindMaybe(1) = 0
        FindMaybe(2) = FindMaybeX  ! from SequoiaSpecs.txt
      else
        FindMaybe(1) = FindMaybeX
        FindMaybe(2) = 0
      endif
    endif

    if (any(FindMaybe==1) .and. CalcLLR==0) then
      if (FindMaybe(2)<1) then
        MaybePO_onlyOH = .TRUE.
      else 
        write(*,*)  "Cannot combine --maybeRel with --noLLR; did you mean --maybePO?"
        stop
      endif
    endif

    if (CalcLLR == 1) then 
      if ((DoDup .or. any(FindMaybe==1) .or. DoPairs) .and. .not. (DoPar .or. DoSibs)) then
        CalcLLR = 0
      endif
    endif

  end subroutine read_args

  !====================================
  ! print overview of command line arguments
  !====================================
  subroutine print_help()
      print '(a, /)', 'command-line options:'
      print '(a)',    ' general:'
      print '(a)',    '  -v, --version       print version information and exit'
      print '(a)',    '  -h, --help          print usage information and exit'
      print '(a)',    '  --dup               check for potential duplicates; default TRUE if --par or --ped'
      print '(a)',    '  --nodup             no check for potential duplicates'
      print '(a)',    '  --par               parentage assignment'
      print '(a)',    '  --nopar             do not read parents.txt prior to pedigree reconstruction' 
      print '(a)',    '  --ped               full pedigree reconstruction: sibship clustering,', &
                      '                        grandparents, and possibly more parents'
      print *,''                
      print '(a)',    ' input/output files: (use NoFile to not use the default)'
      print '(a)',    '  --geno <filename>   file with genotype data, file extension will be added based on',&
                    '                         --informat. Default: Geno'
      print '(a)',    '  --genoformat <x>   SEQ: no header, 0/1/2, missing -9, IDs in column 1; .txt (default)', &  
                    '                       PED: no header, 2 columns per SNP coded 1/2 or A/C/T/G, missing 0,',&
                    '                       IDs in column 2 of 6 non-SNP columns; .ped (+.map)', &
                    '                       RAW: header, 0/1/2/NA, IDs in column 2 of 6 non-SNP columns; .raw', &
                    '                       LMT: no header, 0/1/2 without spacing, IDs in separate file;', &
                    '                        .geno + .id. values 3-9 interpreted as missing'
      print '(a)',    '  --lifehist <filename>    file with lifehistory data, with columns id, sex, birthyear,', &
                      '                      and optionally BY.min, BY.max, YearLast. Default: LifeHist.txt'
      print '(a)',    '  --ageprior <filename>    file with agepriors (columns M-P-FS-MHS-PHS); AgePriors.txt'        
      print '(a)',    '  --pedigreeIN <filename>  starting point for further pedigree reconstruction', & 
                      '                        (<par>, <ped>), or is conditioned on (<maybe>, <pairs>)', & 
                      '                         or for which to calculate parental LLRs (if no other', &
                      '                         options given).'   
      print '(a)',    '  --af <filename>     optional input file with allele frequencies. Either',&
                      '                        1 column and no header, or multiple columns with column',&
                      '                        MAF, AF, or Frequency. E.g. output from plink --freq.'   
      print '(a)',    '  --only <filename>   only find parents for this subset of individuals;',& 
                      '                        their parents in --pedigreeIN are dropped; ',&
                      '                        when combined with --ped other siblings are in DummyParents.txt.',&
                      '                        IDs in column 1, or 2 columns with FID (ignored) + IID'   
      print '(a)',    '  --mt <filename>     optional input file with mitochondrial haplotype sharing'
      print '(a)',    '  --out <filename>    output pedigree file name'  
      print '(a)',    '  --specs <filename>  file with parameter and filename specifications. Default: SequoiaSpecs.txt' 
      
      print *,''
      print '(a)',    ' parameterisation:'
      print '(a)',    '  --complex <option>  breeding system complexity, "mono", "simp" or "full"'
      print '(a)',    '  --herm <option>     hermaphrodites; "no", "A" (distinguish between sex roles)',&
                      '                        or "B" (indiscriminate regarding sex roles)'
      print '(a)',    '  --age <option>      weight of age in assignments: "no", "yes", or "extra"'
      
      print '(a)',    '  --err <value>        presumed genotyping error rate, >0 & <0.5'
      print '(a)',    '  --errV <3 values>   alternative to --err: P(observed|actual) for',&
                      '                        hom|other hom, het|hom, and hom|het' 
      print '(a)',    '  --Tfilter <value>   filtering threshold, <0. Default -2.0' 
      print '(a)',    '  --Tassign <value>   assignment threshold, >0. Default 0.5'          
      print '(a)',    '  --maxsibsize        max size of half-sib group (# offspring per parent)'            
      print *,''
      print '(a)',    ' other:'
      print '(a)',    '  --quiet             suppress (almost) all messages'
      print '(a)',    '  --verbose           print extra many messages'  
      print '(a)',    '  --noLLR             do not calculate parental LLR'
      print '(a)',    '  --dummyprefix <dams> <sires>  prefix for used for dummy parents. Default F,M'  
      print '(a)',    '  --log               write log with assignment steps for each individual'
      print '(a)',    '  --resume <x>        resume pedigree reconstruction at round <x>'
      print '(a)',    '  --maybePO <max>     check for potential (remaining) parent-offspring pairs', &
                      '                        optionally set maximum number of pairs'        
      print '(a)',    '  --maybeRel <max>    check for potential (remaining) relative pairs'
      print '(a)',    '  --pairs <filename>  calculate for each pair LLs for 7 relationships. ', &
                      '                        Can only be combined with --pedigreeIN and --quiet'        
      print '(a)',    ''
      print '(a)',    '  settings specified on the command line will overrule settings in SequoiaSpecs.txt'
      print '(a)',    ''
      print '(a)',    '  When multiple arguments are specified, the execution order is: ', &
                      '  [dup] .. [pedigreeIN].. [par] .. [ped] .. [LLR] .. [maybe]'
  end subroutine print_help
    


  !====================================
  ! initiate all arrays
  !====================================
  subroutine Initiate()
    use sqa_fileIO, ONLY: read_geno
    use sqa_general, ONLY: mk_OcA, InheritanceProbs
    use OHfun, ONLY: init_OH
    use binom, ONLY: calc_MaxMismatch

    integer, allocatable :: BYrange(:,:)                                          
    double precision, allocatable :: AP_IN(:,:), AF(:)
    character(len=nchar_ID), allocatable :: SNP_names(:)  ! not used, output from read_geno()
    integer :: i

    zero = 0.0D0   ! compiler errors when setting it as global parameter

    ! genotype data  ~~~
    if (quiet < 1)  print *, ""
    if (quiet < 1)  call printt("Reading genotypes in "//trim(GenoFileName)//" ... ")
    call read_geno(Geno=Genos, ID=Id, SNP_names=SNP_names, FileName=GenoFileName, FileFormat=GenoFormat)
    nSnp = SIZE(Genos, DIM=1)
    nInd = SIZE(Genos, DIM=2) -1
    allocate(ID_tmp(0:nInd))
    ID_tmp(0) = 'NA'
    ID_tmp(1:nInd) = ID
    call move_alloc(ID_tmp, ID)    

    ! get max ID length (for output files)
    ID_len = MAX(LEN_TRIM(DumPrefix(1)), LEN_TRIM(DumPrefix(2)), LEN_TRIM(DumPrefix(3))) +4   ! min. for writing dummy names
    do i=1,nInd
      if (LEN_TRIM(Id(i)) > ID_len)  ID_len = LEN_TRIM(Id(i))
    enddo

    ! lifehistory data  ~~~
    allocate(BY(nInd))  
    BY=-999
    allocate(BYrange(nInd,2))
    BYrange=-999            
    allocate(Sex(nInd))
    Sex = 3
    allocate(YearLast(nInd))
    YearLast=999   ! NOTE: POSITIVE!

    if (LifehistFileName /= 'NoFile') then
      if (quiet < 1)  call printt("Reading life history data in "//trim(LifehistFileName)//" ... ")
      call ReadLifeHist(LifehistFileName, BYrange)
    endif

    if (any(YearLast /= 999)) then
      AnyYearLast = .TRUE.
    else
      AnyYearLast = .FALSE.
      deallocate(YearLast)
    endif

    if (any(Sex==4) .and. hermaphrodites==0) then
      call Erstop("Please set --herm A or --herm B, found hermaphrodites (Sex=4) in Lifehist file", .FALSE.)
    endif

    ! ageprior  ~~~
    allocate(AP_IN(MaxMaxAgePO, 5))
    AP_IN = 1D0
    if (AgePriorFileName /= 'NoFile') then
      if (quiet < 1)  call printt("Reading age priors in "//trim(AgePriorFileName)//" ... ")
      call ReadAgePrior(AgePriorFileName, AP_IN)
      if (quiet < 1)  call printt("Prepping age data ... ")
    endif
    call PrepAgeData(AP_IN, BYrange)
    deallocate(BYrange)
    deallocate(AP_IN)

    ! general prep  ~~~
    if (quiet < 1)  call printt("generating look-up tables with probabilities ... ")
    AF = getAF(AF_FileName)
    if (erV(1) > TINY(0D0)) then
     OcA = mk_OcA(erV)
    else
      OcA = mk_OcA(Er, ErrFlavour)
    endif
    call InheritanceProbs(nSnp, AF, OcA) 
    if (MaxMismatchDup < 0) MaxMismatchDup = calc_MaxMismatch(0.9999**(1d0/nInd), AF, OcA(:,2:4), 'DUP')
    if (MaxOppHom < 0)      MaxOppHom      = calc_MaxMismatch(0.9999**(1d0/nInd), AF, OcA(:,2:4), 'PO ')
    if (MaxMendelE < 0)     MaxMendelE     = calc_MaxMismatch(0.9999**(1d0/nInd), AF, OcA(:,2:4), 'PPO')
    deallocate(AF)
    call fixbounds_inheritprobs()  ! old (sequoia) vs new (sqa_general) 1:3 vs 0:2 for actual genotype
    call precalc_quicksibs()
    call rchkusr()

    if (quiet < 1)  call printt("Allocating arrays ... ")
    call AllocArrays()

    call init_OH()

    if (mt_FileName /= 'NoFile') then
      DoMtDif = .TRUE.
      call ReadMt(mt_FileName)
    else
      DoMtDif = .FALSE.
    endif

    ! parents ~~~
    ! (pedigree prior or prev. assigned parents)
    allocate(DummyNamesIO(nInd, 2))
    DummyNamesIO = "NA"
    if (trim(PedigreeFileName)/= "NoFile") then
      call ReadPedFile(PedigreeFileName)
      if (quiet<1)  call timestamp()
      if (quiet<1)  print *, " # parents: ", COUNT(Parent /= 0, DIM=1)
      if (quiet<1 .and. any(nC>0))  call timestamp()
      if (quiet<1 .and. any(nC>0))  print *, " # dummies: ", nC
    endif

    ! only-those-individuals  ~~~
    if (trim(OnlyListFileName)/= "NoFile") then
      call ReadOnlyList(OnlyListFileName)
      if (quiet<1)  call timestamp()
      if (quiet<1)  print *, "--only: ", COUNT(.not. skip), " individuals out of ", nInd
    endif
    
    ! helper array to print progress every 5%, 10%, ... of individuals
    viginti =  mk_seq(nInd, 20)
    
  end subroutine Initiate
  
  function getAF(FileName)  result(AF)
    use sqa_fileIO, ONLY: readAF, ilong    
    character(len=*), intent(IN) :: FileName
    double precision, allocatable :: AF(:)
    integer :: l
    
    allocate(AF(nSnp))
    if (FileName == 'NoFile') then
      AF = 1D0
      do l=1,nSnp
        if (ALL(Genos(l,:)==-1)) cycle
        AF(l) = dble(SUM(int(Genos(l,1:),kind=ilong), MASK=Genos(l,1:)/=-1))/(COUNT(Genos(l,1:)/=-1)*2)
      enddo  
    else
      if (quiet ==-1)  call printt("Reading allele frequencies in "//trim(FileName)//" ... ")
      AF = readAF(FileName)
    endif

  end function getAF

end program Main

! ####################################################################

subroutine fixbounds_inheritprobs
use Global, ONLY: nSnp, AHWE, OcA, AKAP, OKAP, AKA2P, OKA2P
implicit none

double precision, allocatable :: TMP2(:,:), TMP3(:,:,:)

allocate(Tmp2(1:3, nSnp))
Tmp2 = AHWE
call move_alloc(Tmp2, AHWE)

allocate(Tmp2(1:3, -1:2))
Tmp2 = OcA
call move_alloc(Tmp2, OcA)

allocate(Tmp3(1:3, 1:3,nSnp))
Tmp3 = AKAP
call move_alloc(Tmp3, AKAP)

allocate(Tmp3(-1:2, 1:3,nSnp))
Tmp3 = OKAP
call move_alloc(Tmp3, OKAP)

allocate(Tmp3(1:3,1:3,1:3))
Tmp3 = AKA2P
call move_alloc(Tmp3, AKA2P)

allocate(Tmp3(-1:2,1:3,1:3))
Tmp3 = OKA2P
call move_alloc(Tmp3, OKA2P)

end subroutine fixbounds_inheritprobs

! ####################################################################

subroutine precalc_quicksibs()
use Global, ONLY: nSnp, PHS, PFS, AHWE, OKAP, OKA2P, OHWE
implicit none 

integer :: l,i,j,m,h
double precision :: Tmp1(3), Tmp2(3,3)

allocate(PHS(-1:2,-1:2,nSnp))
allocate(PFS(-1:2,-1:2,nSnp))
PHS = 1D0
PFS = 1D0
do l=1,nSnp
  do i=0,2  ! obs offspring 1
    do j=0,2    ! obs offspring 2
      Tmp1=0D0
      Tmp2=0D0
      do m=1,3    !act shared parent 
        Tmp1(m) = OKAP(i,m,l) * OKAP(j,m,l) * AHWE(m,l)
        do h=1,3
          Tmp2(m,h) = OKA2P(i,m,h) * OKA2P(j,m,h) * AHWE(m,l) *AHWE(h,l)
        enddo
      enddo
      PHS(i,j,l) = SUM(Tmp1) / (OHWE(i,l) * OHWE(j,l))
      PFS(i,j,l) = SUM(Tmp2) / (OHWE(i,l) * OHWE(j,l))
    enddo
  enddo
enddo

end subroutine precalc_quicksibs

! ####################################################################

subroutine AllocArrays()
use Global
implicit none

integer :: i,l,h

allocate(Parent(nInd,2))
Parent = 0
allocate(Lind(nInd))
Lind = 0D0

nC = 0 
if (DoSibs) then
  allocate(nS(nInd/2,2))
  ns = 0
  allocate(SibID(maxSibSize, nInd/2, 2))
  SibID = 0 
  allocate(GpID(2, nInd/2,2))
  GpID = 0
  allocate(CLL(nInd/2,2))
  CLL = missing
  allocate(IsNewSibship(nInd/2, 2))
  IsNewSibship = .TRUE.
endif
allocate(ToCheck(nInd))
ToCheck = .FALSE.
allocate(SelfedIndiv(nInd))
SelfedIndiv = .FALSE.
allocate(skip(nInd))  ! from --only
skip = .FALSE.

allocate(NFS(nInd))
NFS = 1
allocate(FSID(MaxSibSize+1, nInd))
FSID = 0
FSID(1, :) = (/ (i, i=1, nInd) /)
FSID(MaxSibSize+1, :) = (/ (i, i=1, nInd) /)    ! 'primary' sib

allocate(Mate(nInd))
Mate = 0  
if (DoSibs) then
  allocate(DumMate(nInd/2,2))
  DumMate = 0  
  allocate(DumClone(nInd/2,2))
  DumClone = 0  
endif

allocate(LindX(3,nSnp,nInd))  ! used when missing genotype & at start
if (DoSibs) then
  allocate(DumP(3,nSnp, nInd/2,2))
  allocate(XPr(3,3,nSNP, nInd/2,2))
  XPr = 1D0
endif
do l=1,nSnp
  do h=1,3
    LindX(h,l,:) = AHWE(h,l)  
    if (DoSibs) then
      DumP(h,l,:,:) = AHWE(h,l)
      XPr(2,h,l,:,:) = AHWE(h,l)  ! GP contribution
    endif
  enddo
enddo

end subroutine AllocArrays

! ####################################################################

subroutine CheckMono
use Global
implicit none

integer :: i, j, k, nOff, Off(maxSibSize), sxOff(maxSibSize), s, par
logical :: ParOK

nOff = 0
Off = 0
sxOff = 3
par = 0
ParOK = .TRUE.

do i=1, nInd
  do k=1,2
    call getOff(i, k, .FALSE., nOff, Off, sxOff)
    if (nOff == 0)  cycle
    do j=1, nOff
      if (Parent(Off(j),3-k) /= 0) then
        Mate(i) = Parent(Off(j),3-k)
        exit
      endif
    enddo
    if (Mate(i)==0) then
      call NewSibship(Off(1), 0, 3-k)   ! sets Mate(i) = -nC(3-k)
    endif
    do j=1, nOff
      if (Parent(Off(j),3-k) /= Mate(i) .and. Parent(Off(j), 3-k)/=0) then
        call Erstop("Please change to Complex='simp', assigned parents suggest non-monogamy", .FALSE.)
      else if (Parent(Off(j),3-k) == 0) then
        call ChkValidPar(Off(j), sxOff(j), Mate(i), 3-k, ParOK)
        if (ParOK) then
          call SetPar(Off(j), sxOff(j), Mate(i), 3-k)  ! make all half-sibs full sibs
        else 
          call SetPar(Off(j), sxOff(j), 0, k)
        endif
      endif
    enddo
    ! dummy offspring
    call getOff(i, k, .TRUE., nOff, Off, sxOff)
    do j=1, nOff
      if (Off(j) >0) cycle
      if (GpID(3-k, -Off(j), sxOff(j)) /= Mate(i) .and. GpID(3-k, -Off(j), sxOff(j)) /= 0) then
        call Erstop("Please change to Complex='simp', assigned (grand)parents suggest non-monogamy", .FALSE.)
      endif
      call ChkValidPar(Off(j), sxOff(j), Mate(i), 3-k, ParOK)
      if (ParOK) then
        call setPar(Off(j), sxOff(j), Mate(i), 3-k) 
      else
        call setPar(Off(j), sxOff(j), 0, k)   
      endif
    enddo
  enddo
enddo

if (quiet==-1 .and. any(Mate/=0)) then
  call Rprint("Added dummy parents to ensure monogamy...", (/0/), (/0.0D0/), "NON")
endif     

do k=1,2
  if (nC(k) == 0)  cycle
  do s=1, nC(k)
    call getFSpar(s, k, .TRUE., par)
    if (par /= 0) then
      DumMate(s,k) = par
      if (any(Parent(SibID(1:ns(s,k),s,k), 3-k) == 0)) then
        do j=1, ns(s,k)
          i = SibID(j,s,k)
          call ChkValidPar(i, 3, Par, 3-k, ParOK)
          if (ParOK) then
            call setPar(i, 3, Par, 3-k)
          else
            call setPar(i, 3, 0, k)  ! remove from sibship
          endif
        enddo
      endif
    else if (any(Parent(SibID(1:ns(s,k),s,k), 3-k) /= 0)) then
      call Erstop("Please change to Complex='simp', assigned (dummy) parents suggest non-monogamy", .FALSE.)
    else
      call NewSibship(SibID(1,s,k), SibID(2,s,k), 3-k)   ! also works if ns=1, and SibID(2)=0
      if (ns(s,k) > 2) then
        do j=3, ns(s,k)
          call setPar(SibID(j,s,k), 3, -nC(3-k), 3-k)         
        enddo
      endif
      ToCheck(SibID(1:ns(s,k),s,k)) = .TRUE.
      DumMate(s,k) = Parent(SibID(1,s,k), 3-k)
    endif
    ! dummy offspring
    call getOff(i, k, .TRUE., nOff, Off, sxOff)
    do j=1, nOff
      if (Off(j) >0) cycle
      if (GpID(3-k, -Off(j), sxOff(j)) /= DumMate(s,k) .and. GpID(3-k, -Off(j), sxOff(j)) /= 0) then
        call Erstop("Please change to Complex='simp', assigned (grand)parents suggest non-monogamy", .FALSE.)
      endif
      call ChkValidPar(Off(j), sxOff(j), DumMate(s,k), 3-k, ParOK)
      if (ParOK) then
        call setPar(Off(j), sxOff(j), DumMate(s,k), 3-k) 
      else
        call setPar(Off(j), sxOff(j), 0, k)   
      endif
    enddo  
  enddo
enddo

end subroutine CheckMono

! ####################################################################

subroutine duplicates(dupQuiet)
use Global
implicit none

logical, intent(IN) :: dupQuiet
integer :: i, j, l, CountMismatch, nDupGenoID, nDupGenos, dupGenoIDs(nInd,2), &
  DupGenos(nInd,2), nMismatch(nInd), SnpdBoth(nInd)
double precision :: DupLR(nInd), LLtmp(2), LL(7), LLX(7)
character(len=200) :: HeaderFMT, DataFMT
integer :: IsBothScored(-1:2,-1:2), IsDifferent(-1:2,-1:2), SnpdBoth_ij

nDupGenoID = 0
do i=1,nInd-1
  do j=i+1, nInd
    if (Id(i) == Id(j)) then
      nDupGenoID = nDupGenoID + 1
      dupGenoIDs(nDupGenoID,1) = i
      dupGenoIDs(nDupGenoID,2) = j
    endif
  enddo 
enddo

!====================
! (nearly) identical genotypes?
nDupGenos = 0
DupGenos = -9
nMismatch = -9
SnpdBoth = -9             
DupLR = missing

IsBothScored = 1
IsBothScored(-1,:) = 0
IsBothScored(:,-1) = 0
IsDifferent = 0
IsDifferent(0, 1:2) = 1
IsDifferent(1, (/0,2/)) = 1
IsDifferent(2, 0:1) = 1

LLtmp = missing
LL = missing
LLX = missing
do i=1,nInd-1
  do j=i+1, nInd
    if (skip(i) .and. skip(j))  cycle
    CountMismatch=0
    SnpdBoth_ij = 0
    do l=1, nSnp
      SnpdBoth_ij = SnpdBoth_ij + IsBothScored(Genos(l,i), Genos(l,j))
      CountMismatch = CountMismatch + IsDifferent(Genos(l,i), Genos(l,j))
      if (CountMismatch > MaxMismatchDup)  exit
    enddo
    if (CountMismatch > MaxMismatchDup)  cycle
!    call CalcOppHom(i,j)  ! OH + LLR PO/U  
    LLtmp = missing
    call PairSelf(i, j, LLtmp(1))
    call CheckPair(i, j,3,7,LL, LLX) 
    LLtmp(2) = MaxLL(LL)
    if ((LLtmp(1) - LLtmp(2)) > TF)  then   
      nDupGenos = nDupGenos + 1
      DupGenos(nDupGenos,1) = i
      DupGenos(nDupGenos,2) = j
      nMisMatch(nDupGenos) = CountMismatch
      SnpdBoth(nDupGenos) = SnpdBoth_ij
      DupLR(nDupGenos) = LLtmp(1) - LLtmp(2)
    endif
    if (nDupGenos==nInd) then
      print *, ''
      print *, "reached max for duplicates!"
      print *, ''
      exit
    endif
    if (nDupGenos==nInd) exit
  enddo
  if (nDupGenos==nInd) exit
enddo

if (.not. dupQuiet) then
  call timestamp(.TRUE.)
  write(*, '("Found ", i4  ," pairs of potentially duplicated genotypes")') nDupGenos
  print *, ""
endif

!============================
! write to file
write(HeaderFMT, '( "(a15, 2(a8, 4X, a", I0, ", 4X), 3a10)" )')  ID_len
write(DataFMT,   '( "(a15, 2(i8, 4X, a", I0, ", 4X), 2i10, f10.2)" )')  ID_len

if (nDupGenoID > 0 .or. nDupGenos > 0 .or. .not. dupQuiet) then
  open (unit=201,file="DuplicatesFound.txt",status="replace")
  write (201, HeaderFMT) "Type", "Row1", "ID1", "Row2", "ID2", "nDiffer", "nBoth", "LLR"
  if (nDupGenoID>0) then
      do i=1,nDupGenoID
          write (201, DataFMT) "GenoID", dupGenoIDs(i,1), Id(dupGenoIDs(i,1)), &
          dupGenoIDs(i,2), Id(dupGenoIDs(i,2))
      enddo   
  endif
  if (nDupGenos>0) then
      do i=1,nDupGenos
          write (201, DataFMT) "Genotype", dupGenos(i,1), Id(dupGenos(i,1)), &
          dupGenos(i,2), Id(dupGenos(i,2)), nMismatch(i), SnpdBoth(i), DupLR(i)  
      enddo   
  endif
  close (201)
endif

if (nDupGenoID > 0) then
  call Erstop("Found duplicated IDs in geno file, please fix", .FALSE.)
endif

end subroutine duplicates
 
! ###################################################################### 

subroutine findambig(ParSib, nAmbMax, onlyOH)   
use Global
use OHfun
implicit none

integer, intent(IN) :: ParSib, nAmbMax  ! 1: PO only; 2: all relatives
logical, intent(IN) :: onlyOH
integer :: namb, AmbigID(nAmbMax, 2), ambigrel(nAmbMax,2), BYtmp(2), &
  ambigoh(nAmbMax), ntrio, trioIDs(nInd, 3), trioOH(nInd, 3), ID_len_ambig
double precision :: ambiglr(nAmbMax, 2), trioLR(nInd, 3) 
integer :: i, j, k, x, topX, Anci(2,mxA), Ancj(2,mxA), maybe, Lboth, &
  u,v, ncp, CandPar(mxCP), m, t
double precision :: LL(7), LLtmp(7,3), dLL, LRR(3), LLX(7), LLP(3)
character(len=2) :: RelName(9)
character(len=200) :: HeaderFMT_pairs, DataFMT_pairs, HeaderFMT_triads, DataFMT_triads
logical :: AncOK, ParOK(2)

nAmb = 0
AmbigID = 0
AmbigLR = missing
AmbigOH = -9
AmbigRel = 0   

if (.not. allocated(IndBY)) then   ! when all BY known 
  allocate(IndBY(nYears, nInd, 5))  ! year - indiv - own/wo/w dummy off+par
  IndBY = LOG10(1.0D0/nYears)
endif

t = 1
do i=1,nInd-1
  if (nAmb==nAmbMax)  exit   
  if (MODULO(i,200)==0)   call rchkusr()
  if (quiet==-1 .and. any(viginti==i)) call print_progress(i,t)
  call GetAncest(i,1,Anci)
  
  do j=i+1,nInd 
    if (skip(i) .and. skip(j))  cycle
    Lboth = COUNT(Genos(:,i)/=-1 .and. Genos(:,j)/=-1)  
    if (Lboth < nSnp/2.0)   cycle   ! >1/2th of markers missing
    if (ANY(Parent(i,:)==j) .or. ANY(Parent(j,:)==i)) cycle  ! PO
    if (ALL(Parent(i,:)/=0)) then
      if (Parent(i,1)==Parent(j,1) .and. Parent(i,2)==Parent(j,2)) cycle  ! FS
      if (ANY(Anci(:,3)==j) .and. ANY(Anci(:,4)==j)) cycle  ! double GP
    endif
    if (Parent(j,1)/=0 .and. Parent(j,2)/=0) then
      call GetAncest(j,1,Ancj)
      if (ANY(Ancj(:,3)==i) .and. ANY(Ancj(:,4)==i)) cycle  ! double GP
    endif    
    
    LLtmp = missing
    LL = missing
    topX = 0
    dLL = missing
    ParOK = .FALSE.
    if (ParSib <2 .or. (All(Parent(i,:)/=0) .and. ALL(Parent(j,:)/=0))) then  
    ! check if they're not PO only
      if (calcOH(i,j) > maxOppHom)  cycle      
      if (ParSib==1) then
        if (QLR_PO(i,j) < TA)  cycle
      else if (ParSib==2) then
        if (QLR_PO(i,j) < 2*TA)  cycle
      endif
    endif
    if (ParSib==1) then
      BYtmp(1:2) = BY((/i,j/))
      BY((/i,j/)) = -9
      call ChkValidPar(i,Sex(i), j,Sex(j), ParOK(1))
      call ChkValidPar(j,Sex(j), i,Sex(i), ParOK(2))
      if (.not. ANY(ParOK))  cycle
      if (.not. onlyOH) then
        if (ParOK(1))  call CheckPair(i, j, Sex(j), 1, LLtmp(:,1), LLX)
        if (ParOK(2))  call CheckPair(j, i, Sex(i), 1, LLtmp(:,2), LLX)
        do k=1,7  
          LL(k) = MaxLL(LLtmp(k,1:2)) 
        enddo
        call BestRel2(LL, topX, dLL)
        BY((/i,j/)) = BYtmp(1:2)
        if (topX==6 .or. topX==7) cycle   ! conditionally unrelated
        if (COUNT(Parent == 0) > 0.95*nInd .and. topX>2)  cycle  ! else will exceed nAmbMax
      endif
    
    else if (ParSib == 2) then 
      maybe = 0
      LRR = missing
      LRR(1) = QLR_PO(i,j)
      topX = 0
      do k=1,2 
        if (Parent(i,k)/=0 .and. Parent(i,k)==Parent(j,k)) cycle
        if (Parent(i,k)>0) then
            if (ANY(Parent(Parent(i,k), :)==j)) cycle
        else if (Parent(i,k)<0) then
            if (ANY(GpID(:, -Parent(i,k), k)==j)) cycle
        endif
        if (Parent(j,k)>0) then
            if (ANY(Parent(Parent(j,k), :)==i)) cycle
        else if (Parent(j,k)<0) then
            if (ANY(GpID(:, -Parent(j,k), k)==i)) cycle
        endif
        call PairQFS(i, j, LRR(2)) 
        call PairQHS(i, j, LRR(3))       
        maybe = 0
        do x=1,3
          if (LRR(x) > 2*TA .and. LRR(x) < missing)  maybe=1  
        enddo
        if (maybe==0)  cycle
        if (BY(i) > BY(j)) then
          call CheckPair(i, j, k, 7, LL, LLX)  
        else
          call CheckPair(j, i, k, 7, LL, LLX)
        endif
        call BestRel2(LL, topX, dLL) 
        if (COUNT(Parent == 0) > 0.95*nInd .and. topX>2) then
          cycle  ! else will exceed nAmbMax)
        else if (topX==6 .or. topX==7) then  ! .or. dLL(2)<TA   .or. topX==8
          maybe = 0
          cycle
        else
          exit
        endif
      enddo
      if (maybe==0) cycle
    endif
    
    nAmb = nAmb + 1
    AmbigID(nAmb, :) = (/i, j/)
    AmbigOH(nAmb) = calcOH(i,j)
    if (ParSib==1) then
      AmbigLR(nAmb,1) = QLR_PO(i,j)
      AmbigRel(nAmb,1) = 1
    else if (ParSib==2) then
      AmbigLR(nAmb,1) = MAXVAL(LRR, MASK=LRR<MaybeOtherParent)
      AmbigRel(nAmb,1) = MAXLOC(LRR, MASK=LRR<MaybeOtherParent, DIM=1)
    endif
    if (onlyOH) then
      AmbigRel(nAmb,2) = 8
      AmbigLR(nAmb,2) = 0D0
    else
      AmbigRel(nAmb,2) = TopX
      AmbigLR(nAmb,2) = dLL 
    endif
    if (nAmb==nAmbMax) then
      if (quiet<1) then
        call Rprint("WARNING - reached max for maybe-rel, truncated!",(/0/), (/0.0D0/), "NON")
      endif
      exit
    endif
  enddo
enddo

if(quiet<1)  call timestamp()
if(quiet<1)  print *, "found ", nAmb , " pairs"

!~~~~~~~~~~~~
! triads
ntrio = 0
trioIDs = 0
trioLR = missing 
trioOH = -9  

if (nAmb>1 .and. COUNT(AmbigRel(:,2) == 1) > 1) then
  if (quiet<1)  call printt("Checking for Parent-Parent-Offspring trios ... ") 
    
  do i=1, nInd
    if (MODULO(i,500)==0)   call rchkusr()
    if (ntrio == nInd) exit
    if (ANY(Parent(i,:)/=0) .and. Hermaphrodites==0) cycle
    if (ANY(Parent(i,:)>0) .and. Hermaphrodites>0) cycle  
    ncp = 0
    CandPar = 0  
    if ((COUNT(AmbigID(:,1) == i .and. AmbigRel(:,2) <3) + &   ! PO or FS
      COUNT(AmbigID(:,2) == i .and. AmbigRel(:,2) <3)) < 2) cycle
    do j=1, nAmb
      if (AmbigRel(j,2) >2)  cycle 
      if (.not. ANY(AmbigID(j,:) == i))  cycle
      if (ncp == mxCP) exit
      do m=1,2
        if (AmbigID(j,m) == i) then
          if (AgeDiff(i, AmbigID(j,3-m))<=0)  cycle  ! unknown=999. 
          call ChkAncest(AmbigID(j,3-m), sex(AmbigID(j,3-m)), i, sex(i), AncOK)
          if (.not. AncOK)  cycle
          ncp = ncp + 1
          CandPar(ncp) = AmbigID(j,3-m)
        endif
      enddo
    enddo
    
    if (ncp > 1) then   
      do u=1, ncp-1
        do v=u+1, ncp
          if (Sex(CandPar(u))<3 .and. Sex(CandPar(u))==Sex(CandPar(v))) cycle
          if (Sex(CandPar(u))==1 .or. Sex(CandPar(v))==2) then
            call CheckParentPair(i, Sex(i), CandPar( (/u,v/) ), LLP)
          else
            call CheckParentPair(i, Sex(i), CandPar( (/v,u/) ), LLP)
          endif
          if (LLP(3) < -TA .or. LLP(3)==missing)   cycle
          ntrio = ntrio +1
          trioIDs(ntrio,1) = i
          trioIDs(ntrio, 2:3) = CandPar( (/u,v/) )
          trioLR(ntrio,:) = LLP
          
          trioOH(ntrio,1) = calcOH(i, CandPar(u))
          trioOH(ntrio,2) = calcOH(i, CandPar(v))
          trioOH(ntrio,3) = CalcTrioErr(i, CandPar( (/u,v/) ) )
          if (ntrio == nInd) exit        
        enddo
        if (ntrio == nInd) exit
      enddo
    endif
  enddo
  if (quiet<1)  call timestamp()
  if (quiet<1)  print *, "found ", ntrio, " triads"
endif

  
RelName = (/ "PO", "FS", "HS", "GP", "FA", "HA", "U ", "XX", "X2" /)
ID_len_ambig = MAX(ID_len, 7)
write(HeaderFMT_pairs, '( "(2(a", I0, ", 4X), 2a5, 5a10, a5)" )')  ID_len_ambig
write(DataFMT_pairs, '( "(2(a", I0, ", 4X), 2i5, i10, a10, f10.2, a10, f10.2, i5)" )')  ID_len_ambig
write(HeaderFMT_triads, '( "(3(a", I0, ", 4X), 3a7, 3a12)" )')  ID_len_ambig
write(DataFMT_triads, '( "(3(a", I0, ", 4X), 3i7, 3(3X, f9.2))")')  ID_len_ambig

if (ParSib==1) then
  open (unit=201,file="Unassigned_relatives_par.txt",status="unknown")
else
  open (unit=201,file="Unassigned_relatives_full.txt",status="unknown")
endif
write (201, HeaderFMT_pairs) "ID1", "ID2", "Sex1", "Sex2", &
  "AgeDif", "Top_R_U", "LLR_R_U", "Top_R1_R2", "LLR_R1_R2", "OH"
do x=1, nAmb
  i = AmbigID(x,1)
  j = AmbigID(x, 2) 
  write (201, DataFMT_pairs) Id(i), Id(j), &
    Sex(i), Sex(j), AgeDiff(i,j), RelName(AmbigRel(x, 1)), &
    AmbigLR(x,1), RelName(AmbigRel(x,2)), AmbigLR(x,2),CalcOH(i,j)
enddo  
close(201)

if (ntrio>0) then
  open (unit=202, file="Unassigned_triads.txt",status="unknown")
  write (202, HeaderFMT_triads) "id", "parent1", "parent2", &
   "OH_P1", "OH_P2", "OH_PP", "LLRparent1", "LLRparent2", "LLRpair"
  do i=1, ntrio
     write (202, DataFMT_triads) Id(trioIDs(i,:)), trioOH(i,:), trioLR(i,:)
  enddo
  close(202)
endif

end subroutine findambig

! ######################################################################

subroutine getpairll(PairsFileName, nP)  
use Global
use sqa_fileIO, ONLY: FileNumCol
implicit none

character(len=nchar_filename), intent(IN) :: PairsFileName
integer, intent(IN) :: nP
character(len=8) :: colNamesIN(20), colNamesOUT(9)
character(len=nchar_ID) :: DumC(20), PairNames(nP,2)
character(len=2) :: RelName(9), PairFocalC(nP)
character(len=nchar_filename) :: OutFileName
character(len=200) :: HeaderFMT, DataFMT
integer :: IOerr, nCols, theseCols(9), x, y, pairIDs(nP,2), psex(nP,2), &
  pairAgediff(nP), pairFocal(nP), pairk(nP), pdrop(nP, 2), &
  ij(2), kij(2), Sex_ij(2), a, m, curPar(2,2), top(nP), BYtmp(2), t, viginti_pairs(20)
double precision :: LLpair(nP, 7), LLa(7), dl(nP), IndBYtmp(1:nYears,2,5)

ColNamesOUT = (/ "ID1     ", "ID2     ", "Sex1    ", "Sex2    ", "AgeDif  ", &
  "focal   ", "patmat  ", "dropPar1", "dropPar2" /)

nCols = FileNumCol(trim(PairsFileName))
if (nCols > 20) then
  print *, "WARNING: columns past column 20 are ignored"
  nCols = 20
endif

theseCols = 0
open(unit=103, file=trim(PairsFileName), status="old")
read(103,*) colNamesIN(1:nCols)
do x = 1, 9
  do y = 1, nCols
    if (ColNamesOUT(x) == ColNamesIN(y))  theseCols(x) = y
  enddo
enddo

! defaults, if a column not in file
pSex = 3
PairAgeDiff = 999
pairfocalC = "U "
pairfocal = 7
pairk = 3
pDrop = 0

do x=1,nP
  read(103,*,IOSTAT=IOerr)  DumC(1:nCols)
  if (IOerr > 0) then
    print *, "Wrong input on line ", x
    call ErStop("", .FALSE.)
  else if (IOerr < 0) then  ! EOF
    exit
  else

    if (ANY(theseCols(1:2)==0)) then
      call Erstop(trim(PairsFileName)//" must have at least columns 'ID1' and 'ID2'", .FALSE.)
    else  
      PairNames(x,:) = DumC( theseCols(1:2) )
    endif
    if (theseCols(3)/=0)  read(DumC(theseCols(3))(1:2), '(i2)')  pSex(x,1)
    if (theseCols(4)/=0)  read(DumC(theseCols(4))(1:2), '(i2)')  pSex(x,2)
    if (theseCols(5)/=0)  read(DumC(theseCols(5))(1:4), '(i4)')  pairAgeDiff(x)
    if (theseCols(6)/=0)  read(DumC(theseCols(6))(1:2), '(a2)')  pairfocalC(x)
    if (theseCols(7)/=0)  read(DumC(theseCols(7))(1:2), '(i2)')  pairk(x)
    if (theseCols(8)/=0)  read(DumC(theseCols(8))(1:2), '(i2)')  pDrop(x,1)
    if (theseCols(9)/=0)  read(DumC(theseCols(9))(1:2), '(i2)')  pDrop(x,2)
  end if
enddo
close(103)

! ID name to num
do x=1,nP
  do a=1,2
    do y=1, nInd
      if (PairNames(x,a) == ID(y)) then
        pairIDs(x,a) = y
      endif
    enddo
  enddo
enddo

! focal char to num
RelName = (/ "PO", "FS", "HS", "GP", "FA", "HA", "U ", "XX", "X2" /)
do x=1,nP
  do y=1,7
    if (pairFocalC(x) == RelName(y)) then
      PairFocal(x) = y
    endif
  enddo
enddo

LLpair = Missing
top = 7
dL = -999D0
viginti_pairs =  mk_seq(nP, 20)
t = 1
do x = 1, nP
  if (quiet<1 .and. any(viginti_pairs == x))   call print_progress(x,t)
  ij = pairIDs(x,:)
  
  ! temp. drop parents
  curPar = 0
  do a=1,2
    curPar(a,:) = getPar(ij(a), psex(x,a))
    do m=1,2
      if (pdrop(x,a)==m .or. pdrop(x,a)==3) then 
        call setParTmp(ij(a), psex(x,a), 0, m)
      endif
    enddo
  enddo
  
  ! backup age difference from lifehistory data & sex from pedigree
  Sex_ij = 0
  BYtmp = -9
  IndBYtmp = LOG10(1.0D0/nYears)
  do a=1,2
    if (ij(a) < 0)  cycle
    Sex_ij(a) = sex(ij(a))
    IndBYtmp(:,a,:) = IndBY(:,ij(a),:)
  enddo
  if (all(ij>0))  BYtmp(1:2) = BY(ij)
  
  ! set sex & age difference specified in pairs DF
  do a=1,2
    if (ij(a) < 0)  cycle
    if (Sex(ij(a))==3) then   ! can not change sex if assigned as parent in pedigree
      Sex(ij(a)) = psex(x,a)  
    else  ! for output file only
      psex(x,a) = Sex(ij(a))
    endif
    IndBY(:,ij(a),:) = LOG10(1.0D0/nYears)    
  enddo

  ! set age difference 
  if (theseCols(5)/=0 .and. all(ij>0)) then
    ! TODO
    ! AgeDiff(ij(1), ij(2)) = pairAgeDiff(x) 
    ! if (pairAgeDiff(x) /= 999) then
      ! AgeDiff(ij(2), ij(1)) = -pairAgeDiff(x)
    ! else
      ! AgeDiff(ij(2), ij(1)) = 999
    ! endif
  else  ! for otput file only
    pairAgeDiff(x) = AgeDiff(ij(1), ij(2))
  endif
   
  ! calc likelihoods
  if (all(ij>0) .and. theseCols(7)/=0) then  ! column with 'pairk'
    kij = pairk(x)
  else
    kij = psex(x, :)
  endif
  if (pairAgeDiff(x) < 0 .and. pairfocal(x)==7) then  ! swap, else LL_PO & LL_GP not calculated
    call CheckRel(ij(2), kij(2), ij(1), kij(1), pairfocal(x), LLpair(x,:), LLa)
  else
    call CheckRel(ij(1), kij(1), ij(2), kij(2), pairfocal(x), LLpair(x,:), LLa)
  endif
  call BestRel2(LLpair(x,:), top(x), dL(x))
  
  ! restore sex, agediff, IndBY, Parents
  do a=1,2
    if (ij(a) > 0) then
      sex(ij(a)) = sex_ij(a) 
      IndBY(:,ij(a),:) = IndBYtmp(:,a,:)
    endif
    do m=1,2
      call setParTmp(ij(a), psex(x,a), curPar(a, m), m)
    enddo
  enddo
  if (all(ij>0)) then 
    BY(ij) = BYtmp(1:2)
  endif

enddo

! write output to file
! output file name = input file name + '_LL'
a = index(PairsFileName, ".txt")   ! find location of ".txt" in filename string
OutFileName = PairsFileName(1:(a-1))//"_LL.txt"

write(HeaderFMT,'( "(2(a", I0, ", 4X), 2a5,  3a7,1X, 2a9, 7a10,a7, a10)" )')  ID_len
write(DataFMT,  '( "(2(a", I0, ", 4X), 2i5, i7,a7,i7,1X, 2i9, 7f10.2,a7, f10.2)" )')  ID_len

! 40 = nchar_ID
open (unit=203,file=trim(Outfilename), status="unknown") 
  write (203, HeaderFMT) ColNamesOUT, RelName(1:7), "TopRel", "LLR"
do x=1,nP
  write (203, DataFMT) PairNames(x,:), pSex(x,:), PairAgeDiff(x), PairfocalC(x), &
    pairk(x), pDrop(x,:), LLpair(x,:), RelName(top(x)), dL(x)
enddo   
close (203)

end subroutine getpairll

! ######################################################################

subroutine parents(TotLL, withALog)
use qsort_c_module 
use Global
implicit none

double precision, intent(INOUT) :: TotLL(42)
logical, intent(IN) :: withALog
integer :: i, j, k, Round, isP(2), BYRank(nInd), BYRank_Selfed(nInd)
double precision :: SortBY(nInd), RoundTime(2)
character(len=2) :: RoundC
character(len=200) :: AssignmentLogFile
 
call rchkusr()     
 
AgePhase = 0
  
!============================

call UpdateAllProbs()

if (quiet<1)  call timestamp(.TRUE.)
if (quiet<1)  write(*, '("Initial Total LL   : ", f12.1)')  SUM(Lind)
 
!============================
! get birthyear ranking (increasing)
call getRank_i(BYrank)  ! do earlier birth years before later birth years

if (hermaphrodites/=0) then  ! do selfed individuals first to reduce false positives
  BYRank_Selfed = (/ (i, i=1, nInd, 1) /)
  do i=1, nInd
    call IsSelfed(i, .FALSE., SortBY(i))
    SortBY(i) = -SortBY(i)
  enddo
  call QsortC(SortBy, BYRank_Selfed)
endif

TotLL = 0D0
TotLL(1) = SUM(Lind)
do Round=1,41
  call rchkusr()
  call cpu_time(RoundTime(1))
  write(RoundC, '(i2.2)') Round
  if (withALog) then
    AssignmentLogFile = 'AssignmentLog'//RoundC//'.txt'
  else
    AssignmentLogFile = 'NoLog'
  endif
  if (hermaphrodites/=0) then
    call Parentage(BYRank_Selfed, AssignmentLogFile)
  else
    call Parentage(BYrank, AssignmentLogFile)   
  endif
  call UpdateAllProbs()
  if (any (BY < 0))  call getRank_i(BYrank)
  
  call cpu_time(RoundTime(2))
  if (quiet==-1)  call timestamp(.TRUE.)
  if (quiet==-1)  write(*, '("Round: ", i2, ", Total LL: ", f12.1, ", # parents:", 2i6, ", time (s): ", f7.1)') &
    Round, SUM(Lind), count(Parent/=0, DIM=1), RoundTime(2) - RoundTime(1)
  do i=1,nInd
    if (Sex(i)==3) then
      isP = 0
      do k=1,2
        do j=1,nInd
          if (Parent(j,k) == i) then
            isP(k) = isP(k) + 1
          endif
        enddo
      enddo
      if (isP(1)>0 .and. isP(2)>0) then
 !       call rwarn("individual assigned as both dam & sire")
        print *, "WARNING: individual assigned as both dam & sire: ", i
      else
        do k=1,2
          if (isP(k)>1) then
            Sex(i) = k
          endif
        enddo
      endif
    endif     
  enddo

  TotLL(Round + 1) = SUM(Lind)
  if (TotLL(Round + 1) - TotLL(Round) < ABS(TF))  exit ! convergence
  if (Round==41) then
    call Erstop("parentage not converging - need better SNP data", .FALSE.)
  endif
  
  call writeped("Parents_round"//RoundC//".txt", .FALSE.)
enddo

if (quiet<1)  call timestamp(.TRUE.)
if(quiet<1)   write(*, '("Final Total LL     : ", f12.1, ", # parents:", 2i6)') &
   SUM(Lind), count(Parent/=0, DIM=1)

end subroutine parents

! ####################################################################

! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

! ####################################################################

subroutine sibships(AgeEffect, ResumePed)
use Global
implicit none

integer, intent(IN) :: AgeEffect, ResumePed   ! IN
double precision :: CurLL(8), RoundTime(2), CurTime(0:8)
integer :: MaxRounds, Round, RX, k, PairID(XP*nInd,2), PairType(XP*nInd)
logical :: PairsFileFound 
character(len=2) :: RoundC, StartPedC, RoundXC


MaxRounds = 42  ! fail safe, max. no. iterations
RX = 1  ! no. of initial rounds, pairs-cluster-merge only 
!RoundEA = .FALSE.  ! round with extra ageprior for GGpairs (lags 1 round)
! 0: no ageprior; 1: with + w/o; 2: extra age in last rounds
if (AgeEffect==0) then  !  .or. nYears==1
  AgePhase = 0  ! do not use ageprior
else
  AgePhase = 1  ! use age prior
endif

call UpdateAllProbs()    

if (quiet<1)  call timestamp(.TRUE.)
if (quiet<1)  write(*, '("Initial Total LogLik: ", f12.1, "  # parents:", 2i6)') &
      SUM(Lind), count(Parent/=0, DIM=1)
      
PairsFileFound = .FALSE.

open(unit=202, file="LogLik.txt", status="unknown")
  write(202, '(a6, 8a15)') "round", "start", "cluster", "GGpairs", "merge", &
     "sibPar", "sibGP", "morePar" 
close(202)
  
do Round=1, MaxRounds
  curLL = missing
  write(RoundC, '(i2.2)') Round
  call rchkusr()
    
  !.........................................
  if (ResumePed>=0) then   ! (when debugging) resume at intermediate ped
    if (Round == 1) then   ! StartPed == 0 .and. 
    
      if (ResumePed > 0) then
        write(StartPedC, '(i2.2)') ResumePed
        call ReadPedFile("Pedigree_round"//StartPedC//".txt")
!        do i=1, nInd
!          call CalcFSLik(i)
!        enddo
        call UpdateAllProbs()
        if (quiet<1) call timestamp(.TRUE.)
        if (quiet<1) write(*, '("Total LogLik: ", f12.1, "  # parents:", 2i6)') &
          SUM(Lind), count(Parent/=0, DIM=1)
        do k=1,2
          IsNewSibship(1:nC(k), k) = .FALSE.
        enddo
        if (AgeEffect==2 .and. ResumePed>4)  AgePhase = 2   !!
!        if (AgeEffect==2 .and. ResumePed>5)  RoundEA = .TRUE.  !!
      endif
      
      write(RoundXC, '(i2.2)') ResumePed +1
      call ReadPairs(RoundXC, PairID, PairType, PairsFileFound)      
      
      if (ResumePed>0)  cycle
      
    else if (Round <= ResumePed) then
      cycle
    
    endif
  endif   ! debug
  !.........................................
  
  if(quiet==-1) then 
    write(*,*) ""
    write(*,'(10x, 50a)')  REPEAT('-',29)
    call timestamp()
    write(*, '(" -------   Round ", i3, "   -------")') Round
    write(*,'(10x, 50a)') REPEAT('-',29)
    write(*,*) ""
  endif
  
  open (unit=42,file="log.txt",status="unknown", position="append")
  write (42, *) ""
  write (42, *) "=============================================================="
  write (42, *) "===   ROUND ", Round
  write (42, *) "=============================================================="
  write (42, *) ""
  close(42)
  
  call cpu_time(RoundTime(1))
  CurLL(1) = SUM(Lind)
  call cpu_time(CurTime(0))
    
!  if (nPairs==0) then
  if (ResumePed<0 .or. Round > (ResumePed +1) .or. .not. PairsFileFound) then
    if(quiet==-1)  call Rprint("Find pairs ...", (/0/), (/0.0D0/), "NON")
    call FindPairs(RoundC, PairID, PairType)
  endif
  call cpu_time(CurTime(1))
  if(quiet==-1)  call Rprint("n pairs:", (/npairs/), (/0.0D0/), "INT")  
  if(quiet==-1)  call timestamp(.TRUE.)
  if(quiet==-1)  write(*, '("Total LL : ", f12.1, ", time: ", f5.1, "  ", 2i6)') SUM(Lind), CurTime(1) - CurTime(0), nC
  call rchkusr()
  
  if(quiet==-1)  call Rprint("Clustering ...", (/0/), (/0.0D0/), "NON")
  call Clustering(PairID, PairType) 
  call UpdateAllProbs()  
  CurLL(2) = SUM(Lind)
  call cpu_time(CurTime(2))
  if(quiet==-1)  call timestamp(.TRUE.)
  if(quiet==-1) write(*, '("Total LL : ", f12.1, ", time: ", f5.1, "  ", 2i4, "  ", 2i6)') &
    curLL(2), CurTime(2) - CurTime(1), nC, count(Parent/=0, DIM=1)      
  
  if (Round > RX+1) then
    if(quiet==-1)  call Rprint("Grandparent-grandoffspring pairs ...", (/0/), (/0.0D0/), "NON")  
    call GGpairs()
    call UpdateAllProbs()
  endif
  
  CurLL(3) = SUM(Lind)
  call cpu_time(CurTime(3))
  if (Round > RX+1) then
    if(quiet==-1)  call timestamp(.TRUE.)
    if(quiet==-1) write(*, '("Total LL : ", f12.1, ", time: ", f5.1, "  ", 2i4)') curLL(3), CurTime(3) - CurTime(2), nC
  endif
  call rchkusr()
  
  if(quiet==-1)  call Rprint("Merge clusters ...", (/0/), (/0.0D0/), "NON")
  call Merging()
  call UpdateAllProbs() 
  CurLL(4) = SUM(Lind)
  call cpu_time(CurTime(4))
  if(quiet==-1)  call timestamp(.TRUE.)
  if(quiet==-1) write(*, '("Total LL : ", f12.1, ", time: ", f5.1, "  ", 2i4)') curLL(4), CurTime(4) - CurTime(3), nC  
  call rchkusr()
  
  if(quiet==-1)  call Rprint("Sibship parent replacement...", (/0/), (/0.0D0/), "NON")
  call SibParent()   ! replace dummy parents by indivs
  call UpdateAllProbs()
  CurLL(5) = SUM(Lind)
  call cpu_time(CurTime(5))
  if(quiet==-1)  call timestamp(.TRUE.)
  if(quiet==-1) write(*, '("Total LL : ", f12.1, ", time: ", f5.1, "  ", 2i4)') curLL(5), CurTime(5) - CurTime(4), nC

  if (Round > RX .or. Round==MaxRounds) then  
    if(quiet==-1)  call Rprint("Grandparents of half-sibships ...", (/0/), (/0.0D0/), "NON")
    call SibGrandparents()
    call UpdateAllProbs() 
    CurLL(6) = SUM(Lind)    
    call cpu_time(CurTime(6))
    if(quiet==-1)  call timestamp(.TRUE.)
    if(quiet==-1) write(*, '("Total LL : ", f12.1, ", time: ", f6.1, "  ", 2i4)') curLL(6), CurTime(6) - CurTime(5), nC
   
    if (Round > RX+1) then  
      if(quiet==-1)  call Rprint("Grandparents of full sibships ...", (/0/), (/0.0D0/), "NON")
      call FsibsGPs()
      call UpdateAllProbs() 
    endif
    ! no need to check for GPs of non-changed sibships next round:
    do k=1,2
      IsNewSibship(1:nC(k), k) = .FALSE.
    enddo      
  endif
  CurLL(7) = SUM(Lind)
  call cpu_time(CurTime(7))
  call rchkusr() 
  if (Round > RX+1) then 
    if(quiet==-1)  call timestamp(.TRUE.)
    if(quiet==-1) write(*, '("Total LL : ", f12.1, ", time: ", f6.1, "  ", 2i4, "  ", 2i6)') &
      curLL(7), CurTime(7) - CurTime(6), nC, count(Parent/=0, DIM=1)
  endif

  if(quiet==-1)  call Rprint("Parents & grow clusters...", (/0/), (/0.0D0/), "NON")
  call MoreParent()  !  assign additional parents to singletons 
  call UpdateAllProbs() 
  CurLL(8) = SUM(Lind)
  call cpu_time(CurTime(8))
  if(quiet==-1)  call timestamp(.TRUE.)  
  if(quiet==-1) write(*, '("Total LL : ", f12.1, ", time: ", f6.1, "  ", 2i4, "  ", 2i6)') &
    curLL(8), CurTime(8) - CurTime(7), nC, count(Parent/=0, DIM=1) 
  call rchkusr()
  
  if(quiet<1) then
    call cpu_time(RoundTime(2))
    call timestamp(.TRUE.)
    write(*, '("Round ", i2, ":  Total LogLik: ", f12.1, "  time (sec): ", f6.1, '// &
      '"  # parents:", 2i6)') Round, SUM(Lind), RoundTime(2) - RoundTime(1), count(Parent/=0, DIM=1)
    if(quiet==-1)  call Rprint("No. dummies: ", nC, (/0.0D0/), "INT") 
    if(quiet==-1)  write(*,*) ""
    if(quiet==-1)  write(*,*) ""
  endif

  open(unit=202, file="LogLik.txt", status="unknown", position="append")
  write(202, '(i6, 8f15.1)') Round, curLL
  close(202)    
    
  if (Round == MaxRounds) then
    if (quiet==-1)  write(*,'(10x, 50a)')  REPEAT('-',29)
    exit
  else if (curLL(8) - curLL(1) < 2D0*ABS(TF) .and. &
   ((Round > 2 .and. ResumePed <= 0) .or. (Round > ResumePed+1 .and. ResumePed > 0))) then
    if (AgeEffect==2 .and. AgePhase==1) then
      AgePhase = 2
    else    
      if (quiet==-1)  write(*,'(10x, 50a)')  REPEAT('-',29)
      exit
    endif
  endif
  
  call writeped("Pedigree_round"//RoundC//".txt", .FALSE.)
enddo  

end subroutine Sibships

! ####################################################################

subroutine readpairs(RoundC, PairID, PairType, FileFound)
use Global
implicit none

character(len=2), intent(IN) :: RoundC
integer, intent(OUT) :: PairID(XP*nInd,2), PairType(XP*nInd)
logical, intent(OUT) :: FileFound
double precision :: PairLLR   ! not used by Clustering, only for ranking
integer :: i, IOerr
character(len=nchar_ID) :: DumC(2) 

nPairs = 0  ! nPairs is global
inquire(file="Pairs_"//RoundC//".txt", exist = FileFound)
if (.not. FileFound)  return

if(quiet==-1)  call printt("Reading Pairs... ")
open(unit=103, file="Pairs_"//RoundC//".txt", status="old")
read(103,*)   ! header
do i=1,nInd**2
  read(103,*,IOSTAT=IOerr)  DumC(1:2), PairID(i,1:2), PairType(i), PairLLR
  if (IOerr > 0) then
    print *, "Wrong input in Pairs_"//RoundC//".txt on line ", i
    call Erstop("", .FALSE.)
  else if (IOerr < 0) then
    exit
  else
    nPairs = nPairs +1  
  end if
enddo
close(103)
if(quiet==-1)  call timestamp()
if(quiet==-1)  print *, "Read ", nPairs ," pairs from Pairs_"//RoundC//".txt"

end subroutine readpairs

! ####################################################################

subroutine writeped(filename, CalcLLR)
use Global
use OHfun
implicit none

character(len=*), intent(IN) :: filename
logical, intent(IN) :: CalcLLR
integer :: k, s, i, n, OppHomDF(nInd,3)
double precision :: LLR_Parent(nInd,3), LLR_GP(3, nInd/2, 2)   
character(len=nchar_ID) :: DumName(nInd/2,2), ParentName(nInd,2), GpName(2,nInd/2,2)
character(len=200) :: HeaderFMT, DataFMT
logical :: IsClone2(nInd/2, 2)

LLR_parent = missing
LLR_GP = missing
if (CalcLLR) then
  call CalcParentLLR(LLR_parent, LLR_GP)   
endif

call MakeDumNames(DumName, IsClone2)

ParentName = "NA"
GpName = "NA"

do k=1,2
  do i=1,nInd
    if (Parent(i,k) == 0) cycle
    if (Parent(i,k)>0) then
      ParentName(i,k) = Id(Parent(i,k))
    else if (Parent(i,k)<0) then
      ParentName(i,k) = DumName(-Parent(i,k),k)
    endif
  enddo
  if (nC(k)==0)  cycle
  do s=1,nC(k)
    do n=1,2
      if (GpID(n,s,k)>0) then
        GpName(n,s,k) = Id(GpID(n,s,k))
      else if (GpID(n,s,k)<0) then
        GpName(n,s,k) = DumName(-GpID(n,s,k),n) 
      endif
    enddo
  enddo
enddo

OppHomDF = -9
do i=1,nInd
  if (skip(i)) cycle
  do k=1,2
    if (Parent(i,k) <= 0)  cycle
    OppHomDF(i,k) = calcOH(i, Parent(i,k)) 
  enddo
  if (Parent(i,1)>0 .and. Parent(i,2)>0) then
    OppHomDF(i,3) = CalcTrioErr(i, Parent(i,:))
  endif
enddo

write(HeaderFMT, '( "(3(a", I0, ", 4X), 3a10, 3a8, 4a6)" )')  ID_len
write(DataFMT, '( "(3(a", I0, ", 4X), 3f10.2, 3i8, 4i6)" )')  ID_len

open (unit=201,file=trim(filename), status="unknown") 
  write (201, HeaderFMT) "id", "dam", "sire", "LLRdam", "LLRsire", "LLRpair", &
   "OHdam", "OHsire", "MEpair", "RowO", "RowD", "RowS", "Sex"
  do i=1,nInd
    if (skip(i)) cycle
    write (201,DataFMT) Id(i), ParentName(i,1:2), &
      LLR_parent(i,:), OppHomDF(i,1:3), i, Parent(i,1:2), Sex(i)
  enddo   
  do k=1,2
    if (nC(k)==0)  cycle
    do s=1,nC(k)
      ! TODO: skip clusters if all sibID's are skip(i)=.TRUE. (& not otherwise connected)
      if (IsClone2(s,k))  cycle
      write (201,DataFMT) DumName(s,k), GpName(:, s,k), &
        LLR_GP(:, s, k), -9, -9, -9, -s, GpID(:,s,k), k
    enddo
  enddo
close (201)

!print *, "Assigned ", count(Parent(:,1)/=0), " dams and ", count(Parent(:,2)/=0), " sires."  

end subroutine writeped

! #####################################################################

subroutine MakeDumNames(DumName, IsClone2)
use Global
implicit none

character(len=nchar_ID), intent(OUT) :: DumName(nInd/2,2)
logical, intent(OUT) :: IsClone2(nInd/2, 2)
integer :: k, s
character(len=4) :: DumTmp

DumName = "NA"
IsClone2 = .FALSE.

do k=1,2
  do s=1, nC(k)
    write(DumTmp, '(i4.4)') s
    if (DummyNamesIO(s,k) /= 'NA') then
      DumName(s,k) = DummyNamesIO(s,k)
    else if (hermaphrodites == 0 .or. DumClone(s,k)==0) then
      DumName(s,k) = trim(DumPrefix(k))//DumTmp
    else
      if (k==1) then
        DumName(s,k) = trim(DumPrefix(3))//DumTmp
      else
        write(DumTmp, '(i4.4)') DumClone(s,k)
        DumName(s,k) = trim(DumPrefix(3))//DumTmp
      endif
      if (ns(s,k) < ns(DumClone(s,k), 3-k) .or. &
        (ns(s,k) == ns(DumClone(s,k), 3-k) .and. k==2)) then
        IsClone2(s,k) = .TRUE.   ! print GPs & GP-LLRs for clone w most offspring only
      endif
    endif
  enddo
enddo

end subroutine MakeDumNames

! #####################################################################

! @@@@   SUBROUTINES   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

! #####################################################################


subroutine CalcPO2(A,B,C, LL)   ! LL of A, given B and/or C as parent  
use Global
implicit none

integer, intent(IN) :: A, B, C
double precision, intent(OUT) :: LL
integer :: l, x, y
double precision :: PrL(nSnp), tmp(3,3), PrB(3), PrC(3)

LL = missing
PrL = 0D0
do l=1,nSnp
  if (Genos(l,A)==-1) cycle  
  PrB = LindX(:,l,B) / SUM(LindX(:,l,B))
  PrC = LindX(:,l,C) / SUM(LindX(:,l,C))
  do x=1,3
    do y=1,3
      if (B>0 .and. C>0) then
        tmp(x,y) = OKA2P(Genos(l,A),x,y) * PrB(x) * PrC(y)
      else if (B>0 .and. C==0) then
        tmp(x,y) = OKA2P(Genos(l,A),x,y) * PrB(x) * AHWE(y,l)
      else if (B==0 .and. C>0) then
        tmp(x,y) = OKA2P(Genos(l,A),x,y) * AHWE(x,l) * PrC(y)
      else if (B==0 .and. C==0) then
        tmp(x,y) = OKA2P(Genos(l,A),x,y) * AHWE(x,l) * AHWE(y,l)
      else
        call Erstop("invalid call to CalcPO2", .TRUE.)
      endif
    enddo
  enddo
  PrL(l) = LOG10(sum(tmp))
enddo
LL = SUM(PrL)

end subroutine CalcPO2

! ######################################################################

subroutine CalcP2(A, kA, B, C, kB, LLR)   ! LR of A; B & C both parent vs both U  
use Global
use OHfun
implicit none

integer, intent(IN) :: A, kA, B, C, kB
double precision, intent(OUT) :: LLR
integer :: l, x, y, z, P1, P2
double precision :: PrL(nSnp,2), PrXYZ(3,3,3,2), PrP(3,2), PrA(3)
logical :: Selfed

if (A==0 .or. (B==0 .and. C==0)) then
  LLR = 0D0
  return
endif

if (kB==1) then
  P1 = B
  P2 = C
else
  P2 = B
  P1 = C
endif

LLR = missing
if (A>0 .and. P1>0) then
  if (calcOH(A,P1) > maxOppHom) then
    LLR = impossible
  else
    if (QLR_PO(A,P1) < 5.0*TF)  LLR = impossible
  endif
endif
if (LLR == impossible)  return
if (A>0 .and. P2>0) then
  if (calcOH(A,P2) > maxOppHom) then
    LLR = impossible
  else
    if (QLR_PO(A,P2) < 5.0*TF)  LLR = impossible
  endif
endif
if (LLR == impossible)  return
  
if (A>0 .and. P1>0 .and. P2>0)  then 
  if (CalcTrioErr(A, (/P1,P2/)) > MaxMendelE) then
    LLR = impossible
    return
  endif
endif

! TODO: Qadd if A or P < 0 ?

Selfed = .FALSE.
if (hermaphrodites/=0 .and. A>0) then
  if (SelfedIndiv(A))  Selfed = .TRUE.
endif

PrL = 0D0
do l=1,nSnp
  call ParProb(l, A, kA, -4, 0, PrA)
  call ParProb(l, P1, 1, 0,0, PrP(:,1))
  call ParProb(l, P2, 2, 0,0, PrP(:,2))
  PrXYZ = 0D0
  do x=1,3
    do y=1,3
      if (Selfed) then
        PrXYZ(x,y,y,1) = PrA(x) * AKA2P(x,y,y) * PrP(y,1)
        PrXYZ(x,y,y,2) = PrA(x) * AKA2P(x,y,y) * AHWE(y,l)
      else
        do z=1,3
          PrXYZ(x,y,z,1) = PrA(x) * AKA2P(x,y,z) * PrP(y,1) * PrP(z,2)
          PrXYZ(x,y,z,2) = PrA(x) * AKA2P(x,y,z) * AHWE(y,l) * AHWE(z,l)
        enddo
      endif
    enddo
  enddo
  PrL(l,1) = LOG10(sum(PrXYZ(:,:,:,1)))
  PrL(l,2) = LOG10(sum(PrXYZ(:,:,:,2)))
enddo

if (SUM(PrL(:,1)) < -HUGE(0D0)) then
  LLR = impossible
else
  LLR = SUM(PrL(:,1)) - SUM(PrL(:,2)) 
endif

end subroutine CalcP2

! ######################################################################

subroutine ChkValidPar(A, kA, P, kP, OK)  ! age, ancestor, & OH check
use Global
implicit none

integer, intent(IN) :: A, kA, P, kP
logical, intent(OUT) :: OK
logical :: AncOK, MoreChk
double precision :: ALR, LRQ
integer :: ParA(2), i

if (A==0 .or. P==0) then
  OK = .TRUE.
  return
endif

ParA = getpar(A, kA)
if (kP < 3) then
  if (ParA(kP) == P) then
    OK = .TRUE.
    return
  endif
endif

! if (Complx==0) then
  ! OK = .TRUE.
  ! if (ParA(3-kP) > 0) then
    ! if (Mate(ParA(3-kP))/=0 .and. Mate(ParA(3-kP)) /= P)  OK = .FALSE.
  ! else if (ParA(3-kP) < 0) then
    ! if (DumMate(-ParA(3-kP),3-kP)/=0 .and. DumMate(-ParA(3-kP),3-kP)/=P)  OK = .FALSE.
  ! endif
  ! if (.not. OK)  return
! endif   ! TODO: THIS PART CAUSES R TO CRASH

OK = .FALSE.
AncOK = .FALSE.
ALR = missing
LRQ = missing

call ChkAncest(P, kP, A, kA, AncOK)  ! check that A is not an ancestor of P
if (.not. AncOK)  return

call CalcAgeLR(A,kA, P, kP, 0,1, .TRUE., ALR) 
if (ALR == impossible)  return

if (kP < 3) then
  call CalcP2(A, kA, P, ParA(3-kP), kP, LRQ)  
else
  call CalcP2(A, kA, P, 0, kP, LRQ)  
endif
if (LRQ == impossible)  return

if (A>0 .and. P<0 .and. ALR < -TA) then   ! check age diff with future sibs
  do i=1, ns(-P,kP)
    call CalcAgeLR(A, kA, SibID(i,-P,kP), 3, kP, 3, .TRUE., ALR)
    if (ALR == impossible)  return
  enddo
endif

OK = .TRUE.

MoreChk = .FALSE.
if (any(OKA2P < TINY(0D0))) then  ! some combi of obs. offspring genotype + act parents is impossible
  if (kP < 3) then
    if (ParA(kP)==0)  MoreChk = .TRUE.
  endif
endif

if (MoreChk) then
  if (A > 0  .and. P < 0 .and. ParA(3-kP)<0) then
    Parent(A, kP) = P
    nS(-P,kP) = nS(-P,kP) +1
    SibID(nS(-P,kP), -P,kP) = A  
    call CalcCLL(-P,kP)
    if (CLL(-P,kP) < -HUGE(0D0))  OK = .FALSE.    ! inbreeding / FindEE loop
    SibID(nS(-P,kP), -P,kP) = 0  
    nS(-P,kP) = nS(-P,kP) -1
    Parent(A,kP) = 0
    call CalcCLL(-P,kP)
  
  else if (A < 0) then
      GpID(kP, -A, kA) = P
      call CalcCLL(-A,kA)
      if (CLL(-A,kA) < -HUGE(0D0))  OK = .FALSE.    ! e.g. inbreeding
      GpID(kP, -A, kA) = 0
      call CalcCLL(-A,kA)
  endif
endif

end subroutine ChkValidPar

! ######################################################################

subroutine ChkDoQuick(s,k, DoQuick)
use Global
implicit none

integer, intent(IN) :: s,k
integer, intent(OUT) :: DoQuick
integer :: nOff, Offspr(maxSibSize), sxOff(maxSibSize), i, OpPar, &
  UseEE(ns(s,k)), Sibs(ns(s,k)), MatePar(ns(s,k))

!  1: all parents 3-k >=0, no inbreeding
!  0: some parents 3-k <0, no inbreeding
! -1: inbreeding: Parent(Offspr(i),3-k) == GpID(3-k,s,k)
! -2: all are FS; ns(s,k) = ns(..,3-k)
! -3: s has a dummy clone
!  2: inbreeding: Parent(Bj,3-k) = Offspr(i)
!  3: Parent(Bi,3-k) close relative of Parent(Bj,3-k)

DoQuick = 1
if (ns(s,k)==0)  return

if (any(Parent(SibID(1:ns(s,k),s,k),3-k) < 0)) then   !  .and. ns(s,k)<20
  DoQuick = 0
else
  DoQuick = 1
endif

OpPar = 0         
if (all(Parent(SibID(1:ns(s,k),s,k),3-k) < 0)) then
  call getFSpar(s, k, .TRUE., OpPar)
  if (OpPar < 0) then
    if (ns(-opPar, 3-k) == ns(s,k))  DoQuick = -2
  endif
endif

if (Hermaphrodites/=0) then
  if (DumClone(s,k)/=0)  DoQuick = -3   
endif

nOff = 0
Offspr = 0
sxOff = 3
call getOff(-s,k, .TRUE., nOff, Offspr, sxOff)
if (nOff == 0) then   ! may happen temporarily
  DoQuick = 1
else
  do i=1, nOff
    if (Offspr(i) < 0 .and. sxOff(i)==3-k) then
      if (any(Parent(SibID(1:ns(s,k),s,k),3-k) == Offspr(i))) then
        DoQuick = 2
      endif
    else if (Offspr(i) > 0) then
      if (Parent(Offspr(i),3-k) == GpID(3-k,s,k) .and. GpID(3-k,s,k)/=0) then
        DoQuick = -1
        exit
      endif
    endif
  enddo
endif     

UseEE = 0         
if ((DoQuick==0 .or. DoQuick==1) .and. any(Parent(SibID(1:ns(s,k),s,k),3-k) < 0)) then
  Sibs = SibID(1:ns(s,k), s, k)
  call FindEE(Sibs, ns(s,k), 0, k, UseEE, MatePar)
  if (any(UseEE/=0))  DoQuick = 3
endif

end subroutine ChkDoQuick

! ######################################################################

subroutine CalcPX2(A, kA, P1, P2, LLR)   ! joined LR of A+B+C; B & C as parent vs either or both U  
use Global
implicit none

integer, intent(IN) :: A, kA, P1, P2
double precision, intent(OUT) :: LLR
integer :: x, curPar(2)   ! May or may not be P1, P2
double precision :: LLY(2,2), LLU(4), LLcor(3,2)

curPar = getPar(A, kA)
do x=1,2
  call setParTmp(A, kA, 0, x)
enddo

LLY = missing
call Calc4U((/P1, P2/), 0,0, A,kA, LLU, LLcor)
LLY(1,1) = LLcor(3,1) + LLU(4)  ! no parents. LLU(4) = CLL(-A,kA)
call setParTmp(A,kA,P1,1)
if (Complx/=0) then
  call CalcU(A,kA, P1,1, LLY(2,1))
  LLY(2,1) = LLY(2,1) + LLcor(1,1)   ! only dam
endif
call setParTmp(A,kA,P2,2)        
call CalcU(A,kA, P1,1, LLY(2,2))
LLY(2,2) = LLY(2,2) + LLcor(1,1)   ! dam + sire
if (Complx/=0) then    
  call setParTmp(A,kA,0,1)
  call CalcU(A,kA, P2,2, LLY(1,2))
  LLY(1,2) = LLY(1,2) + LLcor(2,2)   ! only sire
endif

do x = 1,2
  call setParTmp(A, kA, curPar(x), x)
enddo

LLR = LLY(2,2) - MaxLL((/LLY(1,:), LLY(2,1)/))

if (hermaphrodites/=0) then   ! A>0 & A<0 (?)
  if (P1<0 .and. P2<0) then
    if (DumClone(-P1,1) == -P2) then
      LLR = LLY(2,2) - LLY(1,1)
    endif
  else if (P1>0 .and. P2>0) then
    if (P1 == P2) then
      LLR = LLY(2,2) - LLY(1,1)
    endif
  endif
endif

! if (A==2668 .and. P1==2311)  then
  ! open (unit=42,file='log.txt',status="unknown", position="append")
  ! write(42,*) ""
  ! write(42,*) "CalcPX2:"
  ! write(42,*) A, kA, P1, P2
  ! do x=1,2
    ! write(42,'("LLY:   ", 2f8.1)')  LLY(x,:)
  ! enddo
  ! write(42,*) ""
  ! do x=1,3
    ! write(42,'("LLcor: ", 2f8.1)')  LLcor(x,:)
  ! enddo
  ! write(42,*) ""
  ! close(42)
! endif

end subroutine CalcPX2

! ######################################################################

subroutine CheckPedigree(ParOnly, DropSkip)
use Global
implicit none

logical, intent(IN) :: ParOnly, DropSkip  ! T:only parents / F:also dummy parents
integer :: i, k, x, curPar(2), BYrank(nInd), t
logical :: parOK(2), dropS
double precision :: LLRpair!, LL(7,2)

call getRank_i(BYrank)

do x=1, nInd
  if (MODULO(x,100)==0)  call rchkusr()
  if (quiet==-1 .and. any(viginti==x)) call print_progress(x,t)
  
  i = BYRank(x)
  curPar = Parent(i,:)
  
  if (DropSkip .and. any(skip) .and. .not. skip(i)) then  ! --only takes precedent over --pedigreeIN
    do k=1,2
      call setPar(i,Sex(i), 0,k)
      if (curPar(k) < 0) then
        call CheckDropSibship(-curPar(k), k, DropS)
      endif
    enddo
    cycle
  endif

  do k=1,2
    call setParTmp(i,Sex(i), 0,k)  ! also drops i from sibship
  enddo
  
  parOK = .TRUE.
  do k=1,2
    if (curPar(k)>0) then
      if (Sex(curPar(k)) /= k .and. Sex(curPar(k)) < 3) then
        ParOK(k) = .FALSE.    ! male assigned as dam or female as sire.
      endif
    endif
    if (ParOK(k)) then
      call ChkValidPar(i, Sex(i), curPar(k), k, parOK(k))
    endif
  enddo
    
  ! check parent-pair
  if (Parent(i,1)/=0 .and. Parent(i,2)/=0 .and. all(ParOK)) then
    call CalcPX2(i, Sex(i), curPar(1), curPar(2), LLRpair)
    if (LLRpair < TA) ParOK = .FALSE.
  endif
  
  ! do k=1,2
    ! if (ParOK(k) .and. curPar(k)/=0) then
      ! LL = missing
      ! call CheckRel(i,Sex(i),curPar(k),k, 1, LL(:,1), LL(:,2))  
      ! if (LL(1,2) - MaxLL(LL(:,2)) < TF)  ParOK(k) = .FALSE.
    ! endif
  ! enddo 
  
  do k=1,2
    if (ParOK(k)) then
      call setPar(i,Sex(i), curPar(k),k)
    else if (curPar(k) < 0) then
      dropS = .FALSE.
      call CheckDropSibship(-curPar(k), k, DropS)
!      if (hermaphrodites/=0 .and. .not. DropS) then
!        call CheckSelfed(curPar(k),k)
!      endif
    endif
  enddo
  
  if (Parent(i,1)/=curPar(1) .or. Parent(i,2)/=curPar(2)) then
    if (Complx == 0)  call UpdateMate(i, Sex(i), curPar, ParOnly) 
    if (hermaphrodites /= 0) then
      call CheckSelfed(i,Sex(i))
      if (SelfedIndiv(i) .and. all(Parent(i,:)==0)) then
        do k=1,2
          if (Parent(i,k)==0) then
            call NewSibship(i, 0, k)
          endif
        enddo
        do k=1,2
          DumClone( Parent(i,k), k) = Parent(i,3-k)
        enddo
      endif
    endif
  endif
enddo

end subroutine CheckPedigree

! ######################################################################

subroutine Parentage(BYrank, AssignmentLogFile)
use Global
use qsort_c_module
use OHfun
implicit none

integer, intent(IN) :: BYrank(nInd)
character(len=200), intent(IN) :: AssignmentLogFile
integer :: i, j, x, y, k, CandPar(5*mxCP, 2), nCP(2), curPar(2), SexTmp(2), &
  CP_rank(5*mxCP), CP_tmp(5*mxCP), mxxCP, t
double precision :: ALR, LRQ, LLR_CP(5*mxCP)
logical :: AncOK, withLog

mXXCP = 5*mxCP   ! cannot declare as parameter if mxCP is input variable
AncOK = .FALSE.
ALR = missing
LRQ = missing
CP_tmp = 0

if (AssignmentLogFile == 'NoLog') then
  withLog = .FALSE.
else
  withLog = .TRUE.
endif

if (withLog) then
  open(unit=77, file=trim(AssignmentLogFile), status='unknown')
  ! header
  write(77, '(a8, 2X, 3a20, 7a12, 2X, 2a20, 2X, 2a8, 2X, 2a20, 3a12, a14, 2X, 2a20)', advance='no') &
    'G_row', 'id', 'dam_in', 'sire_in', 'n_cand_dam', 'n_cand_sire', &
     'n_pairs_s1', 'n_pairs_s2', 'n_pairs_s3', 'n_pairs_s4', 'mx_LR_pair', 'dam_pair', 'sire_pair', &
      'dup_LR', 'dup_sex', 'dup_id_A', 'dup_id_B', &  
     'n_single_s1', 'n_single_s2', 'n_single_s3', 'mx_LR_single', 'dam_single', 'sire_single'    
  ! FILE STAYS OPEN !
endif

t=1
do x=1, nInd
  if (MOD(x,200)==0) call rchkusr()  
  if (quiet==-1 .and. any(viginti==x)) call print_progress(x,t)
  i = BYRank(x)  
  
  if (skip(i))  cycle
  if (Parent(i,1)>0 .and. Parent(i,2)>0)  cycle
  curPar = Parent(i,:)
  nCP = 0
  CandPar = 0
  SexTmp = 3           
  do k=1,2
    if(Parent(i,k)/=0) then
      nCP(k) = nCP(k) + 1
      CandPar(nCP(k), k) = Parent(i,k)
    endif   
  enddo 

  if (withLog)  write(77, *) ''  ! new line
  if (withLog)  write(77, '(i8, 2X, 3a20)', advance='no')  i, id(i), id(curpar)
  
  do y=1,nInd 
    j = BYRank(y)
    if (i==j) cycle
    if (sex(j)<3) then
      if (nCP(sex(j))==mxxCP)  cycle
    endif
    if (ANY(Parent(j,:)==i)) cycle
    if (ANY(CandPar==j) .and. Sex(j)<3) cycle  ! already included 
    if (CalcOH(i,j) > maxOppHom)  cycle
    if (QLR_PO(i,j) < TF)  cycle
!      if (OppHomM(i,j) > maxOppHom .or. OppHomM(i,j)<0) cycle   ! implied by LLR_O/=missing   
 !   if (LLR_O(i,j) < TF .or. LLR_O(i,j)==missing) cycle   
    if (AgeDiff(i,j) <= 0)  cycle  
!    if (AgeDiff(i,j)==999) then  ! age difference unknown
      call ChkAncest(j,sex(j), i,sex(i), AncOK)  ! check that i is not an ancestor of j
      if (.not. AncOK)  cycle  
!    endif
    call CalcAgeLR(i,sex(i), j,sex(j), 0,1, .TRUE., ALR) 
    if (ALR == impossible)  cycle  

    do k=1,2
      if (Sex(j) <3 .and. Sex(j)/=k) cycle
      if (nCP(k)==mxxCP) cycle
      if (ANY(CandPar(:,k) == j))  cycle
      if (DoMtDif) then
        if (k==1 .and. Sex(j)>2 .and. mtDif(i,j))  cycle    
      endif
      nCP(k) = nCP(k) + 1
      CandPar(nCP(k), k) = j
      if (Complx==0 .and. Mate(j)/=0) then
        if ((.not. any(CandPar == Mate(j))) .and. nCP(3-k)<mxxCP .and. Mate(j)>0) then
          nCP(3-k) = nCP(3-k) +1
          CandPar(nCP(3-k), 3-k) = Mate(j)
        endif
      endif
    enddo
  enddo
  
  if (withLog)  write(77, '(2i12)', advance='no')  ncp
  
  if (ALL(nCP <=1) .and. ALL(candPar(1,:) == Parent(i,:))) cycle
  ! highly unlikely different sire will be assigned in absence of any cand. dam and v.v.
  if (nCP(1)==0 .and. curPar(2)/=0)  cycle  
  if (nCP(2)==0 .and. curPar(1)/=0)  cycle
  
  if (any(nCP > mxCP)) then ! sort by LLR_O
    do k=1,2
      if (nCP(k) <= mxCP)  cycle
      LLR_CP = missing
      do y = 1, nCP(k)
        LLR_CP(y) = QLR_PO(CandPar(y,k), i)
      enddo
      CP_rank = (/ (y, y=1, mxxCP, 1) /)
      LLR_CP = -LLR_CP   ! to sort decreasing
      call QsortC(LLR_CP(1:nCP(k)), CP_rank(1:nCP(k)))
      CP_tmp(1:nCP(k)) = CandPar(1:nCP(k),k)
      CandPar(1:nCP(k),k) = CP_tmp(CP_rank(1:nCP(k)))
      nCP(k) = mxCP
    enddo
  endif  
  
  call SelectParent(i, Sex(i), nCP, CandPar(1:mxCP,:), .TRUE., withLog)  ! does actual assignment

enddo

if (withLog)  close(77)  ! AssignmentLog.txt

end subroutine Parentage

! ######################################################################

subroutine CheckParentPair(A, kA, Par, dLL)
use Global
implicit none

integer, intent(IN) :: A, kA, Par(2)
double precision, intent(OUT) :: dLL(3)  ! 1:dam, 2:sire, 3: both (vs none or either)
integer :: NowPar(2), m
double precision :: LLRP(2), gLL(4,2)

NowPar = getPar(A, kA)
do m=1,2
  call setParTmp(A, kA, 0, m)
enddo

dLL = missing
LLRP = missing
gLL = missing

call CalcP2(A, kA, Par(1), Par(2), 1, LLRP(1))

if (LLRP(1) > 2*TF .and. LLRP(1) /= impossible)  then
  call CalcPX2(A, Sex(A), Par(1), Par(2), LLRP(2))
  
  if (LLRP(2) > TA)  then
    call setParTmp(A, kA, Par(1), 1)
    call CalcPOGPZ(A, kA, Par(2), 2, gLL)
    dLL = gLL(1:3,1)  ! no ageeffect
  endif
endif

do m=1,2
  call setParTmp(A, kA, NowPar(m), m)
enddo

end subroutine CheckParentPair

! ######################################################################

subroutine SelectParent(A, kAIN, nCP, CandPar, ParOnly, withALog)   
   ! assigns parent / grandparent as side effect
use Global
implicit none

integer, intent(IN) :: A, kAIN, nCP(2), CandPar(mxCP, 2)
logical, intent(IN) :: ParOnly, withALog
integer :: m, u, v, best(2), par, AG, kA, curpar(2), nowparUV(2), x, uv(2), &
   i,j, n, fcl, nSingle, dupi(2), dupm
double precision :: LRS, LLRX(nCP(1),nCP(2)), LLRY(nCP(1),nCP(2)), ALR,   & 
  LLRZpair(nCP(1),nCP(2),2), LLRZsingle(2,mxCP,2), gLL(4,2), TAx, dLLrev(2), LLA(7,7,2,2), LRStmp(2,2), DupLR, mxDupLR
logical ::  AgeAUnk, SexUnk(mxCP, 2), MonoPair(nCP(1),nCP(2)), emptySibship(mxCP, 2), &
  PairCand(nCP(1),nCP(2)), AgeUnk(2), maybeRev(2), DoSingle, SingleCand(mxCP,2), &
  DoneSingleCheck(mxCP,2), DoLog
character(len=200) :: logfile

if (ALL(nCP==0))  return

DoLog = .FALSE.
if (A > 0) then
  fcl = 1
!  if (A==91) DoLog = .TRUE.
else !if (A < 0) then
  fcl = 4
! if (any(SibID(:,-A, kAIN) == 91) )  DoLog = .TRUE.
endif
logfile = "log.txt"

!if (DoLog)  print *, A, ' nCP: ', nCP

if (kAIN >2 .or. kAIN==0) then    ! (only?) necessary for checkMaybeRev
  if (Sex(A)<3) then
    kA = Sex(A)
  else
    kA = 1
  endif
else
  kA = kAIN
endif

LRS = missing
if (hermaphrodites/=0) then
  if (A>0)  call IsSelfed(A, .FALSE., LRS)
!  if (A<0)  call IsSelfed(SibID(1,-A,kA), .TRUE., LRS)   ! TODO?
endif

if (DoLog) then
  open (unit=42,file=trim(logfile),status="unknown", position="append")
  write(42,*) ""
  write(42,*)  "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
  write(42,*) "A: ", A, kA
  if (A<0)  write(42, '("sibs: ", 500i6)')  SibID(1:ns(-A,kA),-A,kA)
  do m=1,2
    write(42,*)  m, "CandPar: ", candpar(1:nCP(m), m)
  enddo
  write(42,*) ""
  if (hermaphrodites/=0) then
    write(42,*)  "LR Selfed: ", LRS
    if (A>0)  write(42,*) "SelfedIndiv: ", SelfedIndiv(A)
    if (A<0)  write(42,*) "DumClone: ", DumClone(-A,kA)
  endif
  do m=1,2
    do u=1, nCP(m)
      x=CandPar(u,m)
      if (x<0)  write(42, '(2i6, ": ", 2i6, ", ", 500i6)')  m, x, GpID(:,-x,m), SibID(1:ns(-x,m), -x,m)
    enddo 
    write(42,*) ""
  enddo
  close(42)
endif

curpar = getPar(A,kA)
do m=1,2
  call setParTmp(A, kA, 0, m)   
  call SetEstBY(curPar(m), m)                             
enddo
call SetEstBY(A, kA)   ! conditional on no parents
  
AgeAUnk = .FALSE.                  
if (A < 0) then
  AgeAUnk = .TRUE.
else if (BY(A) < 0) then
  AgeAUnk = .TRUE.
endif
if (AgePhase == 0) then   ! genetics only
  AG = 1
else
  AG = 2
  ! AgePhase = 1 : genetics & [genetics + age] must pass TA
  ! AgePhase = 2 : only [genetics + age] must pass TA
endif

SexUnk = .FALSE.
MonoPair = .FALSE.  
EmptySibship = .FALSE.              
par = 0                  
do m=1,2
  do u = 1, nCP(m)
    if (CandPar(u,m) > 0) then
      if (Sex(CandPar(u,m)) > 2) then
        SexUnk(u,m) = .TRUE.
        if (DoMtDif) then
          if (A>0) then
            if (m==2 .and. mtDif(A, CandPar(u,m)))  SexUnk(u,m) = .FALSE.   ! cannot be maternal relative
          else if (A < 0) then
            if (m==2 .and. mtDif(SibID(1,-A,kA), CandPar(u,m)))  SexUnk(u,m) = .FALSE.   
          endif
        endif
      endif
    else if (CandPar(u,m) < 0) then
      v = -CandPar(u,m)
      if (ns(v,m) == 0 .and. all(GpID(:,v,m)==0)) then
        EmptySibship(u,m) = .TRUE.
      else if (ALL(Parent(SibID(1:nS(v,m), v, m), 3-m) < 0)) then
        call getFSpar(v, m, .TRUE.,par)
        if (par < 0) then
          if (nS(-par, 3-m) == nS(v,m)) then  ! cannot tell if mat or pat
            SexUnk(u,m) = .TRUE.   ! TODO?
            if (DoMtDif .and. A>0) then
              if (m==2 .and. mtDif(A, SibID(1,v,m)))  SexUnk(u,m) = .FALSE. 
            endif                   
            if (ANY(CandPar(:,3-m) == par)) then
              do v=1, nCP(3-m)
                if (CandPar(v,3-m) == par) then
                  if (m==1)  MonoPair(u,v) = .TRUE.  ! monogamous parent pair
                  if (m==2)  MonoPair(v,u) = .TRUE.                                        
                endif
              enddo
            endif
          endif
        endif          
      endif    
    endif
  enddo
enddo

PairCand = .FALSE.   ! check if candidates form an eligible parent-pair
if (ALL(nCP>0)) then   
  do u=1, nCP(1)
    do v=1, nCP(2)
      PairCand(u,v) = .TRUE.
      
      if (Complx==0) then
        if (CandPar(u,1) > 0) then
          if (Mate(CandPar(u,1))/=0 .and. Mate(CandPar(u,1))/= CandPar(v,2))  PairCand(u,v) = .FALSE.
          if (Mate(CandPar(u,1))==CandPar(v,2))  MonoPair(u,v) = .TRUE.   
        else
          if (DumMate(-CandPar(u,1),1)/=0 .and. DumMate(-CandPar(u,1),1)/= CandPar(v,2))  PairCand(u,v) = .FALSE.
          if (DumMate(-CandPar(u,1),1)==CandPar(v,2))  MonoPair(u,v) = .TRUE. 
        endif
        if (CandPar(v,2) > 0) then
          if (Mate(CandPar(v,2))/=0 .and. Mate(CandPar(v,2))/= CandPar(u,1))  PairCand(u,v) = .FALSE.
        else
          if (DumMate(-CandPar(v,2),2)/=0 .and. DumMate(-CandPar(v,2),2)/= CandPar(u,1))  PairCand(u,v) = .FALSE.
        endif
      endif
      
      if (SexUnk(u,1) .and. SexUnk(v,2)) then
        if (hermaphrodites==0 .and. .not. MonoPair(u,v))  PairCand(u,v) = .FALSE.
        if (hermaphrodites==1) then  ! only if selfed
          if (CandPar(u,1) > 0) then
            if (CandPar(u,1)/=CandPar(v,2))  PairCand(u,v) = .FALSE.
          else
            if (DumClone(-CandPar(u,1),1) /= -CandPar(v,2))  PairCand(u,v) = .FALSE.
          endif
        endif
      endif
      
      if (hermaphrodites/=0 .and. A>0) then
        if (LRS > TA .and. CandPar(u,1)/=CandPar(v,2) .and. (CandPar(u,1)>0 .or. CandPar(v,2)>0)) then
          PairCand(u,v) = .FALSE.
        endif
        if (LRS > TA .and. CandPar(u,1)>0) then
          if (Sex(CandPar(u,1)) /= 4) PairCand(u,v) = .FALSE.
        endif
        if (LRS < 5*TF .and. CandPar(u,1)==CandPar(v,2) .and. CandPar(u,1)>0) then
          PairCand(u,v) = .FALSE.
        endif
      endif

      if (hermaphrodites/=0 .and. CandPar(u,1)>0 .and. CandPar(v,2)>0) then    
        if (Hermaphrodites == 2 .and. any(CandPar(1:(u-1),1) ==CandPar(v,2)) .and. &
          any(CandPar(:,2) == CandPar(u,1)))   PairCand(u,v) = .FALSE.
          ! don't care if dam or sire; drop if pair already considered.
        if (Hermaphrodites == 1 .and. any(CandPar(:,1) == CandPar(v,2)) .and. &
          any(CandPar(:,2) == CandPar(u,1)) .and. &   ! can't distinguish dam-sire vs sire-dam
          .not. candPar(u,1) == CandPar(v,2))   PairCand(u,v) = .FALSE.   ! exception: selfing 
      endif
    enddo
  enddo
endif

if (withALog)  write(77, '(2i12)', advance='no') COUNT(PairCand)   ! _s1

LLRX = missing  ! LR(P/U), w/o parents, ignoring inbreeding etc
LLRY = missing  ! LR(P/U), w parents + their relatedness
if (ANY(PairCand)) then   ! find plausible parent-pairs
  do u=1, nCP(1)
    do v=1, nCP(2)
      if (.not. PairCand(u,v))  cycle
      ! calc LLR parent pair / both unrelated --> LLRX(u,v)
      call CalcP2(A, kA, candPar(u,1), CandPar(v,2), 1, LLRX(u,v))   
      if (LLRX(u,v) > 2*TF .and. LLRX(u,v)/=impossible) then
        if (candPar(u,1) < 0) then
          if (ns(-CandPar(u,1),1) <= 1)  cycle   ! CalcPX2 unreliable (?)
        endif
        if (candPar(v,2) < 0) then
          if (ns(-CandPar(v,2),2) <= 1)  cycle  
        endif
        call CalcPX2(A, kA, candPar(u,1), candPar(v,2), LLRY(u,v))
        if (A <0 .and. candPar(u,1) < 0 .and. candPar(v,2) < 0) then
          if (LLRY(u,v) < TF)  PairCand(u,v) = .FALSE.
        else if (LLRY(u,v) < TA) then 
          PairCand(u,v) = .FALSE.
        endif
      else
        PairCand(u,v) = .FALSE.
      endif
    enddo
  enddo
endif

if (withALog)  write(77, '(2i12)', advance='no') COUNT(PairCand)   ! _s2

if (DoLog) then
  open (unit=42,file=trim(logfile),status="unknown", position="append")
  
  ! write(42,*) ""
  ! write(42,*) "SexUnk:"
  ! do m=1,2
    ! write(42,'(500l4)') Sex(CandPar(1:nCP(m),m))
  ! enddo
  
  write(42,*) ""
  write(42,*) "LLRX:"
  do u=1,nCP(1)
    write(42,'(500f8.1)') LLRX(u, 1:nCP(2))
  enddo
  
  write(42,*) ""
  write(42,*) "LLRY:"
  do u=1,nCP(1)
    write(42,'(500f8.1)') LLRY(u, 1:nCP(2))
  enddo
  
  write(42,*) ""
  write(42,*) "PairCand:"
  do u=1,nCP(1)
    write(42,'(500l4)') PairCand(u, 1:nCP(2))
  enddo
  
   write(42,*) ""
 ! write(42,*) "MonoPair:"
 ! do u=1,nCP(1)
   ! write(42,'(500l3)') MonoPair(u, 1:nCP(2))
 ! enddo
  close(42)
endif

  
! age compatibility check
if (AgeAunk .and. ANY(PairCand)) then 
!if (DoLog)  open (unit=42,file=trim(logfile),status="unknown", position="append")                                                                                  
  ALR = missing             
  do u=1, nCP(1)
    call setParTmp(A, kA, CandPar(u,1), 1)
    call SetEstBY(A, kA)
    do v=1, nCP(2)
      if (.not. PairCand(u,v))  cycle
      call CalcAgeLR(A, kA, CandPar(v,2), 2, 0,1, .TRUE., ALR)
!      if (DoLog) write(42, '(2i4, f8.1)') u, v, ALR
      if (ALR == impossible)  PairCand(u,v) = .FALSE.
    enddo
  enddo
  call setParTmp(A, kA, 0, 1)
  call SetEstBY(A, kA)
!  if (DoLog)  close(42)                    
endif


LLRZpair = missing  ! LR(P/next-most-likely)
LLRZsingle = missing
gLL = missing
if (ANY(PairCand)) then  ! possibly parent-pair  
  do u=1, nCP(1)
    do v=1, nCP(2)
      if (.not. PairCand(u,v))  cycle  
      if (Complx==0 .and. A<0 .and. MonoPair(u,v)) then
        call setParTmp(A, kA, 0, 1)
        call CalcPOGPZ(A, kA, CandPar(v,2), 2, gLL) 
        LLRZpair(u,v,:) = gLL(2,:)
        cycle
      endif
      if (Complx==0 .and. A>0) then
        if (EmptySibship(v,2)) then
          call CalcPOGPZ(A, kA, CandPar(u,1), 1, gLL)
        else if (EmptySibship(u,1) .or. (CandPar(u,1) < 0 .and. CandPar(v,2) < 0)) then
          call setParTmp(A, kA, 0, 1)
          call CalcPOGPZ(A, kA, CandPar(v,2), 2, gLL)
        else
          call setParTmp(A, kA, CandPar(u,1), 1)
          call CalcPOGPZ(A, kA, CandPar(v,2), 2, gLL)
        endif  
      else
        call setParTmp(A, kA, CandPar(u,1), 1)
        call CalcPOGPZ(A, kA, CandPar(v,2), 2, gLL)
      endif
!      if (DoLog) then
!        print *, u, v, "; ", Parent(A,:), "; ", CandPar(v,2), "; ", gLL 
!      endif
      LLRZpair(u,v,:) = gLL(3,:)  ! vs one or no parents
      if (gLL(1,AG) < LLRZsingle(AG,u,1)) then
        LLRZsingle(:,u,1) = gLL(1,:)   
      endif
      if (gLL(2,AG) < LLRZsingle(AG,v,2)) then
        LLRZsingle(:,v,2) = gLL(2,:)
      endif
    enddo
  enddo
endif
call setParTmp(A, kA, 0, 1)

TAx = TA
if (A < 0 .and. SUM(nCP)==1) then
 if (ns(-A, kA) == 1)  TAx = 2*TA   ! GGpairs sensitive to false pos.
endif

do u=1, nCP(1)
  do v=1, nCP(2)
    if (.not. PairCand(u,v))  cycle
    do x=1,2
      if (AgePhase == 0 .and. x==2)  cycle   ! only genetics-only needs to pass TA
      if (AgePhase == 2 .and. x==1)  cycle   ! only genetics + age needs to pass TA
      if (LLRZpair(u,v,x) < TAx .or. LLRZpair(u,v,x) >= MaybeOtherParent) then
        PairCand(u,v) = .FALSE.
      endif
    enddo
  enddo
enddo

if (withALog)  write(77, '(2i12)', advance='no') COUNT(PairCand)   ! _s3

if (DoLog) then
  open (unit=42,file=trim(logfile),status="unknown", position="append")
  write(42,*) ""
  write(42,*) "LLRZ-pairs-AG:"
  write(42,'(500i7)') 0, CandPar(1:nCP(2),2)
  do u=1,nCP(1)
    write(42,'(i7, 500f7.2)') CandPar(u,1), LLRZpair(u, 1:nCP(2), AG)
  enddo
  write(42,*) ""
  write(42,*) "LLRZ-pairs-not-AG:"
  write(42,'(500i7)') 0, CandPar(1:nCP(2),2)
  do u=1,nCP(1)
    write(42,'(i7, 500f7.2)') CandPar(u,1), LLRZpair(u, 1:nCP(2), 3-AG)
  enddo
  write(42,*) ""
  write(42,*) "candidate:"
  write(42,'(500i7)') 0, CandPar(1:nCP(2),2)
  do u=1,nCP(1)
    write(42,'(i7, 500l7)') CandPar(u,1), PairCand(u, 1:nCP(2))
  enddo
  write(42,*) ""
  close(42)
endif 


if (ANY(PairCand)) then    ! check if parent-offspring not reversed
  do u=1, nCP(1)
    do v=1, nCP(2)
      if (.not. PairCand(u,v))  cycle
!      if (ParOnly)  cycle    ! reverse check by CheckPair()
      AgeUnk = .FALSE.
      maybeRev = .FALSE.
      dLLrev = missing
      uv = (/u, v/)
      
      do m=1,2
        if (CandPar(uv(m),m) < 0) then
          AgeUnk(m) = .TRUE.
        else if (BY(CandPar(uv(m),m)) < 0) then
          AgeUnk(m) = .TRUE.
        endif
      enddo

      if (AgeAUnk .or. any(AgeUnk)) then                                  
        do m=1,2
          if (CandPar(uv(m),m) < 0) then             
            if (ns(-CandPar(uv(m),m),m) == 0)  cycle   ! A is/was only member of sibship
          endif
          call setParTmp(A, kA, CandPar(uv(3-m), 3-m), 3-m)
          call SetEstBY(A, kA)
          call CheckMaybeRev(A, kA, CandPar(uv(m),m), m, maybeRev(m), dLLrev)
          call setParTmp(A, kA, 0, 3-m)
          call SetEstBY(A, kA)  
!      if (DoLog)  print *, "rev? ", m, maybeRev(m), all(dLLrev==missing), AgePhase, LLRZpair(u,v,AG), dLLrev(AG)
          if (maybeRev(m) .and. ( &
           (all(dLLrev==missing) .or. (AgePhase==1 .and. &
           any(LLRZpair(u,v,:) - dLLrev < 2*TAx)) .or. (AgePhase /= 1 .and. &
           LLRZpair(u,v,AG) - dLLrev(AG) < 2*TAx)))) then
            PairCand(u,v) = .FALSE.
            exit
          endif
        enddo
      endif
    enddo
  enddo
endif

if (withALog)  write(77, '(2i12)', advance='no') COUNT(PairCand)   ! _s4

best = 0
if (ANY(PairCand)) then
  best = MAXLOC(LLRZpair(:,:,AG), MASK=PairCand)
endif 

if (DoLog) then
  open (unit=42,file=trim(logfile),status="unknown", position="append")
  write(42,'("best: ", 2i5, " npairs: ", i4)') best, COUNT(PairCand)
  close(42)
endif  


! check: too small difference in LLR with next-most-likely pair; alternative rels not properly checkable
! NOTE: not entirely sure if this is too conservative
! if (ALL(best > 0) .and. ANY(PairCand)) then
  ! do u=1, nCP(1)
    ! do v=1, nCP(2)
      ! if (.not. PairCand(u,v))  cycle
      ! if (LLRZpair(best(1),best(2),AG) - LLRZpair(u,v,AG) < TAx .and. &
        ! .not. PairCand(best(1), v) .and. .not. PairCand(u, best(2))) then
          ! best = 0
          ! exit
      ! endif
    ! enddo
    ! if (ANY(best == 0))  exit
  ! enddo
  ! if (withALog) then
    ! if (ALL(best >0))  write(77, '(2a12)', advance='no') 'PASS'
    ! if (ANY(best==0))  write(77, '(2a12)', advance='no') 'FAIL'
  ! endif
! else if (withALog) then
  ! write(77, '(2a12)', advance='no') 'NA' 
! endif


if (ALL(best > 0)) then    ! assign (grand)parent pair    
  do m=1,2
    call setParTmp(A, kA, CandPar(best(m), m), m)
  enddo
  nowParUV = best              

  ! ~~~~  check when >1 eligible pair  ~~~~
  if (COUNT(PairCand) > 1) then                                       
    do v=1, nCP(2)
      do u=1, nCP(1)
        if (nCP(1)==0 .and. nCP(2)==2 .and. v==2)  cycle  ! pair already checked
        if (nCP(1)==2 .and. nCP(2)==0 .and. u==2)  cycle
        if (MonoPair(u,v) .and. PairCand(u,v) .and. u/=nowParUV(1)) then   ! MonoPair(NowParUV(1), NowParUV(2)) .and. 
          call CalcPOGPZ(A, kA, CandPar(u,1), 1, gLL)
          LLRZpair(NowParUV(1),NowParUV(2),AG) = MIN(LLRZpair(NowParUV(1),NowParUV(2),AG), gLL(4,AG))                   
          if (gLL(4,AG) < TAx) then  ! nowParUV(1) + (2) not most likely
            call setParTmp(A, kA, CandPar(u,1), 1)
            call setParTmp(A, kA, CandPar(v,2), 2)
            call CalcPOGPZ(A, kA, CandPar(NowParUV(1),1), 1, gLL)
            LLRZpair(u,v,AG) = MIN(LLRZpair(u,v,AG), gLL(4,AG))            
            nowParUV(1) = u
            nowParUV(2) = v
          endif
          cycle
        endif  
        
        if (v/=NowParUV(2) .and. v/=best(2) .and. PairCand(nowParUV(1), v)) then
          call CalcPOGPZ(A, kA, CandPar(v,2), 2, gLL)
          LLRZpair(NowParUV(1),NowParUV(2),AG) = MIN(LLRZpair(NowParUV(1),NowParUV(2),AG), gLL(4,AG))
          LLRZpair(NowParUV(1),v,AG) = MIN(LLRZpair(NowParUV(1),v,AG), gLL(3,AG))
          if (gLL(3,AG) > gLL(4,AG)) then   ! nowParUV(1) + CandPar(v,2) more likely than nowParUV(1) + (2)
            call setParTmp(A, kA, CandPar(v,2), 2)
            nowParUV(2) = v           
          endif  
        endif
        
        if (PairCand(u, nowParUV(2)) .and. u/=nowParUV(1) .and. u/=best(1)) then
          call CalcPOGPZ(A, kA, CandPar(u,1), 1, gLL)         
          LLRZpair(NowParUV(1),NowParUV(2),AG) = MIN(LLRZpair(NowParUV(1),NowParUV(2),AG), gLL(4,AG))
          ! TODO: shouldn't use MIN() if newer based on trio while old based on duo. 
          ! - reset all LLRZpair at start of this loop?
          LLRZpair(u,NowParUV(2),AG) = MIN(LLRZpair(u,NowParUV(2),AG), gLL(3,AG))
          if (gLL(3,AG) > gLL(4,AG)) then
            call setParTmp(A, kA, CandPar(u,1), 1)
            nowParUV(1) = u
          endif  
        endif        
      enddo
    enddo
  endif

  if (LLRZpair(NowParUV(1),NowParUV(2),AG) < TAx) then 
    ! newly assigned pair not convincingly more likely than alternatives
    do m=1,2                                                                          
      call setParTmp(A,kA, 0,m)          
    enddo
  endif
endif

if (withALog) then
  if (any(LLRZpair(:,:,AG) < 222.0)) then
    write(77, '(f12.1)', advance='no')  MAXVAL(LLRZpair(:,:,AG), MASK=LLRZpair(:,:,AG) < 222.0)
  else
    write(77, '(f12.1)', advance='no')  -Missing
  endif
endif

curPar = getPar(A,kA)
if (withALog)  write(77, '(2X, 2a20)', advance='no')  id(curpar)

if (ALL(curPar/=0)) then  ! make it official
  do m=1,2
    call setPar(A,kA, curPar(m),m)
  enddo
  if (DoLog) then
    open (unit=42,file=trim(logfile),status="unknown", position="append")
    write(42,*) ""
    write(42,*) "Assigned: ", getPar(A,kA)
    write(42,*) ""
    close(42)
  endif
  return
else
  do m=1,2
    call setParTmp(A,kA, 0,m)
  enddo
endif


DoSingle = .TRUE.
if (hermaphrodites/=0 .and. LRS > TAx) then  ! selfed
  DoSingle = .FALSE.
else if (Complx==0 .and. A<0) then
  DoSingle = .FALSE.
  do m=1,2
    do u = 1, nCP(m)
      if (CandPar(u,m) > 0) then  
        if (Mate(CandPar(u,m)) == 0) then   ! exception: no genotyped offspring --> no mate yet
          DoSingle = .TRUE.
          exit
        endif
      endif  ! CandPar(u,m) < 0 w/o DumMate shouldn't exist
    enddo
    if (DoSingle)  exit
  enddo
endif

!if (withALog)  write(77, '(2X, l10)', advance='no')  DoSingle

if (.not. DoSingle) then
  do m=1,2
    call setPar(A,kA, 0,m)
  enddo
  if (DoLog) then
    open (unit=42,file=trim(logfile),status="unknown", position="append")
    write(42,*) ""
    write(42,*) "Assigned: ", getPar(A,kA)
    write(42,*) ""
    close(42)
  endif
  return
endif


! ~~~~  single parent  ~~~~
! if (A < 0) then
  ! if (ns(-A, kA) == 1)  TAx = 3*TA   ! GGpairs sensitive to false pos.
! endif

! duplicate check among candidate parents.
if (withALog) then
  mxDupLR = -Missing
  dupi = 0
  dupm = 0
  do m=1,2
    if (nCP(m) < 2)  cycle  ! assume sex is correct or missing
    do u=1,nCP(m)-1
      do v=u+1, nCP(m)
        call CheckMaybeDup(CandPar(u,m), CandPar(v,m), DupLR)
        if (DupLR > mxDupLR) then
          mxDupLR = DupLR
          dupi = (/CandPar(u,m), CandPar(v,m)/)
          dupm = m
        endif
      enddo
    enddo
  enddo
  write(77, '(f8.1, i8, 2X, 2a20)', advance='no')  mxDupLR, dupm, ID(dupi)
endif

SingleCand = .FALSE.
do m=1,2
  do u=1, nCP(m)
    if (SexUnk(u,m) .and. hermaphrodites/=2)  cycle
    if (hermaphrodites==2 .and. m==2) then
      if (CandPar(u,m)>0) then
        if (Sex(CandPar(u,m))==4)  cycle  ! assign as dam
      else if (any(CandPar(1:nCP(1),1) == -DumClone(-CandPar(u,2),2))) then
        cycle ! assign as dam
      endif
    endif
    
    if (ALL(LLRZsingle(:,u,m) < TAx))  cycle

    if (Complx==0) then
      if (CandPar(u,m) < 0) then
        if (DumMate(-candpar(u,m),m)/=0)  cycle
      else
        if (Mate(CandPar(u,m))/=0)  cycle
      endif
    endif

    SingleCand(u,m) = .TRUE.
  enddo
enddo

if (DoLog) then
  open (unit=42,file=trim(logfile),status="unknown", position="append")
  write(42,*) ""
  write(42,*) "LLR-singles-1"
  do m=1,2
    do u=1,nCP(m)
      write(42,'(i3, i6, 2f8.2, l5)') m, CandPar(u,m), LLRZsingle(:,u,m), SingleCand(u,m)
    enddo
  enddo
  write(42,*) ""
  close(42)
endif  

!if (DoLog)  print *, 'nCand: ', SUM(nCP), COUNT(SingleCand)
if (withALog)  write(77, '(i12)', advance='no')  COUNT(SingleCand)   ! _s1

nSingle = COUNT(SingleCand)  ! SingleCand edited inside loop

if (nSingle >= 1) then
  DoneSingleCheck = .TRUE.
  do m=1,2
    do u=1,nCP(m)
      if (SingleCand(u,m))  DoneSingleCheck(u, m) = .FALSE.
    enddo
  enddo

  if (DoLog)  open (unit=42,file=trim(logfile),status="unknown", position="append")
  if (DoLog)  write(42,*) ""
  
  do m=1,2
    do u=1,nCP(m)
      if (.not. SingleCand(u,m))  cycle
      i = (m-1)*mxCP +u    
      LRStmp = Missing      
      
      if (nSingle==1 .or. Complx==0) then   ! do not also condition on other CP being parent
        call CalcCandParLL(A, kA, CandPar(u,m), m, LLA)
        call LLA2LLR(LLA, fcl, LRStmp)
        if (DoLog)  write(42,'("single?? ", i3, i5, 2f7.2)') m, CandPar(u,m), LRStmp(:,1)                              
        if (LRStmp(AG,1) < LLRZsingle(AG,u,m))  LLRZsingle(:,u,m) = LRStmp(:,1)
        if (ALL(LLRZsingle(:,u,m) < TAx))  SingleCand(u,m) = .FALSE.    
      else
        do n=1,2 
          do v=1,nCP(n)
            j = (n-1)*mxCP + v
            if (j>=i) exit  ! only do lower triangle of all singles X all singles matrix
            if (.not. SingleCand(v,n) .and. &
           .not. (COUNT(SingleCand)==1 .and. .not. ALL(DoneSingleCheck)))  cycle                       
            if (m==1 .and. n==2) then
              if (PairCand(u,v))  cycle  ! already done
            endif
            if (m==2 .and. n==1) then
              if (PairCand(v,u))  cycle  ! already done
            endif
            if (CandPar(u,m) == CandPar(v,n) .and. CandPar(u,m)>0)  cycle  ! when sex=4 
            call setParTmp(A, kA, CandPar(v,n), n) 
            call CalcCandParLL(A, kA, CandPar(u,m), m, LLA)
            call LLA2LLR(LLA, fcl, LRStmp)   
            
            if (DoLog)  write(42,'("single? ", i3, i6, 2f7.2, "; ", i3, i6, 2f7.2)') m, CandPar(u,m), LRStmp(:,1), &
              n, CandPar(v,n), LRStmp(:,2)       
            
            if (LRStmp(AG,1) < LLRZsingle(AG,u,m))  LLRZsingle(:,u,m) = LRStmp(:,1)
            if (LRStmp(AG,2) < LLRZsingle(AG,v,n))  LLRZsingle(:,v,n) = LRStmp(:,2)
            if (ALL(LLRZsingle(:,u,m) < TAx))  SingleCand(u,m) = .FALSE.
            if (ALL(LLRZsingle(:,v,n) < TAx))  SingleCand(v,n) = .FALSE.
          enddo  ! v     
          call setParTmp(A, kA, 0, n)
        enddo  ! n
      endif
      
    enddo  ! u
  enddo  ! m
  
  if (DoLog) write(42,*) ""
  if (DoLog) close(42)
endif  


do m=1,2
  do u=1,nCP(m)
    do x=1,2
      if (AgePhase == 0 .and. x==2)  cycle   ! only genetics-only needs to pass TA
      if (AgePhase == 2 .and. x==1)  cycle   ! only genetics + age needs to pass TA 
      if (LLRZsingle(x,u,m) < TAx .or. LLRZsingle(x,u,m) >= 222D0) then
        SingleCand(u,m) = .FALSE.
      endif
    enddo
  enddo
enddo
     
if (withALog)  write(77, '(i12)', advance='no')  COUNT(SingleCand)   ! _s2
     
      
! check if PO pair could be flipped, or relies strongly on uncertain age estimate
dLLrev = missing 
do m=1,2
  do u=1,nCP(m)
    if (.not. SingleCand(u,m))  cycle
    call CheckMaybeRev(A, kA,CandPar(u,m), m, maybeRev(1), dLLrev)  ! includes age check
    if (DoLog) then
      open (unit=42,file=trim(logfile),status="unknown", position="append")
        write(42,'(2i5, " rev? ", l4, 2f7.2)') CandPar(u,m), m, maybeRev(1), dLLrev
      close(42)
    endif  
    if (maybeRev(1)) then
      if (all(dLLrev==missing) .or. (AgePhase==1 .and. &    
       any(dLLrev - LLRZsingle(:,u,m) > TF)) .or. (AgePhase /= 1 .and. &
       dLLrev(AG) - LLRZsingle(AG,u,m) > TF) .or. & ! prone to reversal --> strict threshold. 
       ALL(dLLrev > 2*TA)) then     ! new 2023-02-14 
        SingleCand(u,m) = .FALSE.
        cycle
      endif
    endif
  enddo
enddo

if (withALog)  write(77, '(i12)', advance='no')  COUNT(SingleCand)   ! _s3

if (DoLog) then
  open (unit=42,file=trim(logfile),status="unknown", position="append")
  write(42,*) ""
  write(42,*) "LLR-singles"
  do m=1,2
    do u=1,nCP(m)
      write(42,'(i3, i6, 2f8.2, l5)') m, CandPar(u,m), LLRZsingle(:,u,m), SingleCand(u,m)
    enddo
  enddo
  write(42,*) ""
  close(42)
endif  

best=0
best = MAXLOC(LLRZsingle(AG,:,:), MASK=SingleCand)  ! u,m

!if (DoLog)  print *, 'count: ', COUNT(SingleCand), ', best: ', best(1)

if (COUNT(SingleCand) == 1) then
  call setPar(A, kA, CandPar(best(1), best(2)), best(2))
  if (Complx == 0 .and. A>0 .and. .not. ParOnly) then  ! create dummy mate if assigned parent has no mate yet
    call NewSibship(A, 0, 3-best(2))
  endif
else
  do m=1,2
    call setParTmp(A,kA, 0,m)   ! do no funny bussiness w temp assigned ex-cand-parents
  enddo
  do m=1,2
    call setPar(A,kA, 0,m)
  enddo 
endif

if (withALog) then
  if (any(LLRZsingle(AG,:,:)<222.0)) then
    write(77, '(f14.1)', advance='no')  MAXVAL(LLRZsingle(AG,:,:), MASK=LLRZsingle(AG,:,:)<222.0) 
  else
    write(77, '(f14.1)', advance='no')  -Missing
  endif
  write(77, '(2X, 2a20)', advance='no')  id(getPar(A,kA))
endif

if (DoLog) then
  open (unit=42,file=trim(logfile),status="unknown", position="append")
  write(42,*) ""
  write(42,*) "Assigned: ", getPar(A,kA)
  write(42,*) ""
  close(42)
endif

end subroutine SelectParent

! ######################################################################

subroutine CalcPOGPZ(A, kA, B, kB, pLLR)
use Global
implicit none

integer, intent(IN) :: A, kA, B, kB
double precision, intent(OUT) :: pLLR(4,2)  ! dam - sire - new pair - old pair , w/o - w age
integer :: AG, fcl, notfcl(6), x, curpar(2)
double precision :: LLA(7,7,2,2), LLS, LLP(4,2,2) 
logical :: ParClone

pLLR = missing
curpar = getPar(A, kA)
if (curPar(kB) == B)  return

ParClone = .FALSE.
if (Hermaphrodites/=0 .and. A>0 .and. curPar(kB)==0 .and. &
 B>0 .and. curpar(3-kB) == B) then  
  ParClone = .TRUE.
  call setParTmp(A, kA, 0, 3-kB)
endif

LLA = missing            
call CalcCandParLL(A, kA, B, kB, LLA) 

if (A > 0) then
  fcl = 1
  notfcl = (/2,3,4,5,6,7/)
else 
  fcl = 4
  notfcl = (/1,2,3,5,6,7/)
endif

if (ParClone) then  ! A>0 .and. B >0
  call setParTmp(A, kA, curPar(3-kB), 3-kB)
  call PairPO(A, B, kB, 1, LLS)   
  LLA(fcl,fcl,3-kB,:) = LLS   ! TODO: ALR?
endif

LLP = missing
do AG=1,2
  if (curPar(kB)==0) then   ! --> LLA(:,:,kB,:) empty
    ! B + curPar(3-kB)
    LLP(3   ,1,AG) = LLA(fcl,fcl,3-kB,AG)
    LLA(fcl,fcl,3-kB,AG) = 555D0
    LLP(3   ,2,AG) = MaxLL(RESHAPE(LLA(:,:,3-kB,AG), (/7*7/)))
    ! only B
    LLP(kB  ,1,AG) = MaxLL(LLA(fcl,notfcl,3-kB,AG)) 
    LLP(kB  ,2,AG) = MaxLL(RESHAPE(LLA(notfcl,:,3-kB,AG), (/6*7/))) 
    ! only curPar(3-kB) 
    LLP(3-kB,1,AG) = MaxLL(LLA(notfcl,fcl,3-kB,AG))
    LLP(3-kB,2,AG) = MaxLL(LLA(fcl,:,3-kB,AG)) 
  else if (all(curPar /= 0)) then
    ! B + curPar(3-kB)
    LLP(3   ,1,AG) = MaxLL((/LLA(fcl,fcl,3-kB,AG), LLA(fcl,:,kB,AG)/))
    LLP(3   ,2,AG) = MaxLL((/RESHAPE(LLA(notfcl,:,:,AG), (/6*7*2/)), &
      LLA(fcl,notfcl,3-kB,AG)/))
    ! curPar pair                             
    LLP(4   ,1,AG) = MaxLL(RESHAPE(LLA(notfcl,fcl,:,AG), (/6*2/)))
    LLP(4   ,2,AG) = MaxLL((/RESHAPE(LLA(:,notfcl,:,AG),(/7*6*2/)), &
      LLA(fcl,fcl,3-kB,AG)/))
    LLA(fcl,fcl,:,AG) = 555D0     
    ! only B
    LLP(kB  ,1,AG) = MaxLL(LLA(fcl,notfcl,3-kB,AG))
    LLP(kB  ,2,AG) = MaxLL((/RESHAPE(LLA(notfcl,:,:,AG),(/2*6*7/)),&
      LLA(fcl,:,kB,AG), LLA(fcl,fcl,3-kB,AG)/)) 
    ! only curPar(3-kB)
    LLP(3-kB,1,AG) = MaxLL(RESHAPE(LLA(notfcl,notfcl,kB,AG), (/6*6/)))
    LLP(3-kB,2,AG) = MaxLL((/RESHAPE(LLA(fcl,:,:,AG), (/2*7/)), &
      RESHAPE(LLA(:,fcl,:,AG), (/2*7/))/))
  else 
    LLP(kB  ,1,AG) = LLA(fcl,7,3-kB,AG)
    LLP(kB  ,2,AG) = MaxLL(LLA(notfcl,7,3-kB,AG))
  endif
enddo

pLLR = 555D0
do x=1,4
  if (all(LLP(x,:,:) < 0))  pLLR(x,:) = LLP(x,1,:) - LLP(x,2,:)
enddo  

if (Complx==0 .and. A>0 .and. all(curPar==0)) then  ! monogamous
  pLLR(3,:) = pLLR(kB,:) 
endif


! if (A==2969) then
! if (A<0) then
  ! if (ANY(SibID(:,-A,kA)==91)) then
  ! open (unit=69,file="log-CalcCPL.txt",status="unknown", position="append")
  ! write(69,*) ">>>>>>>>>>>>>>>>>>"
  ! write(69, '("CalcPOGPZ: ", 2i6, " (", 2i6, ") + ", 3i6)') A, kA, curpar, B, kB, fcl
  ! do x=1,2
    ! write(69, '("LLP :", 4f9.1)') LLP(:,x,1)
  ! enddo
  ! write(69, '("pLLR:", 4f9.1)') pLLR(:,1)
  ! write(69,*) ">>>>>>>>>>>>>>>>>>"
  ! write(69,*) ""
  ! close(69)
! endif
! endif

end subroutine CalcPOGPZ

! ######################################################################

subroutine LLA2LLR(LLA, fcl, LLR)   ! output array from CalcCandParLL -> LR(PO/not)
use Global
implicit none

double precision, intent(IN) :: LLA(7,7,2,2)
integer, intent(IN) :: fcl
double precision, intent(OUT) :: LLR(2,2)  ! UseAge; vs not-PO / vs U  ;B, curPar (kB or 3-kB)
integer :: notfcl(6), a
double precision :: LLPO

if (fcl==1) then   ! A>0
  notfcl = (/2,3,4,5,6,7/)
else !if (fcl==4) then  ! A<0
  notfcl = (/1,2,3,5,6,7/)
endif

LLR = Missing
do a=1,2
  if (ALL(LLA(fcl,7,:,a) > 0)) then   ! maybeOtherParent / Impossible w/o co-parent
    LLR(a,1) = -MINVAL(LLA(fcl,7,:,a))
  else
    LLPO = MaxLL( RESHAPE(LLA(fcl,:,:,a), (/2*7/)) )   ! B parent
    LLR(a,1) = LLPO - MaxLL( RESHAPE(LLA(notfcl,:,:,a), (/2*6*7/)) )
  endif

  if (ALL(LLA(7,fcl,:,a) > 0)) then  
      LLR(a,2) = -MINVAL(LLA(7,fcl,:,a))
  else
    LLPO = MaxLL( RESHAPE(LLA(:,fcl,:,a), (/2*7/)) ) ! curPar parent
    LLR(a,2) = LLPO - MaxLL( RESHAPE(LLA(:,notfcl,:,a), (/2*6*7/)) )
  endif
enddo

end subroutine LLA2LLR

! ######################################################################

subroutine CalcCandParLL(A, kA, B, kB, LLA)
use Global
implicit none

! calc LL over A, candidate parent B, + current parent of A, under a set of diff relationships
! combination of prev. CalcPOZ & CalcGPZ

integer, intent(IN) :: A, kA, B, kB
double precision, intent(OUT) :: LLA(7,7,2,2)  ! B, curPar(3-n), curPar(n), Age
integer :: focal, m, curPar(2), mid(5), parB(2), r, fclx(3), z, curGP(2,2)
double precision :: LLcp(3,2), LLU(4), U, ALR(3), ALRtmp(2), LLtmp, LLtrio(6), &
   LLX(2,2), ALRtrio(2,6)
logical :: ParOK!, DoTrios

LLA = missing
curPar = getPar(A, kA)
curGP=0       
if (curPar(kB) == B)  return

if (A > 0) then
  focal = 1
  mid = (/2,3,4,5,6/)
else if (A < 0) then
  focal = 4
  mid = (/1,2,3,5,6/)
else
  return
endif

call CalcU(A, kA, 0, 0, U)

if (ALL(curPar==0)) then
  call CheckRel(A, kA, B, kB, focal, LLA(:,7,3-kB,1), LLA(:,7,3-kB,2))

else
  LLCP = 0D0
  LLU = missing
  ALR = 0D0
  ParOK = .TRUE.  
  do m=1,2
    curGP(:,m) = getPar(curPar(m), m)
  enddo
  parB = getPar(B, kB)

  call Calc4U(curPar, B,kB, A,kA, LLU, LLCP)
  ! LLU: curPar(1), curPar(2), B, A  ; 1-3 w/o A if <0
  ! LLCP: correction factors                          
  ! LLoverlap: likelihoods of overlaps/due to connections                                                               
  
  do m=1,2
    call CalcAgeLR(A,kA, curPar(m), m, 0,1, .FALSE., ALR(m))
    if (A>0 .and. curPar(m)<0) then  
      fclx(m) = 6
    else
      fclx(m) = focal
    endif
  enddo 
  call CalcAgeLR(A,kA, B, kB, 0,1, .FALSE., ALR(3))
  if (A>0 .and. B<0) then  
    fclx(3) = 6
  else 
    fclx(3) = focal
  endif 
  
  do m=1,2 ! sex currently assigned par
    if (curPar(m)==0) cycle 
    call checkRel(A, kA, B, kB, fclx(3), LLA(:,focal,m,1), LLA(:,focal,m,2))   ! curPar(m)=GP + A_7 

    call setParTmp(A, kA, 0, m)
    call checkRel(A, kA, B, kB, fclx(3), LLA(:,7,m,1), LLA(:,7,m,2))   ! A_7
    call checkRel(A, kA, curPar(m), m, fclx(m), LLA(7,:,m,1), LLA(7,:,m,2))  ! curPar(m)_7
    if (curPar(m) < 0 .and. focal==4) then
      if (m/=kB .and. ANY(Parent(SibID(1:nS(-curPar(m),m), -curPar(m),m), kB)==B)) then
        call PairUA(A, curPar(m), kA, m, LLA(7,focal,m,1))
        LLA(7,focal,m,2) = LLA(7,focal,m,1) 
      endif
    endif
    
    call ChkValidPar(A, kA, B, kB, ParOK)
    if (ParOK) then 
      call setParTmp(A, kA, B, kB)
      call checkRel(A, kA, curPar(m),m, fclx(m), LLA(focal,:,m,1), LLA(focal,:,m,2)) 
    endif

    call setParTmp(A, kA, curPar(kB), kB)
    call setParTmp(A, kA, curPar(3-kB), 3-kB)
      
    ! fix: B or curpar parent as side-effect of curpar resp B as HS    
    if (A>0 .and. curPar(m)/=0) then
      if (curPar(m) == parB(m))  LLA(2:3,7,m,:) = impossible  ! HS implies curPar = Par
      if (curGP(kB,m) == B)      LLA(7,2:3,m,:) = impossible  ! HS implies B = Par
      if (curPar(3-m)/=0) then     
        if (curPar(3-m) == parB(3-m))  LLA(2,7,m,:) = impossible
        if (curGP(kB,3-m) == B)        LLA(7,2,m,:) = impossible
      endif
    endif
    if (Complx==0) then
      LLA(focal,7,m,:) = impossible
      LLA(7,focal,m,:) = impossible
    endif

    do z=1,2
      WHERE (LLA(mid,focal,m,z)<0) LLA(mid,focal,m,z) = LLA(mid,focal,m,z) +LLcp(3,m)
      WHERE (LLA(mid,7,m,z)<0)     LLA(mid,7,m,z) = LLA(mid,7,m,z) + LLcp(3,m)
      WHERE (LLA(focal,:,m,z)<0)   LLA(focal,:,m,z) = LLA(focal,:,m,z) + LLcp(m,m)
      WHERE (LLA(7,:,m,z)<0)       LLA(7,:,m,z) = LLA(7,:,m,z) + LLcp(m,m)
    enddo
    
    WHERE (LLA(mid,focal,m,2)<0) LLA(mid,focal,m,2) = LLA(mid,focal,m,2) +ALR(m)  
    WHERE (LLA(focal,:,m,2)<0)   LLA(focal,:,m,2) = LLA(focal,:,m,2) + ALR(3)
    WHERE (LLA(7,:,m,2)<0)       LLA(7,:,m,2) = LLA(7,:,m,2) + ALR(3-m) 
    
    if (m == kB) then
      WHERE (LLA(focal,mid,m,2)<0) LLA(focal,mid,m,2) = LLA(focal,mid,m,2) + ALR(3-m)
      LLA(focal,focal,m,:) = impossible ! cannot have 2 same-sex parents
    endif
    
    if (A<0 .and. ParB(m)==curPar(m)) then  ! then 'B FS of A' identical with & without curpar
      LLA(5,7,m,:) = impossible
    endif

    LLtrio = Missing
    ALRtrio = 0D0
!    DoTrios = .FALSE.    
    if (B/=curPar(m) .and. parB(m)/=curPar(m) .and. curGP(m,m)/=B .and. &
     .not. any(ParB == curGP(:,m) .and. parB/=0)) then
      ! DoTrios = .TRUE.      
      ! ! only do trio calcs if either candidate has LR(fcl/not) > 0 
        ! this is slower than just calculating the trioLL's  ??
      ! LLP(1,1) = MaxLL( RESHAPE(LLA(focal,:,:,2), (/2*7/)) )   ! B parent
      ! LLP(1,2) = MaxLL( RESHAPE(LLA((/mid,7/),:,:,2), (/2*6*7/)) )  ! B not parent
      ! LLP(2,1) = MaxLL( RESHAPE(LLA(:,focal,:,2), (/2*7/)) )   ! curPar parent
      ! LLP(2,2) = MaxLL( RESHAPE(LLA(:,(/mid,7/),:,2), (/2*6*7/)) )  ! curpar not parent
      ! if ((LLP(1,1) < LLP(1,2)) .and. (LLP(2,1) < LLP(2,2))) then
        ! DoTrios = .FALSE.
      ! else
        ! DoTrios = .TRUE.
      ! endif
    ! endif    
    ! if (DoTrios) then
      ! check if B & curPar(m) both FS / both GP of A, + HS & FA if all >0
      call setParTmp(A, kA, 0, m)  ! for ageLR 
      do r=2,6
        if (r==2 .and. A<0 .and. (any(parB/=0) .or. any(curGP(:,m)/=0) .or. curPar(3-m)/=0))  cycle
        if (r==3 .and. (Complx==0 .or. any(parB/=0) .or. any(curGP(:,m)/=0) .or. curPar(3-m)/=0))  cycle
        if (r==4 .and. A<0)  cycle
        if (r==5 .and. (B<0 .or. curPar(m)<0))  cycle      
        if (r/=6) then
          call CalcAgeLR(A,kA, B,kB, kB,r, .FALSE., ALRtrio(1,r))
          call CalcAgeLR(A,kA, CurPar(m),m, kB,r, .FALSE., ALRtrio(2,r))
          if (any(ALRtrio(:,r)==impossible .or. ALRtrio(:,r) < 5*TF)) cycle
        endif
        if (r==2)  call trioFS(A,kA, B,kB, curPar(m),m, LLtrio(2))
        if (r==3)  call trioHS(A,kA, B,kB, curPar(m),m, LLtrio(3))  ! B & curPar HS to eachother or not.
        if (r==4)  call trioGP(A,kA, B,kB, curPar(m),m, LLtrio(4))  ! A>0 only
        if (r==5)  call trioFA(A,kA, B,curPar(m), LLtrio(5))  ! B & C >0 only
        if (r==6) then
          if (A<0) then
            call CalcAgeLR(A,kA, B,kB, kB,4, .FALSE., ALRtrio(1,r))
            call CalcAgeLR(A,kA, CurPar(m),m, kB,4, .FALSE., ALRtrio(2,r))
            if (any(ALRtrio(:,r)==impossible .or. ALRtrio(:,r) < TF)) then
              ALRtrio(:,r) = 0D0   ! use as proxy for other 3rd degree rels
            endif
            call trioGP(A,kA, B,kB, curPar(m),m, LLtrio(6))           
          else if (B>0 .and. curPar(m)>0) then  
            if (all(ParB/=0) .or. all(curGP(:,m)/=0)) then
              call trioGGP(A,kA, B,curPar(m), LLtrio(6))
              call CalcAgeLR(A,kA, B,kB, kB,4, .FALSE., ALRtrio(1,r))
              call CalcAgeLR(A,kA, CurPar(m),m, kB,4, .FALSE., ALRtrio(2,r))
              do z=1,2  ! TODO: ALR for GGP 
                if (ALRtrio(z,r)==impossible .or. ALRtrio(z,r) < 3*TF)  ALRtrio(z,r)=0.0D0
              enddo
            else
              call CalcAgeLR(A,kA, B,kB, kB,5, .FALSE., ALRtrio(1,r))
              call CalcAgeLR(A,kA, CurPar(m),m, kB,5, .FALSE., ALRtrio(2,r)) 
              if (any(ALRtrio(:,r)==impossible .or. ALRtrio(:,r) < 5*TF)) cycle              
              call trioHA(A,kA, B,curPar(m), LLtrio(6))                         
            endif 
          else
            LLtrio(6) = NotImplemented   ! (yet)
            cycle
          endif
        endif  
        if (LLtrio(r) < -HUGE(0D0) .or. LLtrio(r)>=NotImplemented)  cycle
        LLA(r,r,m,:) = LLtrio(r) + LLU(3-m)! &  ! curPar(3-m)
 !       + SUM(LLoverlap)/2   ! trio subroutines assume all 3 are independent. BREAKS EVERYTHING
 !         + SUM(LLoverlap(:,1))   ! overlap with A only?
          LLA(r,r,m,2) = LLA(r,r,m,2) + ALRtrio(1,r) + ALRtrio(2,r)
      enddo
      ! FS+HS
      ! if (all(parB==0) .and. all(curGP(:,m)==0) .and. curPar(3-m)==0) then
        ! call trioFHS(A,kA, B,kB, curPar(m),m, LLtmp)  ! B FS, curPar(m) HS
        ! call CalcAgeLR(A,kA, B,kB, kB,2, .FALSE., ALRtmp(1))
        ! call CalcAgeLR(A,kA, CurPar(m),m, kB,3, .FALSE., ALRtmp(2))
        ! LLA(2,3,m,1) = LLtmp + LLU(3-m)
        ! LLA(2,3,m,2) = LLtmp + LLU(3-m) + ALRtmp(1) + ALRtmp(2)
        
        ! call trioFHS(A,kA, curPar(m),m, B,kB, LLtmp)  ! B FS, curPar(m) HS
        ! call CalcAgeLR(A,kA, B,kB, kB,3, .FALSE., ALRtmp(1))
        ! call CalcAgeLR(A,kA, CurPar(m),m, kB,2, .FALSE., ALRtmp(2))
        ! LLA(3,2,m,1) = LLtmp + LLU(3-m)
        ! LLA(3,2,m,2) = LLtmp + LLU(3-m) + ALRtmp(1) + ALRtmp(2)     
      ! endif
      
      ! HS + GP
      LLX = Missing
      if (curPar(3-m)==0 .and. Complx>0) then   ! A<0  ! TODO? curPar(3-m)/=0
        do z=1,2
!          if (B>0) then
            call CalcAgeLR(A,kA, B,kB, z,3, .FALSE., ALRtmp(1))
            call CalcAgeLR(A,kA, CurPar(m),m, z,4, .FALSE., ALRtmp(2)) 
            if (.not. any(ALRtmp==impossible .or. ALRtmp < 5*TF)) then
              call trioHSGP(A,kA, B,kB, CurPar(m),m, z, LLX(1,z))
!              call trioHSGP2(A,kA, B,kB, CurPar(m),m, z, LLX(3,z))
            endif
!          endif
!          if (curPar(m)>0) then
            call CalcAgeLR(A,kA, CurPar(m),m, z,3, .FALSE., ALRtmp(1))
            call CalcAgeLR(A,kA, B,kB, z,4, .FALSE., ALRtmp(2))  
            if (.not. any(ALRtmp==impossible .or. ALRtmp < 5*TF)) then
              call trioHSGP(A,kA, CurPar(m),m, B,kB, z, LLX(2,z))
!              call trioHSGP2(A,kA, CurPar(m),m, B,kB, z, LLX(4,z))
            endif
!          endif         
        enddo
        LLA(6,5,m,:) = MaxLL((/LLX(:,1), LLX(:,2)/))                                              
      endif      
      call setParTmp(A, kA, curPar(m), m)     
    endif
        
    if (kA>2)  cycle
    ! check if B is offspring & curpar(m) parent
    if ((m/=kB .or. any(curpar==0)) .and. parB(kA)==0 .and. (A>0 .or. B>0)) then  ! A>0 .and. B>0 .and. curPar(m)<0 .and.              
!     call setParTmp(A, kA, curPar(m), m)  ! already is parent      
      call ChkValidPar(B, kB, A, kA, ParOK)
      if (ParOK) then 
        if (A>0) then
          call CalcU(B,kB, curPar(m),m, LLX(1,1)) 
        else
          call CalcU(A,kA, curPar(m),m, LLX(1,1))
        endif
        call setParTmp(B, kB, A, kA)
        if (A>0) then
          call CalcU(B,kB, curPar(m),m, LLX(2,1)) 
        else
          call CalcU(A,kA, curPar(m),m, LLX(2,1))
        endif
        LLtmp = LLX(2,1) - LLX(1,1) + LLA(7,1,m,1)   ! change in LL relative to curpar parent, B unrelated
        if (LLtmp > LLA(6,1,m,1) .or. LLA(6,1,m,1)>0) then
          LLA(6,1,m,1) = LLtmp   ! B not 3rd degree rel, but convenient spot
          call CalcAgeLR(B,kB, A,kA, 0,1, .FALSE., ALRtmp(1)) 
          LLA(6,1,m,2) = LLA(6,1,m,1) + ALRtmp(1) + ALR(m)
        endif
        call setParTmp(B, kB, 0, kA)
      endif
    endif
    
    ! check if B is parent & curpar(m) offspring
    if ((m/=kB .or. any(curpar==0)) .and. curGP(kA,m)==0 .and. (A>0 .or. B>0)) then  ! A>0 .and. B<0 .and. curPar(m)>0 .and. 
      ! curpar(m) - A - B
      call setParTmp(A, kA, 0, m)
      call setParTmp(A, kA, B, kB)
      call ChkValidPar(curPar(m),m, A,kA, ParOK)
      if (ParOK) then              
        if (A>0) then
          call CalcU(B,kB, curPar(m),m, LLX(1,1)) 
        else
          call CalcU(A,kA, curPar(m),m, LLX(1,1))
        endif
        call setParTmp(curPar(m), m, A, kA)  
        if (A>0) then
          call CalcU(B,kB, curPar(m),m, LLX(2,1)) 
        else
          call CalcU(A,kA, curPar(m),m, LLX(2,1))
        endif    
        LLtmp = LLX(2,1) - LLX(1,1) + LLA(1,7,m,1)   ! change in LL relative to B parent, curpar unrelated        
        if (LLtmp > LLA(1,6,m,1) .or. LLA(1,6,m,1)>0) then
          LLA(1,6,m,1) = LLtmp   ! B not 3rd degree rel, but convenient spot
          call CalcAgeLR(curPar(m),m, A,kA, 0,1, .FALSE., ALRtmp(1)) 
          LLA(1,6,m,2) = LLA(1,6,m,1) + ALRtmp(1) + ALR(3)
        endif
        call setParTmp(curPar(m), m, 0, kA)
      endif
      call setParTmp(A, kA, curPar(kB), kB)
      call setParTmp(A, kA, curPar(m), m)
    endif

  enddo  ! m
endif


! if (A<0) then
 ! if (any(SibID(:,-A,kA)==91)) then
! !if (A==2969) then
  ! open (unit=69,file="log-CalcCPL.txt",status="unknown", position="append")
    ! write(69,*) ""
    ! write(69,*)  "A: ", A, kA, " (", curPar, ")   B: ", B, kB
    ! write(69,'("ParB: ", 2i6, " curGP-1: ", 2i6, " curGP-2: ", 2i6)')  ParB, curGP(:,1), curGP(:,2) 
    ! write(69,'("LL U: ", 4f8.1)') LLU
    ! do z=1,2
      ! write(69,'("LLCor: ", 3f8.1)') LLCP(:,z)
    ! enddo
    ! write(69,'("ALR: ", 3f8.1)') ALR
    ! do m=1,2
      ! write(69,*) "m=",m, ',a=1'
      ! do r=1,7
        ! write(69,'(7f8.1)') LLA(r,:,m,1)
      ! enddo
    ! enddo
    ! ! write(69,*) ""
    ! ! do m=1,2
      ! ! write(69,*) "m=",m, ',a=2'
      ! ! do r=1,7
        ! ! write(69,'(7f8.1)') LLA(r,:,m,2)
      ! ! enddo
    ! ! enddo
    ! write(69,*) ""
    ! write(69,'("LLtrio: ", 6f8.1)')  LLtrio
    ! write(69,'("ALRtrio-1: ", 6f8.1)')  ALRtrio(1,:)
    ! write(69,'("ALRtrio-2: ", 6f8.1)')  ALRtrio(2,:)
    ! write(69,'("LLX: ", 8f8.1)')  LLX(:,1), LLX(:,2)
  ! close(69)
! endif
! endif

end subroutine CalcCandParLL

! #####################################################################

subroutine CheckMaybeRev(A, kA, candP, kP, maybe, dLL)   ! check if A could be parent of candP instead.
use Global
implicit none

integer, intent(IN) :: A, kA, candP, kP
logical, intent(OUT) :: maybe
double precision, intent(OUT) :: dLL(2)
integer :: ParCP(2), ParA(2), fcl, n, notfcl(6), fclx
double precision :: ALR(2), LLrev(7,2), LLtmp(2), LLrevX(7,2)
logical :: SexUnk, ParOK                   

dLL = missing
maybe = .TRUE.

if (CandP == 0)  return

SexUnk = .FALSE.
if (A > 0) then
  if (Sex(A)>3)  SexUnk = .TRUE.
endif                              
ParCP = getPar(candP, kP)
ParA = getPar(A,kA)
if (ALL(ParCP > 0)) then
  if (hermaphrodites/=0 .and. candP>0 .and. A>0) then
    if (SelfedIndiv(CandP)) then
      maybe = .TRUE.
    else
      maybe = .FALSE.
    endif
  else
    maybe = .FALSE.
  endif
  if (.not. maybe)  return
else if (Complx==0 .and. candP < 0 .and. ParA(3-kP)/=0) then
  if (DumMate(-CandP, kP) == ParA(3-kP)) then
    maybe = .FALSE.
    return
  endif
else if (A>0 .and. .not. SexUnk) then
  if (ParCP(kA) > 0 .and. Sex(A)<3) then
    maybe = .FALSE.
    return
  endif
endif

ParOK = .TRUE.              
call ChkValidPar(candP, kP, A, kA, ParOK)
if (.not. ParOK) then
  maybe = .FALSE.
  return
endif

ALR = missing             
call CalcAgeLR(A, kA, candP, kP, 0,1, .TRUE., ALR(1))  
call CalcAgeLR(candP, kP, A, kA, 0,1, .TRUE., ALR(2))

if (ALR(2)==impossible .or. ALR(1)-ALR(2) > 2.0*ABS(TF)) then
  maybe = .FALSE.
  return
else if (ALL(ABS(ALR) < 0.01)) then   
  if (all(ParA==0) .and. all(ParCP==0))  return  ! No way to tell which is parent & which offspring
endif

LLrev = missing
if (hermaphrodites/=0 .and. ALL(ParCP > 0) .and. ParCP(1)==ParCP(2)) then
  call ChkValidPar(ParCP(1), 1, A, kA, ParOK)
  if (ParOK) then
    call CheckRel(ParCP(1), 1, A, kA, 1, LLrev(:,1), LLrev(:,2)) 
    notfcl = (/ (n, n = 2, 7) /)    
    do n=1,2
      if (LLrev(1,n) < 0) then
        maybe = .TRUE.
        dLL(n) = LLrev(1, n) - MaxLL(LLrev(notfcl, n))
      else
        maybe = .FALSE.
      endif
    enddo
  else
    maybe = .FALSE.
  endif
  return
endif

!if (A>0 .and. candP>0 .and. all(ParA >=0)) then
!  return
!endif

fclx = 7        
if (candP>0) then
  fcl = 1
  notfcl = (/ (n, n = 2, 7) /)
  if (A>0)  fclx = 1
else !if (candP<0) then
  fcl = 4
  notfcl = (/1,2,3,5,6,7/)
endif

LLrev = missing
LLtmp = missing               
call setParTmp(candP, kP, 0, kA)
call ChkValidPar(candP, kP, A, kA, ParOK)
if (ParOK) then
  call CheckRel(candP, kP, A, kA, fclx, LLrev(:,1), LLrev(:,2))
  if (ParCP(3-kA) < 0) then  ! include changes in CLL
    call CalcU(candP, kP, parCP(3-kA),3-kA, LLtmp(1))
    call setParTmp(candP, kP, A, kA)
    call CalcU(candP, kP, parCP(3-kA),3-kA, LLtmp(2))
    LLrev(fcl,1) = LLrev(fcl,1) + (LLtmp(2) - LLtmp(1))
  endif
endif
call setParTmp(candP, kP, parCP(kA), kA)

LLrevX = missing
if (SexUnk .and. ParOK) then  ! implies A>0; also consider as parent of other sex
  call ChkValidPar(candP, kP, A, 3-kA, ParOK)
  if (ParOK) then
    if (candP > 0) then
      call CheckRel(candP, 3-kA, A, 3-kA, fclx, LLrevX(:,1), LLrevX(:,2))
    else
      call CheckRel(candP, kP, A, 3-kA, fclx, LLrevX(:,1), LLrevX(:,2))
    endif
    if (ParCP(kA) < 0) then  ! include changes in CLL
      call CalcU(candP, kP, parCP(kA),kA, LLtmp(1))
      call setParTmp(candP, kP, A, 3-kA)
      call CalcU(candP, kP, parCP(kA),kA, LLtmp(2))
      LLrevX(fcl,1) = LLrevX(fcl,1) + (LLtmp(2) - LLtmp(1))
      call setParTmp(candP, kP, parCP(3-kA), 3-kA)
    endif
  endif
endif

do n=1,2
  if (LLrev(fcl,n) < 0) then
    maybe = .TRUE.
    dLL(n) = LLrev(fcl, n) - MaxLL(LLrev(notfcl, n))
    if (LLrevX(fcl,n) < 0) then
      dLL(n) = MAX(dLL(n), LLrevX(fcl, n) - MaxLL(LLrevX(notfcl, n)))
    endif
  else
    maybe = .FALSE.
  endif
enddo

end subroutine CheckMaybeRev

! #####################################################################

subroutine CheckMaybeDup(i,j, dLL)   ! check if A and B might be duplicate samples from same indiv
use Global
implicit none

integer, intent(IN) :: i,j
double precision, intent(OUT) :: dLL
integer :: l, CountMismatch, IsBothScored(-1:2,-1:2), IsDifferent(-1:2,-1:2), SnpdBoth_ij
double precision :: LLdup, LL(3)

dLL = -999D0            

IsBothScored = 1
IsBothScored(-1,:) = 0
IsBothScored(:,-1) = 0
IsDifferent = 0
IsDifferent(0, 1:2) = 1
IsDifferent(1, (/0,2/)) = 1
IsDifferent(2, 0:1) = 1

CountMismatch=0
SnpdBoth_ij = 0
do l=1, nSnp
  SnpdBoth_ij = SnpdBoth_ij + IsBothScored(Genos(l,i), Genos(l,j))
  CountMismatch = CountMismatch + IsDifferent(Genos(l,i), Genos(l,j))
  if (CountMismatch > MaxMismatchDup)  exit
enddo
if (CountMismatch > MaxMismatchDup)  return

LLdup = missing
LL = missing
call PairSelf(i, j, LLdup)
call PairPO(i,j, Sex(j),7, LL(1))
call PairFullsib(i,j, LL(2))
call CalcU(i,3,j,3,LL(3))

dLL = LLdup - MaxLL(LL)

end subroutine CheckMaybeDup

! #####################################################################

subroutine FindPairs(RoundC, PairID, PairType)
use Global
use qsort_c_module
implicit none

character(len=2), intent(IN) :: RoundC   ! for output file name
integer, intent(OUT) :: PairID(XP*nInd,2), PairType(XP*nInd)
logical :: UseAge, cPair, matpat(2)
integer :: k, i, j, top, PairTypeTmp(XP*nInd), PairIDtmp(XP*nInd,2), x, t
double precision :: dLL, PairLLRtmp(XP*nInd), LL(7), LLg(7), LRS(2), PairLLR(XP*nInd)
integer, allocatable, dimension(:) :: Rank
double precision, allocatable, dimension(:) :: SortBy
character(len=200) :: HeaderFMT, DataFMT

nPairs = 0
PairID = -9
PairLLR = missing
PairType = 0
UseAge = AgePhase > 0

do i=1,  nInd-1  
  if (MODULO(i,100)==0) call rchkusr()
  if (quiet==-1 .and. any(viginti==i)) call print_progress(i,t)
  if (ALL(Parent(i,:)/=0)) cycle
  do j=i+1,nInd
    if (skip(i) .and. skip(j))  cycle
    if (hermaphrodites==1 .and. ((ANY(parent(i,:)/=0) .and. ALL(parent(j,:)==0)) .or. &
     (ALL(parent(i,:)==0) .and. ANY(parent(j,:)/=0))))  cycle  
    LRS = 0D0
    matpat = .FALSE.
    do k=1,2
      if (Parent(i,k)/=0 .or. Parent(j,k)/=0) cycle
      if (Parent(i,k)==j .or. Parent(j,k)==i) cycle
      if (UseAge .and. getAP(AgeDiff(i,j), 3, 0, k, Impossible) == Impossible)  cycle
      matpat(k) = .TRUE.
    enddo
    if (DoMtDif) then
      if (mtDif(i,j))  matpat(1) = .FALSE.   ! not sharing same mt haplotype --> not mat sibs    
    endif
    if (Complx==0) then
      call PairQFS(i, j, LRS(2))  ! quick check
      if (LRS(2) < 2*TF) cycle  
    else
      call PairQHS(i, j, LRS(1)) 
      if (LRS(1) < 2*TF) cycle  
    endif
    if (Complx>0 .and. ((ALL(Parent(i,:)==0) .and. ALL(Parent(j,:)==0) .and. &
      UseAge .and. ALL(matpat)) .or. &
       (Hermaphrodites/=0 .and. (ALL(Parent(i,:)==0) .or. ALL(Parent(j,:)==0))))) then
      call PairQFS(i, j, LRS(2)) 
    endif
        
    cPair = .FALSE.    
    do k=1,2
      if ((.not. matpat(k)) .or. cPair)  cycle
      if (Complx==0 .and. k==2)  cycle
      if (k==2 .and. matpat(1)) then
        if (.not. UseAge) cycle
        if (LRS(2) < TF) cycle  
      endif    
      if (hermaphrodites==1 .and. (ALL(Parent(i,:)==0) .or. ALL(Parent(j,:)==0))) then
        if (LRS(2) < TF) cycle 
      endif  
      if (Complx==0) then
        x = 2
      else
        x = 3
      endif
      LLg = missing
      LL = missing
      if (AgeDiff(i,j)>=0) then
        call CheckPair(i, j, k, x, LLg, LL)
      else
        call CheckPair(j, i, k, x, LLg, LL)
      endif
      top = 0
      dLL = missing
      if (UseAge)  call BestRel(LL, 3, top, dLL)
      if (.not. UseAge)  call BestRel(LLg, 3, top, dLL)
      if (hermaphrodites==1 .and. (ALL(Parent(i,:)==0) .or. &
       ALL(Parent(j,:)==0)) .and. top/=2)  cycle 
      if (top==2 .or. top==3) then  
        if (nPairs >= XP*nInd) cycle  ! do in next round
        nPairs = nPairs+1
        PairID(nPairs, :) = (/ i, j /)
        PairLLR(nPairs) = dLL
        if (k==1 .and. matpat(2)) then
          pairType(nPairs) = 3
          cPair = .TRUE.
        else
          PairType(nPairs) = k
        endif
      endif
    enddo
  enddo
enddo

! sort by decreasing dLL
PairIDtmp = 0
PairLLRtmp = 0D0
allocate(Rank(nPairs))
allocate(SortBy(nPairs))
Rank = (/ (i, i=1, nPairs, 1) /)
SortBy = PairLLR(1:nPairs)
 
 call QsortC(SortBy, Rank(1:nPairs))
do i=1,nPairs
  PairTypeTmp(i) = PairType(Rank(nPairs-i+1))  ! decreasing order
  PairIDtmp(i,1:2) = PairID(Rank(nPairs-i+1), 1:2)  
  PairLLRtmp(i) = PairLLR(Rank(nPairs-i+1)) 
enddo 

PairType = PairTypeTmp
PairID = PairIDtmp 
PairLLR = PairLLRtmp
deallocate(Rank)
deallocate(SortBy)


write(HeaderFMT, '( "(2(a", I0, ", 4X), 3a8, a10)" )')  ID_len
write(DataFMT, '( "(2(a", I0, ", 4X), 3i8, f10.4)" )')  ID_len

open(unit=405, file="Pairs_"//RoundC//".txt", status="replace")
  write(405, HeaderFMT) "ID1", "ID2", "ID1_i", "ID2_i", "PatMat", "dLLR"
  do i=1,nPairs
    write(405, DataFMT) ID(PairID(i,1:2)), PairID(i,1:2), PairType(i), PairLLR(i)
  enddo
close(405)

end subroutine FindPairs

! #####################################################################

subroutine PairQFS(A, B, LR)  !quick check, not conditioning on parents.
use Global
implicit none

integer, intent(IN) :: A, B
double precision, intent(OUT) :: LR
integer :: l
double precision :: PrL(nSnp)

PrL = 0D0
do l=1,nSnp
  PrL(l) = LOG10(PFS(Genos(l,A), Genos(l,B), l))  
enddo
LR = SUM(PrL)

end subroutine PairQFS

! #####################################################################

subroutine PairQHS(A, B, LR)  !quick check, not conditioning on parents.
use Global
implicit none

integer, intent(IN) :: A, B
double precision, intent(OUT) :: LR
integer :: l
double precision :: PrL(nSnp)

PrL = 0D0
do l=1,nSnp
  PrL(l) = LOG10(PHS(Genos(l,A), Genos(l,B), l)) ! note: >0 for FS 
enddo
LR = SUM(PrL)

end subroutine PairQHS

! #####################################################################

subroutine IsSelfed(A, withFS, LR)  ! A is product of selfing, not conditioning on parents.
use Global
use CalcLik
implicit none

integer, intent(IN) :: A
logical, intent(IN) :: withFS
double precision, intent(OUT) :: LR
integer :: l, x,u, z, v
double precision :: PrX(3), PrXY(3,3,3), PrL(nSnp,4), LLtmp(4), PrZV(3,3), PrA(3,3)

PrL = 0D0
do l=1,nSnp
  do x=1,3
    if (withFS) then
      PrA = FSLik(l,A)
    else
      PrA = OKA2P(Genos(l,A), :, :)
    endif
    PrX(x) = PrA(x, x) * AHWE(x,l)  ! selfed
    PrXY(x,:,1) = PrA(x, :) * AHWE(x,l) * AHWE(:,l)   ! parents U
    PrXY(x,:,2) = PrA(x, :) * AKAP(x,:,l) * AHWE(:,l)   ! parents PO
    do z=1,3
      do v=1,3
        PrZV(z,v) = SUM(PrA(x, :) * AKA2P(x,z,v) * AKA2P(:,z,v) * &
          AHWE(z,l) * AHWE(v,l))    ! parents FS
      enddo
    enddo
    PrXY(x,:,3) = SUM(PrZV)
  enddo
  PrL(l,1) = LOG10(SUM(PrX))
  do u=1,3
    PrL(l,u+1) = LOG10(SUM(PrXY(:,:,u)))
  enddo
enddo
LLtmp = SUM(PrL,DIM=1)
LR = LLtmp(1) - MAXVAL(LLtmp(2:4))

end subroutine IsSelfed

! #####################################################################

subroutine CheckPair(A, B, kIN, focal, LLg, LL) 
! joined LL A,B under each hypothesis
use Global
implicit none

integer, intent(IN) :: A,B,kIN, focal
double precision, intent(OUT) :: Llg(7), LL(7)  ! PO,FS,HS,GG,FAU,HAU,U                                                                      
integer :: x, cgp(2), k, AB(2), i, z
double precision :: ALR(7), LLGR(3), ALRgr(3), LLGGP(5), ALRgg(5), &
  LLtmpAU(2,3), ALRAU(2,3), LLCC, LLFC, ALRp(2,2), LLX(4), LLPA(3), LLFAx(2), &
  LLZ(7), LLC(7,2), LLHH(2,2,3),  LLP(2,2), LLFA, LLPS(2,2), LLU, ALRf, ALRs, LLS!, ALRr
logical :: fclsib

LLg = missing
ALR = missing
LL = missing
  
if (kIN>2) then
  k = 1
else
  k = kIN
endif     
AB = (/A, B/)               

if (focal==1 .and. Sex(B)/=k .and. Sex(B)<3) then 
  LLg(1) = impossible
  LL(1) = impossible
  return
endif

call CalcU(A,k,B,k, LLg(7))
LL(7) = LLg(7)

if (hermaphrodites >0) then
  if (any(Parent(B,:) == A)) then
    LLg(2) = impossible
    if (focal==1 .or. focal==4) then  ! reverse (B-A) is possible
      LLg(1) = impossible
      LLg(4) = impossible
    endif
  else if (any(Parent(A,:) == B)) then
    LLg(2) = impossible
  endif
endif
LL = LLg
if (LLg(focal) == impossible)  return

if (focal==2 .or. focal==3) then
  fclsib = .TRUE.
else
  fclsib = .FALSE.
endif
 
do x=2,4   ! not x=1: done by ChkValidPar()
  if (focal == x) then
    call CalcAgeLR(A,Sex(A), B, Sex(B), k, x, .TRUE., ALR(x)) 
    if (ALR(x)==impossible .or. ALR(x)<5*TF) then
      LLg(x) = impossible    
    else
      if (fclsib) then 
        if (Parent(A, 3-k)/=0 .and. Parent(A,3-k)==Parent(B,3-k)) then
          call PairFullSib(A, B, LLg(2)) 
          LLg(3) = impossible
        else if (focal == 2) then
          call PairFullSib(A, B, LLg(2))
          if (Complx > 0)  call PairHalfSib(A, B, k, LLg(3))          
        else
          call PairHalfSib(A, B, k, LLg(3)) 
        endif
      endif
      if (focal==4)  call PairGP(A, B, k, focal, LLg(4))
    endif
  endif
enddo

if (LLg(focal)==impossible .and. .not. (focal==3 .and. LLg(2)<0D0)) then
  LL = LLg
  return
endif

if (focal>1 .and. focal <= 4 .and. LLg(focal) - LLg(7) < TA .and. &
   .not. (focal==3 .and. LLg(2)<0D0 .and. LLg(2) - LLg(7) > TA)) then
  LL = LLg
  return
endif

do x=1,4
  if (ALR(x) /= missing) cycle
  call CalcAgeLR(A,Sex(A), B, Sex(B), k, x, .TRUE., ALR(x)) 
  if (ALR(x) == impossible) then
    LLg(x) = impossible
    LL(x) = impossible
  endif
enddo

!~~~  parent  ~~~~~~~
if (LLg(1)==missing) then
  call PairPO(A, B, k, focal, LLg(1)) 
endif

!~~~  sibs  ~~~~~~~
if (LLg(2)==missing)  call PairFullSib(A, B, LLg(2)) 
if (LLg(3)==missing .and. Complx>0) then
  if (Parent(A, 3-k)/=0 .and. Parent(A,3-k)==Parent(B,3-k)) then
    LLg(3) = impossible
  else
    call PairHalfSib(A, B, k, LLg(3)) 
  endif
endif

do x=1,3
  LL(x) = addALR(LLg(x), ALR(x)) 
enddo

!~~~  GP  ~~~~~~~
if (LLg(4)==missing .and. ALR(4)/=impossible) then  ! GP?   
  call PairGP(A, B, k, focal, LLg(4))
endif
LL(4) = addALR(LLg(4), ALR(4))

!~~~  GGP  ~~~~~~~
LLGGP = missing
ALRgg = 0D0        
if (ALR(4)/=impossible) then   ! no ageprior for GGP
  if (AgeDiff(A,B)>=3) then
    call PairGGP(A, B, k, focal, LLGGP(1))
    call PairGGP(A,B,3-k, focal, LLGGP(2))    ! TODO drop?
  endif
  do x=1,3  ! hf
    call PairGA(A, B, k, x, LLGGP(2+x))   ! HS of GP 4th degree rel, but giving false pos. 
  enddo
endif
if (LLGGP(1) <0D0) then
  if (Parent(A,k)/=0) then
    call CalcAgeLR(Parent(A,k),k, B,Sex(B), k, 4, .TRUE., ALRgg(1))
    ALRgg(2) = ALRgg(1)
    call CalcAgeLR(Parent(A,k),k, B,Sex(B), k, 6, .TRUE., ALRgg(3))
    ALRgg(4) = ALRgg(3)
    call CalcAgeLR(Parent(A,k),k, B,Sex(B), k, 5, .TRUE., ALRgg(5))
  else
    ! TODO ageprior
  endif
endif

!~~~  FA/HA  ~~~~~~~
LLtmpAU = missing
ALRAU = missing
LLFAx = missing
! call CalcAgeLR(A,k, B,k, 0, 6, .TRUE., ALR(6))
do i=1,2   ! A, B
  do x=1,3  ! mat, pat, FS
    if (x<3 .and. Complx>0) then
      call CalcAgeLR(AB(i),k, AB(3-i),x, k, 6, .TRUE., ALRAU(i,x))
    else
      call CalcAgeLR(AB(i),k, AB(3-i),3, k, 5, .TRUE., ALRAU(i,x))
    endif
    if (ALRAU(i,x)/=impossible .and. .not. (Complx==0 .and. x<3)) then
      call PairUA(AB(i), AB(3-i),k, x, LLtmpAU(i,x))
      if (Complx/=1)  call FAx(AB(i), AB(3-i), LLFAx(i))  ! A inbred
    endif
  enddo
enddo

LLg(5) = MaxLL((/LLtmpAU(:,3), LLFAx/))
if (Complx>0)  LLg(6) = MaxLL(RESHAPE(LLtmpAU(:,1:2), (/2*2/) ))
do i=1,2   
  do x=1,3
    LLtmpAU(i,x) = addALR(LLtmpAU(i,x), ALRAU(i,x))
  enddo
  LLFAx(i) = addALR(LLFAx(i), ALRAU(i,3))
enddo
LL(5) =  MaxLL((/LLtmpAU(:,3), LLFAx/))
if (Complx>0)  LL(6) = MaxLL(RESHAPE(LLtmpAU(:,1:2), (/2*2/) ))
ALR(6) = MAXVAL((/ALRAU(1,1:2), ALRAU(2,1:2)/))

LLCC = missing
call PairCC(A, B, k, LLCC)   ! full 1st cousins
! call pairDHC(A,k,B, .FALSE., LLDHC)   ! TODO? test first                                                          

LLg(6) = MaxLL( (/LLg(6), LLGGP, LLCC/) )  ! most likely 3rd degree
do x=1,5
  LLGGP(x) = addALR(LLGGP(x), ALRgg(x))
enddo
LL(6) = MaxLL( (/LL(6), LLGGP, LLCC/) ) 

!~~~ FS + HC ~~~
LLFC = missing
if (Complx==2 .and. LLg(2)<0D0 .and. Parent(A,3-k)==Parent(B,3-k) .and. &
  (MaxLL(LLtmpAU(:,3)) - MaxLL(LLg(2:3)) > -TA)) then
  call FSHC(A, B, k, LLFC)
  if (LLFC > LLg(2) .and. LLFC<0D0) then
    LLg(2) = LLFC   ! note: no ageprior for FC
    LL(2) = addALR(LLg(2), ALR(2))
  endif
endif    

!~~~  Parent 3-k / B-A  ~~~~~~~  
LLP = missing
LLPS = missing
ALRp = missing 
LLU = missing 
if (all(Parent(A,:)/=B) .and. focal/=7 .and. .not. &
  (focal==1 .and. all(Parent(B,:)==0))) then   ! if no parents, LL indistinguishable.
  ! NOT: .and. all(Parent(A,:)==0) --> inconsistency when called by CalcCandPar
  do x=1,2
    do i=1,2
      if (i==1 .and. x==k .and. hermaphrodites==0)  cycle  ! 'default'  
      if (AgeDiff(AB(i),AB(3-i)) <=0 .or. Sex(AB(3-i))==3-x) then
        LLP(i,x) = impossible  ! known, invalid agedif / sex
        cycle
      endif 
      call CalcAgeLR(AB(i),Sex(AB(i)), AB(3-i),Sex(AB(3-i)), x, 1, .TRUE., ALRp(i,x))  
      if (ALRp(i,x) == impossible) then
        LLP(i,x) = impossible
        cycle
      endif
      if (Parent(AB(i),x) == 0) then
        call PairPO(AB(i), AB(3-i), x, 0, LLP(i,x))   ! WAS: FOCAL=0 ?  
        if (hermaphrodites/=0) then
          call PairPOX(AB(i), AB(3-i), x, focal, LLPS(i,x))    ! AB(3-i) result of selfing
        endif  
      else if (Parent(AB(i),x) < 0) then
        call AddParent(AB(3-i), -Parent(AB(i),x), x, LLP(i,x))
        call CalcU(AB(3-i), Sex(AB(3-i)), Parent(AB(i),x), x, LLU)
        if (LLP(i,x) < 0)  LLP(i,x) = LLP(i,x) - LLU + LLg(7)
      endif
    enddo
  enddo   
  
  if (focal==1) then
    z=6   ! obv. not 3rd degree rel, but most convenient place to store w/o having to call CheckPair again
  else
    z=1
  endif
  LLg(z) = MaxLL((/LLg(z), LLP(:,1), LLP(:,2), LLPS(:,1), LLPS(:,2)/))
  do x=1,2
    do i=1,2
      LLP(i,x) = addALR(LLP(i,x), ALRp(i,x))
      LLPS(i,x) = addALR(LLPS(i,x), ALRp(i,x))
    enddo
  enddo
  LL(z) = MaxLL((/LL(z), LLP(:,1), LLP(:,2), LLPS(:,1), LLPS(:,2)/))
endif

!~~~  Sibs 3-k  ~~~~~~~
LLS = missing
ALRs = missing
if ((ANY(Parent(A,:)/=0) .or. ANY(Parent(B,:)/=0)) &  ! else: could be both. 
  .and. .not. (Parent(A,k) == Parent(B,k))) then
  call CalcAgeLR(A,Sex(A), B, Sex(B), 3-k, 3, .TRUE., ALRs)
  if (ALRs /= impossible) then
    call PairHalfSib(A, B, 3-k, LLS)
    if (LLS<0D0 .and. (LLS - LLg(3) > TA .or. LLg(3)>0D0)) then
      if (focal == 3) then
        LLg(3) = MaybeOtherParent
        LL(3) = MaybeOtherParent
      else
        LLg(3) = LLS
        LL(3) = addALR(LLS, ALRs)
      endif
    endif
  endif
endif

!~~~  GP 3-k  ~~~~~~~
LLGR = missing 
ALRgr = missing    
if (((focal==3 .and. (MaxLL(LLg(2:3))>LLg(4) .or. LLg(4)>0D0)) .or. &
   (focal==4 .and. (LLg(4) > MaxLL(LLg(2:3))) .or. MaxLL(LLg(2:3))>0D0) .or. &
   (focal==1 .and. (LLg(1) > LLg(4) .or. LLg(4)>0))) .and. &
  Sex(B)<3 .and. .not. (Parent(A,3-k)==Parent(B,3-k) .and. Parent(A,3-k)/=0)) then
  cgp = getPar(Parent(A,3-k), 3-k)
  if (cgp(Sex(B)) == 0) then
    call CalcAgeLR(A,Sex(A),B,sex(B),3-k,4,.TRUE.,ALRgr(1))
    if (ALRgr(1)/=impossible) then
      call PairGP(A, B, 3-k, focal, LLGR(1))
    endif
    if (focal==3 .or. focal==1) then
      LLg(4) = MaxLL((/ LLg(4), LLGR(1) /)) 
      LL(4) = MaxLL((/LL(4), addALR(LLGR(1), ALRgr(1))/))
    else if (focal==4 .and. hermaphrodites/=2) then
      if (LLGR(1)<0D0 .and. (LLg(4) - LLGR(1)) < TA) then
        LLg(4) = MaybeOtherParent
        LL(4) = MaybeOtherParent
      endif
    endif
  endif
endif

!~~~  various double rel  ~~~~~~~
LLPA = missing
!call CalcAgeLR(B,Sex(B), A,Sex(A), 0, 1, .TRUE., ALRr)  
if ((Complx==2 .and. AgeDiff(A,B)>0 .and. ALR(1)/=impossible) .and. &  ! .and. (ALRr==impossible .or. ALRr < 3*TF))  
 .not. (focal==1 .and. Parent(A,3-k)==0 .and. Parent(B,3-k)/=0)) then  ! else messes up calccandpar combi's
! .and. LLg(1)<0 .and.  ((any(LLz < 0 .and. LLz > LLg(7)) .or. (any(LLHH < 0 .and. LLHH > LLg(7)))))
  do x=1,3
    if (ALRAU(1,x)/=impossible .and. LLg(1)/=impossible) then
      call PairPOHA(A, B, k, x, LLPA(x))
    endif
  enddo
  if (any(LLPA < 0) .and. MaxLL(LLPA) > LLg(1)) then 
    LLg(1) = MaxLL(LLPA)
    do x=1,3
      LLPA(x) = addALR(LLPA(x), ALR(1))
      LLPA(x) = addALR(LLPA(x), ALRAU(1, x))
    enddo
    LL(1) = MaxLL(LLPA)
  endif
endif

LLZ = missing
LLC = missing   
LLHH = missing 
LLX = missing             
if (Complx==2 .and. LL(2)<0 .and. .not. (SelfedIndiv(A) .or. SelfedIndiv(B))) then    
  if (LL(2) - LL(7) > TA) then  ! check if inbred FS (& boost LL) 
    call PairFSHA(A, B, k, LLZ(1))
    if (hermaphrodites/=2) then
      call PairFSHA(A, B, 3-k, LLZ(2))
    else 
      call PairFSSelfed(A, B, LLZ(2))
    endif
    if (MaxLL(LLZ(1:2)) > LLg(2) .and. ANY(LLZ(1:2)<0D0)) then
      LLg(2) = MaxLL(LLZ(1:2))
      LL(2) = addALR(LLg(2), ALR(2))
    endif
  endif
  if (MaxLL(LL(1:4)) - MaxLL(LL) > -TA) then
    if ((Parent(A,3-k)<=0 .or. Parent(B,3-k)<=0) .and. ALR(6)/=impossible .and. &
      .not. ALR(6) < -HUGE(0.0D0)) then
      do i=1,2
        if (focal==1 .and. Parent(A,i)==0 .and. Parent(B,i)/=0)  cycle   ! else indirect assignment of parent(B,i) to A
        do x=1,3  
          call PairHSHA(A, B, i, x, LLHH(i, 1,x), .FALSE.)
          call PairHSHA(B, A, i, x, LLHH(i, 2,x), .FALSE.)
        enddo
      enddo
      if (ANY(LLHH<0D0)) then
        ! if (LLg(2) - MaxLL(RESHAPE(LLHH, (/2*2*3/))) < TA .and. fclsib .and. hermaphrodites/=2) then
          ! LL(2) = MaybeOtherParent
        ! endif
        if (MaxLL((/LLHH(k,1,:), LLHH(k,2,:)/)) > LLg(3) .and. &
          MaxLL((/LLHH(k,1,:), LLHH(k,2,:)/)) < 0D0) then  ! .and. fclsib ?
          LLg(3) = MaxLL((/LLHH(k,1,:), LLHH(k,2,:)/))
          LL(3) = addALR( addALR(LLg(3), ALR(3)), ALR(6))   ! TODO: fix ageprior
        endif 
        if (MaxLL((/LLHH(3-k,1,:), LLHH(3-k,2,:)/)) > LLg(6) .and. focal/=2) then
          LLg(6) = MaxLL((/LLHH(3-k,1,:), LLHH(3-k,2,:)/))
          LL(6) = addALR( addALR(LLg(6), ALR(3)), ALR(6))
        endif
      endif
      if (LL(2)/=MaybeOtherParent .and. fclsib .and. hermaphrodites/=2 .and. &
       LLg(1)/=impossible) then 
         if (AgeDiff(A,B)>0 .and. Sex(B)/=k) call PairHSPO(A,B,LLX(1))  
         if (AgeDiff(B,A)>0 .and. Sex(A)/=k) call PairHSPO(B,A,LLX(2))
        if (any(LLX(1:2)<0D0) .and. (LLg(2) - MaxLL(LLX(1:2))) < TA) then 
          LL(2) = MaybeOtherParent
        endif
      endif
    endif
    if (LLg(3)/=impossible .and. LLg(4)/=impossible .and. focal/=1 .and. &
     (Parent(A,3-k)==0 .or. Parent(B,3-k)==0)) then               
      call PairHSGP(A, B,k, LLX(3))
      call PairHSGP(B, A,k, LLX(4))
      if (any(LLX(3:4)<0D0)) then
        if ((LLg(2) -MaxLL(LLX(3:4)))<TA  .and. fclsib .and. hermaphrodites/=2) then
          LL(2) = MaybeOtherParent
        endif
        do x=3,4
          if (MaxLL(LLX(3:4)) > LLg(x) .and. MaxLL(LLX(3:4)) < 0D0) then
            LLg(x) = MaxLL(LLX(3:4))
            LL(x) = addALR(LLg(x), ALR(x))
          endif
        enddo
      endif
    endif
    if (ANY(LL(2:3)/=MaybeOtherParent) .and. fclsib .and. hermaphrodites/=2) then
      if (ANY(Parent(AB,3-k)==0) .and. ANY(Parent(AB,3-k)<0)) then   ! check if add to opp. sibship
        if (Parent(A,3-k) < 0 .and. Parent(B,3-k)==0) then  
          call CheckAdd(B, -Parent(A,3-k), 3-k, 3, LLC(:,1), LLC(:,2))  !! DANGER
        else if (Parent(A,3-k)==0 .and. Parent(B,3-k)<0) then
          call CheckAdd(A, -Parent(B,3-k), 3-k, 3, LLC(:,1), LLC(:,2))  !! DANGER 
        endif
        if (LLC(2,1)<0D0 .and. (LLC(2,1) - MaxLL(LLC((/1,3,4,5,6,7/),1))) < TA) then     ! TODO: use as check for GP ?
          LL(2) = MaybeOtherParent
        endif
        if (LLC(3,1)<0D0 .and. (MaxLL(LLC((/1,2,4,5,6,7/),1)) - LLC(3,1)) > TA) then  
          LL(3) = MaybeOtherParent
        endif
      endif
    endif  
    if ((LL(2)/=MaybeOtherParent .and. ALL(Parent(A,:)==0) .and. ALL(Parent(B,:)==0)) .or. &
      (focal==1 .and. LLg(1) - MaxLL(LLg) > -TA .and. ANY(ALRAU /= impossible)))  then
      call pairFAHA(A, B, .FALSE., LLZ(5))
      call pairFAHA(B, A, .FALSE., LLZ(6))
      if (ANY(LLZ(5:6) < 0D0)) then
        if ((LLg(2) - MaxLL(LLZ(5:6))) < TA  .and. fclsib) then
          LL(2) = MaybeOtherParent
        endif
        if (MaxLL(LLZ(5:6)) > LLg(5)) then
          LLg(5) = MaxLL(LLZ(5:6))
          LL(5) = addALR(LLg(5), ALR(6))
        endif
      endif
    endif
    if (LL(2)/=MaybeOtherParent .and. fclsib) then  ! check if GG in any way. can't be FS and GP
      if(LLGR(1)==missing)  call PairGP(A, B, 3-k, focal, LLGR(1))
      if (AgeDiff(A,B)==missing) then
        call PairGP(B, A, k, focal, LLGR(2))
        call PairGP(B, A, 3-k, focal, LLGR(3))
      endif
      if (MaxLL(LLGR)<0D0 .and. (LLg(4) - MaxLL(LLGR)) <TA) then
        LLg(4) = MaxLL(LLGR)
        if (ALRgr(1)==missing) then
          call CalcAgeLR(A,Sex(A),B,sex(B),3-k,4,.TRUE.,ALRgr(1))
          LLGR(1) = addALR(LLGR(1), ALRgr(1))
        endif
        do x=1,2
          call CalcAgeLR(B,Sex(B),A,Sex(A),x,4,.TRUE.,ALRgr(x+1))
           LLGR(x+1) = addALR(LLGR(x+1), ALRgr(x+1))
        enddo
        LL(4) = MaxLL(LLGR) 
      endif
    endif
  endif
  if (LL(2)==MaybeOtherParent)  LLg(2) = MaybeOtherParent
endif

if (Complx==2 .and. (MaxLL(LL)>=LL(3) .or. MaxLL(LL)==LL(2))) then  
  call pairHSHAI(A, B, k, LLZ(3)) ! HS + inbr HA
  call pairHSHAI(B, A, k, LLZ(4))
  if (MaxLL(LLZ(3:4)) < 0D0 .and. MaxLL(LLZ(3:4)) > LLg(3)) then
    LLg(3) = MaxLL(LLZ(3:4))
    LL(3) = addALR(MaxLL(LLZ(3:4)), ALR(3))
  endif
  
  call PairHSCC(A,B, LLZ(7))
  if (LLZ(7) < 0 .and. LLZ(7) > LLg(3)) then
    LLg(3) = LLZ(7)
    LL(3) = addALR(LLZ(7), ALR(3))
  endif  
endif

LLFA = missing
ALRf = missing            
if (focal==1) then  ! check if FS of other parent   !  .and. LL(5) < LL(7)
  call CalcAgeLR(A,3-k, B,0, 3-k, 5, .TRUE., ALRf)
  if (ALRf /= impossible) then  
    call PairUA(A, B, 3-k, 3, LLFA)
  endif
  if (LLFA<0D0 .and. ((LL(5) - addALR(LLFA, ALRf)) < TA .or. LL(5)>0)) then
    LLg(5) = LLFA
    LL(5) = addALR(LLFA, ALRf)
  endif
endif


! if ( (A==1840 .and. B==1841) .or. (A==1841 .and. B==1840) ) then    
  ! open (unit=42,file="log.txt",status="unknown", position="append")
    ! write (42, *) ""
    ! write (42, '("pair?", 2i6, "; parents: ", 2i6, ", ", 2i6)') A, B, Parent(A,:), Parent(B,:)
    ! write (42, '("LL  ", 7f9.2)') LL
    ! write (42, '("LLG ", 7f9.2, "  ", 2i3)') LLg, k, focal
    ! write (42, '("ALR ", 7f9.2, "  ", i3)') ALR
    ! write (42, '("LLHH_k ", 6f9.2)')  LLHH(k,1,:), LLHH(k,2,:)
    ! write (42, '("LLHH_!k ", 6f9.2)')  LLHH(3-k,1,:), LLHH(3-k,2,:)
    ! write (42, '("LLX ", 5f8.1)') LLX
    ! write (42, '("LLZ ", 7f8.1)') LLZ
    ! write (42, '("LL PO-HA ", 3f8.1, "  FAx: ", 2f8.1)') LLPA, LLFAx
    ! write (42, '("LLUA ", 3f8.1, "; ", 3f8.1, "; ", 5f8.1)') LLtmpAU(1,:), LLtmpAU(2,:), LLCC
    ! write (42, '("ALR-AU ", 3f8.1, "; ", 3f8.1)') ALRAU(1,:), ALRAU(2,:)
    ! write (42, '("LLGGP ", 5f8.1)') LLGGP 
    ! write (42, '("LL sib ", f8.1, ",", f8.1)') LLS, ALRs
    ! write (42, '("LLGR ", 3f8.1, "; ", 3f8.1)') LLGR, ALRgr
    ! write (42, '("LLC ", 7f8.1, ", LLFC: ", f8.1)') LLC(:,1), LLFC
    ! write (42, '("LLFA ", 2f8.1)') LLFA, ALRf
    ! write (42, '("LLP ", 4f8.1, ", ", 4f8.1)') LLP(1,:), LLP(2,:), ALRp(1,:), ALRP(2,:) ! LLPS(1,:), LLPS(2,:)
    ! write (42, '("LU ", f8.1, "; ", 3f8.1)') LLg(7), Lind(A) + Lind(B), Lind(A), Lind(B) 
    ! write (42, *) "" 
  ! close(42) 
! endif

end subroutine CheckPair

! #####################################################################

subroutine PairSelf(A, B, LL)  ! A==B; currently only called w/o parents
use Global
implicit none

integer, intent(IN) :: A, B
double precision, intent(OUT) :: LL
integer :: l
double precision :: PrL(nSnp)

forall (l=1:nSnp)  PrL(l) = LOG10(SUM( LindX(:,l,A) * OcA(:,Genos(l,B)) ))
LL = SUM(PrL)

end subroutine PairSelf

! #####################################################################

subroutine GetPOconfigs(A, B, k, focal, Maybe)  ! TODO: check for PairPOX
use Global
implicit none

integer, intent(IN) :: A, B, k, focal
logical, intent(OUT) :: Maybe(5)
integer :: GA(2), PAB, AncB(2,mxA)

Maybe = .FALSE.  ! 1: non-inbred, 2: B PO & GP; 3: B PO & HS, 
                 ! 4: Parent(A,3-k) ancestor of B, 5: B selfing
! B HS of Parent(A,3-k) : see PairPOHA
Maybe(1) = .TRUE.  
GA = getPar(Parent(A,3-k), 3-k)   
if (GA(k)==B .or. (Complx==2 .and. GA(k) == 0)) then
  Maybe(2) = .TRUE.
endif

if ((Complx==2 .and. (Parent(A,3-k)==Parent(B,3-k) .or. Parent(A,3-k)==0 &
  .or. Parent(B,3-k)==0)) .or. (Parent(A,3-k)==Parent(B,3-k) .and. focal==7)) then
    if ((focal==1 .or. focal==7) .and. Parent(A,3-k)==0 .and. Parent(B,3-k)/=0) then
      Maybe(3) = .FALSE.
    else if (Parent(A,3-k)/=0) then
      if (ANY(GA == B)) then
      Maybe(3) = .FALSE.
    else
      Maybe(3) = .TRUE. 
    endif
  else
    Maybe(3) = .TRUE.
  endif
endif

if (Parent(A,3-k)==Parent(B,3-k) .and. Parent(A,3-k)/=0) then
  Maybe(1) = .FALSE.  ! becomes config (3)
else if (Parent(A,3-k) > 0) then
  if (any(GA == B)) then
    Maybe(1) = .FALSE.  ! becomes config (2).
  endif
else if (Parent(A,3-k) < 0) then
  if (any (GpID(:, -Parent(A,3-k), 3-k) == B)) then
    Maybe(1) = .FALSE.
  endif
endif

AncB = 0
call getAncest(B, k, AncB)   
if (ANY(AncB(3-k, 3:4) == Parent(A,3-k))) then
  Maybe(4) = .TRUE.  ! calc at end
endif

if (hermaphrodites>0) then
  if (Parent(A,3-k) == B) then
    Maybe(1:4) = .FALSE.
    Maybe(5) = .TRUE.
  else if (focal/=1 .and. focal/=7 .and. Parent(A,3-k)<=0 .and. Sex(B)==4) then
    Maybe(5) = .TRUE.  
  else
    Maybe(5) = .FALSE.     ! else single/double parent same LL   
  endif
  if (SelfedIndiv(B)) then
    Maybe(3) = .FALSE.
  endif
endif

PAB = Parent(A,3-k)
if(Maybe(3) .and. Parent(A,3-k)==0) then  ! B PO & HS
  PAB = Parent(B,3-k)
  if (PAB < 0) then
    if (any(parent(SibID(1:ns(-PAB,3-k),-PAB,3-k),k) == A)) then
      Maybe(3) = .FALSE.  ! not implemented
    endif
  endif
endif       

end subroutine GetPOconfigs

! #####################################################################

subroutine PairPO(A, B, k, focal, LL)
use Global
use CalcLik
implicit none

integer, intent(IN) :: A, B, k, focal
double precision, intent(OUT) :: LL
integer :: l, x, y,m, curPar(2), AncB(2,mxA), PAB, GA(2), Gtmp(2)
double precision :: PrL(nSnp,5), PrX(3,3,5), PrPA(3), PrB(3),PrPB(3,2),&
  LLtmp(5), PrPAB(3), PrPAX(3), PrG(3), LLX(4)
logical :: Maybe(5), ParOK                

LL = missing
if (Sex(B)<3 .and. Sex(B)/=k) then
  LL = impossible
else if(Parent(A,k)>0) then  ! allow dummy to be replaced (need for AddFS)
  if (Parent(A,k)==B) then
    LL = AlreadyAss
  else if (focal==1) then    ! else do consider (in case current parent wrong)
    LL = impossible  
  endif
else if (Parent(A,k)<0) then
  if (any(SibID(:,-parent(A,k),k) == B)) then
    LL = impossible
  endif
endif
if (LL/=missing) return

ParOK = .TRUE.
call ChkValidPar(A, Sex(A), B, k, ParOK)
if (.not. ParOK .and. focal/=7) then
  LL = impossible
  return
endif

Maybe = .TRUE.
call GetPOconfigs(A, B, k, focal, Maybe)
! 1: non-inbred, 2: B PO & GP; 3: B PO & HS, 
! 4: Parent(A,3-k) ancestor of B, 5: B selfing
GA = getPar(Parent(A,3-k), 3-k)  
PAB = Parent(A,3-k)
if(Maybe(3) .and. Parent(A,3-k)==0) then  ! B PO & HS
  PAB = Parent(B,3-k)
endif    

PrL = 0D0
LLtmp = missing                             
do l=1,nSnp    
  call ParProb(l, Parent(A,3-k), 3-k, A, 0, PrPA)
  PrB = LindX(:,l,B)                  
  if (Maybe(2)) then
    call ParProb(l, Parent(A,3-k), 3-k, A, -4, PrPAX)
    if (GA(3-k)==B) then  ! hermaphrodite  
      m = k
    else
      m = 3-k
    endif
    if (Parent(A,3-k)>0) then
      call ParProb(l, GA(m), m, Parent(A,3-k), 0, PrG)  
    else 
      call ParProb(l, GA(m), m, 0, 0, PrG)       
    endif
  endif
  if (Maybe(3)) then
    do m=1,2  
      call ParProb(l, Parent(B,m), m, B, 0, PrPB(:,m))
    enddo
    call ParProb(l, PAB, 3-k, A, B, PrPAB)
  endif

  do x=1,3  ! B
    do y=1,3  ! parent(A,3-k)
      PrX(x,y,:) = OKA2P(Genos(l, A), x, y)
      if (Maybe(1))  PrX(x,y,1) = PrX(x,y,1) * PrB(x) * PrPA(y)
      if (Maybe(2)) then  ! B PO & GP
        PrX(x,y,2) = PrX(x,y,2) * PrB(x)* PrPAX(y) * SUM(AKA2P(y,x,:) * PrG)
      endif
      if (Maybe(3)) then  ! B PO & HS
        PrX(x,y,3) = PrX(x,y,3) * PrPAB(y) * SUM(AKA2P(x,y,:) * PrPB(:,k)) * OcA(x,Genos(l,B))
      endif
      if (Maybe(5)) then  ! B selfing
        if (x/=y)  PrX(x,y,5) = 0D0
        if (x==y)  PrX(x,y,5) = PrX(x,y,5) * PrB(x)
      endif
    enddo
  enddo
  do x=1,5
    if (Maybe(x))  PrL(l,x) = LOG10(SUM(PrX(:,:,x)))
  enddo 
enddo
LLtmp = SUM(PrL, DIM=1)
if (Parent(A,3-k) > 0 .and. Maybe(2)) then
  LLtmp(2) = LLtmp(2) - Lind(Parent(A,3-k))
endif
do x=1,5
  if (.not. Maybe(x) .or. LLtmp(x)>=0)  LLtmp(x) = impossible
enddo

if (Parent(A,3-k)<0 .and. ANY(Maybe(2:4)) .and. ParOK) then 
  curPar = Parent(A, :)
  LLX = missing                     
  call setParTmp(A, Sex(A), 0, k)     ! B vs none.    
  call CalcU(A,3-k, B,k, LLX(1))
  call CalcU(Parent(A,3-k),3-k, B,k, LLX(2))

  call setParTmp(A, Sex(A), B, k)  
  call CalcU(Parent(A,3-k),3-k, B,k, LLX(3))
  LLtmp(1) = LLX(1) + (LLX(3) - LLX(2))
  
  if (Maybe(2) .and. Parent(A,3-k)<0 .and. .not. ANY(GA == B)) then    ! B PO & GP ?
    AncB = 0                                                                                
    call GetAncest(B, Sex(B), AncB)
    call ChkValidPar(Parent(A,3-k), 3-k, B, k, ParOK)
    if (ParOK .and. .not. ANY(AncB(3-k,:)==Parent(A,3-k))) then
      Gtmp = getPar(Parent(A,3-k), 3-k)
      call setParTmp(Parent(A,3-k), 3-k, B, k)
      call CalcU(Parent(A,3-k),3-k, B,k, LLX(4))
      LLtmp(2) = LLX(1) + (LLX(4) - LLX(2))
      call setParTmp(Parent(A,3-k), 3-k, Gtmp(k), k)
    endif
  endif
  
  call setParTmp(A, Sex(A), curPar(k), k)  ! restore
endif                          

LL = MaxLL(LLtmp)

! if (A==518 .and. B==1580) then
   ! open (unit=42,file="log.txt",status="unknown", position="append")
  ! write (42, '("PairPO : ", 3i6, " + ", i3, 3i6, 5l3, 5f8.1, i3, i6)')  A, Parent(A,:), &
    ! k, B, Parent(B,:), Maybe, LLtmp, focal, PAB
    ! close(42)
! endif

end subroutine PairPO

! #####################################################################

subroutine PairPOX(A, B, k, focal, LL)    ! B parent of A; B result of selfing
use Global
implicit none

integer, intent(IN) :: A, B, k, focal
double precision, intent(OUT) :: LL
integer :: l, x, y, z, GA(2), PAB, m
double precision :: PrL(nSnp,5), PrPB(3), PrPA(3), PrXYZ(3,3,3,5), LLtmp(5), &
  PrG(3), PrPAX(3)
logical :: ParOK, Maybe(5)

LL = missing
if (hermaphrodites==0) then
  LL = impossible
else if (Parent(A,k)/=0) then
  LL = impossible
else if (Parent(B,1) >0) then
  if (Parent(B,1)/=Parent(B,2) .and. Parent(B,2)/=0) then
    LL = impossible
  else if (Sex(Parent(B,1))/=4) then
    LL = impossible
  endif
endif
if (LL == impossible)  return

ParOK = .TRUE.                             
call ChkValidPar(A, Sex(A), B, k, ParOK)
if (.not. ParOK .and. focal/=7) then
  LL = impossible
  return
endif

Maybe = .TRUE.     
call GetPOconfigs(A, B, k, focal, Maybe)
Maybe(4) = .FALSE.   ! not implemented with B selfed
GA = getPar(Parent(A,3-k), 3-k)  
PAB = Parent(A,3-k)
if(Maybe(3) .and. Parent(A,3-k)==0) then  ! B PO & HS
  PAB = Parent(B,3-k)
endif

PrL = 0D0 
PrPAX = missing          
do l=1,nSnp
  call ParProb(l, Parent(B,1), 1, B, 0, PrPB)
  call ParProb(l, Parent(A,3-k), 3-k, A, 0, PrPA)
  if (Maybe(2)) then
    call ParProb(l, Parent(A,3-k), 3-k, A, -4, PrPAX)
    if (GA(3-k)==B) then  ! hermaphrodite    
      m = k
    else
      m = 3-k
    endif
    if (Parent(A,3-k)>0) then
      call ParProb(l, GA(m), m, Parent(A,3-k), 0, PrG)  
    else 
      call ParProb(l, GA(m), m, 0, 0, PrG)    
    endif
  endif
  
  do x=1,3  ! B
    do y=1,3  ! parent A 3-k
      do z=1,3  ! parent B
        PrXYZ(x,y,z,:) = OKA2P(Genos(l, A), x, y) * OcA(x,Genos(l,B)) * AKA2P(x,z,z) * PrPB(z)
        if (Maybe(1))  PrXYZ(x,y,z,1) = PrXYZ(x,y,z,1) * PrPA(y)
        if (Maybe(2))  PrXYZ(x,y,z,2) = PrXYZ(x,y,z,2) * SUM(PrPAX(y) * AKA2P(y,x,:) * PrG)
        if (Maybe(3) .and. y/=z)  PrXYZ(x,y,z,3) = 0D0   ! else as is
        if (Maybe(5) .and. x/=y)  PrXYZ(x,y,z,5) = 0D0   ! else as is
      enddo
    enddo
  enddo
  do x=1,5
    if (Maybe(x))  PrL(l,x) = LOG10(SUM(PrXYZ(:,:,:,x)))
  enddo 
enddo
LLtmp = SUM(PrL, DIM=1)
if (Parent(A,3-k) > 0 .and. Maybe(2)) then
  LLtmp(2) = LLtmp(2) - Lind(Parent(A,3-k))
endif
do x=1,5
  if (.not. Maybe(x) .or. LLtmp(x)>=0 .or. LLtmp(x) < -HUGE(0D0))  LLtmp(x) = impossible
enddo

LL = MaxLL(LLtmp)

! if (A==156 .and. B==106 .and. k==1) then
  ! write (*, '("PairPOX: ", 3i6, " + ", 2i6, 5l3, 5f8.1, i3)')  A, Parent(A,:), &
    ! k, B, Maybe, LLtmp, focal
! endif

end subroutine PairPOX

! #####################################################################

subroutine PairFullSib(A, B, LL)
use Global
use CalcLik
implicit none

integer, intent(IN) :: A, B
double precision, intent(OUT) :: LL
integer :: x, y, l,k, Par(2), ix, i, AlreadyHS, Ei, DoQuick, AB(2), ParA(2)
double precision :: PrL(nSnp), PrXY(3,3), Px(3,2), LUX(2), LLtmp, &
  dx(maxSibSize), LRQ(2), PrE(3), PrXYb(3,3,2), PrFS(3,3)
logical :: AncOK(2)

LL = missing
Par = 0  ! joined parents of A & B
if (Parent(A,1)==Parent(B,1) .and. Parent(A,1)/=0 .and. &
  Parent(A,2)==Parent(B,2) .and. Parent(A,2)/=0) then ! already FS
  LL = AlreadyAss
  return
else 
  do k=1,2
    if (Parent(A,k) == B .or. Parent(B,k) == A) then
      LL = impossible
      return
    else if (Parent(A,k)/=Parent(B,k) .and. .not. (Parent(A,k)==0 .or. &
      Parent(B,k)==0)) then
      LL = impossible
      return
    else if (Parent(A,k)/=0) then
      Par(k) = Parent(A,k)
    else
      Par(k) = Parent(B,k)
    endif        
  enddo
endif  

LRQ = missing
call CalcP2(A, Sex(A), Par(1), Par(2), 1, LRQ(1))
call CalcP2(B, Sex(B), Par(1), Par(2), 1, LRQ(2))
if (any(LRQ == impossible)) then
  LL = impossible
  return
endif

AncOK = .TRUE.
call ChkAncest(A,0,B,0, AncOK(1))
call ChkAncest(B,0,A,0, AncOK(2))
if (any(.not. AncOK)) then  ! cannot be both ancestors of eachother & full siblings
  LL = impossible
  return
endif

AlreadyHS = 0
do k=1,2
  if (Parent(A,k) == Parent(B,k) .and. Parent(A,k)<0) then
    AlreadyHS = k
  endif
enddo

PrL = 0D0 
LUX = 0D0
LLtmp = missing                    
dx = missing
ParA=0              
if ((Par(1) < 0 .or. Par(2)<0) .and. (AlreadyHS==0 .or. (Par(1)/=0 .and. Par(2)/=0))) then
  AB = (/A, B/)
  do k=1,2 
    if (Par(k) >= 0)  cycle
    if (Parent(A,3-k)/=0) then  ! set parA=0 to be consistent with addFS
      ParA = getPar(A,1)
      call setParTmp(A,1,0,3-k)
    endif
    call CalcU(A, 1, B, 1, LUX(1))
    if (ParA(3-k)/=0)  call setParTmp(A,1,ParA(3-k),3-k) 
    do x=1,2
      if (Parent(AB(x),k)==Par(k) .and. Parent(AB(3-x),k)==0) then
        call addFS(AB(3-x), -Par(k), k, 0, k, LLtmp, ix, dx)   ! call AddFS
        do i=1, nS(-Par(k),k)
          if (SibID(i,-Par(k),k) == AB(x)) then
            if (dx(i) < impossible) then
              LL = LUX(1) + dx(i)
            else
              LL = dx(i)
            endif
            return
          endif
        enddo
      endif
    enddo
  enddo
  
else if (AlreadyHS/=0 .and. (Par(1) < 0 .or. Par(2)<0)) then
  k = AlreadyHS
  DoQuick = 1             
  call ChkDoQuick(-Parent(A,k),k,DoQuick) 
  if (DoQuick>1 .or. DoQuick==-1 .or. DoQuick==-3) then
    LL = NotImplemented
    Return
  endif

  do l=1, nSnp  
    call ParProb(l, Par(k), k, -1, 0, Px(:,k))
    call ParProb(l, Par(3-k), 3-k, 0, 0, Px(:,3-k))   ! other parent ==0
    do x=1,3  ! already shared parent
      do y=1,3  ! other parent
        PrXYb(x,y,:) = Px(x,k) * Px(y,3-k)
        do i=1, ns(-Par(k),k)
          Ei = SibID(i,-Par(k), k)
          if (nFS(Ei)==0 .or. Ei==A .or. Ei==B)  cycle
          call ParProb(l, Parent(Ei, 3-k), 3-k, Ei,-1, PrE)
          PrFS = FSLik(l,Ei)
          PrE = PrE * PrFS(:,x)
          if (.not. ALL(PrE==1D0))  PrXYb(x,y,:) = PrXYb(x,y,:) * SUM(PrE)
        enddo
        PrXYb(x,y,2) = PrXYb(x,y,2) * OKA2P(Genos(l,A), x, y) * OKA2P(Genos(l,B), x, y)
      enddo
    enddo
    PrL(l) = LOG10(SUM(PrXYb(:,:,2))) - LOG10(SUM(PrXYb(:,:,1)))
  enddo
  LL = SUM(PrL) 

else  
  do l=1, nSnp
    do k=1,2
      if (Parent(A,k)==Parent(B,k)) then   
        call ParProb(l, Par(k), k, A, B, Px(:,k))
      else if (Parent(A,k)==Par(k)) then
        call ParProb(l, Par(k), k, A, 0, Px(:,k))
      else if (Parent(B,k)==Par(k)) then
        call ParProb(l, Par(k), k, B, 0, Px(:,k))
      else
        call ParProb(l, Par(k), k, 0, 0, Px(:,k))
      endif       
    enddo 
  
    do x=1,3
      do y=1,3
        PrXY(x,y) = Px(x,1) * Px(y,2) * OKA2P(Genos(l,A), x, y) * OKA2P(Genos(l,B), x, y)
      enddo
    enddo
    PrL(l) = LOG10(SUM(PrXY))
  enddo
  LL = SUM(PrL) 
endif

!if (A==77 .and. B==59) print *, "PairFS: ", A, B, Par, ", ", AlreadyHS, LL

end subroutine PairFullSib

! #####################################################################

subroutine PairHalfSib(A, B, k, LL)
use Global
implicit none

integer, intent(IN) :: A, B, k
double precision, intent(OUT) :: LL
integer :: x,y, l, Par, Inbr, AB(2), i, PP
double precision :: PrL(nSnp), PrX(3), PrPx(3,2), PrXY(3,3), LLtmp(2)
logical :: AncOK(2)

LL = missing
Par = 0  ! parent K
if (Parent(A,k)/=0) then
  if (Parent(A,k)/=Parent(B,k) .and. Parent(B,k)/=0) then
    LL = impossible ! mismatch
  else if (Parent(A,k)==Parent(B,k)) then
    LL = AlreadyAss ! already captured under H0
  else
    Par = Parent(A,k)
    if (Par>0) then
      if (AgeDiff(B, Par) <= 0) then  ! Par(k) younger than B
        LL = impossible
      endif
    endif
  endif
else if (Parent(B,k)/=0) then
  Par = Parent(B,k)
  if (Par>0) then
    if (AgeDiff(A, Par) <= 0) then  ! Par(k) younger than A
      LL = impossible
    endif
  endif
endif
if (LL/=missing) return

AncOK = .TRUE.                       
call ChkAncest(Parent(A,k),k,B,0, AncOK(1))
call ChkAncest(Parent(B,k),k,A,0, AncOK(2))
if (any(.not. AncOK)) then
  LL = impossible
  return
endif

AB = (/ A, B /)
LLtmp = missing
if (Par < 0 .and. (all(Parent(A,:)>=0) .or. all(Parent(B,:)>=0))) then
  do i=1,2
    if (Parent(AB(i),k) == Par) then
      call CalcU(AB(3-i), 3, Par, k, LLtmp(1))
      call addSib(AB(3-i), -Par, k, LLtmp(2))
      call CalcU(AB(1),3, AB(2),3, LL)
      LL = LL + (LLtmp(2) - LLtmp(1))
      return
    endif
  enddo

else if (Par < 0 .and. any(Parent(AB,3-k) < 0)) then
  PP = 0
  do i=1,2
    if (Parent(AB(i),3-k) < 0)  PP = Parent(AB(i),3-k)
  enddo
  do i=1,2
    if (Parent(AB(i),k) == Par .and. parent(AB(3-i),k)==0) then
      call CalcU(Par, k, PP, 3-k, LLtmp(1))
      call SetParTmp(AB(3-i),3, Par, k)
      call CalcU(Par, k, PP, 3-k, LLtmp(2))
      call SetParTmp(AB(3-i),3, 0, k)
      call CalcCLL(-Par,k)
      call CalcCLL(-PP, 3-k)
      call CalcU(AB(1),3, AB(2),3, LL)
      LL = LL + (LLtmp(2) - LLtmp(1))
      return
    endif
  enddo
endif

Inbr = 0
do i=1,2
  if (Parent(AB(i),3-k) == AB(3-i)) Inbr = i
enddo

PrL = 0D0
do l=1,nSnp
  if (Par==Parent(A,k) .and. Par/=0) then
    call ParProb(l, Par, k, A, 0, PrX)
  else if (Par==Parent(B,k) .and. Par/=0) then
    call ParProb(l, Par, k, B, 0, PrX)
  else
    call ParProb(l, Par, k, 0, 0, PrX)    
  endif
  call ParProb(l, Parent(A,3-k), 3-k, A, 0, PrPx(:,1))
  call ParProb(l, Parent(B,3-k), 3-k, B, 0, PrPx(:,2))
  if (Inbr==0) then
    do x=1,3
      do y=1,3 
        PrXY(x,y) = PrX(x) * PrPX(y,1) * OKA2P(Genos(l,A),x,y) * &
          SUM(OKA2P(Genos(l,B),x,:) * PrPX(:,2))
      enddo
    enddo
    PrL(l) = LOG10(SUM(PrXY))
  else 
    do x=1,3
      do y=1,3
        PrXY(x,y)=PrX(x) * SUM(AKA2P(y, x, :) * PrPx(:,3-Inbr)) * &
          OKA2P(Genos(l,AB(Inbr)), x, y) * OcA(y,Genos(l, AB(3-Inbr)))
      enddo
    enddo
    PrL(l) = LOG10(SUM(PrXY))
  endif
enddo
LL = SUM(PrL)

!if (A==343 .and. B==190) print *, "PairHS: ", A, B, Par, ", ", LL

end subroutine PairHalfSib

! #####################################################################

subroutine pairHSHA(A, B, k, hf, LL, withFS)  !HS via k, & parent A is HS of B via 3-k
! hf 1: HSHA, 2: HSFA, 3: A inbred
use Global
implicit none

integer, intent(IN) :: A,B, k, hf
logical, intent(IN) :: withFS                             
double precision, intent(OUT) :: LL
integer :: l, x, y, z, PAB, Ai, Bj, i, j, exclFS
double precision :: PrL(nSnp), PrXYZ(3,3,3), PrY(3), PrZ(3), LLX

LL = missing
PAB = 0
if (Parent(A,3-k)/=0) then
  LL = NotImplemented
  return
else if (Parent(A,k)/=Parent(B,k)) then
  if (Parent(A,k)/=0) then
    if(Parent(B,k)/=0) then
      LL = impossible
      return
    else
      PAB = Parent(A,k)
    endif
  else
    PAB = Parent(B,k)
  endif
endif 

if (PAB < 0 .and. hermaphrodites/=0) then
  if (DumClone(-PAB,k)/=0) then
    LL = NotImplemented
    Return
  endif
endif   

Ai = FSID(maxSibSize+1, A)
Bj = FSID(maxSibSize+1, B)             
if (withFS) then
  exclFS = -1
else
  exclFS = 0
endif

PrL = 0D0
do l=1, nSnp
  call ParProb(l, Parent(B,3-k), 3-k, B, 0, PrY)
  if (Parent(A,k)==Parent(B,k)) then
    call ParProb(l, PAB, k, A, B, PrZ)
  else if (Parent(A,k) == PAB) then
    call ParProb(l, Parent(A,k), k, A, exclFS, PrZ)
  else if (Parent(B,k) == PAB) then
    call ParProb(l, Parent(B,k), k, B, exclFS, PrZ)
  endif
  
  do x=1,3
    do y=1,3    
      do z=1,3
        if (hf==1) then
          PrXYZ(x,y,z) = PrY(y) * PrZ(z) * AKAP(x, y, l)
        else if (hf==2) then
          PrXYZ(x,y,z) = PrY(y) * PrZ(z) * AKA2P(x, y, z)
        else if (hf==3) then
          PrXYZ(x,y,z) = PrY(y) * PrZ(z) * AKAP(x, z, l)
        endif      
        do i=1, nFS(Ai)
          if (FSID(i,Ai)/= A .and. .not. withFS)  cycle
          PrXYZ(x,y,z) = PrXYZ(x,y,z) * OKA2P(Genos(l,FSID(i,Ai)), x,z)
        enddo
        do j=1, nFS(Bj)
          if (FSID(j,Bj) /= B .and. .not. withFS)  cycle
          PrXYZ(x,y,z) = PrXYZ(x,y,z) * OKA2P(Genos(l,FSID(j,Bj)), y,z)
        enddo
      enddo
    enddo
  enddo
  PrL(l) = LOG10(SUM(PrXYZ))
enddo

if (.not. withFS) then
  LL = SUM(PrL)
else
  LLX = 0.0D0
  do i=1, nFS(Ai)
    if (FSID(i,Ai) /= A) then
      LLX = LLX + Lind(FSID(i,Ai))
    endif
  enddo
  do j=1, nFS(Bj)
    if (FSID(j,Bj) /= B) then
      LLX = LLX + Lind(FSID(j,Bj))
    endif
  enddo
  LL =  SUM(PrL) - LLX   !! ????
endif

if (LL < -HUGE(0D0))  LL = impossible

end subroutine pairHSHA

! #####################################################################

subroutine pairHSHAI(A, B, k, LL)  !HS via k, & A inbred
use Global
implicit none

integer, intent(IN) :: A, B, k
double precision, intent(OUT) :: LL
integer :: l, x, y, f,g
double precision :: PrL(nSnp,2,2), PrXY(3,3,2,2), PrPA(3), PrPAX(3), PrPB(3), LLU, LLtmp(2)
logical :: AncOK

LL = Missing
if (Parent(A,k)>0 .or. Parent(B,k)>0) then
  LL = NotImplemented
endif 
if (Parent(A,k) <0) then
  if (ns(-parent(A,k),k)/=1 .or. any(GpID(:,-parent(A,k),k)/=0))  LL = NotImplemented
endif
if (Parent(B,k) <0) then
  if (ns(-parent(B,k),k)/=1 .or. any(GpID(:,-parent(B,k),k)/=0))  LL = NotImplemented
endif
if (Parent(A, 3-k)==B) LL = impossible
if (LL/=Missing) return

call ChkAncest(A,0,B,0, AncOK)
if (.not. AncOK)  then
  LL = impossible
  return
endif
if (Parent(A,3-k)==Parent(B,3-k) .and. Parent(A,3-k)/=0) then
    LL = NotImplemented  ! likely picked up elsewhere
    return
endif

PrL = 0D0
do l=1, nSnp
  call ParProb(l, Parent(A,3-k), 3-k, A, 0, PrPA)
  call ParProb(l, Parent(A,3-k), 3-k, A, -4, PrPAX)  ! OcA if Parent(A,3-k)>0
  call ParProb(l, Parent(B,3-k), 3-k, B, 0, PrPB)
  do x=1,3
    do y=1,3    
      PrXY(x,y,1,:) = AKAP(x,y,l) * PrPAX(x) * AHWE(y,l) 
      PrXY(x,y,2,:) = AKAP(y,x,l) * PrPA(x)
      PrXY(x,y,:,:) = PrXY(x,y,:,:) * OKA2P(Genos(l,A), x,y)
      PrXY(x,y,:, 1) = PrXY(x,y,:, 1) * SUM(OKA2P(Genos(l,B), y,:)*PrPB)
      PrXY(x,y,:, 2) = PrXY(x,y,:, 2) * SUM(OKAP(Genos(l,B),:,l)*PrPB)
    enddo
  enddo
  do f=1,2
    do g=1,2
      PrL(l,f,g) = LOG10(SUM(PrXY(:,:,f,g)))
    enddo
  enddo
enddo

LLU = missing
LLtmp = missing
call CalcU(A,k,B,k, LLU)            
do f=1,2
  if (SUM(PrL(:,f,2)) > LLU) then
    LLtmp(f) = SUM(PrL(:,f,1)) - SUM(PrL(:,f,2)) + LLU
  else
    LLtmp(f) = SUM(PrL(:,f,1))
  endif
enddo
LL = MaxLL(LLtmp)

end subroutine pairHSHAI

! #####################################################################

subroutine pairFAHA(A, B, withFS, LL)  !B FA via k & HA via 3-k ; A inbred. 
use Global
implicit none

integer, intent(IN) :: A, B
logical, intent(IN) :: withFS
double precision, intent(OUT) :: LL
integer :: l, k, x, y, z, v, i, AA(maxSibSize), BB(maxSibSize), nA, nB
double precision :: PrL(nSnp, 4), PrXY(3,3,3,3, 4), PrPB(3,2), LLU, LLtmp(4)

LL = missing
if (ANY(Parent(A,:)>0)) then  
  LL = NotImplemented   
  return
endif

AA = 0
BB = 0
if (withFS) then
  nA = nFS(FSID(maxSibSize+1, A))  
  AA(1:nA) = FSID(1:nA, FSID(maxSibSize+1, A))
  nB = nFS(FSID(maxSibSize+1, B))
  BB(1:nB) = FSID(1:nB, FSID(maxSibSize+1, B))
else
  nA = 1
  AA(1) = A
  nB = 1
  BB(1) = B
endif

PrL = 0D0
do l=1, nSnp
  do k=1,2
    if (withFS) then
      call ParProb(l, Parent(B,k), k, B, -1, PrPB(:,k))
    else
      call ParProb(l, Parent(B,k), k, B, 0, PrPB(:,k))
    endif
  enddo
  do x=1,3  ! Par A, FS of B
    do y=1,3   ! par B, double GP of A 
      do z=1,3  ! par B, GP of A
        do v=1,3  ! Par A, HS of B
          PrXY(x,y,z,v,1) = PrPB(y,1) * PrPB(z,2) * AKAP(v,y,l) * AKA2P(x,y,z)
          PrXY(x,y,z,v,2) = PrPB(y,2) * PrPB(z,1) * AKAP(v,y,l) * AKA2P(x,y,z)
          PrXY(x,y,z,v,3) = PrPB(y,1) * PrPB(z,2) * AHWE(v,l) * AHWE(x,l)  ! A, B unrelated
          PrXY(x,y,z,v,4) = PrPB(y,1) * PrPB(z,2) * SUM(AKAP(v,:,l) * AKAP(x,:,l) * AHWE(:,l))
            ! A, B unrelated; A inbred
          do i=1, nA
            PrXY(x,y,z,v,:) = PrXY(x,y,z,v,:) * OKA2P(Genos(l,AA(i)), x, v)
          enddo
          do i=1, nB
            PrXY(x,y,z,v,:) = PrXY(x,y,z,v,:) * OKA2P(Genos(l,BB(i)), y, z)
          enddo
        enddo
      enddo
    enddo 
  enddo
  do i=1,4
    PrL(l,i) = LOG10(SUM(PrXY(:,:,:,:,i)))
  enddo
enddo
LLtmp = SUM(PrL, DIM=1)

if (.not. withFS) then
  LL = MAXVAL(LLtmp(1:2))
else
  LLU = missing
  call CalcU(A, 0, B, 0, LLU)
  LL = MAXVAL(LLtmp(1:2)) - MAXVAL(LLtmp(3:4)) + LLU
endif
end subroutine pairFAHA

! #####################################################################

subroutine pairHSPO(A, B, LL)   ! HS via k, & PO via 3-k
use Global
implicit none

integer, intent(IN) :: A,B
double precision, intent(OUT) :: LL
integer :: l, x, y
double precision :: PrL(nSnp), PrXY(3,3)

if (ANY(Parent(A,:)/=0) .or. ANY(Parent(B,:)/=0)) then
  LL = impossible
  return   ! else not necessary.
endif  

PrL = 0D0
do l=1, nSnp
  do x=1,3 
    do y=1,3    ! B
      PrXY(x,y) = AHWE(x,l) * AKAP(y,x,l) * OKA2P(Genos(l,A), x, y) * OcA(y,Genos(l,B)) 
    enddo
  enddo
  PrL(l) = LOG10(SUM(PrXY))
enddo
LL = SUM(PrL)

end subroutine pairHSPO

! #####################################################################

subroutine PairPOHA(A, B, k, hf, LL)  ! B parent of A via k, and 'hf' sib of A's 3-k parent
use Global
implicit none

integer, intent(IN) :: A, B, k, hf
double precision, intent(OUT) :: LL
integer :: l, x, y,z,m, ParB(2), GA(2), GG(2)
double precision :: PrL(nSnp), PrXYZ(3,3,3), PrPAX(3), PrGG(3,2), PrPB(3), PrGA(3)

ParB = Parent(B,:)
GA = getPar(Parent(A,3-k), 3-k)
GG = 0
LL = missing

if (Parent(A,3-k) == ParB(3-k) .and. ParB(3-k)/=0) then
  LL = impossible
  return
endif

do m=1,2
  if (m/=hf .and. hf/=3)  cycle
  if (GA(m) == B) then
    LL = impossible
  else if (ParB(m) == GA(m) .or. GA(m) == 0) then
    GG(m) = ParB(m)
  else if (ParB(m) == 0) then
    GG(m) = GA(m)
  else
    LL = impossible
  endif
  if (GG(m) < 0 .and. hermaphrodites/=0) then
    if (DumClone(-GG(m),m)/=0) then
      LL = NotImplemented
    endif
  endif
enddo
if (LL /= Missing)  return

PrL = 0D0  
PrPB = missing
PrGA = missing
do l=1,nSnp       
  call ParProb(l, Parent(A,3-k), 3-k, A, -4, PrPAX)  ! offspring only, except A
  do m=1,2
    if (m/=hf .and. hf/=3)  cycle
    if (Parent(A,3-k)>0) then
      call ParProb(l, GG(m), m, B, Parent(A,3-k), PrGG(:,m))
    else
      call ParProb(l, GG(m), m, B, 0, PrGG(:,m))
    endif
  enddo
  if (hf < 3) then
    call ParProb(l, Parent(B,3-hf), 3-hf, B, 0, PrPB)
    if (Parent(A,3-k)>0) then
      call ParProb(l, GA(3-hf), 3-hf, Parent(A,3-k), 0, PrGA)
    else
      call ParProb(l, GA(3-hf), 3-hf, 0, 0, PrGA)
    endif
  endif
  
  PrXYZ = 0D0
  do x=1,3  ! B
    do y=1,3  ! parent(A,3-k)
      do z=1,3  ! double grandparent (dam if hf=3)
        if (hf < 3) then
          PrXYZ(x,y,z) = SUM(AKA2P(x,z,:) * PrPB) * &
            PrPAX(y) * SUM(AKA2P(y,z,:) * PrGA) * PrGG(z, hf)
        else
          PrXYZ(x,y,z) = SUM(AKA2P(x,z,:) * PrPAX(y) * AKA2P(y,z,:) * &
            PrGG(z, 1) * PrGG(:,2))
        endif
        PrXYZ(x,y,z) = PrXYZ(x,y,z) * OKA2P(Genos(l, A), x, y) * OcA(x,Genos(l, B))
      enddo
    enddo
  enddo
  PrL(l) = LOG10(SUM(PrXYZ))
enddo

LL = SUM(PrL)
if (LL < -HUGE(0D0))  LL = impossible

if (Parent(A,3-k)>0 .and. LL/=impossible) then
  LL = LL - Lind(Parent(A,3-k))
endif

! if ((A==80 .and. B==51) .or. (A==51 .and. B==80)) then
  ! write (*, '("PairPOHA : ", 3i6, " + ", 2i3, 3i6, f8.1, 4i6)')  A, Parent(A,:), &
    ! k, hf, B, Parent(B,:), LL, GA, GG
! endif

end subroutine PairPOHA

! #####################################################################

subroutine clustHSHA(SA, SB, k, LL)   ! HS via 3-k, & SB parent of SA; SA,SB FS
use Global
implicit none

integer, intent(IN) :: SA,SB, k
double precision, intent(OUT) :: LL
integer :: l, x, y, z,i, Par(2), GC(2), u
double precision :: PrL(nSnp), PrXYZ(3,3,3), PrGA(3),PrGC(3,2),PrUZ(3,3)

! all checks done by CheckMerge.

! grandparents of opp. parent
 call getFSpar(SA, k, .TRUE., Par(1))  ! TODO: or strict=.FALSE.?
 call getFSpar(SB, k, .TRUE., Par(2))
GC = 0
do i=1,2
  if (Par(1)<0) then
    GC(i) = GpID(i, -Par(1),3-k)
    if (GpID(i, -Par(1),3-k)/=0) then
      if (Par(2) < 0) then
       if (GpID(i,-Par(2),3-k)/=GC(i) .and. GpID(i,-Par(2),3-k)/=0) then
          GC(i) = 0   ! shouldn't happen
        else if (GC(i)==0 .and. GpID(i, -Par(2),3-k)/=0) then
          GC(i) = GpID(i, -Par(2),3-k)
        endif
      endif
    endif
  endif
enddo

PrL = 0D0
do l=1, nSnp
  call ParProb(l, GpID(3-k,SA,k), 3-k, 0, 0, PrGA)
  do i=1,2
    call ParProb(l, GC(i), i, 0, 0, PrGC(:,i))
  enddo
  do z=1,3
    do u=1,3
        PrUZ(u,z) = SUM(AKA2P(z,u,:) * PrGC(u,1) * PrGC(:,2))
    enddo
    do x=1,3    
      do y=1,3
        PrXYZ(x,y,z) = SUM(AKA2P(x,y,:) * PrGA) * XPr(2,y,l,SB,k) *&
          SUM(PrUZ(:,z))
        do i=1,nS(SA,k)
          PrXYZ(x,y,z) = PrXYZ(x,y,z) *OKA2P(Genos(l,SibID(i,SA,k)),x,z)
        enddo 
        do i=1,nS(SB,k)
          PrXYZ(:,y,z) =PrXYZ(:,y,z) *OKA2P(Genos(l, SibID(i,SB,k)),y,z)
        enddo
      enddo
    enddo
  enddo
  PrL(l) = LOG10(SUM(PrXYZ))
enddo
LL = SUM(PrL)

end subroutine clustHSHA

! #####################################################################

subroutine FSHC(A, B, k, LL)  ! FS + parents are HS; B may be neg
use Global
implicit none

integer, intent(IN) :: A,B, k
double precision, intent(OUT) :: LL
integer :: l, x, y, z, Par(2), m, GG(2,2), kG, i, PM(2)
double precision :: PrL(nSnp), PrG(3,2), PrXYZ(3,3,3), PrZ(3)

LL = missing
Par = 0            
if (B < 0 .and. A>0) then
  Par(k) = B
  if (Parent(A,k)/=0) then
    LL = impossible
    return
  endif
  if (ALL(Parent(SibID(1:nS(-B,k), -B, k), 3-k) == 0)) then
    Par(3-k) = Parent(A, 3-k)
  else
    call getFSpar(-B, k, .TRUE., Par(3-k))
    if (Par(3-k)==0 .or. (Parent(A,3-k)/=Par(3-k) .and. Parent(A, 3-k)/=0)) then
      LL = impossible
      return
    endif
  endif
else if (B > 0 .and. A>0) then
  do m=1,2
    if (Parent(B,m)==0) then
      Par(m) = Parent(A,m)
    else if (Parent(B,m) /= Parent(A,m) .and. Parent(A,m)/=0) then
      LL = impossible
      return
    else
      Par(m) = Parent(B,m)
    endif
  enddo
else if (B<0 .and. A<0) then
  if (ANY(GpID(:,-B,k)/=0) .or. ANY(GpID(:,-A,k)/=0)) then
    LL = NotImplemented
    return
  else
    Par = 0
  endif
endif

if (hermaphrodites/=0) then
  LL = NotImplemented
  return
endif

GG = 0
kG = 0
PM = 0
do m=1,2
  GG(:, m) = getPar(Par(m), m)
enddo
do m=1,2
    if (GG(m,1)==0 .or. GG(m,2)==0) then  ! GG(m,1)==GG(m,2) not needed
      kG = m
    endif
    if (Par(m)>0) then
      PM(m) = Par(m)
    endif
enddo
if (kG==0) then
    LL = AlreadyAss
    return
endif

PrL = 0D0
do l=1, nSnp
  do m=1,2
    call ParProb(l, GG(3-kG, m), 3-kG, PM(m),0, PrG(:,m))
  enddo
  if (GG(kG,1)/=0) then
    call ParProb(l, GG(kG,1), kG, PM(1),0, PrZ)
  else
    call ParProb(l, GG(kG,2), kG, PM(2),0, PrZ)
  endif
  do x=1,3  ! Par(1)
    do y=1,3  !Par(2)
      do z=1,3  ! GG(kG)
        PrXYZ(x,y,z) = PrZ(z) * SUM(AKA2P(x, z, :) * PrG(:,1)) * &
          SUM(AKA2P(y, z, :) * PrG(:,2))
        if (A>0) then
          PrXYZ(x,y,z) = PrXYZ(x,y,z) * OKA2P(Genos(l,A), x, y)
        else
          do i=1, nS(-A,k)
            PrXYZ(x,y,z) = PrXYZ(x,y,z) * OKA2P(Genos(l,SibID(i,-A,k)), x, y)
          enddo
        endif
        if (B>0) then
          PrXYZ(x,y,z) = PrXYZ(x,y,z) * OKA2P(Genos(l,B), x, y)
        else
          do i=1, nS(-B,k)
            PrXYZ(x,y,z) = PrXYZ(x,y,z) * OKA2P(Genos(l,SibID(i,-B,k)), x, y)
          enddo
        endif
      enddo
    enddo
  enddo
  PrL(l) = LOG10(SUM(PrXYZ))
enddo
LL = SUM(PrL)

end subroutine FSHC

! #####################################################################

subroutine pairFSHA(A, B, k, LL) !inbred FS: par k offspring of par 3-k
use Global
implicit none

integer, intent(IN) :: A,B,k
double precision, intent(OUT) :: LL
integer :: l, x, y
double precision :: PrL(nSnp), PrXY(3,3), PrY(3)

if (Parent(A,k)/=0 .or. Parent(B,k)/=0) then
  LL = NotImplemented 
  return
endif  

PrL = 0D0
do l=1, nSnp
  call ParProb(l, Parent(A,3-k), 3-k, -1,0, PrY) 
  do x=1,3
    do y=1,3    
      PrXY(x,y) = PrY(y) * AKAP(x, y, l) * OKA2P(Genos(l,B), x, y) * OKA2P(Genos(l,A), x, y)
    enddo
  enddo
  PrL(l) = LOG10(SUM(PrXY))
enddo
LL = SUM(PrL)

end subroutine pairFSHA

! #####################################################################

subroutine pairFSselfed(A, B, LL) !A & B both product of selfing by same parent
use Global
implicit none

integer, intent(IN) :: A,B
double precision, intent(OUT) :: LL
integer :: l, x
double precision :: PrL(nSnp), PrX(3)

if (any(Parent(A,:)/=0) .or. any(Parent(B,:)/=0)) then
  LL = NotImplemented 
  return
endif  

PrL = 0D0
do l=1, nSnp
  do x=1,3
    PrX(x) = OKA2P(Genos(l,B), x, x) * OKA2P(Genos(l,A), x, x) * AHWE(x,l)
  enddo
  PrL(l) = LOG10(SUM(PrX))
enddo
LL = SUM(PrL)

end subroutine pairFSselfed

! #####################################################################

subroutine trioGP(A,kA, B, kB, C, kC, LL)  ! B & C both GP of A, A>0 or A<0
use Global
implicit none

integer, intent(IN) :: A,B,C, kA,kB, kC
double precision, intent(OUT) :: LL
integer :: l, x, y, v, BC(2), kBC(2), oppar_BC(2), w, i
double precision :: PrL(nSnp,3), PrXY(3,3), PrXYV(3,3,3), PrB(3), PrC(3),PrA(3), &
  LLBC, PrMB(3), PrMC(3), PrX(3,3), PrY(3,3), LLtmp(3)
logical :: withMate(2)

BC = (/B,C/)
kBC = (/kB, kC/)
withMate = .FALSE.
if (B<0 .and. C<0) then  !  with mate of B & C
  do x=1,2
    call getFSpar(-BC(x),kBC(x), .TRUE., opPar_BC(x))
    if (opPar_BC(x)/=0 .or. all(parent(SibID(1:ns(-BC(x),kBC(x)),-BC(x),kBC(x)), 3-kBC(x)) == 0)) then
      withMate(x) = .TRUE.
    endif
  enddo
endif

! parents are -not- ignored, but assumes Par(A,:)=0 
PrL = 0D0
do l=1, nSnp
  PrXY = 0D0
  call OffProb(l, A, kA, PrA)  ! OcA if >0, XPr(1,:,) if <0
  call ParProb(l, B, kB, 0, 0, PrB)
  call ParProb(l, C, kC, 0, 0, PrC)
  do x=1,3
    do y=1,3
      ! one mat GP, other pat GP 
      PrXY(x,y) = SUM(PrA * AKA2P(:, x, y)) * SUM(AKAP(x,:,l) * PrB) * &
        SUM(AKAP(y,:,l) * PrC)
      ! both mat or both pat GP
      do v=1,3
        PrXYV(x,y,v) = SUM(PrA * AKA2P(:, x, y)) * AHWE(y,l) * SUM(AKA2P(x,v,:) * PrB(v) * PrC)
      enddo
    enddo
  enddo
  PrL(l,1) = LOG10(SUM(PrXY))
  PrL(l,2) = LOG10(SUM(PrXYV))
enddo

if (all(withMate)) then   ! TODO?: any   ! B<0 & C<0
  do l=1, nSnp
    PrXY = 0D0
    call OffProb(l, A, kA, PrA) 
    call ParProb(l, B, kB, -1, 0, PrB)
    call ParProb(l, opPar_BC(1), 3-kB, -1, 0, PrMB)
    call ParProb(l, C, kC, -1, 0, PrC)
    call ParProb(l, opPar_BC(2), 3-kC, -1, 0, PrMC)
    do x=1,3
      do y=1,3
        PrXY(x,y) = SUM(PrA * AKA2P(:, x, y))
        do v=1,3
          do w=1,3
            PrX(v,w) = AKA2P(x,v,w) * PrB(v) * PrMB(w)
            do i=1,ns(-B,kB)
              PrX(v,w) = PrX(v,w) * OKA2P(Genos(l,SibID(i,-B,kB)),v,w)
            enddo
          enddo
        enddo
        do v=1,3
          do w=1,3
            PrY(v,w) = AKA2P(y,v,w) * PrC(v) * PrMC(w)
            do i=1,ns(-C,kC)
              PrX(v,w) = PrX(v,w) * OKA2P(Genos(l,SibID(i,-C,kC)),v,w)
            enddo
          enddo
        enddo
        PrXY(x,y) = PrXY(x,y) * SUM(PrX) * SUM(PrY)
      enddo
    enddo
    PrL(l,3) = LOG10(SUM(PrXY))
  enddo
endif

LLtmp = SUM(PrL,DIM=1)
if (.not. all(withMate))  LLtmp(3) = NotImplemented
LLBC = missing
call CalcU(B, kB, C, kC, LLBC)

! if (A==-208 .and. B==4846 .and. C==3920)  then
  ! write(*,'("TrioGP: ", 3i5, 2l3, 4f8.1, 2i5)') A, B, C, withMate, LLtmp, LLBC, opPar_BC
! endif

LLtmp(1:2) = LLtmp(1:2) + LLBC

LL = MAXLL(LLtmp)
if (LL < -HUGE(0D0))  LL = impossible

end subroutine trioGP

! #####################################################################

subroutine trioFS(A,kA, B,kB, C,kC, LL) 
use Global
implicit none

integer, intent(IN) :: A,kA, B,kB, C,kC
double precision, intent(OUT) :: LL
integer :: l, x, y, iX(3), kX(3), s, ParT(2), kP, ParX(2), i_excl(2)
double precision :: PrL(nSnp), PrXY(3,3), PrS(3,3), PrP(3,2)

iX = (/A, B, C/)
kX = (/kA, kB, kC/)

if (Hermaphrodites/=0) then
  do s=1,3
    if (iX(s) >=0)  cycle
    if (DumClone(iX(s), kX(s))/=0) then
      LL = NotImplemented
      Return
    endif
  enddo
endif

ParT = 0  ! trio parents
i_excl = 0  ! sibship members to exclude when calling ParProb()
do s=1,3  
  ParX = getPar(iX(s), kX(s))
  do kP=1,2
    if (ParX(kP)/= ParT(kP) .and. ParT(kP)/=0) then
      LL = Impossible
      Return
    else
      ParT(kP) = ParX(kP)
      if (iX(s) > 0 .and. ParT(kP) < 0) then
        i_excl(kP) = iX(s)   ! TODO?: 2 already have ParT as parent
      endif
    endif
  enddo
enddo

do kP=1,2
  do s=1,3  
    if (kX(s)/=kP)  cycle
    if (ParT(kP) == iX(s)) then
      LL = Impossible
      Return 
    endif
  enddo
enddo

PrL = 0D0
do l=1, nSnp
  do kP=1,2
    call ParProb(l,ParT(kP),kP, i_excl(kP),0, PrP(:,kP))
  enddo
  do s=1,3
    call OffProb(l,iX(s), kX(s), PrS(:,s))
  enddo
  do x=1,3
    do y=1,3
      PrXY(x,y) = PrP(x,1) * PrP(y,2) 
      do s=1,3       
        PrXY(x,y) = PrXY(x,y) * SUM(PrS(:,s) * AKA2P(:,x,y))
      enddo   
    enddo
  enddo
  PrL(l) = LOG10(SUM(PrXY))
enddo

LL = SUM(PrL)

end subroutine trioFS

! #####################################################################

subroutine trioHS(A,kA, B,kB, C,kC, LL) 
use Global
implicit none

integer, intent(IN) :: A,kA, B,kB, C,kC
double precision, intent(OUT) :: LL
integer :: l, x, y, iX(3), kX(3), s
double precision :: PrL(nSnp, 2), PrX(3), PrXY(3,3), PrS(3,3)

iX = (/A, B, C/)
kX = (/kA, kB, kC/)

PrL = 0D0
do l=1, nSnp
  do s=1,3
    call OffProb(l,iX(s), kX(s), PrS(:,s))    ! OcA if >0, XPr(1,:,) if <0
  enddo

  ! all 3 are halfsibs
  do x=1,3
    PrX(x) = AHWE(x,l)
    do s=1,3
      PrX(x) = PrX(x) * SUM(PrS(:,s) * AKAP(:,x,l))
    enddo
  enddo
  PrL(l,1) = LOG10(SUM(PrX))
  
  ! A-B HS and A-C HS, but not B-C
  do x=1,3
    do y=1,3
      PrXY(x,y) = SUM(PrS(:,1) * AKA2P(:,x,y)) * SUM(PrS(:,2) * AKAP(:,x,l)) * &
        SUM(PrS(:,3) * AKAP(:,y,l)) * AHWE(x,l) * AHWE(y,l)
    enddo
  enddo
  PrL(l,2) = LOG10(SUM(PrXY))      
enddo

LL = MAXVAL(SUM(PrL, DIM=1))

end subroutine trioHS

! #####################################################################

! subroutine trioFHS(A,kA, B,kB, C,kC, LL)   ! A-B FS, all 3 HS 
! use Global
! implicit none

! integer, intent(IN) :: A,kA, B,kB, C,kC
! double precision, intent(OUT) :: LL
! integer :: l, x, y, iX(3), kX(3), s
! double precision :: PrL(nSnp), PrXY(3,3), PrS(3,3)

! iX = (/A, B, C/)
! kX = (/kA, kB, kC/)

! if (Hermaphrodites/=0) then
  ! do s=1,3
    ! if (iX(s) >=0)  cycle
    ! if (DumClone(iX(s), kX(s))/=0) then
      ! LL = NotImplemented
      ! Return
    ! endif
  ! enddo
! endif

! PrL = 0D0
! do l=1, nSnp
  ! do s=1,3
    ! call OffProb(l,iX(s), kX(s), PrS(:,s))
  ! enddo
  ! do x=1,3
    ! do y=1,3
      ! PrXY(x,y) = AHWE(x,l) * AHWE(y,l)
      ! do s=1,2  ! A, B are FS   
        ! PrXY(x,y) = PrXY(x,y) * SUM(PrS(:,s) * AKA2P(:,x,y))
      ! enddo   
      ! PrXY(x,y) = PrXY(x,y) * SUM(PrS(:,3) * AKAP(:,x,l))
    ! enddo
  ! enddo
  ! PrL(l) = LOG10(SUM(PrXY))
! enddo

! LL = SUM(PrL)

! end subroutine trioFHS

! #####################################################################

subroutine trioFA(A,kA, B,C, LL)   ! B & C >0
use Global
implicit none

integer, intent(IN) :: A,kA, B, C
double precision, intent(OUT) :: LL
integer :: l, k, x, y, z, v,u,w
double precision :: PrL(nSnp), PrXY(3,3), PrPB(3,2), PrPC(3,2), PrUV(3,3), &
   PrWZ(3,3), PrA(3)

if ((all(parent(B,:)==0) .or. all(Parent(C,:)==0)) .and. A>0) then
  LL = NotImplemented   ! LL would be identical to GP+GP or HS+HS 
  return
endif

PrL = 0D0
do l=1, nSnp
  call OffProb(l, A, kA, PrA)  ! OcA if >0, XPr(1,:,) if <0
  do k=1,2
    call ParProb(l,parent(B,k),k,-1,0,PrPB(:,k))
    call ParProb(l,parent(C,k),k,-1,0,PrPC(:,k))
  enddo

  do x=1,3
    do y=1,3     
      do u=1,3
        do v=1,3
          PrUV(u,v) = PrPB(u,1) * PrPB(v,2) * AKA2P(x,u,v) * OKA2P(Genos(l,B),u,v)      
        enddo
      enddo
      do w=1,3
        do z=1,3
          PrWZ(w,z) = PrPC(w,1) * PrPC(z,2) * AKA2P(y,w,z) * OKA2P(Genos(l,C),w,z) 
        enddo
      enddo
      PrXY(x,y) = SUM(PrA * AKA2P(:,x,y)) * SUM(PrUV) * SUM(PrWZ)
    enddo
  enddo
  PrL(l) = LOG10(SUM(PrXY))
enddo

LL = SUM(PrL)

end subroutine trioFA

! #####################################################################

subroutine trioHA(A,kA, B,C, LL)  ! B & C > 0
use Global
implicit none

integer, intent(IN) :: A,kA, B, C
double precision, intent(OUT) :: LL
integer :: l, x, y, z, v, PB, PC, kB, kC
double precision :: PrL(nSnp,2), PrXYZ(3,3,3,2), PrV(3), PrA(3), PrPB(3), PrPC(3)

PB = 0
PC = 0
if (all(Parent(B,:)/=0) .or. all(Parent(C,:)/=0)) then
  LL = NotImplemented   ! TODO? 
  return
else
  if (Parent(B,1)/=0) then
    PB = Parent(B,1)
    kB=1
  else
    PB = Parent(B,2)
    kB = 2
  endif
  if (Parent(C,1)/=0) then
    PC = Parent(C,1)
    kC=1
  else
    PC = Parent(C,2)
    kC = 2
  endif
endif

PrL = 0D0
do l=1, nSnp
  call OffProb(l, A, kA, PrA)  ! OcA if >0, XPr(1,:,) if <0
  call ParProb(l, PB,kB, B, 0, PrPB)
  call ParProb(l, PC,kC, C, 0, PrPC)
  do x=1,3
    do y=1,3
      do z=1,3
        PrXYZ(x,y,z,:) = SUM(PrA * AKA2P(:,x,y)) * AKAP(x,z,l) * AHWE(z,l) * &
          SUM(OKA2P(Genos(l,B),z,:) * PrPB)
        ! B+C HS
        PrXYZ(x,y,z,1) = PrXYZ(x,y,z,1) * AHWE(y,l) * SUM(OKA2P(Genos(l,C),z,:) * PrPC)         
        ! B+C U
        do v=1,3
          PrV(v) = AKAP(y,v,l) * AHWE(v,l) * SUM(OKA2P(Genos(l,C),v,:) * PrPC)   
        enddo
        PrXYZ(x,y,z,2) = PrXYZ(x,y,z,2) * SUM(PrV)
      enddo
    enddo
  enddo
  do v=1,2
    PrL(l,v) = LOG10(SUM(PrXYZ(:,:,:,v)))
  enddo
enddo

LL = MAXVAL(SUM(PrL, DIM=1))

end subroutine trioHA

! #####################################################################

subroutine trioGGP(A,kA, B,C, LL)  ! B & C > 0
use Global
implicit none

integer, intent(IN) :: A,kA, B, C
double precision, intent(OUT) :: LL
integer :: l, x, y, z, v
double precision :: PrL(nSnp), PrXV(3,3,3,3), PrA(3), PrB(3), PrC(3)

PrL = 0D0
do l=1, nSnp
  call OffProb(l, A, kA, PrA)  ! OcA if >0, XPr(1,:,) if <0
  call ParProb(l, B,1, 0,0, PrB)
  call ParProb(l, C,1, 0,0, PrC)
  do x=1,3
    do y=1,3
      do z=1,3
        do v=1,3
          PrXV(x,y,z,v) = SUM(PrA * AKA2P(:,x,y)) * AKAP(x,z,l) * AKAP(y,v,l) * &
            SUM(AKAP(z,:,l) * PrB) * SUM(AKAP(v,:,l) * PrC)
        enddo
      enddo
    enddo
  enddo
  PrL(l) = LOG10(SUM(PrXV))
enddo

LL = SUM(PrL, DIM=1) + Lind(B) + Lind(C)

end subroutine trioGGP

! #####################################################################

subroutine trioHSGP(A,kA, B,kB, C,kC, k, LL)  ! A & B HS, C GP of both, via k
use Global
implicit none

integer, intent(IN) :: A,kA, B,kB, C,kC, k
double precision, intent(OUT) :: LL
integer :: l, x, PB(2), GB(2), BB, y,w, PC(2)
double precision :: PrL(nSnp), PrX(3), PrA(3), PrPB(3,2), PrC(3), PrMC(3), PrY(3),&
  PrB(3), PrW(3)

PB = getPar(B,kB)
GB = getPar(PB(k),k)
if (GB(kC)/=0 .and. GB(kC)/=C) then
  LL = impossible
  return
endif

PC = getPar(C, kC)
if (PC(kB) == B .or. PB(kC)==C) then
  LL = NotImplemented
  return
endif

PrL = 0D0
if (B>0) then
  BB = B
else
  BB = 0
endif
do l=1, nSnp
  call OffProb(l, A, kA, PrA)  ! OcA if >0, XPr(1,:,) if <0
  call OffProb(l, B, kB, PrB)
  call ParProb(l, C, kC, 0, 0, PrC)
  call ParProb(l, PB(3-k),3-k,BB,0,PrPB(:,3-k))
  call ParProb(l, PB(k),k, BB, -4, PrPB(:,k))
  if (PB(k)>0) then
    call ParProb(l,GB(3-kC),3-kC,PB(k),0,PrMC)
  else
    call ParProb(l,GB(3-kC),3-kC,0,0,PrMC)
  endif
  do x=1,3
    do y=1,3
      PrY(y) = SUM(AKA2P(x,y,:) * PrC(y) * PrMC)
    enddo
    do w=1,3
      PrW(w) = SUM(PrB(w) * AKA2P(w,x,:) * PrPB(:,3-k))
    enddo
    PrX(x) = SUM(PrA * AKAP(:,x,l)) * SUM(PrW) * PrPB(x,k) * SUM(PrY)
  enddo
  PrL(l) = LOG10(SUM(PrX))
enddo

if (C > 0) then
  LL = SUM(PrL) + Lind(C)
else
  LL = SUM(PrL) + CLL(-C,kC)
endif
if (PB(k)>0)  LL = LL - Lind(PB(k))

end subroutine trioHSGP

! #####################################################################

! subroutine trio_OP(A,kA, B,kB, C,kC, k, LL)  ! B is parent of A, C is offspring of A
! use Global
! implicit none

! integer, intent(IN) :: A,kA, B,kB, C,kC, k
! double precision, intent(OUT) :: LL
! double precision :: PrL(nSnp), PrX(3), PrA(3), PrB(3), PrC(3) 
! integer :: PC(2), x, l
! logical :: AncOK

! PC = getPar(C, kC)
! if (any(PC /= 0)) then
  ! LL = NotImplemented
  ! return
! endif

! call ChkAncest(B, kB, C, kC, AncOK)  ! check that C is not an ancestor of B
! if (.not. AncOK) then
  ! LL = Impossible
  ! return
! endif

! PrL = 0D0
! do l=1, nSnp
  ! call OffProb(l, A, kA, PrA)   ! OcA if >0, XPr(1,:,) if <0
  ! call OffProb(l, C, kC, PrC)  
  ! call ParProb(l, B, kB, 0, 0, PrB)
  ! do x=1,3
    ! PrX(x) = SUM(PrC * AKAP(:,x,l)) * PrA(x) * SUM(AKAP(x,:,l) * PrB)
  ! enddo
  ! PrL(l) = LOG10(SUM(PrX))
! enddo

! if (B > 0) then
  ! LL = SUM(PrL) + Lind(B)
! else
  ! LL = SUM(PrL) + CLL(-B,kB)
! endif

! end subroutine trio_OP

! #####################################################################

! subroutine trioHSGP2(A,kA, B,kB, C,kC, k, LL)  ! A & B HS, C GP of only A, via k
! use Global
! implicit none

! integer, intent(IN) :: A,kA, B,kB, C,kC, k
! double precision, intent(OUT) :: LL
! integer :: l, x,y,z, PB(2), BB, w
! double precision :: PrL(nSnp), PrXYZ(3,3,3), PrA(3), PrB(3), PrPB(3,2), PrC(3),PrW(3)

! PB = getPar(B,kB)

! PrL = 0D0
! if (B>0) then
  ! BB = B
! else
  ! BB = 0
! endif
! do l=1, nSnp
  ! call OffProb(l, A, kA, PrA)  ! OcA if >0, XPr(1,:,) if <0
  ! call OffProb(l, B, kB, PrB)
  ! call ParProb(l, C, kC, 0, 0, PrC)
  ! call ParProb(l, PB(3-k),3-k,BB,0,PrPB(:,3-k))
  ! call ParProb(l, PB(k),k, BB, 0, PrPB(:,k))
  ! do x=1,3  ! parent of A, offspring of C
    ! do y=1,3  ! C
      ! do z=1,3  ! shared parent of A & B
        ! do w=1,3
          ! PrW(w) = SUM(PrB(w) * AKA2P(w,:,z) * PrPB(:,k))
        ! enddo
        ! PrXYZ(x,y,z) = SUM(PrA * AKA2P(:,x,z)) * PrPB(z,3-k) * AKAP(x,y,l) * &
          ! PrC(y) * SUM(PrW)          
      ! enddo
    ! enddo
  ! enddo
  ! PrL(l) = LOG10(SUM(PrXYZ))
! enddo

! LL = SUM(PrL)
! if (C > 0) then
  ! LL = LL + Lind(C)
! else
  ! LL = LL + CLL(-C,kC)
! endif

! end subroutine trioHSGP2

! #####################################################################

subroutine pairHSGP(A, B,k, LL)   ! HS via k, B is GP of A via 3-k
use Global
implicit none

integer, intent(IN) :: A,B,k
double precision, intent(OUT) :: LL
integer :: l, x, y, z
double precision :: PrL(nSnp), PrXYZ(3,3,3), PrPB(3)

if (Parent(A,3-k)/=0) then
  LL = NotImplemented
  return
endif 

PrL = 0D0
do l=1, nSnp
  call ParProb(l, Parent(B,3-k), 3-k, B, 0, PrPB)
  do x=1,3  ! parent 3-k of A, offspring of B
    do y=1,3  ! shared parent k 
      do z=1,3  ! B
        PrXYZ(x,y,z) =AKAP(x,z,l)*AHWE(y,l)*SUM(AKA2P(z,y,:)*PrPB) * &
          OKA2P(Genos(l,A),x, y) * OcA(z,Genos(l,B))
      enddo
    enddo
  enddo
  PrL(l) = LOG10(SUM(PrXYZ))
enddo
LL = SUM(PrL)

end subroutine pairHSGP

! #####################################################################

subroutine PairGP(A, B, k, focal, LL)  
! calculates LL that B is maternal(k=1) or paternal(k=2) gp of A
use Global
implicit none

integer, intent(IN) :: A,B,k, focal
double precision, intent(OUT) :: LL
integer :: l, x, y, curGP(2,2), m, z, v, i, ParAB
double precision :: PrL(nSnp,8), PrPA(3,2), PrG(3), LLtmp(8),&
   PrXZ(3,3,3,8), PrB(3), PrGx(3), PrPB(3,2), PrV(3), PrPAX(3,2), ALR
logical :: cat(8), AncOK, Aselfed
 
if (AgeDiff(A, B) <= 0) then  ! B younger than A
  LL = impossible
  return
endif            

LL = missing  
AncOK = .TRUE.
call ChkAncest(B, Sex(B), A, Sex(A), AncOK)
if (.not. AncOK)  then
  LL = impossible
else 
  call ChkAncest(B, Sex(B), Parent(A,k), k, AncOK)
  if (.not. AncOK)  LL = impossible
endif
if (LL/=missing) return

if (Parent(A,k)<0) then
  if (any(SibID(:,-parent(A,k),k) == B)) then
    LL = impossible
    return
  endif
  if (ns(-Parent(A,k),k) > 1) then
    call AddGP(B, -Parent(A,k), k, LLtmp(1))
    if (LLtmp(1) < 0) then
      LL = LLtmp(1) - CLL(-Parent(A,k), k) + Lind(A)
    else
      LL = LLtmp(1)
    endif
    return
  endif
endif

curGP = 0 
do z=1,2
  curGP(:,z) = getPar(Parent(A,z), z)
enddo

if (Sex(B)<3) then
  m = Sex(B)
else if (curGP(1,k) == 0) then
  m=1    ! doesn't really matter.
else
  m = 2
endif

if (ANY(curGP(:,k) == B)) then
  LL = AlreadyAss 
else if (focal == 4) then
  if (curGP(m,k) /=0)  LL = impossible
else
  if (curGP(m,k) > 0)  LL = impossible
endif
if (LL/=missing) return

if (Complx==0 .and. curGP(3-m,k)==0 .and. Parent(A,3-m)==0) then
  curGP(3-m,k) = Mate(B)
endif

LLtmp = missing
if (Parent(A,k)/=0) then
  call CalcAgeLR(Parent(A,k),k, B, Sex(B), m, 1, .FALSE., ALR) 
  if (ALR == impossible) then   ! B too young/old to be parent of Parent(A,k)
    LL = impossible
  else if (Parent(A,k)>0) then
    call PairPO(Parent(A,k), B, m, focal, LLtmp(1))
    if (LLtmp(1) > 0) then    ! impossible
      LL = impossible
    else 
      call CalcU(Parent(A,k), k, B,k, LLtmp(2))
      if (LLtmp(1) - LLtmp(2) < TA) then
        LL = impossible  
      endif
    endif
  endif
endif
if (LL/=missing) return

Aselfed = .FALSE.
if (hermaphrodites/=0) then
  if (all(Parent(A,:) > 0)) then
    if (Parent(A,k) == Parent(A,3-k))  Aselfed = .TRUE.
  else if (all(Parent(A,:) < 0)) then
    if (DumClone(-Parent(A,k), k) == -Parent(A,3-k))  Aselfed = .TRUE. 
  endif
endif

! cat: 1: non-inbred, 2: double GP, 3: GP+HS, 4: GP+PO, 5 & 6: P-O mating
ParAB = 0   ! for cat(3)
if (complx < 2 .or. Aselfed) then
  cat = .FALSE.
  cat(1) = .TRUE.
else
  cat = .TRUE.   
  if (Parent(B,3-k)==Parent(A,3-k) .and. Parent(A,3-k)/=0) then
    cat(1) = .FALSE.
  endif
  if ((focal==1 .or. focal==4) .and. (all(parent(A,:)==0) .or. all(parent(B,:)==0))) then
  ! if no parents assigned yet, double GP indistinguishable from PO
    cat(2) = .FALSE.    ! cat(2:3) ?
  endif
  if (Parent(A,3-k)/=0) then
    if (curGP(m, 3-k) == B) then
      cat = .FALSE.
      cat(2) = .TRUE.    ! assignment would make it double GP
    else 
      call ChkAncest(B, Sex(B), Parent(A,3-k), 3-k, AncOK)
      if (.not. AncOK) then
        cat(2) = .FALSE.
      else
        call CalcAgeLR(Parent(A,3-k), 3-k, B,k, 0, 1, .TRUE., ALR)
        if (ALR == impossible .or. ALR < 3.0*TF)  cat(2) = .FALSE.
      endif
    endif
  endif
  if (focal==3) then  
    cat(3) = .FALSE.
  else if (Parent(B,3-k)==0 .or. Parent(A,3-k)==0 .or. & 
   Parent(B,3-k)==Parent(A,3-k)) then
    if (Parent(A,3-k)/=0 .and. Parent(B,3-k)==0) then
      call ChkAncest(Parent(A,3-k), 3-k, B, 0, AncOK)
      if (.not. AncOK) then
        cat(3) = .FALSE.
      else
        call CalcAgeLR(B,k, Parent(A,3-k), 3-k, 0, 1, .TRUE., ALR)
        if (ALR == impossible .or. ALR < 3.0*TF)  cat(3) = .FALSE.
      endif
    endif
    if (Parent(A,3-k)==0) then
      ParAB = Parent(B,3-k)
    else
      ParAB = Parent(A,3-k)
    endif
  else
    cat(3) = .FALSE.
  endif
  if (cat(3)) then
    call CalcAgeLR(A,0, B,0, 3-k, 3, .TRUE., ALR) 
    if (ALR == impossible .or. ALR < 3.0*TF)  cat(3) = .FALSE.
  endif
  
  if (Parent(A,3-k)==0 .and. Sex(B)==3-k .and. focal/=1 .and. focal/=7) then 
    cat(4) = .FALSE. !.TRUE.   ! not implemented. is already implemented in PairPO?
  else
    cat(4) = .FALSE.
  endif
  if (any(Parent(A,:) /= 0)) then
    cat(5:6) = .FALSE.   ! possible but not necessary to consider (?)
  endif
endif

cat(7:8) = .FALSE.
if (Sex(B)==4 .and. focal/=7) then
!  if (Parent(A,3-k)/=B)  cat(7) = .TRUE.   ! in-between parent is selfed; identical to LL(PO)
  if (all(parent(A,:)==0) .or. Aselfed)  cat(8) = .TRUE.   ! A also selfed
endif
 
PrL = 0D0
do l=1,nSnp
  call ParProb(l, curGP(3-m,k), 3-m, 0,0, PrG) 
  call ParProb(l, Parent(A,k), k, A, -4, PrPA(:, k)) 
  if (Aselfed) then
    PrPA(:,3-k) = 1D0
  else if (Parent(A,3-k)==Parent(B,3-k) .and. Parent(A,3-k)/=0) then
    call ParProb(l, Parent(A,3-k), 3-k, A, B, PrPA(:, 3-k))
  else
    call ParProb(l, Parent(A,3-k), 3-k, A, 0, PrPA(:, 3-k))
  endif
  call ParProb(l, B, 0, 0, 0, PrB)
  if (cat(3)) then
    call ParProb(l, Parent(B,k), k, B, 0, PrPB(:,k))
    call ParProb(l, ParAB, 3-k, A, B, PrPB(:,3-k))
  endif
  if (cat(2) .or. cat(6)) then
    call ParProb(l, Parent(A,3-k), 3-k, A, -4, PrPAX(:, 3-k))                                                         
    if (Parent(A,3-k) > 0) then
      call ParProb(l, curGP(3-k,3-m), 3-m, Parent(A,3-k), 0, PrGx)
    else if (Parent(A,3-k) < 0) then
      call ParProb(l, curGP(3-k,3-m), 3-m, 0, 0, PrGx)   
    else
      PrGx = AHWE(:,l)
    endif
  endif

  PrXZ = 0D0
  do x=1,3  ! PA(k)     
    do y=1,3  ! PA(3-k)
      if (Aselfed .and. x/=y)  cycle
      do z=1,3  !  PrG(3-m)
        PrXZ(x,y,z,:) = OKA2P(Genos(l,A),x,y) * PrPA(x,k) * PrG(z)
        if (cat(1)) then
          PrXZ(x,y,z,1) = PrXZ(x,y,z,1) * PrPA(y,3-k) *&
           SUM(AKA2P(x, :, z) * PrB)  !non-inbred
        endif
        if (cat(2)) then   !inbreeding loop; B double gp
          do v=1,3
            PrV(v) = SUM(AKA2P(y,v,:) * PrGx) * AKA2P(x,v,z) * PrB(v)
          enddo
          PrXZ(x,y,z,2) =PrXZ(x,y,z,2) * PrPAX(y,3-k) * SUM(PrV)
        endif
        if (cat(3)) then  !B GP and HS of A
          do v=1,3
            PrV(v) = AKA2P(x, v,z) * SUM(AKA2P(v,y,:)*PrPB(:,k))*PrPB(y,3-k) * OcA(v,Genos(l,B))
          enddo
          PrXZ(x,y,z,3) =  PrXZ(x,y,z,3) * SUM(PrV)
        endif
        if(cat(4)) then
        ! TODO?
        endif
        if (cat(5)) then  ! only when Parent(A,:) == 0
          if (x /= z) then
            PrXZ(x,y,z,5) = 0D0
          else
            PrXZ(x,y,z,5) = SUM(OKA2P(Genos(l,A),x,y) * AKA2P(x,y,:) * PrB * AHWE(y,l))
          endif
        endif
        if (cat(6)) then  ! only when Parent(A,:) == 0
          do v=1,3
            PrV(v) = SUM(AKA2P(y,x,:) * PrGx) * AKA2P(x,v,z) * PrB(v)
          enddo
          PrXZ(x,y,z,6) = PrXZ(x,y,z,6) * SUM(PrV)
        endif
        if (cat(7)) then  ! B double GP, A's parent is selfed
          PrXZ(x,y,z,7) = OKA2P(Genos(l,A),x,y) * PrPA(x,k) * PrPA(y,3-k) * &
            AKA2P(x,z,z) * PrB(z)
        endif
        if (cat(8)) then  ! as 7, A also selfed
          if (x/=y) then
            PrXZ(x,y,:,8) = 0D0
          else
            PrXZ(x,x,z,8) = OKA2P(Genos(l,A),x,x) * AKA2P(x,z,z) * PrB(z)
          endif
        endif
      enddo
    enddo
  enddo
  do i=1,8  ! inbred/non-inbred  
    if (cat(i))   PrL(l,i) = LOG10(SUM(PrXZ(:,:,:,i)))
  enddo
enddo

LLtmp = SUM(PrL, DIM=1)
WHERE(LLtmp((/1,2,5,6,7,8/)) <0)  LLtmp((/1,2,5,6,7,8/)) = LLtmp((/1,2,5,6,7,8/)) + Lind(B)
if (Parent(A,k)>0) then
  WHERE(LLtmp <0)  LLtmp = LLtmp - Lind(Parent(A,k))
endif
if (Parent(A,3-k)>0 .and. cat(2)) then
  LLtmp(2) = LLtmp(2) - Lind(Parent(A,3-k))
endif

LL = MaxLL(LLtmp)
if (LL >= 0) then
  LL = impossible
endif

! if (A==1764 .and. B==92) then
  ! write (*, '("PairGP : ", 3i6, ", ", i3, 3i6, 8f8.1)')  A, Parent(A,:), k, B, Parent(B,:), LLtmp
! endif

end subroutine PairGP

! #####################################################################

subroutine LRGG(A,k, B,kB, LR)
use Global
implicit none

integer, intent(IN) :: A,k,B,kB
double precision, intent(OUT) :: LR
integer :: x, y, l
double precision :: PrXY(3,3,2), PrL(nSnp), PrPA(3), PrB(3)

PrL = 0D0
do l=1,nSnp
  call ParProb(l, Parent(A,3-k), 3-k, A, 0, PrPA)
  call ParProb(l, B, kB, 0, 0, PrB)
  PrXY = 1D0
  do x=1,3  ! PA(k)
    do y=1,3  ! B
      PrXY(x,y,:) = SUM(OKA2P(Genos(l,A),x,:) * PrPA) * PrB(y)
      PrXY(x,y,1) = PrXY(x,y,1) * AKAP(x,y,l)
      PrXY(x,y,2) = PrXY(x,y,2) * AHWE(x,l)
    enddo
  enddo
  PrL(l) = LOG10(SUM(PrXY(:,:,1))) - LOG10(SUM(PrXY(:,:,2)))
enddo
LR = SUM(PrL)

end subroutine LRGG

! #####################################################################

subroutine PairGGP(A, B, k, fcl, LL)   
! calculates LL that B is maternal(k=1) paternal(k=2), or double(k3) ggp
use Global
implicit none

integer, intent(IN) :: A,B,k, fcl
double precision, intent(OUT) :: LL
integer :: l, x, y,z,w, m, n, AncA(2,mxA)
double precision :: PrL(nSnp,4), PrXY(3,3), PrXZ(3,3,3), PrXW(3,3,3,3,2), LLtmp(4),&
  PrPA(3),PrB(3), PrG(3,2), PrPAX(3,2), ALR
logical :: MaybeLoop(2), AncOK(2)
  
LL = missing
LLtmp = missing
AncOK = .TRUE.              
if (AgeDiff(A, B) <= 0) then  ! B younger than A
  LL = impossible
else if (B==Parent(A,k)) then
  LL = impossible
else
  call ChkAncest(B,k, A, 0, AncOK(1))
  if (.not. AncOK(1)) then
    LL = impossible
    return
  endif
  if (Parent(A,3-k)/=0 .and. Parent(A,3-k)==Parent(B,3-k)) then
    LL = NotImplemented   ! not impossible, just unlikely.
    return
  endif
  if (Parent(A,k)/=0) then
    call ChkAncest(B,k,Parent(A,k),k, AncOK(2))
    if (.not. AncOK(2))   LL = impossible
  endif
endif
if (LL==impossible .or. LL==NotImplemented) return

if (Parent(A,k)>0) then 
  if (ANY(Parent(Parent(A,k), :)/=0)) then
    LL = NotImplemented    ! should be picked up elsewere
  else
    call PairGP(Parent(A,k), B, k, 4, LLtmp(1))
    if (LLtmp(1) > 0) then    
      LL = LLtmp(1)
    endif
  endif
else if (Parent(A,k)<0) then
  if (ANY(GpID(:,-Parent(A,k),k)/=0)) LL = NotImplemented
endif
if (LL/=missing) return

MaybeLoop = .FALSE.
if (fcl/=4) then  ! double GGP indistinguishable from GP
  if (ALL(Parent(A,:)==0)) then
    MaybeLoop = .TRUE.
  else
    do m=1,2
      if (ALL(getPar(Parent(A,m),m)==0)) then
        call CalcAgeLR(Parent(A,m),m,B,Sex(B),k,4, .TRUE., ALR)
        if (ALR/=impossible .and. ALR>TF) then
          MaybeLoop(m) = .TRUE.
        endif
      endif
    enddo
  endif
endif     

AncA = 0       
call getAncest(A, k, AncA)

PrL = 0D0    
do l=1,nSnp
  call ParProb(l, B, 0, 0, 0, PrB)
  call ParProb(l, Parent(A,3-k), 3-k, A, 0, PrPA)  ! with GP con
  PrPAX = 1D0
  do n=1,2
    if (n/=k .and. .not. MaybeLoop(n)) cycle
    if (Parent(A,n)/=0) then
      call ParProb(l, Parent(A,n), n, A, -4, PrPAX(:,n))  !=OcA if >0
    endif
    if (ALL(MaybeLoop)) then
      PrG(:,n) = AHWE(:,l)
      do m=1,2
        if (AncA(m, n+2)/=0) then ! either GP ==0
          if (Parent(A,n)>0) then
            call ParProb(l,AncA(m, n+2), m,Parent(A,n),0, PrG(:,n))   
          else
            call ParProb(l,AncA(m, n+2), m,0,0, PrG(:,n))    
          endif
        endif
      enddo
    endif
  enddo
  
  PrXY = 0D0
  PrXZ = 0D0
  PrXW = 0D0
  do x=1,3  
    do y=1,3 
      PrXY(x,y) = SUM(OKA2P(Genos(l,A),x,:) * PrPA) * PrPAX(x,k) * &
       AKAP(x,y,l) * SUM(AKAP(y, :, l) * PrB)
      if (ANY(MaybeLoop)) then
        do z=1,3  !consider double GGP (k & 3-k, or 2x k)
          if (ALL(MaybeLoop)) then
            do w=1,3
              PrXW(x,y,z,w,:) = OKA2P(Genos(l,A),x,z) *PrPAX(x,k) *&
               PrPAX(z,3-k) * SUM(AKA2P(z,w,:)* PrG(:,3-k))* &
               SUM(AKA2P(x,y,:)*PrG(:,k))
              PrXW(x,y,z,w,1) = PrXW(x,y,z,w,1) * SUM(AKAP(y,:,l) * AKAP(w,:,l) * PrB)
              PrXW(x,y,z,w,2) = PrXW(x,y,z,w,2) * SUM(AKAP(y,:,l) * AKAP(w,:,l) * AHWE(:,l))
            enddo
          endif
          if (MaybeLoop(k)) then
            PrXZ(x,y,z) = SUM(OKA2P(Genos(l,A),x,:) *PrPA) * PrPAX(x,k) *&
             AKA2P(x,y,z) * SUM(AKAP(y,:,l) * AKAP(z,:,l) * PrB)
          endif
        enddo
      endif
    enddo
  enddo
  PrL(l,1) = LOG10(SUM(PrXY))
  if (ALL(MaybeLoop)) then
    PrL(l,2) = LOG10(SUM(PrXW(:,:,:,:,1)))
    PrL(l,4) = LOG10(SUM(PrXW(:,:,:,:,2)))
  endif
  if (MaybeLoop(k))    PrL(l,3) = LOG10(SUM(PrXZ))
enddo

LLtmp = SUM(PrL, DIM=1)   
do n=1,2
  if (Parent(A,n)>0) then
    if (n==k) then
      LLtmp(1) = LLtmp(1) - Lind(Parent(A,n))
      if (MaybeLoop(k)) then
        LLtmp(3) = LLtmp(3) - Lind(Parent(A,n))
      endif
    endif
    if (ALL(MaybeLoop)) then
      LLtmp(2) = LLtmp(2) - Lind(Parent(A,n))
      LLtmp(4) = LLtmp(4) - Lind(Parent(A,n))
    endif
  endif
enddo
if (ALL(MaybeLoop)) then  ! compare to inbreeding loop w/o B at the top. 
  if (LLtmp(4) > Lind(A)) then
    LLtmp(2) = LLtmp(2) - LLtmp(4) + Lind(A)
  endif
endif

LL = MaxLL(LLtmp(1:3)) + Lind(B)

! if (A==345 .and. B==27) then
  ! write(*,'("PairGGP: ", 3i6, " + ", i6, 2i3, 2l4, 4f8.2, ", ", f8.2)') &
    ! A, Parent(A,:), B, k, fcl, MaybeLoop, LLtmp, LL
  ! write(*,'("SUM PrL: ", 4f8.2)')  SUM(PrL, DIM=1)
! endif

end subroutine PairGGP

! #####################################################################

 subroutine PairGA(A, B, k, hf, LL)   ! B FS/HS of GP
use Global
implicit none

integer, intent(IN) :: A,B,k, hf
double precision, intent(OUT) :: LL
integer :: l, x, y,v,w, m, n
double precision :: PrL(nSnp), PrX(3,3,3,3), PrPA(3), PrGG(3, 2) 
logical :: AncOK

LL = missing
if (AgeDiff(A, B) <= 0) then  ! B younger than A
  LL = impossible
else if (ANY(Parent(A,:)==B)) then
  LL = NotImplemented
else
  AncOK = .TRUE.              
  call ChkAncest(B,0,A,0,AncOK)
  if (.not. AncOK) then
    LL = impossible
  else if (Parent(A,3-k)/=0 .and. Parent(A,3-k)==Parent(B,3-k)) then
    LL = NotImplemented   ! not impossible, just unlikely.
  endif
endif
if (LL/=missing) return

m = 3-k  ! most neutral, doesn't matter in most cases
if (Parent(A,k)/=0) then
  LL = NotImplemented
  return
endif

PrL = 0D0    
do l=1,nSnp
  PrX = 0D0
  call ParProb(l, Parent(A,3-k), 3-k, A, 0, PrPA) 
  do n=1,2
    call ParProb(l, Parent(B, n), n, B, -1, PrGG(:,n))
  enddo

  do v=1,3
    do w=1,3
      do x=1,3  
        do y=1,3
          PrX(x,y,v,w) = AKAP(x,y,l) * PrGG(v,1) * PrGG(w,2)
          if (hf==3) then
            PrX(x,y,v,w) = PrX(x,y,v,w) * AKA2P(y, v, w)
          else if (hf==1) then
            PrX(x,y,v,w) = PrX(x,y,v,w) * AKAP(y, v, l)
          else if (hf==2) then
            PrX(x,y,v,w) = PrX(x,y,v,w) * AKAP(y, w, l)
          endif
          PrX(x,y,v,w) = PrX(x,y,v,w) * SUM(OKA2P(Genos(l,A),x,:) * PrPA)
        enddo
      enddo
      PrX(:,:,v,w) = PrX(:,:,v,w) * OKA2P(Genos(l,B), v, w)
    enddo
  enddo
  PrL(l) = LOG10(SUM(PrX))
enddo

LL = SUM(PrL)

! if ((A==6024 .and. B==2096)) then
  ! write (*, '("PairGA : ", 3i6, ", ", i3, 3i6, 2f8.1)')  A, Parent(A,:), k, B, Parent(B,:), SUM(PrL(:,2)), SUM(PrL(:,1))
  ! print *, 'nFS B: ', nFSB
! endif

end subroutine PairGA

! #####################################################################

subroutine PairUA(A, B, kA, kB, LL)
! B half sib or full sib (kB=3) of parent kA of A?
use Global
use CalcLik
use qsort_c_module
implicit none

integer, intent(IN) :: A,B,kA, kB  ! kB=3 : full sibs
double precision, intent(OUT) :: LL
integer :: AA(maxSibSize), PA, nA, GA(2), BB(maxSibSize), nB, PB(2), BBx(maxSibSize, 2), &
  nBx(2), Mates(maxSibSize, 2), GG(2), AncG(2, 2,mxA), AncA(2,mxA), AB(2*maxSibSize), &
  catA(maxSibSize), catG(2), catB(maxSibSize), GGP, GGG(2), doneB(maxSibSize), &
  Bj, w, l, x, g, y, z, i, r,u,j,e,Ei,m, PAx, DoQuickB
double precision :: PrL(nSnp), PrG(3,2), PrXYZ(3,3,3), PrPA(3), PrA(3),&! PrLx(nSnp, 2), & 
  PrPB(3), PrGA(3), PrAB(3,3,3,2), PrE(3), PrH(3), PrGG(3), PrEW(3,3), PrW(3), PrXY(3,3), PrFS(3,3)                                         
integer, allocatable, dimension(:) :: UseEE, MateABpar, TypeEE
double precision, allocatable, dimension(:,:,:) :: PrEE
logical :: PAselfed, SIMPL, AncOK, MateLoop(maxSibSize,2), DoAZ

AA = 0
if (A>0) then  
  PA = Parent(A, kA)
  if (PA<0) then
    nA = ns(-PA,kA)
    AA(1:nA) = SibID(1:nA, -PA, kA)
  else
    nA = 1
    AA(1) = A
  endif
else
  nA = nS(-A, kA)
  AA(1:nA) = SibID(1:nA, -A, kA)
  PA = A
endif
GA = getPar(PA, kA)
  
if (nA==0) then
  LL = NotImplemented
  return
endif

BB = 0
nB = 0
PB = 0
BBx = 0
nBx = 0
Mates = 0
if (B > 0) then
  nB = 1
  BB(1) = B  ! for cat checks
  PB = Parent(B,:)
  do m=1,2
    if (kB<3 .and. m/=kB)  cycle
    if (Parent(B, m) >=0) then
      nBx(m) = 1
      BBx(1, m) = B
    else 
      nBx(m) = nS(-Parent(B, m), m)
      BBx(1:nBx(m), m) = SibID(1:nBx(m), -Parent(B, m), m)  ! half sibs
    endif
    do j=1,nBx(m)
      Mates(j,m) = Parent(BBx(j, m), 3-m)
    enddo
  enddo
else !if (B < 0) then
  nB = nS(-B, kB)
  BB(1:nB) = SibID(1:nB, -B, kB)
  PB(kB) = B
  PB(3-kB) = 0
  do j=1,nB
    Mates(j,kB) = Parent(BB(j), 3-kB)
  enddo
endif

if (PB(kA) == PA .and. PA/=0 .and. (kB==kA .or. kB==3)) then
  LL = impossible
  return
endif

LL = missing
if (kB < 3) then
  if (B < 0 .and. GA(kB) == B) then
    LL = AlreadyAss
  else if (B > 0 .and. GA(kB) == PB(kB) .and. PB(kB)/=0)  then
    LL = AlreadyAss
  else if (GA(3-kB)==PB(3-kB) .and. PB(3-kB)/=0) then ! B>0; FA not HA
    LL = impossible
  endif
else if (GA(1)==PB(1) .and. GA(2)==PB(2) .and. GA(1)/=0 .and. &
  GA(2)/=0) then  ! kB==3
  LL = AlreadyAss
endif
if (LL /= missing) return

GG = 0  ! parent of B, GP of A
AncG = 0
do x=1,2
  if (x/=kB .and. kB/=3) cycle
  if (GA(x)==0) then
    GG(x) = PB(x)
  else if (GA(x)/=PB(x) .and. PB(x)/=0) then
    LL = impossible
  else
    GG(x) = GA(x)
  endif
  if (ANY(AA(1:nA)==GG(x))) then
    LL = impossible
  endif
  do j=1,nBx(x)
    if (PA<0 .and. Parent(BBx(j,x), 3-x)<0) then
      if (GpID(kA, -Parent(BBx(j,x), 3-x), 3-x) == PA) then
        LL = NotImplemented
      endif
    endif
  enddo
enddo
if (LL /= missing) return

do x=1,2
  call GetAncest(GG(x), x, AncG(x, :, :))
  if (PA/=0 .and. (x==kB .or. kB==3)) then
     if (ANY(AncG(x,kA,:) == PA)) then
      LL = impossible
      return
    endif
  endif
enddo

AncA = 0               
if (A > 0) then
  if (ANY(AncG == A)) then
    LL = impossible
  endif
else if (A < 0) then
  if (ANY(AncG(:, kA, 2:mxA) == A)) then
    LL = impossible
  endif
endif
if (B > 0) then
  if (ANY(AncG == B)) then
    LL = impossible
  endif
else if (B < 0) then
  if (ANY(AncG(:, kB, 3:mxA) == B)) then
    LL = impossible
  endif
endif
if (kB<3 .and. any(GA < 0)) then
  call GetAncest(A, kA, AncA)   
  if (B<0) then
    if (ANY(AncA(kB,3:mxA)==B))  LL = NotImplemented  ! B is GGP; 
  else if (B>0) then
    if (ANY(AncA(:,3:mxA)==B))  LL = NotImplemented
  endif
endif
if (LL /= missing) return

do x=2,mxA
  do y=1,2
    do g=1,2
      if (AncG(g,y,x) > 0) then
        if (A > 0) then
          if (AgeDiff(A, AncG(g,y,x)) < 0)  LL = impossible  ! A older than putative ancestor
        else if (A<0) then
          do i=1, ns(-A,kA)
            if (AgeDiff(SibID(i,-A,kA),AncG(g,y,x)) < 0)  LL = impossible
          enddo       
        endif
        if (x==2) cycle 
        if (B > 0) then
          if (AgeDiff(B, AncG(g,y,x)) < 0)  LL = impossible  
        else if (B<0) then
          do i=1, ns(-B,kB)
            if(AgeDiff(SibID(i,-B,kB),AncG(g,y,x)) < 0)  LL = impossible
          enddo       
        endif
      endif
    enddo
  enddo
enddo
if (LL /= missing) return

PAselfed = .FALSE.
if (hermaphrodites/=0) then
  if (kB==3) then
    if (all(PB < 0)) then
      if (DumClone(-PB(1),1) == -PB(2))  PAselfed = .TRUE.
    endif
  else if (GA(3-kB) < 0 .and. PB(kB)<0) then
    if (DumClone(-GA(3-kB),3-kB) == -PB(kB))  PAselfed = .TRUE.
  endif
endif

!==============================================

allocate(UseEE(nA+nB))
allocate(TypeEE(nA+nB))
allocate(PrEE(3,3, nA+nB))
allocate(MateABpar(nA+nB))
UseEE = 0
TypeEE = 0
PrEE = 0D0          
MateABpar = 0
AB = 0

if (kA==kB) then
  AB(1:nA) = AA(1:nA)
  AB((nA+1):(nA+nB)) = BB(1:nB)
  call FindEE(AB(1:(nA+nB)), nA, nB, kA, UseEE, MateABpar)  ! may reorder AA, BB
  AA(1:nA) = AB(1:nA)
  BB(1:nB) = AB((nA+1):(nA+nB))
  TypeEE = 3-kA
else if (kA/=kB  .and. kB/=3) then
  call FindEE(AA(1:nA), nA, 0, kA, UseEE(1:nA), MateABpar(1:nA)) 
  call FindEE(BB(1:nB), nB, 0, kB, UseEE((nA+1):(nA+nB)), MateABpar((nA+1):(nA+nB)))
  do i=1, nB
    if (UseEE(nA+i)/=0) then
      UseEE(nA+i) = nA + UseEE(nA+i)
    endif
  enddo
  TypeEE(1:nA) = 3-kA
  TypeEE((nA+1):(nA+nB)) = 3-kB
endif

! TODO: kB==3, BBx i.o. BB
!==============================================

PrL = 0D0
PAx = 0
if (nA>0)  PAx = Parent(AA(1),3-kA)

if (A>0 .and. B>0 .and. PA>=0 .and. PAx>=0 .and. ALL(PB>=0) .and. all(GG>=0)) then  ! quicker.
  do l=1, nSnp
    call ParProb(l, PAx, 3-kA, A, 0, PrPA)
    if (kB == 3) then
      do g=1,2
        call ParProb(l, GG(g), g, 0, 0, PrG(:,g))  ! >=0
      enddo        
    else
      call ParProb(l, GG(kB), kB, 0, 0, PrG(:,kB))  ! >=0
      if (PA>0) then
        call ParProb(l, GA(3-kB), 3-kB, PA, 0, PrGA)
      else
        PrGA = AHWE(:,l)
      endif
      call ParProb(l, PB(3-kB), 3-kB, B, 0, PrPB)
    endif
    
    PrXYZ = 0D0
    do z=1,3
      do y=1,3
        do x=1,3
          if (kB == 3) then
            PrXYZ(x,y,z) = AKA2P(x,y,z)*PrG(y,1)*PrG(z,2)
          else
            PrXYZ(x,y,z) = AKA2P(x,y,z)*PrG(y,kB)*PrGA(z)   
          endif         
          PrXYZ(x,y,z) = PrXYZ(x,y,z) *SUM(OKA2P(Genos(l,A), x, :) * PrPA)
          if (kB==3) then
            PrXYZ(x,y,z) = PrXYZ(x,y,z) *OKA2P(Genos(l,B),y,z)
          else
            PrXYZ(x,y,z) = PrXYZ(x,y,z) *SUM(OKA2P(Genos(l,B),y,:) * PrPB)
          endif
          if (PA>0) then
            PrXYZ(x,y,z) = PrXYZ(x,y,z) * OcA(x,Genos(l,PA))
          endif
        enddo
      enddo
    enddo
    PrL(l) = LOG10(SUM(PrXYZ))
  enddo
  
  LL = SUM(PrL)
  if (PA>0) then
    LL = LL - Lind(PA)
  endif
  deallocate(UseEE)
  deallocate(PrEE)
  deallocate(MateABpar)
  deallocate(TypeEE)                 
  return
endif

!==============================================
if (A>0 .and. B<0 .and. kB/=3) then
  SIMPL = .TRUE.
  AncOK = .TRUE.
  if(ANY(Parent(A,:)<0) .or. Parent(A,kA)/=0) then   ! ANY(Parent(A,:)/=0)
    SIMPL = .FALSE.
  endif
  do j=1,nB
    if (nFS(BB(j))==0)  cycle
    Bj = BB(j)
    call ChkAncest(Bj, 0, A, Sex(A), AncOK)
    if (.not. AncOK)  SIMPL = .FALSE.
    call ChkAncest(A, Sex(A), Bj, 0, AncOK)
    if (.not. AncOK)  SIMPL = .FALSE.
    if (.not. SIMPL)  exit
  enddo
  if (SIMPL) then
    call ChkDoQuick(-B, kB, DoQuickB)
    do l=1, nSnp
      call ParProb(l,Parent(A,3-kA),3-kA,0,0,PrE)
      if (DoQuickB == -2) then  ! all are FS; ns(s,k) = ns(..,3-k)
        call ParProb(l, Parent(Bj,3-kB),3-kB,-1,0,PrH)
      endif
      PrFS = FSLik(l,Bj)
      do y=1,3
        do x=1,3
          if (DoQuickB == -2) then  
            do z=1,3
              PrW(z) = PrFS(z,y) * XPr(2,y,l,-B,kB) * PrH(z) * AKAP(x,y,l) * &
               SUM(OKA2P(Genos(l,A),x,:) * PrE)
            enddo
            PrXY(x,y) = SUM(PrW)         
          else
            PrXY(x,y) = XPr(3,y,l, -B,kB) * AKAP(x,y,l) * SUM(OKA2P(Genos(l,A),x,:) * PrE)
          endif
        enddo
      enddo
      PrL(l) = LOG10(SUM(PrXY))
    enddo
    LL = SUM(PrL)
    deallocate(UseEE)
    deallocate(PrEE)
    deallocate(MateABpar)
    deallocate(TypeEE)
    return
  endif
endif 

!==============================================

if (A<0 .and. B>0 .and. all(UseEE==0)) then
  SIMPL = .TRUE.
  AncOK = .TRUE.
  if (any(Parent(B,:) < 0)) then
    SIMPL = .FALSE.
  else if (kB < 3) then
    if (Parent(B,kB)/=0)  SIMPL = .FALSE.
  else
    if (any(Parent(B,:)/=0))  SIMPL = .FALSE.
  endif
  if (any(GpID(:,-A,kA)/=0))  SIMPL = .FALSE.
  if (any(Parent(AA(1:nA), 3-kA) < 0))  SIMPL = .FALSE.
  if (SIMPL) then
    do j=1,nA
      call ChkAncest(AA(j), 0, B, Sex(B), AncOK)
      if (.not. AncOK)  SIMPL = .FALSE.
      call ChkAncest(B, Sex(B), AA(j), 0, AncOK)
      if (.not. AncOK)  SIMPL = .FALSE.
    enddo
  endif
  if (SIMPL) then
    do l=1, nSnp
      if (kB < 3)  call ParProb(l,Parent(B,3-kB),3-kB,0,0,PrE)
      do x=1,3
        do y=1,3
          PrXY(x,y) = XPr(1,x,l, -A,kA) * AHWE(y,l)
          if (kB<3) then
            PrXY(x,y) = PrXY(x,y) * AKAP(x,y,l) * SUM(OKA2P(Genos(l,B),y,:) * PrE)
          else
            PrXY(x,y) = PrXY(x,y) * SUM(AKA2P(x,y,:) * OKA2P(Genos(l,B),y,:) * AHWE(:,l))
          endif
        enddo
      enddo
      PrL(l) = LOG10(SUM(PrXY))
    enddo
    LL = SUM(PrL)
    deallocate(UseEE)
    deallocate(PrEE)
    deallocate(MateABpar)
    deallocate(TypeEE) 
    return
  endif
endif 

!============================================

catA=0
catG = 0
catB = 0
GGP = 0
GGG = 0       
do i = 1, nA
  if (Parent(AA(i),3-kA)==0) cycle
  if (hermaphrodites/=0 .and. PA<0) then
    if (Parent(AA(i), 3-kA) == -DumClone(-PA,kA)) then
      catA(i) = 12
      cycle
    endif
  endif
  if (kA/=kB .and. Parent(AA(i), 3-kA) == GG(3-kA)) then  !incl. kB=3
    catA(i) = 1  
  else if (kA==kB .and. Parent(AA(i), 3-kA)==GA(3-kA)) then
    catA(i) = 2
    UseEE(i) = 0
  else 
    if (Parent(AA(i), 3-kA)<0) then
      if (GpID(kA,-Parent(AA(i), 3-kA),3-kA) == PA .and. PA/=0) then
        catA(i) = 7  ! Ai inbred
      endif
    endif
    if (kA==kB) then
      do j=1, MAX(nB, nBx(kB))
        if (B>0) then
          Bj = BBx(j,kB)
        else
          Bj = BB(j)
        endif
        if (AA(i) == Bj) cycle
        if (Parent(AA(i), 3-kA) == Parent(Bj, 3-kA) .and. Parent(Bj, 3-kA) < 0) then 
          if (B<0)  catA(i) = 3   
        else if (Parent(AA(i), 3-kA) == Bj) then
          catA(i) = -j
        endif
      enddo
    endif
  endif
  do g=1,2
    if (kB/=g .and. kB/=3) cycle
    if (Parent(AA(i), 3-kA) < 0) then
      if (GpID(g,-Parent(AA(i), 3-kA),3-kA) == GG(g) .and. GG(g)/=0) then
        if (g==kB .or. (kB==3 .and. g==3-kA)) then
          if (ALL(GpID(:,-Parent(AA(i), 3-kA),3-kA) == GG) .and. ALL(GG/=0)) then
            catA(i) = 10
          else
            catA(i) = 8  ! via y
          endif
        else
          catA(i) = 9  ! via z
        endif
      endif
    endif
    GGG = getPar(GG(g), g)
    if (Parent(AA(i),3-kA) == GGG(3-kA)) then
      catA(i) = 5  ! TODO? 4+g when kB==3
      catG(g) = 2
      GGP = GGG(kA)
      UseEE(i) = 0
    ! TODO: parent(parent(A,3-kA),kB) == B), B<0
    endif
  enddo
enddo
if (kB/=3) then   ! TODO: for kB==3
  do j=1, MAX(nB, nBx(1), nBx(2))
    if (B>0) then
      Bj = BBx(j,kB)
    else
      Bj = BB(j)
    endif
    if (Parent(Bj,3-kB)==0) cycle
    if (hermaphrodites/=0 .and. PB(kB) < 0) then
      if (Parent(Bj, 3-kB) == -DumClone(-PB(kB),kB) .and. DumClone(-PB(kB),kB)/=0) then
        catB(j) = 12
        cycle
      endif
    endif
    if (Parent(Bj, 3-kB) == GA(3-kB) .and. GA(3-kB)/=0) then   
      catB(j) = 2
      if (B<0) UseEE(nA+j) = 0
      if (B>0) UseEE(nA+1) = 0
    else if (B<0 .and. Parent(Bj,3-kB)<0) then
      if (GpID(kB, -Parent(Bj,3-kB),3-kB) == GG(kB) .and. GG(kB)/=0) then
        catB(j) = 7
      else if (GpID(kA, -Parent(Bj,3-kB),3-kB) == PA .and. PA/=0) then
        catB(j) = 8
      endif
    endif
    do g=1,2
      GGG = getPar(GG(g), g)
      if (Parent(Bj,3-kB) == GGG(3-kB)) then
        catB(j) = 5  
        if(catG(g)==0)  catG(g) = 3
        GGP = GGG(kB)
        if (B<0)  UseEE(nA+j) = 0  ! ??
        if (B>0) UseEE(nA+1) = 0
      endif
    enddo
    if (ANY(catA == 8) .and. kB/=3 .and. catB(j)==0) then
      do i=1,nA
        if (PA<0 .and. NFS(AA(i))==0) cycle
        if (Parent(AA(i), 3-kA)>=0) cycle
        if (GpID(3-kB,-Parent(AA(i), 3-kA),3-kA) == Parent(Bj,3-kB)) then
          catB(j) = -i
        endif
      enddo
    endif       
  enddo
endif

if (Complx<2 .and. (ALL(catA(1:nA)==1) .or. ALL(catA(1:nA)==3))) then  
  LL = NotImplemented  ! explicit consideration of close inbreeding
  deallocate(UseEE)
  deallocate(PrEE)
  deallocate(MateABpar)
  deallocate(TypeEE) 
  return
endif

if (kB/=3) then
  GGG = getPar(GG(kB), kB)
  if (GGG(3-kB) == GA(3-kB) .and. GA(3-kB)/=0) then
    catG(kB) = 1
    GGP = GGG(kB)
  endif
endif

MateLoop = .FALSE.
if (B>0) then    ! TODO: B<0   superseded by UseEE -- TODO convert mateloop > UseEE
  do m=1,2
    if (kB/=3 .and. m/=kB)  cycle
    do j=1, nBx(m)
      Bj = BBx(j, m)
      if (nFS(Bj)==0) cycle  
      if (kB==3 .and. Parent(Bj,1)==GG(1) .and.  Parent(Bj,2)==GG(2))  cycle
      if (Parent(Bj,m)<0 .and. Parent(Bj,3-m)<0) then
        do g=1, nS(-Parent(Bj, 3-m),3-m)
          Ei = SibID(g,-Parent(Bj,3-m),3-m)
          if (Parent(Ei,m)>=0 .or. Parent(Ei,m)==Parent(Bj,m)) cycle
          if (ANY(Mates(:,3-m) == Parent(Ei, m))) then
            MateLoop(j,m) = .TRUE.
          endif
        enddo
      endif
    enddo
  enddo
endif

!==============================================

PrL = 0D0
!PrLx = 0D0          
DoneB = 0
SIMPL = ALL(catA==0) .and. ALL(catG==0) .and. ALL(catB==0) .and. .not. ANY(MateLoop) .and. &
  ALL(GG >=0) .and. all(Parent(AA(1:nA),3-kA) >=0) .and. .not. PAselfed .and. &
  .not. (all(PB < 0) .and. kB==3)  
if (SIMPL .and. ANY(UseEE /= 0)) then
  if (ALL(PB >= 0)) then
    SIMPL = .TRUE.
  else
    SIMPL = .FALSE.
  endif
endif

DoAZ = any(catA==2 .or. catA==9 .or. catA==10) .or. catG(kA)==2
if (nA==1 .and. catA(1)==8 .and. PA==0)  catA(1) = 0   ! can't find the bug...

do l=1,nSnp
  do g=1,2
    if (g/=kB .and. kB/=3) cycle
    if (catG(g)==0) then
      if (SIMPL .and. B>0) then
        call ParProb(l, GG(g), g, B, -1, PrG(:,g)) 
      else
        call ParProb(l, GG(g), g, -1, 0, PrG(:,g))
      endif
    else
      if (GG(g) > 0) then 
        PrG(:,g) = OcA(:,Genos(l,GG(g)))
        if (catG(g)==1) then
          call ParProb(l, GGP, kB, GG(g), 0, PrGG) 
        else if (catG(g)==2) then
          call ParProb(l, GGP, kA, GG(g), 0, PrGG)
        else if (catG(g)==3) then
          call ParProb(l, GGP, kB, GG(g), 0, PrGG)
        endif
      else if (GG(g) < 0) then
        PrG(:,g) = 1D0
        if (catG(g)==1) then
          call ParProb(l, GGP, kB, 0, 0, PrGG)
        else if (catG(g)==2) then
          call ParProb(l, GGP, kA, 0, 0, PrGG)
        else if (catG(g)==3) then
          call ParProb(l, GGP, kB, 0, 0, PrGG)
        endif
      else 
          PrG(:,g) = AHWE(:,l)
      endif
    endif
  enddo
  if (kB/=3) then 
    if (ANY(catA==2) .or. ANY(catB==2)) then
      call ParProb(l, GA(3-kB), 3-kB, -1, 0, PrGA)
    else if (catG(kB)==1 .and. GG(kB)>0) then
      call ParProb(l, GA(3-kB), 3-kB, GG(kB), 0, PrGA)
    else if (PA>0) then
      call ParProb(l, GA(3-kB), 3-kB, PA, 0, PrGA)
    else
      call ParProb(l, GA(3-kB), 3-kB, 0, 0, PrGA)
    endif
    if (B>0) then
      call ParProb(l, PB(3-kB), 3-kB, B, 0, PrPB)   ! -1?
    endif
  endif

  ! === 
  
  PrA = 1D0
  if (SIMPL) then
    if (A < 0) then
      PrA = Xpr(1,:,l, -A,kA)
    else if (A>0) then
      call ParProb(l, Parent(A,3-kA), 3-kA, A, 0, PrPA)
      call Parprob(l, Parent(A,kA), kA, A, -4, PrH)  ! exclude both GPs & A
      do x=1,3 
        PrA(x) = SUM(OKA2P(Genos(l,A),x,:) * PrPA * PrH(x))
      enddo
    endif
  
    do x=1,3  ! PA, kA
      do y=1,3  ! PrG, kB
        do z=1,3  ! PrGA, 3-kB / PrG, 3-kB
          if (kB==3 .and. B>0) then  ! SA/PA FS of B; 
            PrXYZ(x,y,z) = PrA(x) * AKA2P(x,y,z)*PrG(y,3-kA) * PrG(z,kA) * &
              OKA2P(Genos(l,B), y, z)            
          else 
            PrXYZ(x,y,z) = PrA(x) * AKA2P(x,y,z) * PrGA(z)
            if (B>0) then
              PrXYZ(x,y,z) = PrXYZ(x,y,z) * PrG(y, kB) *SUM(OKA2P(Genos(l,B),y,:) * PrPB)
            else if (B<0) then
              PrXYZ(x,y,z) =PrXYZ(x,y,z) *XPr(3,y,l,-B,kB)
            endif
          endif     
        enddo
      enddo
    enddo
    PrL(l) = LOG10(SUM(PrXYZ))   

  else 

    PrAB = 0D0 
    do y=1,3  ! PrG, kB
      PrEE = 0D0
      do x=1,3  ! PA, kA
        if (PAselfed) then
          if (kB==3) then
            PrAB(x,y,y,:) = AKA2P(x,y,y) * PrG(y,kA)
          else
            PrAB(x,y,y,:) = AKA2P(x,y,y) * PrG(y,kB)
          endif
        else
          do z=1,3
            if (kB==3) then
              PrAB(x,y,z,:) = AKA2P(x,y,z) *PrG(z,kA) *PrG(y,3-kA)
            else if (catG(kB)==1) then
              PrAB(x,y,z,:) = AKA2P(x,y,z)*PrGA(z)*&
                SUM(AKA2P(y,z,:) * PrGG) * PrG(y, kB) 
            else if (catG(kB)==2) then
              PrAB(x,y,z,:) = AKA2P(x,y,z)*PrGA(z)
            else
              PrAB(x,y,z,:) = AKA2P(x,y,z) *PrGA(z) * PrG(y,kB) 
            endif
          enddo
        endif
        if (PA>0) then
          PrAB(x,y,:,:) = PrAB(x,y,:,:) * OcA(x,Genos(l,PA)) 
        endif
      enddo   
        
      do x=1,3
       DoneB = 0
       do z=1,3
        if (z>1 .and. .not. DoAZ)  cycle
        do r=1, nA
          if (PA<0 .and. NFS(AA(r))==0) cycle
          if (catA(r)==7) then
            call ParProb(l, GpID(3-kA,-Parent(AA(r),3-kA),3-kA), 3-kA, 0,0,PrH)
            do e=1,3
              PrE(e) = SUM(AKA2P(e,x,:) * PrH)
            enddo
          else if (catA(r)>=8 .and. catA(r)<=10) then
            if (kB < 3) then
              if (catA(r)==8) then
                g=kB
              else if (catA(r)==9) then
                g=3-kB
              endif
            else 
              if (catA(r)==8) then
                g=3-kA
              else if (catA(r)==9) then
                g=kA
              endif
            endif
            if (catA(r) < 10) then
              call ParProb(l, GpID(3-g,-Parent(AA(r),3-kA),3-kA), 3-g, 0,0,PrH)
            endif
            do e=1,3
              if (catA(r)==8) then
                PrE(e) = SUM(AKA2P(e,y,:) * PrH)
              else if (catA(r)==9) then
                PrE(e) = SUM(AKA2P(e,:,z) * PrH)
              else if (catA(r)==10) then
                PrE(e) = AKA2P(e,y,z)
              endif
            enddo  
            PrE = PrE/SUM(PrE)
          else if (catA(r)==42) then
            cycle ! do with B  (catB(j) = -i)
          else if (UseEE(r)/=0) then
            call ParProb(l, MateABpar(r), 3-TypeEE(r), 0,0,PrH)
            do e=1,3
              do u=1, 3
                PrW(u) = SUM(AKA2P(e,u,:) * PrEE(u,x,UseEE(r)) * PrH)
              enddo
              PrE(e) = SUM(PrW)
            enddo
            PrE = PrE/SUM(PrE)
          else if (catA(r) < 0) then
            if (kB==3) then
              PrH = PrG(:,kA)  
            else 
              if (B>0) then
                Bj = BBx(-catA(r),kA)
              else
                Bj = BB(-catA(r))
              endif
              call ParProb(l, Parent(Bj,3-kB), 3-kB, Bj,0,PrH)
            endif
            do e=1,3
              PrE(e) = SUM(AKA2P(e,y,:) * PrH)
            enddo
            PrE = PrE * OcA(:,Genos(l,Bj))
            PrE = PrE/SUM(PrE)
          else if (catA(r)==0 .or. (catA(r)>2 .and. catA(r)<7)) then
            call ParProb(l,Parent(AA(r),3-kA),3-kA,-1,0,PrE)
          else
            PrE = 1D0
          endif

          if (Parent(AA(r), 3-kA) < 0 .and. Parent(AA(r), 3-kA)/=GG(3-kA)) then
            if (A>0) then
              do i=1, MAX(nFS(AA(r)),1)
                if (FSID(i, AA(r))==A ) cycle      
                if (ANY(GG == FSID(i, AA(r)))) cycle
                PrE=PrE*OKA2P(Genos(l,FSID(i,AA(r))),x,:)  ! FS of A
              enddo
            endif
            
            do e=1,3
              if (catA(r)==1 .and. e/=y)  cycle
              if (catA(r)==2 .and. e/=z)  cycle
              if (catA(r)==12 .and. e/=x)  cycle
              do g=1, nS(-Parent(AA(r), 3-kA), 3-kA)
                Ei = SibID(g, -Parent(AA(r), 3-kA),3-kA)
                if (nFS(Ei) == 0) cycle 
                if (nFS(Ei) == 1 .and. (Ei==A .or. Ei==B))  cycle 
                if (Parent(Ei,kA)==PA .and. PA/=0) cycle
                if (kA==kB .and. Parent(Ei,kA)==GG(kA) .and. GG(kA)/=0)  cycle
                if (kB<3) then
                  if (Parent(Ei, kB)== PB(kB) .and. PB(kB)/=0) cycle
                endif                
                call ParProb(l,Parent(Ei,kA),kA,Ei,-1,PrH)
                if (catA(r)==5 .and. Parent(Ei,kA)==GGP .and. GGP/=0) then  ! FS of GG
                  do i=1, nFS(Ei)
                    if (any(GG == FSID(i, Ei))) cycle
                    if (FSID(i, Ei)==A .or. FSID(i, Ei)==B) cycle                                               
                    PrH = PrH * OKA2P(Genos(l,FSID(i,Ei)),:,e)
                  enddo
                  if (catG(kA)==2 .and. kB==3) then 
                    PrE(e) = PrE(e) * SUM(AKA2P(z,e,:) * PrH)  
                  else if (ANY(catG==2)) then
                    PrE(e) = PrE(e) * SUM(AKA2P(y,e,:) * PrH)
                  endif                
                else
                  do i=1, nFS(Ei)
                    if (FSID(i, Ei)==A .or. FSID(i, Ei)==B) cycle  
                    PrH=PrH*OKA2P(Genos(l,FSID(i,Ei)),:,e)
                  enddo
                  if (.not. all(PrH==1D0))  PrE(e) = PrE(e) * SUM(PrH)
                endif
              enddo  ! g
            enddo  ! e
            if (catA(r)==3 .and. B>0) then  ! Parent(AA(r), 3-kA) == Parent(Bj, 3-kA)
              do j=1,nBx(kA)
                if (Parent(BBx(j,kA), 3-kA) /= Parent(AA(r), 3-kA)) cycle
                do i=1, MAX(nFS(BBx(j,kA)),1)
                  if (FSID(i,BBx(j,kA))==B) cycle   
                  PrE = PrE * OKA2P(Genos(l,FSID(i,BBx(j,kA))), y, :)
                enddo
              enddo
            endif   
          endif
          if (catA(r)==5 .and. (Parent(AA(r), 3-kA)>0 .or. GGP==0)) then
            do e=1,3
              if (catG(kA)==2 .and. kB==3) then 
                PrE(e) = PrE(e) * SUM(AKA2P(z,e,:) * PrGG)
              else
                PrE(e) = PrE(e) * SUM(AKA2P(y,e,:) * PrGG)
              endif
            enddo
          endif
          
          if (catA(r)==1) then 
            if (DoAZ)        PrAB(x,y,z,1) = PrAB(x,y,z,1) * PrE(y)
            if (.not. DoAZ)  PrAB(x,y,:,1) = PrAB(x,y,:,1) * PrE(y)
          else if (catA(r)==2) then 
            PrAB(x,y,z,1) = PrAB(x,y,z,1) * PrE(z)
          else if (catA(r)==12) then 
            if (DoAZ)        PrAB(x,y,z,1) = PrAB(x,y,z,1) * PrE(x)
            if (.not. DoAZ)  PrAB(x,y,:,1) = PrAB(x,y,:,1) * PrE(x)
          else if (.not. all(PrE==1D0)) then  
            if (DoAZ)        PrAB(x,y,z,1) = PrAB(x,y,z,1) * SUM(PrE)
            if (.not. DoAZ)  PrAB(x,y,:,1) = PrAB(x,y,:,1) * SUM(PrE)
          endif       

          do i=1, MAX(nFS(AA(r)),1)
            if (A>0 .and. FSID(i, AA(r))/=A .and. .not. ANY(BB==FSID(i, AA(r)))) cycle
            PrE = PrE * OKA2P(Genos(l,FSID(i,AA(r))), x, :)  ! <- A
          enddo

          if ((catA(r)==3 .or. (catA(r)==5 .and. ANY(catB==5)) .or. &
           (catA(r)==2 .and. ANY(catB==2))) .and. kB<3) then 
            do j=1, MAX(nB, nBx(kB))
              if (B>0) then
                Bj = BBx(j,kB)
              else
                Bj = BB(j)
              endif
              if (Parent(Bj,3-kA) /= Parent(AA(r),3-kA)) cycle
              if (ANY(AA == Bj))  cycle
              PrE = PrE * OKA2P(Genos(l,Bj), y, :)
              DoneB(j) = 1
            enddo 
          endif
          
          if (catA(r)==1) then 
            if (DoAZ)        PrAB(x,y,z,2) = PrAB(x,y,z,2) * PrE(y)
            if (.not. DoAZ)  PrAB(x,y,:,2) = PrAB(x,y,:,2) * PrE(y)
          else if (catA(r)==2) then 
            PrAB(x,y,z,2) = PrAB(x,y,z,2) * PrE(z)
          else if (catA(r)==12) then 
            if (DoAZ)        PrAB(x,y,z,2) = PrAB(x,y,z,2) * PrE(x)
            if (.not. DoAZ)  PrAB(x,y,:,2) = PrAB(x,y,:,2) * PrE(x)
          else if (.not. all(PrE==1D0)) then    
            if (DoAZ)        PrAB(x,y,z,2) = PrAB(x,y,z,2) * SUM(PrE)
            if (.not. DoAZ)  PrAB(x,y,:,2) = PrAB(x,y,:,2) * SUM(PrE)
          endif
          PrEE(:,x,r) = PrE
        enddo  ! r
       enddo  ! z (if DoAZ)
      enddo  ! x
      
      do x=1,3
        if (x>1 .and. ALL(UseEE==0) .and. all(catB>=0))  cycle  !  .and. catB(j)/=5
      if (B<0) then            
        do j=1,nB
          if (nFS(BB(j))==0) cycle
          if (DoneB(j)==1) cycle
          if (kA/=kB .and. PA<0 .and. Parent(BB(j), 3-kB)==PA) cycle
!          DoneB(j) = 2  ! for output check only
          if (catB(j)==2 .or. catB(j)==12) then
            PrE = 1D0
          else if (catB(j)==7) then
            call ParProb(l, GpID(3-kB,-Parent(BB(j),3-kB),3-kB), 3-kB, 0,0,PrH)
            do e=1,3
              PrE(e) = SUM(AKA2P(e,y,:) * PrH)
            enddo
          else if (catB(j)==8) then
            call ParProb(l, GpID(3-kA,-Parent(BB(j),3-kB),3-kB), 3-kA, 0,0,PrH)
            do e=1,3
              PrE(e) = SUM(AKA2P(e,x,:) * PrH)
            enddo
          else if (UseEE(nA+j)/=0) then
            call ParProb(l, MateABpar(nA+j), 3-TypeEE(nA+j), 0,0,PrH)
            do e=1,3
              do u=1, 3
                PrW(u) = SUM(AKA2P(e,u,:) * PrEE(u,x,UseEE(nA+j)) * PrH)
              enddo
              PrE(e) = SUM(PrW)
            enddo
            PrE = PrE/SUM(PrE)
          else
            call ParProb(l, Parent(BB(j), 3-kB), 3-kB, -1,0,PrE)
          endif

          if (Parent(BB(j), 3-kB) < 0) then  
            do e=1,3
              do g=1, nS(-Parent(BB(j), 3-kB), 3-kB)
                Ei = SibID(g, -Parent(BB(j), 3-kB),3-kB)
                if (nFS(Ei) == 0) cycle
                if (Parent(Ei, kB) == PB(kB) .and. PB(kB)/=0) cycle  
                if (Parent(Ei, kA)== PA .and. PA/=0) cycle
                call ParProb(l,Parent(Ei,kB),kB,Ei,-1,PrH) 
                do i=1, nFS(Ei)
                  if (FSID(i, Ei)==A) cycle  
                  PrH = PrH * OKA2P(Genos(l,FSID(i,Ei)), :, e)
                enddo
                if (.not. all(PrH == 1D0))  PrE(e) = PrE(e) * SUM(PrH) 
              enddo 
            enddo                   
          endif
          
          if (catB(j)==5) then
            do e=1,3
              PrE(e) = PrE(e) * SUM(AKA2P(y,e,:) * PrGG)  ! TODO: FS of SB
            enddo
          endif 
                    
          if (ANY(UseEE/=0) .or. ANY(catB<0)) then
            if (catB(j)==2) then
              PrAB(x,y,:,1) = PrAB(x,y,:,1) * PrE
            else if (.not. all(PrE==1D0)) then
              PrAB(x,y,:,1) = PrAB(x,y,:,1) * SUM(PrE)
            endif
          else if (catB(j)==2) then
            do z=1,3
              PrAB(:,y,z,1) = PrAB(:,y,z,1) * PrE(z)
            enddo
          else if (catB(j)==12) then
            PrAB(:,y,:,1) = PrAB(:,y,:,1) * PrE(y)
!          else if (catB(j)==5 .and. SUM(PrE)<3) then
!            PrAB(x,y,:,1) = PrAB(x,y,:,1) * SUM(PrE)
          else if (.not. all(PrE==1D0)) then
            PrAB(:,y,:,1) = PrAB(:,y,:,1) * SUM(PrE)
          endif           

          do u=1, nFS(BB(j))
            PrE = PrE * OKA2P(Genos(l,FSID(u,BB(j))), y, :)
          enddo                 
          
          if (ANY(UseEE/=0) .or. ANY(catB<0)) then
            if (catB(j)==2) then
              PrAB(x,y,:,2) = PrAB(x,y,:,2) * PrE
            else if (.not. all(PrE==1D0)) then
              PrAB(x,y,:,2) = PrAB(x,y,:,2) * SUM(PrE)
            endif
            PrEE(:,x,nA+j) = PrE
          else if (catB(j)==2) then 
            do z=1,3
              PrAB(:,y,z,2) = PrAB(:,y,z,2) * PrE(z)
            enddo    
          else if (catB(j)==12) then
            PrAB(:,y,:,2) = PrAB(:,y,:,2) * PrE(y)    
          else if (.not. all(PrE==1D0)) then
            PrAB(:,y,:,2) = PrAB(:,y,:,2) * SUM(PrE) 
          endif  
        enddo   ! j

      else if (B>0) then 
        do z=1,3  
          do m=1,2
            if (kB/=3 .and. m/=kB)  cycle
            do j=1, nBx(m)
              Bj = BBx(j, m)
              if (nFS(Bj)==0 .and. parent(Bj,m)<0) cycle 
              if (ANY(FSID(:,Bj)==B) .and. DoneB(1)==1)  cycle   
              if (kA/=kB .and. PA<0 .and. Parent(Bj, kA)==PA) cycle
              if (kB==3 .and. Parent(Bj, 3-m)==GG(3-m) .and. GG(3-m)/=0) then  ! FS of B
                if (Parent(Bj,1)<0 .and. Parent(Bj,2)<0 .and. m==2) cycle
!                DoneB(1) = 2  ! for output check only
                do u=1, Max(nFS(Bj), 1)
                  if (FSID(u,Bj)==B) cycle  
                  if (ANY(AA(1:nA)==FSID(u,Bj))) cycle
                  if (ALL(UseEE==0)) then
                    PrAB(:,y,z,:) = PrAB(:,y,z,:) * OKA2P(Genos(l,FSID(u,Bj)),y,z) 
                  else
                    PrAB(x,y,z,:) = PrAB(x,y,z,:) * OKA2P(Genos(l,FSID(u,Bj)),y,z)
                  endif
                enddo
              else if (kB==3 .and. Parent(Bj,1)<0 .and. Parent(Bj,2)<0 &
               .and. MateLoop(j,m)) then  
                if (m==2) cycle
                call ParProb(l,Parent(Bj,m),m,-1,0,PrE)
                call ParProb(l,Parent(Bj,3-m),3-m,-1,0,PrW)
                
                do g=1, nS(-Parent(Bj, 3-m),3-m)
                  Ei = SibID(g,-Parent(Bj,3-m),3-m)
                  if (nFS(Ei) == 0) cycle
                  do e=1,3
                    do w=1,3
                      PrEW(e,w) = PrE(e) * PrW(w)
                      if (Parent(Ei,m)==Parent(Bj,m) .and. &
                       Parent(Ei,3-m)==Parent(Bj,3-m)) then
                        do i=1, nFS(Ei)
                          if (FSID(i, Ei)==B) cycle
                          PrEW(e,w) = PrEW(e,w) * OKA2P(Genos(l,FSID(i,Ei)), e, w)
                        enddo
                      else if (Parent(Ei,m)==Parent(Bj,m)) then
                        call ParProb(l, Parent(Ei,3-m),3-m, Ei, -1, PrH)
                        do i=1, nFS(Ei)
                          if (FSID(i, Ei)==B) cycle  
                          PrH = PrH * OKA2P(Genos(l,FSID(i,Ei)), :, e)
                        enddo
                        if (.not. all(PrH==1D0))  PrEW(e,w) = PrEW(e,w) * SUM(PrH)
                      else if (Parent(Ei,3-m)==Parent(Bj,3-m)) then
                        call ParProb(l, Parent(Ei,m),m, Ei, -1, PrH)
                        do i=1, nFS(Ei)
                          if (FSID(i, Ei)==B) cycle  
                          PrH = PrH * OKA2P(Genos(l,FSID(i,Ei)), :, w)
                        enddo
                        if (.not. all(PrH==1D0))  PrEW(e,w) = PrEW(e,w) * SUM(PrH)
                      endif
                    enddo  ! w
                  enddo  ! e
                enddo  ! sib g
                if (ALL(UseEE==0)) then
                  PrAB(:,y,z,:) = PrAB(:,y,z,:) * SUM(PrEW)
                else
                  do e=1,3
                    PrEE(:,x,nA+1) = SUM(PrEW(e,:))
                  enddo
                  PrAB(x,y,z,:) = PrAB(x,y,z,:) * SUM(PrEW)
                endif

              else
                if (kB==3) then
                  call ParProb(l,Parent(Bj,3-m),3-m,-1,0,PrE)
                else if (catB(j)==2 .or. catB(j)==12) then
                  PrE = 1D0
                else if (catB(j)==7) then
                  call ParProb(l, GpID(3-kB,-Parent(Bj,3-kB),3-kB), 3-kB, 0,0,PrH)
                  do e=1,3
                    PrE(e) = SUM(AKA2P(e,y,:) * PrH)
                  enddo              
                else if (catB(j)==8) then
                  call ParProb(l, GpID(3-kA,-Parent(Bj,3-kB),3-kB), 3-kA, 0,0,PrH)
                  do e=1,3
                    PrE(e) = SUM(AKA2P(e,x,:) * PrH)
                  enddo
                else if (UseEE(nA+1)/=0 .and. Bj==B) then  
                  call ParProb(l, MateABpar(nA+1), 3-TypeEE(nA+1), 0,0,PrH)
                  do e=1,3
                    do u=1, 3
                      PrW(u) = SUM(AKA2P(e,u,:) * PrEE(u,x,UseEE(nA+1)) * PrH)
                    enddo
                    PrE(e) = SUM(PrW)
                  enddo
                  PrE = PrE/SUM(PrE)
                else
                  call ParProb(l,Parent(Bj,3-m),3-m,-1,0,PrE)
                endif
                
                if (Parent(Bj,3-m)<0 .and. Parent(Bj,3-m)/=GG(3-m)) then
                  do g=1, nS(-Parent(Bj, 3-m),3-m)
                    Ei = SibID(g,-Parent(Bj,3-m),3-m)
                    if (nFS(Ei) == 0) cycle
                    if (ANY(AA(1:nA)==Ei)) cycle
                    if (Parent(Ei,m)==GG(m) .and. GG(m)/=0) cycle
                    do e=1,3
                      call ParProb(l, Parent(Ei, m), m, Ei, -1, PrH)
                      do i=1, nFS(Ei)
                        if (ANY(AA(1:nA)==FSID(i, Ei))) cycle
                        if (FSID(i, Ei)==B) cycle   
                        if (FSID(i, Ei)==PA) cycle
                        PrH = PrH * OKA2P(Genos(l,FSID(i,Ei)), :, e)
                      enddo
                      if (.not. all(PrH==1D0))  PrE(e) = PrE(e) * SUM(PrH)
                    enddo
                  enddo 
                endif
                do u=1, MAX(nFS(Bj),1)
                  if (FSID(u,Bj)==B) cycle 
                  if (FSID(u, Bj)==PA) cycle  
                  if (ANY(AA(1:nA)==FSID(u,Bj))) cycle
                  if (kB==3 .and. m==kA) then
                    PrE = PrE * OKA2P(Genos(l,FSID(u,Bj)), z, :) 
                  else
                    PrE = PrE * OKA2P(Genos(l,FSID(u,Bj)), y, :) 
                  endif
                enddo
                if (catB(j)==5 .and. kB/=3) then
                  do e=1,3
                    PrE(e) = PrE(e) * SUM(AKA2P(y,e,:) * PrGG)
                  enddo
                endif
                
                if (kB<3) then
                  if (ANY(UseEE/=0) .or. ANY(catB<0)) then
                    if (catB(j)==2) then 
                      PrAB(x,y,z,1) = PrAB(x,y,z,1) * PrE(z)
                    else if (SUM(PrE)<3.0) then
                      PrAB(x,y,z,1) = PrAB(x,y,z,1) * SUM(PrE)
                    endif
                  else if (catB(j)==2) then 
                    PrAB(:,y,z,1) = PrAB(:,y,z,1) * PrE(z)
                  else if (catB(j)==12) then 
                    PrAB(:,y,z,1) = PrAB(:,y,z,1) * PrE(y)
                  else if (SUM(PrE)<3.0) then
                    PrAB(:,y,z,1) = PrAB(:,y,z,1) * SUM(PrE)
                  endif

                  do u=1, MAX(nFS(Bj),1)
                    if (FSID(u,Bj)==B) then
                      PrE = PrE * OKA2P(Genos(l,B), y, :) 
                    endif
                  enddo
                  
                  if (ANY(UseEE/=0) .or. ANY(catB<0)) then
                    if (catB(j)==2) then 
                      PrAB(x,y,z,2) = PrAB(x,y,z,2) * PrE(z)
                    else if (SUM(PrE)<3.0) then
                      PrAB(x,y,z,2) = PrAB(x,y,z,2) * SUM(PrE)
                    endif
                    if (any(FSID(1:nFS(Bj),Bj)==B))  PrEE(:,x,nA+1) = PrE
                  else if (catB(j)==2) then 
                    PrAB(:,y,z,2) = PrAB(:,y,z,2) * PrE(z)
                  else if (catB(j)==12) then 
                    PrAB(:,y,z,2) = PrAB(:,y,z,2) * PrE(y)
                  else if (SUM(PrE)<3.0) then
                    PrAB(:,y,z,2) = PrAB(:,y,z,2) * SUM(PrE)
                  endif
                else if (SUM(PrE)<3.0) then  ! kB==3
                  PrAB(:,y,z,:) = PrAB(:,y,z,:) * SUM(PrE)
                endif
              endif
            enddo  ! j_m
          enddo  ! m
          if (kB==3) then
            PrAB(:,y,z,2) = PrAB(:,y,z,2) * OKA2P(Genos(l,B), y, z)
          endif
        enddo  ! z
      endif  ! B <0 vs B>0
      enddo  ! x (ANY(UseEE/=0)  only)
    enddo  ! y

    PrL(l) = LOG10(SUM(PrAB(:,:,:,2))) - LOG10(SUM(PrAB(:,:,:,1)))
!    PrLx(l,1) = LOG10(SUM(PrAB(:,:,:,1)))  ! for debug only
!    PrLx(l,2) = LOG10(SUM(PrAB(:,:,:,2)))                                                                                                
  endif
enddo
LL = SUM(PrL)


! if (AA(1)==19 .and. kA==1 .and. B==1575 .and. kB==3) then
  ! print *, ""
  ! write(*,'("UA: ", 2i5, ", ",2i5, ", ", 2i5, f9.2, 3i3, 2l2)') kA, A, kB,B, &
    ! GG, LL, catA(nA+1), catG, SIMPL, ALL(UseEE==0)
 ! ! write(*,'("PrLx: ", 2f9.2)')  SUM(PrLx(:,1)), SUM(PrLx(:,2))
  ! print *, "GGP:", GGP
  ! do i=1, nA
    ! write(*,'(i3, 3i7, ", ", 3i6, ", ", 2i4, i7)') i, AA(i), Parent(AA(i), :), &
     ! UseEE(i), TypeEE(i), MateABpar(i), catA(i), NFS(AA(i)), FSID(1, AA(i))
  ! enddo
  ! print *, ""
  ! if (kB<3) then
    ! do i=1, nB
       ! write(*,'(i3, 3i7, ", ", 3i4, ", ", 3i4)') nA+i, BB(i), Parent(BB(i), :),&
        ! UseEE(i+nA), TypeEE(i+nA), MateABpar(i+nA), catB(i), NFS(BB(i)), DoneB(i)
    ! enddo
  ! else
    ! print *, "nBx: ", nBx
    ! do i=1, nBX(1)
      ! print *, BBx(i,1), Mates(i,1), nFS(BBx(i,1)), MateLoop(i,1)
    ! enddo
    ! print *, "."
    ! do i=1, nBX(2)
      ! print *, BBx(i,2), Mates(i,2), nFS(BBx(i,2)), MateLoop(i,2)
    ! enddo
  ! endif
  ! print *, ""
  ! if (A<0) print *, "GP SA: ", GpID(:,-A,kA)
  ! if (A>0) print *, "Par A: ", Parent(A,:)
  ! if (B<0) print *, "GP SB: ", GpID(:,-B,kB), ", GGP: ", GGP
  ! if (B>0) print *, "Par B: ", Parent(B,:)
  ! print *, ""
! endif

 
deallocate(UseEE)
deallocate(PrEE)
deallocate(MateABpar)
deallocate(TypeEE)

end subroutine PairUA 

! #####################################################################

subroutine addFA(A, SB, kB, LL)  ! SB & partner-of-SB are GP's of A
use Global
implicit none

integer, intent(IN) :: A,SB,kB
double precision, intent(OUT) :: LL
integer :: x, y, z, Par, i, l, Bi
double precision :: PrL(nSnp, 2), PrXYZ(3,3,3,2), PrY(3), PrZ(3), PrPA(3), LLU, &
   ALRq, PrPAx(3)
logical :: AncOK, Ainbr            

LL = missing
if (Parent(A,kB) == -SB)  LL = Impossible
if (hermaphrodites/=0)  LL = NotImplemented
if (LL /= missing)  return 

AncOK = .TRUE.
call ChkAncest(-SB,kB, A,0, AncOK)
if (.not. AncOK .or. Parent(A,kB)<0) then 
  LL = NotImplemented
  return
endif

Par = 0
call getFSpar(SB, kB, .TRUE., Par)  ! TODO: strict=FALSE ?  needs check if Parent(B1,3-k)==Par
if (Par == 0 .or. (any(Parent(SibID(1:ns(SB,kB),SB,kB), 3-kB) /= Par .and. &
  Parent(SibID(1:ns(SB,kB),SB,kB), 3-kB) /= 0))) then  ! NOTE: all sibs made FS as side-effect
  LL = impossible
else
  call CalcAgeLR(A,Sex(A), Par,3-kB, kB,4, .TRUE., ALRq) 
  if (ALRq < 2.0*TF .or. ALRq==impossible)  LL = impossible
  if (LL /= impossible) then
    call ChkAncest(Par, 3-kB, A, 0, AncOK)
    if (.not. AncOK)  LL = impossible
  endif
endif
if (LL == impossible)  return

do i=1, ns(SB, kB)
  if (nFS(SibID(i,SB,kB))==0)  cycle
  if (Parent(SibID(i,SB,kB), 3-kB) == Par) then
    Bi = SibID(i,SB,kB)
    exit
  endif
enddo

if (Parent(A,3-kB)/=0 .and. Parent(A,3-kB)== Par) then
  Ainbr = .TRUE.  ! A would become inbred
else
  Ainbr = .FALSE.
endif

PrL = 0D0
do l=1, nSnp
  call ParProb(l, -SB, kB, -1, 0, PrY)
  call ParProb(l, Par, 3-kB, Bi, -1, PrZ)   
  call ParProb(l, Parent(A,3-kB), 3-kB, A, 0, PrPA)
  call ParProb(l, Parent(A,kB), kB, A, -4, PrPAx)
  do x=1,3
    do y=1,3
      do z=1,3
        PrXYZ(x,y,z,1) = PrY(y) * PrZ(z) * PrPAx(x) * AKA2P(x,y,z)  
        PrXYZ(x,y,z,2) = PrY(y) * PrZ(z) * PrPAx(x) * AHWE(x,l) 
        if (Ainbr) then
          PrXYZ(x,y,z,:) = PrXYZ(x,y,z,:) * OKA2P(Genos(l,A), x, z)
        else
          PrXYZ(x,y,z,:) = PrXYZ(x,y,z,:) * SUM(OKA2P(Genos(l,A), x, :) * PrPA)
        endif
        do i=1, ns(SB,kB)
          PrXYZ(x,y,z,:) = PrXYZ(x,y,z,:) * OKA2P(Genos(l, SibID(i,SB,kB)), y,z)
        enddo  ! NOT FSLik: some siblings may have parent(3-kB)=0
      enddo
    enddo
  enddo
  ! needs 2nd dim otherwise problem: Parent(A,k) gets included in LL. 
  PrL(l,1) = LOG10(SUM(PrXYZ(:,:,:,1)))
  PrL(l,2) = LOG10(SUM(PrXYZ(:,:,:,2)))
enddo

LLU = missing                                       
if (Parent(A,kB) /= 0 .or. Ainbr) then
  call CalcU(-SB, kB, A, Sex(A), LLU)
  LL = SUM(PrL(:,1)) - SUM(PrL(:,2)) + LLU 
else
  LL = SUM(PrL(:,1))
endif

end subroutine addFA

! #####################################################################

subroutine FAx(A, B, LL)  
! A result of FS mating, B offspr of A's double GPs  (complex = 'mono')
use Global
implicit none

integer, intent(IN) :: A, B
double precision, intent(OUT) :: LL
integer :: l, x, y, v, w
double precision :: PrL(nSnp), PrXY(3,3,3,3)

LL = Missing
if (any(Parent(A,:)/=0 .or. any(parent(B,:)/=0))) then
  LL = NotImplemented
  return
endif

PrL = 0D0
do l =1, nSnp
  do x=1,3  ! dam of A
    do y=1,3  ! sire of A
      do v=1,3  ! grandmother
        do w=1,3  ! grandfather
          PrXY(x,y,v,w) = OKA2P(Genos(l, A), x, y) * AKA2P(x,v,w) * AKA2P(y,v,w) * &
            OKA2P(Genos(l,B),v,w) * AHWE(w,l) * AHWE(v,l)
        enddo
      enddo
    enddo
  enddo
  PrL(l) = LOG10(SUM(PrXY))
enddo

LL = SUM(PrL)

end subroutine FAx

! #####################################################################

subroutine HSmating(A, kA, B, kB, hf, LL)   ! SB/PB and HS-of-PB are parents of SA/PA
use Global
implicit none

integer, intent(IN) :: A, kA, B, kB, hf
double precision, intent(OUT) :: LL
integer :: PA, PB, GA(2), GB(2), GGA(2), l, x, y, z, v, DoQuick, GGP  ! , m
double precision :: PrL(nSnp,2), PrXZ(3,3,3,3,2), PrA(3), PrPA(3), PrGA(3), &
  PrB(3), PrPB(3), PrGB(3), PrGGA(3), PrGGP(3), LLU

if (A > 0) then
  PA = Parent(A, kA)
else
  PA = A
endif
if (B > 0) then
  PB = Parent(B, kB)
else
  PB = B
endif

GA = getPar(PA, kA)
GGA = getPar(GA(3-kB), 3-kB)
GB = getPar(PB, kB)

LL = missing
if (PA > 0 .or. PB > 0) then
  LL = NotImplemented
else if (GA(kB) == PB .and. PB/=0) then
  LL = AlreadyAss
else if (GA(kB)/=0) then
  LL = impossible
else if (GB(kA) == PA .and. PA/=0) then
  LL = impossible
endif
if (LL /= Missing)  return

if (GGA(hf)/=0) then
  if (GB(hf)==0 .or. GB(hf)==GGA(hf)) then
    GGP = GGA(hf)
  else
    LL = impossible
  endif
else
  GGP = GB(hf)
endif

if (LL /= Missing)  return

if (PA==0 .or. (PB==0 .and. GGP==0) ) then  
  LL = NotImplemented   ! causes too many false neg 
  return                ! (indistinguishable from just PA inbred)
endif

DoQuick = 1
if (PA < 0) then
  call ChkDoQuick(-PA, kA, DoQuick)
  if (DoQuick /=1)  LL = NotImplemented
endif
if (PB < 0) then
  call ChkDoQuick(-PB, kB, DoQuick)
  if (DoQuick /=1)  LL = NotImplemented
endif  
if (LL /= Missing)  return

PrL = 0D0
do l=1, nSnp
  call ParProb(l, GA(3-kB), 3-kB, -4, 0, PrGA)  ! offspring contribution only
  if (GA(3-kB) > 0) then
    if (PB > 0) then
      call ParProb(l, GGP, hf, GA(3-kB), PB, PrGGP) 
    else
      call ParProb(l, GGP, hf, GA(3-kB), 0, PrGGP) 
    endif
  else
    if (PB > 0) then
      call ParProb(l, GGP, hf, 0, PB, PrGGP) 
    else
      call ParProb(l, GGP, hf, 0, 0, PrGGP) 
    endif
  endif
  call ParProb(l, GB(3-hf), 3-hf, 0, 0, PrGB)
  call ParProb(l, GGA(3-hf), 3-hf, 0, 0, PrGGA)
  if (A < 0) then
    PrA = XPr(1,:,l,-A,kA)
  else 
    call ParProb(l, PA, kA, A, -4, PrA) ! exclude both GPs & A
    call ParProb(l, Parent(A,3-kA), 3-kA, A, 0, PrPA)
    do x=1,3
      PrA(x) = PrA(x) * SUM(OKA2P(Genos(l,A), x, :) * PrPA)
    enddo
  endif
  if (B < 0) then
    PrB = XPr(1,:,l,-B,kB)
  else 
    call ParProb(l, PB, kB, B, -4, PrB)
    call ParProb(l, Parent(B,3-kB), 3-kB, B, 0, PrPB)
    do y=1,3
      PrB(y) = PrB(y) * SUM(OKA2P(Genos(l,B), y, :) * PrPB)
    enddo
  endif
  
  PrXZ = 0D0
  do x=1,3  ! SA/PA
    do y=1,3  ! SB
      do z=1,3  ! HS & mate of SB
        do v=1,3   ! parent of SB & mate-of-SB
          PrXZ(x,y,z,v,1) = PrA(x) * PrB(y) * AKA2P(x,y,z) * PrGA(z) * PrGGP(v) * &
            SUM(AKA2P(y,v,:) * PrGB) * SUM(AKA2P(z,v,:) * PrGGA)   
         ! SA similarly inbred, but B unrelated
         if (PB==0) then  
           PrXZ(x,y,z,v,2) = PrA(x) * AKA2P(x,y,z) * PrGA(z) * PrGGP(v) * &
              AKAP(y,v,l) * SUM(AKA2P(z,v,:) * PrGGA)
          endif
         
        enddo
      enddo
    enddo
  enddo
  PrL(l,1) = LOG10(SUM(PrXZ(:,:,:,:,1)))
  PrL(l,2) = LOG10(SUM(PrXZ(:,:,:,:,2)))
enddo

if (PB==0) then  ! B>0
  call CalcU(A, kA, B, kB, LLU)
  LL = SUM(PrL(:,1)) - (SUM(PrL(:,2)) + Lind(B)) + LLU
else
  LL = SUM(PrL(:,1))
endif

! if (A==-75 .and. kA==1 .and. B==6532)  then
  ! write(*, "('HSmating: ',5i6, 3f8.1)")  A, kA, B, kB, hf, SUM(PrL(:,1)), SUM(PrL(:,2)) + Lind(B), LLU
  ! write(*, '(5i6, "; ", 3i6)')  PA, GA, GGA, PB, GB
  ! print *, ''  
! endif


end subroutine HSmating

! #####################################################################

subroutine addHAselfed(SA,kA, SB,kB, LL)  ! SA result of selfing by SB
use Global
implicit none

integer, intent(IN) :: SA, kA, SB, kB
double precision, intent(OUT) :: LL
integer :: l, x, y
double precision :: PrL(nSnp), PrXY(3,3)

if (any(GpID(:,SA,kA) /= 0)) then
  LL = impossible
  return
endif

if (kA/=kB) then
  if (any(parent(SibID(1:ns(SA,kA),SA,kA), 3-kA) == -SB)) then
    LL = NotImplemented
    return
  endif
endif
! TODO: check if otherwise connected

PrL = 0D0
do l=1,nSnp
  do x=1,3
    do y=1,3
      PrXY(x,y) = XPr(1,x,l,SA,kA) * AKA2P(x,y,y) * XPr(3,y,l,SB,kB)
    enddo
  enddo
  PrL(l) = LOG10(SUM(PrXY))
enddo
LL = SUM(PrL)

end subroutine addHAselfed

! #####################################################################

subroutine pairCC(A,B,k, LL)  ! full 1st cousins
use Global
implicit none

integer, intent(IN) :: A,B,k
double precision, intent(OUT) :: LL
integer :: l, x, y, u, v, z
double precision :: PrL(nSnp, 5), PrXY(3,3), PrUV, PrPA(3), PrPB(3), &
  PrC(3,3,5), PrZ(3), PrXYf(3,3,2,2), LLself(2), LLU, LLtmp(3)
logical :: AncOK(2), AreHS, MaybeInbr(2)
  
LL = missing
if (Parent(A,k)/=0 .and. Parent(B,k)/=0) then
  if (Parent(A,k)==Parent(B,k)) then
    LL = impossible
  endif
endif

if (Parent(A,k)/=0 .or. Parent(B,k)/=0) then  ! TODO? See ParentHFS
  LL = NotImplemented
endif
if (LL/=missing) return

AncOK = .TRUE.                     
call ChkAncest(A,0,B,0, AncOK(1))
call ChkAncest(B,0,A,0, AncOK(2))
if (any(.not. AncOK)) then
  LL = impossible
  return
endif

AreHS = .FALSE.
MaybeInbr = .TRUE.                 
if (Parent(A,3-k)/=0 .and. Parent(A,3-k)==Parent(B,3-k)) then
    AreHS = .TRUE.
endif
if (hermaphrodites>0) then
  LLself = missing                        
  if (Parent(A,3-k)>0) then
    call PairSelf(B, Parent(A,3-k), LLself(1))
    call PairFullSib(B, Parent(A,3-k), LLself(2))
    if (LLself(1) > LLself(2)) then
      MaybeInbr(1) = .FALSE.
    endif
  endif
  if (Parent(B,3-k)>0) then
    call PairSelf(A, Parent(B,3-k), LLself(1))
    call PairFullSib(A, Parent(B,3-k), LLself(2))
    if (LLself(1) > LLself(2)) then
      MaybeInbr(2) = .FALSE.
    endif
  endif
endif

PrL = 0D0
PrXYf = 0D0  ! 4D: x, y, A/B inbred, CC/not
do l=1, nSnp
  if (.not. AreHS) then
    call ParProb(l, Parent(A,3-k), 3-k, A, 0, PrPA)
    call ParProb(l, Parent(B,3-k), 3-k, B, 0, PrPB)
  else
    call ParProb(l, Parent(A,3-k), 3-k, A, B, PrPA)
  endif
  
  PrC = 0D0
  do u=1,3  ! GG 3-k
    do v=1,3  ! GG k
      PrUV = AHWE(u,l) * AHWE(v,l)
      do x=1,3  !PA
        do y=1,3    !PB
          PrXY(x,y) = AKA2P(x,u,v) * AKA2P(y,u,v) * PrUV
          if (.not. AreHS) then
            if (any(MaybeInbr)) then
              PrXYf(x,y,1,1) = PrXY(x,y) * PrPA(u) / AHWE(u,l) 
              PrXYf(x,y,2,1) = PrXY(x,y) * PrPB(u) / AHWE(u,l)
              ! A or B inbred, not CC (reference)
              PrXYf(x,y,1,2) =  AKA2P(x,u,v) * PrUV * AHWE(y,l) * PrPA(u) / AHWE(u,l)
              PrXYf(x,y,2,2) =  AKA2P(y,u,v) * PrUV * AHWE(x,l) * PrPB(u) / AHWE(u,l)
            endif
            PrXY(x,y) = PrXY(x,y) * SUM(OKA2P(Genos(l,A), x, :) * PrPA) * &
              SUM(OKA2P(Genos(l,B), y, :) * PrPB)
            if (any(MaybeInbr)) then
              PrXYf(x,y,1,:) = PrXYf(x,y,1,:) * OKA2P(Genos(l,A), x, u) * &
                SUM(OKA2P(Genos(l,B), y, :) * PrPB)  ! A inbred
              PrXYf(x,y,2,:) = PrXYf(x,y,2,:) * SUM(OKA2P(Genos(l,A), x, :) * PrPA) * &
                OKA2P(Genos(l,B), y, u) 
              ! not considered: both A & B inbred. 
            endif
          else if (AreHS) then
            PrZ = PrPA
            do z=1,3
              PrZ(z) =PrZ(z) * OKA2P(Genos(l,A),x,z) * OKA2P(Genos(l,B),y,z)
            enddo
            PrXY(x,y) = PrXY(x,y) * SUM(PrZ)
          endif
        enddo
      enddo
      PrC(u,v,1) = SUM(PrXY)
      if (.not. AreHS) then
        if (MaybeInbr(1))  PrC(u,v,2) = SUM(PrXYf(:,:,1,1))
        if (MaybeInbr(2))  PrC(u,v,3) = SUM(PrXYf(:,:,2,1))
        if (MaybeInbr(1))  PrC(u,v,4) = SUM(PrXYf(:,:,1,2))
        if (MaybeInbr(2))  PrC(u,v,5) = SUM(PrXYf(:,:,2,2))
      endif
    enddo
  enddo 
  do z=1,5
    PrL(l,z) = LOG10(SUM(PrC(:,:,z)))
  enddo
enddo

LLtmp = 999D0
LLtmp(1) = SUM(PrL(:,1))
if (any(MaybeInbr)) then
  call CalcU(A, k, B, k, LLU)
  if (MaybeInbr(1))  LLtmp(2) = SUM(PrL(:,2)) - SUM(PrL(:,4)) + LLU
  if (MaybeInbr(2))  LLtmp(3) = SUM(PrL(:,3)) - SUM(PrL(:,5)) + LLU
endif
LL = MaxLL(LLtmp)

end subroutine pairCC

! #####################################################################

subroutine pairHSCC(A,B,LL)  ! full 1st cousins + HS
use Global
implicit none

integer, intent(IN) :: A, B
double precision, intent(OUT) :: LL
integer :: x, y, z, v,l
double precision :: PrL(nSnp), PrX(3,3,3,3)

if (any(parent(A,:)/=0) .or. any(Parent(B,:)/=0)) then
  LL = NotImplemented
  return
endif

PrL = 0D0
do l=1, nSnp
  do x=1,3
    do y=1,3
      do z=1,3
        do v=1,3
          PrX(x,y,z,v) = SUM(AKA2P(x,v,:) * AKA2P(y,v,:) * AHWE(v,l) * AHWE(:,l)) * &
            AHWE(z,l) * OKA2P(Genos(l,A),x,z) * OKA2P(Genos(l,B),y,z)
        enddo
      enddo
    enddo
  enddo
  PrL(l) = LOG10(SUM(PrX))
enddo
LL = SUM(PrL)

end subroutine pairHSCC

! #####################################################################

subroutine pairDHC(A, kA, B, withFS, LL)   ! B double half cousin of A
use Global
use CalcLik
implicit none

integer, intent(IN) :: A, kA, B
logical, intent(IN) :: withFS
double precision, intent(OUT) :: LL
integer :: x, y, z, v, w, ParB(2), m, l
double precision :: PrL(nSnp), PrX(3,3,3,3,3), PrPA(3), PrB(3,3)

if (B<0) then
  LL = NotImplemented
  Return
endif

parB = getPar(B,0)
LL = Missing
if (Parent(A,kA)/=0 .or. any(parB>0)) then
  LL = NotImplemented
else  
  do m=1,2
    if (ParB(m) < 0) then
      if (any(GpID(:,-ParB(m),m)/=0)) then
        LL = NotImplemented
      endif
    endif
  enddo
endif
if (LL == NotImplemented)  return
! also assumed that B has no additional siblings, asside from full sibs

PrL = 0D0
do l =1, nSnp
  call ParProb(l, Parent(A,3-kA), 3-kA, A, 0, PrPA)
  if (withFS) then
    PrB = FSLik(l,B)
  else
    PrB = OKA2P(Genos(l,B),:,:)
  endif
  do x=1,3  ! parent of A
    do y=1,3  ! grandparent 1
      do z=1,3  ! grandparent 2
        do v=1,3  ! dam of B
          do w=1,3  ! sire of B  
            PrX(x,y,z,v,w) = SUM(OKA2P(Genos(l,A), x, :) * PrPA) * AKA2P(x,y,z) * &
              AKAP(v,y,l) * AKAP(w,z,l) * AHWE(y,l) * AHWE(z,l) * PrB(v,w)
          enddo
        enddo
      enddo
    enddo
  enddo
  PrL(l) = LOG10(SUM(PrX))
enddo

LL = SUM(PrL)

end subroutine pairDHC

! #####################################################################

subroutine Clustering(PairID, PairType)
use Global
implicit none

integer, intent(IN) :: PairID(XP*nInd,2), PairType(XP*nInd)
integer :: k, x, n, m, ij(2), sx(2), topX, u, fcl, Par(2), topFS, topXi,t,viginti_pairs(10)
double precision :: LL(7,2), dLL, LLx(7, 2,2), dLLtmp(maxSibSize), dLLi
logical :: IsPair, FSM, DoLater

decem_pairs = mk_seq(nPairs, 10)      
t=1      
do x=1, nPairs
  if (MODULO(x,200)==0) call rchkusr()
  LL = missing
  ij = PairID(x,:)
  if (quiet==-1 .and. any(decem_pairs==x)) call print_progress(x,t)
    
  do k=1,2
    if (k/=PairType(x) .and. PairType(x)/=3)  cycle
    if (any(Parent(ij,k)>0)) cycle
    if (Parent(ij(1),k) == Parent(ij(2),k) .and. Parent(ij(1),k) /= 0) cycle
    if (hermaphrodites==1 .and. ANY(parent(ij,3-k)>0) .and. ANY(parent(ij,k)==0))  cycle
    fcl = 0
    LLx = missing                                  
    call getFocal(ij(1), ij(2), 0, k, fcl)     
    sx(1) = -Parent(ij(1),k)  
    sx(2) = -Parent(ij(2),k)
    
    if (k==1 .and. any(sx==0) .and. any(parent(ij,k)<0) .and. all(parent(ij,3-k)<0)) then
      DoLater = .FALSE.
      do n=1,2
        if (sx(n)==0)  cycle
        if (all(parent(SibID(1:ns(sx(n),k),sx(n),k), 3-k) < 0)) then
          call getFSpar(sx(n), k, .TRUE.,par(1))
          if (par(1) < 0) then
            if (nS(-par(1), 3-k) == nS(sx(n),k)) DoLater = .TRUE.  !else possibly add on wrong side i.o. merge
          endif 
        endif
      enddo
      if (DoLater)  cycle
    endif
    
    if (sx(1)==0 .and. sx(2)==0) then
      if (AgeDiff(ij(1),ij(2))==missing) then
        call CheckRel(ij(1), k, ij(2), k, fcl, LLx(:,1,1), LLx(:,1,2))
        call CheckRel(ij(2), k, ij(1), k, fcl, LLx(:,2,1), LLx(:,2,2))
        do u=1,7
          do n=1,2
            LL(u,n) = MaxLL(LLx(u,:,n))
          enddo
        enddo
      else if(AgeDiff(ij(1),ij(2))>=0) then
        call CheckRel(ij(1), k, ij(2), k, fcl, LL(:,1), LL(:,2))
      else
        call CheckRel(ij(2), k, ij(1), k, fcl, LL(:,1), LL(:,2))
      endif
      
      IsPair = .TRUE.
      do n=1,2  
        if (AgePhase==0 .and. n==2)  cycle
        if (AgePhase==2 .and. n==1)  cycle      
        if (LL(2,n)>0 .and. LL(3,n)>0) then
          IsPair = .FALSE.
          exit
        endif
        topX = 0
        dLL = 0D0                                 
        call BestRel(LL(:,n), fcl, topX, dLL)         
        if (fcl==3 .and. (topX==2 .or. topX==3)) then
          IsPair = .TRUE.
        else if (fcl==2 .and. topX==2 .and. dLL > TA .and. &
         (LL(2,n) - MaxLL(LL((/1,4,5,6,7/),n)) > 2*TA .or. Complx==0 .or. &
           ALL(Parent(ij, 3-k)/=0))) then    
          IsPair = .TRUE.
        else if (AgePhase==2 .and. n==1 .and. Complx>0) then   
          ! still do a basic check on non-age-dependend
          if (MaxLL(LL((/2,3/),n)) - MaxLL(LL((/1,6,7/),n))>TA .and. &
            MaxLL(LL((/4,5/),n)) - MaxLL(LL((/2,3/),n)) < TA) then
              IsPair = .TRUE.
          else
            IsPair = .FALSE.
            exit
          endif
        else
          IsPair = .FALSE.
          exit
        endif
      enddo
      if (.not. IsPair) cycle
     
      call NewSibship(ij(1), ij(2), k)   ! new sibship (pair)  
!      print *, "new pair: ", k, ij
      if (fcl==2 .or. (topX==2 .and. dLL>2*TA)) then
        if (ALL(Parent(ij, 3-k)==0)) then  ! another new sibship
          call NewSibship(ij(1), ij(2), 3-k)
        else
          do m=1,2
            if (Parent(ij(m),3-k)/=0) then
              call setPar(ij(3-m), 3, Parent(ij(m),3-k), 3-k)
            endif
          enddo
        endif
      endif
      
    else if (sx(1)>0 .and. sx(2)>0 .and. sx(1) /= sx(2)) then
      IsPair = .TRUE.
      FSM = .FALSE.     
      call CheckMerge(sx(1), sx(2), k,k, 1, LL(:,1), LL(:,2), FSM)
      do n=1,2   
        if (AgePhase==0 .and. n==2)  cycle
        if (AgePhase==2 .and. n==1)  cycle    
        topX = 0
        dLL = 0D0    
        call BestRel(LL(:,n), 1, topX, dLL)
        if (topX /=1 .or. dLL < TA * dble(MIN(nS(sx(1),k), nS(sx(2),k))) &
         .or. (fcl==2 .and. (.not. FSM .or. &
          dLL < 2.0*TA * dble(MIN(nS(sx(1),k), nS(sx(2),k)))))) then
          IsPair = .FALSE.
          exit
        endif
      enddo
      if (.not. IsPair) cycle
      
      if (FSM .and. fcl==2) then
        call DoFSmerge(sx(1), sx(2), k)
      else 
        call DoMerge(sx(1), sx(2), k)
      endif
      
    else
      do m=1,2
        if (sx(m)>0 .and. sx(3-m)==0) then
          IsPair = .TRUE.
          FSM = .FALSE.
          call CheckRel(ij(3-m), 0, -sx(m), k, fcl, LL(:,1), LL(:,2)) 
          do n=1,2
            if (AgePhase==0 .and. n==2)  cycle
            if (AgePhase==2 .and. n==1)  cycle
            topX = 0
            dLL = 0D0 
            call BestRel(LL(:,n), fcl, topX, dLL)
            if (.not. (topX==fcl .or. (fcl==3 .and. topX==2))) then
              IsPair = .FALSE.
              exit
            endif
          enddo
          if (.not. IsPair) cycle
          LLx = missing                                                     
          if (Complx>1 .and. ANY(SibID(1:ns(sx(m),k), sx(m), k) == Parent(ij(3-m),3-k))) then  ! inbreeding
            call CheckPair(ij(3-m), Parent(ij(3-m),3-k), k, 3, LLx(:,1,1), LLx(:,1,2))
            call BestRel(LLX(:,1,2), 3, topXi, dLLi)
            if (topXi /= 3) then
              IsPair = .FALSE. 
              cycle
            endif
          endif
          if (topX==2 .and. fcl/=2) then  
            call BestRel(LL(:,1), 2, topX, dLL)  
            if (topX==2 .and. dll > 2*TA) then 
              FSM = .TRUE.
            endif
          endif
          if (fcl==2 .or. FSM) then
            Par = 0
            topFS = 0
            dLLtmp = missing
            call getFSpar(sx(m), k, .TRUE., Par(m))   
            if (Par(m)/=0 .and. Parent(ij(3-m), 3-k) == Par(m)) then
              ! do nothing
            else if (Par(m)/=0 .and. all(Parent(SibID(1:ns(sx(m),k),sx(m),k),3-k)==Par(m))) then
              call setPar(ij(3-m), 3, Par(m), 3-k)
            else 
              call AddFS(ij(3-m), sx(m), k,0,k, LL(2,1), topFS, dLLtmp)
              if (Complx==0) then
                  call setPar(ij(3-m), 3, Parent(topFS, 3-k), 3-k)  
              else if (topFS>0) then
                if (parent(ij(3-m),3-k) == Parent(topFS, 3-k) .and. Parent(topFS, 3-k)/=0) then
                  ! do nothing
                else if (MAXVAL(dLLtmp, mask=dLLtmp<impossible)>2*TA) then
                  call CheckPair(ij(3-m), topFS, k, 2, LL(:,1), LL(:,2))
                  call BestRel(LL(:,2), 2, topX, dLL)
                  if (topX==2 .and. dll > 2*TA) then
                    if (Parent(topFS, 3-k)/=0) then
                      call setPar(ij(3-m), 3, Parent(topFS, 3-k), 3-k)
                    else if (Parent(ij(3-m), 3-k)/=0) then
                      call setPar(topFS, 3, Parent(ij(3-m), 3-k), 3-k)
                    else
                      call NewSibship(ij(3-m), topFS, 3-k)   ! new sibship (pair)      
                    endif
                  endif
                endif
              endif
            endif     
          endif     

          call setPar(ij(3-m), 3, -sx(m), k)
        endif
      enddo
    endif
  enddo 
enddo

end subroutine Clustering

! #####################################################################

subroutine Merging ! check if any of the existing clusters can be merged
use Global
implicit none

integer :: k, s, r, topX, xr, n
double precision :: LLm(7,2), dLL
logical :: FSM, OK

do k=2,1,-1
!  if (Complx==0 .and. k==2) cycle
  do s=1,nC(k)-1
    if (modulo(s,20)==0)  call rchkusr()         
    if (modulo(s,20)==0 .and. quiet==-1)  call Rprint("", (/k, s/), (/0D0/), "INT")
    if (s >= nC(k)) exit
    r = s
    do xr=s+1, nC(k)
      r = r + 1
      if (r > nC(k)) exit   ! possible due to merged sibships    
      LLm = missing
      FSM = .FALSE.      
      call CheckMerge(s, r, k, k, 1, LLm(:,1), LLm(:,2), FSM) 
      if (LLM(1,2) > 0 .or. LLM(1,2) < LLM(7,2)) cycle
      if (.not. FSM .and. (Complx==0 .or. Hermaphrodites==1))  cycle     
      OK = .TRUE.
      topX = 0
      dLL = missing      
      if (all(GpID(:,s,k) /= 0) .and. GpID(1,s,k)==GpID(1,r,k) .and. GpID(2,s,k)==GpID(2,r,k) .and. &
        ABS(LLM(7,1) - LLM(1,1)) < 0.1 .and. LLM(1,2) > LLM(7,2)) then
         OK = .TRUE.      ! identical LL s & r  FS vs. identical. occams razor --> merge
      else     
        do n=1,2    ! UseAge = (/.FALSE., .TRUE./)
          if (AgePhase==0 .and. n==2)  cycle
          if (AgePhase==2 .and. n==1)  cycle          
          call BestRel(LLm(:,n), 1, topX, dLL)
          if (topX /=1 .or. dLL < TA * dble(MIN(nS(s,k), nS(r,k)))) then
            OK = .FALSE.
            exit   
          endif
        enddo
      endif
      if (.not. OK)  cycle
!      write(*,'("Merge: ", 3i4, ": ", 3i6, " + ", 3i6)') k, s, r, SibID(1:3,s,k), SibID(1:3, r, k)
            
      if (FSM .and. (dLL > 2.0*TA * dble(MIN(nS(s,k), nS(r,k))) .or. &
        Complx==0 .or. Hermaphrodites==1)) then
        call DoFSmerge(s, r, k)
      else 
        call DoMerge(s, r, k)
      endif      
      r = r-1  ! otherwise a cluster is skipped
    enddo
  enddo
enddo

end subroutine Merging

! #####################################################################

subroutine SibParent  
! for each sibship, check if a real indiv can replace the dummy parent
use Global
use OHfun
implicit none

integer :: k, s, xs, i, n, topX, CurNumC(2), Par, SClone, &
  j, nCandPar, CandPar(mxCP), h, SibTmp(maxSibSize), nSib, sib1, sxSib(maxSibSize) 
double precision :: LL(7), dLL, LLtmp(7,2), ALR, LLO, LR, LLg(7)
logical :: NeedsOppMerge, Maybe, MaybeOpp, AncOK, FSM

CurNumC = nC
do k=1,2
  s = 0
  do xs=1, CurNumC(k)
    s = s+1
    if (modulo(s,20)==0)  call rchkusr()       
    if (modulo(s,20)==0 .and. quiet==-1)  call Rprint("", (/k,s/), (/0D0/), "INT")
    if (s > nC(k)) exit   
    if (hermaphrodites==1) then
      call getFSpar(s, k, .TRUE., Par)
      if (Par < 0)  cycle
    endif 
    nCandPar = 0
    CandPar = 0
    NeedsOppMerge = .FALSE.                       
        
    do i=1,nInd
      Maybe = .TRUE.
      MaybeOpp = .FALSE.
      if (nCandPar == mxCP) exit  !unlikely
      if (Sex(i)/=k .and. Sex(i)<3) cycle
      if (Parent(i,k)==-s) cycle
      if (ANY(GpID(:,s,k)==i)) cycle
      if (DumClone(s,k)/=0 .and. sex(i)/=4)  cycle
      do n=1,nS(s,k)
        if (AgeDiff(SibID(n,s,k), i) <= 0) then
          maybe = .FALSE.
          exit
        endif
      enddo
      if (.not. Maybe) cycle      
      if (DoMtDif) then
        if (k==1 .and. mtDif(SibID(1,s,k), i))  cycle 
      endif      
      call CalcAgeLR(-s, k, i, k, k, -1, .TRUE., ALR)
      if (ALR==impossible .or. ALR < 3*TF)  cycle
      call ChkAncest(i,k, -s,k, AncOK)
      if (.not. AncOK)  cycle
      do n=1,nS(s,k)
        call CalcP2(SibID(n,s,k), 3, i, Parent(SibID(n,s,k),3-k), k, LR)
        if (LR == impossible .or. LR < 3*TF) then
          maybe = .FALSE.
          exit
        endif 
      enddo
      if (.not. Maybe) cycle
      call CalcP2(i, Sex(i), GpID(1,s,k), GpID(2,s,k), 1, LR)
      if (LR == impossible .or. LR < 3*TF)  cycle
      call QPO(i, s, k, LR)
      if (LR < TF*nS(s,k))  cycle
      LLg = missing
      LL = missing
      topX = 0
      call CheckRel(-s, k, i, k, 1, LLg, LL)
      if (AgePhase <=1) then
        call BestRel(LLg, 1, topX, dLL)
      else
        call BestRel(LL, 1, topX, dLL)
      endif
      if (topX/=1) then
        Maybe = .FALSE.
        cycle
      endif
      if (Sex(i)>2 .and. DumClone(s,k)==0) then  ! check if parent of opposite sex instead
        MaybeOpp = .TRUE.
        Par = 0       
        call getFSpar(s, k, .TRUE., Par)
        if (Par > 0) then
          MaybeOpp = .FALSE.
        else if (Par/=0 .and. Parent(i, 3-k) == Par) then
          MaybeOpp = .FALSE.  ! are HS
        else if (Par==0) then
          if (ANY(Parent(SibID(1:nS(s,k), s,k),3-k)>0)) then
            MaybeOpp = .FALSE.
          else  ! check if could all be FS
            do j=1, ns(s,k)-1
              do h=j+1, ns(s,k)
                call CalcAgeLR(sibID(j,s,k), 0, SibID(h,s,k), 0, 0, 2, .TRUE., ALR)
                if (ALR == impossible .or. ALR < 3.0*TF) then
                  MaybeOpp = .FALSE.
                  exit
                endif                 
              enddo
              if (.not. MaybeOpp) exit
            enddo
          endif
          if (MaybeOpp) then
            call OppMerge(s,k,LLO)
            if (LLO>NotImplemented .or. (CLL(s,k) - LLO) > ns(s,k)*TF) then
              MaybeOpp = .FALSE.
            endif
          endif                
        else if (Par < 0) then
          do n=1,nS(-Par,3-k)
            if (CalcOH(i, SibID(n,-par,3-k)) > maxOppHom) then
              MaybeOpp = .FALSE.
              exit
            endif   
          enddo
          if (MaybeOpp) then
            do n=1,2
              if (GpID(n,-Par,3-k) <= 0) cycle
              if (CalcOH(i, GpID(n,-Par,3-k)) > maxOppHom) then
                MaybeOpp = .FALSE.
                exit
              endif   
            enddo
          endif                        
        endif
        if (MaybeOpp) then
          LLtmp = missing
          topX = 0
          if (Par < 0) then  ! may have more/fewer sibs
            call CheckRel(Par, 3-k, i, 3-k, 1, LLtmp(:,1), LLtmp(:,2))
            if (AgePhase <=1) then
              call BestRel(LLtmp(:,1), 1, topX, dLL)
            else
              call BestRel(LLtmp(:,2), 1, topX, dLL)
            endif
            if (topX/=1)  MaybeOpp = .FALSE.
          else if (Par == 0 .and. ns(s,k)>0) then
            sib1 = SibID(1,s,k)
            call PairPO(sib1, i, 3-k, 0, LLtmp(1,2))   
            call CalcU(sib1, k, i, 3-k, LLtmp(2,2))
            if (LLtmp(1,2)>0 .or. (LL(1)-LL(7)) - (LLtmp(1,2)-LLtmp(2,2)) > &
             TA*ns(s,k))  MaybeOpp = .FALSE.
          endif
        endif
        if (MaybeOpp) cycle
      endif

      if (Complx==0) then  ! ensure monogamous
        Par = 0       
        call getFSpar(s, k, .FALSE., Par)     
        if (Mate(i)/=Par .and. Mate(i)/=0 .and. Par/=0) then
          LLtmp = missing
          if (Mate(i)<0 .and. Par<0) then
            call CheckMerge(-Par, -Mate(i), 3-k, 3-k, 8, LLtmp(:,1), LLtmp(:,2), FSM)
          else if (Mate(i)>0 .and. Par<0) then
            call CheckRel(Par, 3-k, Mate(i), 3-k, 1, LLtmp(:,1), LLtmp(:,2))
          endif
          if (LLtmp(1,2) < 0 .and. (LLtmp(1,2) - MaxLL(LLtmp(2:7,2)) > -TA)) then
            NeedsOppMerge = .TRUE.
          else
            cycle
          endif
        endif
      endif

      if (Maybe) then               
        nCandPar = nCandPar + 1
        CandPar(nCandPar) = i
      endif
    enddo  ! i  
        
    if (nCandPar == 1) then
      Par = 0       
      if (NeedsOppMerge) then
        call getFSpar(s, k, .FALSE., Par)
        if (Mate(CandPar(1)) < 0) then
          call DoMerge(-Mate(CandPar(1)), -Par, 3-k)  ! Par gets removed
        else
          call getOff(Par, 3-k, .TRUE., nSib, SibTmp, sxSib)
          do n=1,nSib
            call setPar(SibTmp(n), sxSib(n), Mate(CandPar(1)), 3-k)
          enddo
          call DoMerge(0, -Par, 3-k)
        endif
      else if (Complx==0 .and. Mate(CandPar(1))==0) then
        Mate(CandPar(1)) = Parent(SibID(1,s,k), 3-k)
      endif
      
      if (hermaphrodites/=0 .and. Sex(CandPar(1))==4) then
        SClone = DumClone(s,k)
        if (SClone /= 0) then
          call getOff(-s, k, .FALSE., nSib, SibTmp, sxSib) 
          do n=1,nSib 
            call setParTmp(SibTmp(n), sxSib(n), CandPar(1), k)  ! else conflict w CheckSelfed()
          enddo 

          call getOff(-SClone, 3-k, .TRUE., nSib, SibTmp, sxSib)
          do n=1,nSib 
            call setPar(SibTmp(n), sxSib(n), CandPar(1), 3-k)
          enddo  
          call DoMerge(0, SClone, 3-k)  !removes Sclone cluster 
        endif
      endif

      call getOff(-s, k, .TRUE., nSib, SibTmp, sxSib)   ! includes dummy sibs
      do n=1,nSib 
        call setPar(SibTmp(n), sxSib(n), CandPar(1), k)
      enddo  
      call DoMerge(0, s, k)  !removes cluster s 
      s = s-1  ! otherwise a cluster is skipped
!    else   
!       TODO?
    endif
  enddo ! s
enddo ! k
  
end subroutine SibParent

! #####################################################################

subroutine MoreParent
! for each individual, check if a parent can be assigned now.
use Global
use OHfun
use CalcLik
implicit none

integer :: x, i, j, k, s, curPar(2), nCP(2), CandPar(mxCP, 2), &
   TopTmp, BYrank(nInd),t
double precision :: LLP(2), ALR(2), LL(7,2), LRQ, dLL, LRFS(mxCP,2,2), LLtmp(7,2)
logical :: DoNewPars, AncOK, DropS, KeepOld

if (hermaphrodites==1 .or. ALL(BY==BY(1))) then
  DoNewPars = .FALSE.
else
  DoNewPars = .TRUE.   ! do check for additional parent-offspring pairs. 
endif

call getRank_i(BYrank)
t=1
do x=1, nInd
  if (MODULO(x,100)==0)  call rchkusr()
  if (quiet==-1 .and. any(viginti==x)) call print_progress(x,t)
  i = BYRank(x)
  if (skip(i))  cycle
  if (ALL(Parent(i,:)/=0) .and. .not. ToCheck(i)) cycle
  call CalcLind(i)
  CurPar = Parent(i,:)
  nCP = 0
  CandPar = 0
  LRFS = missing  
  
  do k=1,2
    if (curPar(k)/=0) then 
      nCP(k) = nCP(k) +1
      CandPar(nCP(k), k) = curPar(k)
    endif
    call setParTmp(i, Sex(i), 0, k)   
    call SetEstBY(i, Sex(i))
    call SetEstBY(curPar(k), k)
  enddo
  call SetEstBY(i, Sex(i))
    
  do k=1,2
    if (Complx==0 .and. k==2 .and. all(parent <= 0)) exit  ! TODO double check all(par..
    if (nC(k)==0)  cycle
    do s=1, nC(k)   
      if (nCP(k) == mxCP)  exit
      if (ANY(CandPar(:,k) == -s))  cycle    
      if (ANY(GpID(:,s,k)==i)) cycle   
      call ChkAncest(-s, k, i, sex(i), AncOK)
      if (.not. AncOK)  cycle
      call CalcAgeLR(i,0,-s,k,0,1, .TRUE.,ALR(1))
      if (ALR(1)==impossible .or. ALR(1) < (3.0*TF))  cycle
      if (DoMtDif) then
        if (k==1 .and. mtDif(i, SibID(1,s,k)))  cycle    
      endif
      if (Complx==0) then
        call QFS(i,s,k,LRQ)       
        if (LRQ < TF) cycle
      else
        call Qadd(i, s, k, LRQ)
        if (LRQ < ns(s,k)*TF) cycle   
      endif
      LL = missing
      LLtmp = missing
      if (DumClone(s,k)==0) then  
        call SibChk(i,s,k,3, 1, LL(:,1))
        if (LL(3,1)>0 .or. MaxLL(LL(2:3,1)) - LL(7,1) < TA .or. & 
         (MaxLL(LL(2:3,1)) - MaxLL(LL(4:6,1)) < TF .and. ANY(LL(4:6,1) < 0))) cycle 
      endif
      
      if (Complx>1 .and. ANY(SibID(1:ns(s,k), s, k) == Parent(i,3-k))) then  ! inbreeding
        call CheckPair(i, Parent(i,3-k), k, 3, LLtmp(:,1), LLtmp(:,2))
        call BestRel(LLtmp(:,2), 3, topTmp, dLL)
        if (topTmp /= 3)  cycle
      endif

      nCP(k) = nCP(k) +1
      CandPar(nCP(k), k) = -s
      if (Complx==0 .and. nCP(3-k)<mxCP) then
        nCP(3-k) = nCP(3-k) +1
        if (DumMate(s,k) == 0) call Erstop("Mono error: no DumMate", .TRUE.)
        CandPar(nCP(3-k), 3-k) = DumMate(s,k)
      endif
    enddo
  enddo

  if (DoNewPars)  then
    do j=1, nInd  ! candidate parent.
      if (i==j) cycle
      if (ANY(CandPar==j) .and. Sex(j)<3) cycle  ! already included 
      if (ANY(Parent(j,:)==i) .or. ANY(Parent(i,:)==j)) cycle
      if (AgeDiff(i,j) <= 0)  cycle  ! note: unknown = missing > 0
      if (DoMtDif) then      
        if (Sex(j)==1 .and. mtDif(i,j))  cycle   
      endif
      call ChkAncest(j, sex(j), i, sex(i), AncOK)
      if (.not. AncOK)  cycle
      if (Sex(j) < 3) then
        if (Parent(i, Sex(j)) < 0) cycle   ! replacings done elsewhere
        if (nCP(Sex(j))==mxCP) cycle
      else 
        if (ANY(nCP == mxCP)) cycle
      endif
      ALR = missing
      LLP = missing 
      if (CalcOH(i,j) > maxOppHom)  cycle
      if (QLR_PO(i,j) < 5.0*TF)  cycle
!      if ((LLR_O(i,j)==missing .or. LLR_O(i,j) < 5.0*TF))  cycle 
      call CalcAgeLR(i,sex(i), j, sex(j), 0, 1, .TRUE., ALR(1))
      if (ALR(1) == impossible .or. ALR(1)<= 3.0*TF)  cycle
      call CalcAgeLR(j, sex(j), i,sex(i), 0, 1, .TRUE., ALR(2))
      if (ALR(2) /= impossible .and. (ALR(1)-ALR(2)) < TF)  cycle 
       if (Sex(j)<3) then
        call PairPO(i, j, sex(j), 0, LLP(1))   
      else
        call PairPO(i, j, 1, 0, LLP(1))
      endif 
      if (LLP(1) > 0) cycle
      call CalcU(i,sex(i), j, sex(j), LLP(2))
      if ((LLP(1) - LLP(2)) < TF) cycle  
      do k=1,2
        if (Sex(j)<3 .and. Sex(j)/= k) cycle
        if (DoMtDif) then
          if (k==1 .and. mtDif(i,j))  cycle  ! when sex(j)>=3 
        endif
        if (nCP(k) < mxCP .and. .not. any(candPar(:,k) == j))  then
          nCP(k) = nCP(k) +1
          CandPar(nCP(k), k) = j
        endif
        if (Complx==0 .and. Mate(j)/=0 .and. nCP(3-k)<mxCP) then
          nCP(3-k) = nCP(3-k) +1
          CandPar(nCP(3-k), 3-k) = Mate(j)
        endif
      enddo
    enddo  ! j
  endif
  
  KeepOld = .FALSE.
  if (ALL(nCP <=1) .and. ALL(candPar(1,:) == curPar)) then
    KeepOld = .TRUE.
    ! do k=1,2
      ! if (curPar(k) < 0) then
        ! if (IsNewSibship(-curPar(k), k))  KeepOld = .FALSE.
      ! endif
    ! enddo
    if (KeepOld) then  
      do k=1,2
        call setParTmp(i, Sex(i), curPar(k), k)  ! restore
        call SetEstBY(i, Sex(i))
        call SetEstBY(curPar(k), k)
      enddo
      call SetEstBY(i, Sex(i))
      ToCheck(i) = .FALSE.
      cycle
    endif
  endif  
  
  call SelectParent(i, Sex(i), nCP, CandPar, .FALSE., .FALSE.)

  if (ANY(Parent(i,:)/=CurPar)) then
    DropS = .FALSE.               
    if (Complx == 0)  call UpdateMate(i, Sex(i), curPar, .FALSE.)
    do k=1,2     
!      if (curPar(k)/=0 .and. parent(i,k)==0)   write(*,'(i5, " Dropped: ", 2i5, " -> ", 2i5)') i, curpar, parent(i,:)
      if (curPar(k)<0 .and. Parent(i,k)/=curPar(k)) then
        call CheckDropSibship(-curPar(k), k, DropS) 
        if (hermaphrodites/=0 .and. .not. DropS) then
          call CheckSelfed(curPar(k),k)
        endif
      endif
    enddo
    call setEstBY(i,sex(i)) 
  else
    ToCheck(i) = .FALSE.
  endif
enddo

end subroutine MoreParent

! #####################################################################

subroutine UpdateMate(A, kA, OldPar, ParOnly)
use Global
implicit none

integer, intent(IN) :: A, kA, OldPar(2)
logical, intent(IN) :: ParOnly  ! T:only parents / F:also dummy parents
integer :: NewPar(2), m, nOff, sxOff(maxSibSize), Off(maxSibSize), x

if (Complx /= 0)  return

NewPar = getPar(A, kA)
if (all(NewPar == OldPar))  return

do m=1,2
  if (oldPar(m) > 0 .and. oldpar(m) /= newpar(m)) then
    call getOff(oldPar(m), m, .TRUE., nOff, Off, sxOff)
    if (nOff == 0) then
      if (any(Mate == OldPar(m))) then
        x = MINLOC(ABS(Mate - oldPar(m)), DIM=1)
        Mate(x) = 0
      else if (any(DumMate(:,3-m) == OldPar(m))) then
        x = MINLOC(ABS(DumMate(:,3-m) - oldPar(m)), DIM=1)
        DumMate(x, 3-m) = 0
      endif
    endif
  ! else if oldpar(m) < 0, CheckDropSibship will take care of it
  endif  
enddo

do m=1,2
  if (newpar(m) > 0) then
    if (Mate(newpar(m)) == 0) then
      Mate(newpar(m)) = newpar(3-m)
    else if (Mate(newpar(m)) /= newpar(3-m)) then
      call Erstop("Something going wrong with Mate", .TRUE.)
    endif
  else if (newpar(m) < 0) then
    if (DumMate(-newpar(m),m) == 0) then
      DumMate(-newpar(m),m) = newpar(3-m)
    else if (DumMate(-newpar(m),m) /= newpar(3-m)) then
      call Erstop("Something going wrong with DumMate", .TRUE.)
    endif
  endif 
enddo

if (A > 0 .and. any(newPar == 0) .and. any(NewPar/=0) .and. .not. ParOnly) then   ! create singleton sibship
  do m=1,2
    if (NewPar(m) == 0) then
      call NewSibship(A, 0, m)   ! takes care of (Dum)Mate
    endif
  enddo
endif

end subroutine UpdateMate

! #####################################################################

subroutine getfocal(A, B, s, k, focal)
use Global
implicit none

integer, intent(IN) :: A, B, s, k
integer, intent(OUT) :: focal
integer :: j, BB(MaxsibSize), nB, opParB, opParA

BB = 0
nB = 0      
if (B/=0) then
  BB(1) = B
  nB = 1
else if (s/=0) then
  BB = SibID(:,s,k)
  nB = nS(s,k)
else
  call ErStop("getFocal: B=0 and s=0", .TRUE.)
endif

focal = 3
if (Complx==0) then
  focal = 2  ! FS
else if (nB==0) then  ! empty sibship ; should only happen w hermaphrodites or Complx=0 
  focal = 3  
else if (Parent(A,k)==0 .and. Parent(A,3-k)/=0 .and. &
  any(Parent(BB(1:nB), 3-k) == Parent(A,3-k))) then
  focal = 2
else if (ALL(Parent(BB(1:nB),k)==0) .and. Parent(A,3-k)/=0 .and. &
  ALL(Parent(BB(1:nB), 3-k) == Parent(A,3-k))) then
  focal = 2
else if (ALL(Parent(A,:)==0) .and. ALL(Parent(BB(1),:)==0) .and. hermaphrodites/=2) then
  focal = 2
else if (all(Parent(A,:)<=0) .and. Parent(BB(1),k)<0 .and. ALL(Parent(BB(1:nB),3-k)<=0)) then
  call getFSpar(-Parent(BB(1),k), k, .TRUE., opParB)
  if (opParB < 0) then
    if (ALL(Parent(A,:)==0)) then
      focal = 2  ! FS add
    else if (Parent(A,k)<0) then
      call getFSpar(-Parent(A,k), k, .TRUE., opParA)
      if (opParA < 0) then
        focal = 2     ! FS merge
      endif 
    endif
  endif
endif

if (focal==2) then  ! exception: cannot be /unlikely FS based on age
  do j=1,nB 
    if (getAP(AgeDiff(A,BB(j)), 2, 0,0, Impossible) == Impossible) then
      focal = 3
      exit
    else if (AgePhase == 2) then
      if (getAP(AgeDiff(A,BB(j)), 3, 0, k, log10(zero)) - &
       MAX(getAP(AgeDiff(A,BB(j)), 2, 0,0, log10(zero)), &
           getAP(AgeDiff(A,BB(j)), 3, 0, 3-k, log10(zero))) > 2.0*ABS(TF)) then
        focal = 3
      endif
    endif
  enddo 
endif

end subroutine getfocal

! #####################################################################

subroutine SibGrandparents 
! for each sibship, find the most likely male and/or female grandparent
use Global
implicit none

integer :: k, s, i, r, m, par, x, nCG(2), candGP(mxCP, 2), curGP(2), &
   ix, not4(5), BYrankC(nInd/2,2), n
double precision :: LRG, ALRtmp(2), LLX(3,2), dx(maxSibSize), LLA(7)
logical :: skipCluster(maxval(nC),2), AncOK, DropS, Maybe

skipCluster = .FALSE.
do k=1,2
  if (nC(k)==0)  cycle
  do s=1, nC(k)
    if (ns(s,k)==1) then
      skipCluster(s,k) = .TRUE.
    else if (ALL(Parent(SibID(1:nS(s,k), s, k), 3-k) < 0)) then
      call getFSpar(s, k, .TRUE.,par)
      if (par < 0) then
        if (nS(-par, 3-k) == nS(s,k) .and. DumClone(s,k)==0) then  !cannot tell if mat or pat
          skipCluster(s,k) = .TRUE.
        endif
      endif          
    endif
  enddo
enddo  
! done in FsibsGPs()                    

do k=1,2
  call getBYrank_c(k, BYrankC(:,k))
enddo

not4 = (/1,2,3,5,6/)
do x=1, MAXVAL(nC)
  if (MODULO(x,10)==0)  call rchkusr()
  if (MAXVAL(nC) > 20) then
    if (MODULO(x,10)==0 .and. quiet==-1)  call Rprint("", (/x/), (/0D0/), "INT") 
  endif
  
  do k=1,2
    if (x > nC(k))  cycle
    s = BYrankC(x,k)    
    if (ALL(GpID(:,s,k)/=0) .and. .not. IsNewSibship(s,k)) cycle  
    if (skipCluster(s,k) .and. .not. (hermaphrodites==2 .and. k==1))  cycle

    nCG = 0 
    CandGP = 0
    CurGP = GpID(:,s,k)
    do m=1,2
      if (GpID(m,s,k)/=0) then
        nCG(m) = 1
        CandGP(1,m) = GpID(m,s,k)
      endif
      call setParTmp(-s, k, 0, m)
    enddo
    call setEstBY(-s,k)
    
    do i=1,nInd
      if (Parent(i,k)==-s) cycle
      if (ANY(CandGP==i)) cycle
      if (Sex(i)<3) then      
        if (nCG(Sex(i))==mxCP)  cycle
        ! if (ANY(curGP/=0) .and. hermaphrodites/=1) then  ! take curGP for gospel
          ! if (curGP(Sex(i)) > 0) cycle
          ! if (curGP(Sex(i))<0) then
            ! if (ns(-curGP(Sex(i)), Sex(i)) > 1) cycle
          ! endif  
        ! endif
      else
        if (ANY(nCG==mxCP))  cycle 
      endif
      if (DoMtDif) then
        if (k==1 .and. Sex(i)==1 .and. mtDif(SibID(1,s,k), i))  cycle   
      endif      
      Maybe = .TRUE.
      do n=1,ns(s,k)
        if (AgeDiff(SibID(n,s,k), i) <= 1) then
          Maybe = .FALSE.
          exit
        endif
      enddo
      if (.not. Maybe)  cycle
      ALRtmp = missing                
      call CalcAgeLR(-s,k, i,Sex(i), 0,1, .TRUE., ALRtmp(1))
      if (ALRtmp(1) == impossible .or. ALRtmp(1) < 3.0*TF) cycle
      call CalcAgeLR(i,Sex(i), -s,k, 0,1, .TRUE., ALRtmp(2))
      if (ALRtmp(2)/=impossible .and. (ALRtmp(1)-ALRtmp(2)) < TF)  cycle
      call ChkAncest(i,0,-s,k, AncOK)
      if (.not. AncOK)  cycle
      call QGP(i, Sex(i), s, k, LRG) 
      if (LRG < TF * MAX(dble(nS(s,k)),2D0))  cycle      ! 2TF for ns=1      
      LLA = missing             
      call GPfilter(i,s,k,LLA)
      if (LLA(4)>0 .or. (LLA(4) - LLA(7) < TA .and. Complx/=0) .or. & 
       (LLA(4) - LLA(7) < -TA .and. Complx==0) .or. &
       (LLA(4) - MaxLL(LLA(not4)) < 2*TF .and. ANY(LLA(not4) < 0))) cycle  
      do m=1,2
        if (Sex(i)<3 .and. Sex(i)/=m) cycle
        if (ncG(m) < mxCP) then  ! arbitrary threshold to limit comp. time
          ncG(m) = nCG(m) + 1
          CandGP(nCG(m), m) = i
        endif
!        if (Complx==0 .and. Mate(i)/=0 .and. ncG(3-m)<mxCP) then
!          ncG(3-m) = ncG(3-m) +1
!          CandGP(nCG(3-m), 3-m) = Mate(i)
!        endif   ! mate may not be valid
      enddo
    enddo
    
    do m=1,2
      do r=1, nC(m) 
        if (ncG(m) == mxCP) exit
        if (m==k .and. s==r) cycle
        if (any(CandGP(:,m) == -r)) cycle  ! current GP
        call ChkAncest(-r,m, -s,k, AncOK)
        if (.not. AncOK)  cycle
        if (nS(r,m)==1 .and. ANY(SibID(1:nS(s,k),s,k) == SibID(1,r,m))) cycle
        if (DoMtDif) then
          if (k==1 .and. m==1 .and. mtDif(SibID(1,s,k), SibID(1,r,m)))  cycle 
        endif
        if (m/=k .and. complx<2) then
          if (ALL(Parent(SibID(1:ns(s,k),s,k),m) == -r))  cycle
          if (ALL(Parent(SibID(1:ns(r,m),r,m),k) == -s))  cycle
        endif
        call CalcAgeLR(-s,k, -r,m, 0,1, .TRUE., ALRtmp(1))
        if (ALRtmp(1) == impossible .or. ALRtmp(1) < 3.0*TF) cycle
        call CalcAgeLR(-r,m, -s,k, 0,1, .TRUE., ALRtmp(2))
        if (ALRtmp(2)/=impossible .and. (ALRtmp(1)-ALRtmp(2)) < TF)  cycle 
        if (ALL(ABS(ALRtmp) < 0.001))  cycle  ! no age info
        if (hermaphrodites==0) then
          call QGP(-r, m, s, k,  LRG) 
          if (LRG < TF*dble(MIN(nS(s,k), nS(r,m)))) cycle  ! conservative.
        endif

        LLX = missing
        call PairUA(-s, -r, k, m, LLX(1,1))
        if (LLX(1,1)>0) cycle
        call CalcU(-s,k, -r,m, LLX(1,2)) 
        if ((LLX(1,1) - LLX(1,2)) < nS(s,k)*TF) cycle
        call addFS(0, r, m, s, k, LLX(2,1), ix, dx) 
        if ((MaxLL(LLX(:,1)) - LLX(1,2)) < TA)  cycle
        if (ncG(m) < mxCP) then
          nCG(m) = nCG(m) + 1
          CandGP(nCG(m), m) = -r
        endif
      enddo
    enddo
           
    if (ALL(nCG <=1) .and. ALL(CandGP(1,:) == curGP) .and. .not. &
      (any(curGP==0) .and. any(curGP<0))) then
      do m=1,2
        call setParTmp(-s, k, curGP(m), m)  
      enddo
      call setEstBY(-s,k)
      cycle
    endif 

    call SelectParent(-s, k, nCG, candGP, .FALSE., .FALSE.)
    
    if (any(GpID(:,s,k)/= curGP)) then
      if (Complx == 0) call UpdateMate(-s, k, curGP, .FALSE.)
      if (any(GpID(:,s,k) /= curGP .and. curGP/=0)) then
        ToCheck(SibID(1:ns(s,k), s, k)) = .TRUE.
      endif
    endif
    
    if (ALL(GpID(:,s,k)==0) .and. nS(s,k)==1) then  ! single sib left; remove sibship 
      i = SibID(1,s,k)
      call CheckDropSibship(s, k, DropS)
      if (DropS) then
        call getBYrank_c(k, BYrankC(:,k))
        do r=s, nC(k)
          skipCluster(r,k) = skipCluster(r+1,k)
        enddo
        if (Parent(i,3-k) < 0) then
          par = Parent(i,3-k)
          call CheckDropSibship(-Parent(i,3-k), 3-k, DropS)
          if (DropS) then
            call getBYrank_c(k, BYrankC(:,3-k))
            do r=-par, nC(3-k)
              skipCluster(r,3-k) = skipCluster(r+1,k)
            enddo
          endif
        endif
      endif
      cycle
    endif
  enddo  ! k
enddo  ! x

end subroutine SibGrandparents

! #####################################################################

subroutine GPfilter(A, SB, k, LLg)
use Global
implicit none

integer, intent(IN) :: A, SB, k
double precision, intent(OUT) :: LLg(7)
integer :: fsi, sibfcl, kA
double precision :: ALR, dx(maxSibSize), LLtmp

LLg = missing
call AddGP(A, SB, k, LLg(4))
if (LLg(4) > 0)  return
! U
call CalcU(A, k, -SB, k, LLg(7))  
if (LLg(4) - LLg(7) < TA)  return
! GGP / 3rd degree rel   ! before or after FA?
call AddGGP(A, SB, k, LLg(6))  
if (LLg(4) - LLg(6) < TF .and. LLg(6)<0)  return
! FA
if (any(Parent(A,:)/=0)) then
  call CalcAgeLR(-SB,k, A,Sex(A), 0,2, .TRUE., ALR)
  if (ALR /= impossible)  call pairUA(-SB, A, k, 3, LLg(5))
  if (LLg(4) - LLg(5) < TF .and. LLg(5)<0) then
    if (Sex(A) < 3) then
      kA = Sex(A)
    else
      kA = 1
    endif    
    if (Parent(A,3-kA)/=0 .and. GpID(3-kA,SB,k)==0) then
      call setParTmp(-SB,k, Parent(A,3-kA),3-kA)
      call AddGP(A, SB, k, LLtmp)
      call setParTmp(-SB,k, 0,3-kA)
      if (LLtmp < 0 .and. LLtmp > LLg(4)) then
        LLg(4) = LLtmp
        if (LLg(4) - LLg(5) < TF)  return
      endif
    else
      return
    endif
  endif
endif

! FS/HS  
call getfocal(A, 0, SB, k, sibfcl)
if (sibfcl == 2) then
  call AddFS(A, SB, k,0,k, LLg(2), fsi, dx)
else if (any(Parent(A,:)/=0)) then
  call AddSib(A, SB, k, LLg(3))
endif 

end subroutine GPfilter

! #####################################################################

subroutine SibChk(A, SB, k, focal, cat, LLg)  ! 1=filter, 2=confirm
use Global
implicit none

integer, intent(IN) :: A, SB, k, focal, cat  
double precision, intent(OUT) :: LLg(7)
integer :: fsi
double precision :: ALR, dx(maxSibSize), Threshold

if (cat==1) then
  Threshold = TF
else if (cat==2) then
  Threshold = TA
else
  Threshold = missing
  call Erstop("SibChk: cat must be 1 or 2", .TRUE.)
endif

LLg = missing
! FS/HS 
if (focal==2 .or. Complx==0) then
  call AddFS(A, SB, k,0,k, LLg(2), fsi, dx)  ! includes age check
else if (focal==3) then
  call AddSib(A, SB, k, LLg(3))
else
  call Erstop("SibChk: focal must be 2 or 3", .TRUE.)
endif

if (all(LLg > 0))  return
! U
if (cat==1) then
  call CalcU(A, k, -SB, k, LLg(7))  
  if (LLg(focal) - LLg(7) < TA) then
    if (focal==3 .and. Parent(A,3-k)==0 .and. LLg(focal)-LLg(7) > TF) then
      call AddFS(A, SB, k,0,k, LLg(2), fsi, dx)  
      if (LLg(2) - LLg(7) < TA)  return
    else
      return
    endif
  endif
endif

if (focal==3 .and. all(Parent(A,:)==0) .and. cat==1) then  ! 2nd rels indistinguishable
  LLg(4:5) = LLg(3)
  return
endif
! FA
call CalcAgeLR(A,Sex(A), -SB,k, 3,4, .TRUE., ALR) 
if (ALR /= impossible)  call addFA(A, SB, k, LLg(5))
if (LLg(focal) - LLg(5) < Threshold .and. LLg(5)<0)  return
! GP
call AddGP(A, SB, k, LLg(4))
if (LLg(focal) - LLg(4) < Threshold .and. LLg(4)<0)  return
! HA / 3rd degree rel
if (cat==2)  call pairUA(A, -SB, k, k, LLg(6))   !  .and. focal==3
!if (LLg(focal) - LLg(6) < Threshold .and. LLg(6)<0)  return

end subroutine SibChk

! #####################################################################

subroutine Calc4U(Par, B, kB,  A, kA, LLU, LLcor)  
use Global
implicit none

integer, intent(IN) :: Par(2), B, kB,  A, kA
double precision, intent(OUT) :: LLU(4), LLcor(3,2)
integer :: m, y, CY(4), kY(4), x, v, ParA(2)
double precision :: LLtmp(3), LLoverlap(4,4)
logical :: ConPar(4,4)

CY = (/ Par, B, A /)
kY = (/ 1, 2, kB, kA /)

ParA = getPar(A, kA)
do m=1,2
  if (ParA(m)==0)  cycle
  call setParTmp(A, kA, 0, m)   
enddo

LLU = 0D0  
! Individual likelihoods, excluding overlap with A if any
LLcor = 0D0
! correction factors: LL(all 4 indiv) = e.g. LL(A+B) from CheckRel + LLcor(3,m)
LLoverlap = 0D0
! likelihoods of overlaps/ change in likelihoods when considering vs ignoring connections
LLtmp = missing   
                
call CalcU(A, kA, 0, 0, LLU(4))
do y=1,3
  if (CY(y)==0) cycle
  call CalcU(CY(y),kY(y), A,kA, LLtmp(1))
  LLU(y) = LLtmp(1) - LLU(4)
enddo  
! pairs likelihoods, if no overlap
do m=1,2
  LLcor(m,m) = LLU(3-m) + LLU(3) 
enddo
LLcor(3,:) = LLU(1) + LLU(2)

ConPar = .FALSE.
if (ANY(CY(1:3)<0)) then
  do m=1,2
    if (Par(m)==0) cycle
    call Connected(Par(m), m, A, kA, ConPar(4,m))
    if (B/=0)  call Connected(Par(m), m, B, kB, ConPar(3,m))
  enddo
  call Connected(Par(1), 1, Par(2), 2, ConPar(2,1))
  if (B/=0)  call Connected(B, kB, A, kA, ConPar(4,3))
endif

if (ANY(ConPar)) then 
  ! likelihoods of overlaps/ change in likelihoods when considering vs ignoring connections
  do y=1,3
    do x=y+1, 4
      if (.not. ConPar(x,y))  cycle
      call CalcU(CY(x),kY(x), CY(y),kY(y), LLtmp(1))
      call CalcU(CY(x),kY(x), 0,0, LLtmp(2))
      call CalcU(CY(y),kY(y), 0,0, LLtmp(3))
      LLoverlap(x,y) = LLtmp(1) - LLtmp(2) - LLtmp(3)
      LLoverlap(y,x) = LLoverlap(x,y)
    enddo
  enddo

  do m=1,2        
    do y=1,3  ! focal
      if (y/=m .and. y/=3) cycle           
      if (y==1) then
        call CalcU(CY(2),kY(2), CY(3),kY(3), LLcor(y,m))
      else if (y==2) then
        call CalcU(CY(1),kY(1), CY(3),kY(3), LLcor(y,m))    
      else if (y==3) then
        call CalcU(CY(1),kY(1), CY(2),kY(2), LLcor(y,m))
      endif
      do x=1,3
        if (ConPar(4,x) .and. x/=y) then                                     
          LLcor(y,m) = LLcor(y,m) + LLoverlap(x,4)
        endif
        do v=1,2
          if (ConPar(x,v) .and. (x==y .or. v==y)) then                                   
            LLcor(y,m) = LLcor(y,m) + LLoverlap(x,v) 
          endif
        enddo
      enddo
    enddo
  enddo
endif

do m=1,2
  call setParTmp(A, kA, parA(m), m)
enddo

! if (A==2668 .and. Par(1)==2311 .and. Par(2)<0)  then
  ! open (unit=42,file='log.txt',status="unknown", position="append")
  ! write(42,*) ""
  ! write(42,*) "Calc4U ConPar:"
  ! write(42,*) A, kA, Par
  ! do x=1,4
    ! write(42,*)  ConPar(x,:)
  ! enddo
  ! write(42,*) ""
  ! close(42)
! endif

end subroutine Calc4U

! #####################################################################

subroutine GGpairs !(ExtraAge)  ! find & assign grandparents of singletons
use Global
implicit none

!logical, intent(IN) :: ExtraAge
integer :: i, j, k, nCG(2,2), CandG(2,mxCP, 2), n, s, BYrank(nInd), x,t
double precision :: LRS, LRG, ALR, ALRx(2), LRx, LLx(7,2), LL(7,2)
logical :: AncOK, MaybePair

call getRank_i(BYrank)
t=1
do x=1, nInd
  if (MODULO(x,200)==0)  call rchkusr()
  if (quiet==-1 .and. any(viginti==x)) call print_progress(x,t)

  i = BYRank(x)
  if (skip(i))  cycle
  if (ALL(Parent(i,:)==0) .and. AgePhase <2 .and. hermaphrodites/=2) cycle  
  ! can't determine if mat or pat GP  ! TODO: more nuanced.
  if (ALL(Parent(i,:)/=0)) cycle
  nCG = 0  
  CandG = 0
  LL = missing
  
  do k=1,2
    if (Parent(i,k)/=0) cycle
    do j=1, nInd
      if (i==j)  cycle
      if (ANY(nCG(k,:)>=mxCP)) cycle
      if (AgeDiff(i,j) <= 1)  cycle
      if (any(parent(j,:) == i))  cycle
      call ChkAncest(j, sex(j), i, sex(i), AncOK)
      if (.not. AncOK)  cycle
      call PairQHS(i, j, LRS) 
      if (LRS < TF)  cycle
      call LRGG(i,k,j,Sex(j),LRG)
      if (LRG < -TA)  cycle
      call CalcAgeLR(i,Sex(i), j,Sex(j), k,4,.TRUE.,ALR)
      if (ALR==impossible .or. ALR < 2.0*TF)  cycle  ! 3*
      if (DoMtDif) then
        if (k==1 .and. Sex(j)==1 .and. mtDif(i,j))  cycle   
      endif
      if (Parent(i,3-k)<= 0 .and. hermaphrodites/=2) then      
        call CalcAgeLR(i,Sex(i), j,Sex(j), 3-k,4,.TRUE.,ALRx(1))
        if (ALRx(1)/=impossible .and. (ALR - ALRx(1)) < TA) then     
          if (Parent(i,3-k)==0) then   !  .and. .not. DoMtDif   !! it's complicated
            cycle   ! unclear if pat. or mat. GP     
          else if (Parent(i,3-k) < 0) then
            call CalcAgeLR(Parent(i,3-k),3-k, j,Sex(j), 0,1, .TRUE., ALRx(2))
            if (ALRx(2)/=impossible .and. (ALR - ALRx(2)) < TA) then            
              call QGP(j, sex(j), -Parent(i,3-k),3-k, LRx)
              if (LRx > TF*ns(-Parent(i,3-k),3-k)) then
                LLx = missing                                       
                call CheckAdd(j, -Parent(i,3-k),3-k, 4, LLx(:,1), LLx(:,2))
                 if ((LLx(4,1)- MaxLL(LLx((/1,2,6,7/),1))) > TF) then
                  cycle   ! plausible that GP of opposing sibship
                endif
              endif
            endif
          endif
        endif
      endif
      
      MaybePair = .TRUE.
      LL = missing                  
      call CheckPair(i,j,k, 4, LL(:,1), LL(:,2))
      do n=1,2
        if (AgePhase==0 .and. n==2)  cycle
        if (AgePhase==2 .and. n==1)  cycle
        if (LL(4,n)<0 .and. (LL(4,n)- MaxLL(LL((/1,2,6,7/),n))) > TF) then   ! TODO: CHECK, WAS n*TF
          MaybePair = .TRUE.
        else
          MaybePair = .FALSE.
          exit
        endif
      enddo
      if (.not. MaybePair) cycle
      do n=1,2
        if (Sex(j)/=n .and. sex(j)<3)  cycle
        if (nCG(k,n) == mxCP)  cycle
        nCG(k, n) = nCG(k, n) +1
        CandG(k, nCG(k,n), n) = j
      enddo
    enddo
  enddo
    
  if (ANY(nCG>0)) then
    do k=1,2
      if (ANY(nCG(k,:)>0)) then
        call NewSibship(i, 0, k)
        s = nC(k)
        call SelectParent(-s, k, nCG(k,:), CandG(k,:,:), .FALSE., .FALSE.)
        if (ALL(GpID(:,s,k)==0)) then  
          call RemoveSib(i, s, k) 
          call DoMerge(0, s, k)
        else if (Complx == 0) then
          call UpdateMate(-s, k, (/0,0/), .FALSE.)
        endif
      endif
    enddo
  endif  
  
  do k=1,2
    if (nc(k)==0)  cycle
    if (nS(nc(k), k) == 0 .or. SibID(1, nC(k), k) == 0) then
      call Erstop("grandparent pairs -- empty sibship!", .TRUE.)
    endif
  enddo 
enddo

end subroutine GGpairs

! ##############################################################################

subroutine FsibsGPs
! assign grandparents to full-sibling clusters
use Global
implicit none

integer :: x, j, fsx(maxSibSize), k,m,r, nCG(2,2), candGP(mxCP,2,2), SAB(2), ix, n,t
logical :: maybeGP(2,2), AncOK
double precision :: ALR, LRG, LLX(2,2), dx(maxSibSize), LRS

if (.not. (DoMtDif .or. any(AgePriorA(:,1,2) /= AgePriorA(:,1,3)) .or. &
  any(AgePriorA(:,2,2) /= AgePriorA(:,2,3)))) then   ! diff ageprior between MGM-PGM or MGP-PGP
!  print *, 'Not doing FsibsGPs'
  return   ! no way to distinguish between maternal and paternal grandparents
endif
t=1
do x=1, nInd
  if (quiet==-1 .and. any(viginti==x)) call print_progress(x,t)
  if (nFS(x) ==0 )  cycle  ! not 'primary' sib of FS cluster
  if (.not. ALL(Parent(x,:) < 0) .or. ALL(parent(x,:)==0))  cycle
  if (all(Parent(x,:)==0) .and. BY(x) < 0)  cycle   ! high risk wrong way around
  SAB = -parent(x,:)
  if (all(SAB/=0)) then
    if (ns(SAB(1),1) /= ns(SAB(2),2))  cycle  ! resolvable via SibGrandparent 
    if (any(GpID(:,SAB(1),1)/=0) .or. any(GpID(:,SAB(2),2)/=0))  cycle   ! resolvable via SibGrandparent 
!    if (ns(SAB(1),1)==1 .and. any(GpID(:,SAB(1),1)/=0) .and. any(GpID(:,SAB(2),2)/=0))  cycle
    if (all(BY(SibID(1:nS(SAB(1),1),SAB(1),1)) < 0))  cycle  ! high risk wrong way around + slow
  endif
 
  if (ALL(parent(x,:)==0)) then
    do k=1,2
      call NewSibship(x, 0, k)
    enddo
    SAB = -parent(x,:)
  endif
  
  fsx = 0
  fsx(1:nFS(x)) = FSID(1:nFS(x),x)
  nCG = 0  
  CandGP = 0
  ! TODO: something with current GPs?
  
  do j=1, nInd
    maybeGP = .TRUE.
    if (ANY(fsx == j))  cycle
    do n=1,nFS(x)
      if (AgeDiff(fsx(n),j) <= 1) then
        maybeGP = .FALSE.
        exit
      endif
    enddo
    if (.not. any(MaybeGP))  cycle
    do m=1,2
      if (m/=Sex(j) .and. Sex(j)<3)  maybeGP(m,:) = .FALSE.
    enddo
    if (Sex(j)/=2 .and. DoMtDif) then
      if (mtDif(x,j))  maybeGP(1,1) = .FALSE.   ! mat grandmother must have same mt haplo
    endif 
    if (.not. ANY(maybeGP))  cycle          
    do k=1,2
      call ChkAncest(j, sex(j), -SAB(k), k, AncOK)
      if (.not. AncOK)  maybeGP(:,k) = .FALSE.
    enddo 
    if (.not. ANY(maybeGP))  cycle 
    do k=1,2
      do m=1,2
        if (.not. maybeGP(m,k))  cycle
        call CalcAgeLR(-SAB(k), k, j,m, k,1,.TRUE.,ALR)
        if (ALR==impossible .or. ALR < 3.0*TF) then
          maybeGP(m,k) = .FALSE.
        endif
      enddo
    enddo
    if (.not. ANY(maybeGP))  cycle    
    if (ns(SAB(1),1)>1) then
      call QFSGP(j, Sex(j), SAB(1), 1, LRG)   ! TODO: sep subroutine for FS?
      if (LRG < TF)  cycle
    else
      call PairQHS(x, j, LRS)    
      if (LRS < TF)  cycle
      call LRGG(x,1,j,Sex(j),LRG)
      if (LRG < TA)  cycle
    endif
    
    do k=1,2
      do m=1,2
        if (nCG(m,k) == mxCP)  cycle
        if (maybeGP(m,k)) then
          nCG(m,k) = nCG(m,k) +1
          CandGP(nCG(m,k),m,k) = j
        endif
      enddo
    enddo
  enddo
  
  ! dummy GPs
  do m=1,2
    if (ns(SAB(m),m)==1)  exit
    do r=1, nC(m) 
      if (r==SAB(m)) cycle
      if (ANY(GpID(:,r,m) == SAB)) cycle
      maybeGP = .FALSE. 
      do k=1,2       
        if (ncG(m,k) == mxCP) cycle
        if (nS(r,m)==1 .and. ANY(SibID(1:nS(SAB(k),k),SAB(k),k) == SibID(1,r,m))) cycle
        if (DoMtDif) then
          if (k==1 .and. m==1 .and. mtDif(SibID(1,SAB(k),k), SibID(1,r,m)))  cycle
        endif
        if (m/=k .and. complx<2) then
          if (ALL(Parent(SibID(1:ns(SAB(k),k),SAB(k),k),m) == -r))  cycle
          if (ALL(Parent(SibID(1:ns(r,m),r,m),k) == -SAB(k)))  cycle
        endif
        call ChkAncest(-r,m, -SAB(k),k, AncOK)
        if (.not. AncOK)  cycle
        call CalcAgeLR(-SAB(k),k, -r,m, 0,1, .TRUE., ALR)
        if (ALR == impossible .or. ALR < 3.0*TF) cycle        
        maybeGP(m,k) = .TRUE.
      enddo
      if (.not. ANY(maybeGP(m,:)))  cycle 
      if (hermaphrodites==0) then
        call QFSGP(-r, m, SAB(1), 1,  LRG) 
        if (LRG < TF*dble(MIN(nS(SAB(1),1), nS(r,m)))) cycle  ! conservative.   
      endif          
      LLX = missing
      call PairUA(-SAB(1), -r, 1, m, LLX(1,1))
      if (LLX(1,1)>0) cycle
      call CalcU(-SAB(1),1, -r,m, LLX(1,2)) 
      if ((LLX(1,1) - LLX(1,2)) < nS(SAB(1),1)*TF) cycle
      call addFS(0, r, m, SAB(1), 1, LLX(2,1), ix, dx) 
      if ((MaxLL(LLX(:,1)) - LLX(1,2)) < TA)  cycle
      
      do k=1,2
        if (nCG(m,k) == mxCP)  cycle
        if (maybeGP(m,k)) then
          nCG(m,k) = nCG(m,k) +1
          CandGP(nCG(m,k),m,k) = -r
        endif
      enddo
    enddo
  enddo
    
  if (all(CandGP(:,1,1) == CandGP(:,1,2)) .and. all(CandGP(:,2,1) == CandGP(:,2,2))) then
    ! identical candidates for mat & pat side of full sibship, incl all 0
    do k=1,2
      if (ns(SAB(k),k)==1 .and. all(GpID(:,SAB(k),k)==0)) then
        call RemoveSib(x, SAB(k),k)
        call DoMerge(0, SAB(k),k)
      endif
    enddo 
    cycle   
  endif

  if (ANY(nCG>0)) then   
    do k=1,2  !2,1,-1
      if (all(nCG(:,k)==0))  cycle
      call SelectParent(-SAB(k), k, nCG(:,k), CandGP(:,:,k), .FALSE., .FALSE.)
    enddo
    call ChkGPs(SAB, CandGP)  ! drops GPs if could be GPs of 3-k
    
    do k=2,1,-1
      if (all(nCG(:,k)==0))  cycle
      if (any(GpID(:,SAB(k),k) /=0))  cycle
      call SelectParent(-SAB(k), k, nCG(:,k), CandGP(:,:,k), .FALSE., .FALSE.)
    enddo
    call ChkGPs(SAB, CandGP)
       
    do k=1,2      
      if (ns(SAB(k),k)==1 .and. all(GpID(:,SAB(k),k)==0)) then
        call RemoveSib(x, SAB(k),k)
        call DoMerge(0, SAB(k),k)
      endif
    enddo
  endif    
enddo

end subroutine FsibsGPs

! ##############################################################################

subroutine ChkGPs(SAB, CandGP)
use Global
implicit none

integer, intent(IN) :: SAB(2), candGP(mxCP,2,2)
integer :: k,m
logical :: MaybeOpp(2,2)

! if both pairs could have been assigned as GP via parent 3-k, drop all
maybeOpp = .TRUE.
do k=1,2      
  do m=1,2
    if (GpID(m,SAB(k),k)/=0) then
      if (ANY(CandGP(:,m,3-k) == GpID(m,SAB(k),k))) then
        maybeOpp(k,m) = .TRUE.
      else
        maybeOpp(k,m) = .FALSE.
      endif
    endif
  enddo
enddo 
do k=1,2
  if (ALL(maybeOpp)) then
    do m=1,2
      call setPar(-SAB(k),k, 0,m)
    enddo
  endif
enddo

end subroutine ChkGPs

! ##############################################################################

subroutine Qadd(A, SB, kB, LR)
use Global
implicit none

integer, intent(IN) :: A, SB, kB
double precision, intent(OUT) :: LR
integer :: l, x
double precision :: PrL(nSnp), PrX(3)

PrL = 0D0
do l=1,nSnp
  do x=1,3
    PrX(x) = OKAP(Genos(l,A), x, l) * DumP(x,l,SB,kB) / AHWE(x,l)
  enddo   ! simple LL identical for HS and GP
  PrL(l) = LOG10(SUM(PrX))
enddo
LR = SUM(PrL)

end subroutine Qadd

! #####################################################################

subroutine QGP(A, kA, SB, kB, LR)  ! A indiv or dummy, GP of SB
use Global
implicit none

integer, intent(IN) :: A, kA, SB, kB
double precision, intent(OUT) :: LR
integer :: l, x
double precision :: PrLR(nSnp), PrX(3,2), PrA(3)

if (ns(SB,kB)==1 .and. A>0) then
  call PairQHS(SibID(1,SB,kB), A, LR)
else
  PrLR = 0D0
  do l=1,nSnp
    call ParProb(l, A, kA, 0, 0, PrA)  ! no effect on time vs. LindG/DumP 1x
    do x=1,3               
      PrX(x,1) =XPr(1,x,l,SB,kB) * SUM(AKAP(x,:,l) * PrA)
      PrX(x,2) =XPr(1,x,l,SB,kB) * AHWE(x,l)
    enddo
    PrLR(l) = LOG10(SUM(PrX(:,1))) - LOG10(SUM(PrX(:,2))) 
  enddo
  LR = SUM(PrLR)
endif

end subroutine QGP

! #####################################################################

subroutine QPO(A, SB, kB, LR)  ! A replaces dummy SB?
use Global
implicit none

integer, intent(IN) :: A, SB, kB
double precision, intent(OUT) :: LR
integer :: l, x, sib1
double precision :: PrLR(nSnp), PrX(3,2), LL(2), PrA(3)

if (ns(SB,kB)==1) then
  sib1 = SibID(1,SB,kB)
  call CalcU(sib1,kB,A,kB, LL(1))
  call PairPO(sib1, A, kB, 1, LL(2))
  LR = LL(2) - LL(1)
else
  PrLR = 0D0
  do l=1,nSnp
    call ParProb(l, A, kB, 0, 0, PrA)
    do x=1,3
      PrX(x,1) = XPr(1,x,l,SB,kB) * XPr(2,x,l,SB,kB)
      PrX(x,2) = XPr(1,x,l,SB,kB) * PrA(x)
    enddo
    PrLR(l) = LOG10(SUM(PrX(:,2))) - LOG10(SUM(PrX(:,1)))
  enddo
  LR = SUM(PrLR)
endif

end subroutine QPO

! #####################################################################

subroutine QFS(A, SB, kB, LR)   ! only when Complx==0
use Global
implicit none

integer, intent(IN) :: A, SB, kB
double precision, intent(OUT) :: LR
integer :: l, x, y
double precision :: PrLR(nSnp), PrXY(3,3,2), PrY(3)

PrLR = 0D0
do l=1,nSnp
  call ParProb(l, Parent(SibID(1,SB,kB), 3-kB), 3-kB, -1, 0, PrY)
  do x=1,3
    do y=1,3
      PrXY(x,y,1) = OKA2P(Genos(l,A),x,y) * DumP(x,l,SB,kB) * PrY(y)
      PrXY(x,y,2) = OKA2P(Genos(l,A),x,y) * AHWE(x,l) * AHWE(y,l)
    enddo
  enddo  
  PrLR(l) = LOG10(SUM(PrXY(:,:,1))) - LOG10(SUM(PrXY(:,:,2)))
enddo
LR = SUM(PrLR)

end subroutine QFS

! #####################################################################

subroutine QFSGP(A, kA, SB, kB, LR)  ! A indiv or dummy, GP of SB, all B's are FS
use Global
use CalcLik
implicit none

integer, intent(IN) :: A, kA, SB, kB
double precision, intent(OUT) :: LR
integer :: l, x,y, i
double precision :: PrLR(nSnp), PrXY(3,3,2), PrA(3), PrI(3,3) !, PrY(3)

i = FSID(maxSibSize+1, SibID(1,SB,kB))

if (ns(SB,kB)==1 .and. A>0) then
  call PairQHS(SibID(1,SB,kB), A, LR)
else
  PrLR = 0D0
  do l=1,nSnp
    call ParProb(l, A, kA, 0, 0, PrA)
    PrI = FSLik(l,i)
!    call ParProb(l, Parent(i,3-kB),3-kB,-1,0, PrY)  ! GPs only
    do x=1,3
      do y=1,3
        PrXY(x,y,1) = PrI(x,y) * SUM(AKAP(x,:,l) * PrA) * AHWE(y,l)
        PrXY(x,y,2) = PrI(x,y) * AHWE(x,l) * AHWE(y,l)
      enddo
    enddo
    PrLR(l) = LOG10(SUM(PrXY(:,:,1))) - LOG10(SUM(PrXY(:,:,2))) 
  enddo
  LR = SUM(PrLR)
endif

end subroutine QFSGP

! #####################################################################

subroutine CheckRel(A, kA, B, kB, focalIN, LLg, LL)
use Global
implicit none

integer, intent(IN) :: A, kA, B, kB, focalIN
double precision, intent(OUT) :: LLg(7), LL(7)
logical:: FSJ  !do separately?
integer :: k, focal

focal = focalIN
FSJ = .FALSE.
LLg = missing
LL = missing
if (A==0 .or. B==0) then
  call Erstop("CheckRel A or B null ", .TRUE.)
else if (A==B .and. (A>0 .or. kA==kB)) then
  call Erstop("CheckRel A==B ", .TRUE.)
else if (A > 0 .and. B > 0) then
  if (kA == 0 .and. kB==0) then
    call Erstop("CheckRel kA == kB == 0!", .TRUE.)
!  else if (kA /= 0 .and. kB/=0 .and. kA/=kB .and. &
!    focalIN/=1 .and. focalIN/=4) then
!    call Erstop("CheckRel kA /= kB!", .TRUE.)
  else if (kB /= 0) then ! .or. focalIN==1 .or. focalIN==4) then
    k = kB
  else if (kA /= 0) then
    k = kA
  endif
  call CheckPair(A, B, k, focal, LLg, LL)  
else if (A > 0 .and. B < 0) then
  if (kB<1 .or. kB>2)  call Erstop( "CheckRel A>0, B<0, invalid kB", .TRUE.)
  if (focal==0)  call Erstop("CheckRel focal == 0!", .TRUE.)
  if (focalIN==1) focal = 3  ! -B parent of A -> B's HS of A
  call CheckAdd(A, -B, kB, focal, LLg, LL)
  if (focalIN==1 .or. focalIN==6) then
    if (Parent(A,3-kB)==0 .and. Complx/=0) then  ! called by CalcCandParLL, want single vs parent-pair
      LLg(2) = 333D0
      LL(2) = 333D0
    endif
    call ReOrderAdd(LLg)
    call ReOrderAdd(LL) 
  endif
else if (A < 0 .and. B > 0) then
  if (kA<1 .or. kA>2)  call Erstop("CheckRel A<0, B>0, invalid kA", .TRUE.)
  call CheckAdd(B, -A, kA, focal, LLg, LL)
else if (A < 0 .and. B < 0) then
  if (kA<1 .or. kA>2)  call Erstop("CheckRel A<0, B<0, invalid kA", .TRUE.)
  if (kB<1 .or. kB>2)  call Erstop( "CheckRel A<0, B<0, invalid kB", .TRUE.)
  ! note: focal=1: merge, focal=4: SB parent of SA. FSJ: full-sib merge
  call CheckMerge(-A, -B, kA, kB, focal, LLg, LL, FSJ)
endif

end subroutine CheckRel

! #####################################################################

subroutine ReOrderAdd(LL)  
! reorder output from CheckAdd for compatibility with CheckPair (for POZ)
use Global
implicit none

double precision, intent(INOUT) :: LL(7)
double precision :: LLtmp(7)

LLtmp = missing
LLtmp(1) = MaxLL(LL(2:3))
LLtmp(2:3) = LL(5:6)
LLtmp(5) = LL(4)   ! ? not technically correct, but ... 
LLtmp(6) = LL(1)
LLtmp(7) = LL(7) 

LL = LLtmp

end subroutine ReOrderAdd

! #####################################################################

subroutine CheckAdd(A, SB, k, focal, LLg, LL)
use Global
use CalcLik
implicit none

integer, intent(IN) :: A, SB, k, focal
double precision, intent(OUT) :: LLg(7), LL(7)
double precision :: LRHS, ALR(7), LLPH(2), ALRH(2), LLAU(2,3), ALRAU(2,3), LLUi, &
  LLC, LLz(7), ALRz(7), LLM(3), LLp(7), LLpg(7), LLFH(3), LLPX(2,2), dx(maxSibSize), &
  ALRq, LLHH(4,2), ALRtmp, LHH(3), LHH2, LLy(2,2), LLpo(ns(SB,k),2), ALRpo(ns(SB,k),2), &
  LLgp(ns(SB,k),3), ALRgp(ns(SB,k),3), LLfs(3,2)!, LLXi!, LLPHS(2,2), ALRPHS(2,2)
integer :: x, y, FSPar, i, ParTmp(2), OpPar(maxSibSize), nop, fsi, ix, m, Bi, sib1, curpar(2)                               
logical :: AncOK, fclsib, MaybeOpp, ParOK 

LL = missing
LLg = missing
ALR = missing            
  
! quick check
LRHS = missing             
call Qadd(A, SB, k, LRHS)  ! 2nd degree relatives vs unrelated 
if (LRHS < MIN(TF*2, TF*nS(SB,k)) .and. (focal/=4 .and. focal/=7 .and. focal/=6)) return
  
if (Sex(A)<3 .and. Sex(A)/=k) then
  LL(1) = impossible
  if (focal==1)  return
endif
  
if (focal==1) then                
  call CalcAgeLR(-SB,k, A,k, 0,-1, .TRUE., ALR(1))
  if (ALR(1)==impossible .or. ALR(1)<5.0*TF) then
    LL(1) = impossible
    return
  endif
else if (focal==2 .or. focal==3) then
  if (all(GpID(:,SB,k)==0) .and. ns(SB,k)>0) then
    call calcALR_addsib(A,SB,k,3,ALR(3))
  else
    call CalcAgeLR(A,Sex(A), -SB,k, 0,1, .TRUE., ALR(3))
  endif
  ALR(2) = ALR(3)   ! updated below if a FS found
  if (ALR(3)==impossible .or. ALR(3)<5.0*TF) then
    LL(2:3) = impossible
    return
  endif
  do Bi=1, ns(SB,k)
    if (ns(SB,k)==0)  exit
    if ((Parent(A,3-k)==Parent(SibID(Bi,SB,k),3-k) .and. Parent(A,3-k)/=0) .or. &
      Complx==0) then
      call PairQFS(A, SibID(Bi,SB,k), LRHS)
    else
      call PairQHS(A, SibID(Bi,SB,k), LRHS)
    endif
    if (LRHS < 4*TF) then
      LL(2:3) = impossible
      return
    endif
  enddo
endif

call ChkAncest(A,k,-SB,k, AncOK)
if (.not. AncOK) then
  LL(1) = impossible
  LL(4) = impossible
  if (focal==1 .or. focal==4)  return
endif
call ChkAncest(-SB,k, A,k, AncOK)
if (.not. AncOK) then
  LL(2:3) = impossible
  if (focal==2 .or. focal==3)  return
endif

! mt haplotype
if (DoMtDif) then
  if (k==1 .and. mtDif(SibID(1,SB,k), A)) then
    LL(1:3) = impossible
    if (Sex(A)==1)  LL(4) = impossible
    if (LL(focal) == impossible)  return 
  endif
endif
 
call CalcU(A,k, -SB, k, LLg(7))   ! unrelated
LL(7) = LLg(7)
 
fsi=0
if (focal < 4 .and. ns(SB,k) > 0) then                  
  if (focal==1)  call AddParent(A, SB, k, LLg(1))
  if (focal==2)  call AddFS(A, SB, k,0,k, LLg(2), fsi, dx)
  if (focal==3 .and. Complx/=0)  call AddSib(A, SB, k, LLg(3))
  do x=1,3
    if (focal==x) then
      if ((LLg(focal) > 0D0 .or. LLg(focal) - LL(7) < TA)) then
        LL(x) = addALR(LLg(x), ALR(x))
        if (focal ==3 .and. (Parent(A,3-k)==0 .or. Complx==0 .or. &
         ANY(Parent(SibID(1:ns(SB,k),SB,k), 3-k) == Parent(A,3-k)))) then
          call AddFS(A, SB, k,0,k, LLg(2), fsi, dx)
          if ((ALL(LLg(2:3) >0D0) .or. MaxLL(LLg(2:3)) - LL(7) < TA))  return
        else
          return
        endif
      endif
    endif
  enddo
endif 

fclsib = focal==2 .or. focal==3 .or. focal==6
call getFSpar(SB, k, .TRUE., FSpar) 

 !=======
if (LL(1)/=impossible) then   
  if (ALR(1)==missing)  call CalcAgeLR(-SB,k, A,k, 0,-1, .TRUE., ALR(1))
  if (LLg(1)==missing .and. ALR(1)/=impossible)  call AddParent(A, SB, k, LLg(1))  ! A parent of SB                                                                                                   
endif
if (LLg(2)==missing .and. ns(SB,k)>0)  call AddFS(A, SB, k,0,k, LLg(2), fsi, dx)
if (fsi/=0)  call CalcAgeLR(A,Sex(A), fsi,k, 0,2, .TRUE., ALR(2))
! TODO? call calcALR_addsib(A,SB,k,2,ALR(2))
if (LLg(3)==missing .and. (Complx>0 .or. ns(SB,k)==0))  call AddSib(A, SB, k, LLg(3))
!if (all(GpID(:,SB,k)==0) .and. ns(SB,k)>0) then
!  call calcALR_addsib(A,SB,k,3,ALR(3))
!else
  call CalcAgeLR(A,Sex(A), -SB,k, 0,1, .TRUE., ALR(3))  ! SB parent of A 
!endif
if (ns(SB,k)==0) then
  LLg(2) = LLg(3)
  ALR(2) = ALR(3)
  if (Complx==0) then
    LLg(3) = missing
  endif
endif

if (nYears>2 .and. LL(4)/=impossible) then
  call CalcAgeLR(-SB,k, A,Sex(A), 0,1, .TRUE., ALR(4))  ! A parent of SB
  ! TODO? analogous to calcALR_addsib()?
  if (ALR(4)/=impossible .and. (ALR(4) > 5*TF .or. focal==4)) then
      call AddGP(A, SB, k, LLg(4))
  endif
endif

do x=1,4
  LL(x) = addALR(LLg(x), ALR(x))
enddo

! monogamous
if (Complx==0 .and. Mate(A)/=0 .and. any(GpID(:,SB,k)==Mate(A))) then
  return
endif

!~~~~~~~~~~~~
if (Hermaphrodites/=0 .and. focal/=7) then 
!  if (focal/=6) then
!    call AddSelfed(A, SB, k, LLS)     ! A selfed   
!    LLg(3) = MaxLL((/LLg(3), LLS(1)/))
!  endif
    
  if (ns(SB,k)==1) then
    do Bi=1, ns(SB,k)
      LLPH = missing
      ALRH = missing
      call CalcAgeLR(A, Sex(A), SibID(Bi,SB,k),k, 0,1, .TRUE., ALRH(1))
      call CalcAgeLR(A, Sex(A), SibID(Bi,SB,k),k, 3,4, .TRUE., ALRH(2))
      if (ALRH(1)/=impossible) then
        call PairPO(A, SibID(Bi,SB,k), k, 0, LLPH(1))
      endif
      if (ALRH(2)/=impossible) then
        call PairGP(A, SibID(Bi,SB,k), k, 0, LLPH(2))
      endif
      do x=1,2
        if (LLPH(x) < 0)  LLPH(x) = LLPH(x) - Lind(SibID(Bi,SB,k)) + CLL(SB,k)
      enddo
      if (focal==1) then
        LLg(6) = MaxLL((/LLg(6), LLPH(1)/))  
        LL(6) = MaxLL((/LL(6), addALR(LLPH(1), ALRH(1))/))      
      else
        if ((LLPH(1) > LLg(1) .and. LLPH(1)<0) .or. LLg(1)>0) then
          LLg(1) = LLPH(1)
          LL(1) = addALR(LLPH(1), ALRH(1))
        endif
      endif
      if ((LLPH(2) > LLg(4) .and. LLPH(2)<0) .or. LLg(4)>0) then
        LLg(4) = LLPH(2)
        LL(4) = addALR(LLPH(2), ALRH(2))
      endif
    enddo
  endif
endif

!~~~~~~~~~~~~
LLAU = missing
ALRAU = missing
! FA 1: A FS of SB
! FA 2: SB GP of A, SB monogamous, SB's partner (thus) also GP of A
call CalcAgeLR(-SB,k, A,Sex(A), 0,2, .TRUE., ALRAU(1,3))
if (ALRAU(1,3)/=impossible .and. (Complx==0 .or. &
  (.not. (focal==4 .and. ALL(Parent(A,:)/=0)) &
   .and. .not. (focal==7 .and. GpID(3-k,SB,k)/=0))) .and. ns(SB,k)>0) then   
    call pairUA(-SB, A, k, 3, LLAU(1,3))
endif

if (Parent(A,k)/=0) then
  call CalcAgeLR(Parent(A,k),k, -SB,k, 0,1, .TRUE., ALRAU(2,3))
else
  call CalcAgeLR(A,Sex(A), -SB,k, k,4, .TRUE., ALRAU(2,3)) 
endif

if (ALRAU(2,3)/=impossible .and. .not. (focal==7 .and. fsi/=0 .and. Parent(A,3-k)==0)) then  
  ! not considered during CalcParLLR: true other-parent is GP when m=1 (single parent LLR)
  if (ns(SB,k)==1) then
    call pairUA(A, SibID(1,SB,k), k, 3, LLAU(2,3))
  else if (Parent(A,k) < 0 .and. FSpar/=0 .and. Parent(A,k)/=FSpar .and. &
    all(parent(SibID(1:ns(SB,k),SB,k), 3-k)==FSpar) .and. ns(SB,k)>0) then
    call pairUA(Parent(A,k), SibID(1,SB,k), k, 3, LLAU(2,3))
    call CalcU(Parent(A,k), k, SibID(1,SB,k), 3, LLUi)
    LLAU(2,3) = LLAU(2,3) - LLUi + LLg(7)
  else
    call addFA(A, SB, k, LLAU(2,3))
  endif
endif

LLg(5) = MaxLL(LLAU(:,3))
LL(5) = MaxLL((/ addALR(LLAU(1,3),ALRAU(1,3)), addALR(LLAU(2,3),ALRAU(2,3)) /))

LLC = missing 
if (complx==2 .and. ns(SB,k)>0 .and. (fclsib .or. focal==7) .and. LL(2)<0D0 .and. &
 Parent(A,3-k)==FSpar .and. (MaxLL(LLAU(:,3)) - MaxLL(LL(2:3)) > -TA)) then
  call FSHC(A, -SB, k, LLC)  ! Full sib & half-cousin
  if (LLC >LLg(2) .and. LLC<0) then
    LLg(2) = LLC
    LL(2) = addALR(LLg(2), ALR(2))  ! no cousin ageprior yet
  endif
endif

!~~~~~~~~~~~~
! LLg(6) HA (other 3rd degree rel: LLz further down)
if (Complx>0 .and. ns(SB,k)>0) then
  do x=1,2
    ! HA 1: A HS of SB:
    call CalcAgeLR(-SB,k, A,Sex(A), x,3, .TRUE., ALRAU(1,x))
    if (ALRAU(1,x)/=impossible .and. .not. (focal==7 .and. x==3-k .and. Parent(A,3-k)==0) .and. &
         .not. (focal==4 .and. Parent(A,x)>0 .and. GpID(x,SB,k)==0)) then  ! else conflict with CalcCandParLL
      call pairUA(-SB, A, k, x, LLAU(1,x))
    endif   
    
    ! HA 2: SB GP of A 
    if (Parent(A,x)/=0) then
      call CalcAgeLR(Parent(A,x),x, -SB,k, 0,1, .TRUE., ALRAU(2,x))
    else
      call CalcAgeLR(A,Sex(A), -SB,k, x,4, .TRUE., ALRAU(2,x))
    endif
    if (ALRAU(2,x)/=impossible  .and. .not. (focal==7 .and. x==3-k)) then    
      call pairUA(A, -SB, x, k, LLAU(2,x)) 
    endif
  enddo
  LLg(6) = MaxLL(RESHAPE(LLAU(:,1:2), (/2*2/) ))
  do x=1,2
    do y=1,2
      LLAU(y,x) = addALR(LLAU(y,x), ALRAU(y,x))
    enddo
  enddo
  LL(6) = MaxLL(RESHAPE(LLAU(:,1:2), (/2*2/) ))
endif   

! curPar = getPar(Parent(A,k),k)
! if (focal==6 .and. Parent(A,k)<0 .and. curPar(k)==-SB .and. Parent(A,3-k)==0) then
  ! ! check if A result of HS mating, but not with SB as GP
  ! call MkSibInbr(A,k,k, LLXi)
  ! ! TODO: adjust LLAU(2,2)
! endif

LLz = missing
ALRz = missing
!LLPHS = missing
!ALRPHS = missing             
if ((LL(focal)<0D0 .and. LL(focal)>=LL(7)) .or. focal==4 .or. LL(6)>0D0 .or. LL(6)<LL(7)) then
    if (any(GpID(:,SB,k)/=0)) then
      do x=1,2
        if (GpID(x,SB,k)==0) then
          call CalcAgeLR(-SB,k, A,Sex(A), x,4, .TRUE., ALRz(1))
        endif
      enddo
    else
      call CalcAgeLR(-SB,k, A,Sex(A), 3,4, .TRUE., ALRz(1))
    endif
    if (ALRz(1)/=impossible .and. ALRz(1)>3*TF) then  
      call AddGGP(A, SB, k, LLz(1))  
    endif
  if (nS(SB,k)>0) then
    do x=1,2
      if (focal==6 .and. parent(A,x)/=0)  cycle   
      call CalcAgeLR(A,k, -SB,k, x,5, .TRUE., ALRz(x+1))     
      if (ALRz(x+1)==impossible .or. ALRz(x+1)<5*TF) then
        LLz(x+1) = impossible
      else
        call ParentHFS(A, 0,x, SB, k,3, LLz(x+1))
      endif
      ! do y=1,2
        ! call CalcAgeLR(A,y, -SB,k, x,6, .TRUE., ALRPHS(y,x))  
        ! if (ALRPHS(y,x)==impossible .or. ALRPHS(y,x)<5*TF) then
          ! LLPHS(y,x) = impossible
        ! else
          ! call ParentHFS(A, 0,x, SB, k,y, LLPHS(y,x))
        ! endif
      ! enddo      
    enddo    
  endif
  if (Complx==2 .or. Complx==0) then
    do x=1,2   ! as checkmerge: full great-uncle  (2x 1/4)
      call CalcAgeLR(-SB,k, A,Sex(A), x,5, .TRUE., ALRz(3+x))
      if (ALRz(3+x) == impossible .or. ALRz(3+x)<5*TF) then
        LLz(3+x) = impossible
      else
        if (GpID(x,SB,k) <0 .and. .not. any(parent(SibID(1:ns(SB,k),SB,k),x) == GpID(x,SB,k))) then
          call PairUA(GpID(x,SB,k), A, x, 3, LLz(3+x))
          if (LLz(3+x) < 0) then
            LLz(3+x) = LLz(3+x) - CLL(-GpID(x,SB,k), x) + CLL(SB,k)  
          endif
        else if (GpID(x,SB,k)==0) then   ! else cond. indep.
          call addGAU(A, SB, k, x, LLz(3+x))    
        endif
      endif
    enddo
  endif
  if (ns(SB,k)>0) then
    sib1 = SibID(1,SB,k)
    call PairCC(A, sib1, k, LLz(6))  ! full cousins
    if (LLz(6) < 0D0) then
      call CalcU(A, k, sib1, k, LLUi)
      LLz(6) = LLz(6) - LLUi + LL(7)
    endif
    if (FSpar<0 .and. parent(A,k)==0 .and. all(GpID(:,SB,k)==0)) then
      if (ns(SB,k) == ns(-FSpar,3-k)) then
        do i=1,ns(SB,k)
          if (nFS(SibID(i,SB,k))==0)  cycle
          call pairDHC(A,k, SibID(i,SB,k), .TRUE., LLz(7))  ! double half cousins
        enddo
      endif
    endif
  endif
  ALRz(6:7) = 0D0  ! no ALR for cousins yet
  LLg(6) = MaxLL((/LLg(6), LLz/))
  do x=1,6
    LLz(x) = addALR(LLz(x), ALRz(x))
  enddo
  LL(6) = MaxLL((/LL(6), LLz/))
endif

LLM = missing    
LLp = missing
LLFH = missing   
LLPX = missing              
MaybeOpp = .FALSE.    
if (complx>0 .and. fclsib .and. hermaphrodites/=2 .and. &
 abs(MaxLL(LL(2:3)) - MaxLL(LL)) < 0.01 .and. Parent(A,3-k)==0 .and. ns(SB,k)>0) then 
  if (abs(MaxLL(LL)-LL(2)) < 0.01 .and. fsi/=0) then
    !fsi = ID of putative full sib of A within SB, returned by AddFS()                                                                                                                      
    call PairFullSib(A, fsi, LLM(1))
    call PairHalfSib(A, fsi, 3-k, LLM(2))     
    call CalcU(A, k, fsi, k, LLM(3)) 
    if ((LLM(2) - LLM(3)) - (LLg(2) - LLg(7)) > TA) then
      LL(2) = MaybeOtherParent   ! more likely to be HS via 3-k
    endif
  endif
  
  MaybeOpp = .TRUE.
  if (FSpar > 0) then
    if (FSpar/=A) then
       call CheckPair(A, FSpar, k, 1, LLpg, LLp)    !! DANGER !!!
      if (LLp(1)<0 .and. (LLp(1) - MaxLL(LLp)) > TF) then  
        LL(2:3) = MaybeOtherParent  ! FSpar plausible parent of A
      endif
    endif
  else if (FSpar==0 .and. ANY(Parent(SibID(1:ns(SB,k), SB, k), 3-k) < 0)) then                                                                             
    ! get unique opposite-sex dummy parents
    OpPar = 0
    nop = 0
    do i=1, nS(SB,k)
      if (ANY(OpPar == Parent(SibID(i,SB,k), 3-k)))  cycle
      nop = nop +1
      OpPar(nop) = Parent(SibID(i,SB,k), 3-k)
    enddo
    if (ANY(OpPar(1:nop) > 0) .or. nop > 2) then
      MaybeOpp = .FALSE.        
    endif           
    if (MaybeOpp .and. nop==2) then
      call CalcU(OpPar(1), 3-k, OpPar(2), 3-k, LLM(1))
      call MergeSibs(-OpPar(1), -OpPar(2), 3-k, LLM(2))
      if ((LLM(2) - LLM(1)) < TF*nS(SB,k))  MaybeOpp = .FALSE.
    endif
    if (nop>0 .and. MaybeOpp) then  ! .not.  ?
      do x=1, nop
        if (OpPar(x) > 0)  cycle
        call CalcU(A, 3-k, OpPar(x), 3-k, LLM(1))
        call AddSib(A, -OpPar(x), 3-k, LLM(2))
        if (LLM(2)<0D0 .and. (LLM(2) - LLM(1)) - (LLg(2) - LLg(7)) > TA*nS(SB,k)) then
          LL(2) = MaybeOtherParent   ! more likely to be added to opposing sibship only. 
        endif
        if (LLM(2)<0D0 .and. (LLM(2) - LLM(1)) - (LLg(3) - LLg(7)) > TA*nS(SB,k)) then
          LL(3) = MaybeOtherParent  
        endif
        if (LL(2)==MaybeOtherParent .and. LL(3)==MaybeOtherParent)  exit 
      enddo
    endif  
  else if (FSpar < 0) then
    call Qadd(A, -FSpar, 3-k, LLM(1))  ! 2nd degree relatives vs unrelated    
    if (LLM(1) < TF*nS(-FSpar,3-k))  MaybeOpp = .FALSE.
    call CalcAgeLR(A,Sex(A), FSpar,3-k, 0,1, .TRUE., ALRq)
    if (ALRq==impossible)  MaybeOpp = .FALSE.
  endif
  if (MaybeOpp .and. FSpar < 0) then
    LLM = missing
    call AddFS(A, -FSpar, 3-k,0,3-k, LLM(1), ix, dx)
    call AddSib(A, -FSpar, 3-k, LLM(2))
    call CalcU(A, 3-k, FSpar, 3-k, LLM(3))
    if (LLM(2) < 0D0) then
      if (complx>0) then
         if ((LLM(2) - LLM(3)) - (LLg(3) - LLg(7)) > TA*dble(MAX(nS(SB,k),nS(-FSpar,3-k)))) then
          LL(3) = MaybeOtherParent  
        endif
      endif
      if (LLM(1) < 0 .and.(LLM(1) - LLM(2)) > 2*TA) then  
        if (Complx==2) then  ! HS + parents FS/PO?
          ! TODO: not if ns(-FSpar) == ns(SB) & no GPs                                            
          curPar = Parent(A,:)                    
          call setParTmp(A, Sex(A), -SB, k)
          call PairUA(A, FSpar, 3-k, 3-k, LLPX(1,1))  ! HS + HA
          call ParentHFS(A, 0,3-k,-FSpar, 3-k,3, LLPX(1,2))  ! HS + FC
          call setParTmp(A, Sex(A), curPar(k), k)
          call setParTmp(A, Sex(A), FSpar, 3-k)     ! check done by PairFullSib
          call PairUA(A, -SB, k, k, LLPX(2,1))  ! HA + HS
          call ParentHFS(A, 0,k,SB, k,3, LLPX(2,2))  ! FC + HS
          call setParTmp(A, Sex(A), curPar(3-k), 3-k)
          if ((LLg(2) - LLg(7)) - (MaxLL(LLPX(1,:)) - LLM(3)) < TA) then
            LL(2) = MaybeOtherParent
          endif
          if ((MaxLL(LLPX(1,:)) - LLM(3)) > (LLg(3) - LLg(7)) .and. &
           (MaxLL(LLPX(1,:)) - MaxLL(LLPX(2,:))) > TA .and. MaxLL(LLPX(1,:))<0D0) then
            LLg(3) = MaxLL(LLPX(1,:)) - LLM(3) + LLg(7)  
            LL(3) = addALR(LLg(3), ALR(3))
          else if (((MaxLL(LLPX(2,:)) - LLM(3)) - (LLg(3) - LLg(7))) > TA .and. &
            MaxLL(LLPX(2,:))<0D0) then  ! MAX(nS(SB,k),nS(-FSpar,3-k))
            LL(3) = MaybeOtherParent
          endif
        else
          LL(2) = LL(2)
        endif
      else if (LLM(3)<0 .and. (LLM(2) -LLM(3)) >2*TA .and. complx>0) then
        LL(2:3) = MaybeOtherParent  ! as likely to be added to opp. parent
      endif
    endif
  endif
  if (FSpar <= 0 .and. LL(2)<0 .and. Complx==2 .and. ns(SB,k)>0) then  
    sib1 = SibID(1,SB,k)
    call calcU(A,k,sib1, k, LLFH(1))
    call pairFAHA(sib1, A, .TRUE., LLFH(2))
    call pairFAHA(A, sib1, .TRUE., LLFH(3)) 
    WHERE(LLFH(2:3)<0)  LLFH(2:3) = LLFH(2:3) - LLFH(1) + LLg(7) 
    if (ANY(LLFH(2:3)<0) .and. MaxLL(LLFH(2:3)) > LLg(5)) then
      LLg(5) = MaxLL(LLFH(2:3))
      LL(5) = MaxLL((/ addALR(LLFH(2),ALRAU(1,3)), addALR(LLFH(3),ALRAU(2,3)) /))          
    endif
  endif
endif

LLHH = missing
if ((MaxLL(LL)==LL(2) .or. MaxLL(LL)==LL(3)) .and. fclsib .and. complx==2 .and. &
  Parent(A,3-k)==0 .and. fsi/=0 .and. hermaphrodites/=2) then
  if (Parent(fsi,3-k)/=0 .or. MaxLL(LL)==LL(2)) then
    call CalcAgeLR(A,Sex(A), fsi,3-k, 3, 6, .TRUE., ALRtmp)
    if (ALRtmp /= impossible) then
      do x=1,3
        call PairHSHA(A, fsi, k, x, LLHH(x,1), .TRUE.)
      enddo    
    endif
  endif 
  if (Parent(fsi,3-k)/=0) then   ! else symmetrical 
    call CalcAgeLR(A,Sex(A), fsi,k, 3, 6, .TRUE., ALRtmp)  
    if (ALRtmp /= impossible) then
      do x=1,3
        call PairHSHA(A, fsi, 3-k, x, LLHH(x,2), .TRUE.)
      enddo
    endif
  endif 
  call CalcU(A, k, fsi, k, LLHH(4,1))
  do m=1,2
    if (MaxLL(LLHH(1:3,m)) <0D0) then
      do y=2,3 
        if (y==3 .and. Parent(fsi,3-k)==0)  cycle  ! else too many false negs
        if ((LLg(y) - LLg(7)) - (MaxLL(LLHH(1:3,m)) - LLHH(4,1)) < TA) then
          LLg(y) = MaybeOtherParent
          LL(y) = MaybeOtherParent
        endif
      enddo
    endif
  enddo
endif

LHH = missing
LHH2 = missing
if (complx==2 .and. nYears>1 .and. ns(SB,k)>0 .and. (fclsib .or. (focal==7 .and. GpID(3-k,SB,k)==0)) &
 .and. MaxLL(LL(2:3))<0D0 .and. MaxLL(LL(2:3))>=LL(7)) then
  call AddSibInbr(A, SB, k, LHH)  
  ! 1: FSpar(Parent(A,3-k),k)=SB, 2: Parent(A,3-k)=GpID(3-k,SB,k), 3: as 1, A FS of B's (PA == DB)
  if (ns(SB,k)==1 .and. all(GpID(:,SB,k)==0)) then
    call pairHSHAI(SibID(1,SB,k),A,k, LHH2)  ! B1 inbred  (needed for symmetry)
    if (LHH2 < 0D0 .and. LHH2 > LHH(2))  LHH(2) = LHH2
  endif
  if (MaxLL(LHH(1:2)) - LLg(3) > 2*TA .and. MaxLL(LHH(1:2))<0D0) then
    LLg(3) = MaxLL(LHH(1:2))
    LL(3) = addALR(LLg(3), ALR(3))
  endif
  if (LHH(3) - LLg(2) > 2*TA .and. LHH(3)<0D0) then  ! MAX(LLg(3), LLg(2))
    LLg(2) = LHH(3)
    LL(2) = addALR(LLg(2), ALR(2))
  else if (LL(3)==impossible .and. LLg(2)<0 .and. MaxLL(LHH)<0D0 .and. &
    MaxLL(LHH) > LLg(2)) then  ! add inbred FS
    LLg(2) = MaxLL(LHH)
    LL(2) = addALR(MaxLL(LHH), ALR(2))
  endif
endif

LLy = missing
if (Complx==2 .and. fclsib .and. ns(SB,k)>0 .and. MaxLL(LL(2:3))<0D0 .and. MaxLL(LLg(2:3))>MaxLL(LLg(5:7))) then
  do x=1,2
    if (ALRAU(1,x)/= impossible) then
      call HSmating(-SB, k, A, k, x, LLy(1,x))
    endif
    if (ALRAU(2,x)/= impossible) then
      call HSmating(A, k, -SB, k, x, LLy(2,x))
    endif
  enddo
  if (any(LLy < 0D0)) then
    LLg(6) = MaxLL((/LLg(6), LLy(1,:), LLy(2,:)/))
    do x=1,2
      do m=1,2
        LLy(m,x) = addALR(LLy(m,x), ALRAU(m,x))
      enddo
    enddo
    LL(6) = MaxLL((/LL(6), LLy(1,:), LLy(2,:)/))
  endif
endif 

LLpo = missing
ALRpo = missing
LLgp = missing
ALRgp = missing
if (focal/=7 .and. (FSpar<0 .or. FSpar==A .or. &
  COUNT(nFS(SibID(1:ns(SB,k),SB,k))>0) <= 5)) then   ! no. full-sib groups
  ! one of Bi parent of A?  (NOT when called by CalcCandParLL() !)
  call CalcAgeLR(A,Sex(A), -SB,k, 3-k,4, .TRUE., ALRq)
  if (ANY(Parent(A,:)<=0) .and. ALRq /=impossible .and. ALRq > 3*TF .and. focal/=6) then 
    ParTmp = Parent(A,:)
    do m=1,2
      if (Parent(A,m)>0) cycle
      do i=1, ns(SB,k)
        Bi = SibID(i,SB,k)
        if (AgeDiff(A, Bi) < 0) cycle
        if (Sex(Bi)/=m .and. Sex(Bi)<3)  cycle
        call ChkValidPar(A, Sex(A), Bi,m, ParOK)
        if (.not. ParOK)  cycle
        call setParTmp(A, Sex(A), Bi, m)
        call CalcU(A, Sex(A), -SB, k, LLPO(i,1))
        call setParTmp(A, Sex(A), ParTmp(m), m)
        call CalcCLL(SB, k)
        call CalcLind(A)
        call CalcAgeLR(A,Sex(A), Bi,m,0,1,.TRUE., ALRpo(i,1))
      enddo
    enddo
  endif

  ! A parent 3-k of one of Bi?
  if (focal/=1 .and. (Sex(A)==3-k .or. Sex(A)>2)) then  
    do i=1, ns(SB,k)
      Bi = SibID(i,SB,k)
      if (AgeDiff(Bi, A) < 0 .or. Parent(Bi,3-k)>0) cycle
      call ChkValidPar(Bi,Sex(Bi), A, 3-k, ParOK)
      if (.not. ParOK)  cycle
      ParTmp = Parent(Bi,:)
      call setParTmp(Bi, Sex(Bi), A, 3-k)
      call CalcU(A, Sex(A), -SB, k, LLPO(i,2))
      call setParTmp(Bi, Sex(Bi), ParTmp(3-k), 3-k) 
      call CalcCLL(SB, k)
      call CalcAgeLR( Bi,Sex(Bi), A,3-k,0,1,.TRUE., ALRpo(i,2))
    enddo
  endif

  if (any(LLpo < 0D0)) then
    LLg(6) = MaxLL((/LLg(6), LLpo(:,1), LLpo(:,2)/))
    do x=1,2
      do i=1, ns(SB,k)
        LLpo(i,x) = addALR(LLpo(i,x), ALRpo(i,x))
      enddo
    enddo
    LL(6) = MaxLL((/LL(6), LLpo(:,1), LLpo(:,2)/))
  endif

  ! one of Bi grandparent of A?
  if (ANY(Parent(A,:)<=0) .and. focal/=6) then
    do i=1, ns(SB,k)
      Bi = SibID(i,SB,k)
      if (AgeDiff(A, Bi) < 1) cycle
      do x=1,2
        if (Parent(A,x) > 0)  cycle
        call CalcAgeLR(A,Sex(A), Bi,Sex(Bi),x,4,.TRUE., ALRgp(i,x))
        if (ALRgp(i,x) == impossible .or. ALRgp(i,x) < 3*TF)  cycle
        call PairGP(A, Bi, x, 4, LLgp(i,x))
        if (LLgp(i,x) < 0D0) then
          call CalcU(A,k, Bi,k, LLUi)
          LLgp(i,x) = LLgp(i,x) - LLUi + LLg(7)
        endif
      enddo
    enddo
  endif
  
  ! A GP of one of Bi via 3-k?
  do i=1, ns(SB,k)
    Bi = SibID(i,SB,k)
    if (nFS(Bi)==0 .or. Parent(Bi,3-k)>0 .or. nFS(Bi)==ns(SB,k)) cycle  ! last case: safety net elsewhere
    if (AgeDiff(Bi, A) < 1) cycle
    call CalcAgeLR(Bi,Sex(Bi), A,Sex(A), 3-k,4,.TRUE., ALRgp(i,3))
    if (ALRgp(i,3) == impossible .or. ALRgp(i,3) < 3*TF)  cycle
    call PairGP(Bi, A, 3-k, 3, LLgp(i,3))
    if (LLgp(i,3) < 0D0) then
      call CalcU(A,k, Bi,k, LLUi)
      LLgp(i,3) = LLgp(i,3) - LLUi + LLg(7)
    endif
  enddo
  
  if (any(LLgp < 0D0)) then
    LLg(6) = MaxLL((/LLg(6), LLgp(:,1), LLgp(:,2), LLgp(:,3)/))
    do x=1,3
      do i=1, ns(SB,k)
        LLgp(i,x) = addALR(LLgp(i,x), ALRgp(i,x))
      enddo
    enddo
    LL(6) = MaxLL((/LL(6), LLgp(:,1), LLgp(:,2), LLgp(:,3)/))
  endif
endif

! one of Bi FA of A?
! LLfa = missing
! if (focal/=7 .and. ns(SB,k)<=4 .and. Parent(A,k)==0 .and. &
 ! ALRAU(2,3)/=impossible .and. LLAU(2,3)>0) then
  ! do i=1, ns(SB,k)
    ! Bi = SibID(i,SB,k)
    ! if (nFS(Bi)==0)  cycle
    ! call pairUA(A, Bi, k, 3, LLfa(i))
    ! if (LLfa(i) < 0) then
      ! call CalcU(A,k, Bi,k, LLUi)
      ! LLfa(i) = LLfa(i) - LLUi + LLg(7)                                 
    ! endif
  ! enddo
  ! if (any(LLfa < 0D0)) then
    ! LLg(6) = MaxLL((/LLg(6), LLfa/))
    ! do i=1, ns(SB,k)
      ! LLfa(i) = addALR(LLfa(i), ALRAU(2,3))
    ! enddo
    ! LL(6) = MaxLL((/LL(6), LLfa/))
  ! endif                                  
! endif

LLfs = Missing
if (.not. fclsib .and. Parent(A,k)==0 .and. Parent(A,3-k)<0 .and. FSpar<0 .and. &
  LL(2)==impossible .and. LL(3)/=impossible .and. ns(SB,k)>0) then
  ! check if FS anyway, by merging parents 3-k
  ! not merged because unclear if pat or mat merge needed.
  ParTmp = Parent(A,:)
  call setParTmp(A,0,0,3-k)
  call AddFS(A, SB, k,0,k, LLfs(1,1), fsi, dx)
  if (fsi/=0)  call CalcAgeLR(A,Sex(A), fsi,k, 0,2, .TRUE., ALR(2))
  if (LLfs(1,1) < 0 .and. ALR(2)/=impossible) then
    call AddSib(A, SB, k, LLfs(2,1))    ! w/o parent(A,3-k)
    call CalcU(A, k, -SB, k, LLfs(3,1)) 
    if (LLfs(1,1) - MAXVAL(LLfs(2:3,1)) > TA) then
      call MergeSibs(-FSpar, -ParTmp(3-k), 3-k, LLfs(1,2))
      call parenthfs(0,-FSpar,3-k, -ParTmp(3-k),3-k, 3, LLfs(2,2))
      call CalcU(FSpar, 3-k, ParTmp(3-k), 3-k, LLfs(3,2))
     if (LLfs(1,2) - MAXVAL(LLfs(2:3,2)) > TA) then
       LLg(2) = LLfs(1,1) - LLfs(3,1) + LLg(7)
       LL(2) = addALR(LLg(2), ALR(2))
     endif
    endif                            
  endif
  call setParTmp(A,0,ParTmp(3-k),3-k)
endif

do x=1,4
  if (LL(x) > 0) then
    LLg(x) = LL(x)
  endif
enddo

! if ( A==2969 .and. (any(SibID(:,SB,k)==2968) .or. any(SibID(:,SB,k)==3255)) ) then
 ! open (unit=42,file="log.txt",status="unknown", position="append")
  ! write (42, *) ""
    ! write (42, '("add?", 3i6, " + ", 2i6, " GPs: ", 2i6)') A, Parent(A,:), SB, k, GpID(:,SB,k)
    ! write (42, '("LL  ", 7f9.2, " ", i4)') LL, DumClone(SB,k)
    ! write (42, '("LLG ", 7f9.2, "  ", i3)') LLg, focal
    ! write (42, '("ALR ", 7f9.2, "  ", i3)') ALR
    ! write (42, '("LLAU ", 6f8.1)') LLAU(1,:), LLAU(2,:)  !, maybeFA  , ", maybeFA: ", l3
    ! write (42, '("ALRau ", 6f8.1)') ALRAU(1,:), ALRAU(2,:)
    ! write (42, '("LLZ ", 7f8.1)') LLz
    ! write (42, '("ALRz: ", 7f8.1, ", ALRq: ", f8.1)')  ALRz, ALRq
! !    write (42, '("LLPHS ", 2f8.1, "; ", 2f8.1)') LLPHS(:,1), LLPHS(:,2)                                                                   
! !    write (42, '("LLgpX ", 6f8.1)') LLgpX
    ! write (42, '("LLU ", f9.2, "; ", 3f9.2)') LLg(7), Lind(A) + CLL(SB,k), Lind(A), CLL(SB,k) 
    ! write (42, '("LLM ", 2L2, i4, "; ", 3f8.1, ", LLPX: ", 4f8.1)') &
     ! MaybeOpp, FSpar, fsi, LLM, LLPX(1,:), LLPX(2,:)
     ! if (ANY(LLP<missing))  write (42, '("LLP ", 7f8.1)') LLp
    ! write (42, '("LLHH: ", 4f8.1, " ; ", 4f8.1)')  LLHH(:,1), LLHH(:,2)
    ! write (42, '("LLC: ", f8.1, ", LHH: ", 3f8.1, ", ", f8.1)') LLC, LHH, LHH2
     ! if (FSpar<0) write(42,*)  "ns FSpar: ", ns(-FSpar,3-k)
     ! write (42, '("LLPO-1: ", 50f8.1)') LLPO(:, 1) 
     ! write (42, '("LLPO-2: ", 50f8.1)') LLPO(:, 2)
     ! write (42, '("LLGP-1: ", 50f8.1)') LLGP(:, 1) 
     ! write (42, '("LLGP-2: ", 50f8.1)') LLGP(:, 2)
     ! write (42, '("LLGP-3: ", 50f8.1)') LLGP(:, 3)
     ! write (42, '("LLFS  : ", 50f8.1)') LLFS
     ! write (42, '("LLy: ", 4f8.1)')  LLy(1,:), LLy(2,:)
! !     write (42, '("LLXi  : ", f8.1)') LLXi
    ! do x=1, nS(SB, k)
      ! write(42,'(i3, " ",a8, 2i6, 2X, f8.1, " fs", i4, 10i6)') &
       ! SB, ID(SibID(x,SB,k)), Parent(SibID(x,SB,k),:), &
          ! Lind(SibID(x,SB,k)), nFS(SibID(x,SB,k)), FSID(1:nFS(SibID(x,SB,k)), SibID(x,SB,k))
    ! enddo
    ! write (42, *) ""
    
    ! ! write(42,*) 'ALR:'
    ! ! do x=1, nS(SB, k)  
      ! ! call CalcAgeLR(A,sex(A), SibID(x,SB,k),3, k,3, .TRUE., ALR(3)) 
      ! ! call CalcAgeLR(A,sex(A), SibID(x,SB,k),k, k,6, .TRUE., ALR(6))       
      ! ! write(42, '(2i6, 2f9.2)')  SibID(x,SB,k), BY(SibID(x,SB,k)), ALR(3), ALR(6)
    ! ! enddo
    
  ! close(42)
! endif

end subroutine CheckAdd 

! #####################################################################

subroutine Qmerge(SA, SB, k, LR)
use Global
implicit none

integer, intent(IN) :: SA, SB, k
double precision, intent(OUT) :: LR
integer :: l, x, y
double precision :: PrL(nSnp), PrX(3), PrXY(3,3)

PrL = 0D0
do l=1,nSnp
  do x=1,3
    PrX(x) = XPr(1,x,l,SA,k)*XPr(1,x,l,SB,k)* AHWE(x,l)
    do y=1,3
      PrXY(x,y) = XPr(1,x,l,SA,k)*XPr(1,y,l,SB,k)* AHWE(x,l) * AHWE(y,l)
    enddo
  enddo
  PrL(l) = LOG10(SUM(PrX)) - LOG10(SUM(PrXY))   ! merge
enddo
LR = SUM(PrL)

end subroutine Qmerge

! #####################################################################

subroutine QFSmerge(SA, SB, k, LR)
use Global
use CalcLik
implicit none

integer, intent(IN) :: SA, SB, k
double precision, intent(OUT) :: LR
integer :: l, x, z, par(2), i, j
double precision :: PrL(nSnp), PrXZ(3,3,2), PrI(3,3), PrJ(3,3)

par = 0       
call getFSpar(SA,k, .TRUE., par(1))
call getFSpar(SB,k, .TRUE., par(2))
if (any(par==0))  return

i = FSID(maxSibSize+1, SibID(1,SA,k))
j = FSID(maxSibSize+1, SibID(1,SB,k))

PrL = 0D0
do l=1,nSnp
  PrI = FSLik(l,i)
  PrJ = FSLik(l,j)
  do x=1,3
    do z=1,3
      PrXZ(x,z,:) = prI(x,z) * AHWE(x,l) * AHWE(z,l)
      PrXZ(x,z,1) = PrXZ(x,z,1) * prJ(x,z)
      PrXZ(x,z,2) = PrXZ(x,z,2) * SUM(PrJ(x,:) * AHWE(:,l)) 
    enddo
  enddo
  PrL(l) = LOG10(SUM(PrXZ(:,:,1))) - LOG10(SUM(PrXZ(:,:,2)))
enddo
LR = SUM(PrL)

end subroutine QFSmerge

! #####################################################################

subroutine CheckMerge(SA, SB, kA, kB, focal, LLg, LL, FSM) 
use Global
implicit none

integer, intent(IN) :: SA, SB, kA, kB, focal
double precision, intent(OUT) :: LLg(7), LL(7)
logical, intent(OUT) :: FSM
double precision ::  ALRtmp, LRHS, ALR(7), LLtmp(3), dx(maxSibSize), LLx(6), ALRx(6), &
  LLz(2,2), LLY(2,3), ALRy(2,3), LLHA(3), dLH(nS(SB,kB)), LLM(5), LLMo(5), LLHHA(2), &
  LLC, LLP, TAx
integer :: i, j, x, Par(2), ix, tmpGP, NSx(2,2), OpPars(maxSibSize, 2), DoQuickA, DoQuickB
logical :: ShareOpp, ShareSib, AncOK(2), ParOK

LLg = missing
LL = missing
FSM = .FALSE.  ! merge both k & 3-k
ShareOpp = .FALSE.
ShareSib = .FALSE.
if (kA /= kB) then
  LL(1) = impossible
  if (focal==1)  return
endif
if (focal==1 .and. kA==1 .and. DoMtDif) then
  if (mtDif(SibID(1,SA,kA), SibID(1,SB,kB))) then
    LL(1) = impossible
    return
  endif
endif
do i=1, nS(SA, kA)
  do x=1,2
    if (SibID(i, SA, kA)==GpID(x,SB,kB)) then
      LL(1) = impossible
      exit
    endif
  enddo
enddo
do j=1, nS(SB, kB)
  do x=1,2
    if (SibID(j, SB, kB)==GpID(x,SA,kA)) then
      LL(1) = impossible
      exit
    endif
  enddo
enddo
do i=1, nS(SA, kA)
  do j=1, nS(SB, kB)
    if (AgeDiff(SibID(i,SA,kA), SibID(j,SB,kB))==missing) cycle
    if (getAP( AgeDiff(SibID(i,SA,kA), SibID(j,SB,kB)), 3, 0, kA, Impossible) == Impossible) then
      LL(1) = impossible
      exit
    endif
    if (LL(1)==impossible) exit
    if (kA/=kB) then
      if (SibID(i, SA, kA)==SibID(j,SB,kB)) then
        ShareSib = .TRUE.
      endif
    endif
  enddo
enddo 
if (LL(1) == impossible .and. focal==1) return

if (focal==1) then
  call ChkAncest(-SA,kA, -SB,kB, AncOK(1))
  call ChkAncest(-SB,kB, -SA,kA, AncOK(2))
  if (any(.not. AncOK)) then
    LL(1) = impossible
    return
  endif
  
  do x=1,2
    if (GpID(x,SB,kB)/=0) then  !  .and. GpID(x,SB,kB)/=GpID(x,SA,kA)
      call CalcAgeLR(-SA,kA, GpID(x,SB,kB),x, 0,1, .TRUE., ALRtmp)
      if (ALRtmp == impossible) then
        LL(1) = impossible
        return
      endif
    endif
    if (GpID(x,SA,kA)/=0) then
      call CalcAgeLR(-SB,kB, GpID(x,SA,kA),x, 0,1, .TRUE., ALRtmp)
      if (ALRtmp == impossible) then
        LL(1) = impossible
        return
      endif
    endif
  enddo
  
  if (hermaphrodites==1 .and. ((DumClone(SA,kA)/=0 .or. DumClone(SB,kB)/=0) .and. &
   .not. (DumClone(SA,kA)/=0 .and. DumClone(SB,kB)/=0)))  then
    LL(1) = MaybeOtherParent      
    return
  endif
endif

OpPars = 0
if (kA==kB .and. Complx/=0) then
  do i=1, nS(SA, kA)
    OpPars(i,1) = Parent(SibID(i,SA,kA), 3-kA)
  enddo
  do j=1, nS(SB,kB)
    OpPars(j,2) = Parent(SibID(j,SB,kB), 3-kB)
    if (OpPars(j,2)/=0 .and. ANY(opPars(1:ns(SA,kA),1) == opPars(j,2))) then
      ShareOpp = .TRUE.
    endif
  enddo
endif

LRHS = missing
if (.not. ShareOpp .and. .not. ShareSib .and. Complx/=0) then
  call ChkDoQuick(SA, kA, DoQuickA)
  call ChkDoQuick(SA, kA, DoQuickB) 
  if (DoQuickA /= 2 .and. DoQuickB /= 2) then
    call Qmerge(SA, SB, kB,  LRHS)
    if (LRHS < 2.0*TF*dble(MAX(nS(SA,kA), nS(SB,kB)))) then
      LL(1) = impossible
    endif
    if (LL(1) == impossible .and. focal==1) return
  endif
endif

if (focal==1) then
  do i=1, ns(SA,kA)
    do j=1, nS(SB,kB)
      if ((OpPars(i,1)==OpPars(j,2) .and. OpPars(i,1)/=0) .or. Complx==0) then 
        call PairQFS(SibID(i,SA,kA), SibID(j,SB,kB), LRHS)
      else
        call PairQHS(SibID(i,SA,kA), SibID(j,SB,kB), LRHS)
      endif
      if (LRHS < 4*TF) then
        LL(1) = impossible
        return
      endif
    enddo
  enddo
endif                  

 call CalcU(-SA,kA, -SB,kB, LLg(7))
 LL(7) = LLg(7)
 
ALR = missing              
if (LL(1)/=impossible)  call CalcALRmerge(SA, SB, kA, ALR(1))              
if (LL(1)/=impossible .and. ALR(1)/=impossible) then
  if (Complx/=0 .or. focal==8) then
    call MergeSibs(SA, SB, kA, LLg(1))   ! SB parent of A's
  else
    Par = 0
    LLM = missing
    call getFSpar(SA, kA, .TRUE., Par(1))
    call getFSpar(SB, kB, .TRUE., Par(2))
    if ((Par(1)==Par(2) .and. Par(1)/=0) .or. all(Parent(SibID(:,SA,kA),3-kA)==0) .or. &
     all(Parent(SibID(:,SB,kB),3-kB)==0)) then
      call MergeSibs(SA, SB, kA, LLg(1))  
    else
      call FSMerge(SA,SB,kA, LLM)
      LLg(1) = LLM(4)   ! merge both k & 3-k
    endif
  endif    
  LL(1) = addALR(LLg(1), ALR(1))
else
  LL(1) = impossible
  LLg(1) = LL(1)                
endif

if (focal==1 .and. (LLg(1) > 0D0 .or. LL(1)==impossible .or. &
  (LL(1) - LL(7) < TA .and. Complx/=0))) return

call CalcAgeLR(-SB,kB, -SA,kA, 0,1, .TRUE., ALR(2))
if (ALR(2) /= impossible) then
  call addFS(0, SA, kA, SB, kB, LLg(2), ix, dx)  ! SB FS with an A
  if(complx>0)  call PairUA(-SB, -SA, kB, kA, LLg(3))  ! SB HS with an A
  do x=2,3
    LL(x) = addALR(LLg(x), ALR(2))
  enddo
else
  LL(2:3) = impossible
endif

if (focal==1 .and. (LL(1) - MaxLL(LL(2:7)) < TA) .and. ANY(LL(2:7) < 0) .and. Complx/=0) return

LLtmp = missing
tmpGP = 0        
call CalcAgeLR(-SA,kA, -SB,kB, 0,1, .TRUE., ALR(4))
if (ALR(4)/=impossible) then
  if (focal/=4 .or. complx==0) then
    call addFS(0, SB, kB, SA, kA, LLtmp(1), ix, dx)  ! SB GP of A's
  endif
  if (focal==4) then  ! allow for replacement
    tmpGP = GpID(kB,SA,kA)
    call setParTmp(-SA, kA, 0, kB)
  endif
  if (complx>0)  call PairUA(-SA, -SB, kA, kB, LLtmp(2))  ! SB GP of A's
  if (hermaphrodites/=0)  call addHAselfed(SA,kA,SB,kB, LLtmp(3))
  if (focal==4)   call setParTmp(-SA, kA, tmpGP, kB)
  LLg(4) = MaxLL(LLtmp)
  LL(4) = addALR(LLg(4), ALR(4))
else
  LL(4) = impossible
endif

if (.not. focal==4 .and. any(GpID(:,SA,kA)==0 .and. GpID(:,SB,kB)/=0)) then  
  ! else GB assigned as GP of SA as side-effect, messes up CalcCandPar
  call CalcAgeLR(-SA,kA, -SB,kB, 0,2, .TRUE., ALR(5))
  if (ALR(5) /= impossible) then 
    if(complx>0)  call ParentHFS(0, SA, kA, SB, kB,3, LLg(5))  ! SB & SA are FS
  ! TODO: PairUA for FS clusters
    LL(5) = addALR(LLg(5), ALR(5))
  else
    LL(5) = impossible
  endif
endif

LLx = missing
ALRx = 0D0
do x=1,4
  if (x==1 .or. x==2) then
    if (complx==0) cycle
    if (focal==4 .and. GpID(x,SB,kB)>0)  cycle  ! else messes up CalcCandPar
    call CalcAgeLR(-SA,kA, -SB,kB, x,3, .TRUE., ALRx(x))
    if (ALRx(x) /= impossible) then
      call ParentHFS(0, SA, kA, SB, kB, x, LLx(x))
    endif
  else if (x==3) then
    call CalcAgeLR(-SA,kA, -SB,kB, 3,4, .TRUE., ALRx(x))
    if (ALRx(x) /= impossible) then
      call dummyGP(SA, SB, kA, kB, focal, LLx(3))  ! SB GGP of A's
    endif      
  else if (x==4) then
    call CalcAgeLR(-SB,kB, -SA,kA, 3,4, .TRUE., ALRx(x))
    if (ALRx(x) /= impossible) then
      call dummyGP(SB, SA, kB, kA, 0, LLx(4))  ! SA GGP of B's
    endif 
  endif
enddo

LLz = missing
do x=1,2
  if (GpID(x, SA, kA) > 0) then   ! TODO: more general
    if (Parent(GpID(x, SA, kA), kB)==-SB) then
      LLz(x,2) = impossible
    else
      do i=1,2
        if (focal==4 .and. Parent(GpID(x, SA, kA), i)/=0) then
          LLz(x,2) = NotImplemented   ! becomes not about SA/SB
        else if (Parent(GpID(x, SA, kA), i)/=0 .and. GpID(i,SB,kB)/=0 .and. &
         Parent(GpID(x, SA, kA), i)/=GpID(i,SB,kB)) then
          LLz(x,2) = Impossible
        else 
          call ChkValidPar(GpID(x, SA, kA), x, GpID(i,SB,kB), i, ParOK)
          if (.not. ParOK)  LLz(x,2) = impossible
        endif
      enddo
    endif
    if (LLz(x,2)==missing) then
      call CalcU(-SB, kB, GpID(x,SA,kA), x, LLz(x,1))
      call PairUA(-SB, GpID(x,SA,kA), kB, 3, LLz(x,2))
      call CalcAgeLR(-SB,kB, GpID(x,SA,kA),x,0,2, .TRUE., ALRx(4+x))
    endif
    if (LLz(x,2) < 0D0 .and. ALRx(4+x)/=impossible) then
      LLx(4+x) = LL(7) + LLz(x,2) - LLz(x,1)
    endif
  endif
enddo

LLg(6) = MaxLL(LLx)  ! most likely 3rd degree relative
do x=1,6
  LLX(x) = addALR(LLX(x), ALRx(x))
enddo
LL(6) = MaxLL(LLx)

! avuncular relationships between dummies?
LLY = missing
ALRy = missing
do x=1,3
  if (x < 3)  call CalcAgeLR(-SA,kA, -SB, kB, x, 6, .TRUE., ALRy(1,x))
  if (x == 3)  call CalcAgeLR(-SA,kA, -SB, kB, 3, 5, .TRUE., ALRy(1,x))   
  if (ALRy(1,x) /= impossible) then
    call dummyHFA(SA,kA, SB,kB, x, LLy(1,x))
  endif
enddo
do x=1,3
  if (x < 3)  call CalcAgeLR(-SB, kB, -SA,kA,  x, 6, .TRUE., ALRy(2,x))
  if (x == 3)  call CalcAgeLR(-SB, kB, -SA,kA, 3, 5, .TRUE., ALRy(2,x))   
  if (ALRy(2,x) /= impossible) then
    call dummyHFA(SB,kB, SA,kA, x, LLy(2,x))
  endif
enddo

if (any(LLY < 0)) then
  LLg(6) = MaxLL((/LLg(6), LLY(1,:), LLY(2,:)/))
  do x=1,3
    do j=1,2
      LLY(j,x) = addALR(LLY(j,x), ALRy(j,x))
    enddo
  enddo
  LL(6) = MaxLL((/LL(6), LLY(1,:), LLY(2,:)/))
endif

LLHA = missing     
dLH = missing           
if (complx>0 .and. LL(4)<0D0 .and. focal/=4 .and. LLtmp(1)<0D0 .and. &
  LLtmp(1)>=LLtmp(2) .and. LLtmp(1) > MaxLL((/LL(1:3), LL(5:7)/))) then   
  do j=1, nS(SB, kB)
    if (GpID(3-kB,SA,kB)==0) then
      shareOpp = .TRUE.
      par(1) = Parent(SibID(j, SB, kB), 3-kB)
    else if (Parent(SibID(j, SB, kB), 3-kB)==GpID(3-kB,SA,kB) .or. &
      Parent(SibID(j, SB, kB), 3-kB)==0) then
      shareOpp = .TRUE.
      par(1) = GpID(3-kB,SA,kB)
    else
      shareOpp = .FALSE.
    endif
    if (shareOpp) then
      call PairUA(-SA, SibID(j, SB, kB), kA, 3, LLHA(1))
      call PairUA(-SA, SibID(j, SB, kB), kA, 3-kB, LLHA(2))
      call CalcU(-SA, kA, SibID(j, SB, kB), kB, LLHA(3))
      if (LLHA(1)<0)  dLH(j) = LLHA(1) - MaxLL(LLHA(2:3))
    endif
  enddo
  if (MAXVAL(dLH, MASK=dLH<missing) < TA) then
    LLg(4) = LLtmp(2)    ! do not use LLtmp(1) from addFS
    LL(4) = addALR(LLg(4), ALR(4))
  endif
endif

LLM = missing
LLMo = missing
Par = 0
NSx = 0
LLHHA = missing
if (kA == kB .and. (Complx==0 .or. (focal==1 .and. (ABS(MaxLL(LL) - LL(1))<TA)))) then
  call FSMerge(SA,SB,kA, LLM) ! 1:not, 2: via k, 3: via 3-k, 4:both, !!5: 3-k + par HS
  if (Complx/=2)  LLM(5) = 555D0   ! merge via 3-k + par HS
  LLM(1) = MaxLL((/LLM(1), LL(7)/))  ! do not merge  
  LLM(2) = MaxLL((/LLM(2), LLg(1), LLM(5)/))  ! merge via k  
  call getFSpar(SA, kA, .FALSE., Par(1))
  call getFSpar(SB, kB, .FALSE., Par(2))
  if (par(1)<0 .and. par(2)<0) then
    NSx(1,1) = nS(SA,kA)
    NSx(2,1) = nS(SB,kB)
    NSx(1,2) = ns(-par(1), 3-kA)
    NSx(2,2) = ns(-par(2), 3-kB)
    if (ANY(NSx(:,1) /= NSx(:,2))) then
      call FSMerge(-par(1),-par(2),3-kA, LLMo)
      if (Complx/=2)  LLMo(5) = 555D0   ! merge via 3-k + par HS
      call CalcU(Par(1), 3-kA, Par(2), 3-kA, LLMo(1))  ! more accurate
      call MergeSibs(-par(1), -par(2), 3-kA, LLMo(2)) 
      if (LLMo(2) <0) then
        LLM(1) = MaxLL((/LLM(1), LLMo(2) - LLMo(1) + LL(7) /))
      endif
    endif
  endif
  
  TAx = TA * dble(MIN(nS(SA,kA), nS(SB,kB)))
  if (Complx>0 .and. LLM(2)<0D0 .and. LLM(1) - LLM(2) > TAx) then
    LL(1) = MaybeOtherParent 
  endif

  if (Complx==2 .and. par(1)<0 .and. par(2)<0 .and. hermaphrodites/=2) then
    if (SUM(NSx(:,1)) == SUM(NSx(:,2))) then
      if (LLM(4)<0D0) then
        call FSHC(-SA,-SB,kA,LLC)
        if (LLC > LLM(4) .and. LLC<0)  LLM(4) = LLC
      endif
      call clustHSHA(SA, SB, kA, LLHHA(1))
      if (LLHHA(1) - LLM(4) > TA)  LL(1) = MaybeOtherParent
      call clustHSHA(SB, SA, kA, LLHHA(2))
      if (LLHHA(2) - LLM(4) > TA)  LL(1) = MaybeOtherParent
    endif
  endif
  
  if (hermaphrodites > 0) then   
    if ((ns(SA,kA)==1 .and. SelfedIndiv(SibID(1,SA,kA)) .and. all(GpID(:,SA,kA)==0)) .or. &
      (ns(SB,kB)==1 .and. SelfedIndiv(SibID(1,SB,kB)) .and. all(GpID(:,SB,kB)==0))) then
      LLM(4) = LLM(2)  ! selfed singletons
      LLM(2:3) = impossible
    endif
  else if (Complx==0) then
    if (Par(1)>0 .and. Par(1)==Par(2)) then   
      LLM(4) = LLM(2)
    endif
  endif
  
  TAx = TA * dble(MIN(nS(SA,kA), nS(SB,kB)))
  if (MaxLL(LLM(1:4))==LLM(4) .and. (LLM(4)-LLM(2) >TAx .and. &
   (ALL(LLMo==missing) .or. LLMo(4)-LLMo(3) >TAx)) .or. &
   ((Complx==0 .or. hermaphrodites>0) .and. LLM(4)<0 .and. LLM(4) - LLM(1) > TAx)) then
    LLg(1) = LLM(4)  ! FS merge most likely - go ahead.
    LL(1) = addALR(LLg(1), ALR(1))  
    FSM = .TRUE.
  else if (Complx>0 .and. LLM(3)<0D0 .and. LLM(3)-LLM(1) > TA ) then   ! .and. hermaphrodites/=2
    if (LLM(3)-LLM(4) > TA) then
      LL(1) = MaybeOtherParent ! likely that opp. parent need to be merged, 
    else if (Par(1) < 0 .and. Par(2)<0) then
      if (NSx(1,1)==NSx(1,2) .and. NSx(2,1)==NSx(2,2)) then   ! 2 FS groups
        LL(1) = MaybeOtherParent
      endif
    endif
  endif
endif

if (focal==1 .and. (ABS(MaxLL(LL) - LL(1))<TA)) then  ! one of Bi is SA, or vv?
  do i=1, ns(SA,kA)
    if (AgeDiff(SibID(i,SA,kA), SibID(1,SB,kB)) < 0) cycle
    LLP = missing               
    call AddParent(SibID(i,SA,kA), SB, kB, LLP)
    if (LLP < 0D0 .and. (LLP + CLL(SA,kA) - Lind(SibID(i,SA,kA))) > LL(1)) then
      LL(1) = MaybeOtherParent
      exit
    endif
  enddo
  if (LL(1)<0D0) then
    do j=1, ns(SB,kB)
      if (AgeDiff(SibID(j,SB,kB), SibID(1,SA,kA)) < 0) cycle
      call AddParent(SibID(j,SB,kB), SA, kA, LLP)
      if (LLP < 0D0 .and. (LLP + CLL(SB,kB) - Lind(SibID(j,SB,kB))) > LL(1)) then
        LL(1) = MaybeOtherParent
      endif
    enddo
  endif 
endif

do x=1,4
  if (LL(x) > 0) then
    LLg(x) = LL(x)
  endif
enddo

! if ((any(SibID(:,SA,kA)==115) .and. ANY(SibID(:,SB,kB)==127)) .or. &
  ! (any(SibID(:,SA,kA)==127) .and. ANY(SibID(:,SB,kB)==115)))  then
! !if (SA==1 .and. kA==2 .and. SB==5 .and. kB==2) then
  ! open (unit=42,file="log.txt",status="unknown", position="append")
    ! write (42, *) ""
    ! write(42,'("merge? ", 2i5," ,",2i5,": ", i3, l3, "  GA: ", 2i5," , GB: ", 2i5)') &
      ! kA, SA, kB, SB, focal, FSM, GpID(:, SA, kA), GpID(:, SB, kB)      
    ! ! write(42,'("SA 3-k:  ", 10i5)')   SibID(1:10, 1, 3-kA)
    ! ! write(42,'("SB 3-k:  ", 10i5)')   SibID(1:10, 6, 3-kA)
    ! write(42,'("LL:  ", 7f8.1)') LL
    ! write(42,'("LLg: ", 7f8.1)') LLg
    ! write (42,'("ALR ", 7f8.1)') ALR
    ! write(42,'("LLtmp ", 3f8.1)') LLtmp
    ! write (42,'("LLx ", 8f8.1)')  LLx
    ! write (42,'("ALRx ", 6f8.1)')  ALRx
    ! write(42,'("LLZ ", 4f8.1)') LLz(1,:), LLz(2,:)
    ! write(42,'("LLY ", 6f8.1)') LLy(1,:), LLy(2,:)
    ! write(42,'("LLM ", 2i5, 5f8.1, " ; ", 2f8.1, 4i3)') Par, LLM, LLM(2)-LL(7), LLM(3)-LLM(1), NSx(:,1), NSx(:,2)
    ! write(42,'("LLMo, LLHHA ", 5f8.1, "; ", 2f8.1)') LLMo, LLHHA
    ! write(42,'("LLU: ", 3f8.1)') CLL(SA,kA), CLL(SB,kB), CLL(SA,kA) + CLL(SB,kB)
    ! do i=1, nS(SA, kA)
    ! write(42,'(i3, " ", a10, 2i6, " fs", 10i5)') SA, ID(SibID(i,SA,kA)), Parent(SibID(i,SA,kA),:), &
      ! nFS(SibID(i,SA,kA)), FSID(1:nFS(SibID(i,SA,kA)), SibID(i,SA,kA))
    ! enddo
    ! write (42, *) ""
    ! do i=1, nS(SB, kB)
        ! write(42,'(i3, " ",a10, 2i6, " fs", 10i5)') SB, ID(SibID(i,SB,kB)), Parent(SibID(i,SB,kB),:), &
          ! nFS(SibID(i,SB,kB)), FSID(1:nFS(SibID(i,SB,kB)), SibID(i,SB,kB))
    ! enddo
    ! write (42, *) ""
    ! close(42)
! endif

end subroutine CheckMerge 

! #####################################################################

subroutine getFSpar(SA, kA, strict, par)  
! all individuals in SA are FS to eachother
use Global
implicit none

integer, intent(IN) :: SA,  kA
logical, intent(IN) :: strict
integer, intent(OUT) :: Par
integer :: i, j, ParV(ns(SA,kA))

Par = 0
if (ns(SA,kA)==0)  return

ParV = 0
do i=1, nS(SA,kA)
  if (Parent(SibID(i,SA,kA), 3-kA)/=0) then
    Par = Parent(SibID(i,SA,kA), 3-kA)
    if (strict) then
      do j= i+1, nS(SA, kA)
        if (Parent(SibID(j,SA,kA), 3-kA) /= Par .and. &
         Parent(SibID(j,SA,kA), 3-kA)/=0) then
          Par = 0
          return
        endif
      enddo
    else 
      ParV(i) = Par
    endif
  endif
enddo

if (.not. strict) then ! > half by same opp. parent?
  Par = 0
  do i=1, nS(SA,kA)
   if (real(COUNT(ParV == ParV(i))) > nS(SA,kA)/2.0) then
      Par = ParV(i)
      return
    else if (real(COUNT(ParV == ParV(i))) == nS(SA,kA)/2.0 .and. ParV(i)<0) then
      Par = ParV(i)
      return
    endif
  enddo
endif

end subroutine getFSpar

! #####################################################################

subroutine OppMerge(SA, k, LL)  ! could opposing parents of SA all be the same dummy parent?
use Global
implicit none

integer, intent(IN) :: SA, k
double precision, intent(OUT) :: LL ! of SA
integer :: i, l, x, y, m,u, opPar(ns(SA,k)), GPY(2)
double precision :: PrL(nSnp), PrSA(3), PrXY(3, 3), PrGY(3,2)

opPar = 0
GPY = 0
do i=1, ns(SA,k)
  opPar(i) = Parent(SibID(i, SA, k), 3-k)
enddo
if (ANY(opPar > 0)) then
  LL = NotImplemented
  return
else
  do i=1, ns(SA,k)
    if (opPar(i) < 0) then
      if (ns(-opPar(i), 3-k) > COUNT(opPar == opPar(i))) then
        LL = missing  ! TODO?  
        return
      else
        do m=1,2
          if (GpID(m,-OpPar(i),3-k) /= GPY(m)) then
            if (GPY(m) == 0) then
              GPY(m) = GpID(m,-OpPar(i),3-k)
            else
              LL = impossible
              return
            endif
          endif
        enddo
      endif
    endif
  enddo
endif

PrL = 0D0
do l=1, nSnp
  call ParProb(l, -SA, k, -1, 0, PrSA)
  do m=1,2
    call ParProb(l, GPY(m), 3-k, 0, 0, PrGY(:,m))  
  enddo
  do x=1,3
    do y=1,3
      do u=1,3
        PrXY(x,y) = PrSA(x) * SUM(AKA2P(y, u, :) * PrGY(u,1) * PrGY(:,2))
      enddo
      do i=1, ns(SA,k)
        PrXY(x,y) = PrXY(x,y) * OKA2P(Genos(l, SibID(i,SA,k)), x, y)
      enddo
    enddo
  enddo
  PrL(l) = LOG10(SUM(PrXY))
enddo
LL = SUM(PrL)

end subroutine OppMerge

! #####################################################################

subroutine FSmerge(SA,SB,k, LL)  
! calc LL if SA and SB merged via both pat & mat
use Global
implicit none

integer, intent(IN) :: SA, SB, k
double precision, intent(OUT) :: LL(5) ! 1:not, 2: via k, 3: via 3-k, 4:both, !!5: 3-k + par HS
integer :: l, x, y, i, u,v, G(2,2),z, m, Par(2), SX(2)
double precision :: ALR, PrL(nSnp,5), PrXY(3,3), PrUV(3,3), PrXV(3,3,3,3,5),&
  PrG(3,2,2), PrX(3,2), PrTmp(3), PrY(3,2)
logical :: DoParHS, MaybeOpp(2), AncOK(2)

! TODO: currently assumes no gps of sibship 3-k, no close inbreeding
LL = missing
! check if all FS
SX = (/SA, SB/)
MaybeOpp = .FALSE.
do i=1,2
  call getFSpar(SX(i), k, .TRUE., Par(i))
  if (Par(i)>0) cycle
  if (Par(i)==0 .and. ANY(Parent(SibID(1:nS(SX(i),k),SX(i),k),3-k)>0)) cycle
  MaybeOpp(i) = .TRUE.
enddo
if (ALL(MaybeOpp)) then   
  if (Par(1)==Par(2) .and. Par(1)/=0) then 
    if (ALL(Parent(SibID(1:nS(SA,k),SA,k),3-k)==Par(1)) .and. &
      ALL(Parent(SibID(1:nS(SB,k),SB,k),3-k)==Par(2))) then
      MaybeOpp = .FALSE.   ! already share same opp parent
    endif
  else if (Par(1)<0 .and. Par(2)<0) then
    ALR = missing             
    call CalcAgeLR(Par(1),3-k, Par(2),3-k, 0,-1, .TRUE., ALR) 
    if (ALR==impossible) MaybeOpp = .FALSE.
    if (DoMtDif) then
      if (k==2 .and. mtDif(SibID(1,-Par(1),3-k), SibID(1,-Par(2),3-k))) then
        MaybeOpp = .FALSE.
      endif
    endif
    if (nS(-Par(1),3-k) > ns(SA,k) .or. nS(-Par(2),3-k) > ns(SB,k)) then
      LL = NotImplemented   ! called separately from CheckMerge on 'other side'
    endif
  endif
endif
if (ANY(.not. MaybeOpp) .or. ALL(LL==NotImplemented)) return

call ChkAncest(Par(1),3-k, Par(2),3-k, AncOK(1))
call ChkAncest(Par(2),3-k, Par(1),3-k, AncOK(2))
if (any(.not. AncOK)) then
  LL = Impossible
  return
endif 

G = 0
do i=1,2
  if (GpID(i,SA,k)/=0) then
    if(GpID(i,SA,k)/=GpID(i,SB,k) .and. GpID(i,SB,k)/=0) then
      G(i,k) = 0  ! shouldn't happen
    else
      G(i,k) = GpID(i,SA,k)
    endif
  else
    G(i,k) = GpID(i,SB,k)
  endif
  if (Par(1)<0) then
    G(i,3-k) = GpID(i, -Par(1),3-k)
    if (Par(2) < 0) then
      if (GpID(i, -Par(2),3-k) /= G(i,3-k) .and. &
       GpID(i, -Par(2),3-k)/=0 .and. G(i,3-k)/=0) then
        LL = Impossible
        return
      else if (G(i,3-k)==0 .and. GpID(i, -Par(2),3-k)/=0) then
        G(i,3-k) = GpID(i, -Par(2),3-k)
      endif
    endif
  endif
enddo

if (ALL(GPID(:,SA,k)==0) .and. ALL(GPID(:,SB,k)==0) .and. Complx==2) then
  DoParHS = .TRUE.
else
  DoParHS = .FALSE.
endif

PrL = 0D0
do l=1,nSnp 
  do m=1,2
    do i=1,2
      call ParProb(l, G(i,m), i, 0, 0, PrG(:,i,m))
    enddo
    do x=1,3
      do z=1,3
        PrTmp(z) = SUM(AKA2P(x,:,z) * PrG(:,1,m) * PrG(z,2,m))
      enddo
      PrX(x,m) = SUM(PrTmp)
    enddo
  enddo
  do i=1,2
    call ParProb(l, Par(i), 3-k, -1, 0, PrY(:,i))
  enddo
  do x=1,3  ! P1
    do y=1,3  ! P2
      PrXY(x,y) = 1D0  ! XPr(2,x,l, sA,k) * AHWE(y,l)
      do i=1,nS(SA,k)
        PrXY(x,y) = PrXY(x,y) * OKA2P(Genos(l, SibID(i,SA,k)), x, y)
      enddo
    enddo
  enddo
  do u=1,3
    do v=1,3
      PrUV(u,v) = 1D0  ! XPr(2,u,l, sB,k) * AHWE(v,l)
      do i=1,nS(SB,k)
        PrUV(u,v) = PrUV(u,v) * OKA2P(Genos(l, SibID(i,SB,k)), u, v)
      enddo
    enddo
  enddo

  PrXV = 0D0
  do x=1,3 
    do y=1,3
      do u=1,3
        do v=1,3
          PrXV(x,y,u,v,1) = PrXY(x,y) * XPr(2,x,l, sA,k) * PrY(y,1) * &
            PrUV(u,v) * XPr(2,u,l, sB,k) * PrY(v,2)
          if (Complx/=0) then
            PrXV(x,y,x,v,2) = PrXY(x,y) * PrX(x,k) * PrY(y,1) * &
              PrUV(x,v) * PrY(v,2)
          endif
        enddo
        if (Complx/=0) then
          PrXV(x,y,u,y,3) = PrXY(x,y) * XPr(2,x,l, sA,k) * PrX(y,3-k) * &
              PrUV(u,y) * XPr(2,u,l, sB,k)
        endif
        if (DoParHS) then
          do z=1,3
            PrTmp(z) = AKAP(x,z,l) * AKAP(u,z,l) * AHWE(z,l)
          enddo
          PrXV(x,y,u,y,5) = PrXY(x,y) * PrX(y,3-k) * PrUV(u,y) * &
            SUM(PrTMP)
        endif
      enddo
      PrXV(x,y,x,y,4) = PrXY(x,y) * PrX(x,k) * PrX(y,3-k) * PrUV(x,y)
    enddo
  enddo
  do x=1,5
    PrL(l,x) = LOG10(SUM(PrXV(:,:,:,:,x)))
  enddo
enddo
LL = SUM(PrL,DIM=1)
if (.not. DoParHS)  LL(5) = impossible
if (Complx==0)  LL(2:3) = impossible                                    

end subroutine FSmerge

! #####################################################################

subroutine NewSibship(A, B, k)  ! make new sibship
use Global
implicit none

integer, intent(IN) :: A, B, k
integer :: s

nC(k) = nC(k) + 1
s = nC(k)
DumBY(:,s,k,1) = LOG10(1.0D0/nYears)
call SetPar(A, Sex(A), -s, k) 
if (B/=0) then
  call SetPar(B, Sex(B), -s, k) 
  if (BY(A) < 0)  call setEstBY(A, Sex(A))
  call UpdateLL(Parent(A,3-k), 3-k)
endif
call CalcCLL(s, k)

IsNewSibship(s,k) = .TRUE.
if (hermaphrodites/=0)  call CheckSelfed(-s,k)

if (Complx==0) then
  if (Parent(A,3-k)/=0)  DumMate(s,k) = Parent(A,3-k)
  if (Parent(A,3-k) < 0) then
    DumMate(-Parent(A,3-k), 3-k) = -s
  else if (Parent(A,3-k) > 0) then
    Mate(Parent(A,3-k)) = -s
  endif
endif

end subroutine NewSibship

! #####################################################################

subroutine CheckDropSibship(s, k, Drop) 
use Global
implicit none

integer, intent(IN) :: s, k
logical, intent(OUT) :: Drop
integer :: i                     

if (s > nC(k))  return  ! already dropped. 

Drop = .FALSE.
if (ns(s, k) == 0) then
  Drop = .TRUE.
else if (ALL(GpID(:,s,k)==0) .and. ns(s,k)==1) then
  if (DumClone(s,k)/=0 .or. Complx==0) then
    Drop = .FALSE.
  else
    Drop = .TRUE.
  endif
endif
if (.not. Drop)  return

if (ns(s,k)==1) then
  i = SibID(1, s, k)
  call RemoveSib(i, s, k)  
endif
call DoMerge(0, s, k)      ! delete sibship      

end subroutine CheckDropSibship

! #####################################################################

recursive subroutine DoMerge(SA, SB, k)  ! if SA=0, delete SB
use Global
implicit none

integer, intent(IN) :: SA, SB, k
integer :: j, n, m, x, y, BB(ns(SB,k)), nB
logical :: valid(2)         

if (SA == SB) return

valid = .TRUE.
if (SA/=0) then
  call ChkAncest(-SA,k, -SB,k, valid(1))
  call ChkAncest(-SB,k, -SA,k, valid(2))
  if (.not. (all(valid))) then
    call Erstop("Pedigree loop created by merge", .TRUE.)
  endif
endif 

if (SA/=0) then
  if (nS(SA,k) + nS(SB,k) >= maxSibSize) then
    call Erstop("Reached Maximum Sibship Size (number of offspring per parent), please increase '--maxsibsize'", .FALSE.)
  endif 
  nB = ns(SB,k)
  BB = SibID(1:ns(SB,k), SB, k)
  do j=1,nB
    if (nB==0)  exit
    call setPar(BB(j), 3, -SA, k)
  enddo
  ! dummy offspring fixed below

  do m=1,2
    if (GpID(m, SA, k)==0 .and. GpID(m, SB, k)/=0) then  ! checked for mismatches earlier
      call setPar(-SA,k, GpID(m,SB,k), m)   ! takes care of updating CLL, SClone, etc.
    endif  ! else keep GpID(i,SA,k)
  enddo
endif

do x=SB, nC(k)  !remove cluster SB, shift all subsequent ones
  SibID(:, x, k) = SibID(:, x+1, k)  
  nS(x, k) = nS(x+1, k)
  GpID(:, x,k) = GpID(:, x+1,k)
  do n=1, nS(x,k)
    Parent(SibID(n,x,k),k) = -x 
  enddo
  CLL(x,k) = CLL(x+1, k)
  XPr(:,:,:,x,k) = XPr(:,:,:,x+1,k)
  DumP(:,:,x,k) = DumP(:,:,x+1,k)
  DumBY(:,x,k,:) = DumBY(:,x+1,k,:)
  IsNewSibship(x,k) = IsNewSibship(x+1, k)
  DumMate(x,k) = DumMate(x+1, k)
  DumClone(x,k) = DumClone(x+1, k)
enddo
SibID(:,nC(k),k) = 0
GpID(:,nC(k),k) = 0
nS(nC(k), k) = 0
DumMate(nC(k), k) = 0  
DumClone(nC(k), k) = 0  

do x=SB, nC(k)
  if (Complx == 0) then
    if (any(Mate == -x .and. Sex==3-k)) then
      y = MINLOC(ABS(Mate + x), DIM=1, MASK = Sex==3-k)
      if (x==-SB) then
        Mate(y) = 0
      else
        Mate(y) = -x+1   ! shift towards zero.
      endif
    else if (any(DumMate(:,3-k) == -x)) then
      y = MINLOC(ABS(DumMate(:,3-k) + x), DIM=1)
      if (x==-SB) then
        DumMate(y, 3-k) = 0
      else
        DumMate(y, 3-k) = -x+1
      endif
    endif
  endif
  if (hermaphrodites/=0) then
    if (any(DumClone(:,3-k) == x)) then
      y = MINLOC(ABS(DumClone(:,3-k) - x), DIM=1)
      if (x==SB) then
        DumClone(y, 3-k) = 0  
      else
        DumClone(y, 3-k) = x-1
      endif
    endif
  endif
enddo

do m=1,2  !fix GPs
  do n=1, nC(m)
    if (GpID(k, n, m) == -SB) then
      GpID(k, n, m) = -SA
      if (all(GpID(:,n,m)==0) .and. ns(n,m)==1) then
        j = SibID(1,n,m)
        call RemoveSib(j,n,m)
        call DoMerge(0, n, m)   ! Recursive
      endif
    endif
    do x=SB+1, nC(k)  
      if (GpID(k, n, m) == -x)  GpID(k, n, m) = -x+1 
    enddo
  enddo
enddo
nC(k) = nC(k) -1

if (SA/=0) then
  IsNewSibship(SA,k) = .TRUE.
  ToCheck(SibID(1:ns(SA,k),SA,k)) = .TRUE. 
  call setEstBY(-SA, k)                       
  if (hermaphrodites/=0) then
    call CheckSelfed(-SA,k) 
    if (all(SelfedIndiv(SibID(1:ns(SA,k),SA,k)))) then
      ToCheck(SibID(1:ns(SA,k),SA,k)) = .FALSE. 
    endif
  endif
endif

end subroutine DoMerge

! #####################################################################

subroutine DoFSMerge(SA, SB, k)   ! merge via k .and. k-3
use Global
implicit none

integer, intent(IN) :: SA, SB, k
integer :: i, ParA, ParB

! assume all checks have been done beforehand
! not implemented yet: Par(1) > 0, Par(2) <= 0 or vv
 call getFSpar(SA, k, .TRUE., ParA)
 call getFSpar(SB, k, .TRUE., ParB)
 
if (ParA==0 .and. any(parent(SibID(1:ns(SA,k),SA,k), 3-k)/=0)) then
  ParA = 9999
else if (ParA/=0 .and. any(parent(SibID(1:ns(SA,k),SA,k), 3-k)==0)) then
  do i=1, nS(SA,k)
    if (parent(SibID(i,SA,k), 3-k)==0) then
      call SetPar(SibID(i, SA, k), 3, parA, 3-k)
    endif
  enddo
endif
if (ParB==0 .and. any(parent(SibID(1:ns(SB,k),SB,k), 3-k)/=0)) then
  ParB = 9999
else if (ParB/=0 .and. any(parent(SibID(1:ns(SB,k),SB,k), 3-k)==0)) then
  do i=1, nS(SB,k)
    if (parent(SibID(i,SB,k), 3-k)==0) then
      call SetPar(SibID(i, SB, k), 3, parB, 3-k)
    endif
  enddo
endif
 
if (ParA < 0 .and. ParB < 0) then
  call DoMerge(-ParA, -ParB, 3-k)
else if (ParA==0 .and. ParB==0) then
  call NewSibship(SibID(1,SA,k), 0, 3-k)
  do i=2, ns(SA,k)
    call setPar(SibID(i,SA,k), 3, -nC(3-k), 3-k)
  enddo
  do i=1, nS(SB, k)
    call setPar(SibID(i,SB,k), 3, -nC(3-k), 3-k)
  enddo
  call CalcCLL(nC(3-k), 3-k)
  call CalcCLL(SA, k)
  call CalcCLL(SB, k)
  call CalcCLL(nC(3-k), 3-k)
else if (ParA < 0 .and. ParB == 0) then
  do i=1, nS(SB, k)
    call setPar(SibID(i,SB,k), 3, ParA, 3-k)
  enddo
else if (ParB < 0 .and. ParA == 0) then
  do i=1, nS(SA, k)
    call setPar(SibID(i,SA,k), 3, ParB, 3-k)
  enddo
! else not implemented yet
endif

 call DoMerge(SA, SB, k)  ! takes care of MakeFS

end subroutine DoFSMerge

! #####################################################################

subroutine getOff(P, kP, dums, nOff, Off, sxOff)  ! list all offspring for parent P
use Global
implicit none

integer, intent(IN) :: P, kP
logical, intent(IN) :: dums  ! include dummy offspring
integer, intent(OUT) :: nOff, Off(maxSibSize), sxOff(maxSibSize)
integer :: i, k, m, s

nOff = 0
Off = 0
sxOff = 3                                                
if (P==0) return                                
do k=1,2
  if (P>0 .and. kP/=1 .and. kP/=2) then
    if (Sex(P)<3 .and. Sex(P)/=k) cycle
  else if (k/=kP) then 
    cycle
  endif
  do i=1, nInd
    if (Parent(i,k) == P) then
      nOff = nOff + 1
      Off(nOff) = i
      sxOff(nOff) = Sex(i)
    endif
    if (nOff == maxSibSize) then
      call Erstop("Reached Maximum Sibship Size (number of offspring per parent), please increase '--maxsibsize'", .FALSE.)
    endif
  enddo
  if (dums) then
    do m=1,2
      do s=1,nC(m)
        if (GpID(k,s,m) == P) then
          nOff = nOff + 1
          Off(nOff) = -s
          sxOff(nOff) = m 
        endif
        if (nOff == maxSibSize) then
          call Erstop("Reached Maximum Sibship Size (number of offspring per parent), please increase '--maxsibsize'", .FALSE.)
        endif
      enddo
    enddo
  endif
enddo

end subroutine getOff

! #####################################################################

subroutine CalcU(A, kAIN, B, kBIN, LL)  ! A, SB, k, SA, LL
use Global
use CalcLik
implicit none

integer, intent(IN) :: A, kAIN, B, kBIN
double precision, intent(OUT) :: LL
integer :: kA, kB, Ai, Bj, SA, SB, cat, m, n, par(2), i, tmpGP
logical :: con, swap, OpG, conP
double precision :: AgeD

LL = missing
con = .FALSE.
kA = 0
kB = 0
if (A>0) then
  call CalcLind(A)
else if (A<0) then
  kA = kAIN     
  call CalcCLL(-A, kA)
endif
if (B>0) then
  call CalcLind(B)
else if (B<0) then
  kB = kBIN
  call CalcCLL(-B, kB)
endif
!==================================

if (A==0) then
  if (B==0) then
    LL = 0D0
  else if (B>0) then
    LL = Lind(B)
  else if (B<0) then
    LL = CLL(-B, kB)
  endif
  return
else if (B==0) then
  if (A>0) then
    LL = Lind(A)
  else if (A<0) then
    LL = CLL(-A,kA)
  endif
  return
else if (A>0 .and. B<0) then
  if (Parent(A,kB)==B) then
    LL = CLL(-B,kB)
    return
  else if (ANY(GpID(:,-B,kB) == A)) then
    LL = CLL(-B,kB) + Lind(A)  ! CLL already conditional on A
    return
  else if (ALL(Parent(A,:)>=0)) then  
    LL = Lind(A) + CLL(-B, kB)
    return
  else
    call Connected(A,1,B,kB, con)
    if (.not. con) then
      LL = Lind(A) + CLL(-B, kB)
      return
    endif
  endif
else if (B>0 .and. A<0) then
  if (Parent(B,kA)==A) then
    LL = CLL(-A, kA)
    return
  else if (ANY(GpID(:,-A,kA) == B)) then
    LL = CLL(-A,kA) + Lind(B)  
    return
  else if (ALL(Parent(B,:)>=0)) then 
    LL = CLL(-A, kA) + Lind(B) 
    return
  else
    call Connected(B,1,A,kA, con)
    if (.not. con) then
      LL = CLL(-A, kA) + Lind(B) 
      return
    endif
  endif
endif

!==================================
! determine relationship between focal individuals/clusters

Ai = 0
Bj = 0
SA = 0
SB = 0
cat = 0
swap = .FALSE.

if (A>0 .and. B>0) then  ! == pairs ==
  do m=1,2
    if (Parent(A,m)/=0 .and. Parent(A,m)==Parent(B,m)) then
       par(m) = Parent(A,m)
    else
      par(m) = 0  ! unknown or unequal
    endif
  enddo

  if (par(1)/=0 .and. par(2)/=0) then
    cat = 2  ! FS
  else if (par(1)/=0 .or. par(2)/=0) then
    cat = 3  ! HS
  else 
    do m=1,2
      if (parent(A,m) < 0) then
        do n=1,2
          if (GpID(n, -parent(A,m), m) == B) then
            cat = 0 !4  ! already conditioned on.
          else
            if (GpID(n, -parent(A,m), m)==Parent(B, n) .and. &
             Parent(B, n)<0) then
              cat = 5
            endif
          endif
        enddo
      else if (parent(B,m) < 0) then
        if (ANY(GpID(:, -parent(B,m), m) == A)) then
          cat = 0 !4
          swap = .TRUE.
        else
          do n=1,2
            if (GpID(n, -parent(B,m), m) == Parent(A, n) .and. &
             Parent(A, n)<0) then
              cat = 5
              swap = .TRUE.
            endif
          enddo
        endif
      endif
    enddo
  endif

  if (cat==0 .or. cat==5) then  ! TODO? cat=5
    LL = Lind(A) + Lind(B)
    return
  else if (cat==2 .and. par(1)<0 .and. par(2)<0) then
    Ai = A
    Bj = B
    SA = -par(1)
    kA = 1
    SB = -par(2)
    kB = 2
    cat = 0
  else
    call Upair(A, B, cat, LL)
    return
  endif

else if (A>0 .and. B<0) then
  SB = -B
  Ai = A
  if (ALL(Parent(A,:) < 0)) then
    SA = -Parent(A,3-kB)
    kA = 3-kB
    do m=1,2
      do i=1,ns(-Parent(A,m),m)
        if (Parent(SibID(i,-Parent(A,m),m), 3-m)/=Parent(A,3-m)) then
          call Connected(SibID(i,-Parent(A,m),m),m,B,kB, conP)
          if (conP) then
            SA = -parent(A,m)
            kA = m
            exit
          endif
        endif
      enddo
    enddo
  else if (Parent(A,3-kB) < 0) then
    SA = -Parent(A,3-kB)
    kA = 3-kB    
  else if (Parent(A,kB) < 0) then
    SA = -Parent(A,kB)
    kA = kB
  endif ! else: Lind + CLL (earlier) 
else if (B>0 .and. A<0) then
  SA = -A
  Bj = B
  if (ALL(Parent(B,:) < 0)) then
    SB = -Parent(B,3-kA)
    kB = 3-kA
    do m=1,2
      do i=1,ns(-Parent(B,m),m)
        if (Parent(SibID(i,-Parent(B,m),m), 3-m)/=Parent(B,3-m)) then
          call Connected(SibID(i,-Parent(B,m),m),m,A,kA, conP)
          if (conP) then
            SB = -parent(B,m)
            kB = m
            exit
          endif
        endif
      enddo
    enddo
  else if (Parent(B,3-kA) < 0) then
    SB = -Parent(B,3-kA)
    kB = 3-kA
  else if (Parent(B,kA) < 0) then
    SB = -Parent(B,kA)
    kB = kA 
  endif
else if (A<0 .and. B<0) then
  SA = -A  
  SB = -B
endif

cat = 0
if (GpID(kB, SA, kA) == -SB) then
  cat = 1  ! PO
else if (GpID(kA, SB, kB) == -SA) then
  cat = 1
  swap = .TRUE.
else 
  do m=1,2
    if (GpID(m, SA, kA)==GpID(m, SB, kB) .and. GpID(m, SA, kA)/=0) then
     if (GpID(3-m,SA,kA)==GpID(3-m,SB,kB) .and. GpID(3-m,SA,kA)/=0) then  
        cat = 2  ! FS
      else
        cat = 3  ! HS
      endif
    else 
      if (GpID(m, SA, kA)<0) then
        if (GpID(kB, -GpID(m, SA, kA), m) == -SB) then
          cat = 4  ! GP
        endif
      endif
      if (GpID(m, SB, kB)<0) then
        if (GpID(kA, -GpID(m, SB, kB),m) == -SA) then
          cat = 4
          swap = .TRUE.
        endif
      endif
    endif
  enddo  ! FA between SA, SB not currently considered.
endif

OpG = .FALSE.
if (con .and. cat==0) then
  if (A<0 .and. B>0) then
    do i=1, ns(-A,kA)
      if (Parent(SibID(i,-A,kA), 3-kA) < 0) then
        if (ANY(GpID(:,-Parent(SibID(i,-A,kA), 3-kA), 3-kA)==B)) then
          SB = -Parent(SibID(i,-A,kA), 3-kA)
          kB = 3-kA
          Bj = 0
          OpG = .TRUE.
        endif
      endif
    enddo
else if (A>0 .and. B<0) then
  do i=1, ns(-B,kB)
    if (Parent(SibID(i,-B,kB), 3-kB) < 0) then
        if (ANY(GpID(:,-Parent(SibID(i,-B,kB), 3-kB), 3-kB)==A)) then
          SA = -Parent(SibID(i,-B,kB), 3-kB)
          kA = 3-kB
          Ai = 0
          OpG = .TRUE.
          swap = .TRUE.
        endif
      endif
    enddo
  endif
endif

if (cat==0) then ! swap if BY(A) < BY(B)   (A older)
  call EstAgeDif(A, kA, B, kB, AgeD)  ! = AgeDiff(A,B) if BY(A) & BY(B) both exactly known
  if (AgeD < 0D0) then
    swap = .TRUE.
  endif
endif

if (con .and. A<0 .and. B<0 .and. kA==kB) then
  do i=1,ns(-B,kB)
    if (Parent(SibID(i,-B,kB), 3-kB) < 0) then
      if (GpID(3-kA, -Parent(SibID(i,-B,kB), 3-kB), 3-kB) < 0) then
        tmpGP = GpID(3-kA, -Parent(SibID(i,-B,kB), 3-kB), 3-kB)
        if (ANY(Parent(SibID(1:ns(-A,kA),-A,kA), 3-kA) == tmpGP) .and. &
         .not. ANY(Parent(SibID(1:ns(-B,kB),-B,kB), 3-kB) == tmpGP)) then
          swap = .TRUE.
        endif
      endif
    endif
  enddo
endif

if (.not. swap) then
  call UClust(-SA, -SB, kA, kB, cat, Ai, Bj, LL)
else
  call UClust(-SB, -SA, kB, kA, cat, Bj, Ai, LL)
endif

if (opG) then
  if (B>0) then
    LL = LL - CLL(SB,kB) + Lind(B)
  else if (A>0) then
    LL = LL - CLL(SA,kA) + Lind(A)
  endif
endif

end subroutine CalcU

! #####################################################################

subroutine Upair(A, B, cat, LL)
use Global
implicit none

integer, intent(IN) :: A, B, cat
double precision, intent(OUT) :: LL
integer :: m, l, n, x, y, par(2)
double precision :: PrL(nSnp), PrP(3,2), PrPA(3), PrPB(3), PrXY(3,3)

LL = missing
do m=1,2
  if (Parent(A,m)/=0 .and. Parent(A,m)==Parent(B,m)) then
     par(m) = Parent(A,m)
  else
    par(m) = 0  ! unknown or unequal
  endif
enddo
      
PrL = 0D0
do l=1, nSnp  
  if (cat==2) then
    do m=1,2
      call ParProb(l, Par(m), m, A, B, PrP(:,m))
    enddo
    do x=1,3
      do y=1,3
        PrXY(x,y) = OKA2P(Genos(l,A),x,y) * OKA2P(Genos(l,B),x,y) * &
          PrP(x,1) * PrP(y,2)
      enddo
    enddo
  else if (cat==3) then  ! HS
    do m=1,2
      if (Par(m)==0) cycle
      call ParProb(l, Par(m), m, A, B, PrP(:,m))
      call ParProb(l, Parent(A, 3-m), 3-m, A, 0, PrPA)
      call ParProb(l, Parent(B, 3-m), 3-m, B, 0, PrPB)
      do x=1,3  ! shared parent
        do y=1,3  ! parent A
          PrXY(x,y) = OKA2P(Genos(l,A),x,y) * PrP(x,m) * PrPA(y) * &
             SUM(OKA2P(Genos(l,B),x,:) * PrPB)
        enddo
      enddo
    enddo
  else if (cat==4) then
    do m=1,2
      if (Parent(A,m)<0) then
        do n=1,2
          if (GpID(n, -parent(A,m), m) == B) then
            call ParProb(l, parent(A,m), m, A, -4, PrP(:,m))  
            call ParProb(l, parent(A,3-m), 3-m, A, 0, PrPA)
            call ParProb(l, GpID(3-n, -parent(A,m), m), 3-n, 0, 0, PrPB)
            call ParProb(l, B, n, 0, 0, PrP(:,3-m))
            do x=1,3  ! in-between parent
              do y=1,3  ! other parent of A
                PrXY(x,y) = SUM(OKA2P(Genos(l,A),x,:) *PrPA) *PrP(x,m)*&
                   SUM(AKA2P(x, y,:) * PrP(y,3-m) * PrPB)
              enddo
            enddo
          endif
        enddo
      endif
    enddo     
  endif
  PrL(l) = LOG10(SUM(PrXY))
enddo
LL = SUM(PrL)

end subroutine Upair

! #####################################################################

subroutine UClust(A, B, kA, kB, cat, Ai, Bj, LL)
use Global
use CalcLik
implicit none

integer, intent(IN) :: A, B, kA, kB, cat, Ai, Bj
double precision, intent(OUT) :: LL
integer :: nA, AA(ns(-A,kA)), GA(2), nB, BB(ns(-B,kB)), GB(2), &
  AB(2*maxSibSize), UseEE(ns(-A,kA)+ns(-B,kB)), TypeEE(ns(-A,kA)+ns(-B,kB)), &
  MateABpar(ns(-A,kA)+ns(-B,kB)), catA(maxSibSize), catB(maxSibSize), &
  l,x,y, v, i, j, z,m, f, e, u, DoneA(maxSibSize), g, Ei
double precision :: PrL(nSnp,2), PrGA(3,2), PrGB(3,2), PrGGP(3), PrFS(3,3), &
  PrUZ(3,3, 3,3,3,3,2), PrE(3), PrH(3), PrW(3), PrEE(3,(ns(-A,kA)+ns(-B,kB)))
logical :: ParAisClone(maxSibSize), ParBisClone(maxSibSize), AisBclone, &
  SIMPL, DoRSibs(maxSibSize, 2)

LL = missing

nA = nS(-A, kA)
AA = SibID(1:nA, -A, kA)
GA = GpID(:, -A, kA)

nB = nS(-B, kB)
BB = SibID(1:nB, -B, kB)
GB = GpID(:, -B, kB)

!============================================

AB = 0         
UseEE = 0
TypeEE = 0
MateABpar = 0
if (kA==kB) then
  AB(1:nB) = BB
  AB((nB+1):(nB+nA)) = AA
  call FindEE(AB(1:(nB+nA)), nB, nA, kB, UseEE, MateABpar) 
  BB = AB(1:nB)
  AA = AB((nB+1):(nB+nA))
  TypeEE = 3-kB
  do i=1, nB  ! safety net
    if (UseEE(i) > nB) then
      UseEE(i) = 0  ! else use before store
    endif
  enddo
else if (kA/=kB) then
  call FindEE(BB, nB, 0, kB, UseEE(1:nB), MateABpar(1:nB))  ! may reorder BB
  call FindEE(AA, nA, 0, kA, UseEE((nB+1):(nB+nA)), MateABpar((nB+1):(nB+nA)))
  do i=1, nA
    if (UseEE(nB+i)/=0) then
      UseEE(nB+i) = nB + UseEE(nB+i)
    endif
  enddo
  TypeEE(1:nB) = 3-kB
  TypeEE((nB+1):(nB+nA)) = 3-kA
endif

!============================================
catA = 0
catB = 0
do i = 1, nA
  if (kA /= kB .and. Parent(AA(i), kB) == B) then
    catA(i) = 1
    UseEE(nB+i) = 0
  else if (GA(3-kA) == Parent(AA(i), 3-kA) .and. GA(3-kA) /= 0) then
    catA(i) = 2
  else if (GB(3-kA) == Parent(AA(i), 3-kA) .and. GB(3-kA) /= 0) then
    catA(i) = 3
  else if (Parent(AA(i), 3-kA) < 0) then
    if (GpID(kA, -Parent(AA(i), 3-kA), 3-kA)==A) then
      catA(i) = 6
      UseEE(nB+i) = 0
    else if (GpID(kB, -Parent(AA(i), 3-kA), 3-kA)==B) then
      catA(i) = 8
      UseEE(nB+i) = 0
    endif
  endif
enddo

do j=1, nB
  if (kA /= kB .and. Parent(BB(j), kA) == A) then
    catB(j) = 1
    UseEE(j) = 0
  else if (GA(3-kB) == Parent(BB(j), 3-kB) .and. GA(3-kB) /= 0) then
    catB(j) = 2
  else if (GB(3-kB) == Parent(BB(j), 3-kB) .and. GB(3-kB) /= 0) then
    catB(j) = 3  
  else if (Parent(BB(j), 3-kB) < 0) then
    if (GpID(kA, -Parent(BB(j),3-kB), 3-kB)==A) then
      catB(j) = 6
      UseEE(j) = 0
    else if (GpID(kB, -Parent(BB(j),3-kB), 3-kB)==B) then
      catB(j) = 8
      UseEE(j) = 0
    endif 
  endif
enddo

do i = 1, nA
  do j = 1, nB
    if (kA == kB) then
      if (Parent(AA(i), 3-kA) == Parent(BB(j), 3-kB) .and. Parent(BB(j), 3-kB)<0) then 
        catA(i) = 7
        if (catB(j)==0)  catB(j) = 7
      endif
    endif
  enddo
enddo

ParAisClone = .FALSE.
ParBisClone = .FALSE.
AisBclone = .FALSE.
if (hermaphrodites /= 0) then
  if (kA/=kB .and. DumClone(-A,kA) == -B) then
    AisBclone = .TRUE.
  else
    if (DumClone(-A,kA)/=0) then
      do i = 1, nA
        if (Parent(AA(i), 3-kA) == -DumClone(-A,kA))   ParAisClone(i) = .TRUE.
      enddo 
    endif
    if (DumClone(-B,kB)/=0) then
      do j=1,nB
        if (Parent(BB(j), 3-kB) == -DumClone(-B,kB))   ParBisClone(j) = .TRUE.
      enddo
    endif
  endif
endif

!==================================
if (cat==0 .and. ALL(catA==0) .and. ALL(CatB==0) .and. ALL(UseEE==0) .and. &
  .not. AisBclone) then
  LL = CLL(-A,kA) + CLL(-B,kB)
  return
endif
!==================================

SIMPL = ALL(catA==0) .and. ALL(catB==0) .and. Ai==0 .and. Bj==0 .and. &
    ALL(UseEE==0) .and. .not. AisBclone .and. .not. any(ParAisClone) .and. &
    .not. any(ParBisClone)

DoRsibs = .TRUE. 
if (.not. SIMPL) then
  call ChkTooManySibs(AA, nA, kA, DoRsibs(:,1))
  call ChkTooManySibs(BB, nB, kB, DoRsibs(:,2))
endif

PrL = 0D0
do l=1, nSnp
  PrUZ = 0D0
  
  do m=1, 2
    if ((ANY(catA==2) .and. m/=kA) .or. (ANY(catB==2) .and. m/=kB)) then
      call ParProb(l, GA(m), m, -1, 0, PrGA(:, m))
    else
      call ParProb(l, GA(m), m, 0, 0, PrGA(:, m))
    endif
    if ((ANY(catA==3) .and. m/=kA) .or. (ANY(catB==3) .and. m/=kB)) then
      call ParProb(l, GB(m), m, -1, 0, PrGB(:, m))
    else
      call ParProb(l, GB(m), m, 0, 0, PrGB(:, m))
    endif
  enddo
  if (cat==4) then
    do m=1,2
      if (GA(m)<0) then
        if (GpID(kB, -GA(m), m) == B) then
          call ParProb(l, GpID(3-kB, -GA(m), m), 3-kB, 0, 0, PrGGP) 
        endif
      endif
    enddo
  endif
  
  ! == grandparents ==
  do x=1,3 
    do y=1,3
      if (AisBclone .and. y/=x)  cycle
      do u=1,3  ! GP A, kB
        do z=1,3  ! GP A, 3-kB
          do v=1,3  ! GP B, kB
            if (cat == 1) then
             if (kA==kB .and. GA(3-kB)==GB(3-kB) .and. GA(3-kB)/=0) then 
                PrUZ(x,y,y,z,v,z,1) = AKA2P(x,y,z) * AKA2P(y,v,z) *&
                 PrGA(z,3-kB) * PrGB(v,kB)
              else
                PrUZ(x,y,y,z,v,:,1) = AKA2P(x,y,z) * AKA2P(y,v,:) * &
                  PrGA(z,3-kB) * PrGB(v,kB) * PrGB(:,3-kB)
              endif
            else if (cat==2) then
              PrUZ(x,y,u,z,u,z,1) = AKA2P(x,u,z) * AKA2P(y,u,z) *&
               PrGA(u,kB) * PrGA(z,3-kB)
            else if (cat==3) then
              do m=1,2
                if (GA(m)/=0 .and. GA(m) == GB(m)) then
                  if (m==kB) then
                    PrUZ(x,y,u,z,u,:,1) = AKA2P(x,u,z) * AKA2P(y,u,:) *&
                       PrGA(u,m) * PrGA(z,3-m) * PrGB(:,3-m)   
                  else
                    PrUZ(x,y,u,z,v,z,1) = AKA2P(x,u,z) * AKA2P(y,v,z) *&
                       PrGA(u,3-m) * PrGA(z,m) * PrGB(v,3-m)
                  endif
                endif
              enddo
            else if (cat==4) then 
              do m=1,2
                if (GA(m)<0) then
                  if (GpID(kB, -GA(m), m) == B) then
                    if (m==kB) then
                      PrUZ(x,y,u,z,v,:,1) =AKA2P(x,u,z) *&
                       SUM(AKA2P(u,y,:) *PrGGP) *PrGA(z,3-kB) *&
                       AKA2P(y,v,:) * PrGB(v,kB) * PrGB(:, 3-kB) 
                    else
                      PrUZ(x,y,u,z,v,:,1) = AKA2P(x,u,z) *&
                     SUM(AKA2P(z,y,:) *PrGGP) *PrGA(u,kB) *AKA2P(y,v,:)&
                     * PrGB(v,kB) * PrGB(:, 3-kB)
                    endif
                  endif
                endif
              enddo
            else
              PrUZ(x,y,u,z,v,:,1) = AKA2P(x,u,z) * AKA2P(y,v,:) * &
                PrGA(u,kB) * PrGA(z,3-kB) * PrGB(v,kB) * PrGB(:, 3-kB)
            endif
            PrUZ(x,y,u,z,v,:,2) = PrUZ(x,y,u,z,v,:,1)
          enddo
        enddo
      enddo
    enddo
  enddo
   
  ! == siblings ==   
  if (SIMPL) then
    do x=1,3 
      do y=1,3
        PrUZ(x,y,:,:,:,:,2) = PrUZ(x,y,:,:,:,:,2) * XPr(1,x,l,-A,kA) *&
         XPr(1,y,l,-B,kB)
      enddo  ! TODO: needs special for cat<4 ?
    enddo
  
  else

  do x=1,3  ! SA
    doneA = 0
    do y=1,3  ! SB
      if (AisBclone .and. y/=x)  cycle
      PrEE = 0D0
      do j=1, nB
        if (nFS(BB(j))==0) cycle
        if (catB(j)==1 .or. catB(j)==2 .or. catB(j)==3 .or. ParBisClone(j)) then
          PrE = 1D0
        else if (catB(j)==6) then  
          call ParProb(l, GpID(3-kA,-Parent(BB(j),3-kB),3-kB),3-kA, 0, 0, PrH)
          do e=1,3
            PrE(e) = SUM(AKA2P(e,x,:) * PrH) 
          enddo
        else if (catB(j)==8) then  
          call ParProb(l, GpID(3-kB,-Parent(BB(j),3-kB),3-kB),3-kB, 0, 0, PrH)
          do e=1,3
            PrE(e) = SUM(AKA2P(e,y,:) * PrH) 
          enddo
        else if (UseEE(j)/=0) then
          call ParProb(l, MateABpar(j), 3-TypeEE(j), 0,0,PrH)
          do e=1,3
            do u=1, 3
              PrW(u) = SUM(AKA2P(e,u,:) * PrEE(u,UseEE(j)) * PrH)
            enddo
            PrE(e) = SUM(PrW)
          enddo
          PrE = PrE/SUM(PrE)
        else if (catB(j)==1 .or. catB(j)==2 .or. catB(j)==3) then
          PrE = 1D0
        else
          if (DoRsibs(j,2)) then
            call ParProb(l, Parent(BB(j),3-kB), 3-kB, -1, 0, PrE)
          else
            call ParProb(l, Parent(BB(j),3-kB), 3-kB, BB(j), -1, PrE)
          endif
        endif
        
        if (Parent(BB(j),3-kB) < 0 .and. catB(j)/=1 .and. DoRsibs(j,2)) then 
          do e=1,3
            do g=1, nS(-Parent(BB(j),3-kB), 3-kB)
              Ei = SibID(g, -Parent(BB(j),3-kB), 3-kB)
              if (nFS(Ei) == 0) cycle
              if (Parent(Ei, kB) == B) cycle  
              if (Parent(Ei, kA) == A) cycle
              call ParProb(l, Parent(Ei, kB), kB, Ei, -1, PrH) 
              PrFS = FSLik(l,Ei)
              PrH = PrH * PrFS(:,e)
              if (.not. all(PrH==1D0))  PrE(e) = PrE(e) * SUM(PrH)
            enddo
          enddo
        endif
        
        if (Bj/=0 .or. catB(j)==7 .or. (catB(j)==1 .and. Ai/=0)) then 
          do f=1, nFS(BB(j))
            if (Bj==0 .or. FSID(f, BB(j))==Bj) cycle
            if (Parent(BB(j),kA)==A .and. (Ai==0 .or. FSID(f, BB(j))==Ai)) cycle
            PrE = PrE * OKA2P(Genos(l,FSID(f,BB(j))), y, :)
          enddo
        endif
        
        if (catB(j)==7 .and. Ai/=0) then 
          do i=1,nA
            if (Parent(AA(i), 3-kB) /= Parent(BB(j), 3-kB)) cycle
            if (Parent(AA(i), kB) == B) cycle
            if (AA(i)==Ai)  cycle
            PrE = PrE * OKA2P(Genos(l,AA(i)), x, :)
          enddo
        endif
        
        if (catB(j)==1) then  ! Parent(BB(j), 3-kB)==PA
          PrUZ(x,y,:,:,:,:,1) = PrUZ(x,y,:,:,:,:,1) * PrE(x)
        else if (CatB(j)==2) then
          do z=1,3
            PrUZ(x,y,:,z,:,:,1) = PrUZ(x,y,:,z,:,:,1) * PrE(z)   
          enddo
        else if (CatB(j)==3) then
          do z=1,3
            PrUZ(x,y,:,:,:,z,1) = PrUZ(x,y,:,:,:,z,1) * PrE(z)   
          enddo
        else if (ParBisClone(j)) then
          PrUZ(x,y,:,:,:,:,1) = PrUZ(x,y,:,:,:,:,1) * PrE(y)  
        else if (.not. all(PrE==1D0)) then
          PrUZ(x,y,:,:,:,:,1) = PrUZ(x,y,:,:,:,:,1) * SUM(PrE)
        endif
        
        do f=1, nFS(BB(j)) ! includes some AA if cat=1 
          if (Bj==0 .or. FSID(f, BB(j))==Bj .or. &
            (Parent(BB(j),kA)==A .and. (Ai==0 .or. FSID(f, BB(j))==Ai))) then
            PrE = PrE * OKA2P(Genos(l,FSID(f, BB(j))), y, :)
          endif
        enddo

        if (catB(j)==7) then 
          do i=1,nA
            if (Parent(AA(i), 3-kB) /= Parent(BB(j), 3-kB)) cycle
            if (Parent(AA(i), kB) == B) cycle
            if (Ai/=0 .and. AA(i)/=Ai)  cycle
            PrE = PrE * OKA2P(Genos(l,AA(i)), x, :)
            DoneA(i) = 1
          enddo
        endif
        
        if (catB(j)==1) then  ! Parent(BB(j), 3-kB)==PA
          PrUZ(x,y,:,:,:,:,2) = PrUZ(x,y,:,:,:,:,2) * PrE(x)
        else if (CatB(j)==2) then
          do z=1,3
            PrUZ(x,y,:,z,:,:,2) = PrUZ(x,y,:,z,:,:,2) * PrE(z)   
          enddo
        else if (CatB(j)==3) then
          do z=1,3
            PrUZ(x,y,:,:,:,z,2) = PrUZ(x,y,:,:,:,z,2) * PrE(z)   
          enddo
        else if (ParBisClone(j)) then
          PrUZ(x,y,:,:,:,:,2) = PrUZ(x,y,:,:,:,:,2) * PrE(y)  
        else if (.not. all(PrE==1D0)) then
          PrUZ(x,y,:,:,:,:,2) = PrUZ(x,y,:,:,:,:,2) * SUM(PrE)
        endif
        PrEE(:,j) = PrE
      enddo  ! B_j
    
      do i=1, nA
        if (DoneA(i)==1) cycle                      
        if (nFS(AA(i))==0) cycle
        if (Parent(AA(i),kB)==B) cycle
        if ((catA(i)>1 .and. catA(i)<4) .or. ParAisClone(i)) then  ! catA==1 already done
          PrE = 1D0
        else if (catA(i)==6) then ! Parent(AA(i), 3-kA) <0
          call ParProb(l, GpID(3-kA,-Parent(AA(i),3-kA),3-kA),3-kA,0,0,PrH)
          do e=1,3
            PrE(e) = SUM(AKA2P(e,x,:) * PrH) 
          enddo
        else if (catA(i)==8) then ! Parent(AA(i), 3-kA) <0
          call ParProb(l, GpID(3-kB,-Parent(AA(i),3-kA),3-kA),3-kB,0,0,PrH)
          do e=1,3
            PrE(e) = SUM(AKA2P(e,y,:) * PrH) 
          enddo
        else if (UseEE(nB+i)/=0) then
          call ParProb(l, MateABpar(nB+i), 3-TypeEE(nB+i), 0,0,PrH)
          do e=1,3
            do u=1, 3
              PrW(u) = SUM(AKA2P(e,u,:) * PrEE(u,UseEE(nB+i)) * PrH)
            enddo
            PrE(e) = SUM(PrW)
          enddo
          PrE = PrE/SUM(PrE)
        else
          if (DoRsibs(i,1)) then
            call ParProb(l, Parent(AA(i), 3-kA), 3-kA, -1, 0, PrE)  
          else
            call ParProb(l, Parent(AA(i), 3-kA), 3-kA, AA(i), -1, PrE)  
          endif
        endif
        
        if (Parent(AA(i), 3-kA) < 0 .and. DoRsibs(i,1)) then 
          do e=1,3
            do g=1, nS(-Parent(AA(i), 3-kA), 3-kA)
              Ei = SibID(g, -Parent(AA(i), 3-kA), 3-kA)
              if (nFS(Ei) == 0) cycle
              if (Parent(Ei, kA) == A) cycle 
              if (Parent(Ei, kB) == B) cycle
              call ParProb(l, Parent(Ei, kA), kA, Ei, -1, PrH) 
              PrFS = FSLik(l,Ei)
              PrH = PrH * PrFS(:,e)
              if (.not. all(PrH==1D0))  PrE(e) = PrE(e) * SUM(PrH)
            enddo
          enddo
        endif
        
        if (Ai/=0) then
          do f=1, nFS(AA(i))
            if (FSID(f, AA(i))==Ai) cycle
            PrE = PrE * OKA2P(Genos(l,FSID(f,AA(i))), x, :)
          enddo
        endif
        
        if (catA(i)==2) then
          do z=1,3
            if (kA==kB) then
              PrUZ(x,y,:,z,:,:,1) = PrUZ(x,y,:,z,:,:,1) * PrE(z)
            else
              PrUZ(x,y,z,:,:,:,1) = PrUZ(x,y,z,:,:,:,1) * PrE(z) 
            endif
          enddo
        else if (catA(i)==3) then  
          do z=1,3
            if (kA==kB) then
              PrUZ(x,y,:,:,:,z,1) = PrUZ(x,y,:,:,:,z,1) * PrE(z)
            else
              PrUZ(x,y,:,:,z,:,1) = PrUZ(x,y,:,:,z,:,1) * PrE(z)
            endif
          enddo
        else if (ParAisClone(i)) then
          PrUZ(x,y,:,:,:,:,1) = PrUZ(x,y,:,:,:,:,1) * PrE(x)
        else if (.not. all(PrE==1D0)) then
          PrUZ(x,y,:,:,:,:,1) = PrUZ(x,y,:,:,:,:,1) * SUM(PrE)    
        endif
        
        do f=1, nFS(AA(i)) 
          if (Ai/=0 .and. FSID(f, AA(i))/=Ai) cycle
 !         DoneA(i)=2    ! for debuging only
          PrE = PrE * OKA2P(Genos(l,FSID(f,AA(i))), x, :)
        enddo
        
        if (catA(i)==2) then
          do z=1,3
            if (kA==kB) then
              PrUZ(x,y,:,z,:,:,2) = PrUZ(x,y,:,z,:,:,2) * PrE(z)
            else
              PrUZ(x,y,z,:,:,:,2) = PrUZ(x,y,z,:,:,:,2) * PrE(z)
            endif
          enddo
        else if (catA(i)==3) then  
          do z=1,3
            if (kA==kB) then
              PrUZ(x,y,:,:,:,z,2) = PrUZ(x,y,:,:,:,z,2) * PrE(z)
            else
              PrUZ(x,y,:,:,z,:,2) = PrUZ(x,y,:,:,z,:,2) * PrE(z)
            endif
          enddo
        else if (ParAisClone(i)) then
          PrUZ(x,y,:,:,:,:,2) = PrUZ(x,y,:,:,:,:,2) * PrE(x)
        else if (.not. all(PrE==1D0)) then
          PrUZ(x,y,:,:,:,:,2) = PrUZ(x,y,:,:,:,:,2) * SUM(PrE)    
        endif
        PrEE(:,nB+i) = PrE
      enddo  ! i
    enddo  ! x
  enddo  ! y
  endif
  do f=1,2
    PrL(l,f) = LOG10(SUM(PrUZ(:,:,:,:,:,:,f)))
  enddo
enddo
LL = SUM(PrL(:,2)) - SUM(PrL(:,1))


! if ((A==-47 .and. kA==1 .and. Ai==399 .and. B==-30 .and. kB==2) .or. &
  ! (B==-47 .and. kB==1 .and. Bj==399 .and. A==-30 .and. kA==2)) then
  ! write(*,'("UCLUST: ", 4i7, 3f8.1)') GpID(:,-A,kA), GpID(:,-B,kB),SUM(PrL(:,1)), SUM(PrL(:,2)), LL
  ! print *, Ai, Bj
  ! do i=1, nA
    ! write(*,'(i3,i5, ", ", i3,i5, 4i4)') nB+i, AA(i), catA(i), Parent(AA(i), 3-kA), &
      ! UseEE(nB+i), TypeEE(nB+i), MateABpar(nB+i), DoneA(i)
  ! enddo
  ! print *, "."
  ! do i=1, nB
    ! write(*,'(i3,i5, ", ", i3,i5, 3i4)') i, BB(i), catB(i), Parent(BB(i),3-kB), &
      ! UseEE(i), TypeEE(i), MateABpar(i)
  ! enddo
  ! print *, ""
! endif

end subroutine UClust

! #####################################################################

subroutine FindEE(AB, nA, nB, k, UseEE, MatePar)  ! find PO pairs among mates
use Global
use qsort_c_module
implicit none

integer, intent(IN) :: nA, nB, k
integer, intent(INOUT) :: AB(nA+nB)
integer, intent(OUT) :: UseEE(nA+nB), MatePar(nA+nB)
integer :: i, j,x, nAB(2), ABM(MAX(nA,nB),2), MateE(MAX(nA,nB), 2), GGK(2), &
  UseM(MAX(nA,nB),2), Order(2*maxSibSize), MateI
logical :: reorder, OrderAgain
double precision :: EEtmp(2*maxSibSize)

UseEE = 0
MatePar = 0

nAB = (/nA, nB/)                
ABM = 0
ABM(1:nA, 1) = AB(1:nA)
ABM(1:nB, 2) = AB((nA+1):(nA+nB))
MateE = 0
do x=1,2
  do i=1, nAB(x)
    if (nFS(ABM(i,x))==0 .and. nAB(x)>1)  cycle
    MateE(i,x) = Parent(ABM(i,x), 3-k)
  enddo
enddo
if ((nAB(1)==1 .and. nAB(2)<2) .or. COUNT(MateE < 0) < 2) return

GGK = 0
do x=1,2
  if (ABM(1,x)==0)  cycle
  if (Parent(ABM(1,x),k) < 0) then  ! else not called?
    GGK(x) = GpID(3-k, -Parent(ABM(1,x),k), k)   
  endif
enddo

! re-order AA and BB, so that PrE calculated before used
UseM = 0
reorder = .FALSE.
do x=1,2
  if (nAB(x)<=1) cycle
  do i=1, nAB(x)
    if (MateE(i,x) < 0) then 
      if (GpID(3-k, -MateE(i,x), 3-k) < 0 .and. &
       .not. ANY(GGK == GpID(3-k, -MateE(i,x), 3-k))) then
        do j=1, nAB(x)
          if (MateE(j,x) == GpID(3-k, -MateE(i,x), 3-k)) then
            UseM(i,x) = j
            if (j > i) reorder = .TRUE.
            exit
          endif
        enddo
      endif
    endif
  enddo
  
  if (reorder) then
    EEtmp(1:nAB(x)) = dble(UseM(1:nAB(x),x))
    Order = (/ (i, i=1, nAB(x), 1) /)
    call QsortC(EEtmp(1:nAB(x)), Order(1:nAB(x)))
    ABM(1:nAB(x),x) = ABM(Order(1:nAB(x)), x)
    UseM(1:nAB(x),x) = UseM(Order(1:nAB(x)),x)
    OrderAgain = .FALSE.
    do i=1, nAB(x)
      if (UseM(i,x) /= 0) then
        do j=1, nAB(x)
          if (UseM(i,x) == Order(j)) then
            UseM(i,x) = j
            if (j>i)  OrderAgain = .TRUE.
            exit
          endif
        enddo
      endif
    enddo
    if (OrderAgain) then
      EEtmp(1:nAB(x)) = dble(UseM(1:nAB(x),x))
      Order = (/ (i, i=1, nAB(x), 1) /)
      call QsortC(EEtmp(1:nAB(x)), Order(1:nAB(x)))
      ABM(1:nAB(x),x) = ABM(Order(1:nAB(x)), x)
    endif
  endif
enddo

AB = 0
AB(1:nA) = ABM(1:nA, 1)
AB((nA+1):(nA+nB)) = ABM(1:nB, 2)

do i=1, nA+nB
  if (nFS(AB(i))==0 .and. ((i<=nA .and. nA>1) .or. (i>nA .and. nB>1)))  cycle
  MateI = Parent(AB(i), 3-k)
  if (MateI < 0) then 
    if (GpID(3-k, -MateI, 3-k) < 0 .and. &
     .not. ANY(GGK == GpID(3-k, -MateI, 3-k))) then
      do j=1, i
      if (nFS(AB(j))==0 .and. ((j<=nA .and. nA>1) .or. (j>nA .and. nB>1)))  cycle
        if (Parent(AB(j), 3-k) == GpID(3-k, -MateI, 3-k)) then
          UseEE(i) = j
          MatePar(i) = GpID(k, -MateI, 3-k)
          exit
        endif
      enddo
    endif
  endif
enddo

end subroutine FindEE

! #####################################################################

subroutine AddSib(A, SB, k, LL)  
use Global
use CalcLik
implicit none

integer, intent(IN) :: A, SB, k
double precision, intent(OUT) :: LL
integer :: l, x, Bj, f, j, y, Ei, z, v, DoQuick
double precision :: PrL(nSnp), PrX(3), PrY(3), LRQ, LLtmp(2), LLU, PrZ(3), &
  PrYZ(3,3), PrXb(3,2), PrE(3), PrB(3,3), PrFS(3,3)
logical :: Inbr, ParOK

LL = missing
if (Parent(A,k)==-SB) then
  LL = AlreadyAss
else if (Parent(A,k)/=0) then
  LL = impossible
endif
if (LL/=missing) return

do f=1, nS(SB,k)
  if (ns(SB,k)==0)  exit
  Bj = SibID(f, SB, k)
  if (Parent(A, 3-k) /= 0) then
    if (Parent(Bj, 3-k) == Parent(A, 3-k)) then
      LL = impossible  ! use addFS() instead
    endif
  endif
  if (getAP (AgeDiff(A, Bj), 3, 0, k, Impossible) == Impossible) then  
    LL=impossible
  endif 
enddo
if (LL/=missing) return

call ChkValidPar(A,Sex(A), -SB,k, ParOK)
if (.not. ParOK) then
  LL = impossible
  return
endif

call Qadd(A, SB, k, LRQ)
if (LRQ < -HUGE(0D0)) then
  LL = impossible
  return
endif

Inbr = .FALSE.
if (Parent(A,3-k) < 0) then
  if (Parent(A,3-k) == GpID(3-k, SB, k)) then
    Inbr = .TRUE.  ! inbreeding loop created
  else if (GpID(k,-Parent(A,3-k),3-k) == -SB) then
    Inbr = .TRUE.
  endif
endif
do f=1, nS(SB,k)
  if (ns(SB,k)==0)  exit
  Bj = SibID(f, SB, k)
  if (Parent(A,3-k) == Bj)  Inbr = .TRUE.
  if (Parent(Bj,3-k) == A)  Inbr = .TRUE.
enddo

call ChkDoQuick(SB,k,DoQuick)

if ((Parent(A,3-k)<0 .and. DoQuick/=-2) .or. Inbr .or. DoQuick>1 .or. DoQuick==-3) then 
  if (Parent(A,3-k) < 0) then
    call CalcU(-SB, k, A, 3-k, LLU)
    call CalcU(-SB,k, Parent(A,3-k),3-k, LLtmp(1))
  else if (GpID(3-k,SB,k) < 0) then
    call CalcU(-SB, k, 0, 0, LLU)
    call CalcU(-SB,k, GpID(3-k,SB,k),3-k, LLtmp(1))
  endif
  call setParTmp(A,0,-SB,k)
  if (Parent(A,3-k) < 0) then
    call CalcU(-SB,k, Parent(A,3-k),3-k, LLtmp(2)) 
    LL = LLU + (LLtmp(2) - LLtmp(1))
  else if (GpID(3-k,SB,k) < 0) then
    call CalcU(-SB,k, GpID(3-k,SB,k),3-k, LLtmp(2))
    LL = LLU + (LLtmp(2) - LLtmp(1))
  else
    LL = CLL(SB,k)
  endif
  call setParTmp(A,0,0,k)
  if (GpID(3-k,SB,k) < 0)  call CalcCLL(-GpID(3-k,SB,k), 3-k)
  
else 
  PrL = 0D0
  
  if (DoQuick == -2) then    ! all FS
    do l=1,nSnp
      call ParProb(l, Parent(A,3-k), 3-k, A, 0, PrY)
      do j=1, nS(SB,k)
        Bj = SibID(j,SB,k)
        if (nFS(Bj)==0) cycle
        call ParProb(l, Parent(Bj,3-k), 3-k, -1, 0, PrZ)
        PrB = FSLik(l,Bj)
        do x=1,3
          do y=1,3
            do z=1,3
              PrYZ(y,z) = PrB(x,z) * PrZ(z) * OKA2P(Genos(l,A), x, y) * PrY(y) 
            enddo
          enddo
          PrX(x) = XPr(2,x,l, SB,k) * SUM(PrYZ)
        enddo
      enddo
      PrL(l) = LOG10(SUM(PrX))   
    enddo
  
  else if (abs(DoQuick) == 1) then
    do l=1,nSnp
      call ParProb(l, Parent(A,3-k), 3-k, A, 0, PrY)
      do x=1,3
        PrX(x) = XPr(3,x,l, SB,k) * SUM(OKA2P(Genos(l,A), x, :) * PrY)
      enddo
      PrL(l) = LOG10(SUM(PrX))   
    enddo

  else 
    do l=1,nSnp 
      do x=1,3
        PrXb(x,:) = XPr(2,x,l,SB,k)  ! GPs
        do j=1, nS(SB,k)
          Bj = SibID(j,SB,k)
          if (nFS(Bj)==0) cycle
          call ParProb(l, Parent(Bj,3-k), 3-k, -1, 0, PrZ)
          do z=1,3
            if (Parent(Bj,3-k)<0) then
              do v=1, nS(-Parent(Bj, 3-k), 3-k)
                Ei = SibID(v, -Parent(Bj, 3-k), 3-k)  
                if (NFS(Ei) == 0) cycle
                if (Parent(Ei, k) == -SB) cycle
                call ParProb(l, Parent(Ei, k), k, Ei,-1, PrE)
                PrFS = FSLik(l,Ei)
                PrE = PrE * PrFS(:,z)
                if (.not. ALL(PrE==1D0))  PrZ(z) = PrZ(z) * SUM(PrE)  
              enddo  
            endif
          enddo
          if (.not. ALL(PrZ==1D0))  PrXb(x,1) = PrXb(x,1) * SUM(PrZ) 
          PrB = FSLik(l,Bj)
          PrZ = PrZ * PrB(:,x)
          PrXb(x,2) = PrXb(x,2) * SUM(PrZ)
        enddo
      enddo

      call ParProb(l, Parent(A,3-k), 3-k, A, 0, PrY)
      do x=1,3
        PrXb(x,2) = PrXb(x,2) * SUM(OKA2P(Genos(l,A), x, :) * PrY)
      enddo
      PrL(l) = LOG10(SUM(PrXb(:,2))) - LOG10(SUM(PrXb(:,1)))  
    enddo
  endif
  
  LL = SUM(PrL)
endif

! if ((A==169 .and. any(SibID(:,SB,k)==168)) .or. (A==168 .and. any(SibID(:,SB,k)==169)) )  then
 ! open (unit=42,file="log.txt",status="unknown", position="append")  
 ! write (42, *) ""
  ! write (42, '("addSib?", 3i6, ", ", l4, i4, f9.2)') A, SB, k, Inbr, DoQuick , LL
  ! close(42)
! endif

end subroutine AddSib

! #####################################################################

subroutine AddSibInbr(A,SB,k,LL)
use Global
use CalcLik
implicit none

integer, intent(IN) :: A, SB, k
double precision, intent(OUT) :: LL(3)
integer :: l, x, y, GA(2), GG, Par, i, u, Bj, j, GX, Ei
double precision :: PrL(nSnp,3), PrXY(3,3), PrZ(3), PrPA(3), LLtmp(3), &
  ALR(3), LLU(4), PrXYU(3,3,3), PrLU(nSnp,3), PrE(3), PrFS(3,3)
logical :: maybe(3)

! 1: Par(Parent(A,3-k),k)=SB 
! 2: Parent(A,3-k)=GpID(3-k,SB,k)
! 3: as 1, A FS of B's (PA == DB)

LL = missing
maybe = .TRUE.
GA = getPar(Parent(A,3-k), 3-k)
if (GA(k) == -SB) then
  maybe(2) = .FALSE.
else if (GA(k) /= -SB .and. GA(k)/=0) then
  maybe(1) = .FALSE.
endif

if (hermaphrodites/=0) then
  LL = NotImplemented
  return
endif

GG = GpID(k,SB,k)
if (GpID(3-k,SB,k)/=Parent(A,3-k) .and. GpID(3-k,SB,k)/=0 .and. Parent(A,3-k)/=0) then
  maybe(2) = .FALSE.
else if (ANY(SibID(1:ns(SB,k),SB,k)==Parent(A,3-k))) then
  maybe(2) = .FALSE.
endif
if (.not. ANY(maybe(1:2))) then
  LL = impossible
  return
endif

Par = 0
if (maybe(1)) then
  call getFSpar(SB, k, .TRUE., Par)
  if (Par==0) then
    maybe(3) = .FALSE.
  else if (Parent(A,3-k)/=Par .and. Parent(A,3-k)/=0) then
    maybe(3) = .FALSE.
  else if (Par>0) then
    if (Parent(Par, k)/=SB .and. Parent(Par, k)/=0) then
      maybe(3) = .FALSE.
    endif
  else if (Par<0) then
    if (GpID(k, -Par, 3-k)/=SB .and. GpID(k, -Par, 3-k)/=0) then
      maybe(3) = .FALSE.
    endif
  endif
  if (maybe(3) .and. GA(3-k)==0) then
    GA = getPar(Par, 3-k)
  endif
else
  maybe(3) = .FALSE.
endif

call CalcAgeLR(Parent(A,3-k),3-k, -SB,k, 0,1, .TRUE., ALR(1))
call CalcAgeLR(-SB,k, Parent(A,3-k),3-k,  0,1, .TRUE., ALR(2))
call CalcAgeLR(Par,3-k, -SB,k, 0,1, .TRUE., ALR(3))
do x=1,3
  if (ALR(x) == impossible .or. ALR(x) < 3.0*TF)  maybe(x) = .FALSE.
enddo

if (.not. ANY(maybe)) then
  LL = impossible
  return
endif

GX = 0
if (maybe(2)) then
  if (Parent(A,3-k)/=0) then
    GX = Parent(A,3-k)
  else
    GX = GpID(3-k, SB, k)
  endif
endif

PrL = 0D0
PrLU = 0D0          
do l=1,nSnp
  if (maybe(1)) then
    if (Parent(A,3-k)>0) then
      call ParProb(l, GA(3-k), 3-k, Parent(A,3-k), 0, PrZ) 
      call ParProb(l, Parent(A,3-k),3-k,0,0,PrPA)
    else
      call ParProb(l, GA(3-k), 3-k, 0, 0, PrZ) 
      call ParProb(l, Parent(A,3-k),3-k,A,-4,PrPA)
    endif
    if (Parent(A,3-k)==0)   PrPA = 1D0
    do x=1,3
      do y=1,3
        PrXY(x,y) = OKA2P(Genos(l,A), x, y) * XPr(3,x,l, SB,k) * PrPA(y) * &
          SUM(AKA2P(y,x,:) * PrZ)
        do u=1,3
          PrXYU(x,y,u) = OKA2P(Genos(l,A), u, y) * XPr(3,x,l, SB,k) * PrPA(y) * &
             SUM(AKA2P(y,u,:) * PrZ) 
        enddo
      enddo
    enddo
    PrL(l,1) = LOG10(SUM(PrXY))   ! Parent(A,3-k) offspring of SB 
    PrLU(l,1) = LOG10(SUM(PrXYU))                          
  endif
  
 !===
  if(maybe(3)) then
    do i=1, ns(SB, k)
      if (nFS(SibID(i,SB,k))==0 .or. Parent(SibID(i,SB,k),3-k)/=Par)  cycle
      call ParProb(l, Par,3-k,SibID(i,SB,k),-5,PrPA)  ! exclude both GPs & Bi & FS of Bi
    enddo
    if (Par < 0) then
      if (Parent(A,3-k)==Par) then
        PrPA = PrPA/OKAP(Genos(l,A),:,l)
        PrPA = PrPA/SUM(PrPA)
      endif
      call ParProb(l, GA(3-k), 3-k, 0, 0, PrZ) 
    else
      call ParProb(l, GA(3-k), 3-k, Par, 0, PrZ) 
    endif
    do x=1,3
      do y=1,3
        PrXY(x,y) = PrPA(y) * SUM(AKA2P(y,x,:) * PrZ) * XPr(2,x,l, SB,k)  !Xpr(2,) = GP's
        do i=1, ns(SB, k)
          PrXY(x,y) = PrXY(x,y) * OKA2P(Genos(l,SibID(i,SB,k)), x, y)
        enddo
        PrXYU(x,y,:) = PrXYU(x,y,:) * OKAP(Genos(l,A), :, y)
        PrXY(x,y) = PrXY(x,y) * OKA2P(Genos(l,A), x, y)
      enddo
    enddo
    PrL(l,3) = LOG10(SUM(PrXY)) 
    PrLU(l,3) = LOG10(SUM(PrXYU))    
  endif
  
  !===
  if (maybe(2)) then
    call ParProb(l, GG, k, 0, 0, PrZ)
    call ParProb(l, GX,3-k, -1,0, PrPA)   
    PrXY = 1D0
    do x=1,3  ! SB
      do y=1,3  ! GX
        if (GX < 0) then
          do i=1, ns(-GX,3-k)
            Ei = SibID(i,-GX,3-k)
            if (Ei==A .or. Parent(Ei,k)==-SB .or. nFS(Ei)==0)  cycle
            call ParProb(l, Parent(Ei,k), k, Ei, -1, PrE)
            PrFS = FSLik(l,Ei)
            PrE = PrE * PrFS(:,y)
            if (.not. all(PrE==1D0))  PrXY(x,y) = PrXY(x,y) * SUM(PrE)
          enddo
        endif
        do j=1, nS(SB,k)
          Bj = SibID(j,SB,k)
          if (nFS(Bj)==0) cycle
          PrFS = FSLik(l,Bj)
          if (Parent(Bj,3-k)==GX .and. GX/=0) then
            PrXY(x,y) = PrXY(x,y) * PrFS(x,y)
          else
            call ParProb(l, Parent(Bj,3-k), 3-k, Bj, -1, PrE)
            PrE = PrE * PrFS(x,:)
            if (.not. all(PrE==1D0))  PrXY(x,y) = PrXY(x,y) * SUM(PrE)
          endif
        enddo

        if (GpID(3-k,SB,k)==0) then
          do u=1,3
            PrXYU(x,y,u) = PrXY(x,y) * SUM(AKAP(x,:,l) * PrZ) * PrPA(y) *  &
             OKA2P(Genos(l,A), u, y) * SUM(AKA2P(u,y,:) * AHWE(:,l)) 
          enddo
          PrXY(x,y) = PrXY(x,y) * SUM(AKA2P(x,y,:) * PrPA(y) * PrZ)
        else
          PrXY(x,y) = PrXY(x,y) * SUM(AKA2P(x,y,:) * PrPA(y) * PrZ)
          do u=1,3
            if (Parent(A,3-k)==GX .and. GX/=0) then
              PrXYU(x,y,u) = PrXY(x,y) * OKA2P(Genos(l,A), u, y) * SUM(AKA2P(u,y,:) * AHWE(:,l)) 
            else !if (Parent(A,3-k)==0) then
              PrXYU(x,y,u) = PrXY(x,y) * SUM(AKAP(u,:,l) * OKA2P(Genos(l,A), u, :) *  AHWE(:,l))
            endif
          enddo
        endif
        PrXY(x,y) = PrXY(x,y) * OKA2P(Genos(l,A), x, y)       
      enddo
    enddo
    PrL(l,2) = LOG10(SUM(PrXY))   ! SB offspring of Parent(A,3-k)
    PrLU(l,2) = LOG10(SUM(PrXYU)) ! A inbred via another parent 
  endif
enddo
LLtmp = SUM(PrL, dim=1)
LLU(1:3) = SUM(PrLU, dim=1)
call CalcU(A,k, -SB, k, LLU(4))   ! unrelated & A non-inbred
if (maybe(1) .and. Parent(A,3-k)>0) then
  LLtmp(1) = LLtmp(1) - Lind(Parent(A,3-k))
endif

LL = impossible
do x=1,3
  if (.not. maybe(x))  cycle
  if (LLU(x) > LLU(4) .or. x==2) then
    LL(x) = LLtmp(x) - LLU(x) + LLU(4)
  else
    LL(x) = LLtmp(x)
  endif
enddo

! if (A==1650 .and.  any(SibID(:,SB,k)==1651)) then
  ! print *, ""
  ! write (*, '("addSibInbr:", 3l4, 3f8.1)') maybe, SUM(PrL, dim=1)
  ! write (*, '("LLU:", 4f8.1, ", par: ", i5)')  LLU, Par
  ! print *, ""
! endif

end subroutine AddSibInbr

! #####################################################################

! subroutine MkSibInbr(A,k,kG,LL)   ! A result of half-sib sib mating  (TODO? FS mating)
! use Global
! implicit none

! integer, intent(IN) :: A,k, kG
! double precision, intent(OUT) :: LL
! integer :: l, x, y, z, GG
! double precision :: PrL(nSnp), PrY(3), PrZ(3), PrG(3), PrXYZ(3,3,3)

! LL = missing
! if (Parent(A,k)>=0 .or. Parent(A,3-k)/=0) then
  ! LL = NotImplemented
  ! return
! endif

! PrL = 0D0
! do l=1,nSnp
  ! call ParProb(l, GG, kG, 0, 0, PrZ)
  ! call ParProb(l, Parent(A,k),k, A, -4, PrY)
  ! call ParProb(l, GpID(3-kG, -Parent(A,k),k), 3-kG, 0, 0, PrG)  
  ! do x=1,3
    ! do y=1,3
      ! do z=1,3
        ! PrXYZ(x,y,z) = OKA2P(Genos(l,A),x,y) * PrY(y) * SUM(AKA2P(y,z,:) * PrG) * &
          ! AKAP(x,z,l) * PrZ(z)
      ! enddo
    ! enddo
  ! enddo
  ! PrL(l) = LOG10(SUM(PrXYZ))
! enddo
! LL = SUM(PrL)

! end subroutine MkSibInbr

! #####################################################################

subroutine MergeSibs(SA, SB, k, LL)  
use Global
use CalcLik
implicit none

integer, intent(IN) :: SA, SB, k
double precision, intent(OUT) :: LL
integer :: G(2), nAB(2), AB(2,maxSibsize), catG, catGG(2), GGP(2), ParPar(2), &
  catA(ns(SA,k)), catB(ns(SB,k)), l,x,y, r,v, Bj, Ai, i,  m, Ei, e,f, j, z
double precision :: PrL(nSnp, 2), PrG(3,2), PrXYZ(3,3,3,2), PrX(3,2), PrE(3), PrH(3), PrFS(3,3)!, PrXY(3,3)
logical :: AncOK(2), ParIsClone(2,maxSibsize)

LL = missing
G = 0  
do m=1,2  
  if (GpID(m,SA,k) /= 0) then
    if (GpID(m,SB,k) /= 0 .and. GpID(m,SA,k)/=GpID(m,SB,k)) then
      LL = impossible  ! incompatible grandparents
    else
      G(m) = GpID(m,SA,k)  ! including if GP(B) is dummy
    endif
  else if (GpID(m,SA,k) == 0) then
    G(m) = GpID(m,SB,k)
  endif
enddo
if (GpID(k, SA,k)==-SB .or. GpID(k, SB, k)==-SA) then
  LL = impossible
endif
if (LL==impossible) return

AncOK = .TRUE.              
call ChkAncest(-SA,k, -SB,k, AncOK(1))
call ChkAncest(-SB,k, -SA,k, AncOK(2))
if (any(.not. AncOK)) then
  LL = impossible
  return
endif

nAB(1) = ns(SA, k)
nAB(2) = ns(SB, k)
AB(1, 1:ns(SA,k)) = SibID(1:ns(SA,k), SA, k)
AB(2, 1:ns(SB,k)) = SibID(1:ns(SB,k), SB, k)

catG = 0
catGG = 0
GGP = 0
if (ANY(G/=0)) then
  do j=1,2 
    do i=1,nAB(j)
      if (nFS(AB(j,i))==0) cycle
      if (Parent(AB(j,i), 3-k)==0) cycle
      if (Parent(AB(j,i), 3-k) == G(3-k)) then
         if (catG==0) then
          catG = AB(j,i)
          exit  
        endif   
      endif
      ParPar = getPar(Parent(AB(j,i), 3-k), 3-k)
      do m=1,2
        if (ParPar(m) == G(m) .and. G(m)/=0) then
          catGG(m) = AB(j,i)
          GGP(3-m) = ParPar(3-m)
        endif
      enddo
    enddo
    if (catG /= 0) exit
  enddo
endif

catA = 0        
catB = 0
do r = 1, nS(SB, k)
  Bj = SibID(r, SB, k) 
  if (nFS(Bj)==0 .or. Parent(Bj,3-k)==0) cycle
  do v=1,nS(SA,k)
    Ai = SibID(v, SA, k)   
    if (nFS(Ai)==0) cycle
    if (Parent(Ai,3-k) == Parent(Bj,3-k) .and. (Parent(Bj,3-k) < 0 .or. catG==Ai)) then  
      catA(v) = r
      catB(r) = 1 
    endif
  enddo
enddo

! quick for simple cases
! if (catG==0 .and. all(catGG==0) .and. all(catA==0) .and. all(catB==0) .and. &
  ! all(G>0) .and. all(nAB == 1) .and. hermaphrodites==0) then
  ! PrL = 0D0
  ! Ai = SibID(1, SA, k) 
  ! Bj = SibID(1, SB, k) 
  ! do l=1,nSnp
    ! do m=1,2
      ! call ParProb(l, G(m), m, 0, 0, PrG(:,m))
    ! enddo
    ! call ParProb(l, Parent(Ai,3-k),3-k,Ai,0,PrE)
    ! call ParProb(l, Parent(Bj,3-k),3-k,Bj,0,PrH)
    ! do x=1,3
      ! PrXY(x,:) = SUM(OKA2P(Genos(l,Ai),x,:) * PrE) * SUM(OKA2P(Genos(l,Bj),x,:) * PrH)
      ! do y=1,3
        ! PrXY(x,y) = PrXY(x,y) * PrG(y,1) * SUM(AKA2P(x,y,:) * PrG(:,2)) 
      ! enddo
    ! enddo
    ! PrL(l,1) = LOG10(SUM(PrXY))
  ! enddo
  ! LL = SUM(PrL(:,1))
  ! return
! endif

ParIsClone = .FALSE.
if (hermaphrodites/=0) then
  if (DumClone(SA,k)/=0) then
    do i=1, nS(SA,k)
      if (Parent(SibID(i,SA,k), 3-k) == -DumClone(SA,k)) then
        ParisClone(1,i) = .TRUE.
      endif
    enddo
  endif
  if (DumClone(SB,k)/=0) then
    do i=1, nS(SB,k)
      if (Parent(SibID(i,SB,k), 3-k) == -DumClone(SB,k)) then
        ParisClone(2,i) = .TRUE.
      endif
    enddo
  endif
endif

PrL = 0D0
do l=1,nSnp
  do m=1,2
    if (m/=k .and. catG/=0) then
      call ParProb(l, G(m), m, -1, 0, PrG(:,m))
    else if (catGG(m)/=0) then  
      if (Parent(catGG(m),3-k) > 0) then
        call ParProb(l, G(m), m, Parent(catGG(m),3-k), 0, PrG(:,m))
      else
        call ParProb(l, G(m), m, 0, 0, PrG(:,m))
      endif
    else
      call ParProb(l, G(m), m, 0, 0, PrG(:,m))
    endif
  enddo
  do x=1,3
    do y=1,3
      do z=1,3
        PrXYZ(x,y,z,:) = AKA2P(x, y, z) * PrG(y,3-k) * PrG(z,k)
      enddo     
    enddo
    PrX(x,1) = SUM(PrXYZ(x,:,:,1))
    PrX(x,2) = PrX(x,1)
  enddo
  
  do z=1,3   ! GP k
    do y=1,3  ! GP 3-k
      if ((y>1 .or. z>1) .and. ALL(catGG==0) .and. catG==0) cycle
  do x=1,3
    do j=1,2
      do r=1, nAB(j)
       if (j==2) then
         if (catB(r)==1) cycle ! done as FS of an A
       endif
        Ai = AB(j,r)
        if (NFS(Ai) == 0) cycle
        if (catG==Ai .or. ParIsClone(j,r)) then
          PrE = 1D0
        else if (ANY(catGG == Ai)) then
          if (ALL(catGG == Ai)) then  ! parent(Ai, 3-k) FS with SB
            PrE = AKA2P(:,z,y)
          else      
            do m=1,2
              if (catGG(m)==Ai) then
                if (Parent(Ai,3-k)>0) then
                  call ParProb(l, GGP(3-m), 3-m, Parent(Ai,3-k),0,PrH)
                else
                  call ParProb(l, GGP(3-m), 3-m, 0,0,PrH)
                endif
                do e=1,3
                  if (m==k) then
                    PrE(e) = SUM(AKA2P(e, z, :) * PrH)
                  else
                    PrE(e) = SUM(AKA2P(e, y, :) * PrH)
                  endif
                enddo
              endif
            enddo
          endif
          if (Parent(Ai,3-k)>0) then
            PrE = PrE * OcA(:,Genos(l, Parent(Ai,3-k)))
          endif
!          PrE = PrE/SUM(PrE)
        else
          call ParProb(l, Parent(Ai,3-k), 3-k, -1, 0, PrE)
        endif

        if (Parent(Ai,3-k) < 0) then 
          do e=1,3
            if (catG==Ai .and. y/=e)  cycle
            do f=1, nS(-Parent(Ai,3-k), 3-k)
              Ei = SibID(f, -Parent(Ai,3-k), 3-k)
              if (nFS(Ei) == 0) cycle
              if (Parent(Ei, k)==-SB .or. Parent(Ei,k)==-SA) cycle  
              if (catGG(3-k)>0) then
                if (Parent(catGG(3-k),3-k) == Ei)  cycle
              endif  
              call ParProb(l, Parent(Ei, k), k, Ei, -1, PrH)
              PrFS = FSLik(l,Ei)
              PrH = PrH * PrFS(:,e)
              if (.not. ALL(PrH==1D0))  PrE(e) = PrE(e) * SUM(PrH)
            enddo
          enddo
        endif

        if (.not. ALL(PrE==1D0)) then
          if (catG==Ai) then
            PrXYZ(x,y,z,1) = PrXYZ(x,y,z,1) * PrE(y)
          else if (ANY(catGG/=0) .or. catG/=0) then
            PrXYZ(x,y,z,1) = PrXYZ(x,y,z,1) * SUM(PrE)
          else if (ParIsClone(j,r)) then
            PrX(x,1) = PrX(x,1) * PrE(x)
          else 
            PrX(x,1) = PrX(x,1) * SUM(PrE)
          endif
        endif
        
        PrFS = FSLik(l,Ai)
        PrE = PrE * PrFS(:,x)

        if (j==1) then
          if (catA(r)/=0) then
            Bj = AB(2, catA(r))
            PrFS = FSLik(l,Bj)
            PrE = PrE * PrFS(:,x)
          endif
        endif

        if (.not. ALL(PrE==1D0)) then
          if (catG==Ai) then
            PrXYZ(x,y,z,2) = PrXYZ(x,y,z,2) * PrE(y)
          else if (ANY(catGG/=0) .or. catG/=0) then
            PrXYZ(x,y,z,2) = PrXYZ(x,y,z,2) * SUM(PrE)
          else if (ParIsClone(j,r)) then
            PrX(x,2) = PrX(x,2) * PrE(x)
          else 
            PrX(x,2) = PrX(x,2) * SUM(PrE)
          endif 
        endif
      enddo  ! r
    enddo  ! j
  enddo  ! x
  enddo  ! y (catGG>0 only)
  enddo  ! z (catGG>0 only)
  do m=1,2
    if (ANY(catGG/=0) .or. CatG/=0) then
      PrL(l,m) = LOG10(SUM(PrXYZ(:,:,:,m)))
    else
      PrL(l,m) = LOG10(SUM(PrX(:,m)))
    endif
  enddo
enddo
LL = SUM(PrL(:,2)) - SUM(PrL(:,1))

! if (any(SibID(:,SA,k)==66) .and. ANY(SibID(:,SB,k)==63)) then
! open (unit=42,file="log.txt",status="unknown", position="append")
    ! write (42, *) ""
 ! write(42,'("MergeSibs: ", 3f9.1)') SUM(PrL(:,2)), SUM(PrL(:,1)), LL
 ! write (42, *) "catG: ", catG, catGG, " G: ", G
  ! write (42, *) "SA:"
  ! do i=1,ns(SA,k)
    ! Ai = SibID(i, SA, k) 
    ! write(42,*)  i, Ai, nFS(Ai), Parent(Ai,3-k), catA(i)
  ! enddo
  ! write(42,*) "SB:"
  ! do j=1, ns(SB,k)
    ! Bj = SibID(j,SB,k)
    ! write(42,*)  j, Bj, nFS(Bj), Parent(Bj,3-k), catB(j)
  ! enddo
  ! write(42,*) ""
  ! close(42)
! endif

end subroutine MergeSibs

! #####################################################################

subroutine AddFS(A, SB, kB, SA, kAx, LL, TopSib, dLL)  ! A/SA FS with any B?
use Global
use CalcLik
implicit none

integer, intent(IN) :: A, SB, kB, SA, kAx
integer, intent(OUT) :: TopSib    ! most likely FS of A within SB
double precision, intent(OUT) :: LL, dLL(maxSibSize)  ! dLL
integer :: PA, kA, AncA(2,mxA), GA(2), InbrX, Par(nS(SB,kB)), MaybeFS(nS(SB,kB)), &
  f, i, Bj, Inbr(nS(SB,kB)),  GB(ns(SB,kB), 2), DoQuick, l,x,y,g,h,Ei,z,j,Rj, ParAtmp                       
double precision :: ALR, LRQ, LLtmp(2), LLUX, PrL(nSnp, nS(SB,kB),2), PrXY(3,3,2), &
  PrG(3), PrX(3,2), PrY(3,2), PrW(3), PrZ(3), PrV(3), PrB(3,3), PrF(3,3)
logical :: ParOK, ParBisClone(ns(SB,kB))

LL = missing
TopSib = 0                                        
dLL = missing

if (nS(SB,kB)==0) then
  LL = impossible
  return   ! nobody to be FS with
endif

PA = 0
if (A /= 0) then
  if (kAx==1 .or. kAx==2) then
    kA = kAx
  else
    kA = 1
  endif
  PA = Parent(A, 3-kB)
  call GetAncest(A, kA, AncA)
else !if (SA /= 0) then   ! TODO: does it matter if kA=kB?
  kA = kAx
  PA = GpID(3-kB, SA, kA)
  call GetAncest(-SA, kA, AncA)
endif
GA = getPar(PA, 3-kB)          

if (A/=0) then
  if (Parent(A,kB)/=0 .and. Parent(A,kB)/=-SB) then
    LL = impossible
  else
    call ChkValidPar(A,Sex(A), -SB,kB, ParOK)
    if (.not. ParOK)  LL = impossible
  endif
else !if (SA/=0) then
  if (GpID(kB, SA, kA)/=0 .and. GpID(kB, SA, kA)/=-SB) then
    LL = impossible
  else
    call ChkValidPar(-SA,kA, -SB,kB, ParOK)
    if (.not. ParOK)  LL = impossible
  endif
endif
if (LL /= missing) return

InbrX = 0
if (ANY(AncA(kB, 3:mxA) == -SB) .and. A < 0) then  
  LL = NotImplemented   ! TODO: check  
else if (A > 0 .and. PA/=0) then
  if (AncA(kB,5-kB)==-SB .and. PA<0) then
    InbrX = -1  ! P-O mating
  else if (ANY(Parent(SibID(1:nS(SB,kB),SB,kB),3-kB)==GA(3-kB)) .and. GA(3-kB)/=0) then
    InbrX = -2
  else if (GpID(3-kB, SB, kB)==PA) then
    InbrX = -3
  else if (ANY(AncA(kB, 3:mxA) == -SB)) then
    if (.not. any(Parent(SibID(1:ns(SB,kB),SB,kB),3-kB)==PA)) then
      LL = NotImplemented
    endif
  endif
endif
if (LL /= missing) return

if (SA/=0 .and. GpID(3-kA, SB, kB)/=0) then
  do i=1, ns(SA,kA)
    if (Parent(SibID(i,SA,kA), 3-kA) == GpID(3-kA, SB, kB)) then
      LL = NotImplemented
      return
    endif
  enddo
endif

Par = 0  ! shared parent 3-kB  (cand. parent(kB) == SB)
MaybeFS = 1                                                                   
if (PA/=0 .and. ANY(Parent(SibID(1:nS(SB,kB),SB,kB), 3-kB) == PA)) then
  MaybeFS = 0
  do f=1, nS(SB,kB)
    if (NFS(SibID(f, SB, kB))==0) then
      MaybeFS(f) = -1
    else if (Parent(SibID(f, SB, kB), 3-kB) == PA) then
      MaybeFS(f) = 1
      Par(f) = PA
    else
      Par(f) = Parent(SibID(f, SB, kB), 3-kB)
    endif
  enddo
  
else    
 do f=1, nS(SB,kB)
  if (NFS(SibID(f, SB, kB))==0) then
    MaybeFS(f) = -1
    cycle
  endif   
  do i=1,nFS(SibID(f, SB, kB))
    Bj = FSID(i, SibID(f, SB, kB))
    if (A == Bj) then
      LL = AlreadyAss
    else if (A >0) then
      if (Parent(A,3-kB) == Bj) then  
         MaybeFS(f) = 0     ! can't be FS with own parent    
      else if (Parent(Bj, 3-kB) == A) then
        MaybeFS(f) = 0
      else 
        call ChkValidPar(A, Sex(A), Parent(Bj,3-kB), 3-kB, ParOK)
        if (.not. ParOK) then
          MaybeFS(f) = 0
        else
          call CalcAgeLR(A, Sex(A), Bj, Sex(Bj), kB, 2, .TRUE., ALR)
          if (ALR==impossible)  MaybeFS(f) = 0
        endif
      endif
      
    else if (SA/=0) then
      if (kA/=kB .and. Parent(Bj, 3-kB) == -SA) then
        MaybeFS(f) = 0  ! cannot be FS with own parent
        LL = NotImplemented   ! TODO: implement. 
        cycle
      else
        call ChkValidPar(-SA,kA, Parent(Bj,3-kB), 3-kB, ParOK)
        if (.not. ParOK) then
          MaybeFS(f) = 0
        else
          call CalcAgeLR(-SA,kA, Bj, Sex(Bj), kB, 2, .TRUE., ALR)
          if (ALR==impossible)  MaybeFS(f) = 0
        endif  
      endif
    endif
    if (Bj == PA .or. (A/=0 .and. A == Parent(Bj, 3-kB))) then
      MaybeFS(f) = 0
      cycle
    endif
    if (PA>0) then
      if (any(Parent(PA,:)==Bj)) then
        MaybeFS(f) = 0
        cycle
      endif
    endif
    
    Par(f) = Parent(Bj, 3-kB)
    if (PA/=0 .and. PA/=Par(f) .and. Par(f)/=0) then
      MaybeFS(f) = 0
    else if (Par(f)==0) then
      call CalcP2(Bj, Sex(Bj), -SB, PA, kB, LRQ)  
      if (LRQ == impossible) then
        MaybeFS(f) = 0
        cycle
      endif 
      Par(f) = PA
    endif
  enddo
 enddo
endif
if (LL /= missing) return
if (ALL(MaybeFS==0 .or. MaybeFS==-1)) then
  LL = impossible
  return
endif                         

Inbr = 0
GB = 0
do f=1, nS(SB,kB)    
  GB(f,:) = getPar(Par(f), 3-kB)
  if (GB(f,kB) == -SB)  Inbr(f) = 1 
  if (Par(f) == GpID(3-kB, SB, kB) .and. Par(f)/=0) then   ! DoQuick = -1
    Inbr(f) = 2
!    if (InbrX == -3)  InbrX = 0
  endif
enddo

do f=1, nS(SB,kB)
  if (nFS(SibID(f, SB, kB))==0) cycle
  if (MaybeFS(f)<1 .or. Par(f)==0 .or. Par(f)==PA) cycle
  if (A/=0)  call ChkValidPar(A, Sex(A), Par(f), 3-kB, ParOK)
  if (SA/=0) call ChkValidPar(-SA, kA, Par(f), 3-kB, ParOK)
  if (.not. ParOK)  MaybeFS(f) = 0  
enddo
if (ALL(MaybeFS==0 .or. MaybeFS==-1)) then
  LL = impossible
  return
endif

call ChkDoQuick(SB,kB,DoQuick)

ParBisClone = .FALSE.
if (hermaphrodites/=0 .and. DoQuick==-3) then
  do f=1, nS(SB,kB)
    if (DumClone(SB,kB) == -Parent(SibID(f,SB,kB), 3-kB)) then
      ParBisClone(f) = .TRUE.
    endif
  enddo
endif

if (A/=0 .and. nYears>1) then    ! check if A is more likely to be parent of sib instead
  do f=1, nS(SB,kB)
    Bj = SibID(f, SB, kB)
    if (MaybeFS(f)<1 .or. Par(f)/=0 .or. Parent(Bj, 3-kB)/=0) cycle  
    if (AgeDiff(Bj, A) <= 0) cycle
    if (Sex(A)<3 .and. Sex(A)/=3-kB)  cycle  ! TODO check both sexes?
    call ChkValidPar(Bj, 3, A, 3-kB, ParOK)
    if (.not. ParOK)  cycle
    call CalcU(-SB, kB, A, kB, LLtmp(1))
    call setParTmp(Bj, 3, A, 3-kB)
    call CalcU(-SB, kB, A, kB, LLtmp(2))
    call setParTmp(Bj, 3, 0, 3-kB)
    call CalcCLL(SB,kB)                   
    if (LLtmp(1) - LLtmp(2) < TA) then
      MaybeFS(f) = 0
    endif
  enddo
endif
if (ALL(MaybeFS==0 .or. MaybeFS==-1)) then
  LL = impossible
  return
endif

LLtmp = missing
if (A>0 .and. PA/=0 .and. any(Parent(SibID(1:ns(SB,kB),SB,kB),3-kB)==PA) .and. &  ! A already HS via 3-k 
   (InbrX==-2 .or. InbrX==-3 .or. DoQuick>1))  then
  call CalcU(-SB, kB, A, Sex(A), LLUX)
  if (Parent(A,3-kB) < 0) then
    call CalcU(-SB, kB, Parent(A,3-kB), 3-kB, LLtmp(1))
  else
    LLtmp(1) = LLUX
  endif
  ParAtmp = Parent(A,kB)   ! 0 or -SB
  call setParTmp(A, sex(A), -SB, kB)
  if (Parent(A,3-kB) < 0) then
    call CalcU(-SB, kB, Parent(A,3-kB), 3-kB, LLtmp(2))
  else
    LLtmp(2) = CLL(SB,kB)
  endif
  call setParTmp(A, sex(A), ParAtmp, kB)  
  call CalcCLL(SB,kB)
  LL = LLUX + LLtmp(2) - LLtmp(1)
   
  do f=1, ns(SB,kB)
    if (nFS(SibID(f, SB, kB))==0) cycle
    if (Parent(SibID(f, SB, kB), 3-kB) == Parent(A, 3-kB)) then
      TopSib = SibID(f, SB, kB)
      dLL(f) = LLtmp(2) - LLtmp(1)
    endif
  enddo
   
else 
  PrL = 0D0
  
  if (A>0 .and. DoQuick==-2) then  ! SB are all FS
    do f=1, nS(SB,kB)
      if (MaybeFS(f) < 1) cycle
      Bj = SibID(f, SB, kB)
      do l=1, nSnp
        PrB = FSLik(l,Bj)
        do x=1,3
          do y=1,3
            PrXY(x,y,:) = XPr(2,x,l,SB,kB) * XPr(2,y,l,-par(f),3-kB) * PrB(x,y)
            PrXY(x,y,1) = PrXY(x,y,1) * OKA2P(Genos(l,A),x,y)
            if (PA/=0) then
              PrXY(x,y,2) = PrXY(x,y,2) * OKAP(Genos(l,A), y, l)
            else
              PrXY(x,y,2) = PrXY(x,y,2) * OHWE(Genos(l,A), l)
            endif
          enddo
        enddo
        do i=1,2
          PrL(l,f,i) = LOG10(SUM(PrXY(:,:,i)))
        enddo
      enddo
    enddo
  
  else
    do l=1,nSnp
      if (any(inbr == 2)) then
        call ParProb(l, GpID(kB,SB,kB), kB, -1,0, PrG)
      endif 
      do f=1, nS(SB,kB)
        if (MaybeFS(f) < 1) cycle
        if (any(inbr == 2)) then
          PrX = 1D0
        else
          PrX(:,1) = XPr(2,:,l, SB, kB)
          PrX(:,2) = PrX(:,1)
        endif
        do x=1,3   ! SB
          do g=1,nS(SB,kB)
            Bj = SibID(g, SB, kB)
            if (NFS(Bj) == 0) cycle
            if (Inbr(g)==1) then
              if (g==f) then
                PrY(:,1) = 1D0
              else
                call ParProb(l, Parent(Bj, 3-kB), 3-kB, Bj,-5, PrY(:,1))
                  ! no GPs & no Bj & no FS of Bj 
              endif
              call ParProb(l, GB(g,3-kB), 3-kB, 0,0,PrW)
              do y=1,3
                PrY(y,1) = PrY(y,1) * SUM(AKA2P(y,:,x) * PrW)
              enddo
            else if (ParBisClone(g)) then
              PrY(:,1) = 1D0  
            else
              call ParProb(l, Par(g), 3-kB, -1,0, PrY(:,1))
              if (Inbr(g)==2) then
                do y=1,3
                  PrY(y,1) = PrY(y,1) * SUM(AKA2P(x,y,:) * PrG)
                enddo
              endif
            endif
            
            PrY(:,2) = PrY(:,1)  ! 1: FS, 2: HS via 3-k or U 
            do y=1,3
              if (ParBisClone(g) .and. y/=x)  cycle
              if (Par(g) < 0) then 
                do h = 1, nS(-Par(g), 3-kB)
                  Ei = SibID(h, -Par(g), 3-kB)
                  if (Parent(Ei, kB) == -SB .or. Ei==A) cycle  
                  if (NFS(Ei) == 0) cycle 
                  if (g==f .and. Parent(Ei,kB)<0) then  
                    call ParProb(l, Parent(Ei,kB), kB, -1,0, PrZ)
                    do z=1,3
                      do j=1, ns(-parent(Ei,kB),kB)
                        Rj = SibID(j,-parent(Ei,kB),kB)
                        if (Parent(Rj,3-kB) == Par(g) .and. Par(g)/=0)  cycle
                        if (nFS(Rj)==0)  cycle
                        call ParProb(l, Parent(Rj,3-kB),3-kB, Rj,-1, PrV)
                        PrF = FSLik(l,Rj)
                        PrV = PrV * PrF(:,z)
                        if (.not. all(PrV==1D0))  PrZ(z) = PrZ(z) * SUM(PrV)
                      enddo
                    enddo
                  else
                    call ParProb(l, Parent(Ei,kB), kB, Ei,-1, PrZ)
                  endif
                  PrF = FSLik(l,Ei)
                  PrZ = PrZ * PrF(:,y)
                  if (.not. all(PrZ==1D0))  PrY(y,:) = PrY(y,:) * SUM(PrZ)
                enddo
              endif
              
              PrB = FSLik(l,Bj)
              PrY(y,:) = PrY(y,:) * PrB(x,y)
                
              if (g==f) then  
                if (A/=0) then
                  PrY(y,1) = PrY(y,1) * OKA2P(Genos(l,A), x, y)
                  PrY(y,2) = PrY(y,2) * OHWE(Genos(l,A), l)   ! vs unrelated !
                else if (SA/=0) then
                  PrY(y,1) = PrY(y,1) * SUM(XPr(1,:,l, SA,kA) *AKA2P(:,x,y))
                  if (PA/=0) then
                    PrY(y,2) =PrY(y,2) *SUM(XPr(1,:,l, SA,kA) *AKAP(:,y, l))
                  else
                    PrY(y,2) = PrY(y,2) * SUM(XPr(1,:,l, SA,kA) * AHWE(:,l))
                  endif
                endif 
              endif          
            enddo  ! y
            do i=1,2
              if (ParBisClone(g)) then
                PrX(x,i) = PrX(x,i) * PrY(x,i)
              else if (.not. all(PrY(:,i)==1D0)) then
                PrX(x,i) = PrX(x,i) * SUM(PrY(:,i))
              endif
            enddo
          enddo  ! g
        enddo  ! x
        PrL(l,f,:) = LOG10(SUM(PrX, DIM=1))
      enddo  ! f   
    enddo  ! l
  endif

  dLL = impossible
  do f = 1, nS(SB, kB)
    if (MaybeFS(f)<1) cycle
    Bj = SibID(f, SB, kB)
    if (NFS(Bj) == 0) then
      cycle
    else if (nFS(Bj) == 1) then
      dLL(f) = SUM(PrL(:, f,1)) - SUM(PrL(:,f,2))
    else
      do g=1, nS(SB,kB)
        if (Parent(SibID(g, SB, kB), 3-kB) == Parent(Bj, 3-kB)) then
          dLL(g) = SUM(PrL(:, f,1)) - SUM(PrL(:,f,2))
        endif
      enddo
    endif
  enddo

  if (A/=0) then
    ! calc LLUX with parA=0, else LL(FS)|PA/=0 /= LL(FS)|PA==0
    if(PA/=0)  call setParTmp(A,kA,0,3-kB)
    call CalcU(A,kA, -SB, kB, LLUX)
    if(PA/=0)  call setParTmp(A,kA,PA,3-kB)  
    LL = MAXVAL(dLL, MASK=dLL/=impossible) + LLUX
    TopSib = MAXLOC(dLL, MASK=dLL/=impossible, DIM=1)
    if(TopSib>0)  TopSib = SibID(TopSib, SB, kB)
  else if (SA/=0) then
    call CalcU(-SA, kA, -SB, kB, LLUX)
    do f = 1, nS(SB, kB)
      if (dLL(f)==impossible) cycle
      if (Par(f)==0 .and. nS(SA,kA)>1) then
        dLL(f) = dLL(f) + LLUX
      else  ! consider changes in SA (e.g. inbreeding loops) 
        call setParTmp(-SA,kA, Par(f), 3-kB)
        call PairUA(-SA, -SB, kA, kB, dLL(f))  
      endif
    enddo
    call setParTmp(-SA,kA, PA, 3-kB)
    call CalcCLL(SA,kA)
    LL = MaxLL(dLL) 
  endif  
endif


! if (A==77 .and. any(SibID(:,SB,kB)==59)) then
! !if (SB==161 .and. kB==1 .and. SA==100 .and. kA==2) then
  ! open (unit=42,file="log.txt",status="unknown", position="append")
  ! write (42, *) ""
  ! write(42, '("Add FS", 2i5, "--", 3i5, "; ", 2i5, " - ", i4, 2f8.1, " ", 2i3)') & 
     ! kB, SB, kA, A, SA, GpID(:,SB,kB), PA, LLUX, LL, InbrX, DoQuick
  ! if (A>0)  write(42, '("nFS A: ", 2i5, ", LLtmp: ", 2f8.1)') nFS(A), FSID(maxSibSize+1, A), LLtmp
    ! do f=1, nS(SB,kB)
      ! Bj = SibID(f, SB, kB)
        ! write(42,'(a8, 3i6, 3f8.1, 2i3, i6)') ID(Bj), Parent(Bj, 3-kB), &
          ! Par(f), MaybeFS(f), dLL(f), SUM(PrL(:, f, 1)), SUM(PrL(:, f, 2)), Inbr(f), &
          ! nFS(Bj), FSID(maxSibSize+1, Bj)
    ! enddo
    ! write (42, *) "TopSib: ", TopSib!, dLLOUT
  ! write (42, *) ""
  ! close(42)
! endif

end subroutine AddFS

! #####################################################################

subroutine AddParent(A, SB, k, LL)  ! is A parent of sibship SB?  (replace dummy)
use Global
use CalcLik
implicit none

integer, intent(IN) :: A, SB, k
double precision, intent(OUT) :: LL
integer :: G(2), m, n, Inbr, DoQuick, l,x,y,Bi, Ei, i,j,z
double precision :: LRQ, PrL(nSnp), PrG(3,2), PrXY(3,3,2), PrZ(3), PrE(3), PrB(3,3)
logical :: AncOK, ParBisClone(ns(SB,k)), Aselfed

LL = missing
do i=1,ns(SB,k)
  if (AgeDiff(SibID(i,SB,k), A) <= 0) then
    LL = impossible
    return
  endif
enddo

call ChkAncest(A,0, -SB,k, AncOK) ! e.g. if age A unknown, or all age B's unknown
if (.not. AncOK) then
  LL = impossible
  return
endif

G = 0
do m=1,2
  if (Parent(A,m)/= 0) then   
    if (GpID(m,SB,k)/= 0 .and. GpID(m,SB,k) /= Parent(A,m)) then
      LL = impossible
      return
    else
      G(m) = Parent(A,m)
    endif
  else if(GpID(m,SB,k)/=0) then
    G(m) = GpID(m,SB,k)
  endif
enddo

LRQ = missing
do n=1, ns(SB,k)
  call CalcP2(SibID(n,SB,k), 3, A, Parent(SibID(n,SB,k), 3-k), k, LRQ)
  if (LRQ == impossible) then
    LL = impossible
    return
  endif
enddo
call CalcP2(A, k, GpID(1,SB,k), GpID(2,SB,k), 1, LRQ)
if (LRQ == impossible) then
  LL = impossible
  return
endif

Inbr = 0
if (G(3-k)/=0) then
  do n=1, nS(SB, k)
    if (nFS(SibID(n,SB,k))==0) cycle
    if (Parent(SibID(n,SB,k), 3-k)==G(3-k)) then
      Inbr = n
    endif
  enddo 
endif

call ChkDoQuick(SB,k,DoQuick)

ParBisClone = .FALSE.
Aselfed = .FALSE.
if (hermaphrodites/=0) then
  if (DoQuick==-3) then
    do i=1, nS(SB,k)
      if (DumClone(SB,k) == -Parent(SibID(i,SB,k), 3-k)) then
        ParBisClone(i) = .TRUE.
      endif
    enddo
  endif
  if (all(G > 0) .and. G(1) == G(2)) then
    Aselfed = .TRUE.
  else if (all(G < 0)) then
    if (DumClone(-G(1),1) == -G(2))  ASelfed = .TRUE.
  endif
endif

PrL = 0D0
do l=1,nSnp
  if (Inbr==0 .and. DoQuick>0) then
    do m=1,2
      call ParProb(l, G(m), m, A,0, PrG(:,m))    
    enddo
    do x=1,3
      do y=1,3
        PrXY(x,y,2) = XPr(1,x,l, SB,k) * &
          SUM(AKA2P(x, y, :) *  PrG(y,3-k) * PrG(:, k))
      enddo
    enddo
    PrL(l) = LOG10(SUM(PrXY(:,:,2)))
  else
    call ParProb(l, G(k), k, A,0, PrG(:,k))
    call ParProb(l, G(3-k), 3-k, -1,0, PrG(:,3-k))
    do x=1,3  ! A
      do y=1,3  ! G 3-k
        if (Aselfed) then
          PrXY(x,y,:) = OcA(x,Genos(l,A)) * AKA2P(x, y, y) *PrG(y, k)
        else
          PrXY(x,y,:) = OcA(x,Genos(l,A)) * SUM(AKA2P(x, y, :) *PrG(y,3-k) *PrG(:, k))
        endif
        do n=1, nS(SB,k)
          Bi = SibID(n, SB, k)
          if (nFS(Bi)==0) cycle
          if (Inbr == n .or. ParBisClone(n)) then
            PrZ = 1D0
          else
            call ParProb(l, Parent(Bi,3-k), 3-k, -1,0, PrZ)
          endif
          
          if (Parent(Bi,3-k) < 0) then
            do z=1,3
              if (Inbr == n .and. z/=y) cycle
              if (ParBisClone(n) .and. z/=x)  cycle 
              do i=1, ns(-Parent(Bi,3-k),3-k)
                Ei = SibID(i,-Parent(Bi,3-k),3-k)
                if (Parent(Ei,k)==-SB .or. nFS(Ei)==0)  cycle
                call ParProb(l, Parent(Ei,k), k, Ei, -1, PrE)
                do j=1, nFS(Ei)
                  if (FSID(j,Ei)==A)  cycle
                  PrE = PrE * OKA2P(Genos(l,FSID(j,Ei)),z,:)
                enddo
                if (.not. ALL(PrE==1D0))  PrZ(z) = PrZ(z) * SUM(PrE)  
              enddo
            enddo
          endif
          
          PrB = FSLik(l,Bi)
          if (Inbr == n) then
            PrXY(x,y,:) = PrXY(x,y,:) * PrZ(y)
            PrXY(x,y,2) = PrXY(x,y,2) * PrB(x,y)
          else if (ParBisClone(n)) then
            PrXY(x,y,:) = PrXY(x,y,:) * PrZ(x)
            PrXY(x,y,2) = PrXY(x,y,2) * PrB(x,x)
          else
            if (.not. all(PrZ==1D0)) PrXY(x,y,1) = PrXY(x,y,1) * SUM(PrZ)
            PrZ = PrZ * PrB(:,x)
            PrXY(x,y,2) = PrXY(x,y,2) * SUM(PrZ)
          endif
        enddo
      enddo
    enddo
    PrL(l) = LOG10(SUM(PrXY(:,:,2))) - LOG10(SUM(PrXY(:,:,1)))
  endif
enddo
LL = SUM(PrL) + Lind(A)
if (LL < -HUGE(0D0))  LL = impossible

end subroutine AddParent

! #####################################################################

subroutine AddGP(A, SB, k, LL)  ! add A as a grandparent to sibship SB
use Global
use CalcLik
implicit none 

integer, intent(IN) :: A, SB, k
double precision, intent(OUT) :: LL
integer :: cat, catG, DoQuick, curGP(2), l, x,y, m, i,  z, g, Bi,  Ei,w,v,j
double precision :: PrL(nSnp), LLtmp(2), LLU, PrG(3), PrPA(3,2), &
  PrXYZ(3,3,3,2), PrP(3), PrW(3), PrE(3), PrF(3,3), PrA(3)
logical :: ParOK, ParBisClone(ns(SB,k))

LL = missing
if (Sex(A)<3) then
  m = Sex(A)
else if (GpID(1,SB,k)==0) then
  m = 1
else if (GpID(2,SB,k)==0) then
  m = 2
else
  LL = impossible
  return   
endif

call ChkValidPar(-SB,k, A, m, ParOK)
if (.not. ParOK)  then
  LL = impossible
  return
endif

cat = 0
catG = 0
if (Parent(A, 3-k)==GpID(3-k,SB,k) .and. Parent(A, 3-k) /= 0) then
  cat = 1
else
  do i=1,nS(SB,k)
    if (nFS(SibID(i,SB,k))==0) cycle
    if (Parent(SibID(i, SB, k), 3-k) == 0) cycle
    if (Parent(SibID(i, SB, k), 3-k) == A) then
      cat = 4
    else if (Parent(SibID(i, SB, k), 3-k) == Parent(A, 3-k)) then
      catG = i
    else if (Parent(SibID(i, SB, k), 3-k) == GpID(3-k,SB,k)) then
      cat = 2
      exit
    else if (Parent(SibID(i,SB,k),3-k) < 0) then
      if (GpID(k, -Parent(SibID(i,SB,k),3-k), 3-k) == -SB) then
        cat = 3
        exit
      endif
    endif     
  enddo
endif

if (Complx < 2 .and. (cat/=0 .or. catG/=0)) then
  LL = NotImplemented
  return
endif

call ChkDoQuick(SB,k,DoQuick)   

if (DoQuick==2) then  ! inbreeding: Parent(Bj,3-k) = Offspr(i); offspr(i) < 0
  Bi = 0
  do i=1, nS(SB,k)
    if (Parent(SibID(i,SB,k), 3-k) < 0) then
      Bi = Parent(SibID(i,SB,k), 3-k)
      if (GpID(k, -Bi, 3-k) == -SB)  exit
    endif
  enddo
endif 

ParBisClone = .FALSE.
if (hermaphrodites/=0 .and. DoQuick==-3) then
  do i=1, nS(SB,k)
    if (DumClone(SB,k) == -Parent(SibID(i,SB,k), 3-k)) then
      ParBisClone(i) = .TRUE.
    endif
  enddo
endif

curGP = GPID(:, SB,k)
if (Complx==0 .and. Mate(A)/=0) then  
  if ((Mate(A) /= curGP(3-m) .and. curGP(3-m)/=0) .or. &
  (k==3-m .and. Mate(A)==-SB)) then          
    LL = impossible
    return
  else if (curGP(3-m) == 0) then
    curGP(3-m) = Mate(A)
  endif
endif

PrL = 0D0
LLU = missing
LLtmp = missing

if (cat/=0 .or. DoQuick>1 .or. (Parent(A,3-k)<0 .and. catG/=0)) then  
   ! inbreeding loop present / will be created:   
  call CalcU(-SB, k, A, 3-k, LLU)
  if (curGP(3-m) < 0) then
    call CalcU(-SB, k, curGP(3-m), 3-m, LLtmp(1))
  else if (Parent(A,3-k)<0) then
    call CalcU(-SB, k, Parent(A,3-k), 3-k, LLtmp(1))
  else if (DoQuick==2) then
    call CalcU(-SB, k, Bi, 3-k, LLtmp(1))  
  endif
  call setParTmp(-SB, k, A, m)
  if (curGP(3-m) < 0) then
    call CalcU(-SB, k, curGP(3-m), 3-m, LLtmp(2))
    LL = LLU + (LLtmp(2) - LLtmp(1))
  else if (Parent(A,3-k)<0) then
    call CalcU(-SB, k, Parent(A,3-k), 3-k, LLtmp(2))
    LL = LLU + (LLtmp(2) - LLtmp(1))
  else if (DoQuick==2) then
    call CalcU(-SB, k, Bi, 3-k, LLtmp(2)) 
    LL = LLU + (LLtmp(2) - LLtmp(1))    
  else
    LL = CLL(SB,k) + Lind(A)
  endif
  call setParTmp(-SB, k, curGP(m), m)
  if (curGP(3-m) < 0)  call CalcCLL(-curGP(3-m), 3-m)
  if (Parent(A,3-k) < 0)  call CalcCLL(-Parent(A,3-k), 3-k)
  
else
  do l=1,nSnp
    call ParProb(l, curGP(3-m), 3-m, 0, 0, PrG)
    if (catG/=0) then  
      call ParProb(l, Parent(A,3-k), 3-k, -1, 0, PrPA(:,3-k))
      call ParProb(l, Parent(A,k), k, -1, 0, PrPA(:,k))
    else
      PrA = LindX(:,l,A)/SUM(LindX(:,l,A))
    endif
    
    PrXYZ = 0D0
    do x=1,3  ! SB
      do y=1,3
        if (catG == 0) then
          PrXYZ(x,y,1,:) = SUM(AKA2P(x, :, y) * PrG * PrA(y))
        else
          do z=1,3  ! Parent(A,3-k) 
            do g=1,3
              PrP(g) = AKA2P(x,g,y) * PrG(g) * OcA(y,Genos(l,A)) * &
               SUM(AKA2P(y,z,:) * PrPA(:,k) * PrPA(z,3-k))  ! approx
            enddo
            PrXYZ(x,y,z,:) = SUM(PrP)
          enddo
        endif
      enddo
      if (DoQuick > 0 .and. catG==0) then
        PrXYZ(x,:,1,2) = PrXYZ(x,:,1,2) * XPr(1,x,l, SB,k)
      else
        do y=1,3  ! A
          do z=1,3  ! Parent(A,3-k)   
            if (z>1 .and. catG==0)  cycle
            do i=1, nS(SB,k)
              Bi = SibID(i, SB, k)
              if (nFS(Bi)==0) cycle
              if (catG == i .or. ParBisClone(i)) then
                PrW = 1D0
              else
                call ParProb(l, Parent(Bi,3-k), 3-k, -1, 0, PrW)
              endif
              if (Parent(Bi,3-k)<0) then
                do w=1,3
                  if (catG==i .and. w/=z)  cycle
                  if (ParBisClone(i) .and. w/=x)  cycle
                  do v=1, nS(-Parent(Bi, 3-k), 3-k)
                    Ei = SibID(v, -Parent(Bi, 3-k), 3-k)  
                    if (NFS(Ei) == 0) cycle
                    if (Parent(Ei, k) == -SB) cycle
                    call ParProb(l, Parent(Ei, k), k, Ei,-1, PrE)
                    if (catG==i) then 
                      do j=1, nFS(Ei)
                        if (FSID(j,Ei)==A)  cycle
                        PrE = PrE * OKA2P(Genos(l,FSID(j,Ei)),w,:)
                      enddo
                    else
                      PrF = FSLik(l,Ei)
                      PrE = PrE * PrF(:,w)
                    endif
                    if (.not. ALL(PrE==1D0))  PrW(w) = PrW(w) * SUM(PrE)  
                  enddo  
                enddo
              endif
              
              if (catG == i) then
                PrXYZ(x,y,z,1) = PrXYZ(x,y,z,1) * PrW(z)
              else if (catG /= 0) then
                if (.not. all(PrW==1D0))  PrXYZ(x,y,z,1) = PrXYZ(x,y,z,1) * SUM(PrW)
              else if (ParBisClone(i)) then
                PrXYZ(x,y,:,1) = PrXYZ(x,y,:,1) * PrW(x)
              else if (.not. all(PrW==1D0)) then
                PrXYZ(x,y,:,1) = PrXYZ(x,y,:,1) * SUM(PrW)
              endif 
              
              PrF = FSLik(l,Bi)
              PrW = PrW * PrF(:,x)
              
              if (catG == i) then
                PrXYZ(x,y,z,2) = PrXYZ(x,y,z,2) * PrW(z)
              else if (catG /= 0) then
                if (.not. all(PrW==1D0))  PrXYZ(x,y,z,2) = PrXYZ(x,y,z,2) * SUM(PrW)
              else if (ParBisClone(i)) then
                PrXYZ(x,y,:,2) = PrXYZ(x,y,:,2) * PrW(x)
              else if (.not. all(PrW==1D0)) then
                PrXYZ(x,y,:,2) = PrXYZ(x,y,:,2) * SUM(PrW)
              endif
            enddo 
          enddo
        enddo
      endif
    enddo
    PrL(l) = LOG10(SUM(PrXYZ(:,:,:,2))) - LOG10(SUM(PrXYZ(:,:,:,1)))
  enddo
  
  LL = SUM(PrL) + Lind(A)
endif

! if (A==156 .and. any(SibID(:,SB,k)==215))  then
  ! write(*,'("AddGP: ", 3i6, " + ", 4i6, "; cats: ", 3i4, f9.2)') &
    ! A, Parent(A,:), k, SB, GpID(:,SB,k), cat, catG, DoQuick, LL
! endif

end subroutine AddGP

! #####################################################################

subroutine AddGGP(A, SB, k, LL)
use Global
use CalcLik
implicit none
! A a GGP of sibship SB? (only calculating over non-gp-assigned)

integer, intent(IN) :: A, SB, k
double precision, intent(OUT) :: LL
integer :: m, n, catG, GG, AncG(2,mxA), DoQuick, l, x, y,z,i, v, Bj, w,Ei
double precision :: PrL(nSnp), PrXYZ(3,3,3,2), PrZ(3),PrA(3),PrP(3),PrV(3), &
  PrW(3), PrE(3), PrG(3), PrF(3,3)

LL = missing
if (GpID(1, SB,k)/=0) then
  if (GpID(2, SB,k)/=0) then  ! should be assigned as parent-of-gp
    LL = impossible   !(or AlreadyAss)
    return      
  else
    m = 2
  endif
else
  m = 1  ! doesn't really matter (?); GpID(m, SB, k) == 0.
endif

if (Sex(A)<3) then
  n = Sex(A)
else
  n = 1
endif

catG =0
GG = GpID(3-m, SB, k)
AncG = 0
if (GG/=0) then
  if (GG==Parent(A,3-m)) then
    catG = 1
  endif
  call GetAncest(GG, 3-m, AncG)  
  if (ANY(AncG == A)) then
    if ((GG>0 .and. AncG(n,2)==A) .or. (GG<0 .and. &
     AncG(n,3-m+2)==A)) then  
      catG = 2  ! already GGP via 3-m; check if double ggp
    else
      LL = NotImplemented  ! possible; not yet implemented
    endif
  else
    do v=1,2
      if (Parent(A,v)==0) cycle
      if ((GG<0 .and. ANY(AncG(v, 2:4)==Parent(A,v))) .or. &
       (GG>0 .and. ANY(AncG(v, 1:2)==Parent(A,v)))) then
        LL = NotImplemented   ! TODO: stricter implementation?
      endif
    enddo
  endif
endif
if (GpID(3-k,SB,k) < 0) then
  do i=1,nS(SB,k)
    if (Parent(SibID(i, SB, k), 3-k) == GpID(3-k,SB,k)) then 
      LL = NotImplemented
      exit
    endif
  enddo
endif
if (Parent(A,3-k)<0) then
  do i=1,nS(SB,k)
    if (Parent(SibID(i, SB, k), 3-k) == Parent(A,3-k)) then 
      LL = NotImplemented
      exit
    endif
  enddo
endif
if (LL==NotImplemented) return

if (catG/=0) then  ! age check
  if (GG>0) then
    if (AgeDiff(GG, A) >= 0 .and. catG==1 .and. AgeDiff(GG, A)/=missing) &
     LL = impossible  ! A older than GG
    if (AgeDiff(GG, A) <= 0 .and. catG==2)  LL = impossible  ! GG older than A
  else if (GG<0) then
    do v=1, nS(-GG, 3-m)  ! TODO? Age, ancestors
      if (AgeDiff(SibID(v,-GG,3-m), A) >= 0 .and. catG==1 .and. &
        AgeDiff(SibID(v,-GG,3-m), A)/=missing)  LL = impossible
      if (AgeDiff(SibID(v,-GG,3-m), A) <= 0 .and. catG==2)  LL = impossible
    enddo
  endif
endif
if (LL == impossible) return

if (Complx==0 .and. Mate(A)/=0) then  
  if (Mate(A) == GPID(3-n, SB,k)) then
    LL = impossible
    return
  else
    catG = 3
  endif
endif

call ChkDoQuick(SB,k,DoQuick)
if (DoQuick == -1)  then
  LL = NotImplemented  ! inbreeding loops, approx. below invalid
  return
endif

PrL = 0D0
do l=1,nSnp
  call ParProb(l, A, 0, 0, 0, PrA)
  if (catG==1) then
    call ParProb(l, GG, 3-m, A, 0, PrZ)
    call ParProb(l, Parent(A,m), m, -1, 0, PrP)
  else if (catG==2) then
    call ParProb(l, GG, 3-m, -4, 0, PrZ)
    if (GG > 0) then
      call ParProb(l, Parent(GG,3-n), 3-n, GG, 0, PrP) 
    else if (GG < 0) then
      call ParProb(l, GpID(3-n, -GG,3-m), 3-n, 0, 0, PrP)
    else
      PrP = AHWE(:,l)
    endif 
  else
    call ParProb(l, GG, 3-m, 0, 0, PrZ)
  endif
  if (catG==3) then
    call ParProb(l, Mate(A), 3-n, 0, 0, PrG)
  endif
  do x=1,3  ! SB
    do y=1,3  ! in between SB and A
      do z=1,3  ! other GP 
        do v=1,3  ! A
          if (catG==1) then
            PrV(v) = AKAP(y, v, l) * SUM(AKA2P(v,z,:) * PrP)
          else if (catG==2) then
            PrV(v) = AKAP(y, v, l) * SUM(AKA2P(z,v,:) * PrP)
          else if (catG==3) then
            PrV(v) = SUM(AKA2P(y,v,:) * PrG)
          else
            PrV(v) = AKAP(y, v, l)
          endif
        enddo
        PrXYZ(x,y,z,:) = AKA2P(x, y, z) * PrZ(z) * SUM(PrV * PrA)
      enddo
    enddo
    if (DoQuick > 0) then
      PrXYZ(x,:,:,2) = PrXYZ(x,:,:,2) * XPr(1,x,l, SB,k)
    else
      do i=1, nS(SB,k)
        Bj = SibID(i, SB, k)
        if (nFS(Bj)==0)  cycle
        call ParProb(l, Parent(Bj,3-k), 3-k, -1, 0, PrW)
        if (Parent(Bj,3-k)<0) then
          do w=1,3
            do v=1, nS(-Parent(Bj, 3-k), 3-k)
              Ei = SibID(v, -Parent(Bj, 3-k), 3-k)  
              if (NFS(Ei) == 0) cycle
              if (Parent(Ei, k) == -SB) cycle
              call ParProb(l, Parent(Ei, k), k, Ei,-1, PrE)
              PrF = FSLik(l,Ei)
              PrE = PrE * PrF(:,w)
              if (.not. ALL(PrE==1D0))  PrW(w) = PrW(w) * SUM(PrE)  
            enddo  
          enddo
          if (.not. all(PrW==1D0))  PrXYZ(x,:,:,1) = PrXYZ(x,:,:,1) * SUM(PrW)
        endif          
        PrF = FSLik(l,Bj)
        PrW = PrW * PrF(:,x)
        if (.not. all(PrW==1D0))  PrXYZ(x,:,:,2) = PrXYZ(x,:,:,2) * SUM(PrW)
      enddo 
    endif
  enddo
  PrL(l) = LOG10(SUM(PrXYZ(:,:,:,2))) - LOG10(SUM(PrXYZ(:,:,:,1)))          
enddo
if (catG==1) then
  LL = SUM(PrL)
!else if (catG==2 .and. GG>0 .and. DoQuick>0) then  ! double ggp
!  LL = SUM(PrL) + Lind(A) - Lind(GG)   ! Lind(GG) included in PrXYZ(..,1)   
else
  LL = SUM(PrL) + Lind(A)
endif

! if (A==502 .and. any(SibID(:,SB,k)==3485) )  then
  ! print *, "AddGP: ", A, Parent(A,:), " + ", SB, k, "; cats:", catG, DoQuick, GG
  ! write (*, '("LL: ", 2f9.2, "; ", 2f9.2)') SUM(PrL), LL, Lind(A)!, Lind(GG)
! endif

end subroutine AddGGP

! #####################################################################

subroutine addGAU(A, SB, k, m, LL)  ! A great-full-avuncular of B's
use Global
use CalcLik
implicit none

integer, intent(IN) :: A, SB, k, m
double precision, intent(OUT) :: LL                                                                 
integer :: nOff, Offspr(maxSibSize), sxOff(maxSibSize), DoQuick, &
  l, x, y, z, v,g, i, Bi, Ei, w,j, AncG(2,mxA)
double precision :: PrL(nSnp), PrGGG(3,2), PrXV(3,3,3,3,2), PrG(3), PrW(3), PrE(3), PrF(3,3)

LL = missing
if (GpID(m,SB,k)/=0) then
  LL = NotImplemented
else if (Parent(A,3-m) == GpID(3-m,SB,k) .and. Parent(A,3-m)/=0) then
  LL = AlreadyAss
else 
  call getOff(-SB,k, .TRUE., nOff, Offspr, sxOff)
  do g=1,2
    if (ANY(Offspr(1:nOff)==Parent(A,g))) then
      LL = impossible  
    else if (g/=k .and. Parent(A,g)/=0 .and. & 
     ANY(Parent(SibID(1:ns(SB,k),SB,k),g)==Parent(A,g))) then
      LL = NotImplemented
    endif
  enddo
endif
if (LL /= missing) return

AncG = 0
call getAncest(GpID(3-m,SB,k),3-m, AncG)
if (ANY(AncG == A)) then
  LL = NotImplemented
  return
endif

call ChkDoQuick(SB,k,DoQuick)
if (DoQuick == -1)  then
  LL = NotImplemented  ! inbreeding loops, approx. below invalid
  return
endif

PrL = 0D0
do l=1,nSnp
  do g=1,2
    call ParProb(l, Parent(A,g), g, A, 0, PrGGG(:,g))
  enddo
  call ParProb(l, GpID(3-m,SB,k),3-m,0,0,PrG)
  do x=1,3  ! sibship parent
    do y=1,3
      do z=1,3
        do v=1,3
          PrXV(x,y,z,v,:) = SUM(AKA2P(x,y,:) *PrG) *AKA2P(y,z,v) *PrGGG(z,1) *PrGGG(v,2)
          PrXV(x,y,z,v,2) = PrXV(x,y,z,v,2) * OKA2P(Genos(l,A), z, v)
          if (DoQuick > 0) then
            PrXV(x,y,z,v,2) = PrXV(x,y,z,v,2) * XPr(1,x,l, SB,k)
          else
            do i=1, nS(SB,k)
              Bi = SibID(i, SB, k)
              if (nFS(Bi)==0)  cycle
              call ParProb(l, Parent(Bi,3-k), 3-k,-1,0, PrW)
              if (Parent(Bi,3-k)<0) then
                do w=1,3
                  do j=1, nS(-Parent(Bi, 3-k), 3-k)
                    Ei = SibID(j, -Parent(Bi, 3-k), 3-k)  
                    if (NFS(Ei) == 0) cycle
                    if (Parent(Ei, k) == -SB) cycle
                    call ParProb(l, Parent(Ei, k), k, Ei,-1, PrE)
                    PrF = FSLik(l,Ei)
                    PrE = PrE * PrF(:,w)
                    if (.not. ALL(PrE==1D0))  PrW(w) = PrW(w) * SUM(PrE)  
                  enddo  
                enddo
                PrXV(x,y,z,v,1) = PrXV(x,y,z,v,1) * SUM(PrW)
              endif
              PrF = FSLik(l,Bi)
              PrW = PrW * PrF(:,x)
              PrXV(x,y,z,v,2) = PrXV(x,y,z,v,2) * SUM(PrW)
            enddo
          endif
        enddo
      enddo
    enddo
  enddo
  PrL(l) = LOG10(SUM(PrXV(:,:,:,:,2))) -  LOG10(SUM(PrXV(:,:,:,:,1))) 
enddo
LL = SUM(PrL)

end subroutine addGAU

! #####################################################################

subroutine ParentHFS(A, SA, kA, SB, kB, hf, LL)  
! parents of SA and SB HS/FS?
use Global
use CalcLik
implicit none

integer, intent(IN) :: A, SA, kA, SB, kB, hf
double precision, intent(OUT) :: LL
integer :: PA, nA, AncA(2,mxA), AA(maxSibSize), G(2), AncB(2,mxA), GA(2), GB(2), &                                                          
  catA(maxSibSize), catB(nS(SB,kB)), catG, GGP(2), DoQuickA, DoQuickB, &
  m, l, x, y, u,v, i, j,z, r, Ei, e, DoneA(MaxSibSize), Bj
double precision :: LLm(2), ALR, PrL(nSnp), PrG(3,2), PrXV(3,3,3,3,3,2), PrPA(3, 2),&
 PrGA(3), PrGB(3), PrE(3), PrH(3), PrGG(3), PrF(3,3)

PA = 0
nA = 0
LLm = missing  
if (A/=0) then
  PA = Parent(A,kA)
  call GetAncest(A, kA, AncA)
else
  PA = -SA
  call GetAncest(-SA, kA, AncA)
endif
if (PA/=0 .and. GpID(kA,SB,kB)==PA) then
  LL = impossible
  return
endif

do m=1,2
  if (m/=hf .and. hf/=3) cycle
  if (PA < 0) then
    if (GpID(m,-PA,kA)/=GpID(m,SB,kB) .and. GpID(m,-PA,kA)/=0 .and. GpID(m,SB,kB)/=0) then
      LLm(m) = impossible
    endif
    nA = nS(-PA,kA)
    AA(1:nA) = SibID(1:nA, -PA, kA)
  else 
    nA = 1
    AA(1) = A
    if (PA > 0) then
      if (Parent(PA,m)/=GpID(m,SB,kB) .and. parent(PA,m)/=0 .and. GpID(m,SB,kB)/=0) then
        LLm(m) = impossible
      endif
    endif
  endif
enddo
if (ALL(LLm == impossible)) then
  LL = impossible
  return
endif

G = 0
call GetAncest(-SB, kB, AncB)
GA = AncA(:,kA+2)
GB = AncB(:,kB+2)
do m=1,2
  if (m/=hf .and. hf/=3) cycle
  if (GA(m)/=0) then
    if (GA(m) == -SB .and. m==kB) then
      LLm(m) = impossible
    else if (GB(m)/=0 .and. GA(m)/=GB(m)) then
      LLm(m) = impossible
    else if (GB(m)==0) then
      G(m) = GA(m)
    else if (GB(m)==GA(m)) then
      G(m) = GA(m)
      LLm(m) = AlreadyAss  ! already are sibs
    else
      LLm(m) = impossible
    endif
  else    
    G(m) = GB(m)
  endif
  if (hf==3) then  ! FS
    if (ANY(AncA(kB, 3:mxA) == -SB)) then
      LLm = impossible
    else if (A>0) then
      if (ANY(AncB == A)) then
        LLm = impossible
      endif
    else if (SA/=0) then
      if (ANY(AncB(kA,3:mxA) == -SA)) then
        LLm = impossible
      endif
    endif
  endif 
  ALR = missing
  if (A>0) then
     if (hf < 3) then
      call CalcAgeLR(A, kA, -SB, kB, m, 6, .TRUE., ALR) 
    else
      call CalcAgeLR(A, kA, -SB, kB, kA, 5, .TRUE., ALR)
    endif
  else if (SA/=0) then
    if (hf < 3) then
      call CalcAgeLR(-SA, kA, -SB, kB, m, 3, .TRUE., ALR)
    else
      call CalcAgeLR(-SA, kA, -SB, kB, 0, 2, .TRUE., ALR)
    endif
  endif
  if (ALR == impossible)  LLm = impossible
enddo
if (ALL(LLm == impossible)) then
  LL = impossible
  return
endif

if (hf==3) then
  if (LLm(1)==impossible .or. LLm(2)==impossible) then
    LL = impossible
  else if (LLm(1)==AlreadyAss .and. LLm(2)==AlreadyAss) then
    LL = AlreadyAss  ! already are FS
  endif
else 
  if (LLm(hf)==impossible) then 
    LL = impossible
  else if (LLm(hf)==AlreadyAss) then
    LL = AlreadyAss
  else if (LLm(3-hf)==AlreadyAss) then
    LL = impossible   ! already HS, would become FS
  endif
endif

if (ANY(AncA(kB, 5:mxA)==-SB)) then
  LL = NotImplemented  ! highly unlikely (but not strictly impossible: TODO)
else if (AncA(kA,2)/=0 .and. ANY(AncB(kA, 5:mxA) == AncA(kA,2))) then 
  LL = impossible
endif
if (LL /=missing) return  
   
catA = 0  
catB = 0
do i=1, nA
  if (hermaphrodites/=0 .and. PA<0) then
    if (Parent(AA(i), 3-kA) == -DumClone(-PA,kA) .and. DumClone(-PA,kA)/=0) then
      catA(i) = 12
      cycle
    endif
  endif 
  if (kA/=kB) then
    if (Parent(AA(i), kB) == AncB(kB, 2) .and. AncB(kB, 2)<0) then
      catA(i) = 1
    endif
  else if (kA == kB .and. Parent(AA(i), 3-kA) /= 0) then  
    do j=1, nS(SB, kB)
      if (Parent(AA(i), 3-kA) == Parent(SibID(j,SB,kB), 3-kB)) then
        catA(i) = 2
        catB(j) = 2
      endif
    enddo
  endif
  if (Parent(AA(i), 3-kA) /= 0) then
    if (G(3-kA) == Parent(AA(i), 3-kA)) then  ! incl. hf==3
      if (kA==kB) then
        catA(i) = 3  ! (u) 3-kA = 3-kB == hf 
      else if (kA/=kB) then
        catA(i) = 4  ! (z)
      endif 
    else if (kA==hf .and. GA(3-kA) == Parent(AA(i), 3-kA)) then
      catA(i) = 4  ! (z)
    else if (kA==hf .and. GB(3-kA) == Parent(AA(i), 3-kA)) then
      catA(i) = 5  ! (v)
    endif
  endif
enddo    

if (Complx<2 .and. (any(catA/=0) .or. any(catB/=0))) then   ! TODO DOUBLE CHECK IF SOME VALID
  LL = NotImplemented
  return
endif

do i=1, nS(SB, kB)
  if (hermaphrodites/=0 .and. DumClone(SB,kB)/=0) then
    if (Parent(SibID(i,SB,kB), 3-kB) == -DumClone(SB,kB)) then
      catB(i) = 12
      cycle
    endif
  endif
  if (kA/=kB) then
    if (Parent(SibID(i,SB,kB), kA) ==AncA(kA,2) .and. AncA(kA,2)<0) then
      catB(i) = 1
    endif
  endif
  if (Parent(SibID(i,SB,kB), 3-kB) /= 0) then
    if (G(3-kB) == Parent(SibID(i,SB,kB), 3-kB)) then
      catB(i) = 3  ! (u)  (for hf<3 .and. hf==3)
    else if (kB==hf .and. GA(3-kB) == Parent(SibID(i,SB,kB), 3-kB)) then
        catB(i) = 4  ! (z) (GA of type 3-kB if hf==kB) 
    else if (kB==hf .and. GB(3-kB) == Parent(SibID(i,SB,kB), 3-kB)) then
      catB(i) = 5  ! (v)
    endif
  endif
enddo 

catG = 0
GGP = 0
if (hf<3) then
  GGP = getPar(G(hf), hf)
  if (GGP(3-hf) == GA(3-hf) .and. GA(3-hf)/=0) then
    catG = 1
  else if (GGP(3-hf) == GB(3-hf) .and. GB(3-hf)/=0) then
    catG = 2
  endif  
  if (catG == 0)  GGP = 0                       
endif

DoQuickA = 1
DoQuickB = 1                               
if (SA/=0)  call ChkDoQuick(SA,kA,DoQuickA)
call ChkDoQuick(SB,kB,DoQuickB)  

PrL = 0D0
do l=1,nSnp
  do m=1,2
    if (m/=hf .and. hf/=3) cycle
    if (ANY(CatA==3) .or. ANY(CatB==3)) then
      call ParProb(l, G(m), m, -1,0, PrG(:,m)) 
    else if (catG/=0) then
      call ParProb(l, G(m), m, -4, 0, PrG(:,m))
      if (G(m) > 0) then
        call ParProb(l, GGP(hf), hf, G(m), 0, PrGG) 
      else
        call ParProb(l, GGP(hf), hf, 0, 0, PrGG)
      endif
    else
      call ParProb(l, G(m), m, 0,0, PrG(:,m)) 
    endif
  enddo
  if (hf < 3) then
    if (ANY(CatA==4) .or. ANY(CatB==4)) then
      call ParProb(l, GA(3-hf), 3-hf, -1,0, PrGA)
    else
      call ParProb(l, GA(3-hf), 3-hf, 0,0, PrGA)
    endif
    if (ANY(CatA==5) .or. ANY(CatB==5)) then
      call ParProb(l, GB(3-hf), 3-hf, -1,0, PrGB)
    else
      call ParProb(l, GB(3-hf), 3-hf, 0,0, PrGB)
    endif
  endif
  if (A>0) then
    call ParProb(l, Parent(A,kA), kA, A,-4, PrPA(:,kA))
    call ParProb(l, Parent(A,3-kA), 3-kA, A,0, PrPA(:,3-kA))
  endif
  
  PrXV = 0D0
  do x=1,3  ! SA/PA
    do y=1,3  ! SB
      do u=1,3  ! G_hf / G_3-kB (hf==3)
        do z=1,3  ! G_A (hf/=3) / G_kB (hf==3)
          do v=1,3 ! G_B (hf/=3)
            if (hf==3) then  ! 0 for z/=v
              PrXV(x,y,u,z,z,:) = AKA2P(x,u,z) * AKA2P(y,u,z) *&
               PrG(u,3-kB) * PrG(z, kB)
            else
              if (GA(3-hf) < 0 .and. GA(3-hf) == -SB .and. hf==3-kB) then
                PrXV(x,y,u,y,v,:) = AKA2P(x,u,y) * AKA2P(y,u,v) *&
                 PrG(u,hf) * PrGB(v)
              else if (GB(3-hf) < 0 .and. GB(3-hf) == -SA .and. hf==3-kA) then
                PrXV(x,y,u,z,x,:) = AKA2P(x,u,z) * AKA2P(y,u,x) *&
                 PrG(u,hf) * PrGA(z)
              else if (catG == 1) then
                PrXV(x,y,u,z,v,:) =AKA2P(x,u,z)*AKA2P(y,u,v)*PrG(u,hf)*&
                  SUM(AKA2P(u,z,:) * PrGG) * PrGA(z) * PrGB(v)
              else if (catG == 2) then
                PrXV(x,y,u,z,v,:) =AKA2P(x,u,z)*AKA2P(y,u,v)*PrG(u,hf)*&
                  SUM(AKA2P(u,v,:) * PrGG) * PrGA(z) * PrGB(v)
              else
                PrXV(x,y,u,z,v,:) = AKA2P(x,u,z)*AKA2P(y,u,v)*PrG(u,hf)&
                * PrGA(z) * PrGB(v)
              endif
            endif
            if (A /=0) then
              if (Parent(A, kA)/=0) then
                PrXV(x,y,u,z,v,:) = PrXV(x,y,u,z,v,:) * PrPA(x, kA)
              endif
            endif
          enddo
        enddo
      enddo
      
      DoneA = 0            
      if ((DoQuickB > 0 .and. ALL(catA==0) .and. ALL(catB==0)) .or. &
        DoQuickB>1 .or. DoQuickA>1) then
        if (SA/=0) then
          PrXV(x,y,:,:,:,2) = PrXV(x,y,:,:,:,2) * XPr(1,x,l, SA,kA) *&
           XPr(1,y,l, SB,kB)
        else if (A>0) then
         PrXV(x,y,:,:,:,2) =PrXV(x,y,:,:,:,2)*SUM(OKA2P(Genos(l,A),x,:)&
          * PrPA(:,3-kA)) * XPr(1,y,l, SB,kB)
        endif            

      else
        do r=1, nS(SB,kB)
          Bj = SibID(r, SB, kB) 
          if (NFS(Bj) == 0) cycle 
          if (catB(r)==0 .or. catB(r)==2) then
            call ParProb(l, Parent(Bj,3-kB), 3-kB, -1, 0, PrE)
          else
            PrE = 1D0
          endif

          if (Parent(Bj,3-kB) <0 .and. CatB(r)/=1) then
            do e=1,3
              if (catB(r)==12 .and. e/=y)  cycle
              do v=1, nS(-Parent(Bj,3-kB), 3-kB)
                Ei = SibID(v, -Parent(Bj,3-kB), 3-kB)
                if (nFS(Ei) == 0) cycle
                if (Parent(Ei, kB) == -SB) cycle
                if (Parent(Ei, kA) == Parent(AA(1),kA) .and. Parent(AA(1),kA)/=0) cycle  ! FS of A if A>0
                call ParProb(l, Parent(Ei, kB), kB, Ei, -1, PrH)
                do i=1,nFS(Ei)
                  if (A>0 .and. (FSID(i,Ei)==A .or. FSID(i,Ei)==Parent(AA(1),kA))) cycle
                  PrH = PrH * OKA2P(Genos(l,FSID(i,Ei)),:,e)
                enddo
                if (.not. all(PrH==1D0))  PrE(e) = PrE(e) * SUM(PrH)
              enddo
            enddo
          endif

          if ((catB(r)==0 .or. catB(r)==2) .and. .not. all(PrE==1D0)) then 
            PrXV(x,y,:,:,:,1) = PrXV(x,y,:,:,:,1) * SUM(PrE)
          else if (catB(r)==1) then  ! Parent(Bj,3-kB) = PA, kA/=kB
            PrXV(x,y,:,:,:,1) = PrXV(x,y,:,:,:,1) * PrE(x)
          else if (catB(r)==3) then  ! hf==kB, Parent(Bj,3-kB) = GA
            do u=1,3
              PrXV(x,y,u,:,:,1) = PrXV(x,y,u,:,:,1) * PrE(u)
            enddo
          else if (catB(r)==4) then  ! Parent(Bj,3-kB) = G
            do e=1,3
              PrXV(x,y,:,e,:,1) = PrXV(x,y,:,e,:,1) * PrE(e)
            enddo
          else if (catB(r)==5) then  ! hf==kB, Parent(Bj,3-kB) = GB
            do v=1,3
              PrXV(x,y,:,:,v,1) = PrXV(x,y,:,:,v,1) * PrE(v)
            enddo
          else if (catB(r)==12) then ! selfing
            PrXV(x,y,:,:,:,1) = PrXV(x,y,:,:,:,1) * PrE(y)
          endif
     
          PrF = FSLik(l,Bj)
          PrE =  PrE * PrF(:,y)

          if (any(catA == 2) .and. Parent(Bj,3-kB)/=0) then  ! kA==kB, share parent 3-kB
            do v = 1, nA
              if (A/=0 .and. AA(v)/=A) cycle
              if (Parent(AA(v), 3-kA)/=Parent(Bj,3-kB)) cycle
              PrE =  PrE * OKA2P(Genos(l,AA(v)), x, :)
              doneA(v) = 1
            enddo
          endif
          
          if ((catB(r)==0 .or. catB(r)==2) .and. .not. all(PrE==1D0)) then
            PrXV(x,y,:,:,:,2) = PrXV(x,y,:,:,:,2) * SUM(PrE)
          else if (catB(r)==1) then  ! Parent(Bj,3-kB) = PA, kA/=kB
            PrXV(x,y,:,:,:,2) = PrXV(x,y,:,:,:,2) * PrE(x)
          else if (catB(r)==3) then  ! hf==kB, Parent(Bj,3-kB) = GA
            do u=1,3
              PrXV(x,y,u,:,:,2) = PrXV(x,y,u,:,:,2) * PrE(u)
            enddo
          else if (catB(r)==4) then  ! Parent(Bj,3-kB) = G
            do e=1,3
              PrXV(x,y,:,e,:,2) = PrXV(x,y,:,e,:,2) * PrE(e)
            enddo
          else if (catB(r)==5) then  ! hf==kB, Parent(Bj,3-kB) = GB
            do v=1,3
              PrXV(x,y,:,:,v,2) = PrXV(x,y,:,:,v,2) * PrE(v)
            enddo
          else if (catB(r)==12) then ! selfing
            PrXV(x,y,:,:,:,2) = PrXV(x,y,:,:,:,2) * PrE(y)
          endif
        enddo
        
        do r = 1, nA
          if (DoneA(r) == 1)  cycle                      
          if (SA/=0 .and. NFS(AA(r)) == 0) cycle         
          if (kA/=kB .and. Parent(AA(r),3-kA)==-SB) cycle  ! done
          if (catA(r)==0) then
            call ParProb(l, Parent(AA(r),3-kA), 3-kA, -1, 0, PrE)
          else
            PrE = 1D0
          endif

          if (Parent(AA(r), 3-kA) < 0 .and. (SA/=0 .or. &
           ANY(FSID(1:nFS(AA(r)), AA(r))==A))) then
            do e=1,3
              if (catA(r)==12 .and. e/=x)  cycle
              do i=1, nS(-Parent(AA(r), 3-kA), 3-kA)
                Ei = SibID(i, -Parent(AA(r), 3-kA), 3-kA)
                if (nFS(Ei) == 0) cycle
                if (Parent(Ei, kB) == -SB) cycle
                if (A>0 .and. Ei==A) cycle
                if (SA/=0 .and. Parent(Ei, kA) == -SA) cycle
                call ParProb(l, Parent(Ei, kA), kA, Ei, -1, PrH)  
                do j=1, nFS(Ei)
                  if (A/=0 .and. FSID(j, Ei)==A) cycle
                  PrH = PrH * OKA2P(Genos(l,FSID(j,Ei)), :, e)
                enddo
                if (.not. all(PrH==1D0))  PrE(e) = PrE(e) * SUM(PrH)
              enddo
            enddo
          endif
          
          if (catA(r)<3 .and. .not. all(PrE==1D0)) then
            PrXV(x,y,:,:,:,1) = PrXV(x,y,:,:,:,1) * SUM(PrE)
          else if (catA(r)==3) then 
            do u=1,3
              PrXV(x,y,u,:,:,1) = PrXV(x,y,u,:,:,1) * PrE(u)
            enddo
          else if (catA(r)==4) then 
            do z=1,3
              PrXV(x,y,:,z,:,1) = PrXV(x,y,:,z,:,1) * PrE(z)
            enddo
          else if (catA(r)==5) then 
            do v=1,3
              PrXV(x,y,:,:,v,1) = PrXV(x,y,:,:,v,1) * PrE(v)
            enddo
          else if (catA(r)==12) then   
            PrXV(x,y,:,:,:,1) = PrXV(x,y,:,:,:,1) * PrE(x)
          endif
          
          do i=1, MAX(1, nFS(AA(r)))
            if (SA/=0 .or. (FSID(i, AA(r))==A .and. catA(r)/=2)) then 
              PrE =  PrE * OKA2P(Genos(l,FSID(i,AA(r))), x, :)
              doneA(r) = 2
            else if (catA(r)/=2) then
              PrE =  PrE * OKA2P(Genos(l,FSID(i,AA(r))), x, :)
            endif
          enddo
          
          if (catA(r)<3 .and. .not. all(PrE==1D0)) then
            PrXV(x,y,:,:,:,2) = PrXV(x,y,:,:,:,2) * SUM(PrE)
           else if (catA(r)==3) then 
            do u=1,3
              PrXV(x,y,u,:,:,2) = PrXV(x,y,u,:,:,2) * PrE(u)
            enddo
          else if (catA(r)==4) then 
            do z=1,3
              PrXV(x,y,:,z,:,2) = PrXV(x,y,:,z,:,2) * PrE(z)
            enddo
          else if (catA(r)==5) then 
            do v=1,3
              PrXV(x,y,:,:,v,2) = PrXV(x,y,:,:,v,2) * PrE(v)
            enddo
          else if (catA(r)==12) then   
            PrXV(x,y,:,:,:,2) = PrXV(x,y,:,:,:,2) * PrE(x)
          endif
        enddo
      endif
    enddo
  enddo
  PrL(l) = LOG10(SUM(PrXV(:,:,:,:,:,2))) - LOG10(SUM(PrXV(:,:,:,:,:,1)))
enddo

LL = SUM(PrL)


! if (A==147 .and. kA==2 .and. SB==12 .and. kB==2 .and. hf==3) then
  ! print *, ""
  ! write (*,'("ParentHFS: ", 6i4, 2i3, f9.2)') A, SA, kA, SB, kB, hf, DoQuickA, DoQuickB, LL
  ! print *, "G: ", G, ", catG: ", catG, "; ", GA, GB
  ! print *, "catA: ", catA(1:nA) , ", nFS: ", nFS(AA(1:nA))
  ! print *, "catB: ", catB
  ! print *, "doneA: ", doneA(1:nA)
! !  print *, "nFS A1: ", nFS(AA(1))
  ! print *, 'Parent A 3-kA: ', Parent(AA(1:nA), 3-kA)
  ! print *, ""
! endif

end subroutine ParentHFS

! #####################################################################

subroutine DummyGP(SA, SB, kA, kB, focal, LL)  
! SB GP of SA? via observed or unobserved
use Global
use CalcLik
implicit none

integer, intent(IN) :: SA, SB, kA, kB, focal
double precision, intent(OUT) :: LL
integer :: i, r, AncB(2,mxA), Bi, AncBi(2, mxA), catA, catB, G(2), catG, &
  DoQuickA, DoQuickB, GGP(2), m, l, x, y, z, v, Ai, ix, e,f, Ei
double precision :: LLGX(2,2), LLtmp(maxSibSize, 2), PrL(nSnp), PrZ(3), &
  PrG(3), PrPG(3), PrXYZ(3,3,3,3,2), PrW(3), dx(maxSibSize), PrH(3), PrF(3,3)
logical :: MaybeGP(maxSibSize), ParAisClone(ns(SA,kA)), ParBisClone(ns(SB,kB))

LL = missing
do i=1, nS(SB,kB)
  Bi = SibID(i, sB, kB)
  if (any(GpID(:,SA,kA) == Bi)) then
    LL = NotImplemented
  else if (kA /= kB) then
    if (Parent(Bi, kA) == -SA) then
      LL = NotImplemented
    endif
  else if (kA == kB) then
    do r= 1, nS(SA, kA)
      if (Parent(Bi, 3-kB)==Parent(SibID(r,SA,kA), 3-kA) &
       .and. Parent(Bi, 3-kB)/=0) then
        LL = NotImplemented   ! TODO
      endif
    enddo
  endif
enddo
if (LL == NotImplemented)  return 

 call GetAncest(-SB, kB, AncB)
if (ANY(AncB(kA,3:mxA) == -SA)) then
  LL = impossible 
  return
else if (ANY(AncB(3-kA, 3:mxA) < 0)) then
  do r= 1, nS(SA, kA)
    if (ANY(AncB(3-kA, 3:mxA)/=0 .and. AncB(3-kA, 3:mxA) == &
     Parent(SibID(r,SA,kA), 3-kA))) then
      LL = NotImplemented
      return
    endif
  enddo
endif

MaybeGP = .TRUE.   ! Bi potential GP of SA?
do r=1, ns(SB,kB)
  Bi = SibID(r, sB, kB)
  if (Parent(Bi, 3-kB)==0) cycle
  call getAncest(Bi, kB, AncBi)
  if (ANY(AncBi(kA,:) == -SA)) then
    MaybeGP(r) = .FALSE.
  endif
enddo

catA = 0
catB = 0
G = GpID(:, SA, kA)
do r = 1, nS(sB,kB)   ! check if overlap
  Bi = SibID(r, sB, kB)
  if (NFS(Bi) == 0) cycle
  if (Parent(Bi, 3-kB) == G(3-kB) .and. G(3-kB)<0) then
    catB = Bi
  endif
enddo
do r = 1, nS(sA,kA)   ! check if inbreeding loop
  Ai = SibID(r,SA,kA)
  if (NFS(Ai)==0) cycle
  if (Parent(Ai, 3-kA) == G(3-kA) .and. G(3-kA)<0) then 
    catA = Ai
  endif
enddo    

catG = 0
do m=1,2
  if (GpID(m, SA, kA) == GpID(m, SB, kB) .and. GpID(m, SA, kA)/=0) then
    catG = m
  endif
enddo

call ChkDoQuick(SA,kA,DoQuickA)
call ChkDoQuick(SB,kB,DoQuickB)

ParAisClone = .FALSE.
ParBisClone = .FALSE.
if (hermaphrodites/=0) then
  if (DumClone(SA,kA)/=0) then   ! DoQuickA==-3
    do i=1, nS(SA,kA)
      if (Parent(SibID(i,SA,kA), 3-kA) == -DumClone(SA,kA)) then
        ParAisClone(i) = .TRUE.
      endif
    enddo
  endif
  if (DumClone(SB,kB)/=0) then   ! DoQuickB==-3
    do i=1, nS(SB,kB)
      if (Parent(SibID(i,SB,kB), 3-kB) == -DumClone(SB,kB)) then
        ParBisClone(i) = .TRUE.
      endif
    enddo
  endif
endif

LLGX = missing  ! D2: 1: a B is GP; 2: an unsampled offspring of SB is GP
LLtmp = missing
GGP = 0
do m=1,2
  if (focal==4 .and. all(GpID(:,SA,kA)==0) .and. complx/=0) then   ! TODO double check  
    LLGX(m,1) = NotImplemented  
    ! else: SB is GP of SA --> Bi is parent of SA --> mate-of-SB is also GP of SA
  else if (m==kB .and. GpID(kB, SA, kA) == -SB) then
    LLGX(m,:) = impossible
    cycle
  else if (G(m) > 0) then
    if (Parent(G(m), kB) /=0) then
      LLGX(m,:) = impossible
    else 
      call AddSib(G(m), SB, kB, LLtmp(1,m))
      call AddFS(G(m), SB, kB,0,m, LLtmp(2,m), ix, dx)
      if (MaxLL(LLtmp(:,m)) < 0D0) then
        LLGX(m,1) = MaxLL(LLtmp(:,m)) + CLL(SA, kA)
        if (Parent(G(m), kB) /= -SB) then
          LLGX(m,1) = LLGX(m,1) - Lind(G(m))
        endif
      else
        LLGX(m,1) = MaxLL(LLtmp(:,m))
      endif
    endif
    cycle
  else if (G(m) == 0) then
    do r=1, nS(sB,kB)
      Bi = SibID(r, sB, kB) 
      if (Sex(Bi)/=m .and. Sex(Bi)<3) cycle
      if (.not. MaybeGP(r)) cycle
      call AddGP(Bi, SA, kA, LLtmp(r,m))
      if (LLtmp(r,m) < 0) then
        LLtmp(r,m) = LLtmp(r,m) - Lind(Bi) + CLL(SB,kB)
      endif
    enddo
    LLGX(m,1) = MaxLL(LLtmp(:,m))
  else if (G(m) < 0) then 
    if (GpID(kB, -G(m), m) /=0) then
      LL = impossible
    else
      GGP(m) = GpID(3-kB, -G(m), m)
    endif
  endif
  
  if (Complx==0 .and. GGP(m)==0) then
    GGP(m) = DumMate(SB,kB)
  endif
  
  PrL = 0D0
  do l=1,nSnp
    if ((catB /= 0 .and. m==kB) .or. catA/=0) then
      call ParProb(l, G(3-m), 3-m, -1, 0, PrZ)
    else
      call ParProb(l, G(3-m), 3-m, 0, 0, PrZ)
    endif
    if (catB/=0 .and. m==3-kB .and. G(m)<0) then
      PrG = 1D0
    else
      call ParProb(l, G(m), m, -4, 0, PrG)  ! G(m)'s offspring if<0, else 1D0
    endif
    call ParProb(l, GGP(m), 3-kB, 0, 0, PrPG)

    PrXYZ = 0D0
    do x=1,3  ! SA
      do y=1,3  ! parent of SA, offspr of SB
        do z=1,3  ! other parent of SA
          PrXYZ(x,y,z,:,:) = AKA2P(x, y, z) * PrG(y) * PrZ(z)
        enddo

        if ((catA==0 .and. DoQuickA>0) .or. DoQuickA>1) then
          PrXYZ(x,y,:,:,2) = PrXYZ(x,y,:,:,2) * XPr(1,x,l, SA,kA)   
        else
          do z=1,3
            do r=1, nS(sA,kA)
              Ai = SibID(r, sA, kA)  
              if (NFS(Ai) == 0) cycle
              if (Ai == G(m)) cycle
              if (catA==Ai .or. ParAisClone(r)) then
                PrW = 1D0
              else
                call ParProb(l, Parent(Ai, 3-kA), 3-kA, -1, 0, PrW)
              endif
              
              if (Parent(Ai,3-kA) < 0) then
                do e=1,3
                  do f=1, nS(-Parent(Ai, 3-kA), 3-kA)
                    Ei = SibID(f, -Parent(Ai, 3-kA),3-kA)
                    if (nFS(Ei) == 0) cycle
                    if (Parent(Ei,kA) == -SA)  cycle
 !                   if (kA==kB .and. Parent(Ei,kB) == -SB)  cycle
                    call ParProb(l,Parent(Ei,kA),kA,Ei,-1,PrH) 
                    PrF = FSLik(l,Ei)
                    PrH = PrH * PrF(:,x)
                    if (.not. all(PrH==1D0))  PrW(e) = PrW(e) * SUM(PrH)
                  enddo
                enddo
              endif
              
              if (catA==Ai) then
                if (m==kA) then
                  PrXYZ(x,y,z,:,1) = PrXYZ(x,y,z,:,1) * PrW(z)
                else if (m/=kA) then
                  PrXYZ(x,y,z,:,1) = PrXYZ(x,y,z,:,1) * PrW(y)
                endif
              else if (ParAisClone(r)) then
                PrXYZ(x,y,z,:,1) = PrXYZ(x,y,z,:,1) * PrW(x)
              else if (.not. all(PrW==1D0)) then
                PrXYZ(x,y,z,:,1) = PrXYZ(x,y,z,:,1) * SUM(PrW)
              endif
              
              PrF = FSLik(l,Ai)
              PrW = PrW * PrF(:,x)
              
              if (catA==Ai) then
                if (m==kA) then
                  PrXYZ(x,y,z,:,2) = PrXYZ(x,y,z,:,2) * PrW(z)
                else if (m/=kA) then
                  PrXYZ(x,y,z,:,2) = PrXYZ(x,y,z,:,2) * PrW(y)
                endif
              else if (ParAisClone(r)) then
                PrXYZ(x,y,z,:,2) = PrXYZ(x,y,z,:,2) * PrW(x)
              else if (.not. all(PrW==1D0)) then
                PrXYZ(x,y,z,:,2) = PrXYZ(x,y,z,:,2) * SUM(PrW)
              endif
            enddo  ! r
          enddo  ! z
        endif
        
        do v=1,3  ! SB
          if (DoQuickB==-2 .or. Complx==0) then  ! all FS
            do r=1, nS(sB,kB)
              Bi = SibID(r, sB, kB)  
              if (NFS(Bi) > 0) exit
            enddo 
            call ParProb(l, Parent(Bi, 3-kB), 3-kB, -1, 0, PrW)
            PrXYZ(x,y,:,v,:) = PrXYZ(x,y,:,v,:) * XPr(2,v,l, SB,kB)  ! GPs of SB
            
            if (Complx==0 .or. Parent(Bi,3-kB)== GGP(m)) then
              PrW = PrW * AKA2P(y, v, :)
              PrXYZ(x,y,:,v,1) = PrXYZ(x,y,:,v,1) * SUM(PrW)
            else 
              PrXYZ(x,y,:,v,:) = PrXYZ(x,y,:,v,:) * SUM(AKA2P(y,v,:) * PrPG) 
            endif
            
            PrF = FSLik(l,Bi)
            PrW = PrW * PrF(:,v)
            PrXYZ(x,y,:,v,2) = PrXYZ(x,y,:,v,2) * SUM(PrW)
        
          else if (catG==0 .and. DoQuickB/=-3 .and. ns(SB,kB)>1) then
            PrXYZ(x,y,:,v,:) = PrXYZ(x,y,:,v,:) * SUM(AKA2P(y, v, :) * PrPG)
            PrXYZ(x,y,:,v,1) = PrXYZ(x,y,:,v,1) * DumP(v,l, SB,kB) 
            PrXYZ(x,y,:,v,2) = PrXYZ(x,y,:,v,2) * XPr(3,v,l, SB,kB)  ! GPs & sibs
          
          else  ! SA inbred or selfing
          
            PrXYZ(x,y,:,v,:) = PrXYZ(x,y,:,v,:) * SUM(AKA2P(y, v, :) * PrPG)   
            call ParProb(l, GpID(m,SB,kB), m, 0, 0, PrH)
            do z=1,3
              PrXYZ(x,y,z,v,:) = PrXYZ(x,y,z,v,:) * SUM(AKA2P(v,z,:) *PrH)
            enddo

            do r=1, nS(sB,kB)
              Bi = SibID(r, sB, kB)  
              if (NFS(Bi) == 0) cycle
              if (catB==Bi .or. ParBisClone(r)) then
                PrW = 1D0
              else
                call ParProb(l, Parent(Bi, 3-kB), 3-kB, -1, 0, PrW)
              endif
              
              if (Parent(Bi,3-kB) < 0) then
                do e=1,3
                  do f=1, nS(-Parent(Bi, 3-kB), 3-kB)
                    Ei = SibID(f, -Parent(Bi, 3-kB),3-kB)
                    if (nFS(Ei) == 0) cycle
                    if (Parent(Ei,kB) == -SB)  cycle
!                    if (kA==kB .and. Parent(Ei,kA) == -SA)  cycle
                    call ParProb(l,Parent(Ei,kB),kB,Ei,-1,PrH) 
                    PrF = FSLik(l,Ei)
                    PrH = PrH * PrF(:,v)
                    if (.not. all(PrH==1D0))  PrW(e) = PrW(e) * SUM(PrH)
                  enddo
                enddo
              endif
              
              if (catB==Bi .and. m==kB) then
                do z=1,3
                  PrXYZ(x,y,z,v,1) = PrXYZ(x,y,z,v,1) * PrW(z)
                enddo
              else if (catB==Bi .and. m/=kB) then
                PrXYZ(x,y,:,v,1) = PrXYZ(x,y,:,v,1) * PrW(y)
              else if (ParBisClone(r)) then
                PrXYZ(x,y,:,v,1) = PrXYZ(x,y,:,v,1) * PrW(v)
              else if (.not. all(PrW==1D0)) then
                PrXYZ(x,y,:,v,1) = PrXYZ(x,y,:,v,1) * SUM(PrW)
              endif
              
              PrF = FSLik(l,Bi)
              PrW = PrW * PrF(:,v) 
              
              if (catB==Bi .and. m==kB) then
                do z=1,3
                  PrXYZ(x,y,z,v,2) = PrXYZ(x,y,z,v,2) * PrW(z)
                enddo
              else if (catB==Bi .and. m/=kB) then
                PrXYZ(x,y,:,v,2) = PrXYZ(x,y,:,v,2) * PrW(y)
              else if (ParBisClone(r)) then
                PrXYZ(x,y,:,v,2) = PrXYZ(x,y,:,v,2) * PrW(v)
              else if (.not. all(PrW==1D0)) then
                PrXYZ(x,y,:,v,2) = PrXYZ(x,y,:,v,2) * SUM(PrW)
              endif
              
            enddo  ! r
          endif
        enddo  ! v
      enddo  ! y
    enddo  ! x
    PrL(l) = LOG10(SUM(PrXYZ(:,:,:,:,2))) - LOG10(SUM(PrXYZ(:,:,:,:,1)))        
  enddo
  LLGX(m,2) = SUM(PrL)
enddo
LL = MaxLL((/LLGX(:,1), LLGX(:,2)/))

! if (SA==30 .and. kA==2 .and. SB==40 .and. kB==1) then
  ! write(*,'("DummyGP: ", 4i6, 4f8.1, 2i6)') SA, kA, SB, kB, LLGX(:,1), LLGX(:,2), GGP
  ! print *, "cats: ", catA, catB, catG, DoQuickA, DoQuickB
! endif

end subroutine DummyGP

! ######################################################################

subroutine dummyHFA(SA,kA,SB,kB, hf, LL)   ! SB (not Bi's) is HA/FA of SA (not Ai's)
use Global
use CalcLik
implicit none

integer, intent(IN) :: SA, kA, SB, kB, hf
double precision, intent(OUT) :: LL
integer :: i, r, AncB(2,mxA), m, DoQuickA, DoQuickB, l, x, y, z, u,v, Ai
double precision :: PrL(nSnp), PrGG(3,2), PrGA(3), PrXY(3,3,3,3,3), PrA(3), PrPA(3), PrF(3,3)

if (nS(SA,kA)==0 .or. ns(SB,kB)==0) then
  LL = NotImplemented
  return
endif

LL = missing
do i=1, nS(SB,kB)
  if (kA /= kB) then
    if (Parent(SibID(i,SB,kB), kA) == -SA) then
      LL = NotImplemented
      return      
    endif
  else if (kA == kB) then
    do r= 1, nS(SA, kA)
      if (Parent(SibID(i,SB,kB), 3-kB)==Parent(SibID(r,SA,kA), 3-kA) &
       .and. Parent(SibID(i,SB,kB), 3-kB)/=0) then
        LL = NotImplemented
        return  
      endif
    enddo
  endif
enddo

call GetAncest(-SB, kB, AncB)
if (ANY(AncB(kA,3:mxA) == -SA)) then
  LL = impossible 
  return
endif

if (all(GpID(:,SA,kA)/=0)) then
  LL = NotImplemented
  return
else 
  do m=1,2
    if (GpID(m,SA,kA)/=0 .and. GpID(m,SA,kA)==GpID(m,SB,kB)) then
      LL = NotImplemented
      return
    endif
  enddo
endif

if (GpID(1,SA,kA)==0) then
  m = 1
else
  m = 2
endif

call ChkDoQuick(SA,kA,DoQuickA)
call ChkDoQuick(SB,kB,DoQuickB)
if (DoQuickA==-1 .or. DoQuickA==-3 .or. DoQuickB==-1 .or. DoQuickB==-3) then
  LL = NotImplemented
  return
endif

PrL = 0D0
do l=1, nSnp
  do x=1,2
    call ParProb(l, GpID(x,SB,kB), x, 0, 0, PrGG(:,x))
  enddo
  call ParProb(l, GpID(3-m,SA,kA), 3-m, 0,0, PrGA)
  if (DoQuickA==-2) then    ! all FS, incl.  ns(SA,kA)==1 
    do z=1,ns(SA,kA)
      if (nFS(SibID(z,SA,kA))==0)  cycle      
      Ai = SibID(z,SA,kA)
      exit
    enddo
    call ParProb(l, Parent(Ai,3-kA),3-kA, Ai, -1, PrPA)   ! excl Ai & all its FS
    PrF = FSLik(l,Ai)
    do x=1,3
      PrA(x) = SUM(PrF(:,x) * PrPA)  
    enddo
  else
    PrA = XPr(1,:,l,SA,kA) 
  endif
  
  do x=1,3  ! SA
    do y=1,3  ! GpID(m,SA,kA)
      do z=1,3  ! SB
        do u=1,3  ! GpID(1,SB,kB)
          do v=1,3  ! GpID(2,SB,kB)
            PrXY(x,y,z,u,v) = PrA(x) * SUM(AKA2P(x,y,:) * PrGA) * &
              XPr(1,z,l,SB,kB) * AKA2P(z,u,v) * PrGG(u,1) * PrGG(v,2)
            if (hf==3) then
              PrXY(x,y,z,u,v) = PrXY(x,y,z,u,v) * AKA2P(y,u,v)
            else if (hf==1) then
              PrXY(x,y,z,u,v) = PrXY(x,y,z,u,v) * AKAP(y,u,l)
            else if (hf==2) then
              PrXY(x,y,z,u,v) = PrXY(x,y,z,u,v) * AKAP(y,v,l)
            endif
          enddo
        enddo
      enddo
    enddo
  enddo
  PrL(l) = LOG10(SUM(PrXY))
enddo

LL = SUM(PrL)
  
end subroutine dummyHFA

! ######################################################################
!       Age priors  
! ######################################################################

subroutine getEstBY(A, kA, lvl, BYLR)
use Global
implicit none
! lvl: 1=own; 2=own + est. from relatives exact BY + parent yearlast; 
! 3= 2 + est. BY from par + GP, 4= 2 + est. BY from offspr, 5 = all (NOT est BY from sibs)

integer, intent(IN) :: A, kA, lvl
double precision, intent(OUT) :: BYLR(nYears)

BYLR = LOG10(zero)
if (A > 0) then
  if (BY(A) > 0) then
    BYLR(BY(A)) = zero
  else
    BYLR = IndBY(:, A, lvl)
  endif
else if (A < 0) then
  BYLR = DumBY(:, -A, kA, lvl)
endif

end subroutine getEstBY

! ######################################################################

subroutine setEstBY(A, k)
use Global
implicit none
   
! birth year probability distribution, based on offspring, parent & GP BY. (& own min/max) 
! updates IndBY(:,A,:) or DumBY(:,-A, k,:) as side effect, with D3 resp D4:
! 1: own exact BY or BYrange
! 2: own + contributions from Par, GP, Off, sibs with exact BY + YearLast from Par + GP
! 3: [2] + contr from est. BY from par + GP
! 4: [2] + contr from est. BY from offspr
! 5: [3] + contr from est. BY from offspr  (= all)

integer, intent(IN) :: A, k
integer :: w
double precision :: LPBY(nYears, 5), BYup(nYears,2), BYdown(nYears,2), BYsibs(nYears)

if (A == 0) then
  return
else if (A > 0) then
  if (BY(A) > 0)  return  ! birth year known exactly
endif

LPBY = zero
if (A > 0) then
  LPBY(:,1) = IndBY(:, A, 1)   ! own exact BY / BYrange only
else !if (A < 0) then
  LPBY(:,1) = LOG10(1.0D0/(nYears-1))   ! dummy --> by definition no own BY info
  LPBY(nYears,1) = LOG10(zero)   ! dummy can't be 'born' in last year
endif

BYup = 0D0
BYdown = 0D0
BYsibs = 0D0
call CalcBYup(A, k, BYup)    ! info from parents + grandparents
call CalcBYdown(A, k, BYdown)  ! info from offspring
call CalcBYsibs(A, k, BYsibs) 

LPBY(:,2) = LPBY(:,1) + BYup(:,1) + BYdown(:,1) + BYsibs  ! relatives with exact BY
LPBY(:,3) = LPBY(:,1) + BYup(:,2) + BYdown(:,1) + BYsibs
LPBY(:,4) = LPBY(:,1) + BYup(:,1) + BYdown(:,2) + BYsibs
LPBY(:,5) = LPBY(:,1) + BYup(:,2) + BYdown(:,2) + BYsibs  ! everything

! scale to sum to unity 
LPBY = 10**LPBY
do w=2,5  
  if (SUM(LPBY(:,w)) > 0D0) then
    LPBY(:,w) = LPBY(:,w) / SUM(LPBY(:,w)) 
  endif
enddo
LPBY = LOG10(LPBY)

do w=2,5
  if (A > 0) then
    IndBY(:, A, w) = LPBY(:,w)  
  else !if (A < 0) then
    DumBY(:, -A, k, w) = LPBY(:,w)
  endif
enddo

! if (ALL(LPBY(:,3:5) < -HUGE(0D0)) .or. any(LPBY(:,5)/=LPBY(:,5))) then
! !if (A==4432) then
  ! open (unit=51,file="BYprobs.txt",status="unknown", position="append")
  ! write (51, *) ""
  ! if (A>0) then
    ! write (51, *) "LPBY ", A, Parent(A,:), GPID(:, 171,1)
  ! else
    ! write (51, *) "LPBY ", A, GpID(:,-A,k) 
    ! write(51,*)  "Off: ", SibID(1:ns(-A,k),-A,k)
  ! endif
  ! write(51,'(5x, 100i9)')  (/ (w, w=BYzero+1, BYzero+nYears) /)
  ! write (51, '(i5, 100f9.5)') 5, 10**LPBY(:,5)
  ! write (51, *) ""
    ! write (51, '("self ", 100f9.5)') 10**LPBY(:,1)
    ! write (51, '("up   ", 100f9.5)') 10**BYup(:,2)
    ! write (51, '("down ", 100f9.5)') 10**BYdown(:,2)
    ! write (51, '("sibs ", 100f9.5)') 10**BYsibs(:)
    ! write (51, *) ""
    ! close(51)
    ! stop
! endif

end subroutine setEstBY

! ######################################################################

subroutine CalcBYup(A, kA, BYA)
! BY probabilities of indiv A based on its parents + grandparents
! NOTE: not scaled!
use Global
implicit none

integer, intent(IN) :: A, kA
double precision, intent(OUT) :: BYA(nYears, 2)   ! D2: exact BY only; all BY
integer :: Par(2), GP(2,2), y, m, g, x, Ylast, tmpBY
double precision :: BYP(nYears, 2), BYG(nYears, 2,2), tmpX(nYears)

Par = getPar(A, kA)
if (ALL(par == 0)) then  
  BYA = zero
  return
endif

GP = 0
GP(:,1) = getPar(Par(1), 1)
GP(:,2) = getPar(Par(2), 2)

! get current value of parent's & grandparents BY estimates
BYP = LOG10(zero)
BYG = LOG10(zero)
do m=1,2
  if (Par(m)==0)  cycle
  call getEstBY(Par(m), m, 3, BYP(:,m))  ! self + exact + parents + GP
  do g=1,2
    if (GP(g,m)==0)  cycle
    call getEstBY(GP(g,m), g, 3, BYG(:,m,g))
  enddo
enddo

! log --> not-log
BYP = 10**BYP
BYG = 10**BYG

! last possible BY based on parents' & GP's last year of reproduction
Ylast = 999
if (AnyYearLast) then
  do m=1,2
    if (Par(m) == 0)  cycle
    if (Par(m) > 0)  Ylast = MIN(Ylast, YearLast( Par(m) ) )
    do g=1,2
      if (GP(g,m) > 0)  Ylast = MIN(Ylast, YearLast( GP(g,m) ) +MaxAgePO )  
    enddo
  enddo
endif

BYA = zero
BYA(1,:) = LOG10(zero)  ! has parents --> cannot be born in year 1 
if (any(GP /= 0))  BYA(2,:) = LOG10(zero)
do m=1,2
  if (Par(m)==0)  cycle
  if (Par(m) < 0 .and. all(Gp(:,m)==0))  cycle 
  do y=2,nYears  
    if (y > Ylast) then
      BYA(y,:) = LOG10(zero)  ! year y after parents last possible year of reprod
    else if (ANY(BYP(:,m)>=1D0)) then  ! parent has exact BY
      tmpBY = MAXLOC(BYP(:,m), DIM=1)   ! Par(m) may be dummy --> cannot use BY(Par(m)) 
      BYA(y,:) = BYA(y,:) + getAP(y - tmpBY, 1, m, 0, LOG10(zero))
    else  
      ! weighed sum over all possible parent birth years      
      tmpX = 0D0
      do x=1, y-1 
        tmpX(x) = BYP(x,m) * 10**getAP(y-x, 1, m, 0, LOG10(zero))  ! parent born in year x
      enddo
      BYA(y,2) = BYA(y,2) + LOG10(SUM(tmpX))
    endif
    if (BYA(y,2) < -HUGE(0.0D0)) exit
    
    ! grandparents
    if (y==2)  cycle 
    do g=1,2
      if (GP(g,m) == 0)  cycle
      if (ANY(BYG(:,m,g)>=1D0)) then  
        tmpBY = MAXLOC(BYG(:,m,g), DIM=1)
        BYA(y,:) = BYA(y,:) + getAP(y - tmpBY, 4, m, g, LOG10(zero))
      else
        tmpX = 0D0
        do x=1, y-2 
          tmpX(x) = BYG(x,m,g) * 10**getAP(y-x, 4, m, g, LOG10(zero))  ! grandparent born in year x
        enddo 
        BYA(y,2) = BYA(y,2) + LOG10(SUM(tmpX))
      endif
    enddo

  enddo
enddo

! scale 
! BYA = 10**BYA
! do m=1,2
  ! if (SUM(BYA(:,m)) > 0D0) then
    ! BYA(:,m) = BYA(:,m) / SUM(BYA(:,m)) 
  ! endif
! enddo
! BYA = LOG10(BYA)

end subroutine CalcBYup

! ######################################################################

subroutine CalcBYdown(A, kA, BYA)
! BY probabilities of indiv A based on its offspring
! NOTE: not scaled!
use Global
implicit none

integer, intent(IN) :: A, kA
double precision, intent(OUT) :: BYA(nYears, 2)   ! D2: exact BY only; all BY
integer :: nOff, Offspr(maxSibSize), sxOff(maxSibSize), y, x, i, tmpBY
double precision :: tmpX(nYears)
double precision, allocatable :: BYO(:,:)

allocate(BYO(nYears, maxSibSize))                                         

call getOff(A,kA, .TRUE., nOff, Offspr, sxOff)
if (nOff == 0) then
  BYA = zero          
  return
endif

BYO = LOG10(zero)  ! number of offspring born in year y 
if (nOff > 0) then
  do i=1, nOff
    call getEstBY(Offspr(i), sxOff(i),4, BYO(:,i))   ! self + exact + offspring
  enddo
endif

BYO = 10**BYO
  
BYA = zero   
BYA(nYears,:) = LOG10(zero)  ! has offspring --> cannot be born in last year      
do y=1, nYears-1   ! A's BY
  do i=1, nOff
    if (ANY(BYO(:,i)>=1D0)) then  ! offspring i has exact BY
      tmpBY = MAXLOC(BYO(:,i), DIM=1)   ! Off(i) may be dummy --> cannot use BY(Off(i))
      BYA(y,:) = BYA(y,:) + getAP(tmpBY - y, 1, kA, 0, LOG10(zero))
    else
      ! weighed sum over all possible offspring birth years x
      tmpX = 0D0
      do x=y+1, nYears  ! offspring BY
        tmpX(x) = BYO(x,i) * 10**getAP(x-y, 1, kA, 0, LOG10(zero))  
      enddo 
      BYA(y,2) = BYA(y,2) + LOG10(SUM(tmpX))
    endif
    if (BYA(y,2) < -HUGE(0.0D0)) exit  ! e.g. i born in/prior to year y - no need to look at other offspr
  enddo
enddo


! scale 
! BYA = 10**BYA
! do i=1,2
  ! if (SUM(BYA(:,i)) > 0D0) then
    ! BYA(:,i) = BYA(:,i) / SUM(BYA(:,i)) 
  ! endif
! enddo
! BYA = LOG10(BYA)

! if (A==4432) then
  ! open (unit=42,file="log.txt",status="unknown", position="append")
  ! write (42, *) ""
  ! write (42, *) "BYdown ", A
  ! write (42, '(i4, " Off: ", 20i6)')  nOff, Offspr(1:nOff)
  ! write (42, '("down1", 100f8.2)') 10**BYA(1:nYears,1)
  ! write (42, '("down2", 100f8.2)') 10**BYA(1:nYears,2)
  ! close(42)
! !  stop
! endif

deallocate(BYO)

end subroutine CalcBYdown

! ######################################################################

subroutine CalcBYsibs(A, kA, BYA)
! BY probabilities of indiv A based on its siblings (mat + pat + full)
! NOTE: not scaled!    
! NOTE2: sibs with exact BY only, to avoid double contributions of shared parents/shared offspring
use Global
implicit none

integer, intent(IN) :: A, kA
double precision, intent(OUT) :: BYA(nYears) 
integer :: Par(2), nSibs(3), Sibs(maxSibSize, 3), sxSibs(maxSibSize), m, y, i, AgeD
double precision :: BYtmp(nYears, 3)

BYA = zero 
Par = getPar(A, kA)
if (ALL(par == 0)) then  ! no parents --> no siblings.  
  return
endif

nSibs = 0
sibs = 0
do m=1,2
  call getOff(Par(m),m, .FALSE., nSibs(m), Sibs(:,m), sxSibs)  ! non-dummy sibs only
enddo

if (ALL(nSibs(1:2) <= 1)) then   ! no siblings (only offspring of parents = focal indiv, if not dummy)
  return
endif

! FS: intersect between mat + pat & exclude self
do m=1,2
  do i=1, nSibs(m)
    if (Sibs(i,m) == A .or. BY(sibs(i,m)) <0) then  ! do not use if unknown/uncertain BY
      Sibs(i,m) = 0
    else if (Parent(Sibs(i,m),1) == Par(1) .and. Parent(Sibs(i,m),2) == Par(2) .and. &
      Par(1)/=0 .and. Par(2)/=0) then
      if (m==1) then
        nSibs(3) = nSibs(3) +1
        Sibs(nSibs(3), 3) = Sibs(i,m) 
      endif
      Sibs(i,m) = 0   ! only half siblings for pat/mat to avoid double counting
    endif
  enddo
enddo

BYtmp = zero                
do m=1,3
  if (nSibs(m) == 0) then
    BYtmp(:,m) = LOG10(1.0D0/nYears)
  else
    do y=2, nYears
      do i=1, nSibs(m)
         if (Sibs(i,m) == 0)  cycle
        if (BY(Sibs(i,m)) > 0) then  ! sibling i has known BY
          AgeD = ABS(BY(Sibs(i,m)) - y)  ! absolute age difference
          if (m < 3) then  ! half sibs
            if (Par(m) < 0) then
              BYtmp(y,m) = BYtmp(y,m) + getAP(AgeD, 3, 0, m, LOG10(zero))
            else ! genotyped parent: sib genotypes not in likelihood --> unbalanced if multiplying
              BYtmp(y,m) = BYtmp(y,m) + 10**getAP(AgeD, 3, 0, m, LOG10(zero))
            endif
          else  ! full sibs
            if (any(Par < 0)) then
              BYtmp(y,m) = BYtmp(y,m) + getAP(AgeD, 2, 0, 0, LOG10(zero))
            else
              BYtmp(y,m) = BYtmp(y,m) + 10**getAP(AgeD, 3, 0, m, LOG10(zero))
            endif
          endif
        endif
      enddo
    enddo
  endif
enddo

! scale 
do m=1,3
  if (m < 3) then
    if (Par(m) <= 0 .or. COUNT(sibs(:,m)/=0)==0)  BYtmp(:,m) = 10**BYtmp(:,m)
  else
    if (any(Par <= 0) .or. COUNT(sibs(:,m)/=0)==0)  BYtmp(:,m) = 10**BYtmp(:,m)
  endif
  if (SUM(BYtmp(:,m)) > 0D0) then
   BYtmp(:,m) = BYtmp(:,m) / SUM(BYtmp(:,m)) 
  endif
enddo
BYtmp = log10(BYtmp)

BYA = BYtmp(:,1) + BYtmp(:,2) + BYtmp(:,3)

BYA = 10**BYA
if (SUM(BYA) > 0D0) then
  BYA = BYA / SUM(BYA) 
endif
BYA = LOG10(BYA)

! if (ALL(BYA < -HUGE(0D0)) .or. any(BYA/=BYA)) then
! !if (A==6145) then
  ! open (unit=51,file="BYprobs.txt",status="unknown", position="append")
  ! write(51, '("A: ", 4i6)') A, kA, getPar(A,kA)
  ! write (51, '("sib-1", 100f9.5)') 10**BYtmp(:,1)
  ! write (51, '("sib-2", 100f9.5)') 10**BYtmp(:,2)
  ! write (51, '("sib-3", 100f9.5)') 10**BYtmp(:,3)
  ! write (51, *) ""
  ! do m=1,3
    ! write (51, '(i3,i4, " sibs: ", 100i7)') m, nSibs(m), Sibs(1:nSibs(m),m)
  ! enddo
  ! write (51, *) ""
  ! close(51)
! endif

end subroutine CalcBYsibs

! #####################################################################

subroutine EstAgeDif(A, kA, B, kB, AgeD)   ! estimate age difference, incl. from BYrange. Only called by CalcU()
use Global
implicit none

integer, intent(IN) :: A, kA, B, kB
double precision, intent(OUT) :: AgeD
integer :: y, x
double precision :: pBY(nYears, 2)
double precision, allocatable :: ADtmp(:,:)

allocate(ADtmp(nYears, nYears))                             

if (A>0 .and. B>0) then
  if (AgeDiff(A,B) < 999) then
    AgeD = REAL(AgeDiff(A,B), 8)
  endif
endif

pBY = LOG10(zero)
call getEstBY(A, kA, 5, pBY(:, 1))  ! all contributions from all relatives
call getEstBY(B, kB, 5, pBY(:, 2))
pBY = 10**pBY  ! log -> regular scale

ADtmp = 0D0
do x=1, nYears  ! A 
  if (pBY(x,1) < TINY(0.0D0)) cycle  ! no chance that A is born in year x
  do y=1,nYears  ! B
    if (pBY(y,2) < TINY(0.0D0)) cycle
    ADtmp(x,y) = pBY(x,1) * pBY(y,2) * (x-y)   ! if B older, than AgeD > 0
  enddo
enddo

AgeD = SUM(ADtmp)

deallocate(ADtmp)                 

end subroutine EstAgeDif

! #####################################################################

subroutine CalcAgeLR(A, kA, B, kB, m, focal, AllDumRel, ALR) ! m: mat/pat relatives
use Global
implicit none

integer, intent(IN) :: A, kA, B, kB, m, focal
logical, intent(IN) :: AllDumRel
double precision, intent(OUT) :: ALR
integer :: AB(2), kAB(2),  x, y, i, n, fcl, YearLast_B
double precision :: BYLR(nYears, 2), ALRm(2)
double precision, allocatable :: ALRtmp(:,:)

if (.not. ANY((/-1,1,2,3,4,5,6/) == focal))  call Erstop('CalcAgeLR: illegal focal', .TRUE.)                                                                                            
                         
AB = (/ A, B /)
kAB = (/ kA, kB /)
ALR = zero
if (A==0 .or. B==0) then
  return  
else if (A>0 .and. B>0) then
  if (focal==1 .and. BY(A)>0 .and. AnyYearLast) then
    if (BY(A) > YearLast(B)) then  ! YearLast: unknown = +999
      ALR = impossible
      return
    endif
  endif
  
  if (m < 3) then  ! incl. m=0
    ALR = getAP(AgeDiff(A,B), focal, kB, m, LOG10(zero))
  else 
    ALRm = LOG10(zero)                                       
    do n=1,2
      ALRm(n) = getAP(AgeDiff(A,B), focal, kB, n, LOG10(zero))
    enddo
    ALR = MAXVAL(ALRm)
  endif
  if (ALR < -HUGE(0.0D0)) then
    ALR = impossible
    return
  endif
 
  if ((focal==2 .or. focal==3) .and. (any(Parent(A,:)/=0) .or. any(Parent(B,:)/=0))) then
    do i=1,2
      do n=1,2
        if (n/=m .and. focal/=2)  cycle
        if (Parent(AB(i),n) > 0) then
          ALR = getAP(AgeDiff(AB(3-i),Parent(AB(i),n)), 1, n, 0, Impossible)
          if (ALR == Impossible)  return   ! NOT else ALR = ALR * ALRp  - keep it simple. 
        endif
      enddo
    enddo
  endif
  if (AgeDiff(A,B) /= 999)  return
endif  

ALR = Missing
fcl = focal
if (focal==1) then   ! short-cut instead of via dummy BYLR, faster
  if (A>0 .and. B<0) then
    if (BY(A)>0 .and. ns(-B,kB)==0) then
      do n=1,2
        if (GpID(n,-B,kB)>0 .and. GpID(3-n,-B,kB)==0) then
          if (BY(GpID(n,-B,kB))>0) then
            ALR = getAP(AgeDiff(A, GpID(n,-B,kB)), 4, n, kB, Impossible)
          endif
        endif
      enddo
    else if (ns(-B,kB)==1 .and. all(GpID(:,-B,kB)==0)) then
      if (BY(A)>0 .and. BY(SibID(1,-B,kB))>0) then
        ALR = getAP(AgeDiff(A,SibID(1,-B,kB)), 3, 0, kB, Impossible)
      else
        AB(2) = SibID(1,-B,kB)
        fcl = 3       
      endif
    endif
  else if (A<0 .and. B>0) then
    if (ns(-A,kA)==1 .and. all(GpID(:,-A,kA)==0)) then
      if (BY(B)>0 .and. BY(SibID(1,-A,kA))>0) then
        ALR = getAP(AgeDiff(SibID(1,-A,kA),B),4,kB,kA, Impossible)
      else
        AB(1) = SibID(1,-A,kA)
        fcl = 3 
      endif   
    endif
  endif  
endif
if (ALR /= Missing)  return

BYLR = LOG10(zero)  ! likelihood ratio to be born in year X
if (AllDumRel) then
  do i=1,2
    call getEstBY(AB(i), kAB(i), 5, BYLR(:, i))  ! all contributions from all relatives
  enddo
else 
  call getEstBY(A, kA, 4, BYLR(:, 1))  ! excl contributions from est. BY from par + GP
  call getEstBY(B, kB, 3, BYLR(:, 2))  ! excl contributions from est. BY from offspring
endif

do i=1,2
  if (ALL(BYLR(:, i) <= LOG10(1.D0/nYears))) then
    ALR = 0D0  ! no age info available
    return
  endif
enddo  

if (fcl==1 .or. fcl==4) then  ! quick check
  do y=2, nYears  ! B 
    if (BYLR(y,2) < -HUGE(0.0D0)) cycle
    ! at oldest possible BY of B:
    if (ALL(BYLR((y-1):nYears, 1) < -HUGE(0.0D0))) then 
      ALR = impossible
      return
    else
      exit
    endif
  enddo
endif

allocate(ALRtmp(nYears, nYears))  
ALRtmp = LOG10(zero)  ! -Inf
ALRm = LOG10(zero)
if (AnyYearLast .and. B>0) then
  YearLast_B = YearLast(B)
else
  YearLast_B = 9999
endif
do n=1,2
  if (m/=n .and. (m==1 .or. m==2 .or. focal>5))  cycle
  if (m==0 .and. n==2)  cycle  ! e.g. for focal=2 (FS)
  ALRm(n) = calc_ALRm(kA,kB, n,fcl)
enddo

ALR = MAXVAL(ALRm)
if (ALR < -HUGE(0.0D0) .or. ALR/=ALR)   ALR = impossible


! if (A==-176 .and. kA==1 .and. B==1820 .and. focal==4) then
  ! open (unit=42,file="log.txt",status="unknown", position="append")
  ! write (42, *) ""
  ! write(42,'("ALR CAU ", 2i5," ,",2i5,": ", f9.3, i3)') kA, A, kB, B, ALR, m
  ! write(42,'("ALRm ", 2f9.3)') ALRm
   ! do i=1,2
     ! write (42, '("BYLR ", i5, 100f10.4)') AB(i), 10**BYLR(:,i)
  ! enddo
  ! write (42, *) ""   
  ! ! do y=1,nYears !-1
    ! ! if (ALL(ALRtmp(:,y) == LOG10(zero)))  cycle
    ! ! write(42, '("ALRtmp ", i6, 100f10.4)') BYzero + y, 10**ALRtmp(:, y) !, COUNT(ALRtmp(:,y) >0)
  ! ! enddo
  ! ! write (42, *) ""   
  ! close(42)
  ! ! stop
! endif

contains
  function calc_ALRm(kA,kB, m,focal)  
    integer, intent(IN) :: kA,kB, m,focal
    double precision :: calc_ALRm

    do y=1,nYears  ! B
      if (BYLR(y,2) < -HUGE(0.0D0)) cycle
      do x=1, nYears  ! A 
        if (BYLR(x,1) < -HUGE(0.0D0)) cycle
        if (focal==1 .and. x > YearLast_B)  cycle  ! maybe-BY for A after B's last year of reproduction
        if (focal==4 .and. x >(YearLast_B + MaxAgePO))  cycle
        if ((x-y) < -MaxAgePO .or. (x-y) > nYears)  cycle
        if (BYLR(y,2) < -HUGE(0.0D0)) cycle
        if (focal==-1) then  ! A==B
          if (x==y)  ALRtmp(x,y) = BYLR(x,1) + BYLR(y,2)
        else if (focal<=5) then
          ALRtmp(x,y) = BYLR(x,1) + BYLR(y,2) + getAP(x-y, focal, kB, m, LOG10(zero))  
        else
          ALRtmp(x,y) = BYLR(x,1) + BYLR(y,2) + getAP(x-y, focal, m, kA, LOG10(zero))
        endif                    
      enddo
    enddo
    calc_ALRm = LOG10(SUM(10**ALRtmp))  ! sum across age differences
    
  end function calc_ALRm

end subroutine CalcAgeLR

! ######################################################################

subroutine CalcALRmerge(SA, SB, k, ALR)  ! change in ALR when merging
use Global
implicit none

integer, intent(IN) :: SA, SB, k
double precision, intent(OUT) :: ALR
double precision :: ALRj
integer :: i, j

ALR = 0D0
do i = 1, nS(SA,k)
  if (BY(SibID(i,SA,k))<0) cycle  ! TODO
  do j=1, nS(SB, k)
    ALRj = getAP(AgeDiff( SibID(i,SA,k), SibID(j,SB,k)), 3, 0, k, Impossible)
    if (ALRj == Impossible) then
      ALR = impossible
      return
    else
      ALR = ALR + ALRj
    endif
  enddo
enddo

ALR = ALR / (ns(SA,k) * ns(SB,k))   ! else not comparable across sibship sizes

end subroutine CalcALRmerge

! ######################################################################
  
subroutine calcALR_addsib(A,SB,k,focal, dALR)   ! change in ALR when adding A to SB
use Global
implicit none
  integer, intent(IN) :: A, SB, k, focal
  double precision, intent(OUT) :: dALR
  integer :: j
  double precision :: ALRj
  
  dALR = 0D0
  do j=1,ns(SB,k)
    call CalcAgeLR(A,3, SibID(j,SB,k),3, k, focal, .TRUE., ALRj)
    if (ALRj == Impossible) then
      dALR = impossible
      return
    else
      dALR = dALR + ALRj
    endif
  enddo
  
  dALR = dALR / ns(SB,k)

end subroutine calcALR_addsib

! #####################################################################

! #####################################################################

subroutine BestRel(LLIN, focal, X, dLL)
use Global
implicit none
! return which relationship is most likely, by threshold TA
! assuming order PO,FS,HS,GG,FAU,HAU,U in LL vector

double precision, intent(IN) :: LLIN(7)
integer, intent(IN) :: focal
integer, intent(OUT) :: X
double precision, intent(OUT) :: dLL   ! diff best vs next best
integer :: i,j, Y(7)
double precision :: LL(7)
logical:: Maybe(6)

X=0
dLL = 0D0
LL = LLIN

if (ALL(LL(1:6) > 0)) then
  X = 8
  return
endif

if (focal==3 .and. LL(3)<0 .and. LL(2)<0) then   ! want sib vs non-sib
  if (focal==3 .and. LL(2) - LL(3) < TA) then   ! FS less likely, or sliiiightly more likely
    LL(2) = 333D0  
  else if (focal /= 2 .and. LL(2)>=LL(3)) then
    LL(3) = 333D0
  endif
endif

if ((LL(7) - MAXVAL(LL(1:6), MASK=LL(1:6)<0)) > TA) then  
  X = 7  ! unrelated
else
  maybe = .TRUE.
  do i=1,6
    if (LL(i)>0) then
      maybe(i) = .FALSE.
    else
      do j=1,7
        if (i==j) cycle
        if (LL(j)>0) cycle
        if ((LL(i) - LL(j)) < TA) then
          maybe(i) = .FALSE.   ! i has no longer highest LL
        endif
      enddo
    endif
  enddo
  if (COUNT(maybe)==0) then
    X = 8  ! unclear
  else if (COUNT(maybe)==1) then
    X = MAXLOC(LL(1:6), MASK=maybe, DIM=1)
  endif       
endif

Y = (/(i, i=1,7, 1)/)
if (X<8 .and. X>0) then
  dLL = LL(X) - MAXVAL(LL, MASK=(LL<0 .and. Y/=X))
endif  
      
end subroutine BestRel

! #####################################################################

subroutine BestRel2(LLIN, X, dLL)
use Global
implicit none
! as BestRel, but no threshold, and consider all 1st & 2nd degree rel

double precision, intent(IN) :: LLIN(7)
integer, intent(OUT) :: X
double precision, intent(OUT) :: dLL !(2)   ! diff best vs next best
integer :: i,j, Y(7)                    
double precision :: LL(7)
logical:: Maybe(6)

X = 0
dLL = 0D0
LL = LLIN

if (MAXVAL(LL(1:6), MASK=LL(1:6)<0) - LL(7) < TA .or. &
  ALL(LL(1:6) > 0)) then  
  X = 7  ! unrelated 
else
  maybe = .TRUE.
  do i=1,6
    if (LL(i)>0) then
      maybe(i) = .FALSE.
    else
      do j=1,7
        if (i==j) cycle
        if (LL(j)>0) cycle
        if ((LL(i) - LL(j)) < 0.01) then
          maybe(i) = .FALSE.   ! i has no longer highest LL
        endif
      enddo
    endif
  enddo
  if (COUNT(maybe)==0) then
    if (ABS(MaxLL(LL(3:5))-MaxLL(LL))<0.01 .and. COUNT(LL(3:5)<0)>1) then
      X = 9  ! any 2nd degree relative
    else
      X = 8  ! unclear
    endif
  else if (COUNT(maybe)==1) then
    X = MAXLOC(LL(1:6), MASK = Maybe, DIM=1)
  else if (COUNT(maybe)>1) then
   if (ABS(MaxLL(LL(3:5))-MaxLL(LL))<0.01 .and. COUNT(LL(3:5)<0)>1) then
      X = 9  ! any 2nd degree relative
    else
      X = 8
    endif
  endif       
endif

Y = (/(i, i=1,7, 1)/)
if (count(LL < 0) < 2) then
  dLL = -777D0   ! NOT positive !!
else if (X<8) then
  dLL = LL(X) - MAXVAL(LL, MASK=(LL<0 .and. Y/=X))
!  dLL(2) = LL(X) - MaxLL(LL(6:7))
else if (X==9) then
  dLL = MaxLL(LL(3:5)) - MaxLL(LL((/1,2,6,7/)))
!  dLL(2) = MaxLL(LL(3:5)) - MaxLL(LL(6:7))
endif  

end subroutine BestRel2

! #####################################################################

subroutine UpdateAllProbs
use Global
use CalcLik
implicit none

integer :: i, k, s, x, y, r, n=30, BYrankI(nInd), BYrankC(nInd/2,2)
double precision :: Lind_IN(nInd)
double precision, allocatable :: XPr_IN(:,:,:,:,:)

if (DoSibs) then
  allocate(XPr_IN(3,3,nSnp,nInd/2,2))          

  do k=1,2
    do s=1,nC(k)
      if (nC(k)==0)  cycle
      if (ALL(GpID(:,s,k)==0) .and. ALL(SibID(:,s,k)==0)) then
        call Erstop("Empty sibship!", .TRUE.)
      endif
    enddo
  enddo
endif

! individuals ~~~~
call getRank_i(BYrankI)

do x=1, nInd
  i = BYRankI(x)
  call CalcLind(i)   
enddo

do r=1,2
  do x=1, nInd
    i = BYRankI(x)
    call setEstBY(i, Sex(i))     
  enddo

  do x=nInd, 1, -1
    i = BYRankI(x)
    call setEstBY(i, Sex(i))     
  enddo
enddo

! sibships ~~~
if (DoSibs) then
  do k=1,2
    call getBYrank_c(k, BYrankC(:,k))
  enddo

  do y=1,n  
    do r=1,n
      XPr_IN = XPr
      do x=1, MAXVAL(nC)
        do k=1,2
          if (x > nC(k))  cycle
          s = BYrankC(x,k)
          call CalcCLL(s,k)
        enddo
      enddo
      if (all(abs(XPr_IN - XPr) < 0.1))  exit
    enddo
    Lind_IN = Lind
    do x=1, nInd
      i = BYRankI(x)
      call CalcLind(i)
    enddo
    if (all(abs(Lind_IN - Lind) < 0.1))   exit
  enddo

!  do i=1, nInd
!    call CalcFSLik(i)
!  enddo

  do x=1, MAXVAL(nC)
    do k=1,2
      if (x > nC(k))  cycle
      s = BYrankC(x,k)
      call setEstBY(-s, k)
    enddo
  enddo

  deallocate(XPr_IN)
endif

end subroutine UpdateAllProbs

! #####################################################################

subroutine CalcCLL(s,k) 
use Global
use CalcLik
implicit none
! returns XPr: likelihood;  DumP: probability, scaled  (no age prior.),
! split into 1: sibs only 2: gp effect only, 3: all

integer, intent(IN) :: s, k ! S: sibship number, k: mat(1),pat(2),unk(3)
integer :: Sibs(ns(s,k)), UseEE(ns(s,k)), MatePar(ns(s,k)),  cat, catG, &
  IsInbr(ns(s,k)), AncR(2,mxA), FSX, &
  l, x, i, Ei, r, y, z, g, Ri, v, e
double precision :: PrL(nSnp), PrY(3), PrYp(3,ns(s,k)), PrGG(3,2),&
 PrZ(3),PrXZ(3,3,2), PrE(3), PrEE(3, ns(s,k)), LPrX(3,2), PrRF(3,3,ns(s,k)), PrF(3,3)
logical :: ParIsClone(ns(s,k)), DoRsibs(maxSibSize)
integer, allocatable :: HasInbr(:,:)

allocate(HasInbr(ns(s,k), ns(s,k)))
                  
if (ALL(GpID(:,s,k)==0) .and. ALL(SibID(:,s,k)==0)) then
  CLL(s,k) = 0D0
  XPr(1,:,:,s,k) = 1D0 
  do l=1,nSnp
    XPr(2,:,l,s,k) = AHWE(:,l)
    XPr(3,:,l,s,k) = AHWE(:,l)
    DumP(:,l,s,k) = AHWE(:,l)
  enddo
  return 
endif

Sibs = SibID(1:ns(s,k), s, k)
UseEE = 0
MatePar = 0            
call FindEE(Sibs, ns(s,k), 0, k, UseEE, MatePar)   ! may shuffle sibs

!================= 
cat = 0
catG = 0
IsInbr = 0
HasInbr = 0
AncR = 0
ParIsClone = .FALSE.
do r=1,nS(s,k)
  Ri = Sibs(r) 
  if (nFS(Ri) > ns(s,k)) then
    print *, ""
    print *, "nFS > nS! ", k, s, " sibs: ", SibID(1:5, s, k)
    print *, Ri, "  FS: ", FSID(1:nFS(Ri), Ri)
    call ErStop("something wrong with sibship cluster", .TRUE.)
  endif
  
  if (Parent(Ri, 3-k)==0) cycle
  if (Parent(Ri, 3-k)==GpID(3-k,s,k) .and. nFS(Ri)/=0) then  
    cat = Ri
    UseEE(r) = 0
  endif
  do v=1, nS(s,k)
    if (r==v) cycle
    if (nFS(Sibs(v))==0) cycle
    do i=1, nFS(Sibs(v))
      if (Parent(Ri, 3-k) == FSID(i, Sibs(v))) then
        IsInbr(r) = FSID(i, Sibs(v))
        HasInbr(v,i) = r !-1
      endif
    enddo
  enddo
  if (IsInbr(r)/=0) cycle
  call GetAncest(Ri,k,AncR)
  if (AncR(k, 5-k) == -s) then 
    IsInbr(r) = Parent(Ri, 3-k)   ! via dummy
  endif
  if (hermaphrodites/=0 .and. DumClone(s,k)/=0) then
    if (DumClone(s,k) == -Parent(Ri,3-k))  ParIsClone(r) = .TRUE.
  endif
enddo  

FSX = 0  
if (ALL(GpID(:,s,k)<0) .and. cat==0) then  ! check if sibship par inbred
  if (GPID(1,s,k) == GPID(1, -GPID(2,s,k),2)) then
    catG = 2
  else if (GPID(2,s,k) == GPID(2, -GPID(1,s,k),1)) then
    catG = 1
  else 
    do i=1, ns(-GpID(1,s,k), 1)
      if (Parent(SibID(i, -GpID(1,s,k), 1), 2) == GpID(2,s,k)) then  ! FS of dummy par
        catG = 3
        if (nFS(SibID(i, -GpID(1,s,k), 1))/=0) then
          FSX = SibID(i, -GpID(1,s,k), 1)
        endif
      endif
    enddo 
  endif
endif

call ChkTooManySibs(Sibs, ns(s,k), k, DoRsibs)  ! prevent numerical issues w huge sibships

PrL = 0D0       
do l=1,nSnp
  do g=1,2   !grandparents
    if (g/=k .and. cat>0) then
      call ParProb(l, GpID(g,s,k), g, -1,0, PrGG(:,g))
    else if (catG==g) then
      PrGG(:,g) = XPr(3,:,l, -GpID(g,s,k),g)
    else if (catG==3) then
      call ParProb(l, GpID(g,s,k), g, FSX,-1, PrGG(:,g))
    else
      call ParProb(l, GpID(g,s,k), g, 0,0, PrGG(:,g))
    endif
  enddo
  
  do x=1,3  ! genotype dummy parent
    do z=1,3
      if (catG==k) then
        PrXZ(x,z,:) = SUM(AKA2P(x, z, :) * PrGG(:,k))
      else if (catG==3-k) then
        PrXZ(x,z,:) = SUM(AKA2P(x, z, :) * PrGG(z,3-k))
      else if (catG==3) then
        PrE = PrGG(:,k)
        do i=1, nFS(FSX)
          PrE = PrE * OKA2P(Genos(l,FSID(i,FSX)), :, z)
        enddo          
        PrXZ(x,z,:) = SUM(AKA2P(x, z, :) * PrGG(z,3-k) * PrE)
      else  ! catG==0
        PrXZ(x,z,:) = SUM(AKA2P(x, z, :) * PrGG(z,3-k) * PrGG(:,k))  ! GPs
      endif
    enddo
  enddo
  if (catG>0) then
    do v=1,2
      PrXZ(:,:,v) = PrXZ(:,:,v)/SUM(PrXZ(:,:,v)) 
    enddo
  endif
  do x=1,3
    XPr(2,x,l, s,k) = SUM(PrXZ(x,:,2))  ! GP 
  enddo
  LPrX = log10(SUM(PrXZ, DIM=2))  ! sum over z

  PrYp = 0D0
  PrRF = 1D0
  do r=1, nS(s,k)
    Ri = Sibs(r)
    if (NFS(Ri) == 0) cycle
    if (IsInbr(r) /= 0 .or. cat==Ri .or. ParIsClone(r)) then
      PrYp(:,r) = 1D0
    else if (DoRsibs(r)) then
      call ParProb(l, Parent(Ri, 3-k), 3-k, -1,0, PrYp(:,r))
    else
      call ParProb(l, Parent(Ri, 3-k), 3-k, Ri, -1, PrYp(:,r)) 
    endif
    if (all(HasInbr(r,:)==0)) PrRF(:,:,r) = FSLik(l,Ri)
  enddo

  do z=1,3
    if (z>1 .and. cat==0) cycle
  do x=1,3
    PrEE = 0D0
    do r=1, nS(s,k)
      Ri = Sibs(r)  ! array with IDs
      if (NFS(Ri) == 0) cycle  ! moved to its FS
      if (IsInbr(r) > 0) then
        cycle
      else if (IsInbr(r) < 0) then
        call ParProb(l, GpID(3-k,-Parent(Ri, 3-k),3-k), 3-k, 0,0, PrZ) 
        do y=1,3
          PrYp(y,r) = SUM(AKA2P(y,x,:) * PrZ)
        enddo
      else if (UseEE(r) /= 0) then
        call ParProb(l, MatePar(r), k, 0,0, PrZ)  !  GpID(k,-Parent(Ri, 3-k),3-k)
        do y=1,3
          do e=1,3
            PrE(e) = SUM(AKA2P(y,e,:) * PrEE(e,UseEE(r)) * PrZ)
          enddo
          PrYp(y,r) = SUM(PrE)
        enddo
        PrYp(:,r) = PrYp(:,r) / SUM(PrYp(:,r))
      endif      
      PrY = PrYp(:,r)
        
      if (Parent(Ri, 3-k)<0 .and. DoRsibs(r)) then 
        do y=1,3   ! parent 3-k 
          do v=1, nS(-Parent(Ri, 3-k), 3-k)
            Ei = SibID(v, -Parent(Ri, 3-k), 3-k)  
            if (NFS(Ei) == 0) cycle
            if (Parent(Ei, k) == -s) cycle
            call ParProb(l, Parent(Ei, k), k, Ei,-1, PrE)
            PrF = FSLik(l,Ei)
            PrE = PrE * PrF(:,y)
            if (.not. ALL(PrE==1D0))  PrY(y) = PrY(y) * SUM(PrE)
          enddo
        enddo
      endif
      
      if (.not. ALL(PrY==1D0)) then     
        if (cat==Ri) then
          PrXZ(x,z,1) = PrXZ(x,z,1) * PrY(z)
        else if (ParIsClone(r)) then
          PrXZ(x,z,1) = PrXZ(x,z,1) * PrY(x)          
        else if (cat/=0) then
          PrXZ(x,z,1) = PrXZ(x,z,1) * SUM(PrY)
        else
          PrXZ(x,:,1) = PrXZ(x,:,1) * SUM(PrY)
        endif
      endif
      if (cat==0 .and. .not. all(PrY==1D0)) then
        if (ParIsClone(r)) then
          LPrX(x,1) = LPrX(x,1) + log10(PrY(x))
        else
          LPrX(x,1) = LPrX(x,1) + log10(SUM(PrY))
        endif
      endif
     
      if (all(HasInbr(r,:)==0)) then
        PrY = PrY * PrRF(:,x,r)
      else
        do i=1, nFS(Ri)  ! default: nFS = 1    
          if (HasInbr(r,i)==0) then
            PrY = PrY * OKA2P(Genos(l,FSID(i,Ri)), x, :)
          else
            do y=1,3
              do e=1,3
                PrE(e) = AKA2P(e, x, y) * OcA(e,Genos(l,FSID(i,Ri)))
                do v=1, nS(s,k)
                  if (IsInbr(v)==FSID(i,Ri)) then  
                    PrE(e) = PrE(e) * OKA2P(Genos(l,Sibs(v)), x, e)
                  endif
                enddo
              enddo
              if (.not. ALL(PrE==1D0))  PrY(y) = PrY(y) * SUM(PrE)
            enddo
          endif
        enddo
      endif        
            
      if (.not. ALL(PrY==1D0)) then        
        if (cat==Ri) then
          PrXZ(x,z,2) = PrXZ(x,z,2) * PrY(z)
        else if (ParIsClone(r)) then
          PrXZ(x,z,2) = PrXZ(x,z,2) * PrY(x)          
        else if (cat/=0) then
          PrXZ(x,z,2) = PrXZ(x,z,2) * SUM(PrY)
        else
          PrXZ(x,:,2) = PrXZ(x,:,2) * SUM(PrY)
        endif
        if (cat==0 .and. .not. all(PrY==1D0)) then
          if (ParIsClone(r)) then
            LPrX(x,2) = LPrX(x,2) + log10(PrY(x))
          else
            LPrX(x,2) = LPrX(x,2) + log10(SUM(PrY))
          endif
        endif
      endif
      PrEE(:,r) = PrY   
    enddo ! r 
  enddo ! x
  enddo ! z (cat/=0 only)
  if (cat==0) then
    XPr(3,:,l, s,k) = 10**LPrX(:,2) / SUM(10**LPrX(:,1))
  else
    do x=1,3  ! account for GP, dumm offspr & connected sibships
      XPr(3,x,l, s,k) = SUM(PrXZ(x,:,2))/ SUM(PrXZ(:,:,1))    ! offspring + grandparents
    enddo
  endif
  WHERE (XPr(3,:,l,s,k) /= XPr(3,:,l,s,k)) XPr(3,:,l,s,k) = 0D0  ! when Err=0 
  do x=1,3
    DumP(x,l, s,k) = XPr(3,x,l, s,k)/ SUM(XPr(3,:,l, s,k))
    XPr(1,x,l, s,k) = XPr(3,x,l, s,k) / XPr(2,x,l, s,k)   ! offspring only 
  enddo 
  PrL(l) = LOG10(SUM(XPr(3,:,l, s,k)))
enddo
CLL(s,k) = SUM(PrL) 
 
WHERE (XPr(1,:,:,s,k) /= XPr(1,:,:,s,k)) XPr(1,:,:,s,k) = 0D0  ! 0/0 when MAF=0
WHERE (DumP(:,:,s,k) /= DumP(:,:,s,k)) DumP(:,:,s,k) = 0D0

 
if (CLL(s,k)> .001 .or. CLL(s,k)/=CLL(s,k)) then   !.or. CLL(s,k)< -HUGE(1D0) 
  call Rprint("Problem: ", (/k, s, ns(s,k)/), (/0.0D0/), "INT")   
  if (hermaphrodites/=0)  print *, "selfed: ", any(ParIsClone)
  call Erstop("Invalid sibship LL - try increasing Err", .FALSE.)
endif

! if (k==2 .and. any(SibID(:,s,k)==840) .and. any(SibID(:,s,k)==841)) then
! print *, ""
 ! print *, "CLL: ", k, s, ", ", GpID(:,s,k), ", ", cat, catG, CLL(s,k)
  ! do i=1,ns(s,k)
    ! write(*,'(3i5, " - ", i3)') SibID(i,s,k), Parent(SibID(i,s,k),:), &
      ! nFS(SibID(i,s,k))
  ! enddo
  ! print *, "."
  ! do i=1,ns(s,k)
    ! write(*,'(i3, i5, " - ", i5, " . ", 8i2 )') i, Sibs(i), UseEE(i), IsInbr(i), HasInbr(i, 1:5) 
  ! enddo
  ! print *, ""
! endif
deallocate(HasInbr)

end subroutine CalcCLL

! #####################################################################

subroutine ChkTooManySibs(Sibs, n, k, DoRsibs)   ! Not s: called after reshuffle by FindEE
use Global
implicit none

integer, intent(IN) :: n
integer, intent(IN) :: sibs(n), k
logical, intent(OUT) :: DoRsibs(maxSibSize)
integer :: r, i

! prevent numerical issues when sibships are very large
DoRsibs = .FALSE.
do r=1,n
  i = Sibs(r) 
  if (nFS(i)==0)  cycle
  if (Parent(i,3-k) >=0) cycle
  if (ns(-parent(i,3-k),3-k) >50 .and. nFS(i) < ns(-parent(i,3-k),3-k)/5) then  ! What thresholds??
    DoRsibs(r) = .FALSE.
  else
    DoRSibs(r) = .TRUE.
  endif
enddo

end subroutine ChkTooManySibs

! #####################################################################

subroutine ParProb(l, i, k, A, B, prob)  
use Global
use CalcLik
implicit none

integer, intent(IN) :: l, i, k, A,B
double precision, intent(OUT) :: prob(3)
integer :: x,j, AB(2), A1, parA
double precision :: PrP(3, 2), PrY(3), PrAF(3,3)
logical :: AllIN

prob = AHWE(:, l)    
A1 = 0

if (i == 0) then  ! no parent
  if (A==-4 .or. B==-4 .or. B==-5) then
    prob = 1D0
  else
    prob = AHWE(:, l)
  endif

else if (i > 0) then  ! real parent
  if (A==-4 .or. B==-4 .or. B==-5) then
    prob = OcA(:,Genos(l,i))
  else
    prob = LindX(:,l,i)  ! =AHWE if Lind(i) not yet calculated  ! unscaled; scaled below.
  endif

else if (i < 0) then  ! dummy parent
  if (A==0) then   ! probability
    prob = DumP(:,l, -i,k)    
  else if (A == -1) then  ! grandparent contribution only
    prob = XPr(2,:,l, -i, k)
  else if (A==-4) then  ! offspring contribution only
    prob = XPr(1,:,l, -i, k)
  else if (A<0) then  ! shouldn't happen
    call Erstop("Invalid call to ParProb!", .TRUE.)
  
  else if (A>0) then   ! exclude indiv A from calc & standardise
    if (Parent(A,k)/=i .and. B<=0) then
      prob = DumP(:,l, -i,k)
    else if (ns(-i,k)<=1 .and. B> -4) then
      prob = XPr(2,:,l, -i, k)  ! grandparent contribution only
    else
     
      AB = (/ A, B /)
      A1 = FSID(maxSibSize+1, A)
      do j=1,2
        if (j==1 .and. (Parent(A,k)/=i)) cycle
        if (j==2) then
          if(B<=0) cycle
          if (Parent(B,k)/=i)  cycle
        endif
        if (Parent(AB(j), 3-k)==0) then  
          PrP(:,j) = AHWE(:,l)
        else if (Parent(AB(j), 3-k)>0) then
          PrP(:,j) = LindX(:,l, Parent(AB(j), 3-k))
          PrP(:,j) = PrP(:,j)/SUM(PrP(:,j))
        else if (Parent(AB(j), 3-k)<0) then  
          PrP(:,j) = DumP(:,l, -Parent(AB(j),3-k), 3-k)
        endif
      enddo
      
      AllIN = .FALSE.
      if (B==-1) then
        parA = Parent(A, 3-k)
        if (ns(-i,k) <= 1) then
          AllIN = .TRUE.
        else if (ParA/=0 .and. all(Parent(SibID(1:ns(-i,k),-i,k),3-k) == ParA)) then 
          AllIN = .TRUE.
        endif
      endif
   
      prob = 0D0
      PrAF = 1D0
      if ((B==-1 .and. .not. AllIN) .or. (B==-5 .and. nFS(A1) /= ns(-i,k))) then
        PrAF = FSLik(l,A1)
      endif
      do x = 1, 3
        if (B>=0) then    ! .or. (B==-1 .and. nFS(A1)<=1)
          prob(x) = XPr(3,x,l, -i, k)
          if (Parent(A,k)==i) then
            PrY = OKA2P(Genos(l,A),x,:) * PrP(:,1)
            if (SUM(PrY) > 0D0)  prob(x) = prob(x) / SUM(PrY)
          endif
          if (B>0) then
            if (Parent(B,k)==i) then
              PrY = OKA2P(Genos(l,B), x, :) * PrP(:,2)
              if (SUM(PrY) > 0D0)  prob(x) = prob(x) / SUM(PrY)  
            endif
          endif
        else if (B==-1) then  ! exclude all FS of A
          if (ALL(Xpr(3,:,l,-i,k)==1D0)) then   ! at initiate 
            prob(x) = AHWE(x, l)
          else if (AllIN) then
            prob(x) = XPr(2,x,l, -i, k)   ! GP only
          else
            PrY = PrP(:,1) * PrAF(:,x)
            if (SUM(PrY) > 0D0)  prob(x) = XPr(3,x,l, -i, k) / SUM(PrY)
          endif     
        else if (B==-4) then ! exclude both GPs & A
          if (ns(-i,k)==1 ) then  !  
            prob(x) = 1D0   ! AHWE(:,l)  !
          else
            PrY = OKA2P(Genos(l,A),x,:)* PrP(:,1)
            if (SUM(PrY) > 0D0)  prob(x) = XPr(1,x,l,-i,k) /SUM(PrY)
          endif     
        else if (B==-5) then ! exclude both GPs & A & FS of A
          if (nFS(A1) == ns(-i,k)) then  ! ns(-i,k)==1 
            prob(x) = 1D0   ! AHWE(:,l)  !
          else
            PrY = PrP(:,1) * PrAF(:,x)
            if (SUM(PrY) > 0D0)  prob(x) = XPr(1,x,l,-i,k) / SUM(PrY)   ! *AHWE(x,l)
          endif
        else
          call Erstop("ParProb: invalid B!", .TRUE.)
        endif
      enddo
    endif
  endif
endif
if (SUM(prob)>0D0 .and. .not. ALL(prob==1D0)) then
  prob = prob/SUM(prob)
endif
 
if (ANY(prob< 0D0) .or. ANY(prob/=prob) .or. ANY(prob>1.01D0)) then
  call Rprint( "Indiv, k,A,B: ", (/i, k, A, B/), (/0.0D0/), "INT") 
  if (A1/=0)  call Rprint("A1, nFS, ns, par: ", (/A1, nFS(A1), ns(-i,k), Parent(A,:)/), (/0.0D0/), "INT") 
  print *, "AllIN: ", AllIN
  call Rprint("prob: ", (/0/), prob, "DBL")  
  call Erstop("Invalid ParProb!", .TRUE.)
endif

end subroutine ParProb

! #####################################################################

subroutine OffProb(l,i,k, prob)
use Global
implicit none

integer, intent(IN) :: l, i, k
double precision, intent(OUT) :: prob(3)

if (i > 0) then
  prob = OcA(:,Genos(l,i))
else if (i < 0) then
  prob = XPr(1,:,l,-i,k)
else
  prob = 1D0
endif

end subroutine OffProb

! #####################################################################

subroutine Connected(A, kA, B, kB, Con)
use Global
implicit none

integer, intent(IN) :: A, kA, B, kB
logical, intent(OUT) :: Con
integer :: i, j, m, nA, nB, AA(maxSibsize), BB(maxSibsize), n

Con = .FALSE.
if (A==0 .or. B==0)  return

AA = 0
BB = 0
if (A>0) then
  nA = 1
  AA(1) = A
else
  nA = nS(-A,kA)
  AA(1:nA) = SibID(1:nA, -A, kA)
endif

if (B>0) then
  nB = 1
  BB(1) = B
else
  nB = nS(-B,kB)
  BB(1:nB) = SibID(1:nB, -B, kB)
endif

do j=1, nB
  do i=1, nA
    do m=1,2  
      if (Parent(AA(i), m) < 0) then
        if (Parent(AA(i),m) == Parent(BB(j),m)) then
          Con = .TRUE.
          return
        else if(ANY(GpID(:,-Parent(AA(i), m),m) == BB(j))) then
          Con = .FALSE.  ! TODO: update Uclust(). already conditioned on?
!          return
        else if (A<0 .and. m==kA) then
          if(ANY(GpID(:,-Parent(AA(i), m),m) < 0)) then
            do n=1,2
              if (GpID(n,-Parent(AA(i), m),m) == Parent(BB(j),n) .and. &
                Parent(BB(j),n)<0) then
                Con = .TRUE.
                return
              endif 
            enddo
          endif
        endif
      endif
      if (Parent(BB(j),m)<0) then
        if (ANY(GpID(:,-Parent(BB(j),m),m) == AA(i))) then
          Con = .FALSE.  ! TODO
!          return
        else if (B<0 .and. m==kB) then
          if(ANY(GpID(:,-Parent(BB(j),m),m) < 0)) then
            do n=1,2
              if (GpID(n,-Parent(BB(j), m),m) == Parent(AA(i),n) .and. &
                Parent(AA(i),n)<0) then
                Con = .TRUE.
                return
              endif 
            enddo
          endif
        endif
      endif
    enddo
  enddo
enddo

end subroutine Connected

! #####################################################################

subroutine GetAncest(A, kIN, Anc)
use Global
implicit none

integer, intent(IN) :: A, kIN
integer, intent(OUT) :: Anc(2, mxA)  ! 32 = 5 generations
integer :: m, j, i, k, Par(2)

! Anc(1,:)  female ancestors
! Anc(2,:)  male ancestors
! Anc(:,2)  parents
! Anc(:,3)  mat. gp
! Anc(:,4)  pat. gp

Anc = 0
if (A==0) return

k = kIN
if (A > 0) then  ! real indiv
  if (kIN < 1 .or. kIN > 2)  k = 1
  Anc(k,1) = A
else !if (A < 0) then  ! dummy indiv
  if (kIN < 1 .or. kIN > 2) then
   call Erstop("getAncest: k must be 1 or 2 if A<0", .TRUE.)
  else
    Anc(k,2) = A  !! 
  endif
endif

if (A<0 .and. kIN/=1 .and. kIN/=2) then
  print *, 'getAncest: invalid kIN!'
  stop
endif

Par = getPar(A,k)
if (ALL(Par == 0))  return
if (A > 0)  Anc(:, 2) = Par

do j = 2, mxA/2  
  do m = 1, 2
    i = 2 * (j-1) + m
    Anc(:,i) = getPar(Anc(m,j), m)
  enddo
  if (j==2 .and. ALL(Anc(:, 3:4) == 0))  return
  if (j==4 .and. ALL(Anc(:, 5:8) == 0))  return
  if (j==8 .and. ALL(Anc(:, 9:16) == 0))  return
  if (j==16 .and. ALL(Anc(:, 17:32) == 0))  return
enddo

if ((A>0 .and. ANY(Anc(:, 2:mxA)==A)) .or. (A<0 .and. ANY(Anc(k,3:mxA)==A))) then
  call Rprint( "Female ancestors: ", Anc(1,1:8), (/0.0D0/), "INT")
  call Rprint( "Male ancestors: ", Anc(2,1:8), (/0.0D0/), "INT")
  call Erstop("An individual is its own ancestor! Need more birth years or better SNP data", .FALSE.)
endif

end subroutine GetAncest

! #####################################################################

subroutine ChkAncest(A, kA, B, kB, OK)  ! check that B is not an ancestor of A
use Global
implicit none

integer, intent(IN) :: A, kA, B, kB
logical, intent(OUT) :: OK
integer :: AncA(2, mxA), j

OK = .TRUE.
if (A==0 .or. B==0)  return
call GetAncest(A, kA, AncA)

if (B > 0) then
  if (ANY(AncA == B))  OK = .FALSE.  
else if (kB == 1 .or. kB==2) then
  if (ANY(AncA(kB,:) == B))  OK = .FALSE.
  if (hermaphrodites/=0 .and. DumClone(-B,kB)/=0) then
    if (ANY(AncA(3-kB,:) == -DumClone(-B,kB)))  OK = .FALSE.
  endif
else
  call ErStop("ChkAncest: kB must be 1 or 2 if B<0", .TRUE.)
endif

if (OK .and. B < 0 .and. A<0) then   ! check 1 extra generation
  if (ns(-B,kB)==0)  return
  do j=1, ns(-B,kB)
    if (ANY(AncA == SibID(j,-B,kB))) then
      OK = .FALSE.
      exit
    endif
  enddo
endif

end subroutine ChkAncest

! #####################################################################

subroutine CalcParentLLR(LLR_parent, LLR_GP)
! Calc parental LLR (vs next most likely relationship)
use Global
implicit none

double precision, intent(OUT) :: LLR_Parent(nInd,3), LLR_GP(3, nInd/2, 2)
integer :: i, CurPar(2), k, m, nonG(6), CurGP(2), s, g, t
double precision :: LLtmp(2,2,2), LLg(7), LLa(7)
logical :: AllSibsSelfed(2), NoGP, FSM

if (quiet < 1)  call printt('updating all probs ... ')
call UpdateAllProbs()
LLR_parent = missing
LLR_GP = missing

if (quiet<1)  call Rprint("Calculating parental LLR ... ",(/0/), (/0.0D0/), "NON")
t = 1
do i = 1, nInd
  if (MODULO(i,100)==0) call rchkusr()
  if (quiet==-1 .and. any(viginti==i)) call print_progress(i,t)
  if (skip(i))  cycle
  if (Parent(i,1)==0 .and. Parent(i,2)==0) cycle

  CurPar = Parent(i,:) 
  do k=1,2  ! remove i from sibgroup
    call setParTmp(i, Sex(i), 0, k)
  enddo
  
  AllSibsSelfed = .FALSE.
  NoGP = .FALSE.
  if (hermaphrodites/=0 .and. all(curPar <0)) then
    if (DumClone(-curPar(1),1) == -curpar(2)) then
      do k=1,2
        if (ns(-curPar(k),k)==0) then
          AllSibsSelfed(k) = .TRUE.
        else if (all(SelfedIndiv(SibID(1:ns(-curPar(k),k),-curPar(k),k)))) then
          AllSibsSelfed(k) = .TRUE.
        endif
      enddo
      if (all(GpID(:,-curpar(1),1)==0))  NoGP = .TRUE.
    endif
  endif
  
  LLtmp = missing
  do m=1,2  ! m=1: no opp. sex parent;  m=2: with opp. sex parent
    if (m==2 .and. (CurPar(1)==0 .or. CurPar(2)==0)) cycle
    do k=1,2  ! mother, father
      if (CurPar(k)==0) cycle
      if (k==2 .and. Complx==0 .and. .not. any(curPar==0))  cycle
      if (AllSibsSelfed(k) .and. .not. (k==1 .and. AllSibsSelfed(2) .and. .not. NoGP))  cycle
      if (m==2) then  ! temp. assign parent 3-k    .and. .not. AllSibsSelfed(3-k)
        call setParTmp(i, Sex(i), curPar(3-k), 3-k)
      endif
            
      if (CurPar(k) > 0) then
        call CheckPair(i, CurPar(k), k, 7, LLg, LLa) 
        LLtmp(1,k,m) = LLg(1)
        LLtmp(2,k,m) = MaxLL(LLg(2:7))
      else if (CurPar(k) < 0) then 
        call CheckAdd(i, -CurPar(k), k, 7, LLg, LLa)
        if (m==1) LLg(2) = 333D0   ! FS does not count here
        LLtmp(1,k,m) =  MaxLL(LLg(2:3))
        LLtmp(2,k,m) =  MaxLL((/LLg(1), LLg(4:7)/))
      endif
            
      if (m==2 .and. k==1) then  
        call setParTmp(i, Sex(i), 0, 3-k)
      endif

!      if (any(LLtmp(:,k,m) >0) .and. Complx>0) then
 !       write (*, '("parent LLtmp > 0: ", i5, a10, 2i3, 2f8.2 )')  i, trim(ID(i)), k, m, LLtmp(:,k,m)
 !     endif 
   
    ! if (i==994) then 
      ! write(*,'("ParLR ", a12, 2i5, 2f8.1, "; ", 7f8.1)') ID(i), m, k, LLtmp(:,k,m), LLg
    ! endif      
     
    enddo
  enddo 
  
  do k = 1,2   ! restore
    if (Parent(i,k)/=curPar(k))  call setParTmp(i, Sex(i), CurPar(k), k) 
  enddo  
  ! curPar(1) mostly restored in round m=2, k=2
  
  if (all(AllSibsSelfed) .and. NoGP) then
    call IsSelfed(i, .FALSE., LLR_parent(i,3)) 
  else
    call ParLLtoLR(LLtmp, LLR_parent(i,:), AllSibsSelfed)
  endif
enddo

!parents of dummies (Sibship GPs)
nonG = (/1,2,3,5,6,7/)
do k = 1,2  
  if (nC(k)==0)  cycle
  if (quiet==-1 .and. sum(nC)>20) then
    call Rprint("Dummies...", (/k/), (/0D0/), "INT")
  endif
  do s=1, nC(k)
    if (MODULO(s,10)==0) call rchkusr()
    CurGP = GpID(:, s, k)
    do g=1,2
      call setParTmp(-s, k, 0, g)
    enddo      
    LLtmp = missing
    do m=1,2
      if (m==2 .and. (CurGP(1)==0 .or. CurGP(2)==0)) cycle
      do g=1,2
        if (CurGP(g) == 0) cycle
        if (g==2 .and. Complx==0 .and. .not. any(curGP==0))  cycle
        if (m==2) then  ! temp. assign GP 3-g
          call setParTmp(-s, k, CurGP(3-g), 3-g)
        endif
        
        if (curGP(g) > 0) then
          call checkAdd(CurGP(g),s,k, 7, LLg, LLa)  ! B=GP + CurGP(m)_7
        else if (curGP(g) < 0) then
          call checkMerge(s, -CurGP(g), k, g, 4, LLg, LLa, FSM)   !TODO: use FSM?
          if (m==1) then
            call PairUA(-s, CurGP(g), k, g, LLg(4))  
          endif
        endif
        LLtmp(1,g,m) = LLg(4)
        LLtmp(2,g,m) = MaxLL(LLg(nonG))
        if (m==2 .and. g==1) then  ! reset to 0
          call setParTmp(-s, k, 0, 3-g)
        endif 
      enddo
    enddo
    
    do g=1,2
      call setParTmp(-s, k, CurGP(g), g)  ! restore
    enddo
    
    call ParLLtoLR(LLtmp, LLR_GP(:,s,k), (/.FALSE.,.FALSE./))
  enddo
enddo

end subroutine CalcParentLLR

! ######################################################################

subroutine ParLLtoLR(LLtmp, LLR, AllSibsSelfed)
use Global
implicit none

double precision, intent(IN) :: LLtmp(2,2,2)
double precision, intent(OUT) :: LLR(3)
logical, intent(IN) :: AllSibsSelfed(2)
integer :: k
double precision :: LLX(2)

! LLtmp dims: focal/not-focal ; dam/sire ; without/with other parent 
LLR = missing
if (Complx > 0) then
  do k=1,2  ! max with - max w/o 
    if (LLtmp(1,k,1) < 0) then
      LLR(k) = LLtmp(1,k,1) - LLtmp(2,k,1)
    else
      LLR(k) = LLtmp(1,k,1)  ! something wrong  / AllSibsSelfed  / monogamous
    endif
  enddo
endif

if (hermaphrodites/=0) then
  do k=1,2
    if (AllSibsSelfed(k)) then
      LLR(k) = LLR(3-k)
    endif
  enddo
endif

do k=1,2
  if (LLtmp(1,k,2) < 0) then
    LLX(k) = LLtmp(1,k,2) - MaxLL((/LLtmp(2,k,2), LLtmp(:,k,1)/))
  else if (Complx==0 .and. LLtmp(1,k,1) < 0) then   ! single grandparent
    LLX(k) = LLtmp(1,k,1) - LLtmp(2,k,1)
  else
    LLX(k) = LLtmp(1,k,2)
  endif
enddo
LLR(3) = MINVAL(LLX)

end subroutine ParLLtoLR

! ######################################################################

subroutine setPar(A, kA, P, kP)    ! Assigns parent P to A, incl. sex & age update
use Global
implicit none

integer, intent(IN) :: A, kA, P, kP
integer :: curPar(2)

if (A==0)  return

curPar = getPar(A, kA)
if (curPar(kP) /= P)  then
  call setParTmp(A, kA, P, kP)
  call SetEstBY(curPar(kP), kP)
endif

call UpdateLL(P, kP)
call UpdateLL(curPar(3-kP), 3-kP)
call UpdateLL(A,kA)

call SetEstBY(A, kA)
call SetEstBY(P, kP)

if (P > 0) then
  if (Sex(P) == 3)   Sex(P) = kP
!else if (P < 0) then
!  IsNewSibship(-P, kP) = .TRUE.
endif

if (A>0 .and. P/=0)  ToCheck(A) = .TRUE.
if (A<0 .and. P/=0) then
  if (ns(-A,kA) <= 3) then 
    ToCheck(SibID(1:ns(-A,kA),-A,kA)) = .TRUE.
  endif  
endif

if (hermaphrodites/=0) then
  if (A>0) then
    call CheckSelfed(A, sex(A))   ! sets Selfed(P,kP) if P<0
  else if (DumClone(-A,kA)/=0) then
    call setParTmp(-DumClone(-A, kA), 3-kA, P, kP)
    call SetEstBY(-DumClone(-A, kA), 3-kA) 
  endif
endif

end subroutine setPar

! ######################################################################

subroutine setParTmp(A, kA, P, kP)   ! Temporary assigns parent P to A
use Global
use CalcLik
implicit none

integer, intent(IN) :: A, kA, P, kP
integer :: curPar(2), nOffP, OffP(maxSibSize), sxOffP(maxSibSize), i                

curPar = getPar(A, kA)

if (A==0)  call Erstop("SetParTmp: A=0", .TRUE.)
if (kP/=1 .and. kP/=2)  call Erstop("SetParTmp: kP must be 1 or 2", .TRUE.)
if (A<0 .and. kA/=1 .and. kA/=2)  call Erstop("SetParTmp: kA must be 1 or 2 if A<0", .TRUE.)
if (P==0 .and. curPar(kP)==0)  return

if (P < 0) then
  if (-P > nC(kP)) then
    print *, ""
    print *, A, kA, P, kP, "; ", nC
    call Erstop("setParTmp: Sibship number out of bounds", .TRUE.)
  endif
endif

! remove old par
if (A > 0) then
  if (curPar(kP) > 0) then
    call RemoveFS(A)
    Parent(A,kP) = 0
    call CalcLind(A)

  else if (curPar(kP) < 0) then
    call RemoveSib(A, -curPar(kP), kP)   ! NOTE: doesn't drop sibship if singleton w/o GP
  endif
  
  if (P > 0 .and. curPar(3-kP)/=0) then  ! check for FS
    nOffP = 0
    if (ANY(Parent(:,kP) == P)) then   ! P already has some offspring
      call getOff(P, kP, .FALSE., nOffP, OffP, sxOffP)   ! currently only >0 FS considered
    endif
    Parent(A, kP) = P 
    call CalcLind(A)
    if (nOffP > 0) then
      do i=1, nOffP
        if (Parent(A,3-kP) == Parent(OffP(i), 3-kP)) then
          call MakeFS(A, OffP(i))
          call CalcLind(OffP(i))
!          call CalcFSLik(A)
!          call CalcFSLik(OffP(i))
        endif
      enddo
    endif

  else if (P > 0) then
    Parent(A, kP) = P
    call CalcLind(A)
    
  else if (P < 0) then
    call DoAdd(A, -P, kP)   
  endif

  if (curPar(3-kP) < 0) then
    call CalcCLL(-curPar(3-kP), 3-kP)
    if (ns(-curpar(3-kP),3-kP) > 0) then
      do i=1, ns(-curpar(3-kP),3-kP)
        call CalcLind(SibID(i, -curpar(3-kP),3-kP)) 
      enddo
    endif
    if (P < 0) then
      call CalcCLL(-P, kP)
      do i=1, ns(-P,kP)
        call CalcLind(SibID(i,-P,kP))
      enddo
    endif
  endif

else  ! A < 0
  GpID(kP, -A, kA) = P
  call CalcCLL(-A,kA)
  if (ns(-A,kA)>0) then
    do i=1, ns(-A,kA)
      call CalcLind(SibID(i,-A,kA))
    enddo
  endif
endif

end subroutine setParTmp

! #####################################################################

subroutine DoAdd(A, s, k)
use Global
use CalcLik
implicit none

integer, intent(IN) :: A, s, k
integer :: i, u

if (nS(s,k) +1 >= maxSibSize) then
  call Erstop("Reached Maximum Sibship Size (number of offspring per parent), please increase '--maxsibsize'", .FALSE.)
endif

Parent(A, k) = -s
if (.not. ANY(SibID(1:nS(s,k),s,k)==A)) then
  SibID(nS(s,k)+1, s, k) = A  ! add A to sibship
  nS(s,k) = nS(s,k) + 1
endif

if (ns(s,k)<=0 .or. SibID(1,s,k)<=0) then
  print *, ""
  print *, A, s, k, ", ", ns(s,k), SibID(1:10, s,k)
  call Erstop("DoAdd: invalid sibship", .TRUE.)
endif

do u=1, nS(s,k)  ! check for FS   
  i = SibID(u,s,k)
  if (i==0) then
   print *, ""
    print *, A, s, k, ", ", u, "; ",  ns(s,k), " --",  SibID(1:10, s,k)
    call Erstop("DoAdd: invalid sibship", .TRUE.)
  endif
  if (i==A .or. nFS(i)==0) cycle 
  if (Parent(A, 3-k)/=0 .and. Parent(A, 3-k)==Parent(i, 3-k)) then
    call MakeFS(A, i)
!    call CalcFSLik(A)
!    call CalcFSLik(i)
  endif
enddo

call calcCLL(s,k)
if (Parent(A,3-k) < 0) call CalcCLL(-Parent(A,3-k), 3-k)
do u=1,ns(s,k)
  call CalcLind(SibID(u,s,k))
enddo

end subroutine DoAdd

! #####################################################################

subroutine RemoveSib(A, s, k)  ! removes individual A from sibship s. Does NOT drop sibship if ns=0
use Global
use CalcLik
implicit none

integer, intent(IN) :: A, s, k
integer :: u, h

call RemoveFS(A)

do u=1,ns(s,k)
  if (SibID(u,s,k)==A) then
    if (u<ns(s,k)) then  ! shift sibs
      do h=u, nS(s, k)-1  ! drop HS
        SibID(h, s, k) = SibID(h+1, s, k)
      enddo
    endif
    SibID(nS(s,k), s, k) = 0
    nS(s,k) = nS(s,k) -1
    exit
  endif
enddo

Parent(A, k) = 0
call CalcCLL(s,k)
if (Parent(A,3-k) < 0) call CalcCLL(-Parent(A,3-k), 3-k)
if (ns(s,k) >0) then
  do u=1,ns(s,k)
    call CalcLind(SibID(u,s,k))
  enddo
endif
call CalcLind(A)

end subroutine RemoveSib

! ######################################################################

subroutine UpdateLL(A, k)   ! update LL of A & if A<0 of its mates
use Global
use CalcLik
implicit none

integer, intent(IN) :: A, k
integer :: Mates(maxSibSize), x, i, j, nOff, sxOff(maxSibSize), Off(maxSibSize)
double precision :: XPR_x(3, nSnp)
logical :: OK

if (A==0)  return

if (A >0) then
  call CalcLind(A)
  return
endif

if (ns(-A,k)==0)  return

Mates = 0
Mates(1:ns(-A,k)) = Parent(SibID(1:ns(-A,k), -A,k), 3-k)

do x=1,10
  OK = .TRUE.
  do i=1, ns(-A,k)
    if (nFS(SibID(i,-A,k)) == 0)  cycle
    if (Mates(i) < 0) then
      if (ANY(Mates == GpID(3-k, -Mates(i),3-k)) .and. GpID(3-k, -Mates(i),3-k) < 0) then
        call CalcCLL(-GpID(3-k, -Mates(i),3-k), 3-k)   ! Used by UseEE
      endif
      if (ns(-Mates(i),3-k) <= 50 .or. nFS(SibID(i,-A,k)) >= ns(-Mates(i),3-k)/5) then
        XPR_x = XPr(3,:,:,-Mates(i), 3-k)
        call CalcCLL(-Mates(i), 3-k)
        if (any(abs(XPR_x - XPr(3,:,:,-Mates(i), 3-k)) > 0.01))  OK = .FALSE.
      endif
    endif
  enddo
  XPR_x = XPr(3,:,:,-A,k)
  call CalcCLL(-A,k)
  if (OK .and. all(abs(XPR_x - XPr(3,:,:,-A,k)) < 0.01))  exit    ! What tolerance?
enddo

call getOff(A, k, .TRUE., nOff, Off, sxOff)  ! includes dummy offspring
do i=1, nOff
  if (nOff==0)  exit
  if (Off(i) > 0) then
    if (nFS(Off(i)) > 0 .and. Mates(i) < 0) then
      if (ns(-Mates(i),3-k) <= 50 .or. nFS(Off(i)) >= ns(-Mates(i),3-k)/5) then
        if (ns(-Mates(i), 3-k) > 0) then
          do j=1, ns(-Mates(i), 3-k)
            call CalcLind(SibID(j, -Mates(i), 3-k))
          enddo
        endif
      endif
    endif
    call CalcLind(Off(i))                     
  else
    call CalcCLL(-Off(i), sxOff(i))
  endif
enddo

!if (nOff > 0)  call CalcCLL(-A,k)

end subroutine UpdateLL

! #####################################################################

subroutine MakeFS(A, B)
use Global
implicit none

integer, intent(IN) :: A,B
integer :: x, i, j, Ai, Bj

if (nFS(A)>0) then
  Ai = A
else
  Ai = FSID(maxSibSize+1, A)
endif
if (nFS(B)>0) then
  Bj = B
else
  Bj = FSID(maxSibSize+1, B)
endif

if (ANY(FSID(1:nFS(Ai),Ai)==B) .or. ANY(FSID(1:nFS(Bj),Bj)==A)) then
  return ! already are FS.
endif

i = MIN(Ai,Bj)
j = MAX(Ai,Bj)
do x=1, nFS(j)   
  FSID(nFS(i)+x, i) = FSID(x, j)
  FSID(maxSibSize+1, FSID(x,j)) = i
enddo
nFS(i) = nFS(i) + nFS(j)
FSID(maxSibSize+1,i) = i    ! 'primary' sib
FSID(:,j) = 0
FSID(1,j) = j
FSID(maxSibSize+1,j) = i
nFS(j) = 0

end subroutine MakeFS

! ######################################################################

subroutine RemoveFS(A)
use Global
implicit none

integer, intent(IN) :: A
integer :: op, np, i, j

if (nFS(A) == 1) then
!  call CalcFSLik(A)
  return
else if (nFS(A) > 1) then
  op = A
  np = MINVAL(FSID(1:nFS(A), A), MASK=(FSID(1:nFS(A),A)/=A)) 
else !if (nFS(A) == 0) then
  op = FSID(maxSibSize+1, A)  ! 'primary' sib
  np = op
endif

i = 2  ! 1st one stays op
do j=1, nFS(op)
  if (FSID(j,op)==A) then
    FSID(j,op) = 0  ! if nFS(A)=0 .and. nFS(op)=2
    cycle
  endif
  if (FSID(j,op)==np) cycle
  FSID(i, np) = FSID(j, op)
  if (op /= np) then
    FSID(maxSibSize+1, FSID(j, op)) = np
  endif
  i = i+1
enddo

nFS(np) = nFS(op)-1
FSID(maxSibSize+1, np) = np
nFS(A) = 1
FSID(:,A) = 0
FSID(1,A) = A
FSID(maxSibSize+1, A) = A

!call CalcFSLik(op)  ! old primary
!call CalcFSLik(np)  ! new primary of no-longer-FS
!call CalcFSLik(A)

end subroutine RemoveFS

! ######################################################################

subroutine CheckSelfed(A, kA)
use Global
implicit none

integer, intent(IN) :: A, kA
integer :: AS(2), m, j, Aj, x
double precision :: LRself

if (hermaphrodites==0)  return

if (A > 0) then
  call IsSelfed(A, .FALSE., LRself)
  if (all(Parent(A,:) == 0)) then
    if (LRself > TA) then   ! threshold? be consistent with selectparent()
      SelfedIndiv(A) = .TRUE.  
    ! else leave as is?
    endif
    return
  else if (any(Parent(A,:) > 0)) then  
    if (Parent(A,1)==Parent(A,2)) then
      if (LRself < 5*TF) then    
        print *, ""
        print *, A, "; ", Parent(A,:), LRself
        call Erstop("CheckSelfed: dam = sire, but LRself < 5*TF", .TRUE.)
      else
        SelfedIndiv(A) = .TRUE.
      endif
    else if (all(Parent(A,:)/=0) .and. Parent(A,1)/=Parent(A,2)) then
      if (LRself > TA) then   
        print *, A, "; ", Parent(A,:), LRself
        call Erstop("CheckSelfed: dam /= sire, but LRself > TA", .TRUE.)
      else
        SelfedIndiv(A) = .FALSE.  
      endif
    ! else assignment in progress, leave as is?
    endif
  endif
  if (all(Parent(A,:) >= 0))  return
endif

AS = 0
if (A > 0)  AS = -Parent(A,:)
if (A < 0)  AS(kA) = -A

do m=1,2
  if (AS(m) <= 0)  cycle
  do j=1, ns(AS(m),m)
    Aj = SibID(j,AS(m),m)
    if (nFS(Aj)==0)  cycle
    call IsSelfed(Aj, .TRUE., LRself)
    if (LRself > TA) then
      if (Parent(Aj,3-m) < 0) then 
        do x = 1, nFS(Aj)
          SelfedIndiv(FSID(x,Aj)) = .TRUE.
        enddo
        DumClone(AS(m), m) = -Parent(Aj,3-m)
        DumClone(-Parent(Aj,3-m), 3-m) = AS(m)
      else if (Parent(Aj,3-m) > 0) then
        print *, Aj, "; ", Parent(Aj,:), LRself, m
        call Erstop("SetPar: parents incompatible with LRself > TA", .TRUE.)
      else   ! transitory during assignment only (?)
        do x = 1, nFS(Aj)
          SelfedIndiv(FSID(x,Aj)) = .TRUE.
        enddo
      endif
    endif
  enddo
enddo
  
end subroutine CheckSelfed

! #####################################################################

! @@@@   INPUT & PRECALC PROB.   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

! #####################################################################


subroutine ReadLifeHist(LifehistFileName, BYrange)
use Global
use sqa_fileIO, ONLY: FileNumCol, FileNumRow, IOstat_handler
implicit none

character(len=nchar_filename), intent(IN) :: LifehistFileName
integer, intent(INOUT) :: BYrange(nInd, 2)
integer :: k,i,m, numcolLH, ios, dumI(5) 
integer, allocatable, dimension(:) :: SexLH, ByLH, LyLH
integer, allocatable, dimension(:,:) :: BYrangeLH
character(len=nchar_ID), allocatable, dimension(:) :: NameLH
character(len=nchar_ID) :: dumC
logical, allocatable :: done(:)
logical :: dupLH

numcolLH = FileNumCol(trim(LifehistFileName))
nIndLH = FileNumRow(trim(LifehistFileName)) -1  ! has header

if (numcolLH/=3 .and. numcolLH/=5 .and. numcolLH/=6) then
  call Erstop("Invalid number of columns in file "//trim(LifehistFileName), .FALSE.)
endif  

allocate(SexLH(nIndLH))
allocate(NameLH(nIndLH))
allocate(ByLH(nIndLH))
allocate(ByRangeLH(nIndLH, 2))
allocate(LyLH(nIndLH))
SexLH = -999
NameLH = ' '
BYLH = -999
ByRangeLH = -999
LyLH = -999
dumI = -999

open(unit=103, file=trim(LifehistFileName), status="old")
  read(103, *)
  do k=1, nIndLH
    if (numcolLH==3) then
      read(103,*,IOSTAT=ios) dumC, dumI(1:2)
    else if (numcolLH==5) then
      read(103,*,IOSTAT=ios) dumC, dumI(1:4)
    else   ! ignore columns 7+
      read(103,*,IOSTAT=ios) dumC, dumI
    endif
    if (.not. any(ID == dumC))  cycle  ! not genotyped
    call IOstat_handler(ios, k, LifehistFileName)
    NameLH(k) = dumC
    SexLH(k) = dumI(1)
    BYLH(k) = dumI(2)
    if (numcolLH>=5) ByRangeLH(k,:) = dumI(3:4)
    if (numcolLH==6) LyLH(k) = dumI(5)
  enddo
close(103)

! rearrange lifehistory info to same order as genotype file
allocate(done(nInd))
done = .FALSE.
dupLH = .FALSE.
do i=1,nInd
  do k=1,nIndLH 
    if(Id(i) == NameLH(k)) then
      if (done(i))  dupLH = .TRUE.
      done(i) = .TRUE.
      if (SexLH(k) >= 1 .and. SexLH(k) <= 4)  Sex(i) = SexLH(k)
      if (BYLH(k) >= 0) then
        BY(i) = BYLH(k)
      else
        do m=1,2
          if (BYRangeLH(k,m) >= 0)  BYRange(i,m) = BYRangeLH(k,m)
        enddo
      endif
      if (LYLH(k) >= 0) then
        YearLast(i) = LYLH(k)
        if (YearLast(i) < BY(i))  YearLast(i) = -999   ! TODO: warning?
      endif
!      exit  ! increases runtime
    endif    
  enddo
enddo

if (dupLH) then
  print *, ""
  print *, "WARNING: Some IDs are duplicated in LifeHistData"
  print *, "Sex and BirthYear from last instance will be used"
  print *, ""
endif

deallocate(ByLH)
deallocate(SexLH)
deallocate(ByRangeLH)
deallocate(NameLH)

end subroutine ReadLifeHist

! #####################################################################

subroutine ReadPedFile(FileName)
use Global
use sqa_fileIO, ONLY: FileNumRow, IOstat_handler
implicit none

character(len=*) :: FileName
integer :: i, j, k, s, ios, nIndP, m, nnd(2), x, nC_orig(2), y, a, b, nchar_dpfx(3), &
  z, NewDum
integer, parameter :: DumMaleClone = -8888
character(len=nchar_ID) :: tmpC(3), NamePed(3,nInd*2), DumName, Navn
character(len=4) :: s_char
character(len=200) :: DataFMT
logical ::  IsClonedDum, RecodeDums, NameOK

Parent = 0  ! reset 
nIndP = 0
NamePed = "NA"
nC = 0
if (DoSibs) then
  ns = 0
  NewDum = -9999
else
  NewDum = 0
endif

nIndP = FileNumRow(trim(FileName)) -1  ! 1st row = header

if (quiet < 1)  call printt("Reading pedigree in "//trim(FileName)//" ... ")
open(unit=103, file=trim(FileName), status="old")
  read(103,*)   ! header                              
  do i=1,nIndP
    read(103, *,IOSTAT=ios)  tmpC
    call IOstat_handler(ios, i, FileName)
    NamePed(:,i) = tmpC 
  enddo
close(103)

! any non-genotyped, non-standard-dummy-codes in input pedigree?
! i.e., need for recoding dummies?
RecodeDums = .FALSE.
do z=1,3
  nchar_dpfx(z) = LEN_TRIM(DumPrefix(z))
enddo
do j = 1, nIndP
  do x=1,3
    if (NamePed(x,j) == 'NA')  cycle
    if (ANY(ID == NamePed(x,j)))  cycle
    Navn = NamePed(x,j)
    NameOK = .FALSE.
    do z=1,3
      if (Navn(1:nchar_dpfx(z)) == DumPrefix(z) .and. &
        Verify(Navn((nchar_dpfx(z)+1):nchar_ID), "0123456789 ")==0) then
        NameOK = .TRUE.
        exit
      endif
    enddo
    if (.not. NameOK) then
      RecodeDums = .TRUE.
      exit
    endif
  enddo
enddo

! parent names to nums & fix order
!if (quiet == -1)  print *, 'changing parent names to genotype file row numbers ...'
nnd = 0
do j = 1, nIndP
  call NameToNum(NamePed(1,j), RecodeDums, i, IsClonedDum)
  
  if (i > 0) then  ! genotyped individual; do non-genotyped below.
    do k = 1,2
      call NameToNum(NamePed(k+1,j), RecodeDums, x, IsClonedDum)
      Parent(i,k) = x
      if (x == NewDum .and. DoSibs) then
        if (.not. any(DummyNamesIO(:,k) == NamePed(k+1,j))) then
          nnd(k) = nnd(k) +1
          DummyNamesIO(nnd(k), k) = NamePed(k+1,j)   ! further processing below
        endif
      else if (IsClonedDum .and. k==2 .and. DoSibs) then   ! hermaphrodite cloned dummy
        Parent(i,k) = DumMaleClone
      else if (x < 0 .and. DoSibs) then
        if (-x > nInd/2 ) call Erstop('Something wrong with pedigree input: number of sibships > nInd/2', .FALSE.)
        if (nC(k) < -x)  nC(k) = -x
        nS(-x,k) = ns(-x,k) +1
        SibID(ns(-x,k), -x, k) = i
      endif
    enddo
  endif
enddo

 ! find numbers for male dummy clones; ensure no gaps & no overlap
if (hermaphrodites /= 0 .and. DoSibs) then   
  k = 2    
  do j = 1, nIndP
    call NameToNum(NamePed(1,j), RecodeDums, i, IsClonedDum)
    if (Parent(i,k) /= DumMaleClone)  cycle
    call NameToNum(NamePed(k+1,j), RecodeDums, x, IsClonedDum)
    s = 0 
    do y=1, nC(k)   
      if (ns(y,k) == 0) then
        s = y
        exit
      endif
    enddo
    if (s==0) then
      nC(k) = nC(k) +1
      s = nC(k)
    endif
    Parent(i,k) = -s
    nS(s,k) = ns(s,k) +1
    SibID(ns(s,k), s, k) = i
    if (Parent(i,3-k) >=0 .or. Parent(i,3-k)/=x) then
      call ErStop("Something went wrong reading dummy clones from pedigree", .TRUE.)
    endif
    DumClone(s, k) = -x
    DumClone(-x, 3-k) = s
    print *, "new dumclone: ", x, s
    ! find other offspring of s
    do b = j+1, nIndP
      if (NamePed(k+1, b) == NamePed(k+1, j)) then
        do a = 1, nInd
          if (NamePed(1,b) /= Id(a)) cycle
          Parent(a,k) = -s
          nS(s,k) = ns(s,k) +1
          SibID(ns(s,k), s, k) = a
        enddo
      endif
    enddo
  enddo
endif


! non-genotyped non-dummy-codes to nums
nC_orig = nC
if (any(nnd > 0)) then
!  if (quiet == -1)  print *, 'changing non-genotyped non-dummy-codes to numbers ...'
  
  if (hermaphrodites /= 0) then
    call Erstop("Non-genotyped non-dummy IDs not implemented for hermaphrodites", .FALSE.)
  endif

!  write(DataFMT, '( "(a", I0, ", 4X, a10)" )')  ID_len
!  open (unit=201,file="Dummified.txt", status="unknown") 
!  write (201, DataFMT) "id", "dummy"
  
  do k=1,2
    do x=1, nnd(k)
      nC(k) = nC(k) + 1
      s = nC(k)
      write(s_char, '(i4.4)') s
      DumName = trim(DumPrefix(k))//s_char
      write(201, DataFmt) DummyNamesIO(x,k) , DumName
      
      do i = 1, nInd
        if (Parent(i,k)/=NewDum)  cycle
        do j = 1, nIndP
          if (NamePed(1,j) /= Id(i)) cycle
          if (NamePed(k+1, j) /= DummyNamesIO(x,k)) cycle
          Parent(i,k) = -s
          nS(s,k) = ns(s,k) +1
          if (ns(s,k) > maxSibSize) then
            print *, ""
            print *, k, x, trim(DummyNamesIO(x,k)), s, "; ", j, i, Id(i)
            call Erstop("Sibship size in pedigree exceeds max, please increase 'MaxSibshipSize'", .FALSE.)
          endif
          SibID(ns(s,k), s, k) = i
        enddo
      enddo
    enddo
  enddo
  
!  close (201)
endif

! dummy's parents (sibship grandparents)     ! TODO: dumclones.
!if (ANY(nc>0) .and. quiet == -1)  print *, 'fixing numbering of sibship grandparents ...'
do k=1,2
  if (nC(k)==0)  cycle
  do s=1, nC(k)
    if (s <= nc_orig(k)) then
      write(s_char, '(i4.4)') s
      DumName = trim(DumPrefix(k))//s_char
    else
      DumName = DummyNamesIO(s - nC_orig(k), k)
    endif
    do j=1, nIndP
      if (NamePed(1,j) == DumName) then
        do m=1,2
          call NameToNum(NamePed(m+1,j), RecodeDums, GpID(m,s,k), IsClonedDum)
          if (GpID(m,s,k) == NewDum) then
            if (any(DummyNamesIO(:,m) == NamePed(m+1,j))) then
              do x=1, nnd(m)
                if (NamePed(m+1,j) == DummyNamesIO(x,m)) then
                  GpID(m,s,k) = -(Nc_orig(m) +x)
                  exit
                endif
              enddo
            else
              GpID(m,s,k) = 0   ! only dummify if >=0 SNPd offspring
            endif
          endif
        enddo
        exit
      endif
    enddo
  enddo
enddo

if (any(Parent /=0)) then
!  if (quiet == -1)  print *, 'inferring parent sex & full siblings ...'
  do i=1, nInd
    if (Sex(i)==3) then
      if (ANY(Parent(:,1) == i)) then
        Sex(i) = 1
      else if (ANY(Parent(:,2) == i)) then
        Sex(i) = 2
      endif                                     
    endif
  enddo    

  ! find current FS 
  do i=1,nInd-1
    do j=i+1,nInd
      if (Parent(i,1)==Parent(j,1) .and. Parent(i,1)/=0 .and. &
        Parent(i,2)==Parent(j,2) .and. Parent(i,2)/=0) then
        call MakeFS(i, j)
      endif
    enddo
  enddo  
endif

if (Complx == 0 .and. any(Parent /=0))  call CheckMono

if (hermaphrodites/=0) then
  do k=1,2
    if (nC(k)==0)  cycle
    do s=1, nC(k)
      call CheckSelfed(-s, k)    
    enddo
  enddo
  do i=1,nInd
    call CheckSelfed(i, Sex(i))
  enddo
endif

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
contains
  subroutine NameToNum(Navn, NewDumNums, Num, IsClonedHerm)
  implicit none
  
  character(len=nchar_ID), intent(IN) :: Navn
  logical, intent(IN) :: NewDumNums
  integer, intent(OUT) :: Num
  logical, intent(OUT) :: IsClonedHerm
  integer :: j, s

  Num = 0
  IsClonedHerm = .FALSE.
  if (Navn == "NA")  return
  
  do j=1, nInd
    if (Navn == Id(j)) then
      Num = j
      exit
    endif
  enddo
  if (Num /= 0 .or. .not. DoSibs) return
  
  if (.not. NewDumNums) then
    do z=1,3
      if (Navn(1:nchar_dpfx(z)) == DumPrefix(z)) then
        read(Navn((nchar_dpfx(z)+1):(nchar_dpfx(z)+5)), '(i4)') s
        Num = -s
        if (z==3)  IsClonedHerm = .TRUE.
      endif
    enddo
  else 
    Num = NewDum
  endif

  end subroutine NameToNum
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
end subroutine ReadPedFile

! ######################################################################

subroutine ReadOnlyList(FileName)
use Global
use sqa_fileIO, ONLY: FileNumCol, FileNumRow, IOstat_handler
implicit none

character(len=nchar_filename), intent(IN) :: FileName
integer :: x, i, nrows, ios, ncol
character(len=nchar_ID) :: tmpC, tmpX

nrows = FileNumRow(trim(FileName))  ! no header
ncol  = FileNumCol(trim(FileName))                                  

skip = .TRUE.

! single column (ignore all other columns)
if (quiet < 1)  call printt("Reading individuals in --only file "//trim(FileName)//" ... ")
open(unit=103, file=trim(FileName), status="old")
  do x=1, nrows
    if (ncol==2) then   ! PLINK format: FID + IID column
      read(103, *,IOSTAT=ios)  tmpX, tmpC
    else
      read(103, *,IOSTAT=ios)  tmpC
    endif
    call IOstat_handler(ios, x, FileName)
    do i=1, nInd
      if (Id(i) == tmpC) skip(i) = .FALSE.
    enddo
  enddo
close(103)

end subroutine ReadOnlyList

! ######################################################################

subroutine ReadMt(FileName)
use Global
use sqa_fileIO, ONLY: FileNumCol, FileNumRow, IOstat_handler
implicit none

character(len=nchar_filename), intent(IN) :: FileName
integer :: nRows, nCols, x, ios, i,j,y, matchedID
integer, allocatable :: mtSame_temp(:,:)
character(len=nchar_id), allocatable, dimension(:) :: rownames, colnames

nrows = FileNumRow(trim(FileName)) -1
ncols = FileNumCol(trim(FileName)) -1

allocate(mtSame_temp(ncols, nrows))   
allocate(rownames(nrows))
allocate(colnames(ncols))
mtSame_temp = 1
rownames = 'NA'
colnames = 'NA'

if (quiet < 1)  call printt("Reading mt similarity in "//trim(FileName)//" ... ")

open(unit=103, file=trim(FileName), status="old")
  read(103, *)  colnames
  do x=1, nrows
    read(103, *,IOSTAT=ios)  rownames(x), mtSame_temp(:,x)
    call IOstat_handler(ios, x, FileName)
  enddo
close(103)

! rearrange to match order in genotype data (mtSame_temp not necessarily square)
allocate(mtDif(0:nInd,0:nInd))  ! different mt haplotypes
mtDif = .FALSE.
matchedID = 0
do i=1, nInd
  do x=1, nrows
    if (rownames(x) /= Id(i))  cycle
    do j=1, nInd
      do y=1, ncols
        if (colnames(y) /= Id(j)) cycle
        matchedID = matchedID +1
        if (mtSame_temp(y,x) == 0) then
          mtDif(i,j) = .TRUE.
          mtDif(j,i) = .TRUE.
        endif
      enddo
    enddo
  enddo
enddo

if (matchedID == 0) then
  call Erstop('IDs in --mt do not match those in genotype file', .FALSE.)
endif

end subroutine ReadMt

! ######################################################################

subroutine ReadAgePrior(AgePriorFileName, AP_TMP)
use Global
use sqa_fileIO, ONLY: FileNumCol, FileNumRow, IOstat_handler
implicit none

character(len=nchar_filename), intent(IN) :: AgePriorFileName
double precision, intent(OUT) :: AP_TMP(MaxMaxAgePO, 5)
integer :: r,x,y, numcol,  ios, numrow
character(len=3) :: headAPfile(9), headAP(5)
double precision :: AP_IN(MaxMaxAgePO, 9)

!=================
headAP = (/"M  ", "P  ", "FS ","MS ", "PS "/)

numcol = FileNumCol(trim(AgePriorFileName))   ! default: "AgePriors.txt"
if (numcol/=8 .and. numcol/=9 .and. numcol/=5) then
  call Erstop("Invalid number of columns in "//trim(AgePriorFileName), .FALSE.)
endif
numrow = FileNumRow(trim(AgePriorFileName)) -1  ! header row
if (numrow > MaxMaxAgePO) then
  call ErStop("Max parent age >99: increase 'MaxMaxAgePO' on line 26 of source code", .TRUE.)
endif

AP_IN = 0.0D0
open(unit=102, file=trim(AgePriorFileName), status="old")
  read(102,*)  headAPfile(1:numcol)
  do y=1, numrow 
    read(102, *, iostat=ios)   AP_IN(y, 1:numcol)
    call IOstat_handler(ios, y, AgePriorFileName)
  enddo
close(102)

! fix order of columns
AP_TMP = 0.0D0
do r=1,5
  if (ANY(headAPfile == headAP(r))) then
    do x=1,numcol
      if (headAPfile(x) == headAP(r)) then
        AP_TMP(:, r) = AP_IN(:, x)
      endif
    enddo
  else
    call Erstop("column missing from ageprior file! "//headAP(r), .FALSE.)
  endif
enddo

end subroutine ReadAgePrior

! ######################################################################
! ######################################################################

subroutine PrepAgeData(AP_IN, BYrange)
use Global
implicit none

double precision, intent(IN) :: AP_IN(MaxMaxAgePO,5)
integer, intent(INOUT) :: BYrange(nInd, 2)   ! OUT not used, but else not modifiable                                         
integer :: i, BYLast, r, x, y, rik, rkj
double precision :: scl
double precision, allocatable, dimension(:,:) :: BYP, APtmp

! allocate(AgeDiff(nInd,nInd))
! AgeDiff = 999
! do i=1, nInd
  ! do j=1, nInd
    ! if (BY(i)>=0 .and. BY(j)>=0) then
      ! AgeDiff(i,j) = BY(i) - BY(j)   ! if >0, then j older than i
    ! endif   
  ! enddo
! enddo

!===  determine first & last possible birth year  ==============   
BYzero = MINVAL(BY, MASK=BY>=0) -1   ! defaults to HUGE(ARRAY)
BYlast = MAXVAL(BY, MASK=BY>=0)      ! defaults to -HUGE(ARRAY)
if (ANY(BYrange >= 0)) then
  BYzero = MIN(BYzero, MINVAL(BYrange(:,1), MASK=BYrange(:,1)>=0) -1)
  BYlast = MAX(BYlast, MAXVAL(BYrange(:,2), MASK=BYrange(:,2)>=0))
endif

!===  determine MaxAgePO  ==============   
maxAgePO = 1  ! maximum PO age difference (needed for dummy parents)
do y = 2, MaxMaxAgePO
  if (ANY(AP_IN(y, 1:2)>TINY(0D0))) then
    maxAgePO = y - 1  ! first y is agediff of 0
  endif
enddo

if (BYzero < 99999) then
  BYzero = BYzero - MaxAgePO +1
  if ((BYlast - BYzero) < 2*MaxAgePO .or. BYlast==0) then
    BYzero = BYzero - MaxAgePO    ! dummy parents + real grandparents w unknown BY
  endif  
  nYears = BYlast - BYzero   ! defines nYears!
else   ! all birth years unknown
  BYzero = 0
  nYears = MaxMaxAgePO
  MaxAgePO = INT( MaxMaxAgePO/2.0 )
endif

if (nYears==1)  nYears = 2   ! fail safe, shouldn't occur. avoids downstream problems

!===  initiate AgePrior array  ==============
allocate(AgePriorA(-MaxAgePO : nYears, 5, 3))
AgePriorA = 0.0D0

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! AgePriorA D2 + D3:
!  1    2     3
!  M   MGM   PGM
!  P   MGP   PGP
! FS   MFA   PFA
! MS  MMHA  PMHA
! PS  MPHA  PPHA
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~

do r=1,5
  AgePriorA(0:MaxAgePO, r, 1) = AP_IN(1:(MaxAgePO+1), r)
enddo
do r = 3,5  ! FS,MS,PS
  do y=1, MaxAgePO
    AgePriorA(-y, r, 1) = AgePriorA(y, r, 1)
  enddo
enddo
AgePriorA(0, 1:2, 1) = 0.0D0   ! PO CANNOT have age diff of 0. 

if (ALL(AP_IN(:,1:2) < 0.001)) then  ! AP_IN with single row for sibs only
  AgePriorA(1:MaxAgePO, 1:2, 1) = 1.0D0
endif

! calc GP & AU (for individual i born in year 0)
allocate(BYP(3, -MaxAgePO : nYears))
BYP = 0D0
scl = 1D0 / (SUM(AgePriorA(:, 1:2, 1)) / 2.0)  ! scaling factor; similar sum as first 5 columns
allocate(APtmp(-MaxAgePO : nYears, -MaxAgePO : nYears))
do rik = 1,2
  BYP(1:2, :) = 0D0
  BYP(1,:) = AgePriorA(:, rik, 1)  ! possible birth years of in-between indiv k
  BYP(1,:) = BYP(1,:) / SUM(BYP(1,:))
  do rkj = 1,5    
    APtmp = 0D0
    do y = -MaxAgePO, nYears  ! birth year of in-between indiv k
      if (BYP(1,y) < 0.001)  cycle  ! possible birth years of indiv j:
      BYP(2, (y-MaxAgePO) :(y+MaxAgePO)  ) = AgePriorA(-MaxAgePO : MaxAgePO, rkj, 1)  
      BYP(2,:) = BYP(2,:) / SUM(BYP(2,:))
      APtmp(y,:) = BYP(1,y) * BYP(2,:) 
    enddo
    if (ALL(AP_IN==0.0 .or. AP_IN==1.0)) then
      WHERE (SUM(APtmp, DIM=1) > 0D0)  AgePriorA(:, rkj, rik+1) = 1D0
    else
      AgePriorA(:, rkj, rik+1) = SUM(APtmp, DIM=1) / scl
    endif
  enddo
enddo
deallocate(BYP)
deallocate(APtmp)

! ! ageprior for GGP
! allocate(AP_GG(-MaxAgePO:nYears, 0:3, 0:3))
! AP_GG = 0D0
! AP_GG(:,1,0) = AgePriorA(:,1,1)  ! dam
! AP_GG(:,0,1) = AgePriorA(:,2,1)  ! sire
! AP_GG(:,1,1) = AgePriorA(:,1,3)  ! PGM / MGF
! AP_GG(:,2,0) = AgePriorA(:,1,2)  ! MGM
! AP_GG(:,0,2) = AgePriorA(:,2,3)  ! PGF
! do i=0,3  ! no. offspring-dam steps
  ! do j=0,3  ! no. offspring-sire steps
    ! if (i+j<=2 .or. i+j>3)  cycle   ! already done / not required
    ! do y = -MaxAgePO, nYears  ! age difference between the two individuals
          
    ! enddo
  ! enddo
! enddo


!===  shift BY  ==============
WHERE (BY >=0) BY = BY - BYzero  
do x=1,2
  WHERE (BYRange(:,x) >=0) BYRange(:,x) = BYRange(:,x) - BYzero
enddo 
if (AnyYearLast) then
  WHERE(YearLast <= 0)  YearLast = 999                                    
  WHERE (YearLast /= 999)  YearLast = YearLast - BYzero

  ! min/max BY based on YearLast
  do i=1, nInd
    if (YearLast(i)==999)  cycle
    if (BYrange(i,1) < 0)  BYrange(i,1) = MAX(YearLast(i) - MaxAgePO, 1)
    if (BYrange(i,2) < 0)  BYrange(i,2) = MIN(YearLast(i), nYears)
  enddo
endif

WHERE (BYRange(:,1) <0) BYrange(:,1) = 1
WHERE (BYRange(:,2) <0) BYrange(:,2) = nYears 

!===  Initiate indiv BY prob distr  ==============
if (any(BY < 0)) then
  allocate(IndBY(nYears, nInd, 5))  ! year - indiv - own/wo/w dummy off+par
  IndBY = LOG10(1.0D0/nYears)
  do i=1, nInd  
    if (BY(i) >0) then
      IndBY(:, i, :) = LOG10(zero)
      IndBY(BY(i), i, :) = zero  
    else 
      if (BYrange(i,1) < 1)  call ErStop('BY.min < 1', .TRUE.)
      if (BYrange(i,2) > nYears)  call ErStop('BY.max > nYears', .TRUE.)                                                                                                                              
      IndBY(:, i, :) = LOG10(zero)
      IndBY(BYrange(i,1) : BYrange(i,2), i, :) = LOG10(1.0D0/(BYrange(i,2) - BYrange(i,1) +1))
    endif
  enddo
endif

if (DoSibs) then
  allocate(DumBY(1:nYears, nInd/2, 2,5)) 
  DumBY = 0D0
endif

call writeAgePrior  ! FOR CHECKING

end subroutine PrepAgeData

! ######################################################################

subroutine EstBYrange(A, k, MCI)  
use Global      
implicit none

integer, intent(IN) :: A, k
integer, intent(OUT) :: MCI(3)  ! mode - 95% lower - 95% upper
integer :: y, mx, CI(2)
double precision :: DBYP(nYears), cumProp, dd(nYears)

MCI = -9
if (A == 0) return
call getEstBY(A,k, 5, DBYP)  ! all contributions from relatives
DBYP = 10**DBYP
mx = MAXLOC(DBYP, DIM=1)
cumProp = DBYP(mx)
CI = mx
do y=1, nYears 
  if (cumProp > 0.95) then
    dd = 0D0
    WHERE (DBYP > 0.0D0)  dd = ABS(DBYP - DBYP(mx))
    if (ANY(dd > 0.001) .or. CI(1)==CI(2)) then  ! else: DBYP is flat between BYmin & BYmax
      MCI(1) = mx
    endif
    MCI(2:3) = CI
    WHERE(MCI>0)  MCI = MCI + BYzero
    exit
  endif  

  if (CI(1) > 1 .and. CI(2) < nYears) then
    if (DBYP(CI(1)-1) > DBYP(CI(2)+1)) then
      CI(1) = CI(1)-1
      cumProp = cumProp + DBYP(CI(1))
    else 
      CI(2) = CI(2)+1
      cumProp = cumProp + DBYP(CI(2))
    endif
  else if (CI(1) > 1) then
    CI(1) = CI(1)-1
    cumProp = cumProp + DBYP(CI(1))
  else if (CI(2) < nYears) then
    CI(2) = CI(2)+1
    cumProp = cumProp + DBYP(CI(2))
  endif
enddo

end subroutine EstBYrange

! #####################################################################

subroutine getRank_i(BYrank)   ! get ranking based on generation number + birthyear (increasing)
use Global
use qsort_c_module
implicit none

integer, intent(OUT) :: BYrank(nInd)
double precision :: SortBY(nInd)
integer :: i, BYx(nInd), Gen(nInd)

Gen = 0
call getGenerations(Gen) 

BYRank = (/ (i, i=1, nInd, 1) /)
SortBY = REAL(Gen, 8)
call QsortC(SortBy, BYRank)

BYx = BY
do i=1, nInd
  if (BY(i) < 0) then
    BYx(i) = MAXLOC(IndBY(:,i,5), DIM=1)  ! incl parents & offspring w unknown BY
  endif
enddo

SortBY = REAL(BYx, 8)
call QsortC(SortBy, BYRank)  ! do earlier birth years before later birth years

end subroutine getRank_i

! #####################################################################

subroutine getBYrank_c(k, BYrank)   ! get birthyear ranking (increasing)
use Global
use qsort_c_module
implicit none

integer, intent(IN) :: k
integer, intent(OUT) :: BYrank(nInd/2)
double precision :: SortBY(nC(k))
integer :: s, BYx(nC(k)), BYrankk(nC(k))

BYx = -9
do s=1, nC(k)
  BYx(s) = MAXLOC(DumBY(:, s, k, 5), DIM=1)  ! with dummy parents & offspring
enddo

SortBY = REAL(BYx, 8)
BYrankk = (/ (s, s=1, nC(k), 1) /)
call QsortC(SortBY, BYRankk)

BYrank(1:nC(k)) = BYrankk

end subroutine getBYrank_c

! #####################################################################

subroutine getGenerations(Gen)
use Global
implicit none

integer, intent(OUT) :: Gen(nInd)
integer :: i, GenPar(2, nInd), g, Prevgen(nInd), m

Gen = -9
GenPar = 0   ! generation number of dam, sire
Prevgen = 0   ! ids of individuals in generations up to g-1

do i=1,nInd
  if (Parent(i,1)==0 .and. Parent(i,2)==0) then
    Gen(i) = 0  ! founder
    Prevgen(i) = i
  endif
enddo

do g = 0, 1000
  do i=1, nInd
    if (Gen(i) >= 0)  cycle
    do m=1,2
      if (Parent(i,m)==0 .or. GenPar(m,i) > 0)  cycle
      if (ANY(Prevgen == Parent(i,m))) then
        GenPar(m,i) = g
      endif
    enddo
    
    if (ALL(GenPar(:,i) <= g)) then  ! including Parent(i,m)==0 --> GenPar(m,i)==0
      Gen(i) = g+1
      PrevGen(i) = i
    endif    
  enddo
  if (.not. ANY(Gen < 0))  exit   ! all done
enddo

end subroutine getGenerations

! ##############################################################################

subroutine WriteDummies
use Global
implicit none

integer :: s, k, n, DumBYCI(3,nInd/2,2)
character(len=3) :: headTmp
character(len=5) :: OffHeader(maxSibSize)
character(len=nchar_ID) :: DumName(nInd/2,2), GpName(2,nInd/2,2)
character(len=200) :: HeaderFMT, DataFMT
logical :: IsClone2(nInd/2, 2)

DumBYCI = -9
do k=1,2
  if (nC(k)==0)  cycle 
  do s=1,nC(k)
    call EstBYrange(-s, k, DumBYCI(:,s,k))
  enddo
enddo

call MakeDumNames(DumName, IsClone2)
GpName = "NA"
do k=1,2
  do s=1,nC(k)
    do n=1,2
      if (GpID(n,s,k)>0) then
        GpName(n,s,k) = Id(GpID(n,s,k))
      else if (GpID(n,s,k)<0) then
        GpName(n,s,k) = DumName(-GpID(n,s,k),n) 
      endif
    enddo
  enddo
enddo
   
do s=1, MAXVAL(ns)
  write(headTmp, '(i3.3)') s
  OffHeader(s) = "O_"//headTmp
enddo

write(HeaderFMT, '( "(3(a", I0, ", 4X), a4, 4a8, 4X, 999(a",I0,", 2X))" )')  ID_len, ID_len
write(DataFMT,   '( "(3(a", I0, ", 4X), i4, 4i8, 4X, 999(a",I0,", 2X))" )')  ID_len, ID_len

open (unit=301,file="DummyParents.txt",status="replace")
  write (301, HeaderFMT)  "id", "dam", "sire", "sex", "est.BY", "min.BY", "max.BY", &
    "NumOff", OffHeader(1:maxval(ns))
  do k=1,2
    do s=1,nC(k)
      write (301, DataFMT) DumName(s,k), GpName(:, s,k), k, DumBYCI(:, s, k), &
          nS(s,k), ID(SibID(1:nS(s,k), s, k))
    enddo
  enddo
close (301)

end subroutine WriteDummies

! #####################################################################

subroutine WriteBYprob
use Global
implicit none

integer :: i, s, k, Years(nYears)
double precision :: BYLR(nYears)
character(len=nchar_ID) :: DumName(nInd/2,2) 
character(len=200) :: HeaderFMT, DataFMT   
logical :: IsClone2(nInd/2, 2)                      

Years = (/ (i, i=BYzero+1, BYzero+nYears) /)  

call MakeDumNames(DumName, IsClone2)

write(HeaderFMT, '( "(a", I0, ", 4X, 2a6, 999i7)" )')  ID_len
write(DataFMT,  '( "(a", I0, ", 4X, 2i6, 999f7.3)" )')  ID_len

open (unit=401, file="BirthYearProbabilities.txt", status="replace")
write (401, HeaderFMT) "id", "rowO", "Sex", Years
do i=1, nInd
  if (BY(i) > 0)  cycle
  call getEstBY(i, 0, 5, BYLR)
  write(401, DataFMT)  Id(i), i, Sex(i), 10**BYLR
enddo
do k=1,2
  do s=1, nC(k)
    write(401, DataFMT) DumName(s,k), -s, k, 10**DumBY(:,s,k,2)
  enddo
enddo  
close(401)

end subroutine WriteBYprob

! #####################################################################

subroutine writeAgePrior
use Global
implicit none

integer :: y
character(len=3) :: headAP(5, 3)  
headAP(:,1) = (/"M  ", "P  ", "FS ","MS ", "PS "/) 
headAP(:,2) = (/"MGM", "MGF", "MFA","MMA", "MPA"/) 
headAP(:,3) = (/"PGM", "PGF", "PFA","PMA", "PPA"/) 

open(unit=103, file="AgePriors_new.txt", status="unknown") 
write(103, '(a4, "  ", 5a7, "  ", 5a7, "  ", 5a7)') "  dA", headAP(:,1), headAP(:,2), headAP(:,3) 
do y = lbound(AgePriorA, DIM=1), ubound(AgePriorA, DIM=1)
  write(103, '(i4, "  ", 5f7.3, "  ", 5f7.3, "  ", 5f7.3)') y, AgePriorA(y,:, 1), &
    AgePriorA(y,:, 2), AgePriorA(y,:, 3)
enddo
close(103)

end subroutine writeAgePrior

! #####################################################################

subroutine deallocall
use Global
use sqa_general
use OHfun
implicit none

! same order as module Global; allocated in PrepData
if (allocated(ToCheck)) deallocate(ToCheck)
if (allocated(SelfedIndiv)) deallocate(SelfedIndiv)
if (allocated(skip)) deallocate(skip)
if (allocated(IsNewSibship)) deallocate(IsNewSibship)
if (allocated(mtDif)) deallocate(mtDif)                                                    

if (allocated(Sex)) deallocate(Sex)
if (allocated(BY)) deallocate(BY)
if (allocated(nFS)) deallocate(nFS)
if (allocated(Mate))  deallocate(Mate)
if (allocated(YearLast)) deallocate(YearLast)

if (allocated(Genos)) deallocate(Genos)
!if (allocated(AgeDiff)) deallocate(AgeDiff)
if (allocated(Parent)) deallocate(Parent)
!if (allocated(OppHomM)) deallocate(OppHomM)
if (allocated(nS)) deallocate(nS)
if (allocated(FSID)) deallocate(FSID)
if (allocated(DumMate))  deallocate(DumMate)
if (allocated(DumClone))  deallocate(DumClone)

if (allocated(SibID)) deallocate(SibID)
if (allocated(GpID)) deallocate(GpID)

if (allocated(Lind)) deallocate(Lind)
!if (allocated(LLR_O)) deallocate(LLR_O)
if (allocated(CLL)) deallocate(CLL)
! sqa_general:
if (allocated(AHWE)) deallocate(AHWE)
if (allocated(OHWE)) deallocate(OHWE)
if (allocated(OcA))  deallocate(OcA)
if (allocated(lOcA)) deallocate(lOcA)
if (allocated(AKAP)) deallocate(AKAP)
if (allocated(OKAP)) deallocate(OKAP)
if (allocated(OKOP)) deallocate(OKOP)
if (allocated(AKA2P)) deallocate(AKA2P)
if (allocated(OKA2P)) deallocate(OKA2P)
if (allocated(lAKA2P)) deallocate(lAKA2P)
  
!if (allocated(LindG)) deallocate(LindG)
if (allocated(PHS)) deallocate(PHS)
if (allocated(PFS)) deallocate(PFS)
if (allocated(PPO)) deallocate(PPO)
if (allocated(LindX)) deallocate(LindX)
if (allocated(IndBY)) deallocate(IndBY)
if (allocated(AgePriorA)) deallocate(AgePriorA)

if (allocated(DumP)) deallocate(DumP)
if (allocated(DumBY)) deallocate(DumBY)
!if (allocated(FSLik)) deallocate(FSLik)
if (allocated(XPr)) deallocate(XPr)

if (allocated(Id)) deallocate(Id)
if (allocated(DummyNamesIO)) deallocate(DummyNamesIO)

end subroutine deallocall

! #####################################################################

! -9   NA
! missing  NA
! AlreadyAss  Already assigned
! impossible  impossible
! NotImplemented  not yet implemented (typically involves inbreeding)
! MaybeOtherParent  as likely to go via opposite parent