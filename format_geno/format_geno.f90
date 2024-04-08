! re-format SNP genotype data between various formats
! Jisca Huisman, jisca.huisman@gmail.com
!
! This code is available under GNU General Public License v3
!
! Compilation: 
! gfortran -O3 sqa_fileIO.f90 format_geno.f90 -o formatter
! debug:
! gfortran -g -fall-intrinsics -Wall -pedantic -fbounds-check -Og sqa_fileIO.f90 format_geno.f90 -o formatter


!===============================================================================

module Global_variables
  implicit none

 ! integer, parameter :: ishort = selected_int_kind(1)
  integer, parameter :: nchar_filename = 2000!, nchar_ID = 40
  character(len=nchar_filename) :: InFile, OutFile, GenoIN, GenoOUT
  character(len=3) :: InFormat, OutFormat
  logical :: quiet, make_map
  
  
contains
   ! option documentation
  subroutine print_help()
    print '(a, /)', ' ~~ Reformat SNP data between various formats ~~'
    print '(a, /)', 'command-line options:'
    print '(a)',    '  --help           print usage information and exit'
    print '(a)',    '  --in, --genoIN   <filename>  Original genotype data file, without file extension.'  
    print '(a)',    '  --out, --genoOUT <filename>   filename for re-formatted genotypes. An existing ', &                 
                    '                       file will be overwritten unless --warn is used'
    print '(a)',    '  --informat <x>   SEQ: no header, 0/1/2, missing -9, IDs in column 1; .txt', &  
                    '                   PED: no header, IDs in column 2 of 6 non-SNP columns, then',&
                    '                     2 columns per SNP coded 11/12/22, missing 00; .ped (+.map)',&
                    '                     input as A/C/T/G requires more memory; output is always 1/2',&
                    '                   RAW: header, 0/1/2/NA, IDs in column 2 of 6 non-SNP columns; .raw', &
                    '                   LMT: no header, 0/1/2 without spacing, IDs in separate file;', &
                    '                      .geno + .id. Missing values not allowed (coded 9 in output)'                     
    print '(a)',    '  --outformat <x>  same options as for --informat' 
    print '(a)',    '  --make-map        only if --outformat is PED: create associated .map file'   
    print '(a)',    '  --warn          warn if output file exist, & prompt for user whether to continue'   
    print '(a)',    '  --quiet         hide messages'
    print '(a)',    ''
  end subroutine print_help
  
end module Global_variables

!===============================================================================

program main
  use sqa_FileIO
  use Global_variables
  implicit none

  character(len=1) :: answer  
  logical :: FileExists, file_warn
  
  call read_args()
  
  ! ~~~~~~~~~~~~~~~~~~~~~~~  
  ! add file extensions
  GenoIN = add_extension(InFile, InFormat)
  GenoOUT = add_extension(OutFile, OutFormat)
  
  ! ~~~~~~~~~~~~~~~~~~~~~~~
  if (make_map .and. OutFormat /= 'PED') then
    print *, 'WARNING: --makemap only used with --outFormat PED'
  endif
  
  ! check if files exist
  call checkFile(GenoIN)
  
  if (file_warn) then
    inquire(file=trim(OutFile), exist = FileExists)
    if (FileExists) then
      print *, 'WARNING: Output file '//trim(OutFile)//' already exists. Press N to abort, or Y to continue'
      read(*,'(a1)') answer
      if (answer == 'N' .or. answer=='n')  stop
    endif
  endif
  
  call Reformat()   

  print *, 'done.'
  
contains
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !  read in command line arguments
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine read_args()
    integer :: nArg, i, x
    character(len=32) :: arg
  
    ! set default values
    InFile = 'geno'
    OutFile = 'geno'
    InFormat = 'ZZZ'
    OutFormat = 'ZZZ'
    file_warn = .FALSE.
    make_map = .FALSE.
    quiet = .FALSE.
    
    nArg = command_argument_count()
    i = 0
    do x = 1, nArg
      i = i+1
      if (i > nArg)  exit
      call get_command_argument(i, arg)
      
      select case (arg) 

        case ('--help')
          call print_help()
          stop    
        
        case ('--in', '--genoIN')  
          i = i+1
          call get_command_argument(i, InFile)
          
        case ('--out', '--genoOUT')  
          i = i+1
          call get_command_argument(i, OutFile)
          
        case ('--informat', '--inFormat')
          i = i+1
          call get_command_argument(i, InFormat)
          if (.not. any(valid_formats == InFormat)) then
            print *, 'ERROR: inFormat must be one of: ', valid_formats
            stop
          endif
          
        case ('--outformat', '--outFormat')
          i = i+1
          call get_command_argument(i, OutFormat)
          if (.not. any(valid_formats == OutFormat)) then
            print *, 'ERROR: outFormat must be one of: ', valid_formats
            stop
          endif
          
        case ('--make-map')
          make_map = .TRUE.
          
        case ('--quiet')
          quiet = .TRUE.
          
        case ('--warn')
          file_warn = .TRUE.
         
        case default
          print '(2a, /)', 'ERROR: Unrecognised command-line option: ', arg
          stop

      end select
    end do  

  end subroutine read_args
  
end program main

!===============================================================================

subroutine Reformat
  use Global_variables
  use sqa_FileIO
  implicit none
  
  integer :: nInd, nSnp, i, l, IOerr
  character(len=nchar_ID), allocatable :: ID(:)
  character(len=nchar_ID), allocatable :: SNP_names(:)   
  character(len=nchar_ID) :: dumC, dumV(4)
  integer(kind=ishort), allocatable :: G_int(:), G_duos(:), GenoM(:,:)
  character(len=2), allocatable :: G_char(:)
  logical :: FileExists, DidWarn

  nInd = FileNumRow(trim(GenoIN))
  if (InFormat == 'RAW') then
    nInd = nInd -1   ! header row
  endif
  
  allocate(Id(nInd))   ! vector only needed when InFormat or OutFormat is LMT
  Id = ''              ! (but could open LMT .id file in/out similtaneously w geno files)

  select case (InFormat)
    case ('SEQ')
      nSnp = FileNumCol(trim(GenoIN)) -1  ! column 1 = IDs
    case('PED')
      nSnp = (FileNumCol(trim(GenoIN)) -6)/2
    case('RAW')
      nSnp = FileNumCol(trim(GenoIN)) -6
    case('LMT')
      nSnp = FileNumCol(trim(GenoIN), sep='')
  end select  
  
  allocate(SNP_names(nSNP))
  SNP_names = ''
  allocate(G_char(nSnp))
  G_char = ''
  
  print *, '# individuals: ', nInd, '; # SNPs: ', nSnp  
           
  ! check if .ped file in A/C/T/G format; if so need to read in one go to ensure
  ! consistency which homozygote is coded '0' and which '2'.
  if (InFormat=='PED') then
    open (unit=101, file=trim(GenoIN), status='old', action='read')
      read (101,*) dumC, dumC, dumV, G_char
    close(101) 
    if (any(G_char=='A') .or. any(G_char=='C') .or. any(G_char=='G') .or. any(G_char=='T')) then
      allocate(GenoM(nSnp,nInd))    
      call read_geno(Geno=GenoM, ID=Id, SNP_names=SNP_names, FileName=InFile, &
        FileFormat = InFormat)   ! this also reads the .map file, if any
      call write_geno(Geno=GenoM, nInd=nInd, nSnp=nSNP, ID=Id, SNP_names=SNP_names, &
        FileName=OutFile, FileFormat=OutFormat, make_map=make_map)  
      deallocate(GenoM)      
      return
    endif
  endif
  

  allocate(G_int(nSnp))
  if (InFormat=='PED' .or. OutFormat=='PED') then
    allocate(G_duos(nSnp*2))   ! 2 columns per SNP format
  endif
   
  ! read in auxiliary files with SNP names (PED) or IDs (LMT)
  if (InFormat=='PED') then   ! .and. (OutFormat=='PED' .or. OutFormat=='RAW')
    inquire(file=trim(InFile)//'.map', exist = FileExists)
    if (FileExists) then
      open (unit=303, file=trim(InFile)//'.map', status='old', action='read')
        do l=1,nSnp
          read (303,*,IOSTAT=IOerr)  dumC, SNP_names(l) ! SNP names in 2nd column
          if (IOerr > 0) then
            print *, "ERROR: Wrong input in file "//trim(InFile)//".map on line ", l
            stop
          else if (IOerr < 0) then
            print *, "ERROR: file "//trim(InFile)//".map has fewer rows than no. SNPs in genotype file"
            stop
          endif
        enddo     
      close(303)
    ! else see below
    endif
  
  else if (InFormat=='LMT') then
    inquire(file=trim(InFile)//'.id', exist = FileExists)
    if (.not. FileExists) then
      print *, 'ID file '//trim(InFile)//'.id not found!'
      stop
    endif
    open (unit=303, file=trim(InFile)//'.id', status='old', action='read')
      do i=1,nInd
        read (303,*,IOSTAT=IOerr)  ID(i) 
        if (IOerr > 0) then
          print *, "ERROR: Wrong input in file "//trim(InFile)//".id on line ", i
          stop
        else if (IOerr < 0) then
          print *, "ERROR: ID file "//trim(InFile)//".id has fewer rows than genotype file"
          stop
        endif
      enddo
    close(303)
  endif
 
  
  if ((OutFormat=='RAW' .or. (OutFormat=='PED' .and. make_map)) .and. &
    SNP_names(1)=='' .and. InFormat /= 'RAW') then
    ! create fake SNP names
    do l=1,nSnp
      write(SNP_names(l), '("SNP",i6.6)') l 
    enddo
  endif
  
  DidWarn = .FALSE.
  
!  if (.not. quiet .and. nInd>5000) call timestamp(.TRUE.)
  open (unit=101, file=trim(GenoIN), status='old', action='read')
  open (unit=202, file=trim(GenoOUT), status='unknown', action='write')
  
  ! header
  if (InFormat=='RAW') then
    read (101,*) dumC, dumC, dumV, SNP_names
  endif
  if (OutFormat=='RAW') then
    write (202,'(a4,2x,a3,37x,3a5, a6, 2x, 200000a25)') 'FID', 'IID', 'PAT', 'MAT', 'SEX', 'PHENO', SNP_names
  endif
  
  do i=1,nInd
    if (nInd > 5000) then
      if (.not. quiet .and. MODULO(i,2000)==0)  write(*,'(i10, 2x)', advance='no') i
    endif
    
    select case (InFormat)
      case ('SEQ')
        read (101,*)  Id(i), G_int
      case('PED')
        read (101,*) dumC, Id(i), dumV, G_duos
        G_int = Two2One(G_duos)
      case('RAW')
        read (101,*) dumC, Id(i), dumV, G_char
        G_int = -1
        WHERE(G_char=='0')  G_int = 0
        WHERE(G_char=='1')  G_int = 1
        WHERE(G_char=='2')  G_int = 2
      case('LMT')
        read (101,'(500000i1)')  G_int
        WHERE(G_int > 2)  G_int = -1
    end select 
    
    ! character vector for output. set missing values first 
    if (OutFormat /= 'PED') then
      select case (OutFormat)
        case ('SEQ')
          G_char = '-1'
        case ('RAW')
          G_char = 'NA'
        case ('LMT')
          G_char = '9 '
      end select
      WHERE(G_int==0) G_char='0 '
      WHERE(G_int==1) G_char='1 '
      WHERE(G_int==2) G_char='2 '
    endif
    
    select case (OutFormat)
      case ('SEQ')
        write(202, '(a40, 100000a3)') Id(i), G_char
      case('PED')
        write(202, '(i3,2x,a40,4i3,2x, 200000i2)') 0, Id(i), 0,0,0,0, One2Two(G_int)
      case('RAW')
        write(202, '(i4,2x,a40,4i5,2x, 200000a3)') 0, Id(i), 0,0,0,0, G_char
      case('LMT')
        if (any(G_int < 0) .and. .not. DidWarn) then  ! warning once is enough
          print *, 'WARNING: LMT does not support missing values! Coded as 9 in output'
          DidWarn = .TRUE.         
        endif
        write(202,'(500000a1)')  G_char
    end select      

  enddo
  
  close(202)
  close (101)
  
  if (nInd > 5000) write(*,*) ''  ! advance to next line
  
  ! write out auxiliary files
  if (OutFormat=='PED' .and. make_map) then
    open(505, file=trim(OutFile)//'.map', status='unknown', action='write')
      ! no header
      do l=1,nSnp
        write(505,'(i2,1x, a25, 2i2)') 0, SNP_names(l), 0, 0
      enddo
    close(505)
  endif
  
  if (OutFormat=='LMT') then
    open(505, file=trim(OutFile)//'.id', status='unknown', action='write')
      do i=1,nInd
        write(505,'(a40)') ID(i)
      enddo
    close(505)
  endif
  
  
  ! deallocate all the stuff 
  ! TODO
  

  contains
    !~~~~~~~~~~~~~~
    ! From 1-column-per-SNP to 2-columns-per-SNP
    function One2Two(G)  
      integer(kind=ishort), intent(IN) :: G(0:(nSnp-1))
      integer(kind=ishort) :: One2Two(0:(nSnp*2-1))
      integer :: l,m
      
      do l=0,(nSNP-1)
        m = 2*l
        select case (G(l))
          case (0)
            One2Two(m:(m+1)) = (/1_ishort,1_ishort/)  
          case (1)
            One2Two(m:(m+1)) = (/1_ishort,2_ishort/)  
          case (2)
            One2Two(m:(m+1)) = (/2_ishort,2_ishort/)   
          case default  ! missing
            One2Two(m:(m+1)) = (/0_ishort,0_ishort/)   
        end select   
      enddo   
    end function One2Two
    
    !~~~~~~~~~~~~~~
    ! From 2-columns-per-SNP to 1-column-per-SNP
    function Two2One(G)  
      integer(kind=ishort), intent(IN) :: G(0:(2*nSnp-1))
      integer(kind=ishort) :: Two2One(0:(nSnp-1))
      integer :: l,m
      character(len=2) :: Gm
      
      do l=0,(nSNP-1)
        m = 2*l
        write(Gm,'(2i1)')  G(m:(m+1))
        select case (Gm)
          case ('11')
            Two2One(l) = 0  
          case ('12')
            Two2One(l) = 1
          case ('21')
            Two2One(l) = 1
          case ('22')
            Two2One(l) = 2  
          case ('00')
            Two2One(l) = -9 
          case default  ! shouldn't happen, but treat as missing
            Two2One(l) = -9    
        end select   
      enddo   
    end function Two2One

  
end subroutine Reformat

!===============================================================================
! end. 