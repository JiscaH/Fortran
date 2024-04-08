! 1. make imputation test data from real data
! set random prop_del proportion of genotypes to missing
! ignoring current genotype, which may already be missing

! 2. compare imputed genotype data to list of edits made

!===============================================================================

module global_variables 
!  use sqa_general, ONLY: ishort, nchar_ID
  implicit none
  
  integer, parameter :: ishort = selected_int_kind(1)
  integer, parameter :: nchar_ID=40
  
  integer :: nInd, nSnp
  integer(kind=ishort), allocatable :: Geno(:,:)
  character(len=nchar_ID), allocatable :: IdV(:), SNP_names(:)

end module global_variables

!===============================================================================

module stats
  implicit none
  
contains
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! sample n random integers between first & last, inclusive
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function runif_int(n, first, last)
    integer, intent(IN) :: n, first, last
    integer :: runif_int(n)    
    real :: r(n)
    
    call random_number(r)
    runif_int = floor((last-first+1)*r) + first
  end function runif_int

end module stats

!===============================================================================

program main
  use sqa_fileIO, ONLY: checkFile, IOstat_handler, FileNumRow, read_geno, write_geno
  use global_variables
  implicit none
  
  character(len=200) :: GenoFile, GenoOutFile, DelFile, ImpuFile, compare_file
  character(len=3) :: GenoInFormat, GenoOutFormat
  real :: prop_del
  logical :: do_del, do_compare
   
  ! read args: del vs compare
  call read_args()
  
  ! read genotype data
  if (do_del .or. (do_compare .and. ImpuFile == 'none')) then
    call read_geno(Geno=Geno, ID=IdV, SNP_names=SNP_names, FileName=GenoFile, &
      FileFormat = GenoInFormat, transp=.FALSE.)
    nInd = SIZE(Geno, DIM=1) -1 ! dim 0:nInd
    nSnp  = SIZE(Geno, DIM=2)
  endif 
  
  if (do_del) then
    call random_seed()
    call del_snps(prop_del, DelFile)
    call write_geno(Geno=transpose(geno(1:,:)), nInd=SIZE(geno,DIM=1)-1, nSnp=SIZE(geno,DIM=2),&
      ID=IdV, FileName=GenoOutFile, FileFormat=GenoOutFormat)
  
  else if (do_compare) then
    if (ImpuFile == 'none') then
      call chk_geno(DelFile, compare_file)
    else
      call chk_impute(DelFile, ImpuFile, compare_file)
    endif
  endif
  
  if (do_del .or. (do_compare .and. ImpuFile == 'none')) then
    deallocate(geno, IdV, SNP_names)
  endif
  
  
contains
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !  read in command line arguments
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine read_args()
    use sqa_fileIO, only: valid_formats
    integer :: nArg, i, x
    character(len=32) :: arg, argOption
    
    ! set defaults
    GenoFile = 'Geno'
    GenoInFormat = 'SEQ'
    GenoOutFile = 'geno_w_dels'
    GenoOutFormat = 'SEQ'
    DelFile = 'deletions.txt'
    ImpuFile = 'imputation_edits.txt'
    compare_file = 'imputation_compare.txt'
    do_del = .FALSE.
    do_compare = .FALSE.
    prop_del = 0
    
     nArg = command_argument_count()
    if (nArg > 0) then
      i = 0
      do x = 1, nArg
        i = i+1
        if (i > nArg)  exit
        call get_command_argument(i, arg)
        
        select case (arg)
          case ('--help')
            call print_help()
            stop   
            
          case ('--geno', '--geno-in')  
            i = i+1
            call get_command_argument(i, GenoFile)
              
          case ('--informat', '--inFormat')
            i = i+1
            call get_command_argument(i, GenoInFormat)
            if (.not. any(valid_formats == GenoInFormat)) then
              print *, 'ERROR: informat must be one of: ', valid_formats
              stop
            endif
            
          case ('--del')
            do_del = .TRUE.
            i = i+1
            call get_command_argument(i, argOption)
            read(argOption, *)  prop_del
            if (prop_del < 0d0 .or. prop_del > 1d0)  stop 'prop_del must be >=0 & <=1'
                   
         case ('--delfile')
            i = i+1
            call get_command_argument(i, DelFile)
         
          case ('--compare')
            do_compare = .TRUE.
            call get_command_argument(i, argOption)
            if (Len_Trim(argOption) > 0 .and. argOption(1:2)/="--") then   ! optional argument
              i = i+1
              compare_file = argOption
            endif
            
          case ('--impufile')
            i = i+1
            call get_command_argument(i, ImpuFile)
                       
          case ('--out', '--geno-out')
            i = i+1
            call get_command_argument(i, GenoOutFile)
            
          case ('--outformat', '--outFormat')
            i = i+1
            call get_command_argument(i, GenoOutFormat)
            if (.not. any(valid_formats == GenoOutFormat)) then
              print *, 'ERROR: outFormat must be one of: ', valid_formats
              stop
            endif

          case default
            print '(2a, /)', 'Unrecognised command-line option: ', arg
            print *, 'see --help for valid options'
            stop

        end select
      end do
    endif
    
    if (.not. do_del .and. .not. do_compare)  stop 'please choose del or compare'
    if (do_del .and. do_compare)  stop 'please choose *either* del or compare'
  
  
  end subroutine read_args
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine print_help()
    print '(a)',  'Make imputation test data or run comparison'
    print '(a, /)', 'command-line options:'
    print '(a)',    '  --help                print usage information and exit'
    print '(a)',    '  --geno <filename>     file with genotype data, file extension will be added based on',&
                    '                         --informat. Default: Geno'
    print '(a)',    '  --informat <x>        SEQ: no header, 0/1/2/-9, IDs in column 1; .txt (default)', &  
                    '                        PED: no header, 11/12/22/00, IDs in column 2 of 6 ', &
                    '                          non-SNP columns; .ped (+.map)', &
                    '                        RAW: header, 0/1/2/NA, IDs in column 2 of 6 non-SNP columns; .raw', &
                    '                        LMT: no header, 0/1/2 without spacing, IDs in separate file;', &
                    '                          .geno + .id'  
    print '(a)',    '  --del <proportion>    Proportion of genotype data to set to missing. Current genotypes',&
                    '                          are ignored, and may already be missing'
    print '(a)',    '  --geno-out <filename> Genotype file with deletions, to run imputation on. Default geno_w_dels'
    print '(a)',    '  --delfile <filename>  Logfile from deletions, storing original genotypes. Output',&
                    '                          if --del, input if --compare. Default deletions.txt'      
    print '(a)',    '  --compare <filename>  Run comparisons between --delfile and --impufile, output to <filename>.',&
                    '                          Default imputation_compare.txt'
    print '(a)',    '  --impufile <filename> Logfile from imputation to use in comparison. If "none", comparison',&
                    '                          is run between --delfile and --geno. Default imputation_edits.txt'
    print '(a)',    '  --outformat <x>       same options as for --informat. Default: SEQ'                    
    
  end subroutine print_help
    
  
end program main

!===============================================================================

subroutine del_snps(prop_del, DelFile)
  use stats
  use global_variables
  implicit none
  real, intent(IN) :: prop_del
  character(len=*), intent(IN) :: DelFile
  integer :: n_del, x,i,l
  integer, allocatable :: del_rows(:), del_cols(:)
  
  n_del = nInt(prop_del * nInd * nSnp)
  
  write(*, '("# deletions: ", f7.4, 2(" x ", i7), " ==> ", i7 )')  prop_del, nInd, nSnp, n_del
  
  allocate(del_rows(n_del), del_cols(n_del))
  del_rows = runif_int(n_del, first=1, last=nInd)  ! TODO: double check orientation
  del_cols = runif_int(n_del, first=1, last=nSnp)
  
  
  open(unit=4, file=trim(DelFile))
  write(4, '(4(a9,2x), 29x,2(a5,2x))') 'snp_index', 'snp_name', 'id_index', 'id_name', 'g_in', 'g_out'
  
  do x=1,n_del
    i = del_rows(x)  ! individual
    l = del_cols(x)  ! snp
    write(4, '(i9,2x,a9,2x,i9,2x,a40,2x,2i5)') l, SNP_names(l), i, IdV(i), Geno(i,l), -1
    Geno(i,l) = -1
  enddo
  
  close(4)

end subroutine del_snps

!===============================================================================

subroutine chk_impute(DelFile, ImpuFile, OutFile)
  use sqa_fileIO
  use global_variables
  implicit none
  
  character(len=*), intent(IN) :: DelFile, ImpuFile, OutFile
  integer :: n_dels, n_impu, d,x,i,l, counts(-1:2, -1:2), dumI, g_act, g_imp
  integer, allocatable :: imputations(:,:), dels(:,:), imp_i(:), imp_g(:), del_i(:), del_g(:)
  character(len=40) :: dumC(2)
  real :: dumR
  
   
  call CheckFile(trim(DelFile)) 
  n_dels = FileNumRow(trim(DelFile)) -1  ! 1st row = header
  
  call CheckFile(trim(ImpuFile)) 
  n_impu = FileNumRow(trim(ImpuFile)) -1  ! 1st row = header
  
  print *, 'n_dels: ', n_dels, '  n_impu: ', n_impu
  
  ! read imputation file
  allocate(imputations(n_impu,4))   ! l, i, g_in, g_out
  open (5, file=trim(ImpuFile), action='read')
    read(5,*)  ! header
    do x = 1,n_impu
      read (5,*) imputations(x,1), dumC(1), dumR, imputations(x,2), dumC(2), imputations(x,3), imputations(x,4)   
    enddo
  close(5)
  
  allocate(dels(n_dels,3))     ! l, i, g_act
  open (5, file=trim(DelFile), action='read')
    read(5,*)  ! header
    do d = 1, n_dels
      read (5,*) dels(d,1), dumC(1), dels(d,2), dumC(2), dels(d,3)
    enddo
  close(5)

  ! compare imputations to deletions file
  counts = 0
  nSnp = MAX(maxval(imputations(:,1)), maxval(dels(:,1)))
  do l=1,nSnp
    if (mod(l,500)==0)   print *, 'SNP ', l 
    imp_i = PACK(imputations(:,2), MASK=imputations(:,1)==l .and. imputations(:,3)==-1)
    imp_g = PACK(imputations(:,4), MASK=imputations(:,1)==l .and. imputations(:,3)==-1)
    del_i = PACK(dels(:,2), MASK=dels(:,1)==l)
    del_g = PACK(dels(:,3), MASK=dels(:,1)==l)
    do d=1,SIZE(del_i)
      do x=1,SIZE(imp_i)
        if (del_i(d)==imp_i(x)) then
          g_act = del_g(d)
          g_imp = imp_g(x)        
          counts(g_act, g_imp) = counts(g_act, g_imp) +1
        endif
      enddo
    enddo
  enddo
  
  print *, 'counts:'
  do x=-1,2
    print *, counts(x,:)
  enddo
  
  open (101, file=trim(OutFile), action='write')
    write(101, '(2a7,15x,a5)')  'g_act', 'g_imp', 'count'
    do i=-1,2
      do x=-1,2
        write(101, '(2i7, i20)') x,i, counts(x,i)
      enddo
    enddo
  close(101)

end subroutine chk_impute

!===============================================================================

subroutine chk_geno(DelFile, OutFile)
  use sqa_fileIO
  use global_variables
  implicit none
  
  character(len=*), intent(IN) :: DelFile, OutFile
  integer :: n_dels, d,x,i,l, counts(-1:2, -1:2), g_act, g_imp
  character(len=40) :: dumC(2)
     
  call CheckFile(trim(DelFile)) 
  n_dels = FileNumRow(trim(DelFile)) -1  ! 1st row = header

  ! compare genotype array to deletions file
  counts = 0
  open (5, file=trim(DelFile), action='read')
    read(5,*)  ! header
    do d = 1, n_dels
      read (5,*) l, dumC(1), i, dumC(2), g_act
      ! find matching entry in imputation genotype data
      g_imp = geno(i,l)
      counts(g_act, g_imp) = counts(g_act, g_imp) +1
    enddo
  close(5)
  
  print *, 'counts:'
  do x=-1,2
    print *, counts(x,:)
  enddo
  
  open (101, file=trim(OutFile), action='write')
    write(101, '(2a7,15x,a5)')  'g_act', 'g_imp', 'count'
    do i=-1,2
      do x=-1,2
        write(101, '(2i7, i20)') x,i, counts(x,i)
      enddo
    enddo
  close(101)

end subroutine chk_geno

!===============================================================================