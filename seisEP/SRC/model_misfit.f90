program model_misfit
    !! yanhuay@princeton.edu
    !! To calculate model misfit

    use seismo_parameters
    implicit none

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
    integer, parameter :: NARGS = 5
    INTEGER :: i,iter,ier
    character(len=MAX_STRING_LEN) :: arg(NARGS)
    character(len=MAX_STRING_LEN) :: dir_target, dir_estimate
    character(len=MAX_STRING_LEN) :: directory
    character(len=MAX_FILENAME_LEN) :: FILENAME
    integer :: ipar,nspec_start,nspec_end
    REAL(kind=CUSTOM_REAL), DIMENSION(:,:,:,:), ALLOCATABLE :: rho, alpha, beta
    real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: rho_target,alpha_target,beta_target
    REAL(kind=CUSTOM_REAL) :: rho_norm=0.0,alpha_norm=0.0,beta_norm=0.0,sum_norm=0.0

    ! parse command line arguments
    if (command_argument_count() /= NARGS) then
        if (DISPLAY_DETAILS .and. myrank==0) then
            print *, 'USAGE:  mpirun -np NPROC bin/model_misfit.exe ...'
            stop ' Please check command line arguments'
        endif
    endif

    do i = 1, NARGS
    call get_command_argument(i,arg(i), status=ier)
    enddo

    read(arg(1),*) nproc
    read(arg(2),*) iter
    dir_target=arg(3)
    dir_estimate=arg(4)
    directory=arg(5)

    ! slice database file
    allocate(nspec_proc(nproc))
    nspec_proc=0  

    do myrank=0,nproc-1
    ! nspec
    write(filename,'(a,i6.6,a)') &
        trim(dir_target)//'/proc',myrank,'_'//trim(IBOOL_NAME) 
    open(IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
        print *,'Error: could not open database file:',trim(filename)
        stop 'Error opening _NSPEC_IBOOL file'
    endif
    read(IIN) nspec_proc(myrank+1)
    close(IIN)

    !if(DISPLAY_DETAILS .and. myrank==0) print*, 'nspec_proc=',nspec_proc(myrank+1)
    enddo

    nspec=sum(nspec_proc)
    !if(DISPLAY_DETAILS) print*,'NGLLX*NGLLY*NGLLZ*NSPEC*nmod:',NGLLX,NGLLY,NGLLZ,NSPEC,nmod

    !! local
    allocate(rho_target(NGLLX,NGLLY,NGLLZ,NSPEC))
    allocate(alpha_target(NGLLX,NGLLY,NGLLZ,NSPEC))
    allocate(beta_target(NGLLX,NGLLY,NGLLZ,NSPEC))
    allocate(rho(NGLLX,NGLLY,NGLLZ,NSPEC))
    allocate(alpha(NGLLX,NGLLY,NGLLZ,NSPEC))
    allocate(beta(NGLLX,NGLLY,NGLLZ,NSPEC))

    rho_target =  0.0_CUSTOM_REAL
    alpha_target =  0.0_CUSTOM_REAL
    beta_target =  0.0_CUSTOM_REAL
    rho =  0.0_CUSTOM_REAL
    alpha =  0.0_CUSTOM_REAL
    beta =  0.0_CUSTOM_REAL

    !! LOAD models
    do myrank=0,nproc-1
    nspec_start=sum(nspec_proc(1:myrank))+1
    nspec_end=sum(nspec_proc(1:myrank))+nspec_proc(myrank+1)

    do ipar=1,MAX_PAR_NUM

    ! load target model parameters (rho, alpha, beta)
    write(filename,'(a,i6.6,a)') &
        trim(dir_target)//'/proc',myrank,&
        '_'//trim(model_list(ipar))//'.bin'
   ! if(DISPLAY_DETAILS .and. myrank==0 .and. ipar==1) print*,'LOAD m_target ...'
   ! if(DISPLAY_DETAILS .and. myrank==0) print*, trim(filename)
    open(IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
        print*, 'Error: could not open m_target file: ',trim(filename)
        stop 'Error: could not open m_target file: '
    endif
    if(trim(model_list(ipar))== 'rho') then
        read(IIN) rho_target(:,:,:,nspec_start:nspec_end)
    elseif(trim(model_list(ipar))== 'vp') then
        read(IIN) alpha_target(:,:,:,nspec_start:nspec_end)
    elseif(trim(model_list(ipar))== 'vs') then
        read(IIN) beta_target(:,:,:,nspec_start:nspec_end)
    else
        print*,'Wrong model_list'
        stop
    endif
    close(IIN)

    ! load current model parameters (rho, alpha, beta)
    write(filename,'(a,i6.6,a)') &
        trim(dir_estimate)//'/proc',myrank,&
        '_'//trim(model_list(ipar))//'.bin'
    !if(DISPLAY_DETAILS .and. myrank==0 .and. ipar==1) print*,'LOAD m_estimate ...'
    !if(DISPLAY_DETAILS .and. myrank==0) print*, trim(filename)
    open(IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
        print*, 'Error: could not open m_estimate file: ',trim(filename)
        stop 'Error: could not open m_estimate file: '
    endif
    if(trim(model_list(ipar))== 'rho') then
        read(IIN) rho(:,:,:,nspec_start:nspec_end)
    elseif(trim(model_list(ipar))== 'vp') then
        read(IIN) alpha(:,:,:,nspec_start:nspec_end)
    elseif(trim(model_list(ipar))== 'vs') then
        read(IIN) beta(:,:,:,nspec_start:nspec_end)
    else
        print*,'Wrong model_list'
        stop
    endif
    close(IIN)

    enddo ! ipar

    enddo ! myrank

    !! norm 
    rho_norm = norm2(reshape(rho_target,(/1,NGLLX*NGLLY*NGLLZ*NSPEC/)) - &
        reshape(rho,(/1,NGLLX*NGLLY*NGLLZ*NSPEC/))) 
    alpha_norm = norm2(reshape(alpha_target,(/1,NGLLX*NGLLY*NGLLZ*NSPEC/)) - &
        reshape(alpha,(/1,NGLLX*NGLLY*NGLLZ*NSPEC/)))
    beta_norm =  norm2(reshape(beta_target,(/1,NGLLX*NGLLY*NGLLZ*NSPEC/)) - &
        reshape(beta,(/1,NGLLX*NGLLY*NGLLZ*NSPEC/)))
    sum_norm = rho_norm + alpha_norm + beta_norm

    !!! misfit evolution at iterations
    write(filename, "(a)") trim(directory)//'/misfit/model_misfit_hist.dat'
    OPEN (IOUT, FILE=trim(filename),status='unknown',POSITION='APPEND')
    write(IOUT,'(i,3e15.8)') iter,rho_norm,alpha_norm,beta_norm
    close(IOUT)

    if(DISPLAY_DETAILS) then
        print*, 'model misfit -- rho_norm=',rho_norm
        print*, 'model misfit -- alpha_norm=',alpha_norm
        print*, 'model misfit -- beta_norm=',beta_norm
        print*, 'model misfit -- sum_norm=',sum_norm
    endif

    deallocate(rho_target)
    deallocate(alpha_target)
    deallocate(beta_target)
    deallocate(rho)
    deallocate(alpha)
    deallocate(beta)
    deallocate(nspec_proc)

end program model_misfit
