program misfit_adjoint
    ! to calculate adjoint source and misfit
    ! yanhuay@princeton.edu

#ifdef USE_MPI
    use mpi
#endif

    use seismo_parameters
    implicit none

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer, parameter :: NARGS = 5
    integer ::  ier,i,icomp,irec,send_data_tag,receive_data_tag
    real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: adj_master 
    character(len=MAX_STRING_LEN) :: data_names(MAX_DATA_NUM)
    character(len=MAX_STRING_LEN) :: data_names_comma_delimited
    character(len=MAX_STRING_LEN) :: arg(NARGS)
    character(len=MAX_STRING_LEN) :: input_dir
    character(len=MAX_FILENAME_LEN) :: filename
    real t1,t2
    character(len=MAX_STRING_LEN) :: measurement_types(MAX_MISFIT_TYPE)
    integer :: itype,ntype
    character, parameter :: delimiter='+'
    logical :: ex

#ifdef USE_MPI
    call MPI_INIT(ier)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ier)
    call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ier)
#else
    nproc = 1
    myrank = 0
#endif

    if (DISPLAY_DETAILS .and. myrank==0)&
        print *,"Running misfit_adjoint.exe on",NPROC,"processors"
    call cpu_time(t1)

    ! parse command line arguments
    if (command_argument_count() /= NARGS) then
        if (myrank == 0) then
            print *, 'USAGE:  mpirun -np NPROC bin/misfit_adjoint.exe ... ' 
            stop ' Please check command line arguments'
        endif
    endif

#ifdef USE_MPI
    call MPI_BARRIER(MPI_COMM_WORLD,ier)
#endif

    ! get input parameters
    do i = 1, NARGS
    call get_command_argument(i,arg(i), status=ier)
    enddo
    read(arg(1),*) compute_adjoint
    data_names_comma_delimited = arg(2)
    measurement_list=arg(3)
    misfit_type_list=arg(4)
    input_dir=arg(5)
    call split_string(data_names_comma_delimited,',',data_names,ndata)

    !! gloabl initialization
    misfit=0.0_CUSTOM_REAL

    ! loop over comp
    do icomp=1,ndata

    ! step 1 -- load data and syn
    call initialize(input_dir,adjustl(data_names(icomp)))

    ! exist data/syn
    if(nrec_proc > 0 ) then  ! non-zero file

        ! step 2 -- preprocessing data and syn 
        if(norm2(seism_obs(:,:))<SMALL_VAL) cycle
        call process_data_all(seism_obs,'obs')
        if(norm2(seism_syn(:,:))<SMALL_VAL) cycle
        call process_data_all(seism_syn,'syn')

        if(compute_adjoint .and. DISPLAY_DETAILS) then
            write(filename,'(a)') &
                trim(input_dir)//'/DATA_obs/'//trim(adjustl(data_names(icomp)))//'_processed.bin'
            print*,'SAVE processed seismic_obs -- ',trim(filename)
            open(unit=IOUT,file=trim(filename),status='unknown',form='unformatted',iostat=ier)
            if (ier /= 0) then
                print*, 'Error: could not open model file: ',trim(filename)
                stop 'Error reading neighbors external mesh file'
            endif
            write(IOUT) seism_obs
            close(IOUT)
            write(filename,'(a)') &
                trim(input_dir)//'/DATA_syn/'//trim(adjustl(data_names(icomp)))//'_processed.bin'
            print*,'SAVE processed seismic_syn -- ',trim(filename)
            open(unit=IOUT,file=trim(filename),status='unknown',form='unformatted',iostat=ier)
            if (ier /= 0) then
                print*, 'Error: could not open model file: ',trim(filename)
                stop 'Error reading neighbors external mesh file'
            endif
            write(IOUT) seism_syn
            close(IOUT)
        endif

        ! step 3 -- misfit-adjoint
        call split_string(measurement_list,delimiter,measurement_types,ntype)
        if(ntype>mtype) then
            print*,'number of measurement_types exceeds number of measurement_weight: ',&
                ntype, mtype
            stop
        endif
        if (ntype<1) then 
            print*, 'measurement type should be at least 1'
            stop
        endif
        ! sum of the weight to be 1
        !measurement_weight=measurement_weight/sum(measurement_weight(1:ntype))

        do itype=1, ntype
       ! initialize 
        seism_adj_AD = 0.0_CUSTOM_REAL
        seism_adj_DD = 0.0_CUSTOM_REAL
        num_AD=0
        num_DD=0
        misfit_AD=0.0
        misfit_DD=0.0
        mean_AD=0.0
        mean_DD=0.0
        var_AD=1.0
        var_DD=1.0

        if (measurement_weight(itype) > 1e-8) then
        ! absolute mwasurement  
        if(index(misfit_type_list,'AD')>0) then
            write(filename,'(a,i6.6,a)') &
                trim(input_dir)//'/proc',myrank,'_misfit_'//trim(measurement_types(itype))//'_AD'
            OPEN (IOUT,FILE=trim(filename),status='unknown',iostat=ier) 
            call Absolute_diff(trim(measurement_types(itype)))
            close(IOUT)

            ! distribution estimated variance
            if(uncertainty .and. num_AD>1) then
            ! load measurement_AD
            if(trim(measurement_types(itype))=='EP') then 
                allocate(measurement_AD(2*num_AD))
            else
                allocate(measurement_AD(num_AD))
            endif
            OPEN (UNIT=IIN, FILE=filename,iostat=ier)
            if(ier/=0)  stop 'fail to open DD misfit file'
            read(IIN,*) measurement_AD
            close(IIN)
            ! estimate mean and variance 
            mean_AD=sum(measurement_AD)/num_AD
            ! write/read variance into/from file
            write(filename,'(a,i6.6,a)') &
                    trim(input_dir)//'/proc',myrank,'_variance_'//trim(measurement_types(itype))//'_AD'
                ex=.false.
                inquire (file=trim(filename), exist=ex)
                if(ex) then
                    ! use provided var_AD
                    OPEN (IIN,FILE=trim(filename),status='unknown',iostat=ier)
                    read (IIN,*) var_AD
                    close(IIN)
                else
                    ! save estimated var_AD
                    var_AD=sum((measurement_AD-mean_AD)**2)/(num_AD-1)
                    OPEN (IOUT,FILE=trim(filename),status='unknown',iostat=ier)
                    write(IOUT,*) var_AD
                    close(IOUT)
                endif
            deallocate(measurement_AD)
        endif !! unvertainty
            ! misfit and adj
            print*, 'misfit before weight:',misfit_AD, var_AD, max(num_AD,1),measurement_weight(itype) 
            print*, 'misfit after weight:',misfit_AD/var_AD/max(num_AD,1),misfit_AD/var_AD/max(num_AD,1)*measurement_weight(itype)
            misfit_AD=misfit_AD/var_AD/max(num_AD,1)*measurement_weight(itype) 
            seism_adj_AD=seism_adj_AD/var_AD/max(num_AD,1)*measurement_weight(itype)

            if(DISPLAY_DETAILS) then
                print*
                print*, trim(adjustl(data_names(icomp))),' comp'
                print*, 'misfit_',trim(trim(measurement_types(itype))),'_AD'
                print*, 'measurement_weight : ', measurement_weight(itype)
                print*, 'total number of measurements :', num_AD
                print*, 'misfit_AD :',misfit_AD 
                print*, 'min/max of adj : ',&
                minval(seism_adj_AD(:,:)), maxval(seism_adj_AD(:,:))
            endif
        endif !!AD

        ! differential measurement
        if (nrec_proc>1 .and. index(misfit_type_list,'DD')>0) then
            write(filename,'(a,i6.6,a)') &
                trim(input_dir)//'/proc',myrank,'_misfit_'//trim(measurement_types(itype))//'_DD'
            OPEN (IOUT,FILE=trim(filename),status='unknown',iostat=ier)
            call Relative_diff(input_dir,adjustl(data_names(icomp)),trim(measurement_types(itype)))
            close(IOUT)

            ! distribution estimated variance
            if(uncertainty .and. num_DD>1) then
            ! load measurement_DD
            if(trim(measurement_types(itype))=='EP') then
                allocate(measurement_DD(2*num_DD))
            else
                allocate(measurement_DD(num_DD))
            endif
            OPEN (UNIT=IIN, FILE=filename,iostat=ier)
            if(ier/=0)  stop 'fail to open DD misfit file'
            read(IIN,*) measurement_DD
            close(IIN)
            ! estimate mean and variance 
            mean_DD=sum(measurement_DD)/num_DD
            ! write/read variance into/from file
            write(filename,'(a,i6.6,a)') &
                    trim(input_dir)//'/proc',myrank,'_variance_'//trim(measurement_types(itype))//'_DD'
                ex=.false.
                inquire (file=trim(filename), exist=ex)
                if(ex) then
                    ! use provided var_DD
                    OPEN (IIN,FILE=trim(filename),status='unknown',iostat=ier)
                    read (IIN,*) var_DD
                    close(IIN)
                else
                    ! save estimated var_DD
                    var_DD=sum((measurement_DD-mean_DD)**2)/(num_DD-1)
                    OPEN (IOUT,FILE=trim(filename),status='unknown',iostat=ier)
                    write(IOUT,*) var_DD
                    close(IOUT)
                endif
            deallocate(measurement_DD)    
        endif !! uncertainty
            ! misfit and adj
            misfit_DD=misfit_DD/var_DD/max(num_DD,1)*measurement_weight(itype)
            seism_adj_DD=seism_adj_DD/var_DD/max(num_DD,1)*measurement_weight(itype)

           if(DISPLAY_DETAILS) then
                print*
                print*, trim(adjustl(data_names(icomp))),' comp'
                print*, 'misfit_',trim(trim(measurement_types(itype))),'_DD'
                print*, 'measurement_weight : ', measurement_weight(itype)
                print*, 'total number of measurements :', num_DD
                print*, 'misfit_DD :', misfit_DD
                print*, 'min/max of adj : ',&
                minval(seism_adj_DD(:,:)), maxval(seism_adj_DD(:,:))
            endif
        endif !! DD
    endif !  weight>1e-8

        ! combine misfit and adj for AD and DD
        misfit=misfit+misfit_AD + misfit_DD
        if(compute_adjoint) then 
            seism_adj=seism_adj + seism_adj_AD + seism_adj_DD
            call process_adj_all()
        endif

        enddo ! itype

        ! finalize
        call finalize(input_dir,adjustl(data_names(icomp)))

    endif ! nrec_proc>0
    enddo ! icomp

    ! step 5 -- save misfit
    write(filename,'(a,i6.6,a)') &
        trim(input_dir)//'/proc',myrank,'_misfit.dat'
    OPEN (IOUT, FILE=trim(filename),status='unknown',iostat = ier)
    if(ier/=0) then
        print*,'Error opening data misfit file: ',trim(filename)
        stop
    else
        write(IOUT,*) misfit
    endif
    close(IOUT)
    if(DISPLAY_DETAILS .and. compute_adjoint .and. nrec_proc>0) then
        print*
        print*,'SAVE misfit -- ',trim(filename)
        print*,'myrank=',myrank,' final misfit_',trim(measurement_list), &
            '_',trim(misfit_type_list), '= ',misfit
    endif

    call cpu_time(t2)
    if(DISPLAY_DETAILS .and. compute_adjoint .and. myrank==0) &
        print *,'Computation time with CPU:',t2-t1

#ifdef USE_MPI
    ! stop all the processes and exit
    call MPI_FINALIZE(ier)
#endif

end program misfit_adjoint
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine initialize(directory,data_name)
    use seismo_parameters
    implicit none
    integer :: ier,iker,imod,irec,itime
    integer :: filesize, nrec_obs, nrec_syn
    logical :: ex_obs, ex_syn
    character(len=MAX_FILENAME_LEN) :: char_i, filen, filename_obs, filename_syn
    character(len=MAX_STRING_LEN) :: directory
    character(len=MAX_STRING_LEN) :: data_name
    character(len=MAX_FILENAME_LEN) :: stf_filename,filename
    logical :: ex_stf
    real(kind=CUSTOM_REAL), dimension(:,:),allocatable :: temp
    integer :: irow, ncolum
    REAL(KIND=CUSTOM_REAL), DIMENSION(:,:), ALLOCATABLE :: noise

    !! data file format
    filen='empty'
    select case(solver)
    case('specfem2D')   ! single file
        if(myrank==0 .and. data_name == 'x') filen='Ux_file_single.su'
        if(myrank==0 .and. data_name == 'y') filen='Uy_file_single.su'
        if(myrank==0 .and. data_name == 'z') filen='Uz_file_single.su'
        if(myrank==0 .and. data_name == 'p') filen='Up_file_single.su'

    case('specfem3D')  
        write(char_i, '(I5)') myrank        ! convert integer to char
        if(data_name == 'x') write(filen,'(a)') trim(adjustl(char_i))//'_dx_SU'
        if(data_name == 'y') write(filen,'(a)') trim(adjustl(char_i))//'_dy_SU'
        if(data_name == 'z') write(filen,'(a)') trim(adjustl(char_i))//'_dz_SU'
        if(data_name == 'p') write(filen,'(a)') trim(adjustl(char_i))//'_dp_SU'

    case default
        print*,'Currently only work for specfem2D and specfem3D solver'
        stop
    end select

    !! load data & syn
    nrec_proc=0
    write(filename_obs,'(a)') &
        trim(directory)//'/DATA_obs/'//trim(filen)
    write(filename_syn,'(a)') &
        trim(directory)//'/DATA_syn/'//trim(filen)
    !exist and size
    inquire (file=trim(filename_obs), exist=ex_obs)
    inquire (file=trim(filename_syn), exist=ex_syn)

    if(ex_obs .and. ex_syn) then ! exist
        inquire(file=trim(filename_obs), size=filesize)
        nrec_obs=filesize/(240+4*NSTEP) 
        inquire(file=trim(filename_syn), size=filesize)
        nrec_syn=filesize/(240+4*NSTEP)
        if(nrec_obs == nrec_syn) then
            nrec_proc=nrec_obs
            if(DISPLAY_DETAILS) then
                print*,'myrank=',myrank,' LOAD nrec_proc traces : ',nrec_proc
                print*,'seism_obs -- ',trim(filename_obs)
                print*,'seism_syn -- ',trim(filename_syn)
            endif
        else
            print*,'size for obs and syn file is not the same : ',nrec_obs, nrec_syn
            stop
        endif

        !! allocate 
        allocate(seism_obs(NSTEP,nrec_proc))
        allocate(seism_syn(NSTEP,nrec_proc))
        allocate(seism_adj(NSTEP,nrec_proc))
        allocate(seism_adj_AD(NSTEP,nrec_proc))
        allocate(seism_adj_DD(NSTEP,nrec_proc))
        allocate(noise(NSTEP,nrec_proc))
        allocate(st_xval(nrec_proc))
        allocate(st_yval(nrec_proc))
        allocate(st_zval(nrec_proc))
        allocate(win_start(nrec_proc))
        allocate(win_end(nrec_proc))
        allocate(trace_norm(nrec_proc))
        allocate(time(NSTEP))
        allocate(dis_sr(nrec_proc))
        !allocate(which_proc_receiver(NPROC))

        ! initialization
        seism_obs = 0.0_CUSTOM_REAL
        seism_syn = 0.0_CUSTOM_REAL
        seism_adj = 0.0_CUSTOM_REAL
        seism_adj_AD = 0.0_CUSTOM_REAL
        seism_adj_DD = 0.0_CUSTOM_REAL
        win_start=0.0
        win_end=0.0
        event_norm=1.0
        trace_norm=1.0
        do itime=1,NSTEP
        time(itime)=(itime-1)*deltat+t0
        enddo
        ! which_proc_receiver=0

        ! allocate obs and syn
        call readSU(filename_obs,seism_obs)
        ! generate random noise to obs between 0 and 1
        call random_number(noise)
        ! [-1 1]*noise-level   
        noise=(noise-0.5)*2*maxval(abs(seism_obs))*noise_level
        ! add noise 
        seism_obs=seism_obs+noise
        call readSU(filename_syn,seism_syn)
        do irec=1,nrec_proc
        dis_sr(irec)=sqrt((st_xval(irec)-x_source)**2 &            
            +(st_yval(irec)-y_source)**2 &                                            
            +(st_zval(irec)-z_source)**2)
        enddo

        if(DISPLAY_DETAILS) then
            print*,'Min / Max of seism_obs : ',&
                minval(seism_obs(:,:)),maxval(seism_obs(:,:))
            print*,'Min / Max of seism_syn : ',&
                minval(seism_syn(:,:)),maxval(seism_syn(:,:))
        endif
        deallocate(noise)
    endif ! exist

end subroutine initialize
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine process_data_all(seism,tag)
    use seismo_parameters
    implicit none
    integer :: ier,irec
    real(kind=CUSTOM_REAL) :: seism(NSTEP,nrec_proc)
    character(len=3) :: tag
    real(kind=CUSTOM_REAL) :: trace(NSTEP)
    real(kind=CUSTOM_REAL) :: wtr
    real(kind=CUSTOM_REAL) :: taper_len
    integer :: itime
    integer :: NA
    real(kind=CUSTOM_REAL), dimension(:), allocatable :: tas

    !! Event normalization
    if(EVENT_NORMALIZE) then
        wtr=norm2(seism(:,:))
        if(wtr<SMALL_VAL) stop
        if (DISPLAY_DETAILS) print*,'Event normalization ',tag, wtr
        seism=seism/wtr
        if(tag=='syn') event_norm=wtr
    endif

    !! trace-wise processing
    do irec = 1,nrec_proc 
    trace(:) = seism(:,irec)
    tstart=0.0
    tend=(NSTEP-1)*deltat

    wtr=norm2(trace(1:NSTEP))
    if(wtr>SMALL_VAL) then  ! non-zero trace
        !! trace normalization 
        if(TRACE_NORMALIZE) then
            if (DISPLAY_DETAILS .and. irec==1) print*,'Trace normalization ', tag
            trace=trace/wtr 
            if(tag=='syn') trace_norm(irec)=wtr
        endif
        ! WT filtering
        if( Wscale .gt. 0) then
            call WT(trace,NSTEP,Wscale,NA)
        endif
        !! mute near offset 
        if(MUTE_NEAR .and. dis_sr(irec)<=offset_near ) then
            trace(1:NSTEP) = 0.0
            tstart=0.0
            tend=0.0
        endif
        !! muter far offset
        if(MUTE_FAR .and. dis_sr(irec)>=offset_far) then
            trace(1:NSTEP) = 0.0
            tstart=0.0
            tend=0.0
        endif
        !! laplace damping spatially and temporally 
        if (DAMPING) then
            ! spatial
            trace=trace*exp(-(dis_sr(irec)*lambda_x)**2)
            ! temporal
            trace(1:NSTEP)=trace(1:NSTEP)*exp(-(time(1:NSTEP)*lambda_t)**2)
        endif
        ! time-domain window using slopes 
        if(TIME_WINDOW) then
            tstart=dis_sr(irec)/VEL_TOP
            tend=dis_sr(irec)/VEL_BOT
            if(tend>tstart) then
                taper_len = taper_percentage/2.0 * (tend-tstart)
                tstart=max(tstart-taper_len,0.0)
                tend=min(tend+taper_len,(NSTEP-1)*deltat)
            endif
        endif ! window
        ! save
        win_start(irec)=tstart
        win_end(irec)=tend
    endif ! non-zero trace
    seism(:,irec)=trace(:)
    enddo ! irec

end subroutine process_data_all
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine process_adj_all()
    use seismo_parameters
    implicit none
    integer :: irec
    real(kind=CUSTOM_REAL) :: adj(NSTEP), adj_vel(NSTEP), adj_acc(NSTEP)

    !! post-processing of seism_adj
    do irec=1,nrec_proc
    adj(:)=seism_adj(:,irec)
    if(norm2(adj(:))<SMALL_VAL) cycle
    ! trace normalize
    if(TRACE_NORMALIZE) adj=adj/trace_norm(irec)
    ! seismotype
    if(DISPLAY_DETAILS .and. irec==1) print*, trim(seismotype) ,' adjoint'
    if(trim(seismotype) == "velocity") then
        call compute_vel(adj,NSTEP,deltat,NSTEP,adj_vel)
        seism_adj(:,irec)= - adj_vel(:)
    elseif(trim(seismotype) == "acceleration") then
        call compute_acc(adj,NSTEP,deltat,NSTEP,adj_acc)
        seism_adj(:,irec)=adj_acc(:)
    endif
    enddo
    ! event normalize
    if(EVENT_NORMALIZE) seism_adj=seism_adj/event_norm

end subroutine process_adj_all
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine process_adj_trace(trace,dis)
    use seismo_parameters
    implicit none
    real(kind=CUSTOM_REAL) :: trace(NSTEP)
    real(kind=CUSTOM_REAL) :: dis
    real(kind=CUSTOM_REAL) :: wtr
    integer :: NA
    integer :: itime
    real(kind=CUSTOM_REAL) :: tas(NSTEP)
    real(kind=CUSTOM_REAL), dimension(:),allocatable :: stf_reverse

        !! laplace damping spatially and temporally 
        if (DAMPING) then
            ! spatial
            trace=trace*exp(-(dis*lambda_x)**2)
            ! temporal
            trace(1:NSTEP)=trace(1:NSTEP)*exp(-(time(1:NSTEP)*lambda_t)**2)
        endif
        !! mute near offset 
        if(MUTE_NEAR .and. dis<=offset_near) then
            trace(1:NSTEP) = 0.0
        endif
        !! mute far offset
        if(MUTE_FAR .and. dis>=offset_far) then                                                      
            trace(1:NSTEP) = 0.0
        endif
        ! WT filtering
        if( Wscale .gt. 0) then
            call WT(trace,NSTEP,Wscale,NA)
        endif

end subroutine process_adj_trace
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  `  
subroutine Absolute_diff(measurement_type)
    use seismo_parameters
    implicit none
    character(len=2) :: measurement_type
    integer :: irec
    real(kind=CUSTOM_REAL) :: misfit_trace
    real(kind=CUSTOM_REAL) :: d(NSTEP),s(NSTEP),adj(NSTEP)
    real(kind=CUSTOM_REAL) :: wtr
    integer :: nlen,num,num_measure

    ! initialization
    d = 0.0_CUSTOM_REAL
    s = 0.0_CUSTOM_REAL
    misfit_AD=0.0
    num_measure=0

    do irec=1,nrec_proc
    ! get data 
    d(:)=seism_obs(:,irec)
    s(:)=seism_syn(:,irec)

    !! misfit and adjoint evaluation
    ! initialization
    adj(:) = 0.0
    misfit_trace=0.0
    num = 0

    ! window info
    tstart = win_start(irec)
    tend= win_end(irec)
    if((tend-tstart)<min_window_len) cycle
    
    !!  evaluate misfit and adj 
    call misfit_adj_AD(measurement_type,d,s,NSTEP,&
        deltat,f0,tstart,tend,taper_percentage,taper_type,&
        compute_adjoint, &
        adj,num,misfit_trace)
    num_AD = num_AD + num 
    num_measure=num_measure+1
    misfit_AD = misfit_AD + misfit_trace
    if(DISPLAY_DETAILS) then 
        print*
        print*,'irec=',irec, 'misfit_',measurement_type,'_AD=',misfit_trace
    endif

    if(compute_adjoint) then 
        call process_adj_trace(adj,dis_sr(irec))
        seism_adj_AD(:,irec) = seism_adj_AD(:,irec) + adj(:)
    endif
    enddo ! irec
    if(DISPLAY_DETAILS) print*, 'Total number of AD measurements :', num_measure
end subroutine Absolute_diff
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Relative_diff(input_dir,data_name,measurement_type)
    use seismo_parameters
    implicit none
    integer :: irec,jrec,ier
    character(len=MAX_STRING_LEN) :: input_dir
    character(len=MAX_STRING_LEN) :: data_name
    character(len=MAX_FILENAME_LEN) :: filename
    logical :: ex
    real(kind=CUSTOM_REAL) :: misfit_trace
    real(kind=CUSTOM_REAL) :: d(NSTEP),s(NSTEP),adj(NSTEP)
    real(kind=CUSTOM_REAL) :: d_ref(NSTEP),s_ref(NSTEP),adj_ref(NSTEP)
    real(kind=CUSTOM_REAL) :: dis_sr1, dis_sr2, dis_rr
    real(kind=CUSTOM_REAL) :: cc_max_obs
    integer :: nlen,num,num_measure=0
    character(len=2) :: measurement_type

    ! initialization
    d = 0.0_CUSTOM_REAL
    s = 0.0_CUSTOM_REAL
    d_ref = 0.0_CUSTOM_REAL
    s_ref = 0.0_CUSTOM_REAL
    cc_max_obs = 0.0
    misfit=0.0

    allocate(is_pair(nrec_proc,nrec_proc))
    is_pair=0

    write(filename,'(a)') trim(input_dir)//'DATA_obs/'//trim(data_name)//'.similarity.dat'
    ! if existence, read, otherwise, write 
    inquire (file=trim(filename), exist=ex)
    OPEN (UNIT=IIN, FILE=filename,iostat=ier)
    do while(ier==0)
    read(IIN,*,iostat=ier) irec,jrec,cc_max_obs,is_pair(irec,jrec)
    enddo

    ! loop over master trace
    do irec=1,nrec_proc

    ! get data 
    d(:)=seism_obs(:,irec)
    s(:)=seism_syn(:,irec)

    dis_sr1=dis_sr(irec)

    ! window info
    tstart = win_start(irec)
    tend = win_end(irec)
    if((tend-tstart)<min_window_len) cycle

    ! loop over reference trace
    do jrec=irec+1,nrec_proc
    ! get data 
    d_ref(:)=seism_obs(:,jrec)
    s_ref(:)=seism_syn(:,jrec)

    ! window info
    tstart_ref = win_start(jrec)        
    tend_ref= win_end(jrec)
    if((tend_ref-tstart_ref)<min_window_len) cycle

    if(.not. ex) then
        cc_max_obs=0.0
        dis_rr=sqrt((st_xval(jrec)-st_xval(irec))**2 &
            +(st_yval(jrec)-st_yval(irec))**2 &
            +(st_zval(jrec)-st_zval(irec))**2)
        if(dis_rr<=DD_max .and. dis_rr>=DD_min)&
            call CC_similarity(d,d_ref,NSTEP,deltat,&
            tstart,tend,tstart_ref,tend_ref,&
            taper_percentage,taper_type,&
            cc_max_obs)

        if(cc_max_obs>cc_threshold)  is_pair(irec,jrec) = 1
    endif !! ex

    if(DISPLAY_DETAILS .and. compute_adjoint) then
        print*      
        print*,'pair irec=',irec, 'jrec=',jrec           
        print*,'window -- ',tstart,tend,tstart_ref,tend_ref
        print*, 'rr distance dis_rr, DD_min/DD_max: ',dis_rr, DD_min, DD_max 
        print*, 'cc similarity -- ', cc_max_obs 
        print*,'is_pair : ',is_pair(irec,jrec)
    endif 

    if(is_pair(irec,jrec)==1) then
        ! initialization
        adj = 0.0
        adj_ref = 0.0
        num = 0
        misfit_trace=0.0

        dis_sr2=dis_sr(jrec)

        ! number of double difference measurements
        call misfit_adj_DD(measurement_type,d,d_ref,s,s_ref,NSTEP,deltat,f0,&
            tstart,tend,tstart_ref,tend_ref,taper_percentage,taper_type,&
            compute_adjoint, &
            adj,adj_ref,num,misfit_trace)
            num_DD = num_DD + num
            num_measure=num_measure+1
            misfit_DD = misfit_DD + misfit_trace
        if(DISPLAY_DETAILS) then 
            print* 
            print*,'irec=',irec,'jrec=',jrec, 'misfit_',measurement_type,'_DD=',misfit_DD
        endif
    
        if(compute_adjoint) then
            call process_adj_trace(adj,dis_sr1)
            call process_adj_trace(adj_ref,dis_sr2)
            if(DISPLAY_DETAILS) then
                print*, 'Min/Max of adj :',minval(adj(:)),maxval(adj(:))
                print*, 'Min/Max of adj_ref:',minval(adj_ref(:)),maxval(adj_ref(:))
            endif
        
            ! sum of adjoint source over pair
            seism_adj_DD(:,irec) = seism_adj_DD(:,irec) +  adj(:)
            seism_adj_DD(:,jrec) = seism_adj_DD(:,jrec) + adj_ref(:)

        endif ! compute_adjoint

        !! save waveform similarity
        if(.not. ex) write(IIN,'(2I5,1e15.5,I5)') irec,jrec,cc_max_obs,is_pair(irec,jrec)

    else
        write(IOUT,*) 0.0
    endif  ! is_pair

    enddo  ! jrec 
    enddo ! irec 
    close(IIN) ! close similarity file
    if(DISPLAY_DETAILS) print*, 'Total number of DD measurements :', num_measure
    deallocate(is_pair)

end subroutine Relative_diff
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine finalize(directory,data_name)
    use seismo_parameters
    implicit none
    character(len=MAX_STRING_LEN) :: data_name
    character(len=MAX_FILENAME_LEN) :: char_i,filen,filename
    character(len=MAX_STRING_LEN) :: directory

    !! write SU adj
    if(compute_adjoint) then
        !! adjoint file format
        select case(solver)
        case('specfem2D')
            if(myrank==0 .and. data_name == 'x') filen='Ux_file_single.su.adj'
            if(myrank==0 .and. data_name == 'y') filen='Uy_file_single.su.adj'
            if(myrank==0 .and. data_name == 'z') filen='Uz_file_single.su.adj'
            if(myrank==0 .and. data_name == 'p') filen='Up_file_single.su.adj'

        case('specfem3D')
            write(char_i, '(I5)') myrank        ! convert integer to char
            if(data_name == 'x') write(filen,'(a)') trim(adjustl(char_i))//'_dx_SU.adj'
            if(data_name == 'y') write(filen,'(a)') trim(adjustl(char_i))//'_dy_SU.adj'
            if(data_name == 'z') write(filen,'(a)') trim(adjustl(char_i))//'_dz_SU.adj'
            if(data_name == 'p') write(filen,'(a)') trim(adjustl(char_i))//'_dp_SU.adj'

        case default
            print*,'Currently only work for specfem2D and specfem3D solver'
            stop
        end select

        !! write adjoint source
        write(filename,'(a)') &
            trim(directory)//'/SEM/'//trim(filen)
        if(DISPLAY_DETAILS) then
            print*
            print*,'SAVE seism_adj -- ',trim(filename)
            if(DISPLAY_DETAILS) print*,'Min / Max of final seism_adj : ',&
                minval(seism_adj(:,:)),maxval(seism_adj(:,:))
        endif

        call writeSU(filename,seism_adj)
    endif

    deallocate(seism_obs)
    deallocate(seism_syn)
    deallocate(seism_adj)
    deallocate(seism_adj_AD)
    deallocate(seism_adj_DD)
    deallocate(st_xval)
    deallocate(st_yval)
    deallocate(st_zval)
    deallocate(win_start)
    deallocate(win_end)
    deallocate(trace_norm)
    deallocate(time)
    deallocate(dis_sr)

end subroutine finalize
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine readSU(filename,seism)
    use seismo_parameters
    implicit none
    integer :: ier,irec
    real(kind=CUSTOM_REAL) :: seism(NSTEP,nrec_proc)
    character(len=MAX_FILENAME_LEN) :: filename

    seism = 0.0_CUSTOM_REAL
    ! open(IIN,file=trim(filename),access='direct',recl=240+4*NSTEP,iostat = ier)
    open(IIN,file=trim(filename),status='old',form='unformatted',&
        access='direct',recl=240+4*NSTEP,iostat = ier)

    if (ier /= 0) then
        print*, 'Error: could not open data file: ',trim(filename)
        stop
    endif

    do irec = 1, nrec_proc
    read(IIN,rec=irec,iostat=ier) r4head, seism(:,irec)
    ! header info
    z_source=r4head(13) ! Source depth below surface (sdepth) 
    x_source=r4head(19) ! Source x coord (sx)
    y_source=r4head(20) ! Source y coord  (sy)
    st_zval(irec)=r4head(11) ! Receiver group elevation (gelev)
    st_xval(irec)=r4head(21) ! Receiver x coord (gx)
    st_yval(irec)=r4head(22) ! Receiver y coord (gy)
    ! header2=int(r4head(29), kind=2)
    if (DISPLAY_DETAILS .and. irec==1) print *, 'xs,ys,zs',&
        x_source,y_source,z_source
    if (DISPLAY_DETAILS .and. irec==1) print *, 'xr,yr,zr',&
        st_xval(irec),st_yval(irec),st_zval(irec)

    enddo

    close(IIN)

end subroutine readSU
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine writeSU(filename,seism)
    use seismo_parameters
    implicit none
    integer :: ier,irec,itime
    real(kind=CUSTOM_REAL) :: seism(NSTEP,nrec_proc)
    character(len=MAX_FILENAME_LEN) :: filename
    integer :: deltat_int2

    open(IOUT,file=trim(filename),status='unknown',access='direct',recl=4,iostat=ier)
    if (ier /= 0) then
        print*, 'Error: could not open data file: ',trim(filename)
        stop
    endif

    do irec = 1, nrec_proc
    ! write header
    write(IOUT,rec=(irec-1)*60+(irec-1)*NSTEP+1)  irec !receiver ID
    write(IOUT,rec=(irec-1)*60+(irec-1)*NSTEP+10) NINT(st_xval(irec)-x_source)  ! offset
    write(IOUT,rec=(irec-1)*60+(irec-1)*NSTEP+19) NINT(x_source)                ! source location xs
    write(IOUT,rec=(irec-1)*60+(irec-1)*NSTEP+20) NINT(y_source)                ! source location ys
    write(IOUT,rec=(irec-1)*60+(irec-1)*NSTEP+21) NINT(st_xval(irec))           ! receiver location xr
    write(IOUT,rec=(irec-1)*60+(irec-1)*NSTEP+22) NINT(st_yval(irec))           ! receiver location zr
    if (nrec_proc>1) write(IOUT,rec=(irec-1)*60+(irec-1)*NSTEP+48) SNGL(st_xval(2)-st_xval(1)) ! receiver interval
    header2(1)=0  ! dummy
    header2(2)=int(NSTEP, kind=2)
    write(IOUT,rec=(irec-1)*60+(irec-1)*NSTEP+29) header2
    header2(1)=NINT(deltat*1.0d6, kind=2) ! deltat (unit: 10^{-6} s)
    header2(2)=0  ! dummy
    write(IOUT,rec=(irec-1)*60+(irec-1)*NSTEP+30) header2

    ! the "60" in the following corresponds to 240 bytes header (note the
    ! reclength is 4 bytes)
    do itime = 1, NSTEP
    write(IOUT,rec=irec*60+(irec-1)*NSTEP+itime) sngl(seism_adj(itime,irec))
    enddo
    enddo ! irec

    close(IOUT)

end subroutine writeSU
