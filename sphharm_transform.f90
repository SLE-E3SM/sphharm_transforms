    include 'spharmt.f90' ! Spherical harmonic transform module

    program sphharm_transform
    use spharmt

    character(*), parameter :: inputfolder  = 'OUTPUT_SLM/'
    character(*), parameter :: inputfile  = 'dsl_sphharm16'

    character(*), parameter :: outputfolder  = 'OUTPUT_SLM/'
    character(*), parameter :: outputfile  = 'dsl_grid16'
    type(sphere) :: spheredat                                   ! SH transform data to be passed to subroutines
    real :: radius

    integer, parameter :: nglv = 512             ! Number of GL points in latitude
    integer, parameter :: norder = 512           ! Max spherical harmonic degree/order
    complex, dimension(0:norder,0:norder) :: deltasllm
    real, dimension(nglv,2*nglv) :: deltaslxy

    radius = 6.371E6              ! Radius of the Earth (m)

    call spharmt_init(spheredat, 2*nglv, nglv, norder, radius) ! Initialize spheredat (for SH transform subroutines)
   
    write(*,*) 'Reading Spherical Harmonics'
    open(unit = 201, file = trim(inputfolder)//trim(inputfile), form = 'formatted', &
    & access = 'sequential', status = 'old')
    do m = 0,norder
      do n = 0,norder
        read (201,'(ES14.7E2, ES14.7E2)') vr, vi
        deltasllm(m,n) = complex(vr,vi)
      enddo
    enddo

    close(201)

    !Perform Transform
    call spec2spat(deltaslxy, deltasllm, spheredat)

    write(*,*) 'Writing Grid File'
    open(unit = 201, file = trim(outputfolder)//trim(outputfile), form = 'formatted', &
    & access = 'sequential', status = 'replace')
    write(201,'(ES16.9E2)') deltaslxy
    close(201)

end program sphharm_transform

