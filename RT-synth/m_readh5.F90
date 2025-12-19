module m_readh5
    use hdf5
    use RTsynth_module
    implicit none

    integer :: h5err
    character(len=256) :: in_file_h5
    integer(HID_T) :: file_id, dset_id, filespace
    double precision, allocatable :: temp_real(:,:,:)
contains

  subroutine readpdrfiles_h5
    implicit none

    ! locals
    integer :: i, j, k, isp, index
    integer :: access_flags
    integer(HSIZE_T), dimension(3) :: dims
    integer(HSIZE_T), dimension(4) :: dims3
    integer(HSIZE_T), dimension(1) :: dims_network, dims_cooarray
    character(len=128) :: network_name
    character(len=20), dimension(:), allocatable :: coo_array
    character(len=128), dimension(1) :: string_array
    character(len=20) :: file_name, file_index
    integer :: total_levels
    integer(8) :: str_size
    logical :: file_exists

    close(1);open(unit=1,file='chemfiles/species_'//trim(adjustl(chemsuf)),status='old')
    do i=1,nspec
      read(1,*) index,species(i)
      do j=1,coo
        if (coolant(j)%cname.eq.species(i)) coolant(j)%cspec=i
      enddo
    enddo

    call h5open_f(h5err)

    dims = [ctot, ctot, ctot]
    allocate(temp_real(ctot, ctot, ctot))

    in_file_h5 = trim(adjustl(filein)) // ".pdr.h5"
    inquire(file=trim(in_file_h5), exist=file_exists)
    if (.not. file_exists) then
      print *, "Error: Input file not found: ", trim(in_file_h5)
      stop
    end if

    access_flags = H5F_ACC_RDONLY_F
    call h5fopen_f(in_file_h5, access_flags, file_id, h5err)
    if (h5err /= 0) then
      print *, "Failed to open file: ", trim(in_file_h5)
      stop
    end if

    write(6,*) 'Reading abundances'
    call read_h5_one_r3d(file_id, dims, temp_real, "x");        pdr%x = temp_real
    call read_h5_one_r3d(file_id, dims, temp_real, "y");        pdr%y = temp_real
    call read_h5_one_r3d(file_id, dims, temp_real, "z");        pdr%z = temp_real
    call read_h5_one_r3d(file_id, dims, temp_real, "tgas");     pdr%Tgas = temp_real
    call read_h5_one_r3d(file_id, dims, temp_real, "tdust");    pdr%Tdust = temp_real
    call read_h5_one_r3d(file_id, dims, temp_real, "rho");      pdr%rho = temp_real

    do isp = 1, nspec
      write(file_index, '(I3.3)') isp
      file_name = "abundance" // trim(adjustl(file_index))
      call read_h5_one_r3d(file_id, dims, temp_real, file_name)
      do k=1,ctot
      do j=1,ctot
      do i=1,ctot
        pdr(i,j,k)%abun(isp) = temp_real(i,j,k)
      end do
      end do
      end do
    end do

    call h5fclose_f(file_id, h5err)

#ifdef CRATTENUATION
    in_file_h5 = trim(adjustl(filein)) // ".cr.h5"
    inquire(file=trim(in_file_h5), exist=file_exists)
    if (.not. file_exists) then
      print *, "Error: Input file not found: ", trim(in_file_h5)
      stop
    end if

    if (file_exists) then
      call h5fopen_f(in_file_h5, H5F_ACC_RDONLY_F, file_id, h5err)
      if (h5err /= 0) then
        print *, "Failed to open file: ", trim(in_file_h5)
        stop
      end if

      write(6,*) 'Reading CRattenuation'
      call read_h5_one_r3d(file_id, dims, temp_real, "zetalocal"); pdr%zetalocal = temp_real

      call h5fclose_f(file_id, h5err)
    end if
#endif

step = abs(pdr(1,1,2)%z-pdr(1,1,1)%z)*PC
#ifdef VELOCITY
    in_file_h5 = filevel
    inquire(file=trim(in_file_h5), exist=file_exists)
    if (.not. file_exists) then
      print *, "Error: Input file not found: ", trim(in_file_h5)
      stop
    end if

    if (file_exists) then
      call h5fopen_f(in_file_h5, H5F_ACC_RDONLY_F, file_id, h5err)
      if (h5err /= 0) then
        print *, "Failed to open file: ", trim(in_file_h5)
        stop
      end if

      write(6,*) 'Reading velocity'
      call read_h5_one_r3d(file_id, dims, temp_real, "vx"); pdr%vx = temp_real
      call read_h5_one_r3d(file_id, dims, temp_real, "vy"); pdr%vy = temp_real
      call read_h5_one_r3d(file_id, dims, temp_real, "vz"); pdr%vz = temp_real

      call h5fclose_f(file_id, h5err)
    endif 
#endif

  deallocate(temp_real)
  end subroutine readpdrfiles_h5

  subroutine readspopfile_h5
    implicit none

    ! locals
    integer :: i, j, k, ico, total_levels, kr
    integer(HID_T) :: dset_id, filespace
    integer(HSIZE_T), dimension(4) :: dims3
    double precision, allocatable :: allpop(:,:,:,:)
    logical :: file_exists

    in_file_h5 = trim(adjustl(filein)) // trim(adjustl(".RTspop")) // ".h5"
    inquire(file=trim(in_file_h5), exist=file_exists)
    if (.not. file_exists) then
      print *, "Warning: RTspop file not found: ", trim(in_file_h5)
      stop
    end if    

    call h5fopen_f(in_file_h5, H5F_ACC_RDONLY_F, file_id, h5err)
    if (h5err /= 0) then
      print *, "Error: Failed to open RTspop file: ", trim(in_file_h5)
      stop
    end if

    total_levels = 0
    do ico = 1, coo
      total_levels = total_levels + coolant(ico)%cnlev
    end do

    allocate(allpop(1:ctot, 1:ctot, 1:ctot, 1:total_levels))
    dims3 = (/ctot, ctot, ctot, total_levels/)

    call h5dopen_f(file_id, "allpop", dset_id, h5err)
    if (h5err /= 0) then
      print *, "Error: Dataset 'allpop' not found."
      deallocate(allpop)
      call h5fclose_f(file_id, h5err)
      stop
    end if

    write(6,*) 'Reading level populations'
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, allpop, dims3, h5err)
    call h5dclose_f(dset_id, h5err)
    if (h5err /= 0) then
      print *, "Error: Failed to read 'allpop' dataset."
      deallocate(allpop)
      call h5fclose_f(file_id, h5err)
      stop
    end if

    kr=0
    if(jr.ge.2) kr = sum(coolant(1:jr-1)%cnlev)
    do k = 1, ctot
    do j = 1, ctot
    do i = 1, ctot
      pdr(i,j,k)%tpop(1) = allpop(i,j,k,kr+lr)
      pdr(i,j,k)%tpop(2) = allpop(i,j,k,kr+lr+1)
    end do
    end do
    end do

    deallocate(allpop)
    call h5fclose_f(file_id, h5err)    

  end subroutine readspopfile_h5

  subroutine readspopfile_h5_lr
    implicit none

    integer :: i, j, k, ico, total_levels, kr
    integer(HID_T) :: dset_id, filespace, memspace
    integer(HSIZE_T), dimension(4) :: offset, count, memdims
    double precision, allocatable :: pop_slice(:,:,:,:)
    logical :: file_exists

    in_file_h5 = trim(adjustl(filein)) // trim(adjustl(".RTspop")) // ".h5"
    inquire(file=trim(in_file_h5), exist=file_exists)
    if (.not. file_exists) then
        print *, "Warning: RTspop file not found: ", trim(in_file_h5)
        stop
    end if    
    call h5fopen_f(in_file_h5, H5F_ACC_RDONLY_F, file_id, h5err)

    total_levels = 0
    do ico = 1, coo
        total_levels = total_levels + coolant(ico)%cnlev
    end do

    kr = 0
    if(jr >= 2) kr = sum(coolant(1:jr-1)%cnlev)

    call h5dopen_f(file_id, "allpop", dset_id, h5err)
    call h5dget_space_f(dset_id, filespace, h5err)

    offset = [0, 0, 0, kr + lr - 1]
    count = [ctot, ctot, ctot, 2]
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, count, h5err)

    memdims = [ctot, ctot, ctot, 2]
    call h5screate_simple_f(4, memdims, memspace, h5err)

    allocate(pop_slice(ctot, ctot, ctot, 2))
    write(6,*) 'Reading only two levels: ', kr+lr, ' and ', kr+lr+1
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, pop_slice, memdims, h5err, &
                   file_space_id=filespace, mem_space_id=memspace)

    do k = 1, ctot
    do j = 1, ctot
    do i = 1, ctot
      pdr(i,j,k)%tpop(1) = pop_slice(i,j,k,1)  ! kr+lr
      pdr(i,j,k)%tpop(2) = pop_slice(i,j,k,2)  ! kr+lr+1
    end do
    end do
    end do

    deallocate(pop_slice)
    call h5sclose_f(memspace, h5err)
    call h5sclose_f(filespace, h5err)
    call h5dclose_f(dset_id, h5err)
    call h5fclose_f(file_id, h5err)
  end subroutine readspopfile_h5_lr

  subroutine readspopheader_h5
    integer(HID_T) :: dset_id_network, dtype_id, dspace_id
    integer(HID_T) :: dset_id_cooarray, arrtype_id
    integer(HSIZE_T), dimension(1) :: dims
    character(len=128) :: network_name
    character(len=128), allocatable :: coo_array(:)
    logical :: file_exists

    in_file_h5 = trim(adjustl(filein)) // trim(adjustl(".RTspop")) // ".h5"
    inquire(file=trim(in_file_h5), exist=file_exists)
    if (.not. file_exists) then
      print *, "Warning: RTspop file not found: ", trim(in_file_h5)
      stop
    end if    

    call h5fopen_f(in_file_h5, H5F_ACC_RDONLY_F, file_id, h5err)
    if (h5err /= 0) then
      print *, "Error: Failed to open RTspop file: ", trim(in_file_h5)
      stop
    end if

    call h5dopen_f(file_id, "network", dset_id_network, h5err)
    call h5dget_type_f(dset_id_network, dtype_id, h5err)
    call h5dread_f(dset_id_network, dtype_id, network_name, dims, h5err)

    call h5dclose_f(dset_id_network, h5err)
    call h5tclose_f(dtype_id, h5err)

    network = trim(network_name)
    ! print *, "Network name: ", network

    ! call h5dopen_f(file_id, "coolfiles", dset_id_cooarray, h5err)
    ! call h5dget_space_f(dset_id_cooarray, dspace_id, h5err)
    ! call h5sget_simple_extent_dims_f(dspace_id, dims, dims, h5err)
    ! coo = dims(1)

    ! allocate(coo_array(coo))
    ! call h5dget_type_f(dset_id_cooarray, arrtype_id, h5err)
    ! call h5dread_f(dset_id_cooarray, arrtype_id, coo_array, dims, h5err)

    ! call h5dclose_f(dset_id_cooarray, h5err)
    ! call h5tclose_f(arrtype_id, h5err)

    ! do i = 1, coo
    !   coolfile(i) = trim(coo_array(i))
    !   ! print*, "Coolfile",i,coolfile(i),coo
    ! enddo

    coo = 5
    coolfile(1) = 'chemfiles/12co.dat'
    coolfile(2) = 'chemfiles/12c+.dat'
    coolfile(3) = 'chemfiles/12c.dat'
    coolfile(4) = 'chemfiles/16o.dat'
    coolfile(5) = 'chemfiles/hcop.dat'
    
    ! deallocate(coo_array)
  end subroutine readspopheader_h5

  subroutine read_h5_one_r3d(file_id, dims, temp_real, name)
    implicit none
    integer(HID_T), intent(in) :: file_id
    integer(HSIZE_T), dimension(3), intent(in) :: dims
    double precision, intent(out) :: temp_real(:,:,:)
    character(len=*), intent(in) :: name

    integer(HID_T) :: filespace, dset_id
    integer :: h5err

    call h5dopen_f(file_id, trim(name), dset_id, h5err)
    if (h5err /= 0) then
      print *, "Error: Dataset not found: ", trim(name)
      stop 1
    end if

    call h5dread_f(dset_id, H5T_IEEE_F64LE, temp_real, dims, h5err)
    if (h5err /= 0) then
      print *, "Error: Failed to read dataset: ", trim(name)
      stop 1
    end if

    call h5dclose_f(dset_id, h5err)
  end subroutine read_h5_one_r3d

end module m_readh5

