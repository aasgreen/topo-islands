program defectLoc
    use hdf5 !use the hdf5 dataset
!program that will:
!1. read in hd5 grid representing spin angles of xy model
!2. find and label all defects
!3. output this as another hdf5
!this is a seperate program (though I may make it a module). It would
!be more efficient to calculate this during run time to cut down on IO
!but, I think it is more valuable to have this seperate, so we can run it
! over and over again on the same data set.
implicit none

character(100) :: buffer
real(8), allocatable, dimension(:,:) :: grid, dgrid
real(8),parameter :: pi = 4*atan(1.0_8),  deltime=0.1

real(8) :: windingn, avedelphi, dcount, deltar, dvar


integer :: ii, i,j,t, N, NFiles

!hdf parameters

!character(len=8), parameter :: filename = "dsetf.h5" ! file name
character(100) :: datafileName
character(100) :: datasetFileName
character(100) :: dset_name

integer(hid_t) :: file_id, file_idout  ! file identifier
integer(hid_t) :: dset_id, dset_idout    ! dataset identifier
integer(hid_t) :: dspace_id, dspace_idout! dataspace identifier

integer(hsize_t), dimension(2) :: dims ! dataset dimensions

integer     ::   rank = 2                        ! dataset rank

integer     ::   error ! error flag
type(h5o_info_t), target :: infobuf
INTEGER(HSIZE_T), DIMENSION(2) :: data_dims
character(256), allocatable, dimension(:) :: dsetNames





call getarg(1,buffer)
read(buffer,"(A100)") datafileName

call getarg(2, buffer)
read (buffer, "(A100)") datasetFileName

call getarg(3,buffer)
read(buffer,*) N
data_dims(1) = N
data_dims(2) = N

call getarg(4, buffer)
read(buffer, *) NFiles

allocate(dsetNames(NFiles))
open(unit = 99, file = datasetFileName, status = 'old', action='read')
do i = 1, NFiles
    read(99,*) dsetNames(i)
    write(*,*) dsetNames(i)
enddo
!initialize file to output defect grids using hd5 format
!  ! initialize fortran interface.
!  !
call h5open_f(error)
!
!  !
!  ! create a new file using default properties.
!  !
  call h5fcreate_f('defectGrid.h5', h5f_acc_trunc_f, file_idout, error)
!

allocate(grid(N,N))
allocate(dgrid(N,N))
grid = -1
dgrid = 0
call h5open_f(error)
write(*,*) 'reading in'
!CALL h5dopen_f(datafileName, dsetname, dset_id, error)
write(*,*) dsetNames(10)
!open the hdf5 file
CALL h5fopen_f(datafileName, H5F_ACC_RDWR_F, file_id, error)
write(*,*) 'file open'
do ii = 1, NFiles
    write(*,*) trim(dsetNames(ii))

!iteratre through all the datasets
CALL h5dopen_f(file_id, dsetNames(ii), dset_id, error)

CALL h5dread_f(dset_id, H5T_NATIVE_double, grid, data_dims, error)
CALL h5dclose_f(dset_id, error)

!!hdf handling 
!
!dims = (/n,n/) 
!  !

!
        do i=1,n
            do j=1,n
!!            write(1,'(f10.5)',advance="no") acos(cos(grid(i,j)))
!!            measuredt = (hamxy(i,j,grid(i,j),kappa,mu))/n/n+measuredt
            windingn=(windn(grid(i,modulo(j-2,n)+1)-grid(i,j)))+&
                &(windn(grid(modulo(i,n)+1,modulo(j-2,n)+1)-grid(i,modulo(j-2,n)+1)))+&
                &(windn(grid(modulo(i,n)+1,j)-grid(modulo(i,n)+1,modulo(j-2,n)+1)))+&
                &(windn(grid(i,j)-grid(modulo(i,n)+1,j)))
            if (windingn .ge. 1) then
                dgrid(i,j) =1
            else if (-1*windingn .ge. 1) then
                dgrid(i,j)=-1
            else
                dgrid(i,j)=0
            endif
!!           write(3,'(f10.5)',advance="no") dgrid(i,j)
            enddo
!!            write(1,*)
!!            write(3,*)
        enddo
  !
  ! create the dataspace.
  !
  call h5screate_simple_f(rank, data_dims, dspace_idout, error)

  !
  ! create the dataset with default properties.
  !
  call h5dcreate_f(file_idout, dsetNames(ii), h5t_native_double, dspace_idout, &
       dset_idout, error)

  ! write the data
  call h5dwrite_f(dset_idout, h5t_native_integer, dgrid, dims, error)
  !
  ! end access to the dataset and release resources used by it.
  !
  call h5dclose_f(dset_idout, error)

  !
  ! terminate access to the data space.
  !
  call h5sclose_f(dspace_idout, error)
  enddo

call h5fclose_f(file_idout, error)
CALL h5fclose_f(file_id, error)

  !
  ! close fortran interface.
  !
call h5close_f(error)


write(*,*) 'finished langin'
deallocate(dgrid)
deallocate(grid)
contains


    function xi(grid,beta,n)
        integer :: n
        real(8), dimension(n,n) :: grid
        real(8) :: beta, xi
        xi = sum(grid**2)/n**2
    end function xi

    function hamxy(i,j,theta,kappa,mu)
        integer :: i,j,ii,jj
        real(8):: theta
        real(8) :: kappa,hamxy,mu
        real(8), dimension(3) :: x,y
        !write(*,*) 'hamxy', i,j
        x =(/ grid(modulo(i-2,n)+1,j),theta,grid(modulo(i,n)+1,j)/)
        y =(/ grid(i,modulo(j-2,n)+1),theta,grid(i,modulo(j,n)+1)/)
        hamxy = -kappa*(cos(x(2)-x(1))+cos(x(2)-x(3))+cos(y(2)-y(1))+cos(y(2)-y(3)))-mu*cos( (theta-45/2/pi) )
    end function hamxy



    function angledist(theta1,theta2)
        real(8) :: d1,d2,angledist,theta1,theta2
        angledist = acos(cos(theta1-theta2))
        end function angledist

    function windn(angle)
        real(8) :: angle, windn
        windn = (angle-asin(sin(angle))) / (pi) 
        end function windn



function remove_dups(input)
  integer :: input(100)       ! the input
  integer :: remove_dups(size(input))  ! the output
  integer :: k                   ! the number of unique elements
  integer :: i, j
 
  k = 1
  remove_dups(1) = input(1)
  outer: do i=2,size(input)
     do j=1,k
        if (remove_dups(j) == input(i)) then
           ! found a match so start looking again
           cycle outer
        end if
     end do
     ! no match found so add it to the output
     k = k + 1
     remove_dups(k) = input(i)
  end do outer
  end function remove_dups
    function unique(input)
     !   find "indices", the list of unique numbers in "list"
     integer( kind = 4 ) :: kx, input(100)
     integer( kind = 4 ),allocatable :: unique(:)
     logical :: mask(100)
     mask(1)=.true.
     do kx=100,2,-1
       mask(kx)= .not.(any(input(:kx-1)==input(kx)))
     end do
     unique=pack([(kx,kx=1,100)],mask)
    end function unique
end program defectLoc
