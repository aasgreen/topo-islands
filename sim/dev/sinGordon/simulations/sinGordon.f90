program sinGordon
    use hdf5 !use the hdf5 dataset
    use langan
implicit none
!This program will create a linear region in the center where a non-polar force will be applied. 
!Then, we will watch as the dynamics evolve.
character(100) :: buffer
integer :: n = 200,seedn
real(8) :: beta,mu,measuredt,zeroe, meank,avedelphi
integer, dimension(100) :: tpoints
integer, allocatable, dimension(:) :: logtpoints
real(8), allocatable, dimension(:,:) :: grid,dgrid,gridplusdelta,agrid, hgrid
real(8) :: x,kappa,windingn, lnoise
integer, allocatable :: seed(:)
character(40) :: filenames
character(40) :: filenames_h
character(40) :: dfilenames


integer :: i,j,t,timeprint=1,tt,dcount=0,un=4,seedsize,istat=0, deltat

!hdf parameters

character(len=8), parameter :: filename = "dsetf.h5" ! file name
character(len=8), parameter :: hamfilename = 'hsetf.h5'
character(len=10), parameter :: anchorName = 'anchor.dat'

integer(hid_t) :: file_id       ! file identifier
integer(hid_t) :: dset_id       ! dataset identifier
integer(hid_t) :: dspace_id     ! dataspace identifier

integer(hid_t) :: file_id_h   ! file identifier
integer(hid_t) :: dset_id_h      ! dataset identifier
integer(hid_t) :: dspace_id_h    ! dataspace identifier

integer(hid_t) :: file_id_a   ! file identifier
integer(hid_t) :: dset_id_a      ! dataset identifier
integer(hid_t) :: dspace_id_a    ! dataspace identifier



integer(hsize_t), dimension(2) :: dims ! dataset dimensions

integer     ::   rank = 2                        ! dataset rank

integer     ::   error ! error flag



!read in arguments (g,beta,mu,n,endt)
call getarg(1,buffer)
read(buffer,*) kappa !nn-coupling

call getarg(2,buffer)
read(buffer,*) beta !1/T

write(*,*) 1.0/beta
call getarg(3,buffer)
read(buffer,*) mu !strength of electric field

call getarg(4,buffer)
read(buffer,*) n !size of grid

call getarg(5,buffer)
read(buffer,*) endt !number of iterations

call getarg(6,buffer)
read(buffer,*) seedn !random seed


call random_seed(seedsize)
allocate(seed(seedsize))



  ! first try if the os provides a random number generator
open(newunit=un, file="/dev/urandom", access="stream", &
       form="unformatted", action="read", status="old", iostat=istat)
if (istat == 0) then
    read(un) seed
    close(un)
endif

!calculate zero temperature energy
zeroe = -1*kappa*4.0

!lnoise=sqrt(24*1./beta*deltime)
lnoise = sqrt(12.0*1.0/beta)
!this calculates evenlly space logarythmic points
!it is useful to have if you don't want all the points from the simulation,
!as the events progress on a log time scale. This can save time if you are looking at defect coalenscenece.
do t=1,100
    dt = (log10(real(endt/deltime))/100)
    tpoints(t) =  int(10**(t*dt))
enddo
logtpoints= tpoints(unique(tpoints))


!!initial random/zero grid

allocate(grid(n,n))
allocate(hgrid(n,n))
allocate(gridplusdelta(n,n))
allocate(agrid(n,n))

!this part will create a random grid. However, because we are interested in seeing equilibrium,
!we will instead just normalize this to zero.
call random_seed(put=seed)

!initialize grid to have a soliton on it
call initializeGrid(grid,n)


!gridplusdelta is the grid that evolves in time
gridplusdelta = grid
hgrid = 0
measuredt = 0
do i=1,n
    do j=1,n
        measuredt = (hamxy(i,j,grid(i,j), n,kappa,mu, 0d0))/n/n+measuredt
    enddo
enddo

!the following code deals with the radius of an island.
islandr = 0 ! the beginning island radius
deltar = 1 !the beginning width of the island



write(*,*) 'max t', measuredt-zeroe

!open temperature v time file
open(61,file='tvt.dat',status = 'unknown', position='append')
open(66,file='radius.dat',status = 'unknown', position='append')
open(67, file = 'datanames.dat', status = 'unknown', position = 'append')
open(68, file = anchorName, status = 'unknown', position = 'append')

!hdf handling 

dims = (/n,n/) 
  !
  ! initialize fortran interface.
  !
  call h5open_f(error)

  !
  ! create a new file using default properties.
  !
  call h5fcreate_f(filename, h5f_acc_trunc_f, file_id, error)
  call h5fcreate_f(hamfilename, h5f_acc_trunc_f, file_id_h, error)



!initialize defect grid
allocate(dgrid(n,n))
dgrid=0
endt = logtpoints(size(logtpoints))

call h5screate_simple_f(rank, dims, dspace_id_a, error)
call h5dcreate_f(file_id_a, 'anchor.dat', h5t_native_double, dspace_id_a, &
       dset_id_a, error)

call h5dwrite_f(dset_id_a, h5t_native_double, agrid, dims, error)
call h5dclose_f(dset_id_a, error)
call h5fclose_f(file_id_a, error)


do t=1,int(endt)

    avedelphi = 0
    grid=gridplusdelta

    call update(grid,gridplusdelta, hgrid, lnoise, kappa, mu, n)
 
    meank =.5*sum((grid-gridplusdelta)**2)/n/n/deltime**2
    write(filenames,'(a3,f0.4, a4)') 'out', t*deltime,'.dat'
    write(dfilenames,'(a6,f0.4, a4)') 'defect', t*deltime,'.dat'

    if (logspace(t) .eq. 1) then
        print *, trim(filenames)
    end if
    if (1 .eq. 1) then !this is not log time, so we need all time steps
       write(66,*) islandr, 1.0*real(t)/endt
       write(67,*) filenames


  !
  ! create the dataspace.
  !
  call h5screate_simple_f(rank, dims, dspace_id, error)
  call h5screate_simple_f(rank, dims, dspace_id_h, error)

  !
  ! create the dataset with default properties.
  !
  call h5dcreate_f(file_id, filenames, h5t_native_double, dspace_id, &
       dset_id, error)


  call h5dcreate_f(file_id_h, filenames, h5t_native_double, dspace_id_h, &
       dset_id_h, error)
  ! write the data
  call h5dwrite_f(dset_id, h5t_native_double, grid, dims, error)
  call h5dwrite_f(dset_id_h, h5t_native_double, hgrid, dims, error)
  !
  ! end access to the dataset and release resources used by it.
  !
  call h5dclose_f(dset_id, error)
  call h5dclose_f(dset_id_h, error)

  !
  ! terminate access to the data space.
  !
  call h5sclose_f(dspace_id, error)
  call h5sclose_f(dspace_id_h, error)

    endif

end do
close(61)
close(66)

!write to file

  !
  ! close the file.
  !
  call h5fclose_f(file_id, error)
  call h5fclose_f(file_id_h, error)

  !
  ! close fortran interface.
  !
  call h5close_f(error)


write(*,*) 'finished langin'

deallocate(grid)

contains

    function logspace(t)
        integer :: t, logspace
        if (t .eq. logtpoints(timeprint)) then
            logspace = 1
            timeprint = timeprint +1
        else
            logspace = 0
        endif
    end function logspace


end program sinGordon
