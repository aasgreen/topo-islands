program islandC_np
    use hdf5 !use the hdf5 dataset
implicit none
!This program will create a linear region in the center where a non-polar force will be applied. 
!Then, we will watch as the dynamics evolve.
character(100) :: buffer
integer :: n = 200, endt=1001,seedn
real(8) :: beta,mu,measuredt,zeroe, meank,avedelphi
real(8), allocatable, dimension(:) :: muvector
integer, dimension(100) :: tpoints
integer, allocatable, dimension(:) :: logtpoints
real(8), allocatable, dimension(:,:) :: grid,dgrid,gridplusdelta,agrid, hgrid
real(8) :: x,kappa,windingn, lnoise
real(8),parameter :: pi = 4*atan(1.0_8),  deltime=0.01
integer, allocatable :: seed(:)
character(40) :: filenames
character(40) :: filenames_h
character(40) :: dfilenames


real(8) :: t1=1.2,t2=1.0,dt,dvar=0, islandr, deltar
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

!calculate langenvin noise term (c_l in yurke)
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
do i=1,n
    do j=1,n
!        write(*,*) i,j
        call random_number(x)
        grid(i,j) = 0
    end do
enddo

!gridplusdelta is the grid that evolves in time
gridplusdelta = grid
hgrid = 0
measuredt = 0
do i=1,n
    do j=1,n
        measuredt = (hamxy(i,j,grid(i,j),kappa,mu, 0d0))/n/n+measuredt
    enddo
enddo

!the following code deals with the radius of an island.
islandr = 0 ! the beginning island radius
deltar = 1 !the beginning width of the island
!deltar = 10 !deltar is the width of the permimeter of the circle

!allocate mu vector (this is for the 
allocate(muvector(int(endt/deltime)))
muvector = 5.*kappa
!muvector(1:10)=0
!muvector(int(endt/20/deltime):int(endt/deltime)) = (/ (mu+kappa*i/(endt/deltime-2)*100/2, i= 0,int((endt/deltime-1)/2)) /)
!muvector(int(endt/20/deltime):int(endt/deltime)) = 10*kappa
!write(*,*) muvector



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
write(*,*) endt
!this next section will write out the boundary condition for the island. This is mostly just for piece of mind, so folks
!can visualize what the angular dependence of the anchoring condition is.
call h5fcreate_f('anchor.h5', h5f_acc_trunc_f, file_id_a, error)
do i=1,n
    do j = 1, n
        agrid(i,j) = anchor(i,j,20d0,3,n)
        write(68,'(f10.5)',advance='no') anchor(i,j,20d0,3,n)
    enddo

            write(68,*)
end do

close(68)

call h5screate_simple_f(rank, dims, dspace_id_a, error)
call h5dcreate_f(file_id_a, 'anchor.dat', h5t_native_double, dspace_id_a, &
       dset_id_a, error)

call h5dwrite_f(dset_id_a, h5t_native_double, agrid, dims, error)
call h5dclose_f(dset_id_a, error)
call h5fclose_f(file_id_a, error)


do t=1,int(endt)

    avedelphi = 0
    grid=gridplusdelta
    if (( t .gt. 10 )) then
        call islandradius(islandr,t)
    endif
    call update(grid,gridplusdelta, hgrid,lnoise,n)
 
    meank =.5*sum((grid-gridplusdelta)**2)/n/n/deltime**2
    write(filenames,'(a3,f0.4, a4)') 'out', t*deltime,'.dat'
    write(dfilenames,'(a6,f0.4, a4)') 'defect', t*deltime,'.dat'

    if (logspace(t) .eq. 1) then
        print *, trim(filenames)
    end if
    if (1 .eq. 1) then !this is not log time, so we need all time steps
       write(66,*) islandr, 1.0*real(t)/endt
       write(67,*) filenames
       !print *, trim(filenames)

!        open(1,file=filenames)
!        open(3,file=dfilenames)

        ! calculate defects and average energy, and write grid to file
        measuredt = 0.

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


!        do i=1,n
!            do j=1,n
!            write(1,'(f10.5)',advance="no") acos(cos(grid(i,j)))
!            measuredt = (hamxy(i,j,grid(i,j),kappa,mu))/n/n+measuredt
!            windingn=(windn(grid(i,modulo(j-2,n)+1)-grid(i,j)))+&
!                &(windn(grid(modulo(i,n)+1,modulo(j-2,n)+1)-grid(i,modulo(j-2,n)+1)))+&
!                &(windn(grid(modulo(i,n)+1,j)-grid(modulo(i,n)+1,modulo(j-2,n)+1)))+&
!                &(windn(grid(i,j)-grid(modulo(i,n)+1,j)))
!            if (windingn .ge. 1) then
!                dgrid(i,j) =1
!            else if (-1*windingn .ge. 1) then
!                dgrid(i,j)=-1
!            else
!                dgrid(i,j)=0
!            endif
!           write(3,'(f10.5)',advance="no") dgrid(i,j)
!            enddo
!            write(1,*)
!            write(3,*)
!        enddo
!        measuredt = (measuredt-zeroe)
!        write(61,'(f12.2,a, f10.5)') t*deltime,',', measuredt, meank
        !write(*,*) t*deltime

!        close(1)
!        close(3)
    endif
!write(*,*) 'closing'
end do
!write(*,*) dvar/dcount/2
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

    subroutine islandradius(islandr,t)
        !calculate the radius of the island at a given time
        integer :: t
        real(8) :: islandr
        ! case 2: create island out of nothing at .2 percent
        if ( (t .ge. 100) .and. (t .le. 400 ) ) then
            islandr = 20 !this will put the line in the middle of the simulation
        else
            islandr = 0
        end if
    end subroutine islandradius

    function logspace(t)
        integer :: t, logspace
        if (t .eq. logtpoints(timeprint)) then
            logspace = 1
            timeprint = timeprint +1
        else
            logspace = 0
        endif
    end function logspace

    function torque(i,j,theta,kappa,mu,fieldphi)
        integer :: i,j,ii,jj
        real(8):: theta,torque,kappa,fieldphi
        real(8) :: g,hamxy,mu
        real(8), dimension(3) :: x,y
        !write(*,*) 'hamxy', i,j
        x =(/ grid(modulo(i-2,n)+1,j),theta,grid(modulo(i,n)+1,j)/)
        y =(/ grid(i,modulo(j-2,n)+1),theta,grid(i,modulo(j,n)+1)/)
        torque = kappa*(sin(x(2)-x(1))+sin(x(2)-x(3))+sin(y(2)-y(1))+sin(y(2)-y(3)))+4*mu*cos( (theta-fieldphi))* &
            sin((theta-fieldphi) ) 
    end function torque

    function xi(grid,beta,n)
        integer :: n
        real(8), dimension(n,n) :: grid
        real(8) :: beta, xi
        xi = sum(grid**2)/n**2
    end function xi

    function hamxy(i,j,theta,kappa,mu,fieldphi)
        integer :: i,j,ii,jj
        real(8):: theta, fieldphi
        real(8) :: kappa,hamxy,mu
        real(8), dimension(3) :: x,y
        !write(*,*) 'hamxy', i,j
        x =(/ grid(modulo(i-2,n)+1,j),theta,grid(modulo(i,n)+1,j)/)
        y =(/ grid(i,modulo(j-2,n)+1),theta,grid(i,modulo(j,n)+1)/)
        hamxy = -kappa*(cos(x(2)-x(1))+cos(x(2)-x(3))+cos(y(2)-y(1))+cos(y(2)-y(3)))-mu*cos( (theta-fieldphi) )
    end function hamxy
    
   function anchor(ii,jj,ir,dr,n)
      integer :: n, ii, jj, dr
     real(8) :: anchor,ir
     write(*,*) ii,jj, ir, dr, n
    if ( (sqrt(real((ii-n/2)**2 + (jj-n/2)**2)) .lt. ir) .and. &
                    (sqrt(real((ii-n/2)**2+(jj-n/2)**2)) .gt. (ir-dr))) then

                    anchor = atan2(real(jj-n/2),real(ii-n/2))+pi/2
                else
                    anchor = 0.d8
                endif

        end function anchor


    subroutine update(state,newstate,hgrid,lnoise,n)
        real(8) :: theta, thetap,x,dele,e1,e2,lnoise,fluce,u,force,fieldphi,mutemp
        integer :: n,ii,jj,i,j
        real(8), dimension(n,n) :: state,newstate,hgrid
        do i=1,n
            do j=1,n
                call random_number(x)
                ii = 1+floor((n)*x) !this will give periodic boundaries
                call random_number(x)
                jj = 1+floor((n)*x)
                !see if we are within the line boundary
                 if ( (sqrt(real((ii-n/2)**2 + (jj-n/2)**2)) .lt. (islandr+real(deltar)/2.) ) .and. &
                    (sqrt(real((ii-n/2)**2+(jj-n/2)**2)) .gt. (islandr-real(deltar)/2.))) then

                    mutemp = muvector(t) !using global variable laser_noise
                    fieldphi = atan2(real(jj-n/2),real(ii-n/2))+pi/2
                else
                    mutemp = mu
                    fieldphi = 0
                endif

                theta = grid(ii,jj)
                force = torque(ii,jj,theta,kappa,mutemp,fieldphi)
                call random_number(fluce)
                fluce = fluce-.5
                newstate(ii,jj)=theta-deltime*(fluce*lnoise+force)
                hgrid(ii,jj) = hamxy(ii, jj, newstate(ii,jj), kappa, mu, fieldphi)+4*kappa !this will normalize so zero is actually the
        !energy
                avedelphi = avedelphi + deltime*fluce*lnoise
                dvar = dvar+(fluce*lnoise)**2
                dcount = dcount+1
            end do
        end do
    end subroutine update
            

    function dtrack(dgrid,grid)
        real(8), dimension(n,n) :: dgrid
        real(8), dimension(3,3) :: grid
        real(8)                 :: angle=0
        integer                 :: dtrack,i,j
        do i=1,3
            !angle=angle+grid(i,j)-grid(2,2)/9., j=1,3)
            !write(*,*) (grid(i,j), j=1,3)
        enddo
        if (angle .ge. 2*pi) then
            dtrack=1

        else  
            dtrack=0
        endif

    endfunction dtrack

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
end program islandC_np
