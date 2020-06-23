program langan
    use hdf5 !use the hdf5 dataset
implicit none

character(100) :: buffer
integer :: n = 200, endt=1001,seedn
real(8) :: beta,mu,measuredt,zeroe, meank,avedelphi
real(8), allocatable, dimension(:) :: muvector
integer, dimension(100) :: tpoints
integer, allocatable, dimension(:) :: logtpoints
real(8), allocatable, dimension(:,:) :: grid,dgrid,gridplusdelta
real(8) :: x,kappa,windingn,lnoise
real(8),parameter :: pi = 4*atan(1.0_8),  deltime=0.1
integer, allocatable :: seed(:)
character(40) :: filenames
character(40) :: dfilenames


real(8) :: t1=1.2,t2=1.0,dt,dvar=0,islandr
integer :: i,j,t,timeprint=1,tt,dcount=0,un=4,seedsize,istat=0,deltar

!hdf parameters

character(len=8), parameter :: filename = "dsetf.h5" ! file name

integer(hid_t) :: file_id       ! file identifier
integer(hid_t) :: dset_id       ! dataset identifier
integer(hid_t) :: dspace_id     ! dataspace identifier

integer(hsize_t), dimension(2) :: dims ! dataset dimensions

integer     ::   rank = 2                        ! dataset rank

integer     ::   error ! error flag



!read in arguments (g,beta,mu,n,endt)
call getarg(1,buffer)
read(buffer,*) kappa

call getarg(2,buffer)
read(buffer,*) beta

write(*,*) 1.0/beta
call getarg(3,buffer)
read(buffer,*) mu

call getarg(4,buffer)
read(buffer,*) n

call getarg(5,buffer)
read(buffer,*) endt

call getarg(6,buffer)
read(buffer,*) seedn

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
lnoise = sqrt(12*1/beta)
do t=1,100
    dt = (log10(real(endt/deltime))/100)
    tpoints(t) =  int(10**(t*dt))
enddo
logtpoints= tpoints(unique(tpoints))
!write(*,*) logtpoints


!initial random grid
!write(*,*) 'initialize grid'
allocate(grid(n,n))
allocate(gridplusdelta(n,n))
!write(*,*) grid(n,n)




call random_seed(put=seed)
do i=1,n
    do j=1,n
!        write(*,*) i,j
        call random_number(x)
        grid(i,j) = 0
    end do
enddo
gridplusdelta = grid
measuredt = 0
do i=1,n
    do j=1,n
        measuredt = (hamxy(i,j,grid(i,j),kappa,mu))/n/n+measuredt
    enddo
enddo

!calculate the indices where the boundary exists. these indicies will meet the
!equation i^2+j^2 = r^2, where r is the radius of the circle. we will take some 
!delta, so (r-del)^2 .ge. i^2+j^2 .ge. (r)^2
islandr =10
deltar = 3

!allocate mu vector
allocate(muvector(int(endt/deltime)))
muvector = 11.*kappa
!muvector(1:10)=0
!muvector(int(endt/20/deltime):int(endt/deltime)) = (/ (mu+kappa*i/(endt/deltime-2)*100/2, i= 0,int((endt/deltime-1)/2)) /)
!muvector(int(endt/20/deltime):int(endt/deltime)) = 10*kappa
write(*,*) muvector



write(*,*) 'max t', measuredt-zeroe

!open temperature v time file
open(61,file='tvt.dat',status = 'unknown', position='append')
open(66,file='radius.dat',status = 'unknown', position='append')
open(67, file = 'datanames.dat', status = 'unknown', position = 'append')

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

  !
  ! create the dataspace.
  !
!  call h5screate_simple_f(rank, dims, dspace_id, error)

  !
  ! create the dataset with default properties.
  !
!  call h5dcreate_f(file_id, dsetname, h5t_native_integer, dspace_id, &
!       dset_id, error)

  !
  ! end access to the dataset and release resources used by it.
  !
!  call h5dclose_f(dset_id, error)

  !
  ! terminate access to the data space.
  !
!  call h5sclose_f(dspace_id, error)

  !
  ! close the file.
  !
!  call h5fclose_f(file_id, error)

  !
  ! close fortran interface.
  !
!  call h5close_f(error)


!initialize defect grid
allocate(dgrid(n,n))
dgrid=0
endt = logtpoints(size(logtpoints))
write(*,*) endt
!write(*,*) 'begining langan processing...'
!preform langan with euler update
do t=1,int(endt)

    avedelphi = 0
    grid=gridplusdelta
    if (( t .gt. 50 )) then
        call islandradius(islandr,t)
    endif
    call update(grid,gridplusdelta, lnoise,n)
    meank =.5*sum((grid-gridplusdelta)**2)/n/n/deltime**2
    write(filenames,'(a3,f0.4, a4)') 'out', t*deltime,'.dat'
    write(dfilenames,'(a6,f0.4, a4)') 'defect', t*deltime,'.dat'
    write(67,*) filenames

    if (logspace(t) .eq. 1) then
        print *, trim(filenames)
    end if
    if (1 .eq. 1) then !this is not log time, so we need all time steps
        write(66,*) islandr
       !print *, trim(filenames)

!        open(1,file=filenames)
!        open(3,file=dfilenames)

        ! calculate defects and average energy, and write grid to file
        measuredt = 0.

  !
  ! create the dataspace.
  !
  call h5screate_simple_f(rank, dims, dspace_id, error)

  !
  ! create the dataset with default properties.
  !
  call h5dcreate_f(file_id, filenames, h5t_native_double, dspace_id, &
       dset_id, error)

  ! write the data
  call h5dwrite_f(dset_id, h5t_native_double, grid, dims, error)
  !
  ! end access to the dataset and release resources used by it.
  !
  call h5dclose_f(dset_id, error)

  !
  ! terminate access to the data space.
  !
  call h5sclose_f(dspace_id, error)


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
        ! case 2: grow, hit 50 and then shrink
        if (1.0*t/endt .lt. .5) then
            islandr = islandr + 50./(endt/2-1)
        else
            islandr = islandr - 50./(endt/2-1)

        end if
        islandr = 10
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
        torque = kappa*(sin(x(2)-x(1))+sin(x(2)-x(3))+sin(y(2)-y(1))+sin(y(2)-y(3)))-mu*sin( (theta-fieldphi) )
    end function torque

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


    subroutine update(state,newstate,lnoise,n)
        real(8) :: theta, thetap,x,dele,e1,e2,lnoise,fluce,u,force,fieldphi,mutemp
        integer :: n,ii,jj,i,j
        real(8), dimension(n,n) :: state,newstate
        do i=1,n
            do j=1,n
                call random_number(x)
                ii = 1+floor((n)*x)
                call random_number(x)
                jj = 1+floor((n)*x)
                !see if we are on the boundary
                if ( (sqrt(real((ii-n/2)**2 + (jj-n/2)**2)) .lt. islandr) .and. &
                    (sqrt(real((ii-n/2)**2+(jj-n/2)**2)) .gt. (islandr-deltar))) then
                    mutemp = muvector(t)
                    fieldphi = atan2(real(jj-n/2),real(ii-n/2))+pi/2
                   ! write(*,*) ii,jj, fieldphi, mutemp 
                else
                    mutemp = mu
                endif


                theta = grid(ii,jj)
                force = torque(ii,jj,theta,kappa,mutemp,fieldphi)
                call random_number(fluce)
                fluce = fluce-.5
                newstate(ii,jj)=theta-deltime*(fluce*lnoise+force)
                avedelphi = avedelphi + deltime*fluce*lnoise
                dvar = dvar+(fluce*lnoise)**3
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
end program langan
