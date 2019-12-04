program lanGan
implicit none

character(100) :: buffer
integer :: N = 200, endT=1001,seedN
real(8) :: beta,mu,measuredT,zeroE, meanK,avedelPhi
real(8), allocatable, dimension(:) :: muVector
integer, dimension(100) :: tPoints
integer, allocatable, dimension(:) :: logTPoints
real(8), allocatable, dimension(:,:) :: grid,dgrid,gridPlusDelta
real(8) :: x,kappa,windingN,lNoise
real(8),parameter :: PI = 4*atan(1.0_8),  delTime=0.1
integer, allocatable :: seed(:)
character(40) :: fileNames
character(40) :: dfileNames

real(8) :: t1=1.2,t2=1.0,dT,dVar=0,islandR
integer :: i,j,t,timePrint=1,tt,dCount=0,un=4,seedSize,istat=0,deltaR

!read in arguments (g,beta,mu,N,endT)
call getarg(1,buffer)
read(buffer,*) kappa

call getarg(2,buffer)
read(buffer,*) beta

call getarg(3,buffer)
read(buffer,*) mu

call getarg(4,buffer)
read(buffer,*) N

call getarg(5,buffer)
read(buffer,*) endT

call getarg(6,buffer)
read(buffer,*) seedN

call random_seed(seedSize)
allocate(seed(seedSize))

  ! First try if the OS provides a random number generator
open(newunit=un, file="/dev/urandom", access="stream", &
       form="unformatted", action="read", status="old", iostat=istat)
if (istat == 0) then
    read(un) seed
    close(un)
endif

!calculate zero temperature energy
zeroE = -1*kappa*4.0

!calculate langenvin noise term (c_L in Yurke)
!lNoise=sqrt(24*1./beta*delTime)
lNoise = sqrt(12*1/beta)
do t=1,100
    dT = (LOG10(REAL(endT/delTime))/100)
    tPoints(t) =  INT(10**(t*dT))
enddo
logTPoints= tPoints(unique(tPoints))
!write(*,*) logTPoints


!Initial Random Grid
!write(*,*) 'initialize grid'
allocate(grid(N,N))
allocate(gridPlusDelta(N,N))
!write(*,*) grid(N,N)




call random_seed(put=seed)
do i=1,N
    do j=1,N
!        write(*,*) i,j
        call random_number(x)
        grid(i,j) = 0
    end do
enddo
gridPlusDelta = grid
measuredT = 0
do i=1,N
    do j=1,N
        measuredT = (hamXY(i,j,grid(i,j),kappa,mu))/N/N+measuredT
    enddo
enddo

!calculate the indices where the boundary exists. These indicies will meet the
!equation i^2+j^2 = R^2, where R is the radius of the circle. We will take some 
!delta, so (R-del)^2 .ge. i^2+j^2 .ge. (R)^2
islandR =6
deltaR = 3

!Allocate mu vector
allocate(muVector(int(endT/delTime)))
muVector = 10*kappa
muVector(1:10)=0
!muVector(int(endT/20/delTime):int(endT/delTime)) = (/ (mu+kappa*i/(endT/delTime-2)*100/2, i= 0,int((endT/delTime-1)/2)) /)
muVector(int(endT/20/delTime):int(endT/delTime)) = 10*kappa
write(*,*) muVector



write(*,*) 'max T', measuredT-zeroE

!open temperature V time file
open(61,file='tVT.dat',status = 'unknown', position='append')

!initialize defect grid
allocate(dgrid(N,N))
dgrid=0
endT = logTPoints(size(logTPoints))
write(*,*) endT
!write(*,*) 'BEGINING LANGAN PROCESSING...'
!preform langan with euler update
do t=1,int(endT)

    avedelPhi = 0
    grid=gridPlusDelta
    if (( t .gt. 50 )) then
        call islandRadius(islandR,t)
    endif
    call update(grid,gridPlusDelta, lNoise,N)
    meanK =.5*SUM((grid-gridPlusDelta)**2)/N/N/delTime**2
    write(fileNames,'(A3,F0.4, A4)') 'out', t*delTime,'.dat'
    write(dfileNames,'(A6,F0.4, A4)') 'defect', t*delTime,'.dat'

    if (logSpace(t) .eq. 1) then
        print *, trim(fileNames)
    end if
    if (1 .eq. 1) then !this is not log time, so we need all time steps
       !print *, trim(fileNames)

        open(1,file=fileNames)
        open(3,file=dfileNames)

        ! calculate defects and average energy, and write grid to file
        measuredT = 0.

        do i=1,N
            do j=1,N
            write(1,'(F10.5)',advance="no") grid(i,j)
            measuredT = (hamXY(i,j,grid(i,j),kappa,mu))/N/N+measuredT
            windingN=(windN(grid(i,modulo(j-2,N)+1)-grid(i,j)))+&
                &(windN(grid(modulo(i,N)+1,modulo(j-2,N)+1)-grid(i,modulo(j-2,N)+1)))+&
                &(windN(grid(modulo(i,N)+1,j)-grid(modulo(i,N)+1,modulo(j-2,N)+1)))+&
                &(windN(grid(i,j)-grid(modulo(i,N)+1,j)))
            if (windingN .ge. 1) then
                dgrid(i,j) =1
            else if (-1*windingN .ge. 1) then
                dgrid(i,j)=-1
            else
                dgrid(i,j)=0
            endif
           write(3,'(F10.5)',advance="no") dgrid(i,j)
            enddo
            write(1,*)
            write(3,*)
        enddo
        measuredT = (measuredT-zeroE)
        write(61,'(F12.2,A, F10.5)') t*delTime,',', measuredT, meanK
        !write(*,*) t*delTime

        close(1)
        close(3)
    endif
!write(*,*) 'closing'
end do
write(*,*) dVar/dCount/2
close(61)
!Write to File
write(*,*) 'finished langin'

deallocate(grid)
contains

    subroutine islandRadius(islandR,t)
        !calculate the radius of the island at a given time
        integer :: t
        real(8) :: islandR
        ! case 2: grow, hit 50 and then shrink
        if (1.0*t/endT .lt. .5) then
            islandR = islandR + 50./(endT/2-1)
        else
            islandR = islandR - 50./(endT/2-1)

        end if
    end subroutine islandRadius

    function logSpace(t)
        integer :: t, logSpace
        if (t .eq. logTPoints(timePrint)) then
            logSpace = 1
            timePrint = timePrint +1
        else
            logSpace = 0
        endif
    end function logSpace

    function TORQUE(i,j,theta,kappa,mu,fieldPhi)
        integer :: i,j,ii,jj
        real(8):: theta,torque,kappa,fieldPhi
        real(8) :: g,hamXY,mu
        real(8), dimension(3) :: x,y
        !write(*,*) 'hamxy', i,j
        x =(/ grid(modulo(i-2,N)+1,j),theta,grid(modulo(i,N)+1,j)/)
        y =(/ grid(i,modulo(j-2,N)+1),theta,grid(i,modulo(j,N)+1)/)
        torque = kappa*(sin(x(2)-x(1))+sin(x(2)-x(3))+sin(y(2)-y(1))+sin(y(2)-y(3)))-mu*sin( (theta-fieldPhi) )
    end function torque

    function xi(grid,beta,N)
        integer :: N
        real(8), dimension(N,N) :: grid
        real(8) :: beta, xi
        xi = sum(grid**2)/N**2
    end function xi

    function hamXY(i,j,theta,kappa,mu)
        integer :: i,j,ii,jj
        real(8):: theta
        real(8) :: kappa,hamXY,mu
        real(8), dimension(3) :: x,y
        !write(*,*) 'hamxy', i,j
        x =(/ grid(modulo(i-2,N)+1,j),theta,grid(modulo(i,N)+1,j)/)
        y =(/ grid(i,modulo(j-2,N)+1),theta,grid(i,modulo(j,N)+1)/)
        hamXY = -kappa*(cos(x(2)-x(1))+cos(x(2)-x(3))+cos(y(2)-y(1))+cos(y(2)-y(3)))-mu*cos( (theta-45/2/PI) )
    end function hamXY


    subroutine update(state,newstate,lNoise,N)
        real(8) :: theta, thetaP,x,delE,E1,E2,lNoise,flucE,u,force,fieldPhi,muTemp
        integer :: N,ii,jj,i,j
        real(8), dimension(N,N) :: state,newstate
        do i=1,N
            do j=1,N
                call random_number(x)
                ii = 1+floor((N)*x)
                call random_number(x)
                jj = 1+floor((N)*x)
                !see if we are on the boundary
                if ( (sqrt(real((ii-N/2)**2 + (jj-N/2)**2)) .lt. islandR) .and. &
                    (sqrt(real((ii-N/2)**2+(jj-N/2)**2)) .gt. (islandR-deltaR))) then
                    muTemp = muVector(t)
                    FieldPhi = atan2(real(jj-N/2),real(ii-N/2))+PI/2
                   ! write(*,*) ii,jj, fieldPhi, muTemp 
                else
                    muTemp = mu
                endif


                theta = grid(ii,jj)
                force = torque(ii,jj,theta,kappa,muTemp,FieldPhi)
                call random_number(flucE)
                flucE = flucE-.5
                newstate(ii,jj)=theta-delTime*(flucE*lNoise+force)
                avedelPhi = avedelPhi + delTime*flucE*lNoise
                dVar = dVar+(flucE*lNoise)**2
                dCount = dCount+1
            end do
        end do
    end subroutine update
            

    function dtrack(dgrid,grid)
        real(8), dimension(N,N) :: dgrid
        real(8), dimension(3,3) :: grid
        real(8)                 :: angle=0
        integer                 :: dtrack,i,j
        do i=1,3
            !angle=angle+grid(i,j)-grid(2,2)/9., j=1,3)
            !write(*,*) (grid(i,j), j=1,3)
        enddo
        if (angle .ge. 2*PI) then
            dtrack=1

        else  
            dtrack=0
        endif

    endfunction dtrack

    function angleDist(theta1,theta2)
        real(8) :: d1,d2,angleDist,theta1,theta2
        angleDist = acos(cos(theta1-theta2))
        end function angleDist

    function windN(angle)
        real(8) :: angle, windN
        windN = (angle-asin(sin(angle))) / (PI) 
        end function windN



function remove_dups(input)
  integer :: input(100)       ! The input
  integer :: remove_dups(size(input))  ! The output
  integer :: k                   ! The number of unique elements
  integer :: i, j
 
  k = 1
  remove_dups(1) = input(1)
  outer: do i=2,size(input)
     do j=1,k
        if (remove_dups(j) == input(i)) then
           ! Found a match so start looking again
           cycle outer
        end if
     end do
     ! No match found so add it to the output
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
end program lanGan
