module langan
!functions and classes for running the langegin-ginzburg simulations for
!defects. I think seperating this from the main program will be nice
!as it will allow these functions to be called directly using f2py when 
!we need to test them.

    !define global constants
    real(8),parameter :: pi = 4*atan(1.0_8),  deltime=0.01
    real(8) :: t1=1.2,t2=1.0,dt,dvar=0, islandr=20, deltar=.5

    integer :: endt=100



    contains


    subroutine initializeGrid(grid,n)
        !initialize the grid of spins

        !Args: 
        !   grid: the input grid (already allocated
        !   n: the size of the grid (assume square)

        !Returns: 
        !   The new grid, overwriting the old grid

        real(8), dimension(n,n), intent(inout) :: grid
        integer, intent(in) :: n
        real(8) :: xi = 10.d0 ! the initial correlation length

        integer :: initScheme

        initScheme = 1

        select case (initScheme)
        case DEFAULT
            !===================================
            ! create a grid of spins all pointing to the right

            grid = 0

        case (1)
        !============================Case 1=======================
        !Use the solution to the sinGordon equation (tahn) to create a 1D domain wall in the spin

            do i = 1,n
                do j = 1,n
                    grid(i,j) = 2*atan( exp(-real(i-n/2)/xi) ) ! pi -> 0 with kink in the middle

                end do
            end do
        end select

    end subroutine initializegrid

    function muvector(t, kappa)
        real(8) :: muvector
        real(8) :: kappa
        intent(in) :: kappa, t

        muvector = .5*kappa
    end function muvector

    subroutine islandradius(islandr,t)
        !calculate the radius of the island at a given time
        integer :: t, islandScheme
        real(8) :: islandr

        intent(in) :: t
        intent(out) :: islandr

        !first, assign how you want this function to be called
        islandScheme = 2
        ! case 2: create island out of nothing at .2 percent
        select case (islandScheme)

        !=================Default=======================
        ! create an island and let it sit
        case DEFAULT
            if ( (t .ge. 100) ) then
                islandr = 20 !this will put the line in the middle of the simulation
            else
                islandr = 0
            end if

        !==================Case 1==========================
        ! creating an island out of nowhere and disappearing it
        case (1)
            if ( (t .ge. 100) .and. (t .le. 400 ) ) then
                islandr = 20 !this will put the line in the middle of the simulation
            else
                islandr = 0
            end if

        !=================Case 2 ===========================
        !grow an island from nothing, then shrink it
        case (2) 
            if (1.0*t/endt .lt. .3) then
                islandr = islandr + 30./(endt/3)
            else if ( (1.0*t/endt .ge. .3) .and. (1.0*t/endt .le. .6) ) then
                islandr = islandr - 200./(endt/3)
            else
                islandr = -20

            end if
        end select

    end subroutine islandradius


    function torque(i,j,grid,n,kappa,mu,fieldphi)
        !calculate the torque for each grid position
       ! Args: 
        !   i: row co-ordinate
       !    j: column co-ordinate
    !       theta: spin at point grid(i,j)
   !        kappa: elastic constant
  !         mu: electric field strength
 !          fieldphi: the anchoring condition for the electric field 

 !      Returns:
 !          the force (torque) on the spin at site i,j
        integer :: i,j,n
        real(8):: torque,kappa,fieldphi
        real(8) :: mu
        real(8), dimension(3) :: x,y
        real(8), dimension(n,n) :: grid
        intent(in) :: i, j, n, grid, kappa, mu, fieldphi


        x =(/ grid(modulo(i-2,n)+1,j),grid(i,j),grid(modulo(i,n)+1,j)/)
        y =(/ grid(i,modulo(j-2,n)+1),grid(i,j),grid(i,modulo(j,n)+1)/)
        torque = kappa/8*(sin(x(2)-x(1))+sin(x(2)-x(3))+sin(y(2)-y(1))+& 
            sin(y(2)-y(3)))+4*mu*cos( (grid(i,j)-fieldphi))* &
            sin((grid(i,j)-fieldphi) ) 
    end function torque

    function torqueLandau(i,j,grid,n,kappa,mu,fieldphi)
        !calculate the torqueLandau for each grid position
       ! Args: 
        !   i: row co-ordinate
       !    j: column co-ordinate
    !       theta: spin at point grid(i,j)
   !        kappa: elastic constant
  !         mu: electric field strength
 !          fieldphi: the anchoring condition for the electric field 

 !      Returns:
 !          the force (torqueLandau) on the spin at site i,j
        integer :: i,j,n
        real(8):: torqueLandau,kappa,fieldphi
        real(8) :: mu, alpha = .1
        real(8), dimension(3) :: x,y
        real(8), dimension(n,n) :: grid
        intent(in) :: i, j, n, grid, kappa, mu, fieldphi


        x =(/ grid(modulo(i-2,n)+1,j),grid(i,j),grid(modulo(i,n)+1,j)/)
        y =(/ grid(i,modulo(j-2,n)+1),grid(i,j),grid(i,modulo(j,n)+1)/)
        energy = kinetic(i,j,grid,n,kappa)
        torqueLandau = kappa/8*(1-alpha*energy/kappa)*(sin(x(2)-x(1))+ &
            sin(x(2)-x(3))+sin(y(2)-y(1))+sin(y(2)-y(3)))+4*mu*cos( (grid(i,j)-fieldphi))* &
            sin((grid(i,j)-fieldphi) )
    end function torqueLandau


    function kinetic(i,j,grid,n,kappa)
        integer :: i,j,n
        real(8) :: kappa, kinectic
        real(8), dimension(3) :: x,y
        real(8), dimension(n,n) :: grid
        intent(in) :: i, j, grid, n, kappa

        !write(*,*) 'hamxy', i,j
        x =(/ grid(modulo(i-2,n)+1,j),grid(i,j),grid(modulo(i,n)+1,j)/)
        y =(/ grid(i,modulo(j-2,n)+1),grid(i,j),grid(i,modulo(j,n)+1)/)
        kinetic = -kappa/8*(-4+(cos(x(2)-x(1))+cos(x(2)-x(3))+cos(y(2)-y(1))+cos(y(2)-y(3)))) !shift so energy =0 is min.
    end function kinetic
     
    function hamxy(i,j,grid,n,kappa,mu,fieldphi)
        integer :: i,j,n
        real(8):: fieldphi
        real(8) :: kappa,hamxy,mu
        real(8), dimension(3) :: x,y
        real(8), dimension(n,n) :: grid
        intent(in) :: i, j, grid, n, kappa, mu, fieldphi

        !write(*,*) 'hamxy', i,j
        x =(/ grid(modulo(i-2,n)+1,j),grid(i,j),grid(modulo(i,n)+1,j)/)
        y =(/ grid(i,modulo(j-2,n)+1),grid(i,j),grid(i,modulo(j,n)+1)/)
        hamxy = -kappa*(cos(x(2)-x(1))+cos(x(2)-x(3))+cos(y(2)-y(1))+cos(y(2)-y(3)))-mu*cos( (grid(i,j)-fieldphi) )
    end function hamxy
     
   function anchor(ii,jj,ir,dr,n)
      integer :: n, ii, jj, dr
     real(8) :: anchor,ir

     intent(in) ii, jj, ir, dr, n

    if ( (sqrt(real((ii-n/2)**2 + (jj-n/2)**2)) .lt. ir) .and. &
                    (sqrt(real((ii-n/2)**2+(jj-n/2)**2)) .gt. (ir-dr))) then

                    anchor = atan2(real(jj-n/2),real(ii-n/2))+pi/2
                else
                    anchor = 0.d8
                endif

        end function anchor



    subroutine update(state,newstate,hgrid,lnoise,kappa, mu, n)
        real(8) :: theta, x,lnoise,fluce,force,fieldphi,mutemp, kappa, mu
        integer :: n,ii,jj,i,j
        real(8), dimension(n,n) :: state,newstate,hgrid

        intent(in) :: state, n, lnoise, kappa, mu
        intent(out) newstate, hgrid

        do i=1,n
            do j=1,n
                call random_number(x)
                ii = 1+floor((n)*x) !this will give periodic boundaries
                call random_number(x)
                jj = 1+floor((n)*x)
                !see if we are within the line boundary
                mutemp = mu
                fieldphi = 0

                theta = state(ii,jj)
                force = torqueLandau(ii,jj,state, n ,kappa,mutemp,fieldphi)
                call random_number(fluce)
                fluce = fluce-.5
                newstate(ii,jj)=theta-deltime*(fluce*lnoise+force)
                hgrid(ii,jj) = hamxy(ii, jj, newstate, n , kappa, mu, fieldphi)+4*kappa !this will normalize so zero is actually the
        !energy
                avedelphi = avedelphi + deltime*fluce*lnoise
                dvar = dvar+(fluce*lnoise)**2
                dcount = dcount+1
            end do
        end do
    end subroutine update
            

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


    end module langan
