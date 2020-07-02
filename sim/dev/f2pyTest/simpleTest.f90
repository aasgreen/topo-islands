module defects

contains
subroutine simpleTest(A,n, n2, n4)
    ! testing out the f2py funcitonality
    ! 1. Create an nxn array and return it
    ! 2. iterate on the array using update

    implicit none
    integer, intent(in) :: n, n2, n4
    real(8), intent(out), dimension(n,n) :: A

    A = 1
    return
    end subroutine simpleTest

    subroutine update(newstate,state,hgrid,lnoise,n)
        implicit none
    real(8) :: theta, thetap,x,dele,e1,e2,lnoise,fluce,u,force,fieldphi,mutemp
    integer :: n,ii,jj,i,j
    real(8), dimension(n,n) :: state,newstate,hgrid

    intent(in) :: state, lnoise, n
    intent(out) :: newstate, hgrid

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
 
end module defects



