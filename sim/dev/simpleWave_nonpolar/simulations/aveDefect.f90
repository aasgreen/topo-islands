program aveDefect

    character(100) :: buffer
    real, allocatable, dimension(:,:) :: x
    real(8), dimension(10) :: xtest
    complex, dimension(300*300) :: defect
    real :: aveDNN=0,temp=0,aveD=0
    real , dimension(4):: ll
    real(8), allocatable, dimension(:) :: defectD
    integer :: time,row, col,indx=1, defectN=0,N

    character(100) :: fileNames


    call getarg(1,buffer)
    read(buffer,"(A100)"), fileNames

    call getarg(2,buffer)
    read(buffer,*), N

    call getarg(3,buffer)
    read(buffer,*), time
   ! write(*,*) fileNames
    allocate(x(N,N))
    open(unit=99, file=fileNames,status='old',action='read')
    do row = 1, N
        read(99,*) (x(row,col), col=1,N)
    enddo
    do row=1,N
        do col=1,N
            if (x(row,col) .eq. 1 .or. x(row,col) .eq. -1) then
                defect(indx) = cmplx(row,col)
!                write(*,*) defect(indx), x(row,col), indx
                indx=indx+1
            endif
        enddo
    enddo
    defectN=indx-1
    indx=1
    !write(*,*) defectN
    if (defectN .lt. 3) then
        write(*,*)time, N, defectN
        stop
    endif
    

    !write(*,*) 'allocate'
    allocate(defectD(defectN))
    do row=1,defectN
        do col = 1,defectN
            !write(*,*) row,col
            if (row .ne. col) then
                defectD(col)=modulo(cabs(defect(row)-defect(col))-1,real(N/2))+1
                if (int(defectD(col)) == 0) then
                    write(*,*) "Distance error", row, defect(row),col,defect(col)
                    stop
                endif
            endif
        enddo
        ll=sortit(defectD)
        !ll=0
    !    write(*,*) 'smallest ds: ', ll
        aveDNN=aveDNN+ll(1)
        aveD = aveD+SUM(defectD)/defectN
        defectD=0
    enddo
    !write(*,*) 'average distance', aveD/defectN
    !write(*,*) 'defect number', defectN
    write(*,*)time, aveDNN/real(defectN), aveD/real(defectN), defectN



contains
    function dist(a,b)
        real ::a,b
        integer :: i,j
        dist = sqrt(real(a)**2+real(b)**2 )
    end function dist

    function sortit(array)
        implicit none
        real(8) :: array(:)
        integer,dimension(4) :: kL=-1
        real(8), dimension(4) :: sortit
        real(8) :: smallVal=huge(0.0),bigVal=huge(0.0)
        integer ::j,k,l,ksmol=0,rep=0
        sortit(1) = smallVal
        do k =1,ubound(array,1)
            if(array(k) < sortit(1) .and. array(k) > 0) then
                sortit(1) =array(k)
                kL(1) = k
            endif
        enddo 

        do j=2, 4
            smallVal = bigVal
            do k=1, ubound(array,1)
                if(array(k) < smallVal .and. array(k) .ge. sortit(j-1)) then
                    rep=0
        
                    do l=1,j-1
                        if (kL(l) ==k) then
                           rep = 1
                       endif
                    enddo
                    if (rep ==1) then
                        cycle

                    endif
                    !if not duplicate k, then update ksmol
                    sortit(j) =array(k)
                    smallVal=array(k)
                    kL(j)=k
                endif
            enddo
        enddo
        end function sortit
end program aveDefect
