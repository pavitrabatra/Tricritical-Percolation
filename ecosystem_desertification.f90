!using case function 
!using do cycle instead of go to 
!working perfectly

program ecosystem
implicit none
real*8::p,q                                     !control parameters
real*8::r1,r2,r3,r4,r5,r6                       !random numbers
integer::L,S,n,li                               !
integer::i,j,m,k,mm                             !counters
integer,dimension(:,:),allocatable:: matrix     !lattice
character(len=256)::fname                       !file name variation
integer::un                                     !file unit

print*, 'Enter the size of the lattice'
read*, L

print*, 'Environmental Factor, p'
read*,p

print*, 'Feedback Factor, q'
read*,q

print*, 'Number of time steps, n'
read*, n

li = L - 1

allocate(matrix(0:li, 0:li))


! Getting a Microstate (uniform distribution)
open(newunit=un, file='microstate.dat')
do i = 0, li
do j = 0, li
    call random_number(r1)
    matrix(i, j) = floor(2*r1)
end do
end do
do k = 0, li
    write (un, *) matrix(k, :)
end do
close(un)

!Time Steps
do m = 0, n-1
print *, 'time', m
!N^2 Iteration
do mm = 1, L**2
print *, mm
!Choosing a random site
do
    !Check if everyone died
    do i = 0, li
    do j = 0, li
        if (matrix(i, j) == 1) goto 10
    end do
    end do
    print *, 'everything 0 at time step', m
    call writeoutput(m)
    stop

10  call random_number(r2)
    i=floor(r2*L)
    call random_number(r3)
    j=floor(r3*L)
    if (matrix(i,j)==0) cycle
                
    call random_number(r4)
    select case(floor(4*r4))
    !top neighbour selected
    case(0)
        if (matrix(i, modulo(j+1, L)) == 0) then
            call random_number(r5)
            if (r5<p) then 
                matrix(i, modulo(j+1, L)) = 1
            else
                matrix(i,j) = 0
                cycle
            end if
        else
            call random_number(r5)
            if (r5<q) then
                call random_number(r6)
                select case(floor(6*r6))
                !choosing top neigbhour of the pair     
                case(0)
                    matrix(i, modulo(j+2, L)) = 1
                !choosing right top neighbour of the pair
                case(1)
                    matrix(modulo(i+1,L), modulo(j+1,L)) = 1
                !choosing right bottom neighbour of the pair
                case(2)
                    matrix(modulo(i+1,L), j)=1
                !choosing bottom neighbour of the pair
                case(3)
                    matrix(i, modulo(j-1,L))=1
                !choosing left bottom neighbour of the pair
                case(4)
                    matrix(modulo(i-1,L),modulo(j+1,L))=1
                !choosing left top neighbour of the pair
                case(5)
                    matrix(modulo(i-1,L),modulo(j+1,L))=1
                end select
            else
                matrix(i,j)=0
            end if
        end if
    !bottom neighbour selected
    case(1)
        if (matrix(i,modulo(j-1,L)) == 0) then
            call random_number(r5)
            if (r5<p) then 
                matrix(i,modulo(j-1,L))=1
            else
                matrix(i,j)=0
                cycle
            end if
        else
            call random_number(r5)
            if (r5<q) then
                call random_number(r6)
                select case(floor(6*r6))
                !choosing top neigbhour of the pair     
                case(0)
                    matrix(i,modulo(j+1,L)) = 1
                !choosing right top neighbour of the pair
                case(1)
                    matrix(modulo(i+1,L),j)=1
                !choosing right bottom neighbour of the pair
                case(2)
                    matrix(modulo(i+1,L),modulo(j-1,L))=1
                !choosing bottom neighbour of the pair
                case(3)
                    matrix(i,modulo(j-2,L))=1
                !choosing left bottom neighbour of the pair
                case(4)
                    matrix(modulo(i-1,L),modulo(j-1,L))=1
                !choosing left top neighbour of the pair
                case(5)
                    matrix(modulo(i-1,L),j)=1
                end select
            else
                matrix(i,j)=0
            end if
        end if
    !left neighbour selected
    case(2)
        if (matrix(modulo(i-1,L),j) == 0) then
            call random_number(r5)
            if (r5<p) then 
                matrix(modulo(i-1,L),j)=1
            else 
                matrix(i,j)=0
                cycle
            end if
        else
            call random_number(r5)
            if (r5<q) then
                call random_number(r6)
                select case(floor(6*r6))
                !choosing left neigbhour of the pair     
                case(0)
                    matrix(modulo(i-2,L),j) = 1
                !choosing top left neighbour of the pair
                case(1)
                    matrix(modulo(i-1,L),modulo(j+1,L))=1
                !choosing top right neighbour of the pair
                case(2)
                    matrix(i,modulo(j+1,L))=1
                !choosing right neighbour of the pair
                case(3)
                    matrix(modulo(i+1,L),j)=1
                !choosing bottom left neighbour of the pair
                case(4)
                    matrix(modulo(i-1,L),modulo(j-1,L))=1
                !choosing bottom right neighbour of the pair
                case(5)
                    matrix(i,modulo(j-1,L))=1
                end select
            else
                matrix(i,j)=0
            end if
        end if
    !right neighbour selected
    case(3)
        if (matrix(modulo(i+1,L),j) == 0) then
            call random_number(r5)
            if (r5<p) then 
                matrix(i,modulo(j+1,L))=1
            else
                matrix(i,j)=0
                cycle
            end if
        else
            call random_number(r5)
            if (r5<q) then
                call random_number(r6)
                select case(floor(6*r6))
                !choosing left neigbhour of the pair     
                case(0)
                    matrix(modulo(i-1,L),j) = 1
                !choosing top left neighbour of the pair
                case(1)
                    matrix(i,modulo(j+1,L))=1
                !choosing top right neighbour of the pair
                case(2)
                    matrix(modulo(i+1,L),modulo(j+1,L))=1
                !choosing right neighbour of the pair
                case(3)
                    matrix(modulo(i+2,L),j)=1
                !choosing bottom left neighbour of the pair
                case(4)
                    matrix(i,modulo(j-1,L))=1
                !choosing bottom right neighbour of the pair
                case(5)
                    matrix(modulo(i+1,L),modulo(j-1,L))=1
                end select
            else
                matrix(i,j)=0
            end if
        end if
    end select
    exit
end do
end do 
        call writeoutput(m)
end do

contains
        subroutine writeoutput(t)
                integer, intent(in) :: t
                integer :: i

                write(fname, "(a, i0, a)") 'update_', t, '.dat'
                open(newunit=un, file=fname)
                do i = 0, li
                    write (un, *) matrix(i, :)
                end do
                close(un)
        end subroutine writeoutput
end program ecosystem
