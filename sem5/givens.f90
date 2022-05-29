subroutine compute_params(aii, aji, c, s)
    real, intent(in):: aii, aji
    real::c, s
    c = aii / sqrt(aii * aii + aji * aji)
    s = -aji / sqrt(aii * aii + aji * aji)
end subroutine compute_params

subroutine rotate(xi, xj, c, s)
    real:: xi, xj, c, s
    real:: xi_, xj_
    xi_ = xi * c - xj * s
    xj_ = xi * s + xj * c
    xi = xi_
    xj = xj_
end subroutine rotate

subroutine qr(a, q, size)
    real, dimension(:,:):: a
    real, allocatable, dimension(:), intent(out):: q
    integer:: size
    integer:: iq, icol, irow, k
    real:: c, s
    allocate(q(size*(size - 1)))
    iq = 1
    icol = 1
    do while (icol < size)
        irow = icol + 1
        do while(irow < size + 1)
            call compute_params(a(icol, icol), a(irow, icol), c, s)
            q(iq) = c
            q(iq+1) = -s
            iq = iq + 2
            k = icol
            do while (k < size + 1)
                call rotate(a(icol, k), a(irow, k), c, s)
                k = k + 1
            end do
            irow = irow + 1
        end do
        icol = icol + 1
    end do
end subroutine qr


program givens
    implicit none
    integer, dimension(4):: sizes = [256, 512, 1024, 2048]
    integer:: i, n
    real, allocatable, dimension(:,:):: a
    real, allocatable, dimension(:):: q
    real(4):: start_time, finish_time

    interface
        subroutine qr(a, q, size)
            real, dimension(:,:):: a
            real, allocatable, dimension(:), intent(out):: q
            integer:: size
        end subroutine qr
    end interface
    
    i = 1
    do while(i < 5)
        n = sizes(i)
        allocate(a(n, n))
        call random_number(a)
        call cpu_time(start_time)
        call qr(a, q, n)
        call cpu_time(finish_time)
        deallocate(a)
        print *, finish_time - start_time
        i = i + 1
    end do
end program givens