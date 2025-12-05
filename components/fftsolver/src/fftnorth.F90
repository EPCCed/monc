! Module that provides wrapper functions for FFTPACK calls

module fftnorth_mod
    use fftpack, only: zffti, zfftf, zfftb
    use fftpack_kind, only: rk

    implicit none

    ! The array that will have the in-place FFT applied to it
    complex(rk), allocatable, dimension(:) :: data

    real(rk), allocatable, dimension(:) :: wsave

    contains
    
    ! Initialises the FFT routines for the real problem size (n)
    subroutine fftn_init(n)
        integer, intent(in) :: n

        allocate(data(n))
        allocate(wsave(4*n+15))
        call zffti(n, wsave)
    end subroutine
    
    ! Computes a real-to-complex (e.g. forward) FFT
    ! in  : double precision real array of size n (input)
    ! out : double precision complex array of size n/2+1 (output)
    ! n   : integer - size of in (input)
    subroutine fftn_r2c(in, out, n)
        integer, intent(in) :: n
        double precision, intent(in) :: in(n)
        complex*16, intent(out) :: out(n/2+1)

        integer :: i

        !copy real input into complex "data" array
        do i=1,n
            data(i) = dcmplx(in(i), 0.d0)
        enddo

        call zfftf(n,data,wsave)

        !extract the first n/2+1 terms from the FFT and return them as output
        !(The last n/2-1 terms of r2c are complex conjugates of the previous 
        ! ones and so are redundant)
        out(:) = data(1:n/2+1)

    end subroutine
    
    ! Computes a complex-to-real (e.g. inverse) FFT
    ! in  : double precision complex array of size n/2+1 (input)
    ! out : double precision real array of size n (output)
    ! n   : integer - size of out (input)
    subroutine fftn_c2r(in,out,n)
        integer, intent(in) :: n
        complex*16, intent(in) :: in(n/2+1)
        double precision, intent(out) :: out(n)

        integer :: i

        !copy the first n/2+1 terms
        data(1:n/2+1) = in(:)
        !set the remaining entries to the complex conjugates of the previous ones
        do i=n/2+2,n
            data(i) = dconjg(in(n-i+2))
        enddo
        
        call zfftb(n,data,wsave)
        !extract the real part of data and place it in out
        do i=1,n
            out(i) = real(data(i))/n
        enddo
    end subroutine
    
    ! Finalises the FFT routines
    subroutine fftn_finalise()
        ! Deallocate the work arrays
        deallocate(data)
        deallocate(wsave)
    end subroutine

end module fftnorth_mod
