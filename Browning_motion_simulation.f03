! Simulation of Browning motion.

! The MSD(t) is defined as:
! MSD(t) = <|r(t) - r(0)|^2> .

! There should be a relationship like:
! MSD(t) / t = slope = const .

! The coefficient of diffusivity:
! D = slope / (2 * ndim)

subroutine random_init_seed()
    implicit none
    integer(kind=4) :: rand_seed_size
    integer(kind=4), allocatable :: rand_seed(:)
    integer(kind=4) :: time_c

    call system_clock(time_c)
    call random_seed(size = rand_seed_size)
    allocate(rand_seed(rand_seed_size))
    call random_seed(get = rand_seed)
    rand_seed = rand_seed + time_c
    call random_seed(put = rand_seed)
    deallocate(rand_seed)

    return
end subroutine random_init_seed

subroutine std_norm_dist(value, method)
    implicit none
    real(kind=8), intent(out) :: value
    integer(kind=4), intent(in), optional :: method
    integer(kind=4) :: internal_method
    real(kind=8) :: u, v
    real(kind=8) :: s
    real(kind=8) :: tmp

    if (present(method)) then
        internal_method = method
    else
        internal_method = 1
    end if
    do while (.true.)
        call random_number(u)
        call random_number(v)
        u = u * 2.0D0 - 1.0D0
        v = v * 2.0D0 - 1.0D0
        s = u ** 2 + v ** 2
        if (s .lt. 1.0D0) exit
    end do
    tmp = dsqrt(- 2.0D0 * dlog(s) / s)
    if (mod(internal_method, 2) .eq. 1) then
        value = u * tmp
    else
        value = v * tmp
    end if

    return
end subroutine std_norm_dist

subroutine std_norm_dist_useboth(value)
    implicit none
    real(kind=8), intent(out) :: value
    logical(kind=1), save :: status_std_norm_dist = .true.
    real(kind=8), save :: old_std_norm_dist
    real(kind=8) :: u, v
    real(kind=8) :: s
    real(kind=8) :: tmp

    status_std_norm_dist = .not. status_std_norm_dist
    if (status_std_norm_dist) then
        do while (.true.)
            call random_number(u)
            call random_number(v)
            u = u * 2.0D0 - 1.0D0
            v = v * 2.0D0 - 1.0D0
            s = u ** 2 + v ** 2
            if (s .lt. 1.0D0) exit
        end do
        tmp = dsqrt(- 2.0D0 * dlog(s) / s)
        value = u * tmp
        old_std_norm_dist = v * tmp
    else
        value = old_std_norm_dist
    end if

    return
end subroutine std_norm_dist_useboth

real(kind=8) function get_sd(ndim, nsize, arr, pos, interval)
    implicit none
    integer(kind=4), intent(in) :: ndim
    integer(kind=4), intent(in) :: nsize
    real(kind=8), intent(in) :: arr(ndim, nsize)
    integer(kind=4), intent(in) :: pos, interval

    get_sd = sum((arr(:, pos + interval) - arr(:, pos)) ** 2)

    return
end function get_sd

real(kind=8) function get_msd(ndim, nsize, arr, interval)
    ! a better solution is to use the fft method.
    implicit none
    interface
        real(kind=8) function get_sd(ndim, nsize, arr, pos, interval)
            implicit none
            integer(kind=4), intent(in) :: ndim
            integer(kind=4), intent(in) :: nsize
            real(kind=8), intent(in) :: arr(ndim, nsize)
            integer(kind=4), intent(in) :: pos, interval
        end function get_sd
    end interface
    integer(kind=4), intent(in) :: ndim
    integer(kind=4), intent(in) :: nsize
    real(kind=8), intent(in) :: arr(ndim, nsize)
    integer(kind=4), intent(in) :: interval
    integer(kind=4) :: pos

    get_msd = 0.0D0

    do pos = 1, nsize - interval
        get_msd = get_msd + get_sd(ndim, nsize, arr, pos, interval)
    end do
    get_msd = get_msd / (nsize - interval)

    return
end function get_msd

subroutine get_all_msd(ndim, nsize, arr, ret)
    implicit none
    interface
        real(kind=8) function get_msd(ndim, nsize, arr, interval)
            integer(kind=4), intent(in) :: ndim
            integer(kind=4), intent(in) :: nsize
            real(kind=8), intent(in) :: arr(ndim, nsize)
            integer(kind=4), intent(in) :: interval
        end function get_msd
    end interface
    integer(kind=4), intent(in) :: ndim
    integer(kind=4), intent(in) :: nsize
    real(kind=8), intent(in) :: arr(ndim, nsize)
    real(kind=8), intent(out) :: ret(nsize)
    integer(kind=4) :: j

    if (nsize .ge. 1) ret(1) = 0.0D0
    do j = 2, nsize
        ret(j) = get_msd(ndim, nsize, arr, j - 1)
    end do

    return
end subroutine get_all_msd

subroutine get_all_msd_fft(ndim, nsize, arr, ret, d, s1, s2, &
            tmp_fft_inp, tmp_fft_out, tmp_fft_tmp, plan_forward, plan_backward)
    use, intrinsic :: iso_c_binding
    implicit none
    include "fftw3.f03"
    integer(kind=4), intent(in) :: ndim
    integer(kind=4), intent(in) :: nsize
    real(kind=8), intent(in) :: arr(ndim, nsize)
    real(kind=8), intent(out) :: ret(nsize)
    real(kind=8) :: q
    integer(kind=4) :: i
    integer(kind=4) :: j
    real(kind=8), intent(inout) :: d(nsize)
    real(kind=8), intent(inout) :: s1(nsize)
    real(kind=8), intent(inout) :: s2(nsize)
    real(C_DOUBLE) :: tmp_fft_inp(2 * nsize)
    real(C_DOUBLE) :: tmp_fft_out(2 * nsize)
    complex(C_DOUBLE_COMPLEX) :: tmp_fft_tmp(2 * nsize)
    type(C_PTR) :: plan_forward
    type(C_PTR) :: plan_backward

    d = sum(arr ** 2, dim = 1)
    s2 = 0.0D0
    do i = 1, ndim
        tmp_fft_inp(nsize + 1:2 * nsize) = 0.0D0
        tmp_fft_inp(:nsize) = arr(i, :)
        call fftw_execute_dft_r2c(plan_forward, tmp_fft_inp, tmp_fft_tmp)
        tmp_fft_tmp = tmp_fft_tmp * conjg(tmp_fft_tmp)
        call fftw_execute_dft_c2r(plan_backward, tmp_fft_tmp, tmp_fft_out)
        do j = 1, nsize
            tmp_fft_out(j) = tmp_fft_out(j) / dble(nsize + 1 - j)
        end do
        tmp_fft_out = tmp_fft_out / (2 * nsize)
        s2 = s2 + tmp_fft_out(:nsize)
    end do
    q = sum(d) * 2.0D0
    do j = 1, nsize
        if (j /= 1) q = q - (d(j - 1) + d(nsize + 2 - j))
        s1(j) = q / dble(nsize + 1 - j)
    end do
    ret = s1 - 2.0D0 * s2

    return
end subroutine get_all_msd_fft

subroutine linear_through_origin_fit(arr_size, arr, slope, cod, begin, end)
    implicit none
    integer(kind=4), intent(in) :: arr_size
    real(kind=8), intent(in) :: arr(arr_size)
    real(kind=8), intent(out) :: slope, cod
    integer(kind=4), intent(in), optional :: begin, end
    integer(kind=4) :: fit_b, fit_e
    integer(kind=4) :: i
    real(kind=8), allocatable :: var(:)
    real(kind=8) :: arr_average

    fit_b = nint(dble(arr_size) * 0.1)
    fit_e = nint(dble(arr_size) * 0.9)
    if (present(begin)) fit_b = begin
    if (present(end))   fit_e = end
    allocate(var(fit_b:fit_e))

    var = (/(i, i = fit_b, fit_e)/)
    slope = sum(arr(fit_b:fit_e) * var) / sum(var ** 2)
    arr_average = sum(arr(fit_b:fit_e)) / dble(fit_e - fit_b + 1)

    cod = 1.0D0 - sum((arr(fit_b:fit_e) - var * slope) ** 2) / sum((arr(fit_b:fit_e) - arr_average) ** 2)    

    deallocate(var)

    return
end subroutine linear_through_origin_fit


program main
    use, intrinsic :: iso_c_binding
    implicit none
    include "fftw3.f03"

    interface
        real(kind=8) function get_msd(ndim, nsize, arr, interval)
            implicit none
            integer(kind=4), intent(in) :: ndim
            integer(kind=4), intent(in) :: nsize
            real(kind=8), intent(in) :: arr(ndim, nsize)
            integer(kind=4), intent(in) :: interval
        end function get_msd
        subroutine std_norm_dist(value, method)
            implicit none
            real(kind=8), intent(out) :: value
            integer(kind=4), intent(in), optional :: method
        end subroutine std_norm_dist
        subroutine linear_through_origin_fit(arr_size, arr, slope, cod, begin, end)
            implicit none
            integer(kind=4), intent(in) :: arr_size
            real(kind=8), intent(in) :: arr(arr_size)
            real(kind=8), intent(out) :: slope, cod
            integer(kind=4), intent(in), optional :: begin, end
        end subroutine linear_through_origin_fit
        subroutine get_all_msd(ndim, nsize, arr, ret)
            implicit none
            integer(kind=4), intent(in) :: ndim
            integer(kind=4), intent(in) :: nsize
            real(kind=8), intent(in) :: arr(ndim, nsize)
            real(kind=8), intent(out) :: ret(nsize)
        end subroutine get_all_msd
        subroutine get_all_msd_fft(ndim, nsize, arr, ret, d, s1, s2, &
            tmp_fft_inp, tmp_fft_out, tmp_fft_tmp, plan_forward, plan_backward)
            use, intrinsic :: iso_c_binding
            implicit none
            integer(kind=4), intent(in) :: ndim
            integer(kind=4), intent(in) :: nsize            
            real(kind=8), intent(in) :: arr(ndim, nsize)
            real(kind=8), intent(out) :: ret(nsize)
            real(kind=8), intent(inout) :: d(nsize)
            real(kind=8), intent(inout) :: s1(nsize)
            real(kind=8), intent(inout) :: s2(nsize)
            real(C_DOUBLE) :: tmp_fft_inp(2 * nsize)
            real(C_DOUBLE) :: tmp_fft_out(2 * nsize)
            complex(C_DOUBLE_COMPLEX) :: tmp_fft_tmp(2 * nsize)
            type(C_PTR) :: plan_forward
            type(C_PTR) :: plan_backward
        end subroutine get_all_msd_fft
    end interface

    integer(kind=4) :: argc
    integer(kind=4) :: iarg
    integer(kind=4) :: arg_status
    character(kind=1,len=512) :: argv0
    character(kind=1,len=256), allocatable :: argv(:)

    integer(kind=4) :: nsize    =  1000
    integer(kind=4) :: ndim     =     2
    integer(kind=4) :: stepsize =    10
    integer(kind=4) :: nsimu    = 10000
    real(kind=8), allocatable :: x(:, :) ! starts from 0
    real(kind=8), allocatable :: tmp_arr(:) ! starts from 0
    integer(kind=4) :: t, istep, dt
    integer(kind=4), parameter :: write_to_unit = 11
    character(kind=1, len=26), parameter :: write_to_name = "Browning_motion_result.txt"
    integer(kind=4) :: idim
    integer(kind=4) :: isimu
    real(kind=8), allocatable :: msd(:) ! starts from 0
    real(kind=8), allocatable :: tmp_msd_cur(:) ! starts from 0
    real(kind=8) :: slope, cod, co_diffu

    integer(kind=4) :: time_0
    integer(kind=4) :: time_t

    real(kind=8), allocatable :: d(:)
    real(kind=8), allocatable :: s1(:)
    real(kind=8), allocatable :: s2(:)
    real(C_DOUBLE), pointer :: tmp_fft_inp(:)
    real(C_DOUBLE), pointer :: tmp_fft_out(:)
    complex(C_DOUBLE_COMPLEX), pointer :: tmp_fft_tmp(:)
    type(C_PTR) :: plan_forward
    type(C_PTR) :: plan_backward
    type(C_PTR) :: tmp_fft_inp_ptr
    type(C_PTR) :: tmp_fft_out_ptr
    type(C_PTR) :: tmp_fft_tmp_ptr

    call system_clock(time_0)

    argc = command_argument_count()
    call get_command_argument(0, argv0, arg_status)
    if (argc .gt. 0) then
        allocate(argv(argc))
        do iarg = 1, argc
            call get_command_argument(iarg, argv(iarg), arg_status)
        end do
        if (trim(argv(1)) .eq. "--help") then
            write(*, "(a)") "Usage: "
            write(*, "(4x,a,1x,a)") trim(argv0), "[--help]"
            write(*, "(a)") "    or"
            write(*, "(4x,a,1x,a)") trim(argv0), "[--nsize nsize] [--ndim ndim] " // &
                                                 "[--stepsize stepsize] [--nsimu nsimu]"
            write(*, "()")
            write(*, "(a)")       "--help    : print this help message and exit"
            write(*, "(a,1x,i8)") "--nsize   : steps recorded by the simulation, default to", nsize
            write(*, "(a,1x,i8)") "--ndim    : amount of dimensions,             default to", ndim
            write(*, "(a,1x,i8)") "--stepsize: steps between each recored steps, default to", stepsize
            write(*, "(a,1x,i8)") "--nsimu   : amount of simulations,            default to", nsimu
            stop "Help message printed"
        end if
        iarg = 1
        do while(iarg .le. argc)
            if (trim(argv(iarg)) .eq. "--nsize") then
                iarg = iarg + 1
                read(argv(iarg), *) nsize
            else if (trim(argv(iarg)) .eq. "--ndim") then
                iarg = iarg + 1
                read(argv(iarg), *) ndim
            else if (trim(argv(iarg)) .eq. "--stepsize") then
                iarg = iarg + 1
                read(argv(iarg), *) stepsize
            else if (trim(argv(iarg)) .eq. "--nsimu") then
                iarg = iarg + 1
                read(argv(iarg), *) nsimu
            else
                write(*, "(a,1x,a)") "Error: unknown argument:", argv(iarg)
                stop "Unknown option"
            end if
            iarg = iarg + 1
        end do
        deallocate(argv)
    end if
    write(*, "(4(a,1x,i8))") "nsize =", nsize, ", ndim =", ndim, &
                             ", stepsize =", stepsize, ", nsimu =", nsimu

    ! call random_seed()
    call random_init_seed()

    allocate(x(0:ndim - 1, 0:nsize - 1))
    allocate(tmp_arr(0:ndim - 1))
    allocate(msd(0:nsize - 1))
    allocate(tmp_msd_cur(0:nsize - 1))

    open(write_to_unit, file = "Browning_motion_result.txt", status = "replace", action = "write")

    x = 0.0D0
    msd = 0.0D0

    write(*, "(a)") "Browning motion simulation"
    write(write_to_unit, "(a)") "# Browning motion simulation"

    write(write_to_unit, "(a,i8)") "# length of the array:   ", nsize
    write(write_to_unit, "(a,i8)") "# amount of dimensions:  ", ndim
    write(write_to_unit, "(a,i8)") "# stepsize in each step: ", stepsize
    write(write_to_unit, "(a,i8)") "# amount of simulations: ", nsimu
    write(*, "()")

    allocate(d(nsize))
    allocate(s1(nsize))
    allocate(s2(nsize))

    tmp_fft_inp_ptr = fftw_alloc_real(int(2 * nsize, C_SIZE_T))
    call c_f_pointer(tmp_fft_inp_ptr, tmp_fft_inp, (/2 * nsize/))
    tmp_fft_out_ptr = fftw_alloc_real(int(2 * nsize, C_SIZE_T))
    call c_f_pointer(tmp_fft_out_ptr, tmp_fft_out, (/2 * nsize/))
    tmp_fft_tmp_ptr = fftw_alloc_complex(int(2 * nsize, C_SIZE_T))
    call c_f_pointer(tmp_fft_tmp_ptr, tmp_fft_tmp, (/2 * nsize/))

    ! FFTW_ESTIMATE may be used to replace FFTW_MEASURE
    plan_forward  = fftw_plan_dft_r2c_1d(2 * nsize, tmp_fft_inp, tmp_fft_tmp, FFTW_MEASURE)
    plan_backward = fftw_plan_dft_c2r_1d(2 * nsize, tmp_fft_tmp, tmp_fft_out, FFTW_MEASURE)

    do isimu = 1, nsimu
        write(*, "(a,1x,i8,a,i8,a,a)", advance = "no") "simulation", isimu, "/", nsimu, "...", &
        C_CARRIAGE_RETURN
        ! this block makes a cumulatve sum up on x.
        do t = 1, nsize - 1
            x(:, t) = x(:, t - 1)
            do istep = 1, stepsize
                do idim = 0, ndim - 1
                    call std_norm_dist_useboth(tmp_arr(idim))
                end do
                x(:, t) = x(:, t) + tmp_arr
            end do
        end do

        call get_all_msd_fft(ndim, nsize, x, tmp_msd_cur, d, s1, s2, &
            tmp_fft_inp, tmp_fft_out, tmp_fft_tmp, plan_forward, plan_backward)
        ! call get_all_msd(ndim, nsize, x, tmp_msd_cur)
        msd = msd + tmp_msd_cur
    end do

    call fftw_destroy_plan(plan_forward)
    call fftw_destroy_plan(plan_backward)

    deallocate(d)
    deallocate(s1)
    deallocate(s2)
    call fftw_free(tmp_fft_inp_ptr)
    call fftw_free(tmp_fft_out_ptr)
    call fftw_free(tmp_fft_tmp_ptr)

    write(*, "(a,a)", advance = "no") repeat(" ", 32), C_CARRIAGE_RETURN
    msd = msd / dble(nsimu)
    write(write_to_unit, "(a)") "#"
    write(write_to_unit, "(a)") "#  dt         MSD(dt):  "
    do dt = 0, nsize - 1
        write(write_to_unit, "(i8,4x,f18.7)") dt, msd(dt)
    end do

    ! statistic of MSD
    call linear_through_origin_fit(nsize, msd, slope, cod)
    co_diffu = slope / (2 * ndim)
    write(*, "(a)") "statistic is from 10% to 90%."
    write(write_to_unit, "(a)") "# statistic is from 10% to 90%."
    write(*, "(a,1x,e11.4e3)") "slope:", slope
    write(*, "(a,1x,f8.4,4x,a,1x,f7.4)") "diffusivity =", co_diffu, "R**2 =", cod
    write(write_to_unit, "(a,1x,f7.4,4x,a,1x,f8.4)") "# diffusivity =", co_diffu, "R**2 =", cod

    deallocate(tmp_arr)
    deallocate(tmp_msd_cur)
    deallocate(x)
    deallocate(msd)

    close(write_to_unit)

    write(*, "(a,' ""',a,'"".')") "result wrote to", write_to_name

    call system_clock(time_t)
    write(*, "(a,1x,f6.1,1x,a)") "time elapsed:", dble(time_t - time_0) / 1.0D3, "s"

    stop
end program main

