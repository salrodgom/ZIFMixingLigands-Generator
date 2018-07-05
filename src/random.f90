module mod_random
! module for pseudo random numbers
 implicit none
 private
 public init_random_seed, randint, r4_uniform
 contains
 subroutine init_random_seed( )
  use iso_fortran_env, only: int64
  implicit none
  integer, allocatable :: seed(:)
  integer              :: i, n, un, istat, dt(8), pid
  integer(int64)       :: t
  call random_seed(size = n)
  allocate(seed(n))
  ! First try if the OS provides a random number generator
  open(newunit=un, file="/dev/urandom", access="stream", &
       form="unformatted", action="read", status="old", iostat=istat)
  if (istat == 0) then
     read(un) seed
     close(un)
  else
     ! Fallback to XOR:ing the current time and pid. The PID is
     ! useful in case one launches multiple instances of the same
     ! program in parallel.
     call system_clock(t)
     if (t == 0) then
        call date_and_time(values=dt)
        t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
             + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
             + dt(3) * 24_int64 * 60 * 60 * 1000 &
             + dt(5) * 60 * 60 * 1000 &
             + dt(6) * 60 * 1000 + dt(7) * 1000 &
             + dt(8)
     end if
     pid = getpid()
     t = ieor(t, int(pid, kind(t)))
     do i = 1, n
        seed(i) = lcg(t)
     end do
  end if
  call random_seed(put=seed)
 contains
  ! This simple PRNG might not be good enough for real work, but is
  ! sufficient for seeding a better PRNG.
  function lcg(s)
    integer :: lcg
    integer(int64) :: s
    if (s == 0) then
       s = 104729
    else
       s = mod(s, 4294967296_int64)
    end if
    s = mod(s * 279470273_int64, 4294967291_int64)
    lcg = int(mod(s, int(huge(0), int64)), kind(0))
  end function lcg
 end subroutine init_random_seed
 integer function randint(i,j)
  integer,intent(in)    :: i,j
  real                  :: r
  call random_number(r)
  randint = i + floor((j+1-i)*r)
 end function randint
! 
 real function r4_uniform(x,y)
  implicit none
  real,intent(in)       :: x,y
  real                  :: r
  call random_number(r)
  r4_uniform=(r*(y-x))+x
  return
 end function r4_uniform
end module mod_random

program cosa
 use mod_random
 implicit none
 integer :: i
 real    :: r1, r2
 call init_random_seed()
 r1 = r4_uniform(0.0,1.0)
 r2 = r4_uniform(0.0,1.0)
 write(6,*) r1,r2
 do i=1,1000000
  r2=r1
  r1 = r4_uniform(0.0,1.0)
  write(6,*) r1,r2
 end do
end program cosa
