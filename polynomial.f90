module poly 
  implicit none

  type polynomial
     real, allocatable, dimension(:) :: coeff
  end type polynomial

  interface operator(+)
     module procedure poly_plus
  end interface operator(+)

  interface operator(-)
     module procedure poly_minus
  end interface operator(-)

  interface operator(*)
     module procedure poly_mul
  end interface operator(*)
contains 
  function degree(p)
    integer :: degree
    type(polynomial) :: p
    degree = size(p%coeff)-1
  end function degree

  function poly_plus(p, q)
    type(polynomial) :: poly_plus
    type(polynomial), intent(in) :: p, q
    integer :: m, m1, m2
    m1 = degree(p)
    m2 = degree(q)
    m = max(m1, m2)
    allocate(poly_plus%coeff(0:m))
    m = min(m1, m2)
    poly_plus%coeff(:m) = p%coeff(:m) + q%coeff(:m)
    if(m1 > m) poly_plus%coeff(m+1:) = p%coeff(m+1:)
    if(m2 > m) poly_plus%coeff(m+1:) = q%coeff(m+1:)
  end function poly_plus

  function poly_minus(p, q)
    type(polynomial) :: poly_minus
    type(polynomial), intent(in) :: p, q
    integer :: m, m1, m2
    m1 = degree(p)
    m2 = degree(q)
    m = max(m1, m2)
    allocate(poly_minus%coeff(0:m))
    m = min(m1, m2)
    poly_minus%coeff(:m) = p%coeff(:m) - q%coeff(:m)
    if(m1 > m) poly_minus%coeff(m+1:) = p%coeff(m+1:)
    if(m2 > m) poly_minus%coeff(m+1:) = q%coeff(m+1:)
  end function poly_minus

    
  function poly_mul(p, q)
    type(polynomial) :: poly_mul
    type(polynomial), intent(in) :: p, q
    real, dimension(:), allocatable :: pt, qt
    integer :: deg, d, d1, d2, n, k 
    d1 = degree(p)
    d2 = degree(q)
    d = max(d1, d2)
    deg = d1 + d2
    allocate(pt(0:d))
    pt = 0
    pt(0:d1) = p%coeff(0:d1)
    allocate(qt(0:d))
    qt = 0
    qt(0:d2) = q%coeff(0:d2)    
    allocate(poly_mul%coeff(0:deg))
    poly_mul%coeff = 0
    
    do n = 0, deg
       do k = 0, n
          poly_mul%coeff(n) =     poly_mul%coeff(n) + pt(k)*qt(n-k)
       end do
    end do
    deallocate(pt)
    deallocate(qt)
  end function poly_mul

  subroutine poly_print(p)
    type(polynomial), intent(in) :: p
    integer :: i

    do i = size(p%coeff)-1, 0, -1
       if ( i > 0 ) then
          write(*, '(F0.2,"x^",I0," + ")', advance="no") p%coeff(i), i
       else
          write(*, '(F0.2)') p%coeff(i)
       end if
    end do
  end subroutine poly_print

  subroutine poly_div(w, p, q, r)
    type(polynomial), intent(in) :: w, p
    type(polynomial), intent(out) :: q 
    type(polynomial), intent(out) :: r
    real, dimension(:), allocatable :: wt
    integer :: i, dw, dp
    dw = degree(w)
    dp = degree(p)
    if(dw < dp) then
       allocate(q%coeff(0:0))
       q%coeff(0) = 0
       allocate(r%coeff(0:dw))
       r%coeff = w%coeff
    else
       allocate(q%coeff(0:dw-dp))
       allocate(wt(0:dw))
       allocate(r%coeff(0:dp-1))
       wt = w%coeff
       do i = dw, dp, -1
          q%coeff(i-dp) = wt(i)/p%coeff(dp)
          wt(i-dp:i) = wt(i-dp:i) -  q%coeff(i-dp) * p%coeff(0:dp)
       end do
       r%coeff(0:dp-1) = wt(0:dp-1)
       deallocate(wt)
    end if
  end subroutine poly_div
  
end module poly
