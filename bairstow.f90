program bairstow
  use poly
  use blas95
  use lapack95

  type(polynomial) :: w, q, p
  allocate(w%coeff(0:5))
  w%coeff = ((/ 6, 11, -33, -33, 11, 6 /))
  call poly_print(w)
  p = bairstow_find_div(w, q, 0.0001)
  print *, "found factor"
  call poly_print(p)
  print *, "remainder"
  call poly_print(q)
  do 
     if(degree(q) < 3) exit
     deallocate(w%coeff)
     allocate(w%coeff(0:degree(q)))
     w%coeff = q%coeff
     deallocate(q%coeff)
     deallocate(p%coeff)
     p = bairstow_find_div(w, q, 0.0001)
     print *, "found factor"
     call poly_print(p)
     print *, "remainder"
     call poly_print(q)
  end do
contains
  function bairstow_iteration(w, p)
    type(polynomial) :: bairstow_iteration
    type(polynomial), intent(in) :: w, p
    type(polynomial) :: q1, q2, r1, r2
    integer, dimension(:), allocatable :: ipiv
    real :: u, v, c, d, g, h
    real, dimension(2, 2) :: a
    real, dimension(0:1) :: y, cd
    u = p%coeff(1)
    v = p%coeff(0)
    call poly_div(w, p, q1, r1)
    c = r1%coeff(1)
    d = r1%coeff(0)
    !print *, "c = ", c, "d = ", d
    call poly_div(q1, p, q2, r2)
    g = r2%coeff(1)
    h = r2%coeff(0)
    !print *, "g = ", g, "h = ", h
    a(1,1) = g*u - h
    a(1,2) = -g
    a(2,1) = g*v
    a(2,2) = -h
    allocate(ipiv(2))
    call getrf(a, ipiv) ! rozklad LU
    call getri(a, ipiv) ! odwrotnosc
    !print *, "a = ", a
    deallocate(ipiv)
    cd(0) = r1%coeff(1)
    cd(1) = r1%coeff(0)
    call gemv(a, cd, y) ! y = A^-1 [c, d]
    !print *, "y = ", y
    cd(0) = y(1)
    cd(1) = y(0)
    allocate(bairstow_iteration%coeff(0:2))
    bairstow_iteration%coeff(2) = 1
    bairstow_iteration%coeff(0:1) = p%coeff(0:1) - cd
  end function bairstow_iteration

  function bairstow_find_div(w, q, eps)
    type(polynomial) :: bairstow_find_div
    type(polynomial), intent(in) :: w
    type(polynomial), intent(out) :: q
    real, intent(in) :: eps
    type(polynomial) :: p, pp, r
    integer :: d
    d = degree(w)
    if(d < 3) then
       allocate(bairstow_find_div%coeff(0:d))
       bairstow_find_div%coeff = w%coeff
    else
       allocate(bairstow_find_div%coeff(0:2))
       allocate(p%coeff(0:2))
       p%coeff = (/ w%coeff(d-1)/w%coeff(d), w%coeff(d-1)/w%coeff(d), 1 /)
       !call poly_print(p)
       pp = bairstow_iteration(w, p)
       do
          call poly_div(w, pp, q, r)
          !call poly_print(pp)
          if(nrm2(r%coeff) < eps) exit
          p%coeff = pp%coeff
          deallocate(pp%coeff)
          deallocate(q%coeff)
          pp = bairstow_iteration(w, p)
       end do
       bairstow_find_div%coeff = pp%coeff
       deallocate(pp%coeff)
       deallocate(p%coeff)
    end if    
  end function bairstow_find_div
end program bairstow
