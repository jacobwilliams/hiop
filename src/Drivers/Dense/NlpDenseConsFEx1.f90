! =============================================================================
!
!     This example is modified from NlpSparseEx1
!
!
!     min   sum 1/4* { (x_{i}-1)^4 : i=1,...,n}
!     s.t.     4*x_1 + 2*x_2                   == 10
!          5<= 2*x_1         + x_3
!          1<= 2*x_1               + 0.5*x_i   <= 2*n, for i=4,...,n
!          x_1 free
!          0.0 <= x_2
!          1.5 <= x_3 <= 10
!          x_i >=0.5, i=4,...,n
!
! =============================================================================

program example1

   use iso_c_binding, only: ip => c_int, rp => c_double, c_ptr, c_associated
   use iso_fortran_env, only: error_unit, output_unit

   implicit none

   ! =============================================================================
   ! The following should be moved into an hiop_fortran_interface.f90 module
   ! =============================================================================
   abstract interface
      ! Function pointer types for the Fortran callback functions

      subroutine f_eval_f_cb(n, x, new_x, obj) bind(c)
         import
         implicit none
         integer(ip) :: n
         real(rp), dimension(n) :: x
         integer(ip) :: new_x
         real(rp) :: obj
      end subroutine f_eval_f_cb

      subroutine f_eval_grad_cb(n, x, new_x, gradf) bind(c)
         import
         implicit none
         integer(ip) :: n
         real(rp), dimension(n) :: x
         integer(ip) :: new_x
         real(rp), dimension(n) :: gradf
      end subroutine f_eval_grad_cb

      subroutine f_eval_c_cb(n, m, x, new_x, cons) bind(c)
        import
        implicit none
        integer(ip) :: n
        integer(ip) :: m
        real(rp), dimension(n) :: x
        integer(ip) :: new_x
        real(rp), dimension(m) :: cons
    end subroutine f_eval_c_cb

    subroutine f_eval_jac_dense_cb(n, m, x, new_x, mjac) bind(c)
        import
        implicit none
        integer(ip) :: n
        integer(ip) :: m
        real(rp), dimension(n) :: x
        integer(ip) :: new_x
        real(rp), dimension(n*m) :: mjac   ! JW: this should really be an (M,N) matrix instead ...
    end subroutine f_eval_jac_dense_cb

   end interface

   interface
      ! The main procedures in the HIOP library.

      type(c_ptr) function hiopdenseprob(n, m, xlow, xupp, clow, cupp, x0, &
                                         f_eval_f, f_eval_c, f_eval_grad, f_eval_jac) bind(c)
         import
         implicit none
         integer(ip) :: n
         integer(ip) :: m
         real(rp), dimension(*) :: xlow
         real(rp), dimension(*) :: xupp
         real(rp), dimension(*) :: clow
         real(rp), dimension(*) :: cupp
         real(rp), dimension(*) :: x0
         procedure(f_eval_f_cb) :: f_eval_f
         procedure(f_eval_c_cb) :: f_eval_c
         procedure(f_eval_grad_cb) :: f_eval_grad
         procedure(f_eval_jac_dense_cb) :: f_eval_jac
      end function hiopdenseprob

      subroutine hiopdensesolve(f_prob_in, obj, sol) bind(c)
         import
         implicit none
         type(c_ptr) :: f_prob_in
         real(rp) :: obj
         real(rp), dimension(*) :: sol
      end subroutine hiopdensesolve

      subroutine deletehiopdenseprob(f_prob_in) bind(c)
         import
         implicit none
         type(c_ptr) :: f_prob_in
      end subroutine deletehiopdenseprob

   end interface
   ! =============================================================================

   integer, parameter :: n = 50
   integer, parameter :: m = 49
   real(rp), parameter :: obj_saved = 1.10351564683176e-1_rp

   real(rp) :: x(n)
   real(rp) :: x_l(n), x_u(n) 
   real(rp) :: g_l(m), g_u(m)
   real(rp) :: oobj
   type(c_ptr) :: hiopproblem
   integer(ip) :: i

   ! Set initial point and bounds:
   do i = 1, n
      x(i) = 0.0_rp
   end do

   x_l(1) = -1.0e20_rp
   x_u(1) = 1.0e20_rp
   x_l(2) = 0.0_rp
   x_u(2) = 1.0e20_rp
   x_l(3) = 1.5_rp
   x_u(3) = 10.0_rp

   do i = 4, n
      x_l(i) = 0.5_rp
      x_u(i) = 1.0e20_rp
   end do

   ! Set bounds for the constraints
   g_l(1) = 10.0_rp
   g_u(1) = 10.0_rp
   g_l(2) = 5.0_rp
   g_u(2) = 1.0e20_rp

   do i = 3, m
      g_l(i) = 1.0_rp
      g_u(i) = 2.0_rp*n
   end do

   ! create hiop sparse problem
   hiopproblem = hiopdenseprob(n, m, x_l, x_u, g_l, g_u, x, &
                               eval_f, eval_con, eval_grad, eval_jac)
   if (.not. c_associated(hiopproblem)) then
      error stop 'Error creating an HIOP Problem handle.'
   end if

   ! hiop solve
   call hiopdensesolve(hiopproblem, oobj, x)

   write (output_unit, *) ''
   write (output_unit, *) 'The optimal solution is:'
   write (output_unit, *) ''
   write (output_unit, *) 'Optimal Objective = ', oobj
   write (output_unit, *) ''

   if (abs(oobj - obj_saved) > 1e-6_rp) then
      write (error_unit, *) 'Obj mismatches SparseEx1 with 50 variables.'
      write (error_unit, *) 'Saved Obj = ', obj_saved
      error stop -1
   end if

   ! Clean up
   call deletehiopdenseprob(hiopproblem)

contains

! =============================================================================
!                    Computation of objective function
! =============================================================================
   subroutine eval_f(n, x, new_x, obj) bind(c)

      integer(ip) :: n
      real(rp), dimension(n) :: x
      integer(ip) :: new_x
      real(rp) :: obj

      obj = 0.25_rp*sum((x - 1.0_rp)**4)

   end subroutine eval_f

! =============================================================================
!                Computation of gradient of objective function
! =============================================================================
   subroutine eval_grad(n, x, new_x, grad) bind(c)

      integer(ip) :: n
      integer(ip) :: new_x
      real(rp), dimension(n) :: grad
      real(rp), dimension(n) :: x

      grad = (x - 1.0_rp)**3

   end subroutine eval_grad

! =============================================================================
!                     Computation of == equality constraints
! =============================================================================
   subroutine eval_con(n, m, x, new_x, c) bind(c)

      integer(ip) :: n
      integer(ip) :: m
      real(rp), dimension(n) :: x
      integer(ip) :: new_x
      real(rp), dimension(m) :: c

      integer(ip) :: i

      c(1) = 4.0_rp*x(1) + 2.0_rp*x(2)
      c(2) = 2.0_rp*x(1) + 1.0_rp*x(3)
      do i = 3, m
         c(i) = 2.0_rp*x(1) + 0.5_rp*x(i + 1)
      end do

   end subroutine eval_con

! =============================================================================
!                Computation of Jacobian of == equality constraints
! =============================================================================
   subroutine eval_jac(n, m, x, new_x, a) bind(c)

      integer(ip) :: n
      integer(ip) :: m
      real(rp), dimension(n) :: x
      integer(ip) :: new_x
      real(rp), dimension(n*m) :: a

      a(i:n*m) = 0.0_rp

      ! constraint 1 body ---> 4*x_1 + 2*x_2 == 10
      ! constraint 2 body ---> 2*x_1 + x_3
      ! constraint 3 body ---> 2*x_1 + 0.5*x_i, for i>=4
      a(1) = 4.0_rp
      a(2) = 2.0_rp
      a(n + 1) = 2.0_rp
      a(n + 3) = 1.0_rp

      do i = 3, m
         a((i - 1)*n + 1) = 2.0_rp
         a((i - 1)*n + i + 1) = 0.5_rp
      end do

   end subroutine eval_jac

end program example1

