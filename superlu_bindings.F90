module superlu_bindings
  use iso_c_binding
  implicit none

  interface
    subroutine f_pdgsmv(n, A_handle, grid_handle, x, y) bind(C, name="f_pdgsmv")
        use iso_c_binding
        implicit none
        integer(c_int64_t), value :: n
        integer(c_int64_t), value :: A_handle
        integer(c_int64_t), value :: grid_handle
        real(c_double) :: x(*)
        real(c_double) :: y(*)
    end subroutine f_pdgsmv
  end interface

end module superlu_bindings
