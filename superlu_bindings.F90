module superlu_bindings

  interface
    subroutine f_pdgsmv(n, A_handle, grid_handle, x, y) bind(C, name="f_pdgsmv")
        use iso_c_binding
        implicit none
        integer(c_int64_t), value :: n
        integer(superlu_ptr), value :: A_handle
        integer(superlu_ptr), value :: grid_handle
        real(c_double) :: x(*)
        real(c_double) :: y(*)
    end subroutine f_pdgsmv
end interface

endmodule superlu_bindings

