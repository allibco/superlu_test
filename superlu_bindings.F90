module superlu_bindings

  interface
    subroutine f_pdgsmv(n, A_handle, grid_handle, x, y) bind(C, name="f_pdgsmv_handle")
        use iso_c_binding
        implicit none
        integer(c_int64_t), value :: n
        integer(c_int64_t), value :: A_handle
        integer(c_int64_t), value :: grid_handle
        real(c_double) :: x(*)
        real(c_double) :: y(*)
    end subroutine f_pdgsmv_handle
end interface

endmodule superlu_bindings

