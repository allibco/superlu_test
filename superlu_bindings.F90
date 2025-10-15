module superlu_bindings

  interface
     subroutine f_pdgsmv(n, A, grid, x, y) bind(C, name="f_pdgsmv")
       use iso_c_binding
       implicit none
       integer(c_int64_t), value :: n   ! matches C int_t
       type(c_ptr), value :: A          ! pointer to SuperMatrix
       type(c_ptr), value :: grid       ! pointer to gridinfo_t
       real(c_double) :: x(*)           ! input vector
       real(c_double) :: y(*)           ! output vector
     end subroutine f_pdgsmv
  end interface

  
endmodule superlu_bindings

