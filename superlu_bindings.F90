module superlu_bindings

  interface
     subroutine f_pdgsmv(options, A, x, y, grid) bind(c, name='f_pdgsmv')
       use iso_c_binding
       type(c_ptr), value :: options, A, grid
       real(c_double)     :: x(*), y(*)
     end subroutine f_pdgsmv
  end interface

endmodule superlu_bindings
