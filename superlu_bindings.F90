module superlu_bindings
  use iso_c_binding
  implicit none

  interface


   function f_pdgsmv_comm_create() bind(C, name="f_pdgsmv_comm_create")
        use iso_c_binding
        implicit none
        integer(c_int64_t) :: f_pdgsmv_comm_create
    end function f_pdgsmv_comm_create
    
    subroutine f_pdgsmv_comm_destroy(gsmv_comm_handle) &
        bind(C, name="f_pdgsmv_comm_destroy")
        use iso_c_binding
        implicit none
        integer(c_int64_t), value :: gsmv_comm_handle
    end subroutine f_pdgsmv_comm_destroy
    
    subroutine f_pdgsmv_init(A_handle, row_to_proc_handle, grid_handle, &
                             gsmv_comm_handle) bind(C, name="f_pdgsmv_init")
        use iso_c_binding
        implicit none
        integer(c_int64_t), value :: A_handle
        integer(c_int64_t), value :: row_to_proc_handle
        integer(c_int64_t), value :: grid_handle
        integer(c_int64_t), value :: gsmv_comm_handle
    end subroutine f_pdgsmv_init

    
    subroutine f_pdgsmv(abs_flag, A_handle, grid_handle, gsmv_comm_handle, x, ax) &
        bind(C, name="f_pdgsmv")
        use iso_c_binding
        implicit none
        integer(c_int64_t), value :: abs_flag
        integer(c_int64_t), value :: A_handle
        integer(c_int64_t), value :: grid_handle
        integer(c_int64_t), value :: gsmv_comm_handle
        real(c_double) :: x(*)
        real(c_double) :: ax(*)
    end subroutine f_pdgsmv
  end interface
end module superlu_bindings

