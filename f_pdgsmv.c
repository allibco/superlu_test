#include "superlu_ddefs.h"
#include <stdint.h>  // for int64_t

void f_pdgsmv(int64_t abs_flag, int64_t A_handle, int64_t grid_handle,
              int64_t gsmv_comm_handle, double *x, double *ax)
{
    SuperMatrix *A = (SuperMatrix *)A_handle;
    gridinfo_t *grid = (gridinfo_t *)grid_handle;
    pdgsmv_comm_t *gsmv_comm = (pdgsmv_comm_t *)gsmv_comm_handle;
    
    int_t abs_int = (int_t)abs_flag;
    
    pdgsmv(abs_int, A, grid, gsmv_comm, x, ax);
}

void f_pdgsmv_init(int64_t A_handle, int64_t row_to_proc_handle,
                   int64_t grid_handle, int64_t gsmv_comm_handle)
{
    SuperMatrix *A = (SuperMatrix *)A_handle;
    int_t *row_to_proc = (int_t *)row_to_proc_handle;
    gridinfo_t *grid = (gridinfo_t *)grid_handle;
    pdgsmv_comm_t *gsmv_comm = (pdgsmv_comm_t *)gsmv_comm_handle;
    
    pdgsmv_init(A, row_to_proc, grid, gsmv_comm);
}

int64_t f_pdgsmv_comm_create(void)
{
    pdgsmv_comm_t *gsmv_comm = (pdgsmv_comm_t *)malloc(sizeof(pdgsmv_comm_t));
    return (int64_t)gsmv_comm;
}

// Free pdgsmv_comm_t structure
void f_pdgsmv_comm_destroy(int64_t gsmv_comm_handle)
{
    pdgsmv_comm_t *gsmv_comm = (pdgsmv_comm_t *)gsmv_comm_handle;
    if (gsmv_comm != NULL) {
        free(gsmv_comm);
    }
}
