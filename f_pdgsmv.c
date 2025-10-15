#include "superlu_ddefs.h"
#include <stdint.h>  // for int64_t

// Fortran-callable wrapper using integer handles
void f_pdgsmv(int64_t n, int64_t A_handle, int64_t grid_handle,
                     double *x, double *y)
{
    SuperMatrix *A = (SuperMatrix *)A_handle;
    gridinfo_t *grid = (gridinfo_t *)grid_handle;

    // Call public pdgsmv routine; NULL uses default internal comm
    pdgsmv(n, A, grid, NULL, x, y);
}

