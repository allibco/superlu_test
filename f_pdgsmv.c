#include "superlu_ddefs.h"

void f_pdgsmv(int_t *n, void *A, void *grid, double *x, double *y)
{
    pdgsmv_comm_t gsmv_comm;

    dInit_GSMV_comm((SuperMatrix *)A, (gridinfo_t *)grid, &gsmv_comm);

    pdgsmv(*n, (SuperMatrix *)A, (gridinfo_t *)grid, &gsmv_comm, x, y);

    dDestroy_GSMV_comm(&gsmv_comm);
}
