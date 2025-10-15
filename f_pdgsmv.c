#include "superlu_ddefs.h"

void f_pdgsmv(void *opt, void *A, double *x, double *y, void *grid)
{
    pdgsmv((superlu_dist_options_t *)opt,
           (SuperMatrix *)A, x, y,
           (gridinfo_t *)grid);
}
