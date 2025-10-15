#include "superlu_ddefs.h"

// Fortran-callable wrapper
void f_pdgsmv(int_t *n, void *A, void *grid, double *x, double *y)
{
    SuperMatrix *pA = (SuperMatrix *)A;
    gridinfo_t *pgrid = (gridinfo_t *)grid;

    // Direct call to pdgsmv
    // pdgsmv is in libsuperlu_dist and can be linked
    pdgsmv(*n, pA, pgrid, NULL, x, y);
}
