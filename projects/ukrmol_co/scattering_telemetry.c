/* Read-only EPSSolve interposer for the experimental native scattering audit.
 * Build against the headers/libraries of the exact image being measured.
 * It records the actual operator and residuals; it never sets solver options.
 */
#define _GNU_SOURCE
#include <dlfcn.h>
#include <stdio.h>
#include <stdlib.h>
#include <slepceps.h>

PetscErrorCode EPSSolve(EPS eps)
{
    PetscErrorCode (*solve)(EPS) = dlsym(RTLD_NEXT, "EPSSolve");
    if (!solve) abort();
    PetscErrorCode result = solve(eps);
    if (result) return result;

    Mat a;
    MatInfo info;
    MatType matrix_type;
    EPSType solver_type;
    EPSWhich which;
    EPSConvergedReason reason;
    PetscInt n, m, requested, ncv, mpd, converged, iterations;
    int rank;
    PetscErrorCode error;
#define CHECK(call) do { error = (call); if (error) return error; } while (0)
    CHECK(EPSGetOperators(eps, &a, NULL));
    CHECK(MatGetType(a, &matrix_type));
    CHECK(MatGetSize(a, &n, &m));
    CHECK(MatGetInfo(a, MAT_GLOBAL_SUM, &info));
    CHECK(EPSGetType(eps, &solver_type));
    CHECK(EPSGetWhichEigenpairs(eps, &which));
    CHECK(EPSGetDimensions(eps, &requested, &ncv, &mpd));
    CHECK(EPSGetConverged(eps, &converged));
    CHECK(EPSGetIterationNumber(eps, &iterations));
    CHECK(EPSGetConvergedReason(eps, &reason));
    MPI_Comm_rank(PetscObjectComm((PetscObject)eps), &rank);
    const char *path = getenv("UKRMOL_EPS_AUDIT");
    FILE *out = rank == 0 && path ? fopen(path, "ax") : NULL;
    if (rank == 0 && !out) abort();
    if (out) fprintf(out,
        "{\"matrix_type\":\"%s\",\"solver_type\":\"%s\","
        "\"dimension\":%lld,\"requested\":%lld,\"ncv\":%lld,\"mpd\":%lld,"
        "\"converged\":%lld,\"iterations\":%lld,\"reason\":%d,\"which\":%d,"
        "\"nz_used\":%.17g,\"nz_allocated\":%.17g,\"mallocs\":%.17g,"
        "\"petsc_reported_memory_bytes\":%.17g,\"eigenpairs\":[",
        matrix_type, solver_type, (long long)n, (long long)requested,
        (long long)ncv, (long long)mpd, (long long)converged,
        (long long)iterations, (int)reason, (int)which,
        info.nz_used, info.nz_allocated, info.mallocs, info.memory);
    for (PetscInt i = 0; i < converged; ++i) {
        PetscScalar real, imaginary;
        PetscReal absolute;
        CHECK(EPSGetEigenpair(eps, i, &real, &imaginary, NULL, NULL));
        CHECK(EPSComputeError(eps, i, EPS_ERROR_ABSOLUTE, &absolute));
        if (out) fprintf(out,
            "%s{\"value\":%.17g,\"imaginary\":%.17g,\"absolute_residual\":%.17g}",
            i ? "," : "", (double)PetscRealPart(real),
            (double)PetscRealPart(imaginary), (double)absolute);
    }
    if (out) { fprintf(out, "]}\n"); fclose(out); }
    return result;
}
