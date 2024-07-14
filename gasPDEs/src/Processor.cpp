#include "Processor.h"
#include <iostream>
#include "Solver.h"
#include "Inputs.h"

PostProcessor::PostProcessor()
{
   
}

PostProcessor::PostProcessor(Inputs* inputs)
{
    m_inputs = inputs;
}

PostProcessor::~PostProcessor()
{

}

void PostProcessor::initialize()
{

}

void PostProcessor::printIntro() 
{
    printf("Demonstration program for CVODE package - direct linear solvers\n");
    printf("\n\n");
    printf("Problem: Argon Plasma Physics\n");
#if defined(SUNDIALS_EXTENDED_PRECISION)
    printf(" neq = %d,  reltol = %.2Lg,  abstol = %.2Lg", m_inputs->Neq, m_inputs->m_rtol, m_inputs->m_atol);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    printf(" neq = %d,  reltol = %.2g,  abstol = %.2g", m_inputs->Neq, m_inputs->m_rtol, m_inputs->m_atol);
#else
    printf(" neq = %d,  reltol = %.2g,  abstol = %.2g", m_inputs->Neq, m_inputs->m_rtol, m_inputs->m_atol);
#endif
}

void PostProcessor::printHeader() 
{
    printf("\n     t           Yar+              Yar*              Te         qu     hu \n");
}

void PostProcessor::printOutput(sunrealtype t, sunrealtype y0, sunrealtype y1, sunrealtype y2, int qu, sunrealtype hu) 
{
#if defined(SUNDIALS_EXTENDED_PRECISION)
    printf("%10.5Lf    %12.5Le   %12.5Le   %12.5e   %2d    %6.4Le\n", t, y0, y1, y2, qu, hu);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    printf("%10.5f    %12.5e   %12.5e   %12.5e   %2d    %6.4e\n", t, y0, y1, y2, qu, hu);
#else
    printf("%10.5f    %12.5e   %12.5e   %12.5e   %2d    %6.4e\n", t, y0, y1, y2, qu, hu);
#endif
}

void PostProcessor::printErrOutput(sunrealtype tol_factor) 
{
#if defined(SUNDIALS_EXTENDED_PRECISION)
    printf("\n\n Error exceeds %Lg * tolerance \n\n", tol_factor);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    printf("\n\n Error exceeds %g * tolerance \n\n", tol_factor);
#else
    printf("\n\n Error exceeds %g * tolerance \n\n", tol_factor);
#endif
}

void PostProcessor::printFinalStats(void* cvode_mem, sunrealtype ero) 
{
    long int lenrw, leniw, lenrwLS, leniwLS;
    long int nst, nfe, nsetups, nni, ncfn, netf, nje, nfeLS;
    int retval;

    retval = CVodeGetWorkSpace(cvode_mem, &lenrw, &leniw);
    checkReturnValue(&retval, "CVodeGetWorkSpace", 1);
    retval = CVodeGetNumSteps(cvode_mem, &nst);
    checkReturnValue(&retval, "CVodeGetNumSteps", 1);
    retval = CVodeGetNumRhsEvals(cvode_mem, &nfe);
    checkReturnValue(&retval, "CVodeGetNumRhsEvals", 1);
    retval = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
    checkReturnValue(&retval, "CVodeGetNumLinSolvSetups", 1);
    retval = CVodeGetNumErrTestFails(cvode_mem, &netf);
    checkReturnValue(&retval, "CVodeGetNumErrTestFails", 1);
    retval = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
    checkReturnValue(&retval, "CVodeGetNumNonlinSolvIters", 1);
    retval = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
    checkReturnValue(&retval, "CVodeGetNumNonlinSolvConvFails", 1);

    printf("\n Final statistics for this run:\n\n");
    printf(" CVode real workspace length              = %4ld \n", lenrw);
    printf(" CVode integer workspace length           = %4ld \n", leniw);
    printf(" Number of steps                          = %4ld \n", nst);
    printf(" Number of f-s                            = %4ld \n", nfe);
    printf(" Number of setups                         = %4ld \n", nsetups);
    printf(" Number of nonlinear iterations           = %4ld \n", nni);
    printf(" Number of nonlinear convergence failures = %4ld \n", ncfn);
    printf(" Number of error test failures            = %4ld \n\n", netf);

    nje    = nsetups;
    retval = CVodeGetNumRhsEvals(cvode_mem, &nfeLS);
    checkReturnValue(&retval, "CVDiagGetNumRhsEvals", 1);
    retval = CVodeGetWorkSpace(cvode_mem, &lenrwLS, &leniwLS);
    checkReturnValue(&retval, "CVDiagGetWorkSpace", 1);
    
    printf(" Linear solver real workspace length      = %4ld \n", lenrwLS);
    printf(" Linear solver integer workspace length   = %4ld \n", leniwLS);
    printf(" Number of Jacobian evaluations           = %4ld \n", nje);
    printf(" Number of f evals. in linear solver      = %4ld \n\n", nfeLS);

#if defined(SUNDIALS_EXTENDED_PRECISION)
    printf(" Error overrun = %.3Lf \n", ero);
#else
    printf(" Error overrun = %.3f \n", ero);
#endif
}

void PostProcessor::printErrInfo(int nerr) 
{
    printf("\n\n-------------------------------------------------------------");
    printf("\n-------------------------------------------------------------");
    printf("\n\n Number of errors encountered = %d \n", nerr);
}

int checkReturnValue(void* returnvalue, const char* funcname, int opt) 
{
    int* retval;

    if (opt == 0 && returnvalue == NULL) {
        fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n", funcname);
        return 1;
    } else if (opt == 1) {
        retval = (int*)returnvalue;
        if (*retval < 0) {
            fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with retval = %d\n\n", funcname, *retval);
            return 1;
        }
    } else if (opt == 2 && returnvalue == NULL) {
        fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n", funcname);
        return 1;
    }

    return 0;
}
