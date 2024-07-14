#ifndef POSTPROCESSOR_H
#define POSTPROCESSOR_H

#include <sundials/sundials_types.h>
#include <cvode/cvode.h>
#include "Inputs.h"
#include "Solver.h"
#include "Model.h"


#include <cvode/cvode.h>      /* prototypes for CVODE fcts., consts.          */
#include <cvode/cvode_diag.h> /* access to CVDIAG linear solver               */
#include <math.h>
#include <nvector/nvector_serial.h> /* access to serial N_Vector                    */

#include <sundials/sundials_types.h> /* definition of sunrealtype                       */
#include <sunlinsol/sunlinsol_band.h> /* access to band SUNLinearSolver               */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver              */
#include <sunmatrix/sunmatrix_band.h> /* access to band SUNMatrix                     */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix                    */

#include "sunnonlinsol/sunnonlinsol_fixedpoint.h" /* access to the fixed point SUNNonlinearSolver */
#include "sunnonlinsol/sunnonlinsol_newton.h" /* access to the newton SUNNonlinearSolver      */

class PostProcessor 
{
public:
    PostProcessor();
    PostProcessor(Inputs* inputs);
    ~PostProcessor();
    void initialize();
    void printIntro();
    void printHeader();
    void printOutput(sunrealtype t, sunrealtype y0, sunrealtype y1, sunrealtype y2, int qu, sunrealtype hu);
    void printFinalStats(void* cvode_mem, sunrealtype ero);
    void printErrInfo(int nerr);
    void printErrOutput(sunrealtype tol_factor);

    Inputs* m_inputs;
};

sunrealtype MaxError(N_Vector y, sunrealtype t);

extern int checkReturnValue(void* returnvalue, const char* funcname, int opt);

#endif // POSTPROCESSOR_H
