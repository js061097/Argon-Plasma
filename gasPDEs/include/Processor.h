#ifndef PROCESSOR_H
#define PROCESSOR_H

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

class Processor 
{
public:
    Processor();
    Processor(Inputs* inputs);
    ~Processor();
    void initialize();
    void printIntro();
    void printHeader();

    void initializeCSV(const std::string& filename);

    void printOutput(sunrealtype t, double y0, double y1, double y2, double y3, double n1, double n2, double temp, double ydot1, double ydot2, double ydot3, double term4, int qu, sunrealtype hu);
    void appendToCSV(const std::string& filename, double t, double y0, double y1, 
                            double y2, double y3, double n1, double n2, double temp, 
                            double ydot1, double ydot2, double ydot3, double term4, 
                            int qu, sunrealtype hu);
    void createPythonScript();
    void printFinalStats(void* cvode_mem, sunrealtype ero);
    void printErrInfo(int nerr);
    void printErrOutput(sunrealtype tol_factor);

    Inputs* m_inputs;
    bool m_csv_initialized;
};

sunrealtype MaxError(N_Vector y, sunrealtype t);

extern int checkReturnValue(void* returnvalue, const char* funcname, int opt);

#endif // POSTPROCESSOR_H
