#ifndef SOLVER_H
#define SOLVER_H

#include "Inputs.h"
#include "Processor.h"

#include <iostream>

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

/* Linear Solver Option */
#define DENSE_DQ      2

class Model; //Forward declaration - error resulting in not declaring prior to usage during make

class Solver 
{
public:
    Solver();
    ~Solver();
    void initialize(Model* m_model);
    int solve();
    void finalize();
    int PrepareRun(sunindextype mu,
                          sunindextype ml, SUNLinearSolver* LS,
                          SUNNonlinearSolver* NLS);

    Model* m_model;
    SUNContext m_sunctx;
    int m_retval, m_temp_retval, m_nerr, m_iout;
    N_Vector m_y;
    SUNMatrix* m_A;
    SUNLinearSolver m_LS;
    SUNNonlinearSolver m_NLS;
    void* m_cvode_mem;
    sunbooleantype m_firstrun;
    int m_qu;
    sunrealtype m_hu;
    sunrealtype m_tout,m_t;
    sunrealtype m_ero, m_er;
private:
};



#endif // SOLVER_H

