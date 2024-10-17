#include "Solver.h"
#include "Model.h"
#include <iostream>

Solver::Solver()
    : m_cvode_mem(nullptr), m_NLS(nullptr), m_LS(nullptr), m_A(nullptr), m_y(nullptr), m_sunctx(nullptr)
{

}

Solver::~Solver()
{
    std::cout << "Entering Solver Destructor" << std::endl;

    if (m_cvode_mem != nullptr) 
    {
        std::cout << "Freeing CVode memory" << std::endl;
        CVodeFree(&m_cvode_mem);
    }

    if (m_NLS != nullptr) 
    {
        std::cout << "Freeing Nonlinear Solver memory" << std::endl;
        SUNNonlinSolFree(m_NLS);
    }

    if (m_LS != nullptr) 
    {
        std::cout << "Freeing Linear Solver memory" << std::endl;
        SUNLinSolFree(m_LS);
    }

    if (m_A != nullptr) 
    {
        std::cout << "Destroying SUNMatrix" << std::endl;
        SUNMatDestroy(*m_A);
        delete m_A;
    }

    if (m_y != nullptr) 
    {
        std::cout << "Destroying N_Vector" << std::endl;
        N_VDestroy(m_y);
    }

    if (m_sunctx != nullptr) 
    {
        std::cout << "Freeing SUNContext" << std::endl;
        SUNContext_Free(&m_sunctx);
    }

    std::cout << "Exiting Solver Destructor" << std::endl;

}

void Solver::initialize(Model* model)
{
    m_model = model;

    std::cout << "Creating SUNContext" << std::endl;
    m_retval = SUNContext_Create(SUN_COMM_NULL, &m_sunctx); 
    if (checkReturnValue(&m_retval, "SUNContext_Create", 1)) { exit(1); }

    std::cout << "Creating N_Vector" << std::endl;
    m_y = N_VNew_Serial(model->m_inputs->m_Neq, m_sunctx);
    if (checkReturnValue((void*)m_y, "N_VNew_Serial", 0)) { exit(1); }
    
    std::cout << "Creating CVode" << std::endl;
    m_cvode_mem = CVodeCreate(CV_BDF, m_sunctx);
    if (checkReturnValue((void*)m_cvode_mem, "CVodeCreate", 0)) { exit(1); }
    
    NV_Ith_S(m_y, 0) = model->m_inputs->m_Y[1];
    NV_Ith_S(m_y, 1) = model->m_inputs->m_Y[2];
    NV_Ith_S(m_y, 2) = model->m_inputs->m_density_e * model->m_inputs->m_u_e * model->m_inputs->m_volume;
    
    std::cout << "Attaching model pointer to user_data" << std::endl;
    m_retval = CVodeSetUserData(m_cvode_mem, m_model);
    if (checkReturnValue(&m_retval, "SUNContext_Create", 1)) { exit(1); }

    std::cout << "Initializing CVode" << std::endl;
    m_retval = CVodeInit(m_cvode_mem, Model::f1, model->m_inputs->m_tbegin, m_y);
    if (checkReturnValue(&m_retval, "CVodeInit", 1)) { exit(1); }

    /*
    // Set the initial time step for CVode
    std::cout << "Setting initial step size" << std::endl;
    m_retval = CVodeSetInitStep(m_cvode_mem, 1e-8);  // Set small initial step size
    if (checkReturnValue(&m_retval, "CVodeSetInitStep", 1)) { exit(1); }
    
    // Set the maximum order for CVODE
    std::cout << "Setting maximum order to 1" << std::endl;
    m_retval = CVodeSetMaxOrd(m_cvode_mem, 1);
    if (checkReturnValue(&m_retval, "CVodeSetMaxOrd", 1)) { exit(1); }

    // Set the maximum number of steps for CVODE
    std::cout << "Setting the Max Number of Steps" << std::endl;
    CVodeSetMaxNumSteps(m_cvode_mem, 1e5);  // Increase maximum steps
    if (checkReturnValue(&m_retval, "CVodeSetMaxNumSteps", 1)) { exit(1); }
    */
   
    std::cout << "Setting CVode tolerances" << std::endl;
    m_retval = CVodeSStolerances(m_cvode_mem, model->m_inputs->m_rtol, model->m_inputs->m_atol);
    if (checkReturnValue(&m_retval, "CVodeSStolerances", 1)) { exit(1); }
    
}

int Solver::solve() 
{
    std::cout << "Entering solve function" << std::endl;
    m_ero = 0.0;

    std::cout << "Preparing run" << std::endl;
    m_retval = PrepareRun(0, 0, &m_LS, &m_NLS);
    if (checkReturnValue(&m_retval, "PrepareRun", 1)) { return 1; }

    Processor pp;
    pp.printHeader();

    //pp.printOutput(m_t, m_model->m_inputs->m_Y[0], m_model->m_inputs->m_Y[1], m_model->m_inputs->m_Y[2], m_model->m_inputs->m_Y[3], m_model->m_inputs->m_N[1], m_model->m_inputs->m_N[2], m_model->m_inputs->m_Te, m_model->m_inputs->m_ydot[0], m_model->m_inputs->m_ydot[1], m_model->m_inputs->m_ydot[2], m_model->m_inputs->term4, m_qu, m_hu);
        
    for (m_iout = 1, m_tout = m_model->m_inputs->m_t1; m_iout <= m_model->m_inputs->m_nout; m_iout++, m_tout += m_model->m_inputs->m_dt) 
    {
        //std::cout << "Calling CVode at time " << m_tout << std::endl;
        //std::cout << "Initial values: Yar+ = " << NV_Ith_S(m_y, 0) << ", Yar* = " << NV_Ith_S(m_y, 1) << ", Te = " << NV_Ith_S(m_y, 2) << std::endl;
        
        m_retval = CVode(m_cvode_mem, m_tout, m_y, &m_t, CV_NORMAL);
        //std::cout << "CVode return value: " << m_retval << std::endl;
        if (checkReturnValue(&m_retval, "CVode", 1)) { return 1; }
        
        int temp_retval = CVodeGetLastOrder(m_cvode_mem, &m_qu);
        if (checkReturnValue(&temp_retval, "CVodeGetLastOrder", 1)) { ++m_nerr; }
        m_temp_retval = CVodeGetLastStep(m_cvode_mem, &m_hu);
        if (checkReturnValue(&temp_retval, "CVodeGetLastStep", 1)) { ++m_nerr; } 
        
        //std::cout<<m_t<<" "<<m_model->m_inputs->m_u_e<<"\n";


        pp.printOutput(m_t, m_model->m_inputs->m_Y[0], m_model->m_inputs->m_Y[1], m_model->m_inputs->m_Y[2], m_model->m_inputs->m_Y[3], m_model->m_inputs->m_N[1], m_model->m_inputs->m_N[2], m_model->m_inputs->m_Te, m_model->m_inputs->m_ydot[0], m_model->m_inputs->m_ydot[1], m_model->m_inputs->m_ydot[2], m_model->m_inputs->term4, m_qu, m_hu);
        if (m_retval != CV_SUCCESS) {
            m_nerr++;
            break;
        }
        
        
        if (m_iout % 2 == 0) 
        {
            m_er = fabs(NV_Ith_S(m_y, 0)) / m_model->m_inputs->m_atol;
            if (m_er > m_ero) { m_ero = m_er; }
            if (m_er > 1e12) {
                m_nerr++;
                pp.printErrOutput(1e12);
                //P1_TOL_FACTOR
            }
        }
        pp.createPythonScript();
        
    }

    pp.printFinalStats(m_cvode_mem, m_ero);
    pp.printErrInfo(m_nerr);

    std::cout << "Exiting solve function" << std::endl;
    return m_nerr;
}

int Solver::PrepareRun(sunindextype mu, sunindextype ml, SUNLinearSolver* LS, SUNNonlinearSolver* NLS) 
{
    std::cout << "Entering PrepareRun" << std::endl;
    int m_retval = CV_SUCCESS;

    std::cout << "Creating new NLS" << std::endl;
    *NLS = SUNNonlinSol_Newton(m_y, m_sunctx);
    if (checkReturnValue((void*)*NLS, "SUNNonlinSol_Newton", 0)) { return 1; }

    std::cout << "Setting Nonlinear Solver" << std::endl;
    m_retval = CVodeSetNonlinearSolver(m_cvode_mem, *NLS);
    if (checkReturnValue(&m_retval, "CVodeSetNonlinearSolver", 1)) { return 1; }

    std::cout << "Creating new SUNMatrix" << std::endl;
    std::cout << "model->m_inputs->m_Neq: "  << m_model->m_inputs->m_Neq << std::endl;
    std::cout << "SUNContext: " << (void*)m_sunctx << std::endl;

    if (m_model->m_inputs->m_Neq <= 0 || m_sunctx == nullptr) 
    {
        std::cerr << "Invalid NEQ or SUNContext" << std::endl;
        return 1;
    }

    m_A = new SUNMatrix;
    if (!m_A) 
    {
        std::cerr << "Failed to allocate memory for SUNMatrix" << std::endl;
        return 1;
    }

    *m_A = SUNDenseMatrix(m_model->m_inputs->m_Neq, m_model->m_inputs->m_Neq, m_sunctx);
    if (*m_A == nullptr) 
    {
        std::cerr << "Failed to create SUNMatrix" << std::endl;
        delete m_A;
        return 1;
    }
    
    if (checkReturnValue((void*)*m_A, "SUNDenseMatrix", 0)) { return 1; }

    std::cout << "Creating new LS" << std::endl;
    *LS = SUNLinSol_Dense(m_y, *m_A, m_sunctx);
    if (checkReturnValue((void*)*LS, "SUNLinSol_Dense", 0)) { return 1; }

    std::cout << "Setting Linear Solver" << std::endl;
    m_retval = CVodeSetLinearSolver(m_cvode_mem, *LS, *m_A);
    if (checkReturnValue(&m_retval, "CVodeSetLinearSolver", 1)) { return 1; }

    std::cout << "Setting Jacobian function" << std::endl;
    m_retval = CVodeSetJacFn(m_cvode_mem, NULL);
    if (checkReturnValue(&m_retval, "CVodeSetJacFn", 5)) { return 1; }

    std::cout << "Exiting PrepareRun" << std::endl;
    return m_retval;
}
