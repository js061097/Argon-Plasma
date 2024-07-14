#include "Model.h"
#include "Solver.h"
#include <iostream>

Model::Model(Inputs* inputs, Solver* solver)
    : m_inputs(inputs), m_solver(solver)
{

}

Model::~Model()
{
    std::cout << "Model Destructor" << std::endl;
}

void Model::initialize()
{
    std::cout << "Model Initialize" << std::endl;
    m_inputs->initialize();
    std::cout << "Inputs Initialized" << std::endl;
    m_solver->initialize(this);
    std::cout << "Solver Initialized" << std::endl;
}

int Model::solve()
{
    std::cout << "Model Solve" << std::endl;
    return m_solver->solve();
}

void Model::finalize()
{
    std::cout << "Model Finalize" << std::endl;
}

int Model::f1(double t, N_Vector y, N_Vector ydot, void* user_data)
{  
    sunrealtype y0, y1, y2;

    y0 = NV_Ith_S(y, 0);
    y1 = NV_Ith_S(y, 1);
    y2 = NV_Ith_S(y, 2);

    //std::cout << "In f1: t = " << t << ", y0 = " << y0 << ", y1 = " << y1 << ", y2 = " << y2 << std::endl;

    NV_Ith_S(ydot, 0) = -y0/m_inputs->m_residtime/0.0634; //m_inputs->m_residtime; - invalid use of member in static function
    NV_Ith_S(ydot, 1) = -y1/0.0634; //m_inputs->m_residtime; 
    NV_Ith_S(ydot, 2) = 0;   
    
    //std::cout << "In f1: ydot0 = " << NV_Ith_S(ydot, 0) << ", ydot1 = " << NV_Ith_S(ydot, 1) << ", ydot2 = " << NV_Ith_S(ydot, 2) << std::endl;
    
    return 0;
}
