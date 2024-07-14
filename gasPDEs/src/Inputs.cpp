#include "Inputs.h"
#include <iostream>

//  Constructor
Inputs::Inputs() : m_rtol(0.0),m_atol(1.0e-6),m_tbegin(0.0),m_tend(1.0),m_dt(1.0e-6),m_volume(0.016128),m_surface_area(0.1871),m_pressure(0.666612),m_residtime(0.0634),m_Pinp(20),m_A_coeff(0.0),m_beta_coeff(0.0),m_ea_coeff(0.0),m_N_E_init(0.0),m_N_Ar_init(0.0),m_N_Arion_init(0.0),m_N_Armeta_init(0.0),m_Te(300),m_Tgas(0.0)   
{
        
}

Inputs::~Inputs() {
}

void Inputs::initialize()
{
    readInputs();
}

void Inputs::readInputs()
{ 
        //  solver related inputs
        m_rtol = 1.e-5;
        m_atol = 1.0e-12;
        Neq = 3;
        
        //  timings
        m_tbegin = 0.0;
        m_tend = 1.0;
        m_dt = 1.0e-6;
        
        //  reactor parameters
        m_volume = 0.016128;
        m_surface_area = 0.1871;

        //  reactor operating conditions
        m_pressure = 0.666612;  // (units: Pa)
        m_residtime = 0.0634; // (units: sec)
        m_Pinp = 20;  // input electron power (units: W) //Assigned with random initial constant value

        //  gas-chemistry related
        m_A_coeff = 0.0;
        m_beta_coeff = 0.0;
        m_ea_coeff = 0.0;

        //  initial conditions
        m_N_E_init = 0.0;
        m_N_Ar_init = 0.0;
        m_N_Arion_init = 1.0;
        m_N_Armeta_init = 0.0;
        m_Te = 300; //Verify
        m_Tgas = 0.0;
}

void Inputs::finalize()
{

}
