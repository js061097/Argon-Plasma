#include "Model.h"
#include "Solver.h"
#include "Inputs.h"
#include <iostream>
#include <algorithm>

int k=0;

enum Species {
    AR = 0,         // Represents Argon (Ar)
    AR_ION = 1,     // Represents Argon Ion (Ar⁺)
    AR_META = 2,    // Represents Argon Metastable (Ar*)
    ELECTRON = 3,   // Represents Electron (e⁻)
    NUM_SPECIES     // Represents the total number of species
};

Model::Model(Inputs* inputs, Solver* solver)
    : m_inputs(inputs),m_solver(solver)
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
    //std::cout << "Inputs Initialized" << std::endl; //Made inputs as static
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
    //  For minimizing pointer to pointer calls - Think of it as typedef
    Inputs* user_inputs = ((Model*)user_data)->m_inputs;

    sunrealtype y0, y1, y2;

    //  Access each element from y vector : y0 - Mass fraction of Argon ion, y1 - Mass fraction of Argon meta stable,
    //  y2 - product of electron internal energy, electron enthalpy, and volume
    y0 = NV_Ith_S(y, 0);
    y1 = NV_Ith_S(y, 1);
    y2 = NV_Ith_S(y, 2);

    double y_e = (user_inputs->m_W_E/user_inputs->m_W[0])*y0;
    double y_Ar = 1 - y0 - y1 - y_e; 

    //  Updating mass fraction values in inputs object
    user_inputs->m_Y[AR]= y_Ar;
    user_inputs->m_Y[AR_ION] = y0;
    user_inputs->m_Y[AR_META]= y1;
    user_inputs->m_Y[ELECTRON] = y_e;

    user_inputs->m_X[AR] =  user_inputs->m_mole_fract_eval(AR);
    user_inputs->m_X[AR_ION] =  user_inputs->m_mole_fract_eval(AR_ION);
    user_inputs->m_X[AR_META] =  user_inputs->m_mole_fract_eval(AR_META);
    user_inputs->m_X[ELECTRON] =  user_inputs->m_mole_fract_eval(ELECTRON);

    user_inputs->m_N[AR] =  user_inputs->m_numdensity_eval(AR);
    user_inputs->m_N[AR_ION] =  user_inputs->m_numdensity_eval(AR_ION);
    user_inputs->m_N[AR_META] =  user_inputs->m_numdensity_eval(AR_META);
    user_inputs->m_N[ELECTRON] =  user_inputs->m_numdensity_eval(ELECTRON);

    //std::cout<<"\n"<<user_inputs->m_N[AR_ION]<<"\n";

    //  Updating electron properties in input object
    user_inputs->m_density_e = user_inputs->m_density_elec_eval();
    //user_inputs->m_u_e = user_inputs->m_electron_internal_energy_eval(); //This is u_e from prev timestep as the temperature has changed in this timestep
    user_inputs->m_u_e = y2 / user_inputs->m_density_e / user_inputs->m_volume;
    user_inputs->m_Te = user_inputs->m_electron_temperature_eval();
    
    if(user_inputs->m_Y[AR_ION]<0 || user_inputs->m_Y[AR_META]<0)
    {
        user_inputs->m_Y[AR_ION]=0;
        user_inputs->m_N[AR_ION]=0;
        user_inputs->m_Y[ELECTRON]=0;
        user_inputs->m_N[ELECTRON]=0;
        
    }

    user_inputs->m_h_e = user_inputs->m_electron_enthalpy_eval();

    //  Constant value used in the equations denoting: A_DV =  Asurf / (density * Volume)
    double A_DV = user_inputs->m_surface_area / user_inputs->m_density / user_inputs->m_volume;

    // Calculate rate constants k for each reaction
    for (int i = 0; i < 5; i++) 
    {
        user_inputs->m_k[i] = (user_inputs->m_reaction[i]->m_A_coeff * exp(-user_inputs->m_reaction[i]->m_ea_coeff / user_inputs->m_Te));
    }

    user_inputs->m_bohm_speed_ion = user_inputs->m_bohm_speed_ion_eval();   
    user_inputs->m_sdot_Arion = user_inputs->m_surf_ion_prod_rate_eval();   //Molar surf ion prod rate
    user_inputs->m_sdot_Armeta = user_inputs->m_surf_meta_prod_rate_eval(); //Molar surf meta prod rate

    // Calculate production rate for Argon Ion (Ar⁺)
    user_inputs->m_wdot_Arion = user_inputs->m_molar_prod_rate_eval(AR_ION);
    // Calculate production rate for Argon Metastable (Ar*)
    user_inputs->m_wdot_Armeta = user_inputs->m_molar_prod_rate_eval(AR_META);

     // Calculate production rate for Electrons (m_q_r)
    user_inputs->m_q_r = 0.0;
    for (int i : user_inputs->m_reaction_num[ELECTRON]) 
    {
        double concentration_product = 1.0;
        for (size_t j = 0; j < user_inputs->m_reaction[i]->m_reactants.size(); j++)
        {
            concentration_product*= pow(user_inputs->m_N[user_inputs->m_reaction[i]->m_reactants[j]], user_inputs->m_reaction[i]->m_coefficients[j]);
        }
        user_inputs->m_q_r += user_inputs->m_k[i] * concentration_product * user_inputs->m_delta_E[i];
    }

    NV_Ith_S(ydot, 0) = (- y0) / user_inputs->m_residtime + user_inputs->m_wdot_Arion * user_inputs->m_W[AR_ION] / user_inputs->m_density + A_DV * user_inputs->m_sdot_Arion * user_inputs->m_W[AR_ION];
    user_inputs->m_ydot[0]=NV_Ith_S(ydot, 0);

    NV_Ith_S(ydot, 1) = (- y1) / user_inputs->m_residtime + user_inputs->m_wdot_Armeta * user_inputs->m_W[AR_META] / user_inputs->m_density + A_DV * user_inputs->m_sdot_Armeta * user_inputs->m_W[AR_META];
    user_inputs->m_ydot[1]=NV_Ith_S(ydot, 1);

    NV_Ith_S(ydot, 2) = - user_inputs->m_m_dot * y_e * user_inputs->m_h_e + (user_inputs->m_sdot_Arion * user_inputs->m_W_E * (user_inputs->m_h_e + user_inputs->m_h_sh) * user_inputs->m_surface_area) + user_inputs->c3 * user_inputs->m_density_e * (user_inputs->m_Te - user_inputs->m_Tgas) * sqrt(user_inputs->m_Te) - user_inputs->m_volume * user_inputs->m_q_r; 
    /*
    if(t<=1e-9)
        NV_Ith_S(ydot, 2) += user_inputs->m_Pinp * (t/1e-9);
    else
        NV_Ith_S(ydot, 2) += user_inputs->m_Pinp;
    */
    
    if(std::fmod(t,user_inputs->m_time_period)<=user_inputs->m_time_period* user_inputs->m_duty_ratio)
        NV_Ith_S(ydot, 2) += user_inputs->m_Pinp * (1/user_inputs->m_duty_ratio);
    else
        NV_Ith_S(ydot, 2) += 0.1;
    

    user_inputs->term4=user_inputs->m_volume * user_inputs->m_q_r;
    user_inputs->m_ydot[2]=NV_Ith_S(ydot, 2);

    return 0;
    
}