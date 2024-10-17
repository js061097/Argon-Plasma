#include "Inputs.h"
#include <iostream>
#include <cmath>

// Declare the Species enum
enum Species {
    AR = 0,         // Represents Argon (Ar)
    AR_ION = 1,     // Represents Argon Ion (Ar⁺)
    AR_META = 2,    // Represents Argon Metastable (Ar*)
    ELECTRON = 3,   // Represents Electron (e⁻)
    NUM_SPECIES     // Represents the total number of species
};

enum ReactionType {
    REACTION1 = 0,
    REACTION2 = 1,
    REACTION3 = 2,
    REACTION4 = 3,
    REACTION5 = 4,
    NUM_REACTIONS   // Represents the total number of reactions
};

// Constructor
Inputs::Inputs() : m_rtol(0.0), m_atol(1.0e-6), m_tbegin(0.0), m_tend(1.0), m_dt(1.0e-6), 
                   m_volume(0.016128), m_surface_area(0.1871), m_pressure(0.666612), 
                   m_residtime(0.0634), m_Pinp(20), m_A_coeff(0.0), m_beta_coeff(0.0), 
                   m_ea_coeff(0.0), m_Y_E_init(0.0), m_Y_Ar_init(0.0), m_Y_Arion_init(0.0), 
                   m_Y_Armeta_init(0.0), m_Te(300), m_Tgas(0.0)
{

}

Inputs::~Inputs() 
{

}

void Inputs::initialize()
{
    readInputs();

    if(m_Tgas==0){
        std::cerr<<"Error: Tgas cannot be zero - Division by zero";
        exit(1);
    }

    // Set initial mass fractionss
    m_Y = std::vector<double>(NUM_SPECIES, 0.0);
    m_Y[AR] = m_Y_Ar_init;
    m_Y[AR_ION] = m_Y_Arion_init;
    m_Y[AR_META] = m_Y_Armeta_init;
    m_Y[ELECTRON] = m_Y_E_init;

    m_W = std::vector<double>(NUM_SPECIES, 0.0);
    m_W[AR] = m_W_Ar;
    m_W[AR_ION] = m_W_Arion;
    m_W[AR_META] = m_W_Armeta;
    m_W[ELECTRON] = m_W_E;

    m_X = std::vector<double>(NUM_SPECIES, 0.0);
    m_X[AR] =  m_mole_fract_eval(AR);
    m_X[AR_ION] =  m_mole_fract_eval(AR_ION);
    m_X[AR_META] =  m_mole_fract_eval(AR_META);
    m_X[ELECTRON] =  m_mole_fract_eval(ELECTRON);

    m_N = std::vector<double>(NUM_SPECIES, 0.0);
    m_N[AR] = m_numdensity_eval(AR);
    m_N[AR_ION] = m_numdensity_eval(AR_ION);
    m_N[AR_META] = m_numdensity_eval(AR_META);
    m_N[ELECTRON] = m_numdensity_eval(ELECTRON);

    m_moles_total = m_total_moles_eval();

    m_density = m_density_eval();
    m_density_e = m_density_elec_eval();

    // Conversion multiplier // Have to remove this as routines have been included
    m_conc_conv = m_pressure / m_k_B / m_Tgas / m_N_A;    // y to Concentration  //Is /N_A necessary?
    m_numdensity_conv = m_pressure / m_k_B / m_Tgas;      // y to Number Densitys
  
    // Calculate electron internal energy u_e, enthalpy h_e, and enthalpy due to ground sheath h_sh
    m_u_e = m_electron_internal_energy_eval();
    m_h_e = m_electron_enthalpy_eval();
    m_h_sh = m_electron_enthalpy_sh_eval();

    // Calculate molar surface production rate of Arion and Armeta
    m_bohm_speed_ion = m_bohm_speed_ion_eval();
    m_sdot_Arion = m_surf_ion_prod_rate_eval();
    m_sdot_Armeta = m_surf_meta_prod_rate_eval();

    // Calculate rate constants k for each reaction
    for (int i = 0; i < NUM_REACTIONS; i++) {
        m_k.push_back(m_reaction[i]->m_A_coeff * exp(-m_reaction[i]->m_ea_coeff / m_Te));
    }

    // Initialize production rates to zero
    m_wdot_Arion = m_molar_prod_rate_eval(AR_ION);
    m_wdot_Armeta = m_molar_prod_rate_eval(AR_META);
    m_q_r = 0.0;
    
    
    // Calculate production rate for Electrons (m_q_r)
    for (int i : m_reaction_num[ELECTRON]) {
        double concentration_product = 1.0;
        for (size_t j = 0; j < m_reaction[i]->m_reactants.size(); j++) {
            concentration_product *= pow(m_N[m_reaction[i]->m_reactants[j]], m_reaction[i]->m_coefficients[j]);
        }
        m_q_r += m_k[i] * concentration_product * m_delta_E[i];
    }

    //std::cout<<"\n Tot qr*deltaE: "<<m_volume * m_q_r<<"\n";
    // Coefficients for the 3rd Te PDE
    C = m_volume * 3 * m_R_u * m_N_A * m_m_e * m_density / 2 / m_W_E / m_W_Arion;
    //c1 = -m_m_dot * m_h_e; // Gets multiplied with Ye
    c2 = m_W_E * (m_h_e + m_h_sh) * m_surface_area; // Gets multiplied with m_sdot_Arion
    c3 = -(3 * m_volume * m_R_u / m_W_Ar) * (m_pressure / m_k_B / m_Tgas) * sqrt(8 * m_k_B / pi / m_m_e) * m_sigma_en; // Gets multiplied with m_density_e and Te-Tgas and sqrt(Te)
    c4 = -1; // Complete this after getting the qr value

    m_ydot = std::vector<double>(3,0.0);
    term4=m_volume * m_q_r;
    m_ydot[0] = (- m_Y[1]) / m_residtime + m_wdot_Arion * m_W[AR_ION] / m_density + (m_surface_area/m_density/m_volume) * m_sdot_Arion * m_W[AR_ION];
    m_ydot[1] = (- m_Y[2]) / m_residtime + m_wdot_Armeta * m_W[AR_META] / m_density + (m_surface_area/m_density/m_volume) * m_sdot_Armeta * m_W[AR_META];
    m_ydot[2] = - m_m_dot * m_Y[3] * m_h_e + (m_sdot_Arion * m_W_E * (m_h_e + m_h_sh) * m_surface_area) + c3 * m_density_e * (m_Te- m_Tgas) * sqrt(m_Te) - m_volume * m_q_r + m_Pinp; 

    //std::cout<<std::endl<<m_m_e<<std::endl;
    //std::cout<<std::endl<<3*m_volume*(m_R_u/m_W_Ar)*m_density_e*(m_Te - m_Tgas)*(m_pressure/m_k_B/m_Tgas)*sqrt(8*m_k_B*m_Te/pi/m_m_e)*m_sigma_en<<std::endl;
}

void Inputs::readInputs()
{

    // Solver related inputs
    m_rtol = 1.e-5;
    m_atol = 1.0e-10;
    m_Neq = 3;

    // Reactor parameters
    m_volume = 0.016128;
    m_surface_area = 0.1871;
    m_D_eff = 1.04;
    m_delta = 0.04056;

    // Reactor operating conditions
    m_pressure = 0.666612;
    m_residtime = 0.0634;
    m_m_dot = 2.72e-6;
    m_Pinp = 0.3e3;
    m_freq = 50e3;
    m_time_period = 1/m_freq;
    m_duty_ratio = 0.1;

    // Timings
    m_tbegin = 0.0;
    m_tend = 5*m_residtime;
    m_nout = 1e4;
    m_dt = 1e-6;
    //m_dt = m_tend / m_nout;
    m_t1 = 1e-6;

    // Production rates
    m_wdot_Arion = 0.0;
    m_wdot_Armeta = 0.0;
    m_q_r = 0.0;

    // Constants
    m_eta = 0.5;
    m_R_u = 8.314;
    m_W_Ar = 0.039948;
    m_W_Arion = 0.0399474515;
    m_W_Armeta = 0.0399480206;
    m_W_E = 5.4857990962e-7;
    m_e = 1.602176634e-19;
    m_m_e = 9.109e-31;
    m_k_B = 1.380649e-23;
    m_N_A = 6.023e23;

    m_sigma_en = 1e-19;

    // Initial conditions
    
    m_Y_Arion_init = 1e-6;
    //m_Y_E_init = 1e-6;
    m_Y_E_init = (m_W_E/m_W_Arion)*m_Y_Arion_init;
    m_Y_Armeta_init = 1e-6;
    m_Y_Ar_init = 1.0 - m_Y_Arion_init - m_Y_E_init - m_Y_Armeta_init;

    m_Te = 30000; 
    m_Tgas = 300;

    // Setting up the reactions
    m_A_coeff = {3.71e-14, 1.23e-13, 2.05e-13, 6.2e-16, 2.0e-13}; // In SI units - m3/s
    m_ea_coeff = {1.748e5, 2.168e5, 0.575e5, 0.0, 0.0};
    m_delta_E = {1.8496e-18, 2.512e-18, 6.656e-19, 0.0, -1.8496e-18};

    m_reaction.push_back(new Reaction(3.71e-14, 1.748e5, 1.8496e-18, {ELECTRON, AR}, {1, 1}));  // Reaction 1: e + Ar → e + Ar*
    m_reaction.push_back(new Reaction(1.23e-13, 2.168e5, 2.512e-18, {ELECTRON, AR}, {1, 1}));  // Reaction 2: e + Ar → 2e + Ar⁺
    m_reaction.push_back(new Reaction(2.05e-13, 0.575e5, 6.656e-19, {ELECTRON, AR_META}, {1, 1}));  // Reaction 3: e + Ar* → 2e + Ar⁺
    m_reaction.push_back(new Reaction(6.2e-16, 0.0, 0.0, {AR_META, AR_META}, {1, 1}));  // Reaction 4: Ar* + Ar* → e + Ar⁺ + Ar
    m_reaction.push_back(new Reaction(2.0e-13, 0.0, -1.8496e-18, {ELECTRON, AR_META}, {1, 1}));  // Reaction 5: e + Ar* → e + Ar

    // Map reactions to products
    
    m_reaction_num[AR_ION] = {1, 2, 3};  // Reactions producing AR_ION (Reactions 2, 4)
    m_reaction_num[AR_META] = {0};  // Reactions producing AR_META (Reaction 1)
    m_reaction_num[AR] = {3, 4};  // Reactions producing AR (from AR_META, Reaction 5)
    m_reaction_num[ELECTRON] = {0, 1, 2, 3, 4};  // Reactions producing ELECTRON (Reactions 1, 2, 3, 4)
    

   //   Using only Rxn2 - electron ionization reaction
    
   
   m_reaction_num[AR_ION] = {1};
   m_reaction_num[AR_META] = {};
   m_reaction_num[AR] = {};
   m_reaction_num[ELECTRON] = {1};
   


}

// Conversion functions
double Inputs::m_massfract_to_concentration(double massFraction)
{
    return massFraction * m_pressure / m_k_B / m_Tgas;  
}

double Inputs::m_massfract_to_numdensity(double massFraction)
{
    return massFraction * m_pressure / m_k_B / m_Tgas;  
}

double Inputs::m_molefract_to_numdensity(double moleFraction)
{
    return moleFraction * m_pressure / m_k_B / m_Tgas;  
}

double Inputs::m_numdensity_eval(int specie)
{
    return m_X[specie] * m_pressure / m_k_B / m_Tgas;
}

double Inputs::m_density_eval()
{
    return m_pressure * m_W[AR] / m_R_u / m_Tgas;
}

double Inputs::m_density_elec_eval()
{
    return m_N_A * m_m_e * m_density * m_Y[AR_ION] / m_W[AR_ION];
}

double Inputs::m_electron_internal_energy_eval()
{
    return 3 * m_R_u * m_Te / 2 / m_W[ELECTRON];
}

double Inputs::m_electron_enthalpy_eval()
{
    return 5 * m_R_u * m_Te / 2 / m_W[ELECTRON];
}

double Inputs::m_electron_enthalpy_sh_eval()
{
    return m_k_B / 2 / m_m_e / log(m_W[AR_ION] / 2 / pi / m_W[ELECTRON]);
}

double Inputs::m_electron_temperature_eval()
{
    return 2 * m_W[ELECTRON] * m_u_e / 3 / m_R_u;
}

double Inputs::m_bohm_speed_ion_eval()
{
    return sqrt(m_R_u * m_Te / m_W[AR_ION]);
}

double Inputs::m_surf_ion_prod_rate_eval()
{
    return -m_eta * m_bohm_speed_ion * (m_density / m_W[AR_ION]) * m_Y[AR_ION];
}

double Inputs::m_surf_meta_prod_rate_eval()
{
    return - (m_D_eff * m_density / m_delta / m_W[AR_META]) * m_Y[AR_META];
}

double Inputs::m_total_moles_eval()
{
    double total_moles = 0.0;
    for(int i=0;i<m_Y.size();i++)
    {
        total_moles+= m_Y[i]/m_W[i];
    }
    return total_moles;
}

double Inputs::m_mole_fract_eval(int specie)
{   
    return m_Y[specie]/m_W[specie]/m_total_moles_eval();;
}

double Inputs::m_molar_prod_rate_eval(int specie)
{
    double res = 0.0;
    for (int i : m_reaction_num[specie]) {
        double concentration_product = 1.0;
        for (size_t j = 0; j < m_reaction[i]->m_reactants.size(); j++) {
            concentration_product *= pow(m_N[m_reaction[i]->m_reactants[j]], m_reaction[i]->m_coefficients[j]);
        }
        res += m_k[i] * concentration_product;
    }
    return res / m_N_A;
}

void Inputs::finalize()
{
    // Output saying inputs are finalized
    // Display the given inputs
    std::cout << "Initialized input variables " << std::endl;
    std::cout << "Reactor geometry and Problem parameters: " << std::endl;
    std::cout << "Volume = " << m_volume << " (m3)" << std::endl;
    std::cout << "Surface Area = " << m_surface_area << " (m2)" << std::endl;
    std::cout << "Pressure = " << m_pressure << " (Pa)" << std::endl;
    std::cout << "Gas Temperature = " << m_Tgas << " (K)" << std::endl;
    std::cout << "Inlet Mass Flow Rate = " << m_m_dot << " (kg/s)" << std::endl;
    std::cout << "Flow residence time = " << m_residtime << " (s)" << std::endl;
    std::cout << "Effective diffusion coefficient = " << m_D_eff << " (m2/s)" << std::endl;
    std::cout << "Effective diffusion length scale = " << m_delta << " (m)" << std::endl;

    std::cout << "Gas chemistry variables and constants: " << std::endl;
    //  Constants
    std::cout << "Molar mass of Argon = " << m_W_Ar << " (kg)" << std::endl;
    std::cout << "Molar mass of Argon ion = " << m_W_Arion << " (kg)" << std::endl;
    std::cout << "Molar mass of Argon meta = " << m_W_Armeta << " (kg)" << std::endl;
    std::cout << "Molar mass of Electron = " << m_W_E << " (kg)" << std::endl;

    //  Evaluated variables
    std::cout << "Mole fraction of Argon = " << m_X[AR] << std::endl;
    std::cout << "Mole fraction of Argon ion = " << m_X[AR_ION] << std::endl;
    std::cout << "Mole fraction of Argon meta = " << m_X[AR_META] << std::endl;
    std::cout << "Mole fraction of Electron = " << m_X[ELECTRON] << std::endl;

    std::cout << "Molar surface production rate of Argon ion = " << m_sdot_Arion << " (mol/m2)" << std::endl;
    std::cout << "Molar surface production rate of Argon meta = " << m_sdot_Armeta << " (mol/m2)" << std::endl;
    std::cout << "Molar gas phase production rate of Argon ion = " << m_wdot_Arion << " (mol/m3)" << std::endl;
    std::cout << "Molar gas phase production rate of Argon meta = " << m_wdot_Armeta << " (mol/m3)" << std::endl;
    std::cout << "Molar electron energy lost due to collision = " << m_q_r << " (J/s)" << std::endl;
    std::cout << "Mass density = " << m_density << " (kg/m3)" << std::endl;
    std::cout << "Electron mass density = " << m_density_e << " (kg)" << std::endl;
    std::cout << "Bohm speed of Argon ion = " << m_bohm_speed_ion << " (m/s)" << std::endl;
}