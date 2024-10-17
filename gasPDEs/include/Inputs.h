#ifndef INPUTS_H
#define INPUTS_H

#include <iostream>
#include <vector>
#include <unordered_map>
#include <functional>
#include "Reaction.h"

#define pi 3.14

class Reaction;

class Inputs {
    public:
        Inputs();
        ~Inputs();
        int m_Neq;

        void readInputs();
        void initialize();
        void finalize();

        // Conversion functions
        double m_massfract_to_concentration(double massFraction);
        double m_massfract_to_numdensity(double massFraction);
        double m_molefract_to_numdensity(double moleFraction);
        double m_numdensity_eval(int specie);

        // Parameter evaluation functions
        double m_density_eval();
        double m_density_elec_eval();

        double m_electron_internal_energy_eval();
        double m_electron_enthalpy_eval();
        double m_electron_enthalpy_sh_eval();

        double m_bohm_speed_ion_eval();
        double m_surf_ion_prod_rate_eval();
        double m_surf_meta_prod_rate_eval();

        double m_total_moles_eval();
        double m_mole_fract_eval(int specie);

        double m_molar_prod_rate_eval(int specie);
        double m_electron_temperature_eval();


        // Solver-related inputs
        double m_rtol;
        double m_atol;

        // Timings
        double m_tbegin;
        double m_tend;
        double m_dt;
        double m_nout;
        double m_t1;

        // Reactor parameters
        double m_volume;
        double m_surface_area;
        double m_D_eff;   // Effective diffusion coefficient
        double m_delta;   // Effective diffusion length

        // Reactor operating conditions
        double m_pressure;  // (units: Pa)
        double m_residtime; // (units: sec)
        double m_m_dot;     // Feed gas flow rate   (units: kg/s)
        double m_Pinp;      // Input electron power (units: W)
        double m_freq;      //  Frequency of power input (units: Hz or 1/sec)
        double m_time_period;     //  Time period of power input pulse (units: sec)
        double m_duty_ratio = 0.1;

        // Gas-chemistry related - Arrhenius equation coefficients
        std::vector<double> m_A_coeff;
        double m_beta_coeff;
        std::vector<double> m_ea_coeff;
        std::vector<double> m_delta_E;  // Energy lost/electron due to inelastic collision

        std::vector<double> m_k;
        std::vector<Reaction*> m_reaction;
        std::unordered_map<int, std::vector<int>> m_reaction_num;
        std::unordered_map<int, std::reference_wrapper<double>> m_specie;

        // Production rates - surface and molar gas-phase per mass fraction of respective ions
        double m_bohm_speed_ion;
        double m_sdot_Arion;
        double m_sdot_Armeta;
        double m_wdot_Arion;
        double m_wdot_Armeta;
        double m_q_r;

        // Constants
        double m_eta;   // Correction factor
        double m_R_u;   // Universal gas constant
        double m_W_Ar;  // Molar mass of argon 
        double m_W_Arion;   // Molar mass of argon ion
        double m_W_Armeta;  // Molar mass of argon metastable
        double m_W_E;
        std::vector<double> m_W;    //  Vector containing molar masses
        double m_m_e;   // Electron particle mass
        double m_e;     // Electron charge
        double m_density;
        double m_density_e; // Electron mass density
        double m_k_B;   // Boltzmann constant
        double m_N_A;   // Avogadro's number

        double m_u_e;   // Electron internal energy
        double m_h_e;   // Electron enthalpy
        double m_h_sh;  // Electron enthalpy corresponding to potential drops through grounded sheath
        double m_sigma_en;  // Momentum transfer collision between electrons and background neutrals

        // Concentrations
        double m_conc_conv; // For converting mass fraction to concentration

        // Number density
        double m_numdensity_conv;

        // Initial conditions
        double m_Y_E_init;
        double m_Y_Ar_init;
        double m_Y_Arion_init;
        double m_Y_Armeta_init;
        double m_Te;
        double m_Tgas;

        // Current conditions
        std::vector<double> m_Y;
        std::vector<double> m_X;
        std::vector<double> m_N;    //  Number densities
        double m_Y_Ar;
        double m_Y_Arion;
        double m_Y_Armeta;
        double m_Y_E;
        double m_moles_total;

        //  Differential eqtns
        std::vector<double> m_ydot;

        double m_X_Ar;
        double m_X_Arion;
        double m_X_Armeta;
        double m_X_E;

        // Coefficients for the 3rd Te PDE
        double C, c1, c2, c3, c4;
        double term4;

    private:
};

#endif // INPUTS_H
