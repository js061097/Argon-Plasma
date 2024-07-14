#ifndef INPUTS_H
#define INPUTS_H

#include <stdio.h>
#include <stdlib.h>
#include <iostream>


/* Equations and IVs */
#define P1_ETA        SUN_RCONST(3.0)
#define P1_NOUT       4
#define P1_T0         SUN_RCONST(0.0)
#define P1_T1         SUN_RCONST(10.39283880203)
#define P1_DTOUT      SUN_RCONST(2.214773875)
#define P1_TOL_FACTOR SUN_RCONST(1.0e4)



/*Classes*/

//To take in custom inputs and initializing 
class Inputs 
{
    public:
        Inputs();
        ~Inputs();
        int Neq;

        void readInputs();

        void initialize();
        void setEquations(int);
        void setTime(double,double);
        void finalize();

        //  solver related inputs
        double m_rtol;
        double m_atol;
        
        //  timings
        double m_tbegin;
        double m_tend;
        double m_dt;
        
        //  reactor parameters
        double m_volume;
        double m_surface_area;

        //  reactor operating conditions
        double m_pressure;  // (units: Pa)
        double m_residtime; // (units: sec)
        double m_Pinp;  // input electron power (units: W)

        //  gas-chemistry related
        double m_A_coeff;
        double m_beta_coeff;
        double m_ea_coeff;

        //  initial conditions
        double m_N_E_init;
        double m_N_Ar_init;
        double m_N_Arion_init;
        double m_N_Armeta_init;
        double m_Te;
        double m_Tgas;

        //

    /*Variables used in Problem 1 DirectDemo main*/
    //sunrealtype reltol, abstol, t, tout, ero, er;
    //int retval, temp_retval, iout, nerr;
    

    private:
    

};


#endif // INPUTS_H




