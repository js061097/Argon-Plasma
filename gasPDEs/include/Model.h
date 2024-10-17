#ifndef MODEL_H
#define MODEL_H

#include "Inputs.h"
#include "Solver.h" 
#include<vector>

class Solver; //Forward declaration - error resulting in not declaring prior to usage during make


//UserFunction prototype

class Model
{
    public:
        Model(Inputs* inputs, Solver* solver);
        ~Model();
        static int f1(double t, N_Vector y, N_Vector ydot, void* user_data);
        void initialize();
        int solve();
        void finalize();
    
        //std::vector<std::vector<int>> c;
        Inputs* m_inputs;   // local copy
        Solver* m_solver;   // local copy
    private:
};

#endif