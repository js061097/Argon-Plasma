/* -----------------------------------------------------------------
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh and
 *                Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2024, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * Demonstration program for CVODE - direct linear solvers.
 * Two separate problems are solved using both the CV_ADAMS and CV_BDF
 * linear multistep methods in combination with the
 * SUNNONLINSOL_FIXEDPOINT and SUNNONLINSOL_NEWTON nonlinear solver
 * modules:
 * 
 * The NEWTON iteration cases use the following types of Jacobian
 * approximation: (1) dense, user-supplied, (2) dense, difference
 * quotient approximation, (3) diagonal approximation.
 * -----------------------------------------------------------------*/

#include <iostream> 
#include "Inputs.h"
#include "Processor.h"
#include "Model.h"
#include "Solver.h"

int main() 
{
    //  inputs
    Inputs* inputs = new Inputs();
    inputs->initialize();std::cout<<"\n\n"<<"****************"<<"\n\n";

    PostProcessor* pp = new PostProcessor(inputs);  //For handling print statements both intro and end
    //pp->initialize();
    pp->printIntro();

    //  Solver
    Solver* solver = new Solver();

    // model
    Model* model = new Model(inputs, solver);

    //  initialize model
    model->initialize();

    int nerr = model->solve();

    model->finalize();

    std::cout<<"\n\n"<<"****************"<<"\n\n";
    
    //if (inputs.initialize() != 0) return 1; //returns 1 if checkretval at any initialization step returns 1 
    //Solver solver(inputs);  
    //int nerr = solver.solve();
    
    return nerr;
}
