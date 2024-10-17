
#include "Processor.h"
#include <iostream>
#include "Solver.h"
#include "Inputs.h"
#include <fstream>

Processor::Processor()
{
   m_csv_initialized = false;
}

Processor::Processor(Inputs* inputs)
{
    m_inputs = inputs;
    m_csv_initialized = false;
}

Processor::~Processor()
{

}

void Processor::initialize()
{

}

void Processor::printIntro() 
{
    printf("Demonstration program for CVODE package - direct linear solvers\n");
    printf("\n\n");
    printf("Problem: Argon Plasma Physics\n");
#if defined(SUNDIALS_EXTENDED_PRECISION)
    printf(" neq = %d,  reltol = %.2Lg,  abstol = %.2Lg", m_inputs->m_Neq, m_inputs->m_rtol, m_inputs->m_atol);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    printf(" neq = %d,  reltol = %.2g,  abstol = %.2g", m_inputs->m_Neq, m_inputs->m_rtol, m_inputs->m_atol);
#else
    printf(" neq = %d,  reltol = %.2g,  abstol = %.2g", m_inputs->m_Neq, m_inputs->m_rtol, m_inputs->m_atol);
#endif
}

void Processor::printHeader() 
{
    printf("\n     t           Yar              Yar+              Yar*           Ye-           Nar+           Nar*            Te           ydot0           ydot1           ydot2            term4         qu     hu \n");
}

void Processor::printOutput(double t, double y0, double y1, double y2, double y3, double n1, double n2, double temp, double ydot1, double ydot2, double ydot3, double term4, int qu, sunrealtype hu) 
{
    if (!m_csv_initialized) {
            initializeCSV("output.csv");
            m_csv_initialized = true;  // Set the flag to prevent re-initialization
        }
#if defined(SUNDIALS_EXTENDED_PRECISION)
    printf("%10.5Lf    %12.5Le   %12.5Le   %12.5Le   %12.5e   %12.5Le   %12.5e    %12.5e    %12.5e    %12.5e    %12.5e    %12.5e    %2d    %6.4Le\n", t, y0, y1, y2, y3, n1, n2, temp, ydot1, ydot2, ydot3, term4, qu, hu);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    printf("%10.5Lf    %12.5Le   %12.5Le   %12.5Le   %12.5e   %12.5Le   %12.5e    %12.5e    %12.5e    %12.5e    %12.5e    %12.5e    %2d    %6.4Le\n", t, y0, y1, y2, y3, n1, n2, temp, ydot1, ydot2, ydot3, term4, qu, hu);
#else
    printf("%10.5Lf    %12.5Le   %12.5Le   %12.5Le   %12.5e   %12.5Le   %12.5e    %12.5e    %12.5e    %12.5e    %12.5e    %12.5e    %2d    %6.4Le\n", t, y0, y1, y2, y3, n1, n2, temp, ydot1, ydot2, ydot3, term4, qu, hu);
#endif
appendToCSV("output.csv", t, y0, y1, y2, y3, n1, n2, temp, ydot1, ydot2, ydot3, term4, qu, hu);

}

void Processor::initializeCSV(const std::string& filename) {
    std::ofstream file(filename);  // Open in write mode to create or overwrite
    if (file.is_open()) {
        file << "t,Yar,Yar+,Yar*,Ye-,Nar+,Nar*,Te,ydot0,ydot1,ydot2,term4,qu,hu\n";
        file.close();
    } else {
        std::cerr << "Unable to create or open CSV file: " << filename << std::endl;
    }
}

void Processor::printErrOutput(sunrealtype tol_factor) 
{
#if defined(SUNDIALS_EXTENDED_PRECISION)
    printf("\n\n Error exceeds %Lg * tolerance \n\n", tol_factor);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    printf("\n\n Error exceeds %g * tolerance \n\n", tol_factor);
#else
    printf("\n\n Error exceeds %g * tolerance \n\n", tol_factor);
#endif
}

void Processor::appendToCSV(const std::string& filename, double t, double y0, double y1, 
                            double y2, double y3, double n1, double n2, double temp, 
                            double ydot1, double ydot2, double ydot3, double term4, 
                            int qu, sunrealtype hu) {
    std::ofstream file(filename, std::ios::app);  // Open in append mode
    if (file.is_open()) {
        file << std::fixed << std::setprecision(10)
             << t << "," << y0 << "," << y1 << "," << y2 << "," << y3 << ","
             << n1 << "," << n2 << "," << temp << "," << ydot1 << "," << ydot2 << ","
             << ydot3 << "," << term4 << "," << qu << "," << hu << "\n";
        file.close();
    } else {
        std::cerr << "Unable to open CSV file for appending: " << filename << std::endl;
    }
}

#include <fstream>  // Required for file handling

void Processor::createPythonScript() {
    std::ofstream pythonFile("plot_data.py");  // Create the Python file
    if (pythonFile.is_open()) {
        pythonFile << R"(
import pandas as pd
import matplotlib.pyplot as plt
import sys

# Check if the correct number of arguments are provided
if len(sys.argv) != 3:
    print("Usage: python plot_data.py <X-axis-column> <Y-axis-column>")
    sys.exit(1)

# Get column names for X and Y axes from command-line arguments
x_col = sys.argv[1]
y_col = sys.argv[2]

# Read the CSV file
try:
    df = pd.read_csv('output.csv')
except FileNotFoundError:
    print("Error: 'output.csv' not found. Please run your C++ program to generate the CSV.")
    sys.exit(1)
print(df['t'][0]==df['t'][1])
# Check if the provided columns exist in the CSV
if x_col not in df.columns or y_col not in df.columns:
    print(f"Error: Columns '{x_col}' or '{y_col}' not found in the CSV file.")
    print(f"Available columns: {', '.join(df.columns)}")
    sys.exit(1)

# Plot the data
plt.figure(figsize=(8, 5))
plt.plot(df[x_col], df[y_col], linestyle='-')
plt.xlabel(x_col)
plt.ylabel(y_col)
plt.title(f'{y_col} vs {x_col}')
plt.grid(True)
plt.show()
        )";
        pythonFile.close();
        //std::cout << "Python script 'plot_data.py' created successfully." << std::endl;
    } else {
        std::cerr << "Unable to create 'plot_data.py'." << std::endl;
    }
}


void Processor::printFinalStats(void* cvode_mem, sunrealtype ero) 
{
    long int lenrw, leniw, lenrwLS, leniwLS;
    long int nst, nfe, nsetups, nni, ncfn, netf, nje, nfeLS;
    int retval;

    retval = CVodeGetWorkSpace(cvode_mem, &lenrw, &leniw);
    checkReturnValue(&retval, "CVodeGetWorkSpace", 1);
    retval = CVodeGetNumSteps(cvode_mem, &nst);
    checkReturnValue(&retval, "CVodeGetNumSteps", 1);
    retval = CVodeGetNumRhsEvals(cvode_mem, &nfe);
    checkReturnValue(&retval, "CVodeGetNumRhsEvals", 1);
    retval = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
    checkReturnValue(&retval, "CVodeGetNumLinSolvSetups", 1);
    retval = CVodeGetNumErrTestFails(cvode_mem, &netf);
    checkReturnValue(&retval, "CVodeGetNumErrTestFails", 1);
    retval = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
    checkReturnValue(&retval, "CVodeGetNumNonlinSolvIters", 1);
    retval = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
    checkReturnValue(&retval, "CVodeGetNumNonlinSolvConvFails", 1);

    printf("\n Final statistics for this run:\n\n");
    printf(" CVode real workspace length              = %4ld \n", lenrw);
    printf(" CVode integer workspace length           = %4ld \n", leniw);
    printf(" Number of steps                          = %4ld \n", nst);
    printf(" Number of f-s                            = %4ld \n", nfe);
    printf(" Number of setups                         = %4ld \n", nsetups);
    printf(" Number of nonlinear iterations           = %4ld \n", nni);
    printf(" Number of nonlinear convergence failures = %4ld \n", ncfn);
    printf(" Number of error test failures            = %4ld \n\n", netf);

    nje    = nsetups;
    retval = CVodeGetNumRhsEvals(cvode_mem, &nfeLS);
    checkReturnValue(&retval, "CVDiagGetNumRhsEvals", 1);
    retval = CVodeGetWorkSpace(cvode_mem, &lenrwLS, &leniwLS);
    checkReturnValue(&retval, "CVDiagGetWorkSpace", 1);
    
    printf(" Linear solver real workspace length      = %4ld \n", lenrwLS);
    printf(" Linear solver integer workspace length   = %4ld \n", leniwLS);
    printf(" Number of Jacobian evaluations           = %4ld \n", nje);
    printf(" Number of f evals. in linear solver      = %4ld \n\n", nfeLS);

#if defined(SUNDIALS_EXTENDED_PRECISION)
    printf(" Error overrun = %.3Lf \n", ero);
#else
    printf(" Error overrun = %.3f \n", ero);
#endif
}

void Processor::printErrInfo(int nerr) 
{
    printf("\n\n-------------------------------------------------------------");
    printf("\n-------------------------------------------------------------");
    printf("\n\n Number of errors encountered = %d \n", nerr);
}

int checkReturnValue(void* returnvalue, const char* funcname, int opt) 
{
    int* retval;

    if (opt == 0 && returnvalue == NULL) {
        fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n", funcname);
        return 1;
    } else if (opt == 1) {
        retval = (int*)returnvalue;
        if (*retval < 0) {
            fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with retval = %d\n\n", funcname, *retval);
            return 1;
        }
    } else if (opt == 2 && returnvalue == NULL) {
        fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n", funcname);
        return 1;
    }

    return 0;
}