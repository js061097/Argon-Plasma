#ifndef REACTION_H
#define REACTION_H

#include <vector>

class Reaction
{
    public:
        double m_A_coeff;                   // Pre-exponential factor
        double m_ea_coeff;                  // Activation energy
        double m_delta_Ee;                  // Energy lost/electron due to inelastic collision
        double m_k;                         // Rate constant
        std::vector<int> m_reactants;       // Reactants as integers (using enum values)
        std::vector<int> m_coefficients;    // Stoichiometric coefficients of the reactants

        // Default constructor
        Reaction();

        // Parameterized constructor
        Reaction(double A, double ea, double delta_Ee, std::vector<int> reactants, std::vector<int> coefficients);

        // Method to evaluate the rate constant based on temperature Te
        double m_rateCoeff_eval(double Te);

        // Method to initialize the reaction (if needed)
        void initialize();
};

#endif // REACTION_H
