#include "Reaction.h"

Reaction::Reaction() : m_A_coeff(0.0), m_ea_coeff(0.0), m_delta_Ee(0.0)
{

}

Reaction::Reaction(double A, double ea, double delta_Ee, std::vector<int> reactants, std::vector<int> coefficients)
    : m_A_coeff(A), m_ea_coeff(ea), m_delta_Ee(delta_Ee), m_reactants(reactants), m_coefficients(coefficients)
{

}

void Reaction::initialize()
{
    // Initialize reaction-specific properties here
}

double Reaction::m_rateCoeff_eval(double Te)
{
    return m_A_coeff * exp(-m_ea_coeff / Te);
}
