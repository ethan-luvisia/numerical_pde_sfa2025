#ifndef CHEBYSHEV_HPP
#define CHEBYSHEV_HPP

#include "utils.hpp"
#include <vector>
#include <functional>

class ChebyshevSolver {
private:
    OptionParams p;
    std::vector<std::vector<double>> grid;
    int n_basis;  // Number of Chebyshev basis functions
    
    // Chebyshev nodes and basis functions
    std::vector<double> computeNodes(int n, double a = -1.0, double b = 1.0);
    double chebyshevBasis(int i, double x);
    double chebyshevBasisDerivative(int i, double x);
    double chebyshevBasisSecondDerivative(int i, double x);
    
    // Transform between [0, Smax] and [-1, 1]
    double transformToStandard(double S);
    double transformFromStandard(double x);
    
    // Matrix operations
    std::vector<std::vector<double>> createVandermondeMatrix(const std::vector<double>& nodes);
    std::vector<double> solveLinearSystem(std::vector<std::vector<double>>& A, std::vector<double>& b);
    
    // PDE solving
    std::vector<double> fitCoefficients(const std::vector<double>& values, const std::vector<double>& nodes);
    double evaluateApproximation(const std::vector<double>& coeffs, double S);
    double evaluateDerivative(const std::vector<double>& coeffs, double S);
    double evaluateSecondDerivative(const std::vector<double>& coeffs, double S);

public:
    ChebyshevSolver(const OptionParams& params, int basis_functions = 20);
    
    void run();
    void writeCSV(const std::string& filename);
    
    // Greeks computation
    std::vector<std::vector<double>> computeDelta() const;
    std::vector<std::vector<double>> computeGamma() const;
    void writeDeltaCSV(const std::string& filename);
    void writeGammaCSV(const std::string& filename);
    
    // Access to grid for analysis
    const std::vector<std::vector<double>>& getGrid() const { return grid; }
};

#endif