#include "analytic.hpp"
#include "solver.hpp"
#include "utils.hpp"
#include <iostream>

int main() {
    OptionParams p = {200, 100, 0.05, 0.2, 1.0, 100, 1000};

    std::cout << "Computing analytical solutions..." << std::endl;
    
    // Compute analytical option prices
    auto analyticGrid = computeAnalyticGrid(p);
    writeCSV("../data/output_analytic.csv", analyticGrid);
    
    // Compute analytical Greeks
    auto analyticDelta = computeAnalyticDelta(p);
    writeCSV("../data/delta_analytic.csv", analyticDelta);
    
    auto analyticGamma = computeAnalyticGamma(p);
    writeCSV("../data/gamma_analytic.csv", analyticGamma);

    std::cout << "Computing finite difference solutions..." << std::endl;
    
    // Compute finite difference option prices
    CrankNicolsonSolver solver(p);
    solver.run();
    solver.writeCSV("../data/output_fd.csv");
    
    // Compute finite difference Greeks
    solver.writeDeltaCSV("../data/delta_fd.csv");
    solver.writeGammaCSV("../data/gamma_fd.csv");

    std::cout << "All computations complete!" << std::endl;
    std::cout << "Files generated:" << std::endl;
    std::cout << "  - output_analytic.csv / output_fd.csv (option prices)" << std::endl;
    std::cout << "  - delta_analytic.csv / delta_fd.csv (delta values)" << std::endl;
    std::cout << "  - gamma_analytic.csv / gamma_fd.csv (gamma values)" << std::endl;
    
    return 0;
}