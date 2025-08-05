#include "analytic.hpp"
#include "solver.hpp"
#include "chebyshev.hpp"
#include "utils.hpp"
#include <iostream>
#include <chrono>

int main() {
    OptionParams p = {200, 100, 0.05, 0.2, 1.0, 100, 1000};

    std::cout << "=== European Call Option Pricing Comparison ===" << std::endl;
    std::cout << "Parameters:" << std::endl;
    std::cout << "  S_max: " << p.Smax << ", K: " << p.K << std::endl;
    std::cout << "  r: " << p.r << ", σ: " << p.sigma << ", T: " << p.T << std::endl;
    std::cout << "  Grid: " << p.Nx + 1 << " × " << p.Nt + 1 << " points" << std::endl;
    std::cout << std::endl;

    // Timing variables
    auto start_time = std::chrono::high_resolution_clock::now();
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

    // 1. Analytical Solution
    std::cout << "1. Computing analytical (Black-Scholes) solutions..." << std::endl;
    start_time = std::chrono::high_resolution_clock::now();
    
    auto analyticGrid = computeAnalyticGrid(p);
    writeCSV("../data/output_analytic.csv", analyticGrid);
    
    auto analyticDelta = computeAnalyticDelta(p);
    writeCSV("../data/delta_analytic.csv", analyticDelta);
    
    auto analyticGamma = computeAnalyticGamma(p);
    writeCSV("../data/gamma_analytic.csv", analyticGamma);
    
    end_time = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    std::cout << "   Completed in " << duration.count() << " ms" << std::endl;

    // 2. Crank-Nicolson Finite Difference
    std::cout << "2. Computing Crank-Nicolson finite difference solutions..." << std::endl;
    start_time = std::chrono::high_resolution_clock::now();
    
    CrankNicolsonSolver cnSolver(p);
    cnSolver.run();
    cnSolver.writeCSV("../data/output_cn.csv");
    cnSolver.writeDeltaCSV("../data/delta_cn.csv");
    cnSolver.writeGammaCSV("../data/gamma_cn.csv");
    
    end_time = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    std::cout << "   Completed in " << duration.count() << " ms" << std::endl;

    // 3. Chebyshev Polynomial Approximation
    std::cout << "3. Computing Chebyshev polynomial approximation solutions..." << std::endl;
    start_time = std::chrono::high_resolution_clock::now();
    
    ChebyshevSolver chebSolver(p, 20);  // Using 20 basis functions as in the paper
    chebSolver.run();
    chebSolver.writeCSV("../data/output_cheb.csv");
    chebSolver.writeDeltaCSV("../data/delta_cheb.csv");
    chebSolver.writeGammaCSV("../data/gamma_cheb.csv");
    
    end_time = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    std::cout << "   Completed in " << duration.count() << " ms" << std::endl;

    std::cout << std::endl;
    std::cout << "All computations complete!" << std::endl;
    std::cout << std::endl;
    std::cout << "Files generated:" << std::endl;
    std::cout << "Option Prices:" << std::endl;
    std::cout << "  - output_analytic.csv  (Black-Scholes analytical)" << std::endl;
    std::cout << "  - output_cn.csv        (Crank-Nicolson FDM)" << std::endl;
    std::cout << "  - output_cheb.csv      (Chebyshev polynomial)" << std::endl;
    std::cout << std::endl;
    std::cout << "Delta (∂V/∂S):" << std::endl;
    std::cout << "  - delta_analytic.csv   (Black-Scholes analytical)" << std::endl;
    std::cout << "  - delta_cn.csv         (Crank-Nicolson FDM)" << std::endl;
    std::cout << "  - delta_cheb.csv       (Chebyshev polynomial)" << std::endl;
    std::cout << std::endl;
    std::cout << "Gamma (∂²V/∂S²):" << std::endl;
    std::cout << "  - gamma_analytic.csv   (Black-Scholes analytical)" << std::endl;
    std::cout << "  - gamma_cn.csv         (Crank-Nicolson FDM)" << std::endl;
    std::cout << "  - gamma_cheb.csv       (Chebyshev polynomial)" << std::endl;
    std::cout << std::endl;
    std::cout << "Use the Python analysis script to visualize and compare results!" << std::endl;
    
    return 0;
}