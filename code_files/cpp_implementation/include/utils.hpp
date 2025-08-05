#ifndef UTILS_HPP
#define UTILS_HPP

#include <vector>
#include <string>

struct OptionParams {
    double Smax;    // Maximum stock price
    double K;       // Strike price
    double r;       // Risk-free rate
    double sigma;   // Volatility
    double T;       // Time to maturity
    int Nx;         // Number of spatial grid points
    int Nt;         // Number of temporal grid points
};

// CSV writing utility
void writeCSV(const std::string& filename, const std::vector<std::vector<double>>& grid);

#endif