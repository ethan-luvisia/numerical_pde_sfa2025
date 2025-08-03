#include "analytic.hpp"
#include "utils.hpp"
#include <cmath>

// Normal CDF using standard C++ only
double blackScholesFormula(double S, double K, double r, double sigma, double T) {
    if (T <= 0.0) return std::max(S - K, 0.0);
    double d1 = (std::log(S / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T));
    double d2 = d1 - sigma * std::sqrt(T);
    auto N = [](double x) { return 0.5 * std::erfc(-x / std::sqrt(2)); };
    return S * N(d1) - K * std::exp(-r * T) * N(d2);
}


std::vector<std::vector<double>> computeAnalyticGrid(const OptionParams& p) {
    std::vector<std::vector<double>> grid(p.Nt + 1, std::vector<double>(p.Nx + 1));
    for (int n = 0; n <= p.Nt; ++n) {
        double t = p.T * n / p.Nt;
        for (int i = 0; i <= p.Nx; ++i) {
            double S = p.Smax * i / p.Nx;
            grid[n][i] = blackScholesFormula(S, p.K, p.r, p.sigma, p.T - t);
        }
    }
    return grid;
}
