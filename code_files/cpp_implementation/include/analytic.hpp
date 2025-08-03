#ifndef ANALYTIC_HPP
#define ANALYTIC_HPP

#include <vector>

struct OptionParams {
    double Smax, K, r, sigma, T;
    int Nx, Nt;
};

std::vector<std::vector<double>> computeAnalyticGrid(const OptionParams& params);
void writeAnalyticCSV(const std::string& filename, const OptionParams& params);

#endif