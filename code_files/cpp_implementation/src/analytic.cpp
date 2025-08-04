#include "analytic.hpp"
#include "utils.hpp"
#include <cmath>

// Standard normal CDF
double normalCDF(double x) {
    return 0.5 * std::erfc(-x / std::sqrt(2.0));
}

// Standard normal PDF
double normalPDF(double x) {
    return std::exp(-0.5 * x * x) / std::sqrt(2.0 * M_PI);
}

double blackScholesFormula(double S, double K, double r, double sigma, double T) {
    if (T <= 0.0) return std::max(S - K, 0.0);
    
    double d1 = (std::log(S / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T));
    double d2 = d1 - sigma * std::sqrt(T);
    
    return S * normalCDF(d1) - K * std::exp(-r * T) * normalCDF(d2);
}

double blackScholesDelta(double S, double K, double r, double sigma, double T) {
    if (T <= 0.0) return (S > K) ? 1.0 : 0.0;
    
    double d1 = (std::log(S / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T));
    return normalCDF(d1);
}

double blackScholesGamma(double S, double K, double r, double sigma, double T) {
    if (T <= 0.0) return 0.0;
    
    double d1 = (std::log(S / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T));
    return normalPDF(d1) / (S * sigma * std::sqrt(T));
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

std::vector<std::vector<double>> computeAnalyticDelta(const OptionParams& p) {
    std::vector<std::vector<double>> delta(p.Nt + 1, std::vector<double>(p.Nx + 1));
    for (int n = 0; n <= p.Nt; ++n) {
        double t = p.T * n / p.Nt;
        for (int i = 0; i <= p.Nx; ++i) {
            double S = p.Smax * i / p.Nx;
            delta[n][i] = blackScholesDelta(S, p.K, p.r, p.sigma, p.T - t);
        }
    }
    return delta;
}

std::vector<std::vector<double>> computeAnalyticGamma(const OptionParams& p) {
    std::vector<std::vector<double>> gamma(p.Nt + 1, std::vector<double>(p.Nx + 1));
    for (int n = 0; n <= p.Nt; ++n) {
        double t = p.T * n / p.Nt;
        for (int i = 0; i <= p.Nx; ++i) {
            double S = p.Smax * i / p.Nx;
            gamma[n][i] = blackScholesGamma(S, p.K, p.r, p.sigma, p.T - t);
        }
    }
    return gamma;
}