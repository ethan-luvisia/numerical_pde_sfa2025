#pragma once
#include "utils.hpp" // This has OptionsParams struct definition
#include <vector>

double blackScholesFormula(double S, double K, double r, double sigma, double T);
double blackScholesDelta(double S, double K, double r, double sigma, double T);
double blackScholesGamma(double S, double K, double r, double sigma, double T);

std::vector<std::vector<double>> computeAnalyticGrid(const OptionParams& params);
std::vector<std::vector<double>> computeAnalyticDelta(const OptionParams& params);
std::vector<std::vector<double>> computeAnalyticGamma(const OptionParams& params);