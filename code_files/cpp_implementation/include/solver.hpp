#pragma once
#include <vector>
#include <string>
#include "utils.hpp"  // for OptionParams

class CrankNicolsonSolver {
public:
    CrankNicolsonSolver(const OptionParams& p);
    void run();
    void writeCSV(const std::string& filename);

    std::vector<std::vector<double>> computeDelta() const;
    std::vector<std::vector<double>> computeGamma() const;
    void writeDeltaCSV(const std::string& filename);
    void writeGammaCSV(const std::string& filename);
    
private:
    OptionParams p;
    std::vector<std::vector<double>> grid;
    void thomasSolve(std::vector<double>& a, std::vector<double>& b,
                     std::vector<double>& c, std::vector<double>& d,
                     std::vector<double>& x);
};
