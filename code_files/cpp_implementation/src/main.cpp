#include "analytic.hpp"
#include "solver.hpp"
#include "utils.hpp"
#include <iostream>

int main() {
    OptionParams p = {200, 100, 0.05, 0.2, 1.0, 100, 1000};

    auto analyticGrid = computeAnalyticGrid(p);
    std::cout << "Writing analytic grid..." << std::endl;
    writeCSV("../data/output_analytic.csv", analyticGrid);

    CrankNicolsonSolver solver(p);
    solver.run();
    std::cout << "Writing FDM grid..." << std::endl;
    solver.writeCSV("../data/output_fd.csv");

    std::cout << "Done." << std::endl;
    return 0;
}