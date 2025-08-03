#include "solver.hpp"
#include "utils.hpp"
#include <cmath>

CrankNicolsonSolver::CrankNicolsonSolver(const OptionParams& params) : p(params) {
    grid.resize(p.Nt + 1, std::vector<double>(p.Nx + 1, 0.0));
    for (int i = 0; i <= p.Nx; ++i) {
        double S = p.Smax * i / p.Nx;
        grid[p.Nt][i] = std::max(S - p.K, 0.0);  // final condition: payoff
    }
}

void CrankNicolsonSolver::run() {
    double dt = p.T / p.Nt;
    double dS = p.Smax / p.Nx;

    std::vector<double> a(p.Nx - 1), b(p.Nx - 1), c(p.Nx - 1), d(p.Nx - 1), x(p.Nx - 1);

    for (int n = p.Nt - 1; n >= 0; --n) {
        double t = n * dt;

        for (int i = 1; i < p.Nx; ++i) {
            double S = i * dS;
            double sigma2 = p.sigma * p.sigma;

            // Finite difference coefficients
            double alpha = 0.25 * dt * (sigma2 * S * S / (dS * dS) - p.r * S / dS);
            double beta = -0.5 * dt * (sigma2 * S * S / (dS * dS) + p.r);
            double gamma = 0.25 * dt * (sigma2 * S * S / (dS * dS) + p.r * S / dS);

            // LHS matrix coefficients (implicit part)
            a[i - 1] = -alpha;                    // coefficient of V[n][i-1]
            b[i - 1] = 1.0 - beta;               // coefficient of V[n][i]
            c[i - 1] = -gamma;                   // coefficient of V[n][i+1]

            // RHS vector (explicit part)
            d[i - 1] = alpha * grid[n + 1][i - 1] + 
                      (1.0 + beta) * grid[n + 1][i] + 
                      gamma * grid[n + 1][i + 1];
        }

        // Apply boundary conditions to the system
        // Lower boundary: V(0,t) = 0
        d[0] -= a[0] * 0.0;  // grid[n][0] = 0
        
        // Upper boundary: V(Smax,t) = Smax - K*exp(-r*tau)
        double tau = p.T - t;
        double upper_bc = p.Smax - p.K * std::exp(-p.r * tau);
        d[p.Nx - 2] -= c[p.Nx - 2] * upper_bc;

        // Solve tridiagonal system
        thomasSolve(a, b, c, d, x);
        
        // Set interior points
        for (int i = 1; i < p.Nx; ++i) {
            grid[n][i] = x[i - 1];
        }

        // Set boundary conditions
        grid[n][0] = 0.0;
        grid[n][p.Nx] = upper_bc;
    }
}

void CrankNicolsonSolver::thomasSolve(std::vector<double>& a, std::vector<double>& b,
                                      std::vector<double>& c, std::vector<double>& d,
                                      std::vector<double>& x) {
    int n = d.size();
    std::vector<double> c_star(n), d_star(n);

    // Forward elimination
    c_star[0] = c[0] / b[0];
    d_star[0] = d[0] / b[0];
    
    for (int i = 1; i < n; ++i) {
        double m = 1.0 / (b[i] - a[i] * c_star[i - 1]);
        c_star[i] = c[i] * m;
        d_star[i] = (d[i] - a[i] * d_star[i - 1]) * m;
    }

    // Back substitution
    x[n - 1] = d_star[n - 1];
    for (int i = n - 2; i >= 0; --i) {
        x[i] = d_star[i] - c_star[i] * x[i + 1];
    }
}

void CrankNicolsonSolver::writeCSV(const std::string& filename) {
    ::writeCSV(filename, grid);  // uses utils.cpp CSV writer
}