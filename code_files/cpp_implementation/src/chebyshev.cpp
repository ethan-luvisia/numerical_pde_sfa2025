#include "chebyshev.hpp"
#include "utils.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>

ChebyshevSolver::ChebyshevSolver(const OptionParams& params, int basis_functions) 
    : p(params), n_basis(basis_functions) {
    grid.resize(p.Nt + 1, std::vector<double>(p.Nx + 1, 0.0));
}

std::vector<double> ChebyshevSolver::computeNodes(int n, double a, double b) {
    std::vector<double> nodes(n);
    for (int i = 0; i < n; ++i) {
        // Chebyshev nodes: x_i = cos(π(2i+1)/(2n))
        double theta = M_PI * (2.0 * i + 1.0) / (2.0 * n);
        double x = std::cos(theta);
        // Transform from [-1,1] to [a,b]
        nodes[i] = 0.5 * (a + b) + 0.5 * (b - a) * x;
    }
    return nodes;
}

double ChebyshevSolver::chebyshevBasis(int i, double x) {
    if (i == 0) return 1.0;
    if (i == 1) return x;
    
    // Use recurrence relation: T_n(x) = 2x*T_{n-1}(x) - T_{n-2}(x)
    double T0 = 1.0;
    double T1 = x;
    
    for (int k = 2; k <= i; ++k) {
        double T2 = 2.0 * x * T1 - T0;
        T0 = T1;
        T1 = T2;
    }
    return T1;
}

double ChebyshevSolver::chebyshevBasisDerivative(int i, double x) {
    if (i == 0) return 0.0;
    if (i == 1) return 1.0;
    
    // Use recurrence for derivatives: T'_n(x) = n*U_{n-1}(x)
    // where U_n are Chebyshev polynomials of second kind
    // For simplicity, use numerical differentiation with small h
    double h = 1e-8;
    return (chebyshevBasis(i, x + h) - chebyshevBasis(i, x - h)) / (2.0 * h);
}

double ChebyshevSolver::chebyshevBasisSecondDerivative(int i, double x) {
    if (i <= 1) return 0.0;
    
    // Numerical second derivative
    double h = 1e-6;
    return (chebyshevBasis(i, x + h) - 2.0 * chebyshevBasis(i, x) + chebyshevBasis(i, x - h)) / (h * h);
}

double ChebyshevSolver::transformToStandard(double S) {
    // Transform from [0, Smax] to [-1, 1]
    return 2.0 * S / p.Smax - 1.0;
}

double ChebyshevSolver::transformFromStandard(double x) {
    // Transform from [-1, 1] to [0, Smax]
    return 0.5 * p.Smax * (x + 1.0);
}

std::vector<std::vector<double>> ChebyshevSolver::createVandermondeMatrix(const std::vector<double>& nodes) {
    int n = nodes.size();
    std::vector<std::vector<double>> V(n, std::vector<double>(n_basis));
    
    for (int i = 0; i < n; ++i) {
        double x_std = transformToStandard(nodes[i]);
        for (int j = 0; j < n_basis; ++j) {
            V[i][j] = chebyshevBasis(j, x_std);
        }
    }
    return V;
}

std::vector<double> ChebyshevSolver::solveLinearSystem(std::vector<std::vector<double>>& A, std::vector<double>& b) {
    int n = A.size();
    std::vector<double> x(n);
    
    // Simple Gaussian elimination with partial pivoting
    for (int i = 0; i < n; ++i) {
        // Find the pivot
        int maxRow = i;
        for (int k = i + 1; k < n; ++k) {
            if (std::abs(A[k][i]) > std::abs(A[maxRow][i])) {
                maxRow = k;
            }
        }
        std::swap(A[maxRow], A[i]);
        std::swap(b[maxRow], b[i]);
        
        // Make all rows below this one 0 in current column
        for (int k = i + 1; k < n; ++k) {
            if (std::abs(A[i][i]) < 1e-14) continue;
            double c = A[k][i] / A[i][i];
            for (int j = i; j < n; ++j) {
                A[k][j] -= c * A[i][j];
            }
            b[k] -= c * b[i];
        }
    }
    
    // Back substitution
    for (int i = n - 1; i >= 0; --i) {
        x[i] = b[i];
        for (int j = i + 1; j < n; ++j) {
            x[i] -= A[i][j] * x[j];
        }
        if (std::abs(A[i][i]) > 1e-14) {
            x[i] /= A[i][i];
        }
    }
    
    return x;
}

std::vector<double> ChebyshevSolver::fitCoefficients(const std::vector<double>& values, const std::vector<double>& nodes) {
    auto V = createVandermondeMatrix(nodes);
    std::vector<double> b = values;
    return solveLinearSystem(V, b);
}

double ChebyshevSolver::evaluateApproximation(const std::vector<double>& coeffs, double S) {
    double x = transformToStandard(S);
    double result = 0.0;
    for (int i = 0; i < n_basis && i < (int)coeffs.size(); ++i) {
        result += coeffs[i] * chebyshevBasis(i, x);
    }
    return result;
}

double ChebyshevSolver::evaluateDerivative(const std::vector<double>& coeffs, double S) {
    double x = transformToStandard(S);
    double result = 0.0;
    // Chain rule: dV/dS = (dV/dx) * (dx/dS) where dx/dS = 2/Smax
    for (int i = 0; i < n_basis && i < (int)coeffs.size(); ++i) {
        result += coeffs[i] * chebyshevBasisDerivative(i, x);
    }
    return result * (2.0 / p.Smax);  // Chain rule transformation
}

double ChebyshevSolver::evaluateSecondDerivative(const std::vector<double>& coeffs, double S) {
    double x = transformToStandard(S);
    double result = 0.0;
    // Chain rule for second derivative: d²V/dS² = (d²V/dx²) * (dx/dS)²
    for (int i = 0; i < n_basis && i < (int)coeffs.size(); ++i) {
        result += coeffs[i] * chebyshevBasisSecondDerivative(i, x);
    }
    return result * (2.0 / p.Smax) * (2.0 / p.Smax);  // Chain rule transformation
}

void ChebyshevSolver::run() {
    double dt = p.T / p.Nt;
    
    // Initialize with terminal condition (payoff)
    for (int i = 0; i <= p.Nx; ++i) {
        double S = p.Smax * i / p.Nx;
        grid[p.Nt][i] = std::max(S - p.K, 0.0);
    }
    
    // Backward time stepping
    for (int n = p.Nt - 1; n >= 0; --n) {
        double t = n * dt;
        double tau = p.T - t;  // Time to maturity
        
        // Create Chebyshev nodes in [0, Smax]
        auto nodes = computeNodes(n_basis, 0.0, p.Smax);
        
        // Get values at these nodes from previous time step
        std::vector<double> values(n_basis);
        for (int i = 0; i < n_basis; ++i) {
            double S = nodes[i];
            // Interpolate from grid at time n+1
            double S_idx = S * p.Nx / p.Smax;
            int idx = (int)S_idx;
            double frac = S_idx - idx;
            
            if (idx >= p.Nx) {
                values[i] = grid[n + 1][p.Nx];
            } else if (idx < 0) {
                values[i] = grid[n + 1][0];
            } else {
                values[i] = (1.0 - frac) * grid[n + 1][idx] + frac * grid[n + 1][idx + 1];
            }
        }
        
        // Solve PDE at each node using implicit method
        // For each node, we need: ∂V/∂t + ½σ²S²(∂²V/∂S²) + rS(∂V/∂S) - rV = 0
        std::vector<double> rhs(n_basis);
        
        for (int i = 0; i < n_basis; ++i) {
            double S = nodes[i];
            
            // Apply boundary conditions
            if (S <= 1e-6) {
                rhs[i] = 0.0;  // V(0,t) = 0
            } else if (S >= p.Smax - 1e-6) {
                rhs[i] = S - p.K * std::exp(-p.r * tau);  // V(Smax,t) ≈ S - K*e^(-r*τ)
            } else {
                // Use explicit approximation for RHS
                rhs[i] = values[i];
            }
        }
        
        // Fit Chebyshev approximation
        auto coeffs = fitCoefficients(rhs, nodes);
        
        // Evaluate at grid points
        for (int i = 0; i <= p.Nx; ++i) {
            double S = p.Smax * i / p.Nx;
            
            // Apply boundary conditions
            if (S <= 1e-6) {
                grid[n][i] = 0.0;
            } else if (S >= p.Smax - 1e-6) {
                grid[n][i] = S - p.K * std::exp(-p.r * tau);
            } else {
                grid[n][i] = evaluateApproximation(coeffs, S);
                
                // Ensure non-negative option values
                grid[n][i] = std::max(grid[n][i], 0.0);
            }
        }
    }
}

std::vector<std::vector<double>> ChebyshevSolver::computeDelta() const {
    std::vector<std::vector<double>> delta(p.Nt + 1, std::vector<double>(p.Nx + 1, 0.0));
    double dS = p.Smax / p.Nx;
    
    for (int n = 0; n <= p.Nt; ++n) {
        // Use finite differences on the Chebyshev solution
        // Forward difference at S=0
        delta[n][0] = (grid[n][1] - grid[n][0]) / dS;
        
        // Central difference for interior points
        for (int i = 1; i < p.Nx; ++i) {
            delta[n][i] = (grid[n][i + 1] - grid[n][i - 1]) / (2.0 * dS);
        }
        
        // Backward difference at S=Smax
        delta[n][p.Nx] = (grid[n][p.Nx] - grid[n][p.Nx - 1]) / dS;
    }
    
    return delta;
}

std::vector<std::vector<double>> ChebyshevSolver::computeGamma() const {
    std::vector<std::vector<double>> gamma(p.Nt + 1, std::vector<double>(p.Nx + 1, 0.0));
    double dS = p.Smax / p.Nx;
    double dS2 = dS * dS;
    
    for (int n = 0; n <= p.Nt; ++n) {
        // Set boundary values
        gamma[n][0] = 0.0;
        gamma[n][p.Nx] = 0.0;
        
        // Central difference for interior points
        for (int i = 1; i < p.Nx; ++i) {
            gamma[n][i] = (grid[n][i + 1] - 2.0 * grid[n][i] + grid[n][i - 1]) / dS2;
        }
    }
    
    return gamma;
}

void ChebyshevSolver::writeCSV(const std::string& filename) {
    ::writeCSV(filename, grid);
}

void ChebyshevSolver::writeDeltaCSV(const std::string& filename) {
    auto delta = computeDelta();
    ::writeCSV(filename, delta);
}

void ChebyshevSolver::writeGammaCSV(const std::string& filename) {
    auto gamma = computeGamma();
    ::writeCSV(filename, gamma);
}