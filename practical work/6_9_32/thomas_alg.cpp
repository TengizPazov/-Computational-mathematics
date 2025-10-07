//Алгоритм прогонки(Томаса)
#include <vector>
#include <iostream>
std::vector<double> thomas_algorithm(const std::vector<double>& h, const std::vector<double>& u) {
    int N = h.size();
    int n_unknown = N - 1;
    std::vector<double> xi(n_unknown, 0.0);
    std::vector<double> eta(n_unknown, 0.0);
    std::vector<double> c(N + 1, 0.0);
    if (n_unknown >= 1) {
        double beta1 = 2.0;
        double gamma1 = h[1] / (h[0] + h[1]);
        
        xi[0] = -gamma1 / beta1;
        eta[0] = (6.0 * u[0]) / beta1;
    }
    for (int i = 1; i < n_unknown - 1; ++i) {
        double alpha = h[i] / (h[i] + h[i + 1]);
        double beta = 2.0;
        double gamma = h[i + 1] / (h[i] + h[i + 1]);
        double denominator = beta + alpha * xi[i - 1];
        
        xi[i] = -gamma / denominator;
        
        eta[i] = (6.0 * u[i] - alpha * eta[i - 1]) / denominator;
    }

    if (n_unknown >= 2) {
        int i = n_unknown - 1;
        double alpha = h[N - 1] / (h[N - 2] + h[N - 1]);
        double beta = 2.0;
        double denominator = beta + alpha * xi[i - 1];
        eta[i] = (6.0 * u[i] - alpha * eta[i - 1]) / denominator;
    }
    if (n_unknown >= 1) {
        c[N - 1] = eta[n_unknown - 1];
        for (int i = n_unknown - 2; i >= 0; --i) {
            c[i + 1] = xi[i] * c[i + 2] + eta[i];
        }
    }
    c[0] = 0.0;
    c[N] = 0.0;
    
    return c;
}
