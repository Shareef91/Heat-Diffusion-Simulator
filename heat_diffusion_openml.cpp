#include <iostream>
#include <vector>
#include <iomanip>
#include <omp.h>
#include <fstream>

const int N = 100;
const int STEPS = 1000;
const double alpha = 0.1;
const double dt = 0.1;
const double dx = 1.0;

double compute_new_value(double up, double down, double left, double right, double center) {
    return center + alpha * dt / (dx * dx) * (up + down + left + right - 4 * center);
}

int main() {
    std::vector<std::vector<double>> grid(N, std::vector<double>(N, 25.0));
    std::vector<std::vector<double>> new_grid = grid;
    grid[N/2][N/2] = 100.0;

    double start = omp_get_wtime();

    for (int step = 0; step < STEPS; ++step) {
        #pragma omp parallel for collapse(2)
        for (int i = 1; i < N - 1; ++i) {
            for (int j = 1; j < N - 1; ++j) {
                new_grid[i][j] = compute_new_value(
                    grid[i-1][j], grid[i+1][j], grid[i][j-1], grid[i][j+1], grid[i][j]
                );
            }
        }
        grid.swap(new_grid);
    }

    double end = omp_get_wtime();
    std::cout << "Runtime: " << end - start << " seconds\n";

    // --- SAVE FINAL GRID TO CSV ---
    std::ofstream out("heat_output_openmp.csv");
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (j) out << ",";
            out << grid[i][j];
        }
        out << "\n";
    }
    out.close();
    std::cout << "Saved heat_output_openmp.csv (OpenMP)\n";
    return 0;
}

