#include <iostream>
#include <vector>
#include <iomanip>
#include <fstream>
#include <string>
#include <fstream>

const int N = 100;        // Grid size (N x N)
const int STEPS = 1000;   // Number of time steps
const double alpha = 0.1; // Heat transfer rate
const double dt = 0.1;    // Time increment
const double dx = 1.0;    // Distance between cells

double compute_new_value(double up, double down, double left, double right, double center) {
    return center + alpha * dt / (dx * dx) * (up + down + left + right - 4 * center);
}

int main() {
    std::vector<std::vector<double>> grid(N, std::vector<double>(N, 25.0));
    std::vector<std::vector<double>> new_grid = grid;

    // Start with a hot spot in the center
    grid[N/2][N/2] = 100.0;

    for (int step = 0; step < STEPS; ++step) {
        for (int i = 1; i < N - 1; ++i) {
            for (int j = 1; j < N - 1; ++j) {
                new_grid[i][j] = compute_new_value(
                    grid[i-1][j], grid[i+1][j], grid[i][j-1], grid[i][j+1], grid[i][j]
                );
            }
        }
        grid.swap(new_grid);

        // Show progress every 100 steps and write CSV dump
        if (step % 100 == 0) {
            std::cout << "Step " << step << ": Center = " << grid[N/2][N/2] << "\n";
            std::string filename = std::string("heat_output_seq_") + std::to_string(step) + ".csv";
            std::ofstream out(filename);
            for (int i = 0; i < N; ++i) {
                for (int j = 0; j < N; ++j) {
                    if (j) out << ",";
                    out << grid[i][j];
                }
                out << "\n";
            }
            out.close();
        }
    }

    // --- SAVE FINAL GRID TO CSV ---
    std::ofstream out("heat_output_seq.csv");
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (j) out << ",";
            out << grid[i][j];
        }
        out << "\n";
    }
    out.close();
    std::cout << "Saved heat_output_seq.csv (Sequential)\n";

    return 0;
}

