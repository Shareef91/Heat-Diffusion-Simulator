#include <iostream>
#include <vector>
#include <iomanip>
#include <omp.h>
#include <fstream>

const int N = 100; //Used to create the grid size (N x N)
const int STEPS = 1000; //Number of steps to iterate
const double alpha = 0.1; //Constant for Heat Diffusion
const double dt = 0.1; //step for time step
const double dx = 1.0; //distance between the grid cells

//discrete heat diffusion function used to update the temperature for each cell
double compute_new_value(double up, double down, double left, double right, double center) {
    return center + alpha * dt / (dx * dx) * (up + down + left + right - 4 * center);
}

int main() {
    std::vector<std::vector<double>> grid(N, std::vector<double>(N, 25.0)); //each grid is at 25 degrees C (77 degrees F)
    std::vector<std::vector<double>> new_grid = grid;
    grid[N/2][N/2] = 100.0; //Hot spot in the center of the plate 

    double start = omp_get_wtime();// start time for the parallel execution

    for (int step = 0; step < STEPS; ++step) { //main iteration loop
        #pragma omp parallel for collapse(2) //Use openmp to parallel both loops 
        for (int i = 1; i < N - 1; ++i) {
            for (int j = 1; j < N - 1; ++j) {
                new_grid[i][j] = compute_new_value( //Compute the new value
                    grid[i-1][j], grid[i+1][j], grid[i][j-1], grid[i][j+1], grid[i][j]
                );
            }
        }
        grid.swap(new_grid); //Swap old grid for new grid
    }

    double end = omp_get_wtime(); //Timer stopped
    std::cout << "Runtime: " << end - start << " seconds\n"; //print the runtime total

    // --- SAVE FINAL GRID TO CSV ---
    std::ofstream out("heat_output_openmp.csv");
    for (int i = 0; i < N; ++i) { 
        for (int j = 0; j < N; ++j) {
            if (j) out << ",";
            out << grid[i][j];
        }
        out << "\n"; // new line
    }
    out.close();
    std::cout << "Saved heat_output_openmp.csv (OpenMP)\n";
    return 0;
}

