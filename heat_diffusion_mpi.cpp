#include <mpi.h>
#include <iostream>
#include <vector>

const int N = 100;
const int STEPS = 500;
const double alpha = 0.1;
const double dt = 0.1;
const double dx = 1.0;

double compute_new_value(double up, double down, double left, double right, double center) {
    return center + alpha * dt / (dx * dx) * (up + down + left + right - 4 * center);
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int local_rows = N / size;
    std::vector<std::vector<double>> local_grid(local_rows + 2, std::vector<double>(N, 25.0));
    std::vector<std::vector<double>> new_grid = local_grid;

    if (rank == size / 2) local_grid[local_rows / 2][N / 2] = 100.0; // heat source in middle node

    for (int step = 0; step < STEPS; ++step) {
        // exchange boundary rows
        if (rank > 0)
            MPI_Sendrecv(&local_grid[1][0], N, MPI_DOUBLE, rank - 1, 0,
                         &local_grid[0][0], N, MPI_DOUBLE, rank - 1, 0,
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        if (rank < size - 1)
            MPI_Sendrecv(&local_grid[local_rows][0], N, MPI_DOUBLE, rank + 1, 0,
                         &local_grid[local_rows + 1][0], N, MPI_DOUBLE, rank + 1, 0,
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        for (int i = 1; i <= local_rows; ++i) {
            for (int j = 1; j < N - 1; ++j) {
                new_grid[i][j] = compute_new_value(
                    local_grid[i - 1][j], local_grid[i + 1][j],
                    local_grid[i][j - 1], local_grid[i][j + 1], local_grid[i][j]
                );
            }
        }
        local_grid.swap(new_grid);
    }

    MPI_Finalize();
    return 0;
}
