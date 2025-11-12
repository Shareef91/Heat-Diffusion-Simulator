#include <iostream>
#include <vector>
#include <chrono>
#include <fstream>
#include <string>

#if defined(__has_include)
#  if __has_include(<mpi.h>)
#    define USE_MPI 1
#    include <mpi.h>
#  else
#    define USE_MPI 0
#  endif
#else
#  define USE_MPI 0
#endif

const int N = 100;
const int STEPS = 500;
const double alpha = 0.1;
const double dt = 0.1;
const double dx = 1.0;

double compute_new_value(double up, double down, double left, double right, double center) {
    return center + alpha * dt / (dx * dx) * (up + down + left + right - 4 * center);
}

int main(int argc, char** argv) {
#if USE_MPI
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
#else
    int rank = 0, size = 1;
#endif

    int local_rows = N / size;
    std::vector<std::vector<double>> local_grid(local_rows + 2, std::vector<double>(N, 25.0));
    std::vector<std::vector<double>> new_grid = local_grid;

    if (rank == size / 2) local_grid[local_rows / 2][N / 2] = 100.0; // heat source in middle node

#if USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
    double t0 = MPI_Wtime();
#else
    auto t0_clock = std::chrono::steady_clock::now();
#endif

    for (int step = 0; step < STEPS; ++step) {
        // exchange boundary rows
#if USE_MPI
        if (rank > 0)
            MPI_Sendrecv(&local_grid[1][0], N, MPI_DOUBLE, rank - 1, 0,
                         &local_grid[0][0], N, MPI_DOUBLE, rank - 1, 0,
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        if (rank < size - 1)
            MPI_Sendrecv(&local_grid[local_rows][0], N, MPI_DOUBLE, rank + 1, 0,
                         &local_grid[local_rows + 1][0], N, MPI_DOUBLE, rank + 1, 0,
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE);
#endif

        for (int i = 1; i <= local_rows; ++i) {
            for (int j = 1; j < N - 1; ++j) {
                new_grid[i][j] = compute_new_value(
                    local_grid[i - 1][j], local_grid[i + 1][j],
                    local_grid[i][j - 1], local_grid[i][j + 1], local_grid[i][j]
                );
            }
        }
        local_grid.swap(new_grid);

        // periodically gather full grid on rank 0 and write to CSV
        if (step % 100 == 0) {
#if USE_MPI
            // pack local (non-halo) rows into a contiguous buffer
            std::vector<double> local_chunk(local_rows * N);
            for (int ii = 0; ii < local_rows; ++ii) {
                for (int jj = 0; jj < N; ++jj) {
                    local_chunk[ii * N + jj] = local_grid[ii + 1][jj];
                }
            }

            // root allocates full buffer to receive all pieces
            std::vector<double> full_grid;
            if (rank == 0) full_grid.resize(N * N);

            MPI_Gather(local_chunk.data(), local_rows * N, MPI_DOUBLE,
                       rank == 0 ? full_grid.data() : nullptr, local_rows * N, MPI_DOUBLE,
                       0, MPI_COMM_WORLD);

            if (rank == 0) {
                std::string filename = std::string("heat_output_mpi_") + std::to_string(step) + ".csv";
                std::ofstream out(filename);
                for (int i = 0; i < N; ++i) {
                    for (int j = 0; j < N; ++j) {
                        if (j) out << ",";
                        out << full_grid[i * N + j];
                    }
                    out << "\n";
                }
                out.close();
                std::cout << "Saved heat_output_mpi_" << step << ".csv (MPI Rank 0)\n";
            }
#else
            if (rank == 0) {
                std::string filename = std::string("heat_output_mpi_") + std::to_string(step) + ".csv";
                std::ofstream out(filename);
                // local_grid has halo rows: valid data is in local_grid[1..local_rows]
                for (int i = 0; i < N; ++i) {
                    for (int j = 0; j < N; ++j) {
                        if (j) out << ",";
                        out << local_grid[i + 1][j];
                    }
                    out << "\n";
                }
                out.close();
                std::cout << "Saved heat_output_mpi_" << step << ".csv\n";
            }
#endif
        }
    }

#if USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
    double t1 = MPI_Wtime();
    if (rank == 0) {
        std::cout << "Runtime: " << (t1 - t0) << " seconds\n";
    }
    MPI_Finalize();
#else
    auto t1_clock = std::chrono::steady_clock::now();
    double elapsed = std::chrono::duration<double>(t1_clock - t0_clock).count();
    std::cout << "Runtime: " << elapsed << " seconds\n";
#endif

    // --- SAVE FINAL GRID TO CSV (only rank 0) ---
    if (rank == 0) {
        std::ofstream out("heat_output_mpi.csv");
        // local_grid has halo rows: valid data is in local_grid[1..local_rows]
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                if (j) out << ",";
                out << local_grid[i + 1][j];
            }
            out << "\n";
        }
        out.close();
        std::cout << "Saved heat_output_mpi.csv (MPI Rank 0)\n";
    }

    return 0;
}
