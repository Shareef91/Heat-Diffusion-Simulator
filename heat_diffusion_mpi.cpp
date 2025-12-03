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

const int N = 100; // Grid size (N x N)
const int STEPS = 500; // Number of time steps
const double alpha = 0.1; // Diffusion cefficient
const double dt = 0.1; // Time step size
const double dx = 1.0; // Spatial step size

// Computes new temperature at a cell based on a 5 point stencil
double compute_new_value(double up, double down, double left, double right, double center) {
    return center + alpha * dt / (dx * dx) * (up + down + left + right - 4 * center);
}

int main(int argc, char** argv) {
#if USE_MPI
    // Initialize MPI environment
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Rank ID (0.. size-1)
    MPI_Comm_size(MPI_COMM_WORLD, &size); // Total number of processes
#else
    int rank = 0, size = 1; // Fallback single process (No MPI)
#endif

    int local_rows = N / size; // Each process handles a block of rows: N/size rows per process

    // Local grid includes 2 halo rows like top + bottom for boundary exchanges
    // Dimensions: (local_rows + 2) x N
    std::vector<std::vector<double>> local_grid(local_rows + 2, std::vector<double>(N, 25.0));
    std::vector<std::vector<double>> new_grid = local_grid; // Buffer updastes on values

    if (rank == size / 2) local_grid[local_rows / 2][N / 2] = 100.0; // heat source in middle node

#if USE_MPI
    // Synchronize all ranks before starting timing
    MPI_Barrier(MPI_COMM_WORLD);
    double t0 = MPI_Wtime(); // Start time (MPI wall clock)
#else
    auto t0_clock = std::chrono::steady_clock::now(); // Start time
#endif 
    // Main simulation loop
    for (int step = 0; step < STEPS; ++step) {
        // exchange boundary rows
#if USE_MPI
        if (rank > 0)
            MPI_Sendrecv(&local_grid[1][0],  // send buffer (first interior row)
                         N, MPI_DOUBLE, rank - 1, 0,  // send N doubles to rank-1
                         &local_grid[0][0], // receive buffer (top halo row)
                         N, MPI_DOUBLE, rank - 1, 0,   // receive N doubles from rank-1
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE);
         // Send our last owned row to the next rank; receive its first owned row into halo row local_rows+1
        if (rank < size - 1)
            MPI_Sendrecv(&local_grid[local_rows][0], 
                         N, MPI_DOUBLE, rank + 1, 0,
                         &local_grid[local_rows + 1][0],
                         N, MPI_DOUBLE, rank + 1, 0,
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE);
#endif
   // --- Update local interior cells using 5-point stencil ---
        for (int i = 1; i <= local_rows; ++i) {  // skip halo rows (1 .. local_rows)
            for (int j = 1; j < N - 1; ++j) {      // skip left/right boundaries (0 and N-1)
                new_grid[i][j] = compute_new_value(
                    local_grid[i - 1][j], // UP
                    local_grid[i + 1][j], // DOWN
                    local_grid[i][j - 1],  // LEFT
                    local_grid[i][j + 1],  // RIGHT
                    local_grid[i][j]   // CENTER
                );
            }
        }

        // Swaps new and old grids for next iteration
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

            // Gather all local chunkc to rank 0
            MPI_Gather(local_chunk.data(),  // send buffer (each rank)
                       local_rows * N, MPI_DOUBLE,  // number of doubles from each rank
                       rank == 0 ? full_grid.data() // receive buffer (only on rank 0)
                                 : nullptr,
                       local_rows * N, MPI_DOUBLE,    // receive count per rank
                       0, MPI_COMM_WORLD); // root = 0

             // Rank 0 writes full grid to CSV
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
            // Non-MPI fallback: just dump the local grid as if it were the full grid
            if (rank == 0) {
                std::string filename = std::string("heat_output_mpi_") + std::to_string(step) + ".csv";
                std::ofstream out(filename);
                // local_grid has halo rows: valid data is in local_grid[1..local_rows]
                for (int i = 0; i < N; ++i) {
                    for (int j = 0; j < N; ++j) {
                        if (j) out << ",";
                        out << local_grid[i + 1][j]; // +1 to skip halo
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
     // Stop timer and print runtime from rank 0
    MPI_Barrier(MPI_COMM_WORLD);
    double t1 = MPI_Wtime();
    if (rank == 0) {
        std::cout << "Runtime: " << (t1 - t0) << " seconds\n";
    }
    MPI_Finalize(); // Clean up MPI
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
