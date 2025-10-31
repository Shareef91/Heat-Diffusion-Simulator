## Repository snapshot

This is a tiny single-file C++ project implementing a 2D heat diffusion simulation in `heat_diffusion.cpp`.
Keep changes small and self-contained: the program is a single translation unit and is executed as a standalone binary.

## Big picture

- Entry point: `main()` in `heat_diffusion.cpp`.
- Core compute kernel: `compute_new_value(double up, double down, double left, double right, double center)` — this implements the finite-difference update for one cell.
- Configuration is encoded as top-level constants in the file: `N`, `STEPS`, `alpha`, `dt`, `dx`.
- Data layout: `std::vector<std::vector<double>> grid` and `new_grid`. The simulation updates inner cells `[1 .. N-2]` and leaves boundary cells at their initial values (ambient 25.0).

Why things are structured this way
- Single-file simplicity: the project is intended as a demonstration or small simulation; expect changes to either keep it simple or be refactored into multiple files if adding features/tests.
- Explicit constants: many parameters are compile-time constants inside the file — prefer editing these or adding CLI/config parsing if needed.

## Project-specific patterns and gotchas

- Boundary handling: edges are intentionally not updated (loop runs `i=1..N-2`, `j=1..N-2`). Do not change inner-loop ranges without checking boundary semantics.
- Hotspot initialization: `grid[N/2][N/2] = 100.0;` — central cell is seeded hot. Tests and modifications should be aware of this deterministic initial condition.
- Progress output: program prints only summary lines every 100 steps (`Step X: Center = ...`). Any change to logging frequency will change console output used for quick checks.
- Numeric stability: `alpha`, `dt`, and `dx` form the CFL-like factor `alpha*dt/(dx*dx)` used in `compute_new_value`. Keep this in mind when changing time step or spatial resolution.

## Build / run (assumptions)

Assumption: this is a single-file C++ app that can be built with a standard C++ toolchain (g++ / clang++ / MSVC). If your environment uses a different toolchain, tell me and I'll add exact commands.

Recommended examples (PowerShell on Windows):

- With g++ (MinGW / mingw-w64):

  g++ -std=c++17 -O2 heat_diffusion.cpp -o heat_diffusion.exe

- With MSVC Developer Command Prompt (cl):

  cl /EHsc /O2 heat_diffusion.cpp

- Run:

  .\heat_diffusion.exe

If the workspace later adds a build system (CMake, Makefile), prefer using that instead of the raw compiler commands.

## Where to look when changing behavior

- To change simulation parameters: edit the constants near the top of `heat_diffusion.cpp` (`N`, `STEPS`, `alpha`, `dt`, `dx`).
- To change initial conditions: inspect the section that initializes `grid` and the hotspot lines.
- To expose parameters to users: add minimal CLI parsing (e.g., `argc/argv`) or a thin `struct Config` and keep `main()` small.

## Examples to include in PRs

- When making algorithmic changes, include a small table showing a few step outputs (e.g., center temperature at steps 0,100,500,1000) so reviewers can quickly validate behavior.
- If you refactor into multiple files, keep the public kernel function signature comparable to `compute_new_value(...)` so it's easy to test in isolation.

## Tests & verification

- There are no tests present; for quick CI-friendly checks, add a unit that runs a tiny grid (e.g., N=5, STEPS=1) and asserts expected numeric relationships (e.g., center decreases/increases according to alpha). Keep tests deterministic by controlling RNG and initial conditions.

## Questions for maintainers (ask when unclear)

- Preferred compiler/toolchain on CI (MSVC vs GCC/Clang)?
- Expected CLI / input/output format if this binary will be integrated into larger tooling or datasets.

If anything here looks wrong or you'd like the instructions expanded (toolchain specifics, CI, tests), tell me which parts to clarify and I will update the file.
