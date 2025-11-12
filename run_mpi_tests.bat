@echo off
set EXE=heat_diffusion_mpi.exe
if not exist %EXE% (echo %EXE% not found & exit /b 1)
echo Procs,Runtime(s)
for %%P in (1 2 4 8) do (
  for /f "tokens=2 delims= " %%R in ('mpiexec -n %%P %EXE% ^| findstr /i "Runtime:"') do echo %%P,%%R
)
