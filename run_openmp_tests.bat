@echo off
setlocal enabledelayedexpansion
set EXE=heat_diffusion_openml.exe
if not exist %EXE% (
  echo %EXE% not found! Compile it first.
  exit /b 1
)
echo Threads,Runtime(s)
for %%T in (1 2 4 8) do (
  set OMP_NUM_THREADS=%%T
  for /f "tokens=2 delims= " %%R in ('%EXE% ^| findstr /i "Runtime:"') do (
    echo %%T,%%R
  )
)
endlocal
