@echo off
robocopy X:\apps\horses3d-gpu-miguel_dev_cpu C:\Users\mchav\UPMdrive\00_horses3d-gpu-miguel_dev\repo_local /MIR /XF *.o *.mod *.exe *.a *.so *.out *.ns *.mu *.ins *.nssa *.ch *.bin /XD .git
pause