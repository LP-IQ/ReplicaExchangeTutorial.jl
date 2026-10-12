@echo off
rem Abre a GUI de REMD no Windows. Requer Julia no PATH e GLMakie instalado (ver README).
cd /d "%~dp0"
julia -t auto remd_gui.jl
pause
