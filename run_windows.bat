@echo off
rem Abre a GUI de REMD no Windows. Requer Julia no PATH; na primeira vez instala as dependencias.
cd /d "%~dp0"
julia -t auto remd_gui.jl
pause
