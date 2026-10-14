#!/usr/bin/env julia
# Abre a GUI a partir do repositório clonado, sem instalar o pacote:
#     julia -t auto remd_gui.jl
# Com o pacote instalado (pkg> add ReplicaExchangeTutorial), basta:
#     using ReplicaExchangeTutorial; remd_gui()

using Pkg
Pkg.activate(@__DIR__; io=devnull)
Pkg.instantiate()   # na primeira vez baixa o GLMakie (demora alguns minutos)

using ReplicaExchangeTutorial
wait(display(remd_gui()))
