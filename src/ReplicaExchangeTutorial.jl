"""
    ReplicaExchangeTutorial

Tutorial interativo de Replica Exchange MD em temperatura (T-REMD) num fluido de
Lennard-Jones 2D. Abra a interface com `remd_gui()`.
"""
module ReplicaExchangeTutorial

using GLMakie, Random, Printf

export remd_gui

include("simulacao.jl")
include("interface.jl")

end
