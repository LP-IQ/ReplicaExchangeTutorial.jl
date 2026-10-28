using ReplicaExchangeTutorial, Test, Random
using ReplicaExchangeTutorial: Replica, new_sim, prepare!, step_all!, attempt_exchange!,
                               forces!, langevin_step!, HistGR, acumula_gr!, gr_pontos,
                               energia_cinetica, graus_liberdade, kT_inst_teorico

# escada geométrica, como o botão da GUI
escada(a, b, R) = [a * (b / a)^((i - 1) / (R - 1)) for i in 1:R]

@testset "forças = -dU/dx" begin
    s = new_sim([0.6], 20, 30.0, 1.0, 2.0, 0.05, 0.1; seed=1)
    prepare!(s)
    r = s.reps[1]; h = 1e-6
    for i in (1, 7, 20)
        x0 = r.x[i]
        r.x[i] = x0 + h; Up = forces!(r, s.L, s.eps, s.sig)
        r.x[i] = x0 - h; Um = forces!(r, s.L, s.eps, s.sig)
        r.x[i] = x0;     forces!(r, s.L, s.eps, s.sig)
        @test r.fx[i] ≈ -(Up - Um) / (2h) rtol=1e-4 atol=1e-6
    end
end

@testset "termostato: <v²> = kT (gás quase ideal)" begin
    kT = 0.7; n = 10; L = 1.0e4             # caixa enorme: interações desprezíveis
    r = Replica(n, 1); rng = MersenneTwister(3)
    r.x .= L .* rand(rng, n); r.y .= L .* rand(rng, n)
    forces!(r, L, 1.0, 2.0)
    soma = 0.0; nam = 0
    for passo in 1:20_000
        langevin_step!(r, kT, 0.05, 1.0, L, 1.0, 2.0, rng)
        passo > 1000 || continue
        soma += sum(abs2, r.vx) + sum(abs2, r.vy); nam += 2n
    end
    @test soma / nam ≈ kT rtol=0.05
    @test all(0 .<= r.x .< L) && all(0 .<= r.y .< L)   # condição periódica
end

@testset "simulação curta" begin
    kTs = escada(0.55, 0.86, 6)
    s = new_sim(kTs, 30, 40.0, 1.0, 2.0, 0.05, 0.1; seed=42)
    @test [r.walker for r in s.reps] == 1:6
    prepare!(s)
    @test s.prepared && all(isfinite(r.U) for r in s.reps)
    for _ in 1:20
        for _ in 1:10; step_all!(s); end
        for (i, ok, P) in attempt_exchange!(s, MersenneTwister(s.step))
            @test 0 <= P <= 1
        end
    end
    @test s.step == 200
    @test sort([r.walker for r in s.reps]) == 1:6      # trocas só permutam as configurações
    @test all(s.acc .<= s.att)
    @test sum(s.att) == 10 * 3 + 10 * 2                 # 10 rodadas com 3 pares, 10 com 2
end

@testset "critério de troca" begin
    s = new_sim([0.5, 1.0], 5, 20.0, 1.0, 2.0, 0.05, 0.1; seed=7)
    prepare!(s)
    a, b = s.reps
    # configuração fria com U maior que a quente: delta > 0, troca sempre aceita
    a.U = 0.0; b.U = -10.0
    va = copy(a.vx)
    (i, ok, P), = attempt_exchange!(s, MersenneTwister(1))
    @test i == 1 && ok && P == 1.0
    @test s.reps[1] === b && s.reps[2] === a
    @test a.vx ≈ va .* sqrt(1.0 / 0.5)                # v ∝ sqrt(T_novo/T_antigo)
    # caso contrário: P = exp[(1/kT_i - 1/kT_j)(U_i - U_j)]
    s.nattempt = 0                                     # volta para o par 1-2
    s.reps[1].U = -10.0; s.reps[2].U = 0.0
    (_, _, P), = attempt_exchange!(s, MersenneTwister(1))
    @test P ≈ exp((1 / 0.5 - 1 / 1.0) * (-10.0))
end

@testset "g(r) em 2D" begin
    # partículas sem correlação: g(r) = 1 em todas as distâncias
    n = 100; L = 100.0; rng = MersenneTwister(5)
    h = HistGR(10.0, 0.05)
    r = Replica(n, 1)
    for _ in 1:400
        r.x .= L .* rand(rng, n); r.y .= L .* rand(rng, n)
        acumula_gr!(h, r, L)
    end
    pts = gr_pontos(h, n, L)
    g = [p[2] for p in pts if p[1] > 1]
    @test sum(g) / length(g) ≈ 1 rtol=0.02
    # um par a 3,0 de distância, atravessando a borda periódica: um único bin preenchido
    h2 = HistGR(10.0, 0.05); r2 = Replica(2, 1)
    r2.x .= [0.5, L - 2.5]; r2.y .= [10.0, 10.0]
    acumula_gr!(h2, r2, L)
    @test sum(h2.counts) == 1
    @test findfirst(>(0), h2.counts) == floor(Int, 3.0 / 0.05) + 1
    @test length(HistGR(60.0, 0.05).counts) * 0.05 <= 60.0
end

# um estado só (kT = 0,7), n = 40 numa caixa 40 x 40: denso o bastante para haver colisões
function corrida_termostato(termostato; nsteps=60_000, kT=0.7)
    s = new_sim([kT], 40, 40.0, 1.0, 2.0, 0.01, 1.0; seed=3, termostato)   # τ = 1/gamma = 1
    prepare!(s)
    r = s.reps[1]
    Ks = Float64[]
    for passo in 1:nsteps
        step_all!(s)
        passo > nsteps ÷ 5 && passo % 5 == 0 && push!(Ks, energia_cinetica(s.reps[1]))
    end
    nf = graus_liberdade(40, termostato)
    media = sum(Ks) / length(Ks)
    desvio = sqrt(sum(abs2, Ks .- media) / length(Ks))
    return 2 * media / nf, desvio / media, sqrt(2 / nf)   # <kT inst>, desvio relativo de K, valor canônico
end

@testset "termostatos" begin
    for t in (:langevin, :bussi, :nosehoover)
        kTm, rel, rel_canonico = corrida_termostato(t)
        @test kTm ≈ 0.7 rtol=0.05
        @test rel ≈ rel_canonico rtol=0.15            # flutuações de K do ensemble canônico
    end
    kTm, rel, rel_canonico = corrida_termostato(:berendsen)
    @test kTm ≈ 0.7 rtol=0.05
    @test rel < 0.5 * rel_canonico                    # Berendsen: flutuações pequenas demais

    # NVE: energia total e momento total conservados (velocity Verlet)
    s = new_sim([0.7], 40, 40.0, 1.0, 2.0, 0.01, 1.0; seed=3, termostato=:nve)
    prepare!(s)
    r = s.reps[1]
    E1 = r.U + energia_cinetica(r)
    for _ in 1:20_000
        step_all!(s)
    end
    r = s.reps[1]
    @test abs(r.U + energia_cinetica(r) - E1) < 0.005 * 40
    @test abs(sum(r.vx)) < 1e-8 && abs(sum(r.vy)) < 1e-8

    # a distribuição canônica do kT instantâneo integra a 1
    dT = 1e-3
    @test sum(kT_inst_teorico(T, 0.7, 78) for T in dT:dT:3.0) * dT ≈ 1 rtol=1e-3
    @test_throws ErrorException new_sim([0.7], 10, 20.0, 1.0, 2.0, 0.01, 1.0; termostato=:andersen)
end
