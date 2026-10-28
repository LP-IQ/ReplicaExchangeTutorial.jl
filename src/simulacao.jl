# Núcleo da simulação: fluido de Lennard-Jones 2D, Langevin BAOAB e trocas T-REMD.
# Os outros termostatos estão em termostatos.jl.
# Não depende da janela; é o que os testes em test/runtests.jl exercitam.
# Unidades reduzidas: massa = 1, energia em unidades de epsilon, temperatura como kT.

const MAXREP = 12

# uma configuração (walker); o campo walker é a identidade e vai junto nas trocas
mutable struct Replica
    x::Vector{Float64}; y::Vector{Float64}       # posições
    vx::Vector{Float64}; vy::Vector{Float64}     # velocidades
    fx::Vector{Float64}; fy::Vector{Float64}     # forças
    U::Float64                                   # energia potencial (unidades de epsilon)
    walker::Int                                  # identidade da configuração (cor)
end
Replica(n::Int, w::Int) = Replica(zeros(n), zeros(n), zeros(n), zeros(n), zeros(n), zeros(n), 0.0, w)

# reps[k] = configuração que está no estado k (temperatura kTs[k])
# trocar = permutar dois elementos de reps
Base.@kwdef mutable struct Sim
    reps::Vector{Replica}
    kTs::Vector{Float64}
    L::Float64 = 100.0
    eps::Float64 = 1.0
    sig::Float64 = 2.0
    dt::Float64 = 0.05
    gam::Float64 = 0.1
    step::Int = 0
    nattempt::Int = 0
    att::Vector{Int} = Int[]              # tentativas por par
    acc::Vector{Int} = Int[]              # trocas aceitas por par
    last_prob::Vector{Float64} = Float64[]  # P da última tentativa de cada par (NaN = nunca tentado)
    rngs::Vector{MersenneTwister} = MersenneTwister[]
    termostato::Symbol = :langevin         # um de TERMOSTATOS (termostatos.jl)
    vxi::Matrix{Float64} = zeros(0, 0)     # Nosé-Hoover: velocidades da cadeia, uma coluna por estado
    prepared::Bool = false
end

# LJ: U = 4eps[(sig/r)^12 - (sig/r)^6], todos os pares (O(n^2)), sem cutoff
function forces!(r::Replica, L, eps, sig)
    n = length(r.x); fill!(r.fx, 0.0); fill!(r.fy, 0.0); U = 0.0; s2 = sig^2
    @inbounds for i in 1:n-1
        xi = r.x[i]; yi = r.y[i]
        for j in i+1:n
            dx = xi - r.x[j]; dy = yi - r.y[j]
            # imagem mínima
            dx -= L * round(dx / L); dy -= L * round(dy / L)
            r2 = dx * dx + dy * dy
            # (sig/r)^6, sem sqrt
            s6 = (s2 / r2)^3
            U += 4 * eps * (s6 * s6 - s6)
            # |F|/r
            fm = 24 * eps * (2 * s6 * s6 - s6) / r2
            r.fx[i] += fm * dx; r.fy[i] += fm * dy
            r.fx[j] -= fm * dx; r.fy[j] -= fm * dy
        end
    end
    r.U = U
    return U
end

# anda cada partícula na direção da força, passo máx. 0.05 (só para tirar sobreposições)
function minimize!(r::Replica, L, eps, sig; niter=500)
    n = length(r.x)
    for _ in 1:niter
        forces!(r, L, eps, sig)
        @inbounds for i in 1:n
            f = sqrt(r.fx[i]^2 + r.fy[i]^2)
            f < 1e-12 && continue
            d = min(0.05, 0.01 * f)
            r.x[i] = mod(r.x[i] + d * r.fx[i] / f, L)
            r.y[i] = mod(r.y[i] + d * r.fy[i] / f, L)
        end
    end
    return forces!(r, L, eps, sig)
end

# Langevin BAOAB na temperatura kT
# c1 = exp(-gam*dt); c2 escolhido para manter <v^2> = kT (Leimkuhler & Matthews 2013)
function langevin_step!(r::Replica, kT, dt, gam, L, eps, sig, rng)
    n = length(r.x); c1 = exp(-gam * dt); c2 = sqrt(kT * (1 - c1^2)); h = 0.5 * dt
    @inbounds for i in 1:n
        # B
        r.vx[i] += h * r.fx[i]; r.vy[i] += h * r.fy[i]
        # A
        r.x[i] += h * r.vx[i]; r.y[i] += h * r.vy[i]
        # O: atrito + ruído
        r.vx[i] = c1 * r.vx[i] + c2 * randn(rng)
        r.vy[i] = c1 * r.vy[i] + c2 * randn(rng)
        # A
        r.x[i] = mod(r.x[i] + h * r.vx[i], L)
        r.y[i] = mod(r.y[i] + h * r.vy[i], L)
    end
    forces!(r, L, eps, sig)
    @inbounds for i in 1:n
        # B
        r.vx[i] += h * r.fx[i]; r.vy[i] += h * r.fy[i]
    end
    return nothing
end

# entre trocas as réplicas são independentes -> uma por thread
function step_all!(s::Sim)
    Threads.@threads for i in eachindex(s.reps)
        if s.termostato == :langevin
            langevin_step!(s.reps[i], s.kTs[i], s.dt, s.gam, s.L, s.eps, s.sig, s.rngs[i])
        else
            passo_vv!(s, i)
        end
    end
    s.step += 1
    return nothing
end

# P = min{1, exp[(β_i − β_j)(U_i − U_j)]}; pares pares e ímpares alternados
function attempt_exchange!(s::Sim, rng)
    R = length(s.reps); out = Tuple{Int,Bool,Float64}[]
    # nattempt par: pares 1-2, 3-4, ...; ímpar: 2-3, 4-5, ...
    for i in (1 + s.nattempt % 2):2:(R - 1)
        a = s.reps[i]; b = s.reps[i+1]
        delta = (1 / s.kTs[i] - 1 / s.kTs[i+1]) * (a.U - b.U)
        P = min(1.0, exp(delta))
        ok = rand(rng) < P
        s.att[i] += 1
        s.last_prob[i] = P
        if ok
            s.acc[i] += 1
            sc = sqrt(s.kTs[i+1] / s.kTs[i])       # v ∝ sqrt(T_novo/T_antigo)
            a.vx .*= sc; a.vy .*= sc
            b.vx ./= sc; b.vy ./= sc
            # troca de estado (as forças não mudam, só dependem das posições)
            s.reps[i], s.reps[i+1] = b, a
        end
        push!(out, (i, ok, P))
    end
    s.nattempt += 1
    return out
end

# mesma configuração inicial em todas; walker w começa no estado w
function new_sim(kTs, n, L, eps, sig, dt, gam; seed=rand(1:10^6), termostato=:langevin)
    termostato in TERMOSTATOS || error("termostato desconhecido: $termostato")
    R = length(kTs); rng = MersenneTwister(seed)
    x0 = L .* rand(rng, n); y0 = L .* rand(rng, n)
    reps = Replica[]
    for w in 1:R
        r = Replica(n, w); r.x .= x0; r.y .= y0; push!(reps, r)
    end
    return Sim(; reps, kTs=collect(Float64, kTs), L, eps, sig, dt, gam,
               att=zeros(Int, R - 1), acc=zeros(Int, R - 1), last_prob=fill(NaN, R - 1),
               rngs=[MersenneTwister(seed + 17 * i) for i in 1:R],
               termostato, vxi=zeros(NCADEIA, R))
end

# Minimiza uma configuração e copia para todas as réplicas.
# fator_kT: as velocidades iniciais são sorteadas em fator_kT × kT do estado
# (1 = já no alvo; 2 = começa quente, para ver o termostato trazer a temperatura de volta)
function prepare!(s::Sim; fator_kT=1.0)
    R = length(s.reps); r1 = s.reps[1]
    minimize!(r1, s.L, s.eps, s.sig)
    for i in 1:R
        r = s.reps[i]; kT = fator_kT * s.kTs[i]
        r.x .= r1.x; r.y .= r1.y
        # velocidades de Maxwell-Boltzmann (m = 1)
        r.vx .= sqrt(kT) .* randn(s.rngs[i], length(r.x))
        r.vy .= sqrt(kT) .* randn(s.rngs[i], length(r.x))
        if s.termostato != :langevin           # momento total zero: 2n - 2 graus de liberdade
            r.vx .-= sum(r.vx) / length(r.vx); r.vy .-= sum(r.vy) / length(r.vy)
        end
        forces!(r, s.L, s.eps, s.sig)
    end
    fill!(s.vxi, 0.0)
    s.prepared = true
    return nothing
end

# histograma normalizado -> pontos (centro do bin, densidade)
function hist_points(v::Vector{Float64}, lo, hi, nb)
    c = zeros(Int, nb); w = (hi - lo) / nb
    for u in v
        c[clamp(floor(Int, (u - lo) / w) + 1, 1, nb)] += 1
    end
    tot = max(1, length(v))
    return [Point2f(lo + (b - 0.5) * w, c[b] / (tot * w)) for b in 1:nb]
end
