# Termostatos além do Langevin. Todos integram a dinâmica com velocity Verlet e usam
# τ = 1/gamma como tempo de acoplamento (mesma unidade de tempo de dt).
#
#   :bussi       reescalonamento estocástico (Bussi, Donadio e Parrinello, 2007): canônico
#   :berendsen   reescalonamento determinístico para o alvo (Berendsen et al., 1984): NÃO canônico,
#                as flutuações da energia cinética ficam menores que as do ensemble canônico
#   :nosehoover  cadeia de Nosé-Hoover com 3 termostatos (Martyna, Klein e Tuckerman, 1992):
#                canônico se a dinâmica for ergódica
#   :nve         sem termostato: a energia total só muda nas trocas

const TERMOSTATOS = (:langevin, :bussi, :berendsen, :nosehoover, :nve)
const NCADEIA = 3   # termostatos na cadeia de Nosé-Hoover

energia_cinetica(r::Replica) = 0.5 * (sum(abs2, r.vx) + sum(abs2, r.vy))   # massa = 1

# Langevin troca momento com o banho: 2n graus de liberdade.
# Os outros conservam o momento total, que é zerado no início: 2n - 2.
graus_liberdade(n, termostato) = termostato == :langevin ? 2n : 2n - 2

function escala_v!(r::Replica, a)
    r.vx .*= a; r.vy .*= a
    return nothing
end

# velocity Verlet: meio passo em v, passo inteiro em x, forças novas, meio passo em v
function vv_step!(r::Replica, dt, L, eps, sig)
    n = length(r.x); h = 0.5 * dt
    @inbounds for i in 1:n
        r.vx[i] += h * r.fx[i]; r.vy[i] += h * r.fy[i]
        r.x[i] = mod(r.x[i] + dt * r.vx[i], L)
        r.y[i] = mod(r.y[i] + dt * r.vy[i], L)
    end
    forces!(r, L, eps, sig)
    @inbounds for i in 1:n
        r.vx[i] += h * r.fx[i]; r.vy[i] += h * r.fy[i]
    end
    return nothing
end

# Berendsen: v *= λ, λ² = 1 + (dt/τ)(kT/kT_inst - 1)
# puxa a energia cinética para o alvo, mas não gera as flutuações canônicas
function berendsen!(r::Replica, kT, dt, tau, nf)
    kTi = 2 * energia_cinetica(r) / nf
    kTi > 0 || return nothing
    escala_v!(r, sqrt(max(0.0, 1 + dt / tau * (kT / kTi - 1))))
    return nothing
end

# Bussi (CSVR): sorteia a nova energia cinética de modo que P(K) seja a canônica
# K' = K + (1-c)[K̄ (R1² + soma de nf-1 normais²)/nf - K] + 2 R1 sqrt(c(1-c) K K̄/nf),  c = exp(-dt/τ)
function bussi!(r::Replica, kT, dt, tau, nf, rng)
    K = energia_cinetica(r)
    K > 0 || return nothing
    Kalvo = 0.5 * nf * kT                      # K̄ = nf kT/2
    c = exp(-dt / tau)
    R1 = randn(rng)
    soma = 0.0
    for _ in 2:nf
        soma += randn(rng)^2
    end
    Knovo = K + (1 - c) * (Kalvo * (soma + R1^2) / nf - K) + 2 * R1 * sqrt(c * (1 - c) * K * Kalvo / nf)
    escala_v!(r, sqrt(max(0.0, Knovo) / K))
    return nothing
end

# "força" sobre a velocidade do termostato j da cadeia
# j = 1: (2K - nf kT)/Q1;  j > 1: (Q_{j-1} ξ_{j-1}² - kT)/Q_j
@inline function forca_nhc(j, K2, xi, nf, kT, Q1, Qj)
    j == 1 && return (K2 - nf * kT) / Q1
    Qant = j == 2 ? Q1 : Qj
    return (Qant * xi[j-1]^2 - kT) / Qj
end

# meio passo (dt/2) da cadeia de Nosé-Hoover (esquema de Martyna et al., 1996)
# xi: velocidades dos termostatos do ESTADO (não vão junto com a configuração nas trocas).
# Com Q1 = nf kT τ² e Qj = kT τ², a variância de ξ_j, kT/Q_j, não depende de kT; por isso
# ξ não precisa ser ajustado quando chega uma configuração de outra temperatura.
function nhc_meio_passo!(r::Replica, xi::AbstractVector, kT, dt, tau, nf)
    M = length(xi)
    Q1 = nf * kT * tau^2; Qj = kT * tau^2
    dt2 = dt / 2; dt4 = dt / 4; dt8 = dt / 8
    K2 = 2 * energia_cinetica(r)
    xi[M] += forca_nhc(M, K2, xi, nf, kT, Q1, Qj) * dt4
    for j in M-1:-1:1
        xi[j] *= exp(-xi[j+1] * dt8)
        xi[j] += forca_nhc(j, K2, xi, nf, kT, Q1, Qj) * dt4
        xi[j] *= exp(-xi[j+1] * dt8)
    end
    s = exp(-xi[1] * dt2)
    escala_v!(r, s); K2 *= s^2
    for j in 1:M-1
        xi[j] *= exp(-xi[j+1] * dt8)
        xi[j] += forca_nhc(j, K2, xi, nf, kT, Q1, Qj) * dt4
        xi[j] *= exp(-xi[j+1] * dt8)
    end
    xi[M] += forca_nhc(M, K2, xi, nf, kT, Q1, Qj) * dt4
    return nothing
end

# um passo da réplica que está no estado i, com o termostato da simulação
function passo_vv!(s, i)
    r = s.reps[i]; kT = s.kTs[i]; tau = 1 / s.gam
    nf = graus_liberdade(length(r.x), s.termostato)
    nh = s.termostato == :nosehoover
    nh && nhc_meio_passo!(r, view(s.vxi, :, i), kT, s.dt, tau, nf)
    vv_step!(r, s.dt, s.L, s.eps, s.sig)
    if nh
        nhc_meio_passo!(r, view(s.vxi, :, i), kT, s.dt, tau, nf)
    elseif s.termostato == :berendsen
        berendsen!(r, kT, s.dt, tau, nf)
    elseif s.termostato == :bussi
        bussi!(r, kT, s.dt, tau, nf, s.rngs[i])
    end
    return nothing
end

# kT instantâneo = 2K/nf
kT_inst(r::Replica, termostato) = 2 * energia_cinetica(r) / graus_liberdade(length(r.x), termostato)

# densidade canônica do kT instantâneo: com a = nf/2, K tem distribuição gama(a, kT) e
# p(T) = a^a T^(a-1) exp(-aT/kT) / [Γ(a) kT^a]; nf é par, então a é inteiro
function kT_inst_teorico(T, kT, nf)
    a = nf ÷ 2
    T > 0 || return 0.0
    lgam = sum(log(j) for j in 1:a-1; init=0.0)     # log Γ(a) = log (a-1)!
    return exp(a * log(a) + (a - 1) * log(T) - a * T / kT - lgam - a * log(kT))
end
