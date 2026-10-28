# Função de distribuição radial g(r) em 2D, acumulada separadamente para cada estado (kT).
#
# g(r) = (pares com distância em [r, r+dr]) / (pares esperados num gás ideal na mesma casca)
#      = 2 C(r) / [ M n (n-1)/L² π((r+dr)² - r²) ]
# C(r): contagem de pares somada nos M quadros; n-1 em vez de n deixa g = 1 exato
# para partículas sem correlação, mesmo com poucas partículas.

mutable struct HistGR
    dr::Float64               # largura do bin (mesma unidade de L)
    counts::Vector{Float64}   # pares por bin, somados em todos os quadros
    nframes::Int
end

# rmax <= L/2: além disso a casca sai da caixa e a imagem mínima não vê todos os pares
HistGR(rmax, dr) = HistGR(dr, zeros(max(1, floor(Int, rmax / dr))), 0)

function zera!(h::HistGR)
    fill!(h.counts, 0.0)
    h.nframes = 0
    return nothing
end

function acumula_gr!(h::HistGR, r::Replica, L)
    n = length(r.x); nb = length(h.counts)
    rmax2 = (nb * h.dr)^2
    @inbounds for i in 1:n-1
        xi = r.x[i]; yi = r.y[i]
        for j in i+1:n
            dx = xi - r.x[j]; dy = yi - r.y[j]
            dx -= L * round(dx / L); dy -= L * round(dy / L)   # imagem mínima
            d2 = dx * dx + dy * dy
            d2 < rmax2 || continue
            b = floor(Int, sqrt(d2) / h.dr) + 1
            b <= nb && (h.counts[b] += 1.0)
        end
    end
    h.nframes += 1
    return nothing
end

# pontos (r no centro do bin, g) para o gráfico
function gr_pontos(h::HistGR, n, L)
    h.nframes == 0 && return Point2f[]
    pares = n * (n - 1) / L^2            # n(n-1)/L²: densidade de pares do gás ideal
    pts = Point2f[]
    for b in eachindex(h.counts)
        casca = π * ((b * h.dr)^2 - ((b - 1) * h.dr)^2)
        g = 2 * h.counts[b] / (h.nframes * pares * casca)
        push!(pts, Point2f((b - 0.5) * h.dr, g))
    end
    return pts
end
