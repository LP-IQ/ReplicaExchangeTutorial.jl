# Interface gráfica (GLMakie). Cada painel da direita é um ESTADO (kT fixo); a cor é a
# identidade da CONFIGURAÇÃO (walker). Troca aceita = duas configurações trocam de painel.

tbtext(tb) = (d = tb.displayed_string[]; d === nothing ? "" : strip(String(d)))
getf(tb) = parse(Float64, tbtext(tb))
geti(tb) = round(Int, getf(tb))
# Testes do tutorial (docs/tutorial.md, seção 9): nome no menu => valores dos campos.
# Campos omitidos voltam ao padrão (PADRAO). A lista kTs é gerada pela escada geométrica,
# exceto quando "kts" é dado explicitamente.
const PADRAO = Dict("nrep" => "6", "kmin" => "0.55", "kmax" => "0.86", "n" => "100", "L" => "100.0",
                    "dt" => "0.05", "nsteps" => "20000", "nex" => "50", "pause" => "0", "nequil" => "2000",
                    "gam" => "0.1", "eps" => "1.0", "sig" => "2.0", "seed" => "42", "termo" => "Langevin",
                    "fatorkT" => "1.0")
const TESTES = [
    "(escolher)"          => Dict{String,String}(),
    "1: referência"       => Dict{String,String}(),
    "2a: 3 réplicas"      => Dict("nrep" => "3"),
    "2b: 2 réplicas"      => Dict("nrep" => "2"),
    "3: larga, 3 rép."    => Dict("nrep" => "3", "kmin" => "0.35", "kmax" => "1.2", "nequil" => "5000"),
    "4a: larga, 8 rép."   => Dict("nrep" => "8", "kmin" => "0.35", "kmax" => "1.2", "nequil" => "5000"),
    "4b: larga, 10 rép."  => Dict("nrep" => "10", "kmin" => "0.35", "kmax" => "1.2", "nequil" => "5000"),
    "6a: troca a cada 10" => Dict("nex" => "10"),
    "6b: troca a cada 50" => Dict("nex" => "50"),
    "6c: troca a cada 500"=> Dict("nex" => "500"),
    "7a: n = 200"         => Dict("n" => "200", "L" => "141.0"),
    "7b: n = 400"         => Dict("n" => "400", "L" => "200.0"),
    "8: gargalo"          => Dict("kts" => "0.55,0.58,0.61,0.80,0.83,0.86"),
    "9: faixa fria"       => Dict("kmin" => "0.35", "kmax" => "0.45", "nequil" => "5000"),
    "10a: Bussi"          => Dict("termo" => "Bussi (v-rescale)"),
    "10b: Berendsen"      => Dict("termo" => "Berendsen"),
    "10c: Nosé-Hoover"    => Dict("termo" => "Nosé-Hoover"),
    "10d: sem termostato" => Dict("termo" => "Nenhum (NVE)"),
]
# nome no menu => termostato da simulação (termostatos.jl)
const OPCOES_TERMO = ["Langevin" => :langevin, "Bussi (v-rescale)" => :bussi, "Berendsen" => :berendsen,
                      "Nosé-Hoover" => :nosehoover, "Nenhum (NVE)" => :nve]

function settext!(tb, s::String)
    tb.displayed_string[] = s; tb.stored_string[] = s
    return nothing
end

const PALETTES = Dict(
    "Destaque" => [RGBf(0.9, 0.1, 0.1); fill(RGBf(0.65, 0.65, 0.65), MAXREP - 1)],
    "Colorido" => [
        RGBf(0.95, 0.1, 0.1), RGBf(0.1, 0.5, 0.8), RGBf(0.2, 0.7, 0.2), RGBf(0.9, 0.6, 0.1),
        RGBf(0.6, 0.3, 0.8), RGBf(0.2, 0.8, 0.8), RGBf(0.9, 0.3, 0.6), RGBf(0.4, 0.4, 0.4),
        RGBf(0.7, 0.8, 0.2), RGBf(0.6, 0.4, 0.2), RGBf(0.1, 0.4, 0.4), RGBf(0.8, 0.8, 0.3)],
)
# cor por temperatura: mapa 'turbo' (azul = menor kT -> verde -> amarelo -> vermelho = maior kT)
const TURBO = GLMakie.Makie.to_colormap(:turbo)
function tempcolor(k, R)
    t = R > 1 ? (k - 1) / (R - 1) : 0.0
    c = TURBO[clamp(round(Int, 1 + (0.08 + 0.84 * t) * (length(TURBO) - 1)), 1, length(TURBO))]
    return RGBf(c.r, c.g, c.b)
end

"""
    remd_gui()

Abre a interface do tutorial de T-REMD e devolve a `Figure`.

No REPL a janela abre sozinha. Num script, use `wait(display(remd_gui()))` para que o
Julia espere a janela ser fechada. Inicie o Julia com `julia -t auto`: cada réplica roda
numa thread. Figuras e trajetórias são salvas a partir do diretório atual (`pwd()`).

Tutorial: docs/tutorial.md no repositório do pacote.
"""
# tudo da interface fica aqui dentro; as funções internas compartilham as variáveis abaixo
function remd_gui()   # se ainda cortar: redimensione a janela ou edite o tamanho em Figure(size=...)
    Threads.nthreads() == 1 && @warn "Julia com 1 thread: as réplicas vão rodar em série. Inicie com  julia -t auto"
    monitor = GLMakie.GLFW.GetPrimaryMonitor()
    vidmode = GLMakie.GLFW.GetVideoMode(monitor)
    # os gráficos preenchem a janela e acompanham o redimensionamento; altura mínima útil ~650 px
    fig = Figure(size=(vidmode.width * 0.92, vidmode.height * 0.80), fontsize=12)

    left = fig[1, 1] = GridLayout(tellheight=false, valign=:top)   # não fixa a altura da janela
    center = fig[1, 2] = GridLayout()
    right_grid = fig[1, 3] = GridLayout()
    colsize!(fig.layout, 1, Fixed(300))
    colsize!(fig.layout, 3, Relative(0.24))
    # a linha principal ocupa toda a altura disponível, sem perguntar a altura aos painéis.
    # Sem isso, com poucas réplicas (uma linha de painéis à direita) a altura da janela passava
    # a seguir a altura desses painéis, e tudo mudava de tamanho a cada troca de título.
    rowsize!(fig.layout, 1, Auto(false))

    # ---------------- controles ----------------
    # telas baixas (1080p ou menos): controles mais compactos, para a coluna caber na janela
    compacto = vidmode.height < 1200
    fs = compacto ? 12.0 : 13.5                 # fonte dos rótulos
    ht, hm, hb = compacto ? (21, 22, 24) : (25, 27, 29)   # alturas: caixa de texto, menu, botão
    gap = compacto ? 2 : 4
    row = Ref(0)
    function field(label, default)
        row[] += 1
        Label(left[row[], 1], label; halign=:right, fontsize=fs)
        return Textbox(left[row[], 2]; stored_string=default, width=125, height=ht, fontsize=fs, textpadding=(7, 7, 3, 3))
    end
    tb_nrep = field("réplicas (máx $MAXREP)", "6")
    tb_kmin = field("kT mín (em epsilon)", "0.55")
    tb_kmax = field("kT máx (em epsilon)", "0.86")
    row[] += 1
    btn_ladder = Button(left[row[], 1:2]; label="Gerar escada geométrica", height=hm, fontsize=fs)
    tb_kts = field("kTs", "0.550,0.601,0.658,0.719,0.786,0.860")
    tb_n = field("n partículas", "100")
    tb_L = field("L (caixa)", "100.0")
    tb_dt = field("dt", "0.05")
    tb_nsteps = field("nsteps", "20000")
    tb_nex = field("troca a cada (passos)", "50")
    tb_pause = field("pausa se trocou (s)", "0.3")
    tb_nequil = field("equilibração", "2000")
    tb_fkT = field("kT inicial (× alvo)", "1.0")
    tb_gam = field("gamma (1/τ do termostato)", "0.1")
    row[] += 1
    Label(left[row[], 1], "termostato"; halign=:right, fontsize=fs)
    menu_termo = Menu(left[row[], 2]; options=first.(OPCOES_TERMO), default="Langevin", height=hm, fontsize=fs, width=145)
    tb_eps = field("epsilon (poço LJ)", "1.0")
    tb_sig = field("sigma (diâmetro LJ)", "2.0")
    tb_seed = field("semente (0 = aleatória)", "0")
    row[] += 1
    Label(left[row[], 1], "Teste"; halign=:right, fontsize=fs)
    menu_teste = Menu(left[row[], 2]; options=first.(TESTES), default="(escolher)", height=hm, fontsize=fs, width=145)
    row[] += 1
    Label(left[row[], 1], "Cores"; halign=:right, fontsize=fs)
    menu_cores = Menu(left[row[], 2]; options=["Destaque", "Colorido"], default="Destaque", height=hm, fontsize=fs, width=145)
    row[] += 1
    Label(left[row[], 1], "Mapa"; halign=:right, fontsize=fs)
    menu_map = Menu(left[row[], 2]; options=["Ocupação", "Menor energia", "Taxa de troca (kT)", "Trocas entre walkers"], default="Ocupação", height=hm, fontsize=fs, width=145)
    row[] += 1
    Label(left[row[], 1], "Gráfico"; halign=:right, fontsize=fs)
    menu_graf = Menu(left[row[], 2]; options=["P(U)", "g(r)", "kT instantâneo", "kT no tempo"], default="P(U)", height=hm, fontsize=fs, width=145)
    row[] += 1
    Label(left[row[], 1], "menor U no random walk"; halign=:right, fontsize=fs)
    tg_min = Toggle(left[row[], 2]; active=false, halign=:left)
    row[] += 1
    bl = left[row[], 1:2] = GridLayout()
    bk = (; height=hb, fontsize=fs)
    btn_min = Button(bl[1, 1:3]; label="1. Minimizar + Equilibrar", buttoncolor=RGBf(0.55, 0.80, 0.55), bk...)
    btn_run = Button(bl[2, 1]; label="2. Run", buttoncolor=RGBf(1.0, 0.85, 0.2), bk...)
    btn_stop = Button(bl[2, 2]; label="Stop", buttoncolor=RGBf(0.9, 0.45, 0.45), bk...)
    btn_reset = Button(bl[2, 3]; label="Reset", buttoncolor=RGBf(0.85, 0.75, 0.9), bk...)
    btn_save = Button(bl[3, 1:3]; label="Salvar figuras (PNG por gráfico)", buttoncolor=RGBf(0.75, 0.85, 0.95), bk...)
    btn_traj = Button(bl[4, 1:3]; label="Salvar trajetória do menor kT (XYZ)", buttoncolor=RGBf(0.75, 0.85, 0.95), bk...)
    rowgap!(bl, 4); colgap!(bl, 4)
    row[] += 1
    stats = Observable("Taxa de troca: —")
    Label(left[row[], 1:2], stats; halign=:left, valign=:top, justification=:left, tellwidth=false, tellheight=false, fontsize=fs - 1)   # última linha absorve a sobra
    rowgap!(left, gap)

    # ---------------- centro: random walk, energia, histogramas, réplica de referência e mapa ----------------
    ax_W = Axis(center[1, 1:2]; title="Random walk das configurações", xlabel="rodada de tentativas de troca", ylabel="índice do kT (1 = menor)")
    ax_E = Axis(center[2, 1]; title="Energia potencial por réplica (média móvel)", xlabel="passo (< 0 = equilibração)", ylabel="U (unidades de epsilon)")
    ax_H = Axis(center[2, 2]; title="P(U) por réplica (azul = menor kT, vermelho = maior kT)", xlabel="U (unidades de epsilon)", ylabel="P(U)")
    ax_main = Axis(center[3, 1]; aspect=DataAspect(), title="Réplica de referência (menor kT)", xticklabelsvisible=false, yticklabelsvisible=false)
    g_O = center[3, 2] = GridLayout()
    ax_O = Axis(g_O[1, 1]; aspect=DataAspect(), title="Ocupação", xlabel="índice do kT (1 = menor)", ylabel="walker")
    rowsize!(center, 1, Relative(0.27)); rowsize!(center, 2, Relative(0.27))

    # gráficos ligados a Observables: mudou o valor (obs[] = ...) -> redesenha
    pts_main = Observable(Point2f[]); col_main = Observable(RGBf(0, 0, 0))
    scatter!(ax_main, pts_main; color=col_main, markersize=9, strokewidth=1.0, strokecolor=:black)

    epts = [Observable(Point2f[]) for _ in 1:MAXREP]      # energia por estado (kT)
    wpts = [Observable(Point2f[]) for _ in 1:MAXREP]      # random walk por walker
    hpts = [Observable(Point2f[]) for _ in 1:MAXREP]      # histograma por estado (kT)
    gpts = [Observable(Point2f[]) for _ in 1:MAXREP]      # g(r) por estado (kT), no mesmo painel do P(U)
    kpts = [Observable(Point2f[]) for _ in 1:MAXREP]      # histograma do kT instantâneo por estado
    kteo = [Observable(Point2f[]) for _ in 1:MAXREP]      # distribuição canônica esperada (tracejada)
    tpts = [Observable(Point2f[]) for _ in 1:MAXREP]      # kT instantâneo ao longo do tempo, por estado
    talvo = [Observable(Point2f[]) for _ in 1:MAXREP]     # kT alvo de cada estado (tracejado)
    tcols = [Observable(tempcolor(k, MAXREP)) for k in 1:MAXREP]   # cor do estado k
    lw_e = [Observable(k == 1 ? 2.8 : 1.2) for k in 1:MAXREP]
    lw_t = [Observable(1.0) for _ in 1:MAXREP]          # série do kT: linhas finas, ela é ruidosa
    line_cols = [Observable(RGBf(0, 0, 0)) for _ in 1:MAXREP]
    lw_lines = [Observable(1.5) for _ in 1:MAXREP]
    alpha_W = Observable(0.85)                 # opacidade das linhas do random walk
    for k in MAXREP:-1:1
        lines!(ax_E, epts[k]; color=tcols[k], linewidth=lw_e[k])   # menor kT desenhado por último (fica por cima)
        lines!(ax_W, wpts[k]; color=line_cols[k], linewidth=lw_lines[k], alpha=alpha_W)
        lines!(ax_H, hpts[k]; color=tcols[k], linewidth=lw_e[k])
        lines!(ax_H, gpts[k]; color=tcols[k], linewidth=lw_e[k])
        lines!(ax_H, kpts[k]; color=tcols[k], linewidth=lw_e[k])
        lines!(ax_H, kteo[k]; color=tcols[k], linewidth=1.2, linestyle=:dash)
        lines!(ax_H, tpts[k]; color=tcols[k], linewidth=lw_t[k])
        lines!(ax_H, talvo[k]; color=tcols[k], linewidth=1.5, linestyle=:dash)
    end
    # referências da g(r): g = 1 (gás ideal) e o mínimo do LJ, r = 2^(1/6) sigma
    # desenhadas como segmentos, vazios fora do modo g(r); hlines!/vlines! com lista vazia
    # quebram o cálculo dos limites do eixo
    ref_g1 = Observable(Point2f[]); ref_rmin = Observable(Point2f[])
    lines!(ax_H, ref_g1; color=:gray50, linestyle=:dash, linewidth=1)
    lines!(ax_H, ref_rmin; color=:gray50, linestyle=:dot, linewidth=1.5)

    # rastreio da configuração de menor U: em qual kT ela está a cada tentativa de troca (pontos magenta, opcionais)
    minpts = Observable(Point2f[])
    COR_MENOR_U = RGBf(0.85, 0.0, 0.75)   # magenta: cor que não aparece na paleta dos walkers
    scatter!(ax_W, minpts; color=COR_MENOR_U, markersize=8, strokewidth=0, visible=tg_min.active)
    # com a chave ligada, as linhas dos walkers ficam claras para destacar os pontos da menor U
    on(tg_min.active) do ativo
        alpha_W[] = ativo ? 0.2 : 0.85
    end

    # mapa de ocupação: occ[estado, walker] = fração do tempo; NaN = célula não usada
    vlines!(ax_E, [0.0]; color=:gray40, linestyle=:dash, linewidth=1.5)   # fim da equilibração

    occ = Observable(fill(NaN32, MAXREP, MAXREP))
    hm = heatmap!(ax_O, occ; colormap=:Blues, colorrange=(0.0, 0.34))
    occ_pos = [Point2f(k, w) for w in 1:MAXREP for k in 1:MAXREP]
    occ_txt = Observable(fill("", MAXREP * MAXREP))
    text!(ax_O, occ_pos; text=occ_txt, align=(:center, :center), fontsize=12, color=:black)
    cb_O = Colorbar(g_O[1, 2], hm; label="fração", width=14)

    pts = [Observable(Point2f[]) for _ in 1:MAXREP]
    cols = [Observable(RGBf(0, 0, 0)) for _ in 1:MAXREP]
    axs = Axis[]

    logtext = Observable("Pronto. Configure e clique em '1. Minimizar + Equilibrar'.")
    Label(fig[2, 1:3], logtext; halign=:left, tellwidth=false)
    rowsize!(fig.layout, 2, Fixed(30))

    # ---------------- estado ----------------
    sim = Ref{Union{Nothing,Sim}}(nothing)
    running = Ref(false)        # pedido de continuar (Stop põe em false)
    busy = Ref(false)           # há uma tarefa de simulação ativa
    lastflag = fill("", MAXREP)
    xrng = MersenneTwister(1234)
    current_colors = Observable(PALETTES["Destaque"])
    Uh = [Float64[] for _ in 1:MAXREP]     # amostras de U por estado (para P(U))
    Th = [Float64[] for _ in 1:MAXREP]     # amostras do kT instantâneo (2K/nf) por estado
    Tser = [Point2f[] for _ in 1:MAXREP]   # (passo, kT instantâneo) por estado, inclusive na equilibração
    occ_count = zeros(Int, MAXREP, MAXREP)  # [estado, walker]: nº de registros
    swap_count = zeros(Int, MAXREP, MAXREP) # [walker a, walker b]: trocas aceitas entre os dois
    minU_count = zeros(Int, MAXREP, MAXREP) # [estado, walker]: nº de registros em que a menor U estava ali
    # frames do estado de menor kT durante o REMD: (passo, walker, U, x, y)
    traj = Tuple{Int,Int,Float64,Vector{Float32},Vector{Float32}}[]
    map_mode = Ref("Ocupação")
    graf_mode = Ref("P(U)")
    gr_hist = [HistGR(10.0, 0.05) for _ in 1:MAXREP]   # recriadas no reset! com rmax = min(5 sigma, L/2)
    Usm = fill(NaN, MAXREP)                 # média móvel exponencial de U por estado (só para o gráfico)
    rt_stage = zeros(Int, MAXREP)          # 0: ainda não tocou o menor kT; 1: a caminho do maior kT; 2: voltando ao menor kT
    rt_count = zeros(Int, MAXREP)          # round trips (menor kT -> maior kT -> menor kT) por walker

    mean_or_nan(v) = isempty(v) ? NaN : sum(v) / length(v)

    # faixa comum aos histogramas de todos os estados
    function faixa(amostras, R)
        lo = Inf; hi = -Inf
        for k in 1:R
            isempty(amostras[k]) && continue
            lo = min(lo, minimum(amostras[k])); hi = max(hi, maximum(amostras[k]))
        end
        return lo, hi
    end

    # painel P(U), g(r) ou kT instantâneo: só o modo escolhido tem pontos; os outros ficam vazios
    function update_hist()
        s = sim[]; s === nothing && return
        R = length(s.reps); modo = graf_mode[]
        for k in 1:MAXREP
            modo == "P(U)" || (hpts[k][] = Point2f[])
            modo == "g(r)" || (gpts[k][] = Point2f[])
            if modo != "kT instantâneo"
                kpts[k][] = Point2f[]; kteo[k][] = Point2f[]
            end
            if modo != "kT no tempo"
                tpts[k][] = Point2f[]; talvo[k][] = Point2f[]
            end
        end
        modo == "g(r)" || (ref_g1[] = Point2f[]; ref_rmin[] = Point2f[])
        if modo == "kT no tempo"
            x0 = Inf; x1 = -Inf
            for k in 1:R
                isempty(Tser[k]) && continue
                x0 = min(x0, Tser[k][1][1]); x1 = max(x1, Tser[k][end][1])
            end
            isfinite(x0) || return nothing
            for k in 1:R
                tpts[k][] = copy(Tser[k])
                talvo[k][] = [Point2f(x0, s.kTs[k]), Point2f(x1, s.kTs[k])]
            end
            nome = first(OPCOES_TERMO[findfirst(p -> p[2] == s.termostato, OPCOES_TERMO)])
            ax_H.title = "kT instantâneo ao longo do tempo (tracejado: kT alvo)\ntermostato: " * nome * " | passo < 0 = equilibração"
            ax_H.xlabel = "passo"; ax_H.ylabel = "kT instantâneo (unidades de epsilon)"
            autolimits!(ax_H)
            return nothing
        end
        if modo == "kT instantâneo"
            nf = graus_liberdade(length(s.reps[1].x), s.termostato)
            ax_H.title = @sprintf("kT instantâneo = 2K/nf (tracejado: canônico)\nmédia: %.3f no menor kT (alvo %.3f), %.3f no maior (alvo %.3f)",
                                  mean_or_nan(Th[1]), s.kTs[1], mean_or_nan(Th[R]), s.kTs[R])
            ax_H.xlabel = "kT instantâneo (unidades de epsilon)"; ax_H.ylabel = "P"
            lo, hi = faixa(Th, R)
            hi - lo < 1e-9 && return nothing
            xs = range(lo, hi; length=150)
            for k in 1:R
                kpts[k][] = hist_points(Th[k], lo, hi, 28)
                kteo[k][] = [Point2f(x, kT_inst_teorico(x, s.kTs[k], nf)) for x in xs]
            end
            autolimits!(ax_H)
            return nothing
        end
        if modo == "g(r)"
            for k in 1:MAXREP
                gpts[k][] = k <= R ? gr_pontos(gr_hist[k], length(s.reps[k].x), s.L) : Point2f[]
            end
            rmin = 2^(1 / 6) * s.sig
            rmax = length(gr_hist[1].counts) * gr_hist[1].dr
            gmax = maximum((p[2] for k in 1:R for p in gpts[k][]); init=1.0)
            ref_g1[] = [Point2f(0, 1), Point2f(rmax, 1)]
            ref_rmin[] = [Point2f(rmin, 0), Point2f(rmin, gmax)]
            ax_H.title = @sprintf("g(r) por réplica (azul = menor kT)\npontilhado: mínimo do LJ, r = %.2f; tracejado: g = 1", 2^(1 / 6) * s.sig)
            ax_H.xlabel = "r (mesma unidade de L)"; ax_H.ylabel = "g(r)"
            gr_hist[1].nframes > 0 && autolimits!(ax_H)
            return nothing
        end
        ax_H.title = "P(U) por réplica (azul = menor kT, vermelho = maior kT)"
        ax_H.xlabel = "U (unidades de epsilon)"; ax_H.ylabel = "P(U)"
        lo, hi = faixa(Uh, R)
        hi - lo < 1e-9 && return nothing
        for k in 1:R
            hpts[k][] = hist_points(Uh[k], lo, hi, 28)
        end
        autolimits!(ax_H)
        return nothing
    end

    function update_occ()
        s = sim[]; R = length(s.reps)
        m = fill(NaN32, MAXREP, MAXREP); txt = fill("", MAXREP * MAXREP)
        if map_mode[] == "Ocupação"
            # m[estado, walker] = fração do tempo do walker naquele kT (cada linha soma 1)
            for w in 1:R
                tot = max(1, sum(@view occ_count[1:R, w]))
                for k in 1:R
                    f = occ_count[k, w] / tot
                    m[k, w] = f
                    txt[(w - 1) * MAXREP + k] = @sprintf("%.2f", f)
                end
            end
            hm.colorrange = (0.0, 2.0 / R)
            ax_O.xlabel = "índice do kT (1 = menor)"; ax_O.ylabel = "walker"
            ax_O.title = @sprintf("Ocupação: fração do tempo de cada walker em cada kT (ideal = %.2f)", 1 / R)
        elseif map_mode[] == "Menor energia"
            # m[estado, walker] = fração dos registros em que a configuração de MENOR U estava naquele kT e era aquele walker
            tot = max(1, sum(@view minU_count[1:R, 1:R])); mx = 0.0
            for w in 1:R, k in 1:R
                f = minU_count[k, w] / tot
                m[k, w] = f; mx = max(mx, f)
                txt[(w - 1) * MAXREP + k] = @sprintf("%.2f", f)
            end
            pk = [sum(@view minU_count[k, 1:R]) / tot for k in 1:R]
            hm.colorrange = (0.0, max(mx, 1e-3))
            ax_O.xlabel = "índice do kT (1 = menor)"; ax_O.ylabel = "walker"
            ax_O.title = "Onde está a configuração de menor U?\nP por kT: " * join([@sprintf("%.2f", q) for q in pk], "  ")
        elseif map_mode[] == "Taxa de troca (kT)"
            # m[i, j] = P_acc = aceitas/tentadas entre os estados kT_i e kT_j (só vizinhos tentam trocar)
            for i in 1:R-1
                s.att[i] == 0 && continue
                f = s.acc[i] / s.att[i]
                m[i, i+1] = f; m[i+1, i] = f
                txt[i * MAXREP + i] = @sprintf("%.2f", f)             # célula (i, i+1)
                txt[(i - 1) * MAXREP + i + 1] = @sprintf("%.2f", f)   # célula (i+1, i)
            end
            hm.colorrange = (0.0, 1.0)
            ax_O.xlabel = "índice do kT (1 = menor)"; ax_O.ylabel = "índice do kT"
            ax_O.title = "Taxa de troca: P_acc = aceitas/tentadas por par de kT (alvo 0.2-0.3)"
        else
            # m[a, b] = fração de TODAS as trocas aceitas que foram entre os walkers a e b (simétrico)
            tot = max(1, sum(@view swap_count[1:R, 1:R]) ÷ 2)
            npairs = max(1, R * (R - 1) ÷ 2)
            for b in 1:R, a in 1:R
                a == b && continue
                f = swap_count[a, b] / tot
                m[a, b] = f
                txt[(b - 1) * MAXREP + a] = @sprintf("%.3f", f)
            end
            hm.colorrange = (0.0, 2.0 / npairs)
            ax_O.xlabel = "walker"; ax_O.ylabel = "walker"
            ax_O.title = @sprintf("Trocas entre walkers: fração do total de trocas aceitas (ideal = %.3f)", 1 / npairs)
        end
        occ[] = m; occ_txt[] = txt
        return nothing
    end

    function update_plots()
        s = sim[]; s === nothing && return
        R = length(s.reps); pal = current_colors[]
        for k in 1:R
            r = s.reps[k]
            pts[k][] = [Point2f(r.x[i], r.y[i]) for i in eachindex(r.x)]
            cols[k][] = pal[r.walker]
            pstr = (k < R && !isnan(s.last_prob[k])) ? @sprintf(" | P(->)=%.2f", s.last_prob[k]) : ""
            # segunda linha nunca vazia: um título que muda de 2 para 1 linha muda a altura do painel
            axs[k].title = @sprintf("kT=%.2f | U=%.0f%s\n%s", s.kTs[k], r.U, pstr, isempty(lastflag[k]) ? " " : lastflag[k])
        end
        kmin = argmin([r.U for r in s.reps])            # painel com a configuração de menor U: moldura magenta
        for k in 1:R
            c = k == kmin ? COR_MENOR_U : RGBf(0, 0, 0); lw = k == kmin ? 4.0 : 1.0
            ax = axs[k]
            ax.leftspinecolor = c; ax.rightspinecolor = c; ax.topspinecolor = c; ax.bottomspinecolor = c
            ax.spinewidth = lw
        end
        notify(minpts)
        pts_main[] = pts[1][]
        col_main[] = pal[s.reps[1].walker]
        ax_main.title = @sprintf("Réplica de referência (menor kT) | kT=%.2f | walker %d | U=%.0f", s.kTs[1], s.reps[1].walker, s.reps[1].U)
        ax_E.title = @sprintf("Energia potencial (média móvel) | <U>: menor kT = %.1f, maior kT = %.1f", mean_or_nan(Uh[1]), mean_or_nan(Uh[R]))
        for w in 1:MAXREP
            if w <= R
                line_cols[w][] = pal[w]
                lw_lines[w][] = (w == 1) ? 3.5 : 1.2
            end
            notify(epts[w]); notify(wpts[w])
        end
        autolimits!(ax_E); autolimits!(ax_W)
        update_occ()
        if sum(s.att) > 0
            lines_acc = [@sprintf("%.2f-%.2f: %.2f (%d/%d)", s.kTs[i], s.kTs[i+1], s.acc[i] / max(1, s.att[i]), s.acc[i], s.att[i]) for i in 1:R-1]
            stats[] = "Taxa de troca acumulada por par de kT:\n" * join(lines_acc, "\n") *
                      "\nRound trips (kT mín-máx-mín): " * string(sum(rt_count[1:R])) *
                      "\n  por walker: " * join(rt_count[1:R], " ")
        end
        return nothing
    end

    function trim!(v, n)
        length(v) > n && deleteat!(v, 1:(length(v) - n))
        return nothing
    end

    # energia (sempre) e amostras para P(U) (só durante o REMD, collect = true)
    function record_energy(collect::Bool)
        s = sim[]
        for k in eachindex(s.reps)
            U = s.reps[k].U
            Usm[k] = isnan(Usm[k]) ? U : 0.9 * Usm[k] + 0.1 * U     # janela de ~10 registros (100 passos no Run; 200 na equilibração)
            v = epts[k][]
            push!(v, Point2f(s.step, Usm[k])); trim!(v, 4000)
            push!(Tser[k], Point2f(s.step, kT_inst(s.reps[k], s.termostato))); trim!(Tser[k], 4000)
            if collect
                push!(Uh[k], s.reps[k].U); trim!(Uh[k], 20000)
                acumula_gr!(gr_hist[k], s.reps[k], s.L)     # g(r) do estado k (kT fixo), a cada 10 passos
                push!(Th[k], kT_inst(s.reps[k], s.termostato)); trim!(Th[k], 20000)
            end
        end
        if collect
            kmin = argmin([r.U for r in s.reps])
            minU_count[kmin, s.reps[kmin].walker] += 1
            r = s.reps[1]                                   # réplica no menor kT (a cada 10 passos)
            length(traj) < 50000 && push!(traj, (s.step, r.walker, r.U, Float32.(r.x), Float32.(r.y)))
        end
        return nothing
    end

    # round trip: tocou kT mín -> tocou kT máx -> voltou ao kT mín
    function record_walk()
        s = sim[]; R = length(s.reps)
        push!(minpts[], Point2f(s.nattempt, argmin([r.U for r in s.reps])))
        for k in 1:R
            w = s.reps[k].walker
            push!(wpts[w][], Point2f(s.nattempt, k))
            occ_count[k, w] += 1
            if k == 1
                rt_stage[w] == 2 && (rt_count[w] += 1)
                rt_stage[w] = 1
            elseif k == R && rt_stage[w] == 1
                rt_stage[w] = 2
            end
        end
        return nothing
    end

    function clear_series!()
        for k in 1:MAXREP
            pts[k][] = Point2f[]; epts[k][] = Point2f[]; wpts[k][] = Point2f[]; hpts[k][] = Point2f[]
            empty!(Uh[k]); empty!(Th[k]); gpts[k][] = Point2f[]; kpts[k][] = Point2f[]; kteo[k][] = Point2f[]
            empty!(Tser[k]); tpts[k][] = Point2f[]; talvo[k][] = Point2f[]
            zera!(gr_hist[k])
        end
        minpts[] = Point2f[]; empty!(traj); fill!(occ_count, 0); fill!(swap_count, 0); fill!(minU_count, 0); fill!(rt_stage, 0); fill!(rt_count, 0); fill!(lastflag, ""); fill!(Usm, NaN)
        return nothing
    end

    function reset!()
        v = [parse(Float64, strip(a)) for a in split(tbtext(tb_kts), ',') if !isempty(strip(a))]
        sort!(v); R = length(v)
        R < 2 && error("informe ao menos 2 temperaturas")
        R > MAXREP && error("máximo de $MAXREP réplicas")
        L = getf(tb_L); seed = geti(tb_seed)
        seed <= 0 && (seed = rand(1:10^6))
        nome_termo = something(menu_termo.selection[], "Langevin")
        termostato = Dict(OPCOES_TERMO)[nome_termo]
        sim[] = new_sim(v, geti(tb_n), L, getf(tb_eps), getf(tb_sig), getf(tb_dt), getf(tb_gam); seed, termostato)
        Random.seed!(xrng, seed + 99)
        rmax = min(5 * getf(tb_sig), L / 2)              # até ~5 diâmetros, sem passar de L/2
        for k in 1:MAXREP
            gr_hist[k] = HistGR(rmax, 0.05)
        end
        clear_series!()
        if length(axs) != R
            for ax in axs
                delete!(ax)
            end
            empty!(axs)
            for k in 1:R
                ax = Axis(right_grid[(k - 1) ÷ 2 + 1, (k - 1) % 2 + 1]; aspect=DataAspect(), titlesize=12)
                hidedecorations!(ax)
                scatter!(ax, pts[k]; color=cols[k], markersize=6, strokewidth=0.5, strokecolor=:black)
                push!(axs, ax)
            end
        end
        for k in 1:R
            axs[k].title = ""; axs[k].titlecolor = :black
            limits!(axs[k], 0, L, 0, L)
        end
        limits!(ax_main, 0, L, 0, L)
        limits!(ax_O, 0.5, R + 0.5, 0.5, R + 0.5)
        ax_O.xticks = 1:R; ax_O.yticks = 1:R
        ax_W.yticks = 1:R                          # índices inteiros no random walk
        hm.colorrange = (0.0, 2.0 / R)            # ideal (1/R) no meio da escala
        for k in 1:MAXREP
            tcols[k][] = tempcolor(min(k, R), R)
            lw_e[k][] = (k == 1 || k == R) ? 3.0 : 1.3      # extremos (menor e maior kT) em destaque
            lw_t[k][] = (k == 1 || k == R) ? 1.6 : 0.9
        end
        stats[] = "Taxa de troca: —"
        record_walk()
        update_plots()
        logtext[] = "Reset: $R réplicas criadas (semente $seed, termostato $nome_termo). Clique em '1. Minimizar + Equilibrar'."
        return nothing
    end

    # botão 1: minimiza e equilibra sem trocas; a equilibração fica em passos < 0
    function equilibrate()
        s = sim[]
        logtext[] = "Minimizando..."; yield()
        prepare!(s; fator_kT=getf(tb_fkT))
        neq = geti(tb_nequil)
        logtext[] = "Equilibrando cada réplica no seu kT (sem trocas)..."
        running[] = true
        for it in 1:neq
            running[] || break
            step_all!(s)
            if it % 20 == 0
                record_energy(false); update_plots()
                graf_mode[] == "kT no tempo" && update_hist()
                sleep(0.001)
            end
        end
        stopped = !running[]
        running[] = false
        for k in 1:MAXREP                       # equilibração fica em passos negativos; produção começa em 0
            epts[k][] = [Point2f(q[1] - s.step, q[2]) for q in epts[k][]]
            Tser[k] = [Point2f(q[1] - s.step, q[2]) for q in Tser[k]]
        end
        s.step = 0
        minpts[] = Point2f[]; empty!(traj); fill!(occ_count, 0); fill!(swap_count, 0); fill!(minU_count, 0); fill!(rt_stage, 0); fill!(rt_count, 0)
        for k in 1:MAXREP
            wpts[k][] = Point2f[]; zera!(gr_hist[k])
        end
        record_walk()
        update_plots()
        logtext[] = stopped ? "Equilibração interrompida (Stop). O sistema pode não estar equilibrado." :
                              "Equilibração concluída. Clique em '2. Run'."
        return nothing
    end

    # botão 2: laço principal (trocas a cada nex passos, energia a cada 10)
    function run_remd()
        s = sim[]
        s.prepared || error("minimize e equilibre primeiro (botão 1)")
        nsteps = geti(tb_nsteps); nex = max(1, geti(tb_nex)); ptime = getf(tb_pause)
        s.dt = getf(tb_dt); s.gam = getf(tb_gam)
        running[] = true
        for it in 1:nsteps
            running[] || break
            step_all!(s)
            if s.step % nex == 0
                res = attempt_exchange!(s, xrng)
                fill!(lastflag, "")
                for k in eachindex(s.reps)
                    axs[k].titlecolor = :black
                end
                swapped = false
                for (i, ok, P) in res
                    lastflag[i] = ok ? "Trocou" : "Falhou"; lastflag[i+1] = lastflag[i]
                    c = ok ? RGBf(0.1, 0.6, 0.1) : RGBf(0.8, 0.2, 0.2)
                    axs[i].titlecolor = c; axs[i+1].titlecolor = c
                    swapped |= ok
                    if ok
                        wa = s.reps[i].walker; wb = s.reps[i+1].walker
                        swap_count[wa, wb] += 1; swap_count[wb, wa] += 1
                    end
                end
                record_walk()
                update_plots(); update_hist()
                if swapped && ptime > 0
                    sleep(ptime)                 # pausa só quando houve troca aceita
                end
            end
            if s.step % 10 == 0
                record_energy(true); update_plots(); sleep(0.001)
                logtext[] = "Rodando REMD: passo $(s.step) ($it/$nsteps)"
            end
        end
        running[] = false
        update_plots(); update_hist()
        logtext[] = "Parado no passo $(s.step). '2. Run' continua daqui; Reset recomeça."
        return nothing
    end

    # Executa f em tarefa assíncrona: uma por vez, com erro mostrado no rodapé
    # @async para a janela não travar durante a simulação
    function launch(f)
        if busy[]
            logtext[] = "Já há uma simulação em andamento. Use Stop antes."
            return nothing
        end
        busy[] = true
        @async try
            f()
        catch err
            logtext[] = "ERRO: " * sprint(showerror, err)
        finally
            running[] = false; busy[] = false
        end
        return nothing
    end

    # carrega os valores de um teste do tutorial nos campos (não roda nada)
    on(menu_teste.selection) do sel
        (sel === nothing || sel == "(escolher)") && return
        v = merge(PADRAO, Dict(TESTES)[sel])
        for (k, tb) in (("nrep", tb_nrep), ("kmin", tb_kmin), ("kmax", tb_kmax), ("n", tb_n), ("L", tb_L),
                        ("dt", tb_dt), ("nsteps", tb_nsteps), ("nex", tb_nex), ("pause", tb_pause),
                        ("nequil", tb_nequil), ("gam", tb_gam), ("eps", tb_eps), ("sig", tb_sig), ("seed", tb_seed),
                        ("fatorkT", tb_fkT))
            settext!(tb, v[k])
        end
        menu_termo.i_selected[] = findfirst(==(v["termo"]), first.(OPCOES_TERMO))   # vale no botão 1 / Reset
        if haskey(v, "kts")
            settext!(tb_kts, v["kts"])
        else
            R = parse(Int, v["nrep"]); a = parse(Float64, v["kmin"]); b = parse(Float64, v["kmax"])
            settext!(tb_kts, join([@sprintf("%.3f", a * (b / a)^((i - 1) / (R - 1))) for i in 1:R], ","))
        end
        logtext[] = "Teste $sel carregado (semente 42). Clique em '1. Minimizar + Equilibrar' e depois em '2. Run'."
    end
    on(menu_cores.selection) do pal
        pal === nothing && return
        current_colors[] = PALETTES[pal]
        update_plots()
    end
    on(menu_map.selection) do sel
        sel === nothing && return
        map_mode[] = sel
        sim[] === nothing || update_occ()
    end
    on(menu_graf.selection) do sel
        sel === nothing && return
        graf_mode[] = sel
        update_hist()
    end
    on(btn_ladder.clicks) do _
        try
            R = clamp(geti(tb_nrep), 2, MAXREP); a = getf(tb_kmin); b = getf(tb_kmax)
            kTs = [a * (b / a)^((i - 1) / (R - 1)) for i in 1:R]
            settext!(tb_kts, join([@sprintf("%.3f", k) for k in kTs], ","))
            logtext[] = "Escada gerada. Clique em Reset (ou no botão 1) para aplicar."
        catch err
            logtext[] = "ERRO: " * sprint(showerror, err)
        end
    end
    on(btn_reset.clicks) do _
        if busy[]
            running[] = false
            logtext[] = "Simulação parada. Clique em Reset de novo para recriar o sistema."
        else
            try
                reset!()
            catch err
                logtext[] = "ERRO: " * sprint(showerror, err)
            end
        end
    end
    on(btn_min.clicks) do _
        launch() do
            reset!(); equilibrate()
        end
    end
    on(btn_run.clicks) do _
        launch(run_remd)
    end
    # recorta da imagem da janela a região ocupada pelos blocos 'objs' (eixo + título + rótulos) e salva
    # pasta de saída com o nome do teste escolhido no menu (ex.: "teste8_gargalo"); sem teste: "sem_teste"
    function pasta_teste()
        sel = menu_teste.selection[]
        (sel === nothing || sel == "(escolher)") && return "sem_teste"
        nome = Base.Unicode.normalize(lowercase(sel), stripmark=true)   # tira acentos
        return "teste" * strip(replace(nome, r"[^a-z0-9]+" => "_"), '_')
    end

    function crop_save(fn, img, objs; pad=10)
        H, W = size(img); sc = W / size(fig.scene)[1]          # pixels da imagem por unidade da figura
        x0 = Inf; x1 = -Inf; y0 = Inf; y1 = -Inf
        for o in objs
            bb = o.layoutobservables.computedbbox[]; pr = o.layoutobservables.protrusions[]
            x0 = min(x0, bb.origin[1] - pr.left - pad);   x1 = max(x1, bb.origin[1] + bb.widths[1] + pr.right + pad)
            y0 = min(y0, bb.origin[2] - pr.bottom - pad); y1 = max(y1, bb.origin[2] + bb.widths[2] + pr.top + pad)
        end
        # não usar 'cols' aqui: sobrescreveria a variável da interface
        icol = clamp(floor(Int, x0 * sc) + 1, 1, W):clamp(ceil(Int, x1 * sc), 1, W)
        irow = clamp(H - ceil(Int, y1 * sc) + 1, 1, H):clamp(H - floor(Int, y0 * sc), 1, H)   # linha 1 = topo
        save(fn, img[irow, icol])
        return nothing
    end

    on(btn_save.clicks) do _
        try
            s = sim[]; R = length(s.reps)
            dir = joinpath(pasta_teste(), @sprintf("remd_figs_R%d_kT%.2f-%.2f_n%d_troca%d_passo%d", R, s.kTs[1], s.kTs[end],
                           length(s.reps[1].x), geti(tb_nex), s.step))
            mkpath(dir)
            img = GLMakie.Makie.colorbuffer(fig)
            save(joinpath(dir, "janela_completa.png"), img)
            crop_save(joinpath(dir, "random_walk.png"), img, [ax_W])
            crop_save(joinpath(dir, "energia.png"), img, [ax_E])
            oldg = graf_mode[]                                  # os três gráficos do painel, seja qual for o escolhido
            for (mode, tag) in (("P(U)", "P_U"), ("g(r)", "g_r"), ("kT instantâneo", "kT_instantaneo"), ("kT no tempo", "kT_tempo"))
                graf_mode[] = mode; update_hist()
                crop_save(joinpath(dir, tag * ".png"), GLMakie.Makie.colorbuffer(fig), [ax_H])
            end
            graf_mode[] = oldg; update_hist()
            crop_save(joinpath(dir, "replica_referencia.png"), img, [ax_main])
            crop_save(joinpath(dir, "replicas.png"), img, axs)
            old = map_mode[]                                   # um PNG para cada opção do seletor "Mapa"
            for (mode, tag) in (("Ocupação", "ocupacao"), ("Menor energia", "menor_energia"),
                                ("Taxa de troca (kT)", "taxa_de_troca"), ("Trocas entre walkers", "trocas_walkers"))
                map_mode[] = mode; update_occ()
                crop_save(joinpath(dir, "mapa_" * tag * ".png"), GLMakie.Makie.colorbuffer(fig), [ax_O, cb_O])
            end
            map_mode[] = old; update_occ()
            logtext[] = "Figuras salvas em: " * joinpath(pwd(), dir)
        catch err
            logtext[] = "ERRO ao salvar: " * sprint(showerror, err)
        end
    end
    on(btn_traj.clicks) do _
        try
            isempty(traj) && error("sem frames: rode o REMD primeiro (botão 2)")
            s = sim[]; R = length(s.reps); n = length(traj[1][4])
            pasta = pasta_teste(); mkpath(pasta)
            fn = joinpath(pasta, @sprintf("remd_traj_kT%.2f_R%d_n%d_%dframes.xyz", s.kTs[1], R, n, length(traj)))
            open(fn, "w") do io
                for (st, w, U, x, y) in traj
                    println(io, n)
                    @printf(io, "passo=%d walker=%d U=%.4f kT=%.3f L=%.2f\n", st, w, U, s.kTs[1], s.L)
                    for i in 1:n
                        @printf(io, "H %.4f %.4f 0.0000\n", x[i], y[i])   # "H": o VMD nunca liga H com H
                    end
                end
            end
            # script do VMD: sem ligações automáticas, esferas de diâmetro ~sig, caixa periódica, vista ortográfica
            tcl = replace(fn, r"\.xyz$" => ".tcl")
            xyzpath = replace(abspath(fn), '\\' => '/')   # caminho absoluto com '/': o VMD no Windows abre em outra pasta
            open(tcl, "w") do io
                println(io, "mol new {$xyzpath} type xyz autobonds off waitfor all")
                println(io, "mol delrep 0 top")
                println(io, "mol representation VDW $(round(s.sig / 2 / 1.0, digits=3)) 20")   # raio = sig/2 (H: 1.0 A)
                println(io, "mol color ColorID 1")
                println(io, "mol addrep top")
                println(io, "pbc set {$(s.L) $(s.L) 10.0} -all")
                println(io, "pbc box -center origin -shiftcenter {$(s.L / 2) $(s.L / 2) 0}")
                println(io, "display projection Orthographic")
                println(io, "display resetview")
            end
            logtext[] = "Trajetória salva (.xyz + .tcl, $(length(traj)) frames) em $(abspath(pasta)). Abrir com: vmd -e " * tcl * " (Windows: Tk Console do VMD, source " * replace(abspath(tcl), '\\' => '/') * ")"
        catch err
            logtext[] = "ERRO ao salvar: " * sprint(showerror, err)
        end
    end
    on(btn_stop.clicks) do _
        running[] = false
    end

    reset!()
    return fig
end
