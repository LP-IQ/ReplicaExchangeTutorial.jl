# ReplicaExchangeTutorial.jl

[![CI](https://github.com/LP-IQ/ReplicaExchangeTutorial.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/LP-IQ/ReplicaExchangeTutorial.jl/actions/workflows/CI.yml)

**Site:** https://lp-iq.github.io/ReplicaExchangeTutorial.jl/ · tutorial em [português](docs/tutorial.md) e em [English](docs/tutorial_en.md)

> **English:** an interactive Julia GUI and tutorial on temperature replica exchange molecular dynamics (T-REMD). Install with `import Pkg; Pkg.add("ReplicaExchangeTutorial")`, start Julia with `julia -t auto` and run `using ReplicaExchangeTutorial; remd_gui()`. The full tutorial in English is in [docs/tutorial_en.md](docs/tutorial_en.md) and on the [website](https://lp-iq.github.io/ReplicaExchangeTutorial.jl/tutorial-en.html).

Tutorial sobre **Replica Exchange Molecular Dynamics (T-REMD)**: quais são os parâmetros de controle do método, como cada um afeta a amostragem e o que variar quando a amostragem falha.

A ferramenta do tutorial é uma interface gráfica interativa em Julia (`remd_gui()`) que roda um REMD num sistema simples, um fluido de Lennard-Jones em 2D, e mostra o efeito de cada parâmetro nos gráficos. Inspirada na interface do [FundamentosDMC.jl](https://github.com/m3g/FundamentosDMC.jl).

## Conteúdo

| Arquivo | O que é |
| --- | --- |
| [docs/tutorial.md](docs/tutorial.md) | O tutorial completo: teoria, instalação, visita guiada à GUI, nove testes passo a passo, diagnóstico e relação com simulações moleculares |
| [docs/tutorial_en.md](docs/tutorial_en.md) | O mesmo tutorial em inglês |
| `site/` | Página do projeto (GitHub Pages), gerada a partir dos tutoriais por `site/build.sh` |
| `src/simulacao.jl` | Núcleo da simulação: Lennard-Jones 2D, Langevin e critério de troca |
| `src/interface.jl` | A GUI usada nos testes |
| `test/runtests.jl` | Testes do núcleo (forças, termostato, critério de troca) |
| `remd_gui.jl` | Abre a GUI a partir do repositório clonado |
| `run_windows.bat` | Atalho para abrir a GUI no Windows a partir do repositório |

## Como seguir o tutorial

1. Instale e abra a GUI (seções abaixo).
2. Abra [docs/tutorial.md](docs/tutorial.md) e siga as seções na ordem: teoria (1 a 4), primeira simulação (5 a 7), parâmetros (8) e testes (9).
3. Em cada teste, mude um parâmetro por vez, com a mesma semente, e use "Salvar figuras" ao final.
4. Use a tabela de diagnóstico (seção 10) para interpretar o que viu.

## O que a GUI mostra

| Painel | Pergunta que responde |
| --- | --- |
| Painéis de réplicas (direita) | Onde está cada configuração agora? A última troca foi aceita? |
| Random walk | Cada configuração percorre todos os kT? |
| Energia potencial por réplica | O sistema equilibrou? As réplicas se separam por temperatura? |
| P(U) por réplica | As distribuições de energia de kT vizinhos se sobrepõem? |
| Réplica de referência | Qual configuração está no menor kT? |
| Mapa: ocupação | Fração do tempo de cada walker em cada kT |
| Mapa: menor energia | Em que kT está a configuração de menor U |
| Mapa: taxa de troca (kT) | Trocas aceitas / tentadas por par de kT vizinhos |
| Mapa: trocas entre walkers | Há grupos de configurações isolados? |

## Instalação

Requer [Julia](https://julialang.org/downloads/) 1.10 ou mais recente e suporte a OpenGL 3.3 (exigência do GLMakie). No REPL do Julia:

```julia
julia> import Pkg; Pkg.add("ReplicaExchangeTutorial")
```

Ou, no modo de pacotes (aperte `]` no REPL): `pkg> add ReplicaExchangeTutorial`.

## Uso

```bash
julia -t auto
```

```julia
julia> using ReplicaExchangeTutorial
julia> remd_gui()
```

`-t auto` usa uma thread por núcleo; cada réplica é propagada numa thread. A primeira execução é mais lenta, porque o Julia compila o GLMakie. Figuras e trajetórias são salvas na pasta atual do Julia (`pwd()`).

Num só comando: `julia -t auto -e 'using ReplicaExchangeTutorial; wait(display(remd_gui()))'`.

### A partir do repositório

```bash
git clone https://github.com/LP-IQ/ReplicaExchangeTutorial.jl.git
cd ReplicaExchangeTutorial.jl
julia -t auto remd_gui.jl
```

Na primeira vez, o script instala as dependências no ambiente do repositório.

### Windows

1. Instale o Julia ([julialang.org/downloads](https://julialang.org/downloads/) ou `winget install --name Julia --id 9NJNWW8PVKMN -e -s msstore`).
2. No PowerShell, abra o Julia com `julia -t auto` e siga os passos de Instalação e Uso acima.

Sem digitar comandos: baixe o repositório (**Code → Download ZIP**) e dê duplo clique em `run_windows.bat`. Detalhes e problemas comuns na seção 5.5 do [tutorial](docs/tutorial.md).

### Testes

```julia
julia> import Pkg; Pkg.test("ReplicaExchangeTutorial")
```

## Na janela

Sequência básica:

1. Ajuste o número de réplicas e os limites de kT e clique em **Gerar escada geométrica** (ou edite a lista `kTs` à mão).
2. Clique em **1. Minimizar + Equilibrar**.
3. Clique em **2. Run**. Um novo clique continua de onde parou.
4. Use **Salvar figuras** para gravar um PNG de cada gráfico e **Salvar trajetória do menor kT** para gerar um `.xyz` e um script `.tcl` (abrir com `vmd -e arquivo.tcl`).

## Modelo

- Partículas de Lennard-Jones em 2D, caixa quadrada periódica, sem raio de corte (todos os pares, imagem mínima).
- Unidades reduzidas: massa 1, energia em unidades de epsilon, temperatura como kT.
- Dinâmica de Langevin (integrador BAOAB) em cada réplica.
- Tentativas de troca entre kT vizinhos, alternando entre os pares (1-2, 3-4, …) e (2-3, 4-5, …), com o critério

```math
P = \min\left\{1,\ \exp\left[\left(\frac{1}{kT_i}-\frac{1}{kT_j}\right)\left(U_i-U_j\right)\right]\right\}
```

- Após uma troca aceita, as velocidades são reescaladas por √(kT novo / kT antigo).

## Limitações

- O cálculo de forças é O(n²), em que n é o número de partículas; acima de algumas centenas de partículas a simulação fica lenta.
- É uma ferramenta didática, não um código de produção.

## Licença

MIT. Veja [LICENSE](LICENSE).

## Autores

- Lucas Avila Pinheiro
- Eduard D. S. Mourão
