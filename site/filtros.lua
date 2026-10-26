-- Filtros do pandoc para gerar as páginas do tutorial no site.

-- ```math ... ``` (sintaxe do GitHub) -> equação em bloco
function CodeBlock(el)
  if el.classes:includes("math") then
    return pandoc.Div({ pandoc.Para({ pandoc.Math("DisplayMath", el.text) }) }, { class = "equacao" })
  end
end

-- links entre os .md viram links entre as páginas do site
local paginas = { ["tutorial.md"] = "tutorial.html", ["tutorial_en.md"] = "tutorial-en.html" }
function Link(el)
  local alvo, ancora = el.target:match("^([^#]*)(.*)$")
  if paginas[alvo] then el.target = paginas[alvo] .. ancora end
  return el
end

-- imagens sem loading="lazy": com carregamento atrasado, os links para
-- seções (#teste-8 etc.) param no lugar errado quando as figuras de cima carregam

-- o primeiro H1 vira o título da página; o sumário escrito à mão sai,
-- porque o site gera o índice lateral
function Pandoc(doc)
  local blocos, i = doc.blocks, 1
  local saida = {}
  while i <= #blocos do
    local b = blocos[i]
    if b.t == "Header" and b.level == 1 and doc.meta.title == nil then
      doc.meta.title = pandoc.Inlines(b.content)
    elseif b.t == "Header" and b.level == 2 and
           (pandoc.utils.stringify(b) == "Sumário" or pandoc.utils.stringify(b) == "Contents") then
      -- pula o sumário e a regra horizontal que vem depois
      i = i + 1
      while i <= #blocos and blocos[i].t ~= "Header" do i = i + 1 end
      i = i - 1
    else
      table.insert(saida, b)
    end
    i = i + 1
  end
  doc.blocks = saida
  return doc
end

-- tabelas largas rolam na horizontal em telas pequenas
function Table(el)
  return pandoc.Div({ el }, { class = "tabela" })
end
