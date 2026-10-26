// ReplicaExchangeTutorial.jl: comportamento do site (idioma, tema, cópia, índice)
(function () {
  "use strict";
  var raiz = document.documentElement;

  function guardar(chave, valor) { try { localStorage.setItem(chave, valor); } catch (e) {} }

  // ---- idioma ----
  // Página inicial: os dois idiomas estão no mesmo HTML e o botão só troca data-lang.
  // Tutorial: cada idioma é uma página; o seletor é um link e só guarda a escolha.
  function aplicarIdioma(l) {
    raiz.dataset.lang = l;
    raiz.lang = l === "pt" ? "pt-BR" : "en";
    document.querySelectorAll("[data-idioma]").forEach(function (b) {
      b.setAttribute("aria-pressed", String(b.dataset.idioma === l));
    });
    var t = document.querySelector('meta[name="titulo-' + l + '"]');
    if (t) document.title = t.content;
  }
  document.querySelectorAll("button[data-idioma]").forEach(function (b) {
    b.addEventListener("click", function () { guardar("idioma", b.dataset.idioma); aplicarIdioma(b.dataset.idioma); });
  });
  document.querySelectorAll("a[data-idioma]").forEach(function (a) {
    a.addEventListener("click", function () { guardar("idioma", a.dataset.idioma); });
  });
  if (document.querySelector("button[data-idioma]")) aplicarIdioma(raiz.dataset.lang || "pt");

  // ---- tema claro/escuro ----
  var escuroSistema = window.matchMedia("(prefers-color-scheme: dark)");
  function temaAtual() { return raiz.dataset.theme || (escuroSistema.matches ? "dark" : "light"); }
  document.querySelectorAll(".botao-tema").forEach(function (b) {
    b.addEventListener("click", function () {
      var novo = temaAtual() === "dark" ? "light" : "dark";
      raiz.dataset.theme = novo;
      guardar("tema", novo);
    });
  });

  // ---- botões de copiar ----
  function rotulos() { return raiz.dataset.lang === "en" ? ["Copy", "Copied"] : ["Copiar", "Copiado"]; }
  function botaoCopiar(texto) {
    var b = document.createElement("button");
    b.className = "copiar"; b.type = "button"; b.textContent = rotulos()[0];
    b.addEventListener("click", function () {
      var t = typeof texto === "function" ? texto() : texto;
      var feito = function () {
        b.textContent = rotulos()[1]; b.classList.add("ok");
        setTimeout(function () { b.textContent = rotulos()[0]; b.classList.remove("ok"); }, 1500);
      };
      if (navigator.clipboard) navigator.clipboard.writeText(t).then(feito, function () {});
    });
    return b;
  }
  // tira o prompt do REPL (pkg> / julia>) do texto copiado
  function semPrompt(s) { return s.replace(/^(pkg|julia)> /gm, "").replace(/\s+$/, ""); }
  document.querySelectorAll(".comando").forEach(function (c) {
    var code = c.querySelector("code");
    c.appendChild(botaoCopiar(function () { return semPrompt(code.textContent); }));
  });
  document.querySelectorAll(".passo pre, .artigo pre").forEach(function (pre) {
    var code = pre.querySelector("code");
    if (!code) return;
    pre.appendChild(botaoCopiar(function () { return semPrompt(code.textContent); }));
  });
  document.querySelectorAll("[data-idioma]").forEach(function (b) {
    b.addEventListener("click", function () {
      document.querySelectorAll(".copiar").forEach(function (c) { c.textContent = rotulos()[0]; });
    });
  });

  // ---- índice lateral: destaca a seção visível ----
  var links = Array.prototype.slice.call(document.querySelectorAll(".indice a[href^='#']"));
  if (links.length && "IntersectionObserver" in window) {
    var porId = {};
    links.forEach(function (a) { porId[decodeURIComponent(a.getAttribute("href").slice(1))] = a; });
    var titulos = Array.prototype.slice.call(document.querySelectorAll(".artigo h2[id], .artigo h3[id]"));
    var visiveis = new Set();
    var obs = new IntersectionObserver(function (entradas) {
      entradas.forEach(function (e) { if (e.isIntersecting) visiveis.add(e.target); else visiveis.delete(e.target); });
      var atual = titulos.filter(function (h) { return visiveis.has(h); })[0];
      if (!atual) return;
      links.forEach(function (a) { a.classList.remove("ativo"); });
      var a = porId[atual.id];
      if (!a) return;
      a.classList.add("ativo");
      // mantém o item ativo visível dentro do índice lateral, sem rolar a página
      var caixa = document.querySelector(".indice");
      if (caixa && caixa.scrollHeight > caixa.clientHeight) {
        var topoA = a.getBoundingClientRect().top - caixa.getBoundingClientRect().top + caixa.scrollTop;
        if (topoA < caixa.scrollTop + 40 || topoA > caixa.scrollTop + caixa.clientHeight - 60) {
          caixa.scrollTop = topoA - caixa.clientHeight / 3;
        }
      }
    }, { rootMargin: "-80px 0px -70% 0px" });
    titulos.forEach(function (h) { obs.observe(h); });
  }
  // no celular, fechar o índice ao escolher uma seção
  document.querySelectorAll(".indice a").forEach(function (a) {
    a.addEventListener("click", function () {
      var d = a.closest("details");
      if (d && window.matchMedia("(max-width: 980px)").matches) d.open = false;
    });
  });
  var indice = document.querySelector(".indice");
  if (indice && window.matchMedia("(max-width: 980px)").matches) indice.open = false;

  // ---- voltar ao topo ----
  var topo = document.querySelector(".voltar-topo");
  if (topo) {
    window.addEventListener("scroll", function () {
      topo.classList.toggle("visivel", window.scrollY > 900);
    }, { passive: true });
  }

  // ---- animação da escada: para se o usuário pediu menos movimento ----
  if (window.matchMedia("(prefers-reduced-motion: reduce)").matches) {
    document.querySelectorAll("svg.animada").forEach(function (s) { if (s.pauseAnimations) s.pauseAnimations(); });
  }
})();
