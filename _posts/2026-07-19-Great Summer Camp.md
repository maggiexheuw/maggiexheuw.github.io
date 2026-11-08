---
layout: post
title: "Greater Bay Area Summer School — Lecture Series"
subtitle: "BCS, Quantum Geometry, Mesoscopic Physics & Topological Phases"
date: 2025-11-04
author: Maggie
header-img: img/EdWitten.jpg
catalog: true
published: true
permalink: /greater-bay-area-summer-school-lecture-series/
---

<style>
/* ======================================================================
   GREATER BAY AREA SUMMER SCHOOL
   Academic Archive / Editorial Edition
   All selectors are locally scoped so they do not restyle other posts.
   ====================================================================== */
.asl-page {
  --asl-paper: #fbfaf6;
  --asl-card: #fffefa;
  --asl-ink: #202a37;
  --asl-deep: #24394d;
  --asl-subtle: #596777;
  --asl-muted: #6d7885;
  --asl-rule: #dedbd2;
  --asl-gold: #a7834e;
  --asl-wine: #7e4144;
  --asl-teal: #426d69;
  --asl-soft: #f1efe8;
  --asl-shadow: 0 14px 36px rgba(28, 39, 49, .07);
  position: relative;
  left: 50%;
  transform: translateX(-50%);
  width: min(1160px, calc(100vw - 42px));
  max-width: none;
  margin: 36px 0 100px;
  padding: 0;
  background: var(--asl-paper);
  color: var(--asl-ink);
  font-family: Georgia, 'Iowan Old Style', 'Baskerville', 'Palatino Linotype', 'Noto Serif', serif;
  font-size: 20px;
  line-height: 1.82;
  -webkit-font-smoothing: antialiased;
  text-rendering: optimizeLegibility;
  isolation: isolate;
}
.asl-page, .asl-page * { box-sizing: border-box; }
.asl-page a { text-decoration: none; }
.asl-page a:focus-visible { outline: 3px solid var(--asl-wine); outline-offset: 4px; }
.asl-page .asl-sans {
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Arial, 'Noto Sans', sans-serif;
}

/* Main book-like cover. */
.asl-page .asl-cover {
  position: relative;
  overflow: hidden;
  display: grid;
  grid-template-columns: minmax(0, 1.58fr) minmax(240px, .72fr);
  gap: 24px;
  min-height: 475px;
  padding: clamp(30px, 5.4vw, 72px) clamp(28px, 5.8vw, 74px);
  border: 1px solid #dedbd2;
  border-top: 5px solid var(--asl-deep);
  background: linear-gradient(115deg, #fffefa 0%, #f9f7f1 70%, #efede5 100%);
  box-shadow: var(--asl-shadow);
}
.asl-page .asl-cover:before {
  content: '';
  position: absolute;
  inset: 14px;
  border: 1px solid rgba(167,131,78,.24);
  pointer-events: none;
}
.asl-page .asl-cover-copy { position: relative; z-index: 2; align-self: center; }
.asl-page .asl-overline,
.asl-page .asl-cover-stats,
.asl-page .asl-section-label,
.asl-page .asl-card-eyebrow,
.asl-page .asl-nav-title,
.asl-page .asl-card-action,
.asl-page .asl-endmark,
.asl-page .asl-roadmap-index {
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Arial, sans-serif;
}
.asl-page .asl-overline {
  color: var(--asl-wine);
  font-size: 12px;
  font-weight: 750;
  line-height: 1.45;
  letter-spacing: .21em;
  text-transform: uppercase;
}
.asl-page .asl-overline:before {
  content: '';
  display: inline-block;
  width: 27px;
  height: 1px;
  margin: 0 11px 3px 0;
  background: var(--asl-wine);
}
.asl-page h1.asl-cover-title {
  margin: 21px 0 16px !important;
  color: var(--asl-deep) !important;
  font-family: Georgia, 'Iowan Old Style', 'Baskerville', 'Palatino Linotype', serif !important;
  font-size: clamp(43px, 5.5vw, 72px) !important;
  font-weight: 500 !important;
  line-height: 1.08 !important;
  letter-spacing: -.05em;
}
.asl-page .asl-cover-title em { color: var(--asl-wine); font-style: italic; font-weight: 400; }
.asl-page .asl-cover-rule { width: 78px; height: 2px; margin: 24px 0 22px; background: var(--asl-gold); }
.asl-page p.asl-cover-description {
  max-width: 660px;
  margin: 0 !important;
  color: var(--asl-subtle);
  font-size: clamp(19px, 1.8vw, 23px) !important;
  line-height: 1.77;
}
.asl-page .asl-cover-stats {
  display: flex;
  flex-wrap: wrap;
  gap: 12px 24px;
  margin-top: 26px;
  color: var(--asl-deep);
  font-size: 13px;
  letter-spacing: .08em;
  text-transform: uppercase;
}
.asl-page .asl-cover-stats span + span:before {
  content: '•';
  color: var(--asl-gold);
  padding-right: 22px;
}

/* Pure CSS/typographic academic illustration: no third-party assets. */
.asl-page .asl-figure {
  position: relative;
  align-self: center;
  justify-self: center;
  width: min(100%, 310px);
  aspect-ratio: 1;
  display: grid;
  place-items: center;
  opacity: .95;
}
.asl-page .asl-ring,
.asl-page .asl-ring:before,
.asl-page .asl-ring:after {
  position: absolute;
  display: block;
  border: 1px solid rgba(167,131,78,.52);
  border-radius: 50%;
}
.asl-page .asl-ring { width: 82%; aspect-ratio: 1; }
.asl-page .asl-ring:before { content: ''; inset: 17%; border-color: rgba(126,65,68,.46); }
.asl-page .asl-ring:after { content: ''; inset: 33%; border-color: rgba(36,57,77,.49); }
.asl-page .asl-axis { position: absolute; width: 94%; height: 1px; background: rgba(36,57,77,.17); transform: rotate(-32deg); }
.asl-page .asl-axis--b { transform: rotate(57deg); }
.asl-page .asl-diagram-center {
  position: relative;
  display: grid;
  place-items: center;
  width: 72px;
  height: 72px;
  border-radius: 50%;
  background: #f5f0e6;
  color: var(--asl-wine);
  font-size: 55px;
  line-height: 1;
}
.asl-page .asl-symbol {
  position: absolute;
  z-index: 1;
  color: var(--asl-deep);
  font-size: 24px;
  font-style: italic;
}
.asl-page .asl-symbol--a { top: 7%; left: 43%; }
.asl-page .asl-symbol--b { right: 6%; top: 48%; }
.asl-page .asl-symbol--c { bottom: 8%; left: 45%; }
.asl-page .asl-symbol--d { left: 3%; top: 48%; }
.asl-page .asl-figure-label {
  position: absolute;
  bottom: -26px;
  color: var(--asl-muted);
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Arial, sans-serif;
  font-size: 10px;
  letter-spacing: .17em;
  text-align: center;
  text-transform: uppercase;
  white-space: nowrap;
}

/* Navigation and editorial introduction. */
.asl-page .asl-nav {
  display: flex;
  align-items: center;
  gap: 16px 28px;
  flex-wrap: wrap;
  margin: 29px 0 0;
  padding: 19px 25px;
  border-block: 1px solid var(--asl-rule);
}
.asl-page .asl-nav-title {
  flex: 0 0 auto;
  color: var(--asl-wine);
  font-size: 12px;
  font-weight: 800;
  letter-spacing: .17em;
  text-transform: uppercase;
}
.asl-page .asl-nav-links { display: flex; gap: 12px 24px; flex-wrap: wrap; }
.asl-page .asl-nav-links a {
  color: var(--asl-deep) !important;
  font-size: 16px;
  font-weight: 700;
  border-bottom: 1px solid transparent;
  transition: color .15s ease, border-color .15s ease;
}
.asl-page .asl-nav-links a:hover { color: var(--asl-wine) !important; border-color: var(--asl-wine); }
.asl-page .asl-introduction { max-width: 980px; padding: 50px 15px 10px; }
.asl-page h2.asl-introduction-title {
  margin: 0 0 16px !important;
  color: var(--asl-deep) !important;
  font-size: clamp(29px, 3.3vw, 41px) !important;
  font-weight: 500 !important;
  line-height: 1.24 !important;
}
.asl-page .asl-introduction p { margin: 0 0 19px; font-size: 20px; color: var(--asl-subtle); line-height: 1.86; }
.asl-page .asl-smallcaps { color: var(--asl-wine); font-variant: small-caps; letter-spacing: .045em; }

.asl-page .asl-roadmap {
  display: grid;
  grid-template-columns: repeat(5, minmax(0, 1fr));
  margin-top: 25px;
  border: 1px solid var(--asl-rule);
  background: #f7f5ef;
}
.asl-page .asl-roadmap-step { min-width: 0; padding: 25px 19px 27px; }
.asl-page .asl-roadmap-step + .asl-roadmap-step { border-left: 1px solid var(--asl-rule); }
.asl-page .asl-roadmap-index { display: block; color: var(--asl-wine); font-size: 12px; letter-spacing: .16em; font-weight: 750; }
.asl-page .asl-roadmap-step strong { display: block; margin-top: 8px; color: var(--asl-deep); font-size: 19px; line-height: 1.28; font-weight: 600; }
.asl-page .asl-roadmap-step small { display: block; margin-top: 8px; color: var(--asl-muted); font-size: 14px; line-height: 1.47; }

/* Section heading mimics a journal chapter opener. */
.asl-page .asl-section { margin-top: 83px; scroll-margin-top: 28px; }
.asl-page .asl-section-head {
  display: grid;
  grid-template-columns: 95px minmax(0, 1fr);
  gap: 18px 23px;
  align-items: start;
  padding-bottom: 26px;
  border-bottom: 1px solid var(--asl-rule);
}
.asl-page .asl-chapter-num {
  color: #a9987b;
  font-size: clamp(65px, 7.3vw, 90px);
  font-weight: 400;
  font-style: italic;
  line-height: 1;
  letter-spacing: -.075em;
}
.asl-page .asl-section-label { color: var(--asl-wine); font-size: 12px; font-weight: 800; letter-spacing: .17em; text-transform: uppercase; }
.asl-page h2.asl-section-title {
  margin: 7px 0 10px !important;
  color: var(--asl-deep) !important;
  font-size: clamp(31px, 3.3vw, 43px) !important;
  line-height: 1.19 !important;
  font-weight: 500 !important;
  letter-spacing: -.032em;
}
.asl-page .asl-section-intro { margin: 0 !important; max-width: 920px; color: var(--asl-subtle); font-size: 19px; line-height: 1.75; }
.asl-page .asl-inquiry { display: flex; gap: 12px; margin: 17px 0 24px 118px; color: var(--asl-subtle); font-size: 16px; line-height: 1.68; }
.asl-page .asl-inquiry b { color: var(--asl-wine); font-size: 12px; font-weight: 800; font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Arial, sans-serif; white-space: nowrap; text-transform: uppercase; letter-spacing: .09em; padding-top: 4px; }

/* Resource tiles: much larger text and clear note-vs-slide hierarchy. */
.asl-page .asl-cards { display: grid; grid-template-columns: repeat(2, minmax(0, 1fr)); gap: 21px; }
.asl-page .asl-card {
  --asl-topic-color: var(--asl-wine);
  position: relative;
  display: flex;
  flex-direction: column;
  min-width: 0;
  padding: 27px 29px 25px;
  border: 1px solid var(--asl-rule);
  border-top: 3px solid var(--asl-topic-color);
  background: var(--asl-card);
  box-shadow: 0 8px 21px rgba(30,40,51,.035);
  transition: transform .2s ease, box-shadow .2s ease;
}
.asl-page .asl-card:hover { transform: translateY(-3px); box-shadow: var(--asl-shadow); }
.asl-page .asl-card--geometry { --asl-topic-color: var(--asl-teal); }
.asl-page .asl-card--mesoscopic { --asl-topic-color: #7b6549; }
.asl-page .asl-card--topology { --asl-topic-color: #586582; }
.asl-page .asl-card-top { display: flex; justify-content: space-between; align-items: baseline; gap: 10px; margin-bottom: 16px; }
.asl-page .asl-card-eyebrow { color: var(--asl-topic-color); font-size: 12px; font-weight: 800; letter-spacing: .12em; text-transform: uppercase; }
.asl-page .asl-card-count { color: var(--asl-muted); font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Arial, sans-serif; font-size: 12px; letter-spacing: .12em; white-space: nowrap; }
.asl-page h3.asl-card-title {
  margin: 0 0 12px !important;
  color: var(--asl-deep) !important;
  font-size: clamp(28px, 2.5vw, 33px) !important;
  font-weight: 500 !important;
  line-height: 1.24 !important;
  letter-spacing: -.023em;
}
.asl-page .asl-card-topic { margin: 0 0 13px; color: var(--asl-muted); font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Arial, sans-serif; font-size: 13px; letter-spacing: .055em; text-transform: uppercase; }
.asl-page p.asl-card-desc { flex: 1; margin: 0 0 26px !important; color: var(--asl-subtle); font-size: 20px !important; line-height: 1.77 !important; }
.asl-page .asl-card-footer { display: flex; align-items: center; justify-content: space-between; gap: 12px; padding-top: 14px; border-top: 1px solid var(--asl-rule); }
.asl-page .asl-card-kind { color: var(--asl-muted); font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Arial, sans-serif; font-size: 12px; text-transform: uppercase; letter-spacing: .08em; }
.asl-page a.asl-card-action {
  display: inline-flex;
  align-items: center;
  gap: 8px;
  color: var(--asl-topic-color) !important;
  font-size: 15px;
  font-weight: 750;
  letter-spacing: .015em;
  white-space: nowrap;
  border-bottom: 1px solid currentColor;
  transition: gap .18s ease;
}
.asl-page a.asl-card-action:hover { gap: 12px; }
.asl-page .asl-card-action:after { content: '↗'; font-size: 19px; line-height: 1; }

/* Colophon. */
.asl-page .asl-colophon { margin-top: 88px; padding: 27px 10px 0; border-top: 2px solid var(--asl-deep); }
.asl-page .asl-colophon-head { color: var(--asl-wine); font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Arial, sans-serif; font-size: 12px; font-weight: 800; letter-spacing: .17em; text-transform: uppercase; }
.asl-page .asl-colophon p { max-width: 960px; margin: 11px 0 0 !important; color: var(--asl-subtle); font-size: 17px; line-height: 1.75; }
.asl-page .asl-endmark { margin-top: 24px; color: var(--asl-muted); font-size: 12px; letter-spacing: .15em; text-transform: uppercase; }

/* Tablet and narrow content themes. */
@media (max-width: 890px) {
  .asl-page .asl-cover { grid-template-columns: 1fr; }
  .asl-page .asl-figure { display: none; }
  .asl-page .asl-roadmap { grid-template-columns: repeat(3, minmax(0,1fr)); }
  .asl-page .asl-roadmap-step:nth-child(4) { border-left: 0; border-top: 1px solid var(--asl-rule); }
  .asl-page .asl-roadmap-step:nth-child(5) { border-top: 1px solid var(--asl-rule); }
  .asl-page .asl-cards { grid-template-columns: 1fr; gap: 16px; }
  .asl-page .asl-card { padding: 23px; }
}
@media (max-width: 670px) {
  .asl-page { width: calc(100vw - 24px); font-size: 18px; margin-top: 12px; }
  .asl-page .asl-cover { padding: 32px 23px 38px; min-height: unset; }
  .asl-page .asl-cover:before { inset: 9px; }
  .asl-page h1.asl-cover-title { font-size: clamp(40px, 10vw, 56px) !important; }
  .asl-page p.asl-cover-description { font-size: 18.5px !important; }
  .asl-page .asl-cover-stats { gap: 5px 13px; }
  .asl-page .asl-cover-stats span + span:before { padding-right: 10px; }
  .asl-page .asl-nav { align-items: start; padding: 18px 10px; }
  .asl-page .asl-nav-links { gap: 8px 18px; }
  .asl-page .asl-nav-links a { font-size: 14px; }
  .asl-page .asl-introduction { padding: 36px 8px 7px; }
  .asl-page .asl-introduction p { font-size: 18px; }
  .asl-page .asl-roadmap { grid-template-columns: 1fr; }
  .asl-page .asl-roadmap-step { padding: 15px 20px; }
  .asl-page .asl-roadmap-step + .asl-roadmap-step { border-left: 0; border-top: 1px solid var(--asl-rule); }
  .asl-page .asl-roadmap-step strong { margin-top: 2px; }
  .asl-page .asl-roadmap-step small { margin-top: 3px; }
  .asl-page .asl-section { margin-top: 64px; }
  .asl-page .asl-section-head { grid-template-columns: 59px minmax(0, 1fr); gap: 12px; }
  .asl-page .asl-chapter-num { font-size: 55px; }
  .asl-page h2.asl-section-title { font-size: clamp(28px, 7vw, 36px) !important; }
  .asl-page .asl-section-intro { font-size: 17px; }
  .asl-page .asl-inquiry { margin: 17px 2px 23px; display: block; font-size: 16px; }
  .asl-page .asl-inquiry b { display: block; padding-bottom: 4px; }
  .asl-page .asl-cards { grid-template-columns: 1fr; }
  .asl-page .asl-card { padding: 23px; }
  .asl-page h3.asl-card-title { font-size: 27px !important; }
  .asl-page p.asl-card-desc { font-size: 19px !important; }
}
@media (prefers-reduced-motion: reduce) {
  .asl-page *, .asl-page *:before, .asl-page *:after { transition: none !important; }
}
@media print {
  .asl-page { left: auto; transform: none; width: 100%; background: white; margin: 0; }
  .asl-page .asl-card { break-inside: avoid; box-shadow: none; }
  .asl-page .asl-section { break-before: auto; }
  .asl-page .asl-figure { display: none; }
}
</style>

<div class="asl-page">
  <!-- COVER / ACADEMIC ARCHIVE -->
  <header class="asl-cover">
    <div class="asl-cover-copy">
      <div class="asl-overline">Theoretical Physics · Lecture Archive</div>
      <h1 class="asl-cover-title">Greater Bay Area<br><em>Summer School</em></h1>
      <div class="asl-cover-rule" aria-hidden="true"></div>
      <p class="asl-cover-description">
        A guided archive of lecture notes and presentations, following a conceptual path
        from the Cooper instability and BCS theory to quantum geometry, mesoscopic physics,
        and Majorana modes in topological superconductors.
      </p>
      <div class="asl-cover-stats" aria-label="Collection information">
        <span>Summer School 2026</span>
        <span>4 Chapters</span>
        <span>15 PDFs</span>
      </div>
    </div>
    <div class="asl-figure" aria-hidden="true">
      <div class="asl-axis"></div>
      <div class="asl-axis asl-axis--b"></div>
      <div class="asl-ring"></div>
      <div class="asl-diagram-center">Δ</div>
      <span class="asl-symbol asl-symbol--a">ψ</span>
      <span class="asl-symbol asl-symbol--b">g</span>
      <span class="asl-symbol asl-symbol--c">ξ</span>
      <span class="asl-symbol asl-symbol--d">γ</span>
      <div class="asl-figure-label">Pairing · Geometry · Topology</div>
    </div>
  </header>

  <nav class="asl-nav" aria-label="Lecture chapters">
    <span class="asl-nav-title">Contents</span>
    <div class="asl-nav-links">
      <a href="#pairing">01 &nbsp; Pairing &amp; BCS</a>
      <a href="#geometry">02 &nbsp; Quantum Geometry</a>
      <a href="#mesoscopic">03 &nbsp; Mesoscopic Physics</a>
      <a href="#topology">04 &nbsp; Majorana &amp; Topology</a>
    </div>
  </nav>

  <section class="asl-introduction" aria-labelledby="archive-intro-title">
    <div class="asl-overline">Preface · How to read this archive</div>
    <h2 class="asl-introduction-title" id="archive-intro-title">One physical narrative, four connected chapters.</h2>
    <p>
      This collection is organized by <span class="asl-smallcaps">physical questions</span>,
      rather than by file name. We start with why electrons pair, ask how band geometry
      affects pairing, introduce real-space scales and scattering, and finish with
      topological superconductivity and its possible signatures.
    </p>
    <p>
      Where both a written note and a presentation are available, they appear together.
      Read the notes for derivations and use the slides for a faster conceptual overview.
    </p>
  </section>

  <div class="asl-roadmap" aria-label="Suggested conceptual reading path">
    <div class="asl-roadmap-step"><span class="asl-roadmap-index">I</span><strong>Cooper instability</strong><small>Why a pair forms</small></div>
    <div class="asl-roadmap-step"><span class="asl-roadmap-index">II</span><strong>BCS coherence</strong><small>From a pair to a condensate</small></div>
    <div class="asl-roadmap-step"><span class="asl-roadmap-index">III</span><strong>Quantum metric</strong><small>Geometry in flat bands</small></div>
    <div class="asl-roadmap-step"><span class="asl-roadmap-index">IV</span><strong>Interference</strong><small>Length scales and disorder</small></div>
    <div class="asl-roadmap-step"><span class="asl-roadmap-index">V</span><strong>Majorana physics</strong><small>Topology and transport</small></div>
  </div>

  <!-- ========================================================== -->
  <!-- CHAPTER 01 -->
  <!-- ========================================================== -->
  <section class="asl-section" id="pairing" aria-labelledby="pairing-title">
    <div class="asl-section-head">
      <div class="asl-chapter-num" aria-hidden="true">01</div>
      <div>
        <div class="asl-section-label">FOUNDATIONS · THE PAIRING INSTABILITY</div>
        <h2 class="asl-section-title" id="pairing-title">From the Cooper Problem to BCS Theory</h2>
        <p class="asl-section-intro">Begin with the instability of a Fermi sea under an attractive interaction, then move from a single bound pair to the coherent, self-consistent many-body BCS ground state.</p>
      </div>
    </div>
    <p class="asl-inquiry"><b>Guiding question</b> <span>How does a two-electron bound-state instability become a macroscopic superconducting condensate?</span></p>
    <div class="asl-cards">
      <article class="asl-card asl-card--pairing" aria-labelledby="resource-01">
        <div class="asl-card-top">
          <span class="asl-card-eyebrow">LECTURE NOTE</span>
          <span class="asl-card-count">01 / 15</span>
        </div>
        <h3 class="asl-card-title" id="resource-01">The Cooper Problem: Full Derivation</h3>
        <p class="asl-card-topic">Cooper Problem</p>
        <p class="asl-card-desc">In 1956 Leon Cooper showed that the normal Fermi sea of a metal is unstable against the formation of a bound pair of electrons, provided there is any net attractive interaction. This note provides the full derivation of the bound-state energy.</p>
        <div class="asl-card-footer">
          <span class="asl-card-kind">Written notes · PDF</span>
          <a class="asl-card-action" href="https://maggiexheuw.github.io/PPTX/KT-Law1.pdf" target="_blank" rel="noopener noreferrer" aria-label="Open PDF: The Cooper Problem: Full Derivation">Read PDF</a>
        </div>
      </article>
      <article class="asl-card asl-card--pairing" aria-labelledby="resource-02">
        <div class="asl-card-top">
          <span class="asl-card-eyebrow">LECTURE SLIDES</span>
          <span class="asl-card-count">02 / 15</span>
        </div>
        <h3 class="asl-card-title" id="resource-02">The Cooper Problem: Instability of the Fermi Sea</h3>
        <p class="asl-card-topic">Cooper Problem</p>
        <p class="asl-card-desc">Lecture slides covering the physical setup, trial wave function, Schrödinger equation, and Cooper's model interaction leading to the bound-state energy of the pair.</p>
        <div class="asl-card-footer">
          <span class="asl-card-kind">Beamer slides · PDF</span>
          <a class="asl-card-action" href="https://maggiexheuw.github.io/PPTX/KT-Law-beamer1.pdf" target="_blank" rel="noopener noreferrer" aria-label="Open PDF: The Cooper Problem: Instability of the Fermi Sea">Read PDF</a>
        </div>
      </article>
      <article class="asl-card asl-card--pairing" aria-labelledby="resource-03">
        <div class="asl-card-top">
          <span class="asl-card-eyebrow">LECTURE NOTE</span>
          <span class="asl-card-count">03 / 15</span>
        </div>
        <h3 class="asl-card-title" id="resource-03">From Cooper Pair to BCS Ground State</h3>
        <p class="asl-card-topic">BCS Theory</p>
        <p class="asl-card-desc">A document developing both the single-pair instability and the full BCS variational treatment, working in the language of second quantization.</p>
        <div class="asl-card-footer">
          <span class="asl-card-kind">Written notes · PDF</span>
          <a class="asl-card-action" href="https://maggiexheuw.github.io/PPTX/KT-Law2.pdf" target="_blank" rel="noopener noreferrer" aria-label="Open PDF: From Cooper Pair to BCS Ground State">Read PDF</a>
        </div>
      </article>
      <article class="asl-card asl-card--pairing" aria-labelledby="resource-04">
        <div class="asl-card-top">
          <span class="asl-card-eyebrow">LECTURE SLIDES</span>
          <span class="asl-card-count">04 / 15</span>
        </div>
        <h3 class="asl-card-title" id="resource-04">From Cooper Pair to BCS Ground State</h3>
        <p class="asl-card-topic">BCS Theory</p>
        <p class="asl-card-desc">Slides illustrating the promotion of the single-pair instability into a self-consistent many-body ground state of a coherent condensate of pairs.</p>
        <div class="asl-card-footer">
          <span class="asl-card-kind">Beamer slides · PDF</span>
          <a class="asl-card-action" href="https://maggiexheuw.github.io/PPTX/KT-Law-beamer2.pdf" target="_blank" rel="noopener noreferrer" aria-label="Open PDF: From Cooper Pair to BCS Ground State">Read PDF</a>
        </div>
      </article>
    </div>
  </section>

  <!-- ========================================================== -->
  <!-- CHAPTER 02 -->
  <!-- ========================================================== -->
  <section class="asl-section" id="geometry" aria-labelledby="geometry-title">
    <div class="asl-section-head">
      <div class="asl-chapter-num" aria-hidden="true">02</div>
      <div>
        <div class="asl-section-label">BAND GEOMETRY · BEYOND DISPERSION</div>
        <h2 class="asl-section-title" id="geometry-title">Flat Bands &amp; Quantum Geometry</h2>
        <p class="asl-section-intro">The flat-band limit raises a subtle question: when ordinary band velocity is suppressed, what fixes the spatial extent of a Cooper pair? These materials focus on the role of the quantum metric.</p>
      </div>
    </div>
    <p class="asl-inquiry"><b>Guiding question</b> <span>What geometric information survives when a band becomes nearly dispersionless?</span></p>
    <div class="asl-cards">
      <article class="asl-card asl-card--geometry" aria-labelledby="resource-05">
        <div class="asl-card-top">
          <span class="asl-card-eyebrow">LECTURE NOTE</span>
          <span class="asl-card-count">05 / 15</span>
        </div>
        <h3 class="asl-card-title" id="resource-05">Size of a Cooper Pair in a Flat Band</h3>
        <p class="asl-card-topic">Flat-Band</p>
        <p class="asl-card-desc">In a flat band, electrons do not move, and the familiar estimate predicts a vanishing pair size. These notes resolve the paradox by showing a new length scale encoded in the quantum metric \(g_{\mu\nu}(\mathbf{k})\).</p>
        <div class="asl-card-footer">
          <span class="asl-card-kind">Written notes · PDF</span>
          <a class="asl-card-action" href="https://maggiexheuw.github.io/PPTX/KT-Lwa3.pdf" target="_blank" rel="noopener noreferrer" aria-label="Open PDF: Size of a Cooper Pair in a Flat Band">Read PDF</a>
        </div>
      </article>
      <article class="asl-card asl-card--geometry" aria-labelledby="resource-06">
        <div class="asl-card-top">
          <span class="asl-card-eyebrow">LECTURE SLIDES</span>
          <span class="asl-card-count">06 / 15</span>
        </div>
        <h3 class="asl-card-title" id="resource-06">The Size of a Cooper Pair in a Flat Band</h3>
        <p class="asl-card-topic">Flat-Band</p>
        <p class="asl-card-desc">Graduate lecture slides on quantum geometry and superconductivity, demonstrating how the quantum metric governs the physical size of pairs in flat bands.</p>
        <div class="asl-card-footer">
          <span class="asl-card-kind">Beamer slides · PDF</span>
          <a class="asl-card-action" href="https://maggiexheuw.github.io/PPTX/KT-Law-beamer3.pdf" target="_blank" rel="noopener noreferrer" aria-label="Open PDF: The Size of a Cooper Pair in a Flat Band">Read PDF</a>
        </div>
      </article>
    </div>
  </section>

  <!-- ========================================================== -->
  <!-- CHAPTER 03 -->
  <!-- ========================================================== -->
  <section class="asl-section" id="mesoscopic" aria-labelledby="mesoscopic-title">
    <div class="asl-section-head">
      <div class="asl-chapter-num" aria-hidden="true">03</div>
      <div>
        <div class="asl-section-label">REAL SPACE · INTERFERENCE &amp; DISORDER</div>
        <h2 class="asl-section-title" id="mesoscopic-title">Mesoscopic Physics &amp; Localization</h2>
        <p class="asl-section-intro">Develop physical intuition for electronic length scales, then examine how scattering produces Friedel oscillations and how disorder leads to localization in one dimension.</p>
      </div>
    </div>
    <p class="asl-inquiry"><b>Guiding question</b> <span>How do energy scales, scattering, and interference reappear as measurable real-space structures?</span></p>
    <div class="asl-cards">
      <article class="asl-card asl-card--mesoscopic" aria-labelledby="resource-07">
        <div class="asl-card-top">
          <span class="asl-card-eyebrow">LECTURE NOTE</span>
          <span class="asl-card-count">07 / 15</span>
        </div>
        <h3 class="asl-card-title" id="resource-07">Fermi-Velocity-Controlled Length Scales</h3>
        <p class="asl-card-topic">Length Scales</p>
        <p class="asl-card-desc">Every energy scale \(E\) in a metal defines a characteristic length \(L_E \sim \hbar v_F / E\). These notes develop this master relation and tour its realizations in superconductors, Kondo screening, and localization.</p>
        <div class="asl-card-footer">
          <span class="asl-card-kind">Written notes · PDF</span>
          <a class="asl-card-action" href="https://maggiexheuw.github.io/PPTX/KT-Law-beamer4.pdf" target="_blank" rel="noopener noreferrer" aria-label="Open PDF: Fermi-Velocity-Controlled Length Scales">Read PDF</a>
        </div>
      </article>
      <article class="asl-card asl-card--mesoscopic" aria-labelledby="resource-08">
        <div class="asl-card-top">
          <span class="asl-card-eyebrow">LECTURE SLIDES</span>
          <span class="asl-card-count">08 / 15</span>
        </div>
        <h3 class="asl-card-title" id="resource-08">Solid State Length Scales</h3>
        <p class="asl-card-topic">Length Scales</p>
        <p class="asl-card-desc">Lecture slides explaining how \(\hbar v_F \approx 0.66 \text{ eV} \cdot \text{nm}\) converts every energy scale in a metal into a length, unifying Josephson junctions and Kondo physics.</p>
        <div class="asl-card-footer">
          <span class="asl-card-kind">Beamer slides · PDF</span>
          <a class="asl-card-action" href="https://maggiexheuw.github.io/PPTX/L2P1SolidStateLengthScaleBeamerL2P1SolidStateLengthScaleBeamer.pdf" target="_blank" rel="noopener noreferrer" aria-label="Open PDF: Solid State Length Scales">Read PDF</a>
        </div>
      </article>
      <article class="asl-card asl-card--mesoscopic" aria-labelledby="resource-09">
        <div class="asl-card-top">
          <span class="asl-card-eyebrow">LECTURE NOTE</span>
          <span class="asl-card-count">09 / 15</span>
        </div>
        <h3 class="asl-card-title" id="resource-09">Friedel Oscillations in One Dimension</h3>
        <p class="asl-card-topic">Friedel Oscillations</p>
        <p class="asl-card-desc">A self-contained derivation of Friedel oscillations produced by a single delta-function impurity in a 1D gas of non-interacting electrons, using the T-matrix/Green's function approach.</p>
        <div class="asl-card-footer">
          <span class="asl-card-kind">Written notes · PDF</span>
          <a class="asl-card-action" href="https://maggiexheuw.github.io/PPTX/KT-Law-beamer5.pdf" target="_blank" rel="noopener noreferrer" aria-label="Open PDF: Friedel Oscillations in One Dimension">Read PDF</a>
        </div>
      </article>
      <article class="asl-card asl-card--mesoscopic" aria-labelledby="resource-10">
        <div class="asl-card-top">
          <span class="asl-card-eyebrow">LECTURE SLIDES</span>
          <span class="asl-card-count">10 / 15</span>
        </div>
        <h3 class="asl-card-title" id="resource-10">Friedel Oscillations in 1D (L2P2)</h3>
        <p class="asl-card-topic">Friedel Oscillations</p>
        <p class="asl-card-desc">Slides on the Lehmann representation, Dyson equation, and the thermal length cutoff for Friedel oscillations in one dimension.</p>
        <div class="asl-card-footer">
          <span class="asl-card-kind">Beamer slides · PDF</span>
          <a class="asl-card-action" href="https://maggiexheuw.github.io/PPTX/L2P2FriedelOscillationsBeamer.pdf" target="_blank" rel="noopener noreferrer" aria-label="Open PDF: Friedel Oscillations in 1D (L2P2)">Read PDF</a>
        </div>
      </article>
      <article class="asl-card asl-card--mesoscopic" aria-labelledby="resource-11">
        <div class="asl-card-top">
          <span class="asl-card-eyebrow">LECTURE NOTE</span>
          <span class="asl-card-count">11 / 15</span>
        </div>
        <h3 class="asl-card-title" id="resource-11">Anderson Localization in One Dimension</h3>
        <p class="asl-card-topic">Anderson Localization</p>
        <p class="asl-card-desc">An introduction to Anderson localization through the exactly tractable 1D tight-binding model, deriving the transfer-matrix formulation and computing the localization length in the weak-disorder limit.</p>
        <div class="asl-card-footer">
          <span class="asl-card-kind">Written notes · PDF</span>
          <a class="asl-card-action" href="https://maggiexheuw.github.io/PPTX/KT-Law-beamer7.pdf" target="_blank" rel="noopener noreferrer" aria-label="Open PDF: Anderson Localization in One Dimension">Read PDF</a>
        </div>
      </article>
      <article class="asl-card asl-card--mesoscopic" aria-labelledby="resource-12">
        <div class="asl-card-top">
          <span class="asl-card-eyebrow">LECTURE SLIDES</span>
          <span class="asl-card-count">12 / 15</span>
        </div>
        <h3 class="asl-card-title" id="resource-12">Anderson Localization in 1D (L2P3)</h3>
        <p class="asl-card-topic">Anderson Localization</p>
        <p class="asl-card-desc">Lecture slides covering transfer matrices, Lyapunov exponents, numerical QR re-orthogonalization, and the scaling theory of localization.</p>
        <div class="asl-card-footer">
          <span class="asl-card-kind">Beamer slides · PDF</span>
          <a class="asl-card-action" href="https://maggiexheuw.github.io/PPTX/L2P3Andersonlocalization1DBeamer.pdf" target="_blank" rel="noopener noreferrer" aria-label="Open PDF: Anderson Localization in 1D (L2P3)">Read PDF</a>
        </div>
      </article>
    </div>
  </section>

  <!-- ========================================================== -->
  <!-- CHAPTER 04 -->
  <!-- ========================================================== -->
  <section class="asl-section" id="topology" aria-labelledby="topology-title">
    <div class="asl-section-head">
      <div class="asl-chapter-num" aria-hidden="true">04</div>
      <div>
        <div class="asl-section-label">TOPOLOGICAL PHASES · MAJORANA PHYSICS</div>
        <h2 class="asl-section-title" id="topology-title">Topological Superconductivity &amp; Majorana Modes</h2>
        <p class="asl-section-intro">Close the series with the Kitaev-chain viewpoint, Majorana-induced selective equal-spin Andreev reflection, and quantum-anomalous-Hall-based superconducting platforms.</p>
      </div>
    </div>
    <p class="asl-inquiry"><b>Guiding question</b> <span>Which transport signatures and device geometries can probe Majorana physics?</span></p>
    <div class="asl-cards">
      <article class="asl-card asl-card--topology" aria-labelledby="resource-13">
        <div class="asl-card-top">
          <span class="asl-card-eyebrow">LECTURE SLIDES</span>
          <span class="asl-card-count">13 / 15</span>
        </div>
        <h3 class="asl-card-title" id="resource-13">Introduction to Topological Superconductors</h3>
        <p class="asl-card-topic">Topological SC</p>
        <p class="asl-card-desc">An overview of topological superconductivity, covering the 1D Kitaev model, realizing p-wave SC in condensed matter, and detection schemes like selective equal spin Andreev reflections.</p>
        <div class="asl-card-footer">
          <span class="asl-card-kind">Beamer slides · PDF</span>
          <a class="asl-card-action" href="https://maggiexheuw.github.io/PPTX/3b05a110-03e0-49e9-8dba-674de5ac3de2.pdf" target="_blank" rel="noopener noreferrer" aria-label="Open PDF: Introduction to Topological Superconductors">Read PDF</a>
        </div>
      </article>
      <article class="asl-card asl-card--topology" aria-labelledby="resource-14">
        <div class="asl-card-top">
          <span class="asl-card-eyebrow">LECTURE SLIDES</span>
          <span class="asl-card-count">14 / 15</span>
        </div>
        <h3 class="asl-card-title" id="resource-14">Majorana Induced Equal Spin Andreev Reflections</h3>
        <p class="asl-card-topic">Majorana</p>
        <p class="asl-card-desc">A detailed presentation on the selective equal spin Andreev reflections induced by Majorana fermions in topological superconductors, exploring the underlying scattering matrices.</p>
        <div class="asl-card-footer">
          <span class="asl-card-kind">Beamer slides · PDF</span>
          <a class="asl-card-action" href="https://maggiexheuw.github.io/PPTX/b967d75c-aec2-4c62-8266-a8d4335ced0c.pdf" target="_blank" rel="noopener noreferrer" aria-label="Open PDF: Majorana Induced Equal Spin Andreev Reflections">Read PDF</a>
        </div>
      </article>
      <article class="asl-card asl-card--topology" aria-labelledby="resource-15">
        <div class="asl-card-top">
          <span class="asl-card-eyebrow">LECTURE SLIDES</span>
          <span class="asl-card-count">15 / 15</span>
        </div>
        <h3 class="asl-card-title" id="resource-15">Quantum Anomalous Hall Insulators</h3>
        <p class="asl-card-topic">Majorana Platforms</p>
        <p class="asl-card-desc">Exploring Quantum Anomalous Hall Insulators (QAHI) as Majorana platforms. Covers quasi-1D QAHI/SC junctions, devices with cuts, and flat-band topological superconductors.</p>
        <div class="asl-card-footer">
          <span class="asl-card-kind">Beamer slides · PDF</span>
          <a class="asl-card-action" href="https://maggiexheuw.github.io/PPTX/60dfd905-fe36-4173-8bfe-c6be461c113c.pdf" target="_blank" rel="noopener noreferrer" aria-label="Open PDF: Quantum Anomalous Hall Insulators">Read PDF</a>
        </div>
      </article>
    </div>
  </section>

  <footer class="asl-colophon">
    <div class="asl-colophon-head">Colophon · A note on the collection</div>
    <p>
      All fifteen PDF destinations are preserved from the original archive.
      The chapter order is a suggested learning route; individual resources can
      also be read independently. External PDFs open in a new tab.
    </p>
    <div class="asl-endmark">Greater Bay Area Summer School &nbsp; · &nbsp; Theoretical Physics</div>
  </footer>
</div>