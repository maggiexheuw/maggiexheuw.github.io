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

<!-- ========================================================== -->
<!-- Greater Bay Area Summer School · Lectures 01-15            -->
<!-- PDFs are stored in the PPTX/ folder.                       -->
<!-- ========================================================== -->

<style>
  /* ---------- Page scope ---------- */
  .qm-page {
    --qm-ink: #172033;
    --qm-muted: #667085;
    --qm-line: #e7eaf0;
    --qm-surface: #ffffff;
    --qm-soft: #f7f8fb;
    --qm-accent: #4f46e5;
    --qm-accent-2: #7c3aed;
    --qm-shadow: 0 18px 45px rgba(24, 32, 51, 0.08);
    color: var(--qm-ink);
    font-size: 20px;
    line-height: 1.78;
  }

  .qm-page * {
    box-sizing: border-box;
  }

  /* ---------- Hero ---------- */
  .qm-hero {
    position: relative;
    overflow: hidden;
    margin: 1.5rem 0 2.4rem;
    padding: clamp(2rem, 5vw, 4.2rem);
    border: 1px solid rgba(79, 70, 229, 0.13);
    border-radius: 24px;
    background:
      radial-gradient(circle at 88% 14%, rgba(124, 58, 237, 0.16), transparent 30%),
      radial-gradient(circle at 10% 88%, rgba(14, 165, 233, 0.12), transparent 32%),
      linear-gradient(135deg, #fafaff 0%, #f7f8ff 48%, #fbfdff 100%);
    box-shadow: var(--qm-shadow);
  }

  .qm-hero::after {
    content: "";
    position: absolute;
    inset: 0;
    pointer-events: none;
    background-image:
      linear-gradient(rgba(79, 70, 229, 0.035) 1px, transparent 1px),
      linear-gradient(90deg, rgba(79, 70, 229, 0.035) 1px, transparent 1px);
    background-size: 28px 28px;
    -webkit-mask-image: linear-gradient(to bottom right, rgba(0,0,0,.75), transparent 72%);
    mask-image: linear-gradient(to bottom right, rgba(0,0,0,.75), transparent 72%);
  }

  .qm-eyebrow {
    position: relative;
    z-index: 1;
    display: inline-flex;
    align-items: center;
    gap: 0.5rem;
    margin-bottom: 1.1rem;
    padding: 0.45rem 0.8rem;
    border: 1px solid rgba(79, 70, 229, 0.15);
    border-radius: 999px;
    background: rgba(255, 255, 255, 0.74);
    color: #5145cd;
    font-size: 1.02rem;
    font-weight: 800;
    letter-spacing: 0.12em;
    text-transform: uppercase;
  }

  .qm-hero-title {
    position: relative;
    z-index: 1;
    max-width: 900px;
    margin: 0 0 1rem !important;
    color: #111827;
    font-size: clamp(2.45rem, 5.4vw, 4.45rem) !important;
    line-height: 1.06;
    font-weight: 850;
    letter-spacing: -0.045em;
  }

  .qm-hero-lead {
    position: relative;
    z-index: 1;
    max-width: 820px;
    margin: 0;
    color: #596275;
    font-size: clamp(1.22rem, 2.35vw, 1.52rem);
    line-height: 1.75;
  }

  .qm-meta {
    position: relative;
    z-index: 1;
    display: flex;
    flex-wrap: wrap;
    gap: 0.7rem;
    margin-top: 1.7rem;
  }

  .qm-meta span {
    padding: 0.48rem 0.78rem;
    border: 1px solid #e6e8f1;
    border-radius: 999px;
    background: rgba(255, 255, 255, 0.76);
    color: #4b5565;
    font-size: 1.12rem;
    font-weight: 650;
  }

  /* ---------- Section intro ---------- */
  .qm-section-head {
    display: flex;
    align-items: end;
    justify-content: space-between;
    gap: 1rem;
    margin: 3.5rem 0 1.5rem;
    padding-bottom: 0.85rem;
    border-bottom: 1px solid var(--qm-line);
  }

  .qm-section-head h2 {
    margin: 0 !important;
    color: #151a28;
    font-size: clamp(1.95rem, 3.4vw, 2.55rem) !important;
    line-height: 1.2;
    letter-spacing: -0.025em;
  }

  .qm-section-kicker {
    color: #8a93a4;
    font-size: 1.08rem;
    font-weight: 750;
    letter-spacing: 0.08em;
    text-transform: uppercase;
    white-space: nowrap;
  }

  /* ---------- Lecture cards ---------- */
  .qm-grid {
    display: grid;
    grid-template-columns: repeat(2, minmax(0, 1fr));
    gap: 1.5rem;
  }

  .qm-card {
    --card-accent: #4f46e5;
    position: relative;
    overflow: hidden;
    display: flex;
    flex-direction: column;
    min-height: 330px;
    padding: 2rem 2rem 1.8rem;
    border: 1px solid var(--qm-line);
    border-radius: 18px;
    background: linear-gradient(180deg, #ffffff 0%, #fdfdff 100%);
    box-shadow: 0 9px 28px rgba(17, 24, 39, 0.055);
    transition: transform 0.22s ease, box-shadow 0.22s ease, border-color 0.22s ease;
  }

  .qm-card::before {
    content: "";
    position: absolute;
    top: 0;
    left: 0;
    right: 0;
    height: 4px;
    background: var(--card-accent);
  }

  .qm-card:hover {
    transform: translateY(-4px);
    border-color: var(--card-accent);
    box-shadow: 0 18px 38px rgba(17, 24, 39, 0.095);
  }

  .qm-card-top {
    display: flex;
    align-items: center;
    justify-content: space-between;
    gap: 1rem;
    margin-bottom: 1.2rem;
  }

  .qm-number {
    display: inline-flex;
    align-items: center;
    justify-content: center;
    min-width: 3.1rem;
    height: 2.1rem;
    padding: 0 0.8rem;
    border-radius: 999px;
    background: #f5f7ff;
    color: var(--card-accent);
    font-size: 1.05rem;
    font-weight: 850;
    letter-spacing: 0.08em;
  }

  .qm-topic {
    color: #9299a8;
    font-size: 1.02rem;
    font-weight: 700;
    letter-spacing: 0.04em;
    text-transform: uppercase;
  }

  .qm-card h3 {
    margin: 0 0 0.9rem !important;
    color: #1b2232;
    font-size: 1.75rem !important;
    line-height: 1.35;
    font-weight: 800;
    letter-spacing: -0.02em;
  }

  .qm-desc {
    flex: 1;
    margin: 0 0 1.4rem;
    color: #687184;
    font-size: 1.18rem;
    line-height: 1.72;
  }

  .qm-btn {
    display: inline-flex;
    align-items: center;
    justify-content: center;
    align-self: flex-start;
    gap: 0.5rem;
    padding: 0.7rem 1.1rem;
    border: 1px solid #dfe4ee;
    border-radius: 10px;
    background: #f7f8fc;
    color: var(--card-accent) !important;
    text-decoration: none !important;
    font-size: 1.15rem;
    font-weight: 800;
    transition: background 0.18s ease, transform 0.18s ease;
  }

  .qm-btn:hover {
    background: #eef2ff;
    transform: translateY(-1px);
  }

  .qm-btn::after {
    content: "↗";
    font-size: 1rem;
  }

  /* Color accents for cards */
  .qm-card.c1  { --card-accent: #1d4ed8; }
  .qm-card.c2  { --card-accent: #2563eb; }
  .qm-card.c3  { --card-accent: #0f766e; }
  .qm-card.c4  { --card-accent: #0891b2; }
  .qm-card.c5  { --card-accent: #7c3aed; }
  .qm-card.c6  { --card-accent: #a855f7; }
  .qm-card.c7  { --card-accent: #ea580c; }
  .qm-card.c8  { --card-accent: #d97706; }
  .qm-card.c9  { --card-accent: #059669; }
  .qm-card.c10 { --card-accent: #10b981; }
  .qm-card.c11 { --card-accent: #dc2626; }
  .qm-card.c12 { --card-accent: #e11d48; }
  .qm-card.c13 { --card-accent: #ca8a04; }
  .qm-card.c14 { --card-accent: #eab308; }
  .qm-card.c15 { --card-accent: #9333ea; }

  /* ---------- Responsive ---------- */
  @media (max-width: 820px) {
    .qm-grid {
      grid-template-columns: 1fr;
    }

    .qm-card {
      min-height: auto;
    }
  }

  @media (max-width: 600px) {
    .qm-page {
      font-size: 18px;
    }

    .qm-hero {
      margin-top: 1rem;
      padding: 1.55rem;
      border-radius: 18px;
    }

    .qm-section-head {
      align-items: flex-start;
      flex-direction: column;
    }

    .qm-section-kicker {
      white-space: normal;
    }

    .qm-card {
      padding: 1.4rem;
      border-radius: 15px;
    }

    .qm-card h3 {
      font-size: 1.55rem !important;
    }
  }

  /* ---------- Dark-mode friendly ---------- */
  @media (prefers-color-scheme: dark) {
    .qm-page {
      --qm-ink: #eef2ff;
      --qm-muted: #aeb7c8;
      --qm-line: #2b3447;
      --qm-surface: #151b28;
      --qm-soft: #111827;
      --qm-shadow: none;
    }

    .qm-hero {
      border-color: #2b3550;
      background:
        radial-gradient(circle at 88% 14%, rgba(124, 58, 237, 0.20), transparent 30%),
        radial-gradient(circle at 10% 88%, rgba(14, 165, 233, 0.12), transparent 32%),
        linear-gradient(135deg, #121827 0%, #151827 48%, #111827 100%);
    }

    .qm-eyebrow,
    .qm-meta span {
      border-color: #313a50;
      background: rgba(17, 24, 39, 0.72);
    }

    .qm-hero-title,
    .qm-section-head h2,
    .qm-card h3 {
      color: #f5f7ff;
    }

    .qm-hero-lead,
    .qm-desc,
    .qm-intro,
    .qm-path span,
    .qm-footer-note {
      color: #aeb7c8;
    }

    .qm-path span,
    .qm-card {
      background: #151b28;
      box-shadow: none;
    }
  }
</style>

<div class="qm-page">

  <section class="qm-hero">
    <div class="qm-eyebrow">Greater Bay Area Summer School · Theoretical Physics</div>
    <h1 class="qm-hero-title">Greater Bay Area<br>Summer School</h1>
    <p class="qm-hero-lead">
      A curated collection of lecture notes and beamer presentations from the Greater Bay Area Summer School. 
      Delving into the Cooper problem, BCS theory, flat-band superconductivity via quantum geometry, 
      mesoscopic phenomena, and the frontiers of topological superconductivity & Majorana fermions.
    </p>
    <div class="qm-meta">
      <span>04 Topics</span>
      <span>15 Files</span>
      <span>Superconductivity & Majorana</span>
      <span>Summer School 2026</span>
    </div>
  </section>

  <!-- ========================================== -->
  <!-- TOPIC 01 -->
  <!-- ========================================== -->
  <div class="qm-section-head">
    <h2>Superconductivity: From Cooper Pairs to BCS</h2>
    <div class="qm-section-kicker">Topic 01</div>
  </div>

  <div class="qm-grid">

    <article class="qm-card c1">
      <div class="qm-card-top">
        <span class="qm-number">Note</span>
        <span class="qm-topic">Cooper Problem</span>
      </div>
      <h3>The Cooper Problem: Full Derivation</h3>
      <p class="qm-desc">
        In 1956 Leon Cooper showed that the normal Fermi sea of a metal is unstable against the formation of a bound pair of electrons, provided there is any net attractive interaction. This note provides the full derivation of the bound-state energy.
      </p>
      <a href="https://maggiexheuw.github.io/PPTX/KT-Law1.pdf" target="_blank" rel="noopener noreferrer" class="qm-btn">Open PDF</a>
    </article>

    <article class="qm-card c2">
      <div class="qm-card-top">
        <span class="qm-number">Beamer</span>
        <span class="qm-topic">Cooper Problem</span>
      </div>
      <h3>The Cooper Problem: Instability of the Fermi Sea</h3>
      <p class="qm-desc">
        Lecture slides covering the physical setup, trial wave function, Schrödinger equation, and Cooper's model interaction leading to the bound-state energy of the pair.
      </p>
      <a href="https://maggiexheuw.github.io/PPTX/KT-Law-beamer1.pdf" target="_blank" rel="noopener noreferrer" class="qm-btn">Open PDF</a>
    </article>

    <article class="qm-card c3">
      <div class="qm-card-top">
        <span class="qm-number">Note</span>
        <span class="qm-topic">BCS Theory</span>
      </div>
      <h3>From Cooper Pair to BCS Ground State</h3>
      <p class="qm-desc">
        A document developing both the single-pair instability and the full BCS variational treatment, working in the language of second quantization.
      </p>
      <a href="https://maggiexheuw.github.io/PPTX/KT-Law2.pdf" target="_blank" rel="noopener noreferrer" class="qm-btn">Open PDF</a>
    </article>

    <article class="qm-card c4">
      <div class="qm-card-top">
        <span class="qm-number">Beamer</span>
        <span class="qm-topic">BCS Theory</span>
      </div>
      <h3>From Cooper Pair to BCS Ground State</h3>
      <p class="qm-desc">
        Slides illustrating the promotion of the single-pair instability into a self-consistent many-body ground state of a coherent condensate of pairs.
      </p>
      <a href="https://maggiexheuw.github.io/PPTX/KT-Law-beamer2.pdf" target="_blank" rel="noopener noreferrer" class="qm-btn">Open PDF</a>
    </article>

  </div>

  <!-- ========================================== -->
  <!-- TOPIC 02 -->
  <!-- ========================================== -->
  <div class="qm-section-head">
    <h2>Flat-Band &amp; Quantum Geometry</h2>
    <div class="qm-section-kicker">Topic 02</div>
  </div>

  <div class="qm-grid">

    <article class="qm-card c5">
      <div class="qm-card-top">
        <span class="qm-number">Note</span>
        <span class="qm-topic">Flat-Band</span>
      </div>
      <h3>Size of a Cooper Pair in a Flat Band</h3>
      <p class="qm-desc">
        In a flat band, electrons do not move, and the familiar estimate predicts a vanishing pair size. These notes resolve the paradox by showing a new length scale encoded in the quantum metric \(g_{\mu\nu}(\mathbf{k})\).
      </p>
      <a href="https://maggiexheuw.github.io/PPTX/KT-Lwa3.pdf" target="_blank" rel="noopener noreferrer" class="qm-btn">Open PDF</a>
    </article>

    <article class="qm-card c6">
      <div class="qm-card-top">
        <span class="qm-number">Beamer</span>
        <span class="qm-topic">Flat-Band</span>
      </div>
      <h3>The Size of a Cooper Pair in a Flat Band</h3>
      <p class="qm-desc">
        Graduate lecture slides on quantum geometry and superconductivity, demonstrating how the quantum metric governs the physical size of pairs in flat bands.
      </p>
      <a href="https://maggiexheuw.github.io/PPTX/KT-Law-beamer3.pdf" target="_blank" rel="noopener noreferrer" class="qm-btn">Open PDF</a>
    </article>

  </div>

  <!-- ========================================== -->
  <!-- TOPIC 03 -->
  <!-- ========================================== -->
  <div class="qm-section-head">
    <h2>Mesoscopic Physics &amp; Localization</h2>
    <div class="qm-section-kicker">Topic 03</div>
  </div>

  <div class="qm-grid">

    <article class="qm-card c7">
      <div class="qm-card-top">
        <span class="qm-number">Note</span>
        <span class="qm-topic">Length Scales</span>
      </div>
      <h3>Fermi-Velocity-Controlled Length Scales</h3>
      <p class="qm-desc">
        Every energy scale \(E\) in a metal defines a characteristic length \(L_E \sim \hbar v_F / E\). These notes develop this master relation and tour its realizations in superconductors, Kondo screening, and localization.
      </p>
      <a href="https://maggiexheuw.github.io/PPTX/KT-Law-beamer4.pdf" target="_blank" rel="noopener noreferrer" class="qm-btn">Open PDF</a>
    </article>

    <article class="qm-card c8">
      <div class="qm-card-top">
        <span class="qm-number">Beamer</span>
        <span class="qm-topic">Length Scales</span>
      </div>
      <h3>Solid State Length Scales</h3>
      <p class="qm-desc">
        Lecture slides explaining how \(\hbar v_F \approx 0.66 \text{ eV} \cdot \text{nm}\) converts every energy scale in a metal into a length, unifying Josephson junctions and Kondo physics.
      </p>
      <!-- 修正这里的文件名 -->
      <a href="https://maggiexheuw.github.io/PPTX/L2P1SolidStateLengthScaleBeamerL2P1SolidStateLengthScaleBeamer.pdf" target="_blank" rel="noopener noreferrer" class="qm-btn">Open PDF</a>
    </article>

    <article class="qm-card c9">
      <div class="qm-card-top">
        <span class="qm-number">Note</span>
        <span class="qm-topic">Friedel Oscillations</span>
      </div>
      <h3>Friedel Oscillations in One Dimension</h3>
      <p class="qm-desc">
        A self-contained derivation of Friedel oscillations produced by a single delta-function impurity in a 1D gas of non-interacting electrons, using the T-matrix/Green's function approach.
      </p>
      <a href="https://maggiexheuw.github.io/PPTX/KT-Law-beamer5.pdf" target="_blank" rel="noopener noreferrer" class="qm-btn">Open PDF</a>
    </article>

    <article class="qm-card c10">
      <div class="qm-card-top">
        <span class="qm-number">Beamer</span>
        <span class="qm-topic">Friedel Oscillations</span>
      </div>
      <h3>Friedel Oscillations in 1D (L2P2)</h3>
      <p class="qm-desc">
        Slides on the Lehmann representation, Dyson equation, and the thermal length cutoff for Friedel oscillations in one dimension.
      </p>
      <!-- 修正这里的文件名，去掉了非法的 displayName= 尾巴 -->
      <a href="https://maggiexheuw.github.io/PPTX/L2P2FriedelOscillationsBeamer.pdf" target="_blank" rel="noopener noreferrer" class="qm-btn">Open PDF</a>
    </article>

    <article class="qm-card c11">
      <div class="qm-card-top">
        <span class="qm-number">Note</span>
        <span class="qm-topic">Anderson Localization</span>
      </div>
      <h3>Anderson Localization in One Dimension</h3>
      <p class="qm-desc">
        An introduction to Anderson localization through the exactly tractable 1D tight-binding model, deriving the transfer-matrix formulation and computing the localization length in the weak-disorder limit.
      </p>
      <a href="https://maggiexheuw.github.io/PPTX/KT-Law-beamer7.pdf" target="_blank" rel="noopener noreferrer" class="qm-btn">Open PDF</a>
    </article>

    <article class="qm-card c12">
      <div class="qm-card-top">
        <span class="qm-number">Beamer</span>
        <span class="qm-topic">Anderson Localization</span>
      </div>
      <h3>Anderson Localization in 1D (L2P3)</h3>
      <p class="qm-desc">
        Lecture slides covering transfer matrices, Lyapunov exponents, numerical QR re-orthogonalization, and the scaling theory of localization.
      </p>
      <!-- 修正这里的文件名 -->
      <a href="https://maggiexheuw.github.io/PPTX/L2P3Andersonlocalization1DBeamer.pdf" target="_blank" rel="noopener noreferrer" class="qm-btn">Open PDF</a>
    </article>

  </div>

  <!-- ========================================== -->
  <!-- TOPIC 04 -->
  <!-- ========================================== -->
  <div class="qm-section-head">
    <h2>Topological Superconductivity &amp; Majorana</h2>
    <div class="qm-section-kicker">Topic 04</div>
  </div>

  <div class="qm-grid">

    <article class="qm-card c13">
      <div class="qm-card-top">
        <span class="qm-number">Beamer</span>
        <span class="qm-topic">Topological SC</span>
      </div>
      <h3>Introduction to Topological Superconductors</h3>
      <p class="qm-desc">
        An overview of topological superconductivity, covering the 1D Kitaev model, realizing p-wave SC in condensed matter, and detection schemes like selective equal spin Andreev reflections.
      </p>
      <a href="https://maggiexheuw.github.io/PPTX/3b05a110-03e0-49e9-8dba-674de5ac3de2.pdf" target="_blank" rel="noopener noreferrer" class="qm-btn">Open PDF</a>
    </article>

    <article class="qm-card c14">
      <div class="qm-card-top">
        <span class="qm-number">Beamer</span>
        <span class="qm-topic">Majorana</span>
      </div>
      <h3>Majorana Induced Equal Spin Andreev Reflections</h3>
      <p class="qm-desc">
        A detailed presentation on the selective equal spin Andreev reflections induced by Majorana fermions in topological superconductors, exploring the underlying scattering matrices.
      </p>
      <a href="https://maggiexheuw.github.io/PPTX/b967d75c-aec2-4c62-8266-a8d4335ced0c.pdf" target="_blank" rel="noopener noreferrer" class="qm-btn">Open PDF</a>
    </article>

    <article class="qm-card c15">
      <div class="qm-card-top">
        <span class="qm-number">Beamer</span>
        <span class="qm-topic">Majorana Platforms</span>
      </div>
      <h3>Quantum Anomalous Hall Insulators</h3>
      <p class="qm-desc">
        Exploring Quantum Anomalous Hall Insulators (QAHI) as Majorana platforms. Covers quasi-1D QAHI/SC junctions, devices with cuts, and flat-band topological superconductors.
      </p>
      <a href="https://maggiexheuw.github.io/PPTX/60dfd905-fe36-4173-8bfe-c6be461c113c.pdf" target="_blank" rel="noopener noreferrer" class="qm-btn">Open PDF</a>
    </article>

  </div>

</div>