---
layout: post
title: "Advanced Theoretical Physics — Lecture Notes & Summaries"
subtitle: "Quantum Anomalies, Conformal Field Theory & Higher-Form Symmetries"
date: 2026-09-14
author: Maggie
header-img: img/EdWitten.jpg
catalog: true
---

<!-- ========================================================== -->
<!-- Advanced Theoretical Physics · Notes & Lectures            -->
<!-- PDFs are stored in the Wu-note/ folder.                     -->
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
    font-size: 18px;
    line-height: 1.75;
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
    font-size: 0.98rem;
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
    font-size: 1.06rem;
    font-weight: 650;
  }

  /* ---------- Section intro ---------- */
  .qm-section-head {
    display: flex;
    align-items: end;
    justify-content: space-between;
    gap: 1rem;
    margin: 3rem 0 1rem;
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
    font-size: 1.03rem;
    font-weight: 750;
    letter-spacing: 0.08em;
    text-transform: uppercase;
    white-space: nowrap;
  }

  .qm-intro {
    margin: 0 0 2rem;
    padding: 1rem 1.15rem;
    border-left: 3px solid #6366f1;
    border-radius: 0 12px 12px 0;
    background: #f8f9ff;
    color: #566074;
    font-size: 1.2rem;
    line-height: 1.75;
  }

  /* ---------- Topic path ---------- */
  .qm-path {
    display: flex;
    flex-wrap: wrap;
    gap: 0.65rem;
    margin: 0 0 2.2rem;
  }

  .qm-path span {
    display: inline-flex;
    align-items: center;
    gap: 0.45rem;
    padding: 0.48rem 0.78rem;
    border: 1px solid var(--qm-line);
    border-radius: 10px;
    background: var(--qm-surface);
    color: #5e6677;
    font-size: 1.03rem;
    font-weight: 650;
    box-shadow: 0 4px 12px rgba(17, 24, 39, 0.035);
  }

  .qm-path span::before {
    content: "";
    width: 7px;
    height: 7px;
    border-radius: 50%;
    background: linear-gradient(135deg, var(--qm-accent), var(--qm-accent-2));
  }

  /* ---------- Lecture cards ---------- */
  .qm-grid {
    display: grid;
    grid-template-columns: repeat(2, minmax(0, 1fr));
    gap: 1.15rem;
  }

  .qm-card {
    --card-accent: #4f46e5;
    position: relative;
    overflow: hidden;
    display: flex;
    flex-direction: column;
    min-height: 330px;
    padding: 1.85rem 1.85rem 1.7rem;
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
    border-color: color-mix(in srgb, var(--card-accent) 24%, #e7eaf0);
    box-shadow: 0 18px 38px rgba(17, 24, 39, 0.095);
  }

  .qm-card-top {
    display: flex;
    align-items: center;
    justify-content: space-between;
    gap: 1rem;
    margin-bottom: 1.1rem;
  }

  .qm-number {
    display: inline-flex;
    align-items: center;
    justify-content: center;
    min-width: 3.1rem;
    height: 2rem;
    padding: 0 0.7rem;
    border-radius: 999px;
    background: color-mix(in srgb, var(--card-accent) 10%, white);
    color: var(--card-accent);
    font-size: 0.98rem;
    font-weight: 850;
    letter-spacing: 0.08em;
  }

  .qm-topic {
    color: #9299a8;
    font-size: 0.96rem;
    font-weight: 700;
    letter-spacing: 0.04em;
    text-transform: uppercase;
  }

  .qm-card h3 {
    margin: 0 0 0.8rem !important;
    color: #1b2232;
    font-size: 1.66rem !important;
    line-height: 1.35;
    font-weight: 800;
    letter-spacing: -0.02em;
  }

  .qm-desc {
    flex: 1;
    margin: 0 0 1.35rem;
    color: #687184;
    font-size: 1.16rem;
    line-height: 1.72;
  }

  .qm-btn {
    display: inline-flex;
    align-items: center;
    justify-content: center;
    align-self: flex-start;
    gap: 0.5rem;
    padding: 0.68rem 0.95rem;
    border: 1px solid color-mix(in srgb, var(--card-accent) 18%, #e7eaf0);
    border-radius: 10px;
    background: color-mix(in srgb, var(--card-accent) 7%, white);
    color: var(--card-accent) !important;
    text-decoration: none !important;
    font-size: 1.06rem;
    font-weight: 800;
    transition: background 0.18s ease, transform 0.18s ease;
  }

  .qm-btn:hover {
    background: color-mix(in srgb, var(--card-accent) 12%, white);
    transform: translateY(-1px);
  }

  .qm-btn::after {
    content: "↗";
    font-size: 0.95rem;
  }

  .qm-card.c1 { --card-accent: #2563eb; }
  .qm-card.c2 { --card-accent: #0f766e; }
  .qm-card.c3 { --card-accent: #7c3aed; }
  .qm-card.c4 { --card-accent: #dc2626; }

  /* ---------- Compatibility fallback ---------- */
  @supports not (color: color-mix(in srgb, red 50%, white)) {
    .qm-number,
    .qm-btn {
      background: #f6f7fb;
      border-color: #e6e8ef;
    }
  }

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
      font-size: 17px;
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
      padding: 1.3rem;
      border-radius: 15px;
    }

    .qm-card h3 {
      font-size: 1.48rem !important;
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
    .qm-path span {
      color: #aeb7c8;
    }

    .qm-intro {
      background: #141b2a;
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
    <div class="qm-eyebrow">Advanced Theoretical Physics · Selected Notes</div>
    <h1 class="qm-hero-title">From Quantum Anomalies to Conformal &amp; Higher-Form Symmetries</h1>
    <p class="qm-hero-lead">
      A collection of personal study notes focusing on modern theoretical physics: advanced quantum anomalies, quantum anomalous phenomena, conformal field theories (CFT), and generalized higher-form symmetries.
    </p>
    <div class="qm-meta">
      <span>4 specialized files</span>
      <span>PDF notes</span>
      <span>Quantum field theory</span>
      <span>Updated Sep 14, 2026</span>
    </div>
  </section>

  <div class="qm-section-head">
    <h2>Lecture &amp; Study Notes</h2>
    <div class="qm-section-kicker">Core topics in modern QFT &amp; condensed matter theory</div>
  </div>

  <div class="qm-intro">
    This collection brings together key topics bridging modern quantum field theory and theoretical condensed matter: topological and quantum anomalies, conformal field theory frameworks, and the modern perspective of higher-form global symmetries.
  </div>

  <div class="qm-path">
    <span>Quantum Anomalies</span>
    <span>Quantum Anomalous Hall / Phases</span>
    <span>Conformal Field Theory (CFT)</span>
    <span>Higher-Form Symmetries</span>
  </div>

  <div class="qm-grid">

    <article class="qm-card c1">
      <div class="qm-card-top">
        <span class="qm-number">Note 01</span>
        <span class="qm-topic">Advanced QFT</span>
      </div>
      <h3>许霄焱 — 讲义笔记</h3>
      <p class="qm-desc">
        深入探讨现代量子场论与多体物理中的核心理论推导与专题分析，涵盖前沿物理机制与严谨的数学结构。
      </p>
      <a href="https://maggiexheuw.github.io/Wu-note/许霄焱.pdf" target="_blank" rel="noopener noreferrer" class="qm-btn">Open PDF</a>
    </article>

    <article class="qm-card c2">
      <div class="qm-card-top">
        <span class="qm-number">Note 02</span>
        <span class="qm-topic">Topological Physics</span>
      </div>
      <h3>姚元 — 量子反常专题</h3>
      <p class="qm-desc">
        聚焦量子反常（Quantum Anomalies）及相关拓扑物态的物理内涵，讨论微观对称性破坏与宏观输运现象的深刻联系。
      </p>
      <a href="https://maggiexheuw.github.io/Wu-note/姚元-量子反常.pdf" target="_blank" rel="noopener noreferrer" class="qm-btn">Open PDF</a>
    </article>

    <article class="qm-card c3">
      <div class="qm-card-top">
        <span class="qm-number">Note 03</span>
        <span class="qm-topic">CFT</span>
      </div>
      <h3>朱伟 — 共形场论 (CFT)</h3>
      <p class="qm-desc">
        系统梳理共形场论（Conformal Field Theory）的基础与高级内容，包括算符乘积展开（OPE）、共形对称性及关联函数。
      </p>
      <a href="https://maggiexheuw.github.io/Wu-note/朱伟-CFT.pdf" target="_blank" rel="noopener noreferrer" class="qm-btn">Open PDF</a>
    </article>

    <article class="qm-card c4">
      <div class="qm-card-top">
        <span class="qm-number">Note 04</span>
        <span class="qm-topic">Symmetry</span>
      </div>
      <h3>Higher-Form Symmetries</h3>
      <p class="qm-desc">
        现代量子场论中高阶对称性（Higher-Form Symmetries）的理论框架，探索其在规范理论、禁闭及拓扑相分类中的应用。
      </p>
      <a href="https://maggiexheuw.github.io/Wu-note/Higher-form-symmetry.pdf" target="_blank" rel="noopener noreferrer" class="qm-btn">Open PDF</a>
    </article>

  </div>

</div>