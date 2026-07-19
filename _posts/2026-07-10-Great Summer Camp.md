---
layout: post
title: "Topological Moiré Bands & Fractional Chern Insulators – Lecture Notes"
subtitle: "2026 Greater Bay Area Quantum Science Summer School · Fengcheng Wu (Wuhan University)"
date: 2026-07-19
author: Maggie
header-img: img/EdWitten.jpg
catalog: true
---

<!-- ========================================================== -->
<!-- 修改说明：将相对路径 ./PPTX/ 改为绝对 URL，指向 PPTX 文件夹 -->
<!-- 使用与 CCMP 一致的格式：https://maggiexheuw.github.io/PPTX/   -->
<!-- 文件名保留空格（%20 编码），确保链接可直接点击访问          -->
<!-- ========================================================== -->

<style>
  /* Page Introduction Styling */
  .lecture-intro {
    font-size: 1.4rem;          /* 原 1.1rem → 放大 */
    line-height: 1.8;           /* 增加行高更舒适 */
    color: #2d3748;
    margin-bottom: 2.8rem;
  }
  
  /* Grid Layout for Cards */
  .lecture-list {
    display: grid;
    grid-template-columns: 1fr;
    gap: 2rem;                  /* 卡片间距增大 */
  }

  /* Modern Card Design */
  .lecture-card {
    background: #ffffff;
    border: 1px solid #e2e8f0;
    border-radius: 12px;
    padding: 2rem 2.2rem;       /* 内边距加大 */
    box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.05), 0 2px 4px -1px rgba(0, 0, 0, 0.03);
    transition: transform 0.2s ease, box-shadow 0.2s ease;
  }
  
  .lecture-card:hover {
    transform: translateY(-3px);
    box-shadow: 0 12px 20px -6px rgba(0, 0, 0, 0.10);
  }

  /* Subtle Top Border Accents */
  .card-blue   { border-top: 5px solid #3b8bba; }
  .card-orange { border-top: 5px solid #e67e22; }
  .card-purple { border-top: 5px solid #8e44ad; }
  .card-red    { border-top: 5px solid #c0392b; }

  /* Typography */
  .lecture-title {
    margin-top: 0 !important;
    margin-bottom: 0.8rem !important;
    font-size: 1.8rem !important;    /* 原 1.3rem → 放大 */
    color: #1a202c;
    font-weight: 700;
    letter-spacing: -0.01em;
  }
  
  .lecture-desc {
    color: #2d3748;
    font-size: 1.2rem;              /* 原 1rem → 放大 */
    line-height: 1.7;
    margin-bottom: 1.8rem;
  }

  /* Button Styling */
  .download-btn {
    display: inline-block;
    color: #ffffff !important;
    padding: 0.6rem 1.8rem;         /* 按钮内边距加大 */
    border-radius: 8px;
    text-decoration: none;
    font-size: 1.05rem;             /* 原 0.9rem → 放大 */
    font-weight: 600;
    transition: filter 0.2s ease;
  }
  
  .download-btn:hover {
    filter: brightness(110%);
    text-decoration: none;
  }

  .btn-blue   { background-color: #3b8bba; }
  .btn-orange { background-color: #e67e22; }
  .btn-purple { background-color: #8e44ad; }
  .btn-red    { background-color: #c0392b; }

  /* 额外：适配移动端，防止溢出 */
  @media (max-width: 600px) {
    .lecture-title {
      font-size: 1.5rem !important;
    }
    .lecture-desc {
      font-size: 1.05rem;
    }
    .lecture-card {
      padding: 1.5rem;
    }
  }
</style>

<!-- Page Introduction -->
<div class="lecture-intro">
  <p>
    This page collects the PDF lecture notes from Prof. Fengcheng Wu's series of talks at the <strong>2026 Greater Bay Area Quantum Science Summer School</strong>. Topics include topological moiré bands, integer and fractional quantum anomalous Hall effects, collective excitations, and non-Abelian fractional quantum states. All files are stored in the <code>PPTX/</code> folder. Click to download or view online.
  </p>
</div>

<!-- Lecture List (Card Style) -->
<div class="lecture-list">

  <!-- Lecture 2.2 -->
  <div class="lecture-card card-blue">
    <h3 class="lecture-title">📘 Lecture 2.2 – Topological Moiré Bands</h3>
    <p class="lecture-desc">
      Theory of topological bands in moiré superlattices, including continuum models, layer-pseudospin skyrmions, quantum geometry, and Wilson loop methods.
    </p>
    <a href="https://maggiexheuw.github.io/PPTX/Lecture%202.2.pdf" target="_blank" class="download-btn btn-blue">📄 Download PDF</a>
  </div>

  <!-- Lecture 3.1 -->
  <div class="lecture-card card-orange">
    <h3 class="lecture-title">📘 Lecture 3.1 – Integer Chern Insulators &amp; Competing States in tMoTe₂</h3>
    <p class="lecture-desc">
      Quantum anomalous Hall insulators based on the Kane-Mele-Hubbard model, Hartree-Fock approximation, and phase diagrams with competing orders in tMoTe₂.
    </p>
    <a href="https://maggiexheuw.github.io/PPTX/Lecture%203.1.pdf" target="_blank" class="download-btn btn-orange">📄 Download PDF</a>
  </div>

  <!-- Lecture 3.2 -->
  <div class="lecture-card card-purple">
    <h3 class="lecture-title">📘 Lecture 3.2 – Collective Excitations</h3>
    <p class="lecture-desc">
      Collective excitations in integer QAHIs: excitonic optical response (Bethe-Salpeter equation) and topological magnons, plus spin models and domain walls.
    </p>
    <a href="https://maggiexheuw.github.io/PPTX/Lecture%203.2.pdf" target="_blank" class="download-btn btn-purple">📄 Download PDF</a>
  </div>

  <!-- Lecture 4 -->
  <div class="lecture-card card-red">
    <h3 class="lecture-title">📘 Lecture 4 – Abelian &amp; Non-Abelian Fractional Chern Insulators in tMoTe₂</h3>
    <p class="lecture-desc">
      Lattice analogs of fractional quantum Hall effects: skyrmion picture, generalized Landau levels, ideal vs. nonideal quantum geometry, and Abelian (Jain sequences) as well as non-Abelian (ν = 5/2) fractional states.
    </p>
    <a href="https://maggiexheuw.github.io/PPTX/Lecture%204.pdf" target="_blank" class="download-btn btn-red">📄 Download PDF</a>
  </div>

</div>