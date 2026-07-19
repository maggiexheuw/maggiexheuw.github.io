---
layout: post
title: "Additional Lecture Notes – Quantum Monte Carlo, Higher-Form Symmetry, Anomalies and CFT"
subtitle: "Collected from various summer schools and workshops"
date: 2026-07-19
author: Maggie
header-img: img/EdWitten.jpg
catalog: true
---

<!-- ========================================================== -->
<!-- 说明：本页面列出四份补充讲义，PDF 文件均存放在 PPTX/ 文件夹 -->
<!-- 链接为绝对 URL，文件名已 URL 编码（空格用 %20）            -->
<!-- 按主题分为三个 Section，每 Section 配有简介                   -->
<!-- ========================================================== -->

<style>
  /* 与主页面相同的样式，保持统一 */
  .lecture-intro {
    font-size: 1.4rem;
    line-height: 1.8;
    color: #2d3748;
    margin-bottom: 2.8rem;
  }
  
  .lecture-list {
    display: grid;
    grid-template-columns: 1fr;
    gap: 2rem;
  }

  .lecture-card {
    background: #ffffff;
    border: 1px solid #e2e8f0;
    border-radius: 12px;
    padding: 2rem 2.2rem;
    box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.05), 0 2px 4px -1px rgba(0, 0, 0, 0.03);
    transition: transform 0.2s ease, box-shadow 0.2s ease;
  }
  
  .lecture-card:hover {
    transform: translateY(-3px);
    box-shadow: 0 12px 20px -6px rgba(0, 0, 0, 0.10);
  }

  .card-green  { border-top: 5px solid #2ecc71; }
  .card-teal   { border-top: 5px solid #1abc9c; }
  .card-orange { border-top: 5px solid #e67e22; }
  .card-purple { border-top: 5px solid #8e44ad; }

  .section-header {
    font-size: 1.6rem;
    font-weight: 700;
    color: #1a202c;
    margin-top: 2.5rem;
    margin-bottom: 0.8rem;
    padding-bottom: 0.4rem;
    border-bottom: 3px solid #e2e8f0;
    letter-spacing: -0.01em;
  }

  .section-intro {
    font-size: 1.15rem;
    color: #4a5568;
    line-height: 1.7;
    margin-bottom: 1.8rem;
    padding-left: 0.5rem;
    border-left: 4px solid #a0aec0;
    padding-left: 1.2rem;
  }

  .lecture-title {
    margin-top: 0 !important;
    margin-bottom: 0.8rem !important;
    font-size: 1.8rem !important;
    color: #1a202c;
    font-weight: 700;
    letter-spacing: -0.01em;
  }
  
  .lecture-desc {
    color: #2d3748;
    font-size: 1.2rem;
    line-height: 1.7;
    margin-bottom: 1.8rem;
  }

  .download-btn {
    display: inline-block;
    color: #ffffff !important;
    padding: 0.6rem 1.8rem;
    border-radius: 8px;
    text-decoration: none;
    font-size: 1.05rem;
    font-weight: 600;
    transition: filter 0.2s ease;
  }
  
  .download-btn:hover {
    filter: brightness(110%);
    text-decoration: none;
  }

  .btn-green  { background-color: #2ecc71; }
  .btn-teal   { background-color: #1abc9c; }
  .btn-orange { background-color: #e67e22; }
  .btn-purple { background-color: #8e44ad; }

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
    .section-header {
      font-size: 1.3rem;
    }
  }
</style>

<!-- ====== 页面总介绍 ====== -->
<div class="lecture-intro">
  <p>
    This page collects four supplementary lecture notes from various sources, covering topics in <strong>Quantum Monte Carlo sign problems, higher-form symmetries, quantum anomalies, and conformal field theory</strong>. All PDF files are stored in the <code>PPTX/</code> folder. Click to download or view online.
  </p>
</div>

<!-- ============================================================ -->
<!-- Section 1: Quantum Monte Carlo & Sign Problem                  -->
<!-- ============================================================ -->

<h2 class="section-header">I. Quantum Monte Carlo &amp; The Sign Problem</h2>

<div class="section-intro">
  Quantum Monte Carlo (QMC) is one of the most powerful numerical methods for studying strongly correlated quantum many-body systems. However, its application is often severely limited by the infamous <strong>fermion sign problem</strong>, which causes the computational cost to grow exponentially with system size and inverse temperature. This section provides a comprehensive review of the origin, diagnosis, and mitigation of the sign problem across different QMC formulations — from worldlines and stochastic series expansion to determinant quantum Monte Carlo. It also discusses design principles for sign-free models, complexity-theoretic limits, and recent advances in understanding when and why the sign problem can be cured.
</div>

<div class="lecture-list">

  <!-- 1. Sign Problem -->
  <div class="lecture-card card-green">
    <h3 class="lecture-title">📘 The Sign Problem in Quantum Monte Carlo</h3>
    <p class="lecture-desc">
      <strong>Author:</strong> Xiao Yan Xu
      <br><br>
      This lecture note offers a systematic and pedagogical exposition of the fermion sign problem. It begins by tracing the microscopic origin of negative and complex weights in three representative QMC formulations: stochastic series expansion (SSE), world-line methods, and auxiliary-field determinant QMC (DQMC). The concept of the <em>average sign</em> is introduced as a quantitative measure of the problem, and the exponential cost of reweighting is derived from a free-energy difference.
      <br><br>
      The note then surveys structural principles that can remove the sign problem entirely — including Marshall rotations, determinant/Kramers pairing, fermion-bag resummation, meron-cluster algorithms, and Majorana positivity — while emphasizing that <strong>the sign problem is representation-dependent</strong>, not an intrinsic property of the Hamiltonian. Finally, it discusses fundamental limits: generic efficient curing is NP-hard, finding a stoquastic basis can be hard, and certain chiral topological phases possess intrinsic obstructions to local sign-free representations. 
      <br><br>
      <em>This is an excellent entry point for anyone seeking a deep understanding of both the practical challenges and the mathematical structure underlying the sign problem.</em>
    </p>
    <a href="https://maggiexheuw.github.io/PPTX/sign_problem.pdf" target="_blank" class="download-btn btn-green">📄 Download PDF</a>
  </div>

</div>

<!-- ============================================================ -->
<!-- Section 2: Higher-Form Symmetry & Quantum Anomalies           -->
<!-- ============================================================ -->

<h2 class="section-header">II. Higher-Form Symmetry &amp; Quantum Anomalies</h2>

<div class="section-intro">
  Generalized global symmetries — including higher-form symmetries and non-invertible symmetries — have emerged as organizing principles in modern quantum field theory and condensed matter physics. This section covers two complementary perspectives. The first note provides a pedagogical introduction to the <strong>compact boson</strong> in 1+1 dimensions, where the concepts of momentum/winding symmetries, T-duality, discrete gauging, emergent dual symmetries, and non-invertible defects all appear naturally. It then generalizes the discussion to <strong>higher-form symmetries</strong> using Maxwell theory in 3+1 dimensions as the canonical example. The second note offers an elementary treatment of <strong>'t Hooft anomalies</strong> through the lens of a (0+1)-dimensional Dirac fermion, illustrating how the conflict between U(1) gauge symmetry and charge conjugation leads to a mixed anomaly, how anomaly inflow restores symmetry in the bulk, and how anomalies imply ingappability.
</div>

<div class="lecture-list">

  <!-- 2. Higher-Form Symmetry -->
  <div class="lecture-card card-teal">
    <h3 class="lecture-title">📘 Compact Boson, Duality, Gauging, and Anomalies</h3>
    <p class="lecture-desc">
      <strong>Source:</strong> Summer school handout
      <br><br>
      This note is a self-contained pedagogical tour through the <strong>compact boson</strong> in 1+1 dimensions, arguably the simplest quantum field theory that already contains the key structures of modern generalized symmetry. It begins with canonical quantization, introduces the <em>momentum</em> and <em>winding</em> U(1) symmetries, and shows how they have a mixed 't Hooft anomaly: turning on background flux for one symmetry obstructs conservation of the other.
      <br><br>
      The note then explores T-duality, discrete gauging of a ℤ₂ subgroup, and the emergence of a dual ℤ₂ symmetry. It develops the <strong>BF pairing</strong> and the four-sector formulation (charge parity vs. twist/flux) to make Kramers-Wannier duality concrete at the level of partition functions. The discussion extends to gapped phases, showing that gauging exchanges order and disorder, and how ℤ₂ × ℤ₂ SPT phases are related by Kennedy-Tasaki duality.
      <br><br>
      The final part generalizes the framework to <strong>higher-form symmetries</strong> using Maxwell theory in 3+1 dimensions, where electric and magnetic 1-form symmetries are topological, their charged objects are Wilson and 't Hooft lines, and their mixed anomaly is captured by a 4+1-dimensional inflow action. The Coulomb, Higgs, and confined phases are reinterpreted in terms of spontaneous breaking of higher-form symmetries. 
      <br><br>
      <em>A must-read for anyone interested in the modern language of generalized symmetries.</em>
    </p>
    <a href="https://maggiexheuw.github.io/PPTX/higher_form_symmetry.pdf" target="_blank" class="download-btn btn-teal">📄 Download PDF</a>
  </div>

  <!-- 3. Quantum Anomaly -->
  <div class="lecture-card card-orange">
    <h3 class="lecture-title">📘 量子反常简介（'t Hooft Anomaly）</h3>
    <p class="lecture-desc">
      <strong>Author:</strong> 姚元 (Yuan Yao)
      <br><br>
      This note provides an elementary and highly pedagogical introduction to <strong>quantum anomalies</strong>, using a simple (0+1)-dimensional system — a single Dirac fermion on a thermal circle — to illustrate the core concepts without the technical overhead of higher dimensions.
      <br><br>
      The central problem is posed as follows: after introducing a background U(1) gauge field, the system has two symmetries — local U(1) gauge invariance and global charge conjugation (C.C.). The note shows that no regularization scheme can preserve both simultaneously: if C.C. is preserved, a large U(1) gauge transformation changes the partition function by a sign; if U(1) is preserved, C.C. fails. This is a <strong>mixed 't Hooft anomaly</strong>. 
      <br><br>
      The note then demonstrates how the anomaly can be cancelled by coupling the 0+1D system to a 1+1D bulk (anomaly inflow), and proves the important consequence that <strong>anomalous symmetries imply ingappability</strong> — the system cannot have a unique, gapped ground state. The discussion is complemented by explicit Hamiltonian analysis in the 2-dimensional Hilbert space of the Dirac fermion, making the abstract field-theoretic arguments concrete and accessible.
      <br><br>
      <em>This is an excellent bridge between formal anomaly theory and physical intuition.</em>
    </p>
    <a href="https://maggiexheuw.github.io/PPTX/anomaly.pdf" target="_blank" class="download-btn btn-orange">📄 Download PDF</a>
  </div>

</div>

<!-- ============================================================ -->
<!-- Section 3: Conformal Field Theory                             -->
<!-- ============================================================ -->

<h2 class="section-header">III. Conformal Field Theory</h2>

<div class="section-intro">
  Conformal field theory (CFT) is the language of continuous phase transitions and critical phenomena. In two dimensions, conformal symmetry becomes infinite-dimensional, making many theories exactly solvable. This section provides a self-contained introduction to 2D CFT, starting from conformal transformations and the global conformal group, developing the Virasoro algebra as the central extension of the Witt algebra, and culminating in the classification of unitary minimal models via the Kac determinant. The critical Ising model serves as a concrete illustration throughout, connecting the abstract representation theory of the Virasoro algebra to physical critical exponents and fusion rules.
</div>

<div class="lecture-list">

  <!-- 4. CFT -->
  <div class="lecture-card card-purple">
    <h3 class="lecture-title">📘 Conformal Field Theory – From Basics to Minimal Models</h3>
    <p class="lecture-desc">
      <strong>Author:</strong> 朱伟 (Wei Zhu)
      <br><br>
      This lecture note offers a self-contained and pedagogical introduction to <strong>conformal field theory in two dimensions</strong>. It is structured as a gradual ascent from the geometric notion of conformal transformations to the full machinery of the Virasoro algebra and its representations.
      <br><br>
      The note begins with the definition of conformal transformations and the global conformal group in d dimensions, deriving the constraints imposed on two- and three-point correlation functions of primary fields. It then specializes to two dimensions, where conformal transformations are holomorphic (or anti-holomorphic) functions, and the Witt algebra of infinitesimal generators is introduced. The <strong>Virasoro algebra</strong> is derived as the central extension of the Witt algebra, with the central charge c as the key parameter encoding the quantum anomaly of the energy-momentum tensor.
      <br><br>
      The note then develops radial quantization, the state-operator correspondence, and the operator product expansion (OPE). The representation theory of the Virasoro algebra is explored through highest-weight states and Verma modules. The <strong>Kac determinant</strong> is computed, revealing the existence of null states and leading to the classification of unitary minimal models with \( c = 1 - 6/[m(m+1)] \) and conformal weights \( h_{r,s} = [(m+1)r - sm]^2 - 1 / [4m(m+1)] \). 
      <br><br>
      The entire framework is illustrated with the <strong>critical Ising model</strong> (m=3, c=1/2), where the identification of the spin field σ, energy density ε, and identity I is made explicit, and the fusion rules are derived. The note also shows how the conformal data determine critical exponents such as η and ν, bridging the gap between abstract CFT and measurable physics. 
      <br><br>
      <em>An ideal companion for anyone learning CFT for the first time, or seeking a concise reference.</em>
    </p>
    <a href="https://maggiexheuw.github.io/PPTX/cft.pdf" target="_blank" class="download-btn btn-purple">📄 Download PDF</a>
  </div>

</div>