---
layout: post
title: "Lecture Notes – Topological Moiré Bands, Fractional Chern Insulators, and Advanced Topics"
subtitle: "Collected from 2026 Summer School and various workshops"
date: 2026-07-19
author: Maggie
header-img: img/EdWitten.jpg
catalog: true
---

<!-- ========================================================== -->
<!-- 说明：所有 PDF 文件均存放在 PPTX/ 文件夹中，链接为绝对 URL -->
<!-- 涵盖：吴冯成课程 + 补充讲义（QMC/高形式对称性/反常/CFT）  -->
<!--       + K.T. Law 系列（Cooper/BCS/量子几何/长度尺度等）    -->
<!-- ========================================================== -->

<style>
  .lecture-intro {
    font-size: 1.6rem;            /* 增大 */
    line-height: 1.8;
    color: #2d3748;
    margin-bottom: 2rem;
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

  .section-header {
    font-size: 1.8rem;            /* 增大 */
    font-weight: 700;
    color: #1a202c;
    margin-top: 2.5rem;
    margin-bottom: 0.8rem;
    padding-bottom: 0.4rem;
    border-bottom: 3px solid #e2e8f0;
    letter-spacing: -0.01em;
  }

  .subsection-header {
    font-size: 1.6rem;            /* 增大 */
    font-weight: 600;
    color: #1a202c;
    margin-top: 2rem;
    margin-bottom: 0.5rem;
  }

  .section-intro {
    font-size: 1.3rem;            /* 增大 */
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
    font-size: 2.0rem !important;  /* 增大 */
    color: #1a202c;
    font-weight: 700;
    letter-spacing: -0.01em;
  }
  
  .lecture-desc {
    color: #2d3748;
    font-size: 1.3rem;            /* 增大 */
    line-height: 1.7;
    margin-bottom: 1.8rem;
  }

  .download-btn {
    display: inline-block;
    color: #ffffff !important;
    padding: 0.6rem 1.8rem;
    border-radius: 8px;
    text-decoration: none;
    font-size: 1.15rem;           /* 增大 */
    font-weight: 600;
    transition: filter 0.2s ease;
  }
  
  .download-btn:hover {
    filter: brightness(110%);
    text-decoration: none;
  }

  /* 吴冯成课程颜色 */
  .card-blue   { border-top: 5px solid #3b8bba; }
  .card-orange { border-top: 5px solid #e67e22; }
  .card-purple { border-top: 5px solid #8e44ad; }
  .card-red    { border-top: 5px solid #c0392b; }
  /* 补充讲义颜色 */
  .card-green  { border-top: 5px solid #2ecc71; }
  .card-teal   { border-top: 5px solid #1abc9c; }
  .card-gold   { border-top: 5px solid #f39c12; }
  .card-violet { border-top: 5px solid #9b59b6; }
  /* K.T. Law 课程颜色 */
  .card-cooper   { border-top: 5px solid #2980b9; }
  .card-geometry { border-top: 5px solid #8e44ad; }
  .card-length   { border-top: 5px solid #d35400; }
  .card-friedel  { border-top: 5px solid #27ae60; }
  .card-local    { border-top: 5px solid #c0392b; }
  .card-scale    { border-top: 5px solid #16a085; }

  .btn-blue   { background-color: #3b8bba; }
  .btn-orange { background-color: #e67e22; }
  .btn-purple { background-color: #8e44ad; }
  .btn-red    { background-color: #c0392b; }
  .btn-green  { background-color: #2ecc71; }
  .btn-teal   { background-color: #1abc9c; }
  .btn-gold   { background-color: #f39c12; }
  .btn-violet { background-color: #9b59b6; }
  .btn-cooper   { background-color: #2980b9; }
  .btn-geometry { background-color: #8e44ad; }
  .btn-length   { background-color: #d35400; }
  .btn-friedel  { background-color: #27ae60; }
  .btn-local    { background-color: #c0392b; }
  .btn-scale    { background-color: #16a085; }

  @media (max-width: 600px) {
    .lecture-title {
      font-size: 1.6rem !important;
    }
    .lecture-desc {
      font-size: 1.1rem;
    }
    .lecture-card {
      padding: 1.5rem;
    }
    .section-header {
      font-size: 1.5rem;
    }
  }
</style>

<!-- ====== 简短总介绍 ====== -->
<div class="lecture-intro">
  <p>
    This page collects lecture notes from <strong>Prof. Fengcheng Wu</strong> (Wuhan University), <strong>Prof. K.T. Law</strong>, and additional contributors on topics including topological moiré bands, fractional Chern insulators, quantum Monte Carlo sign problems, higher‑form symmetries, quantum anomalies, conformal field theory, BCS superconductivity, quantum geometry, Anderson localization, and Friedel oscillations. All PDFs are stored in the <code>PPTX/</code> folder.
  </p>
</div>

<!-- ============================================================ -->
<!-- Section 1: 吴冯成课程                                         -->
<!-- ============================================================ -->

<h2 class="section-header">I. Topological Moiré Bands &amp; Fractional Chern Insulators</h2>

<div class="section-intro">
  Lecture notes from Prof. Fengcheng Wu's series at the <strong>2026 Greater Bay Area Quantum Science Summer School</strong>. Topics include topological moiré bands, integer and fractional quantum anomalous Hall effects, collective excitations, and non‑Abelian fractional quantum states.
</div>

<div class="lecture-list">

  <div class="lecture-card card-blue">
    <h3 class="lecture-title">📘 Lecture 2.2 – Topological Moiré Bands</h3>
    <p class="lecture-desc">Theory of topological bands in moiré superlattices, including continuum models, layer‑pseudospin skyrmions, quantum geometry, and Wilson loop methods.</p>
    <a href="https://maggiexheuw.github.io/PPTX/Lecture%202.2.pdf" target="_blank" class="download-btn btn-blue">📄 Download PDF</a>
  </div>

  <div class="lecture-card card-orange">
    <h3 class="lecture-title">📘 Lecture 3.1 – Integer Chern Insulators &amp; Competing States in tMoTe₂</h3>
    <p class="lecture-desc">Quantum anomalous Hall insulators based on the Kane‑Mele‑Hubbard model, Hartree‑Fock approximation, and phase diagrams with competing orders in tMoTe₂.</p>
    <a href="https://maggiexheuw.github.io/PPTX/Lecture%203.1.pdf" target="_blank" class="download-btn btn-orange">📄 Download PDF</a>
  </div>

  <div class="lecture-card card-purple">
    <h3 class="lecture-title">📘 Lecture 3.2 – Collective Excitations</h3>
    <p class="lecture-desc">Collective excitations in integer QAHIs: excitonic optical response (Bethe‑Salpeter equation) and topological magnons, plus spin models and domain walls.</p>
    <a href="https://maggiexheuw.github.io/PPTX/Lecture%203.2.pdf" target="_blank" class="download-btn btn-purple">📄 Download PDF</a>
  </div>

  <div class="lecture-card card-red">
    <h3 class="lecture-title">📘 Lecture 4 – Abelian &amp; Non‑Abelian Fractional Chern Insulators in tMoTe₂</h3>
    <p class="lecture-desc">Lattice analogs of fractional quantum Hall effects: skyrmion picture, generalized Landau levels, ideal vs. nonideal quantum geometry, and Abelian (Jain sequences) as well as non‑Abelian (\(\nu = 5/2\)) fractional states.</p>
    <a href="https://maggiexheuw.github.io/PPTX/Lecture%204.pdf" target="_blank" class="download-btn btn-red">📄 Download PDF</a>
  </div>

</div>

<!-- ============================================================ -->
<!-- Section 2: 补充讲义                                          -->
<!-- ============================================================ -->

<h2 class="section-header">II. Supplementary Lecture Notes</h2>

<div class="section-intro">
  Additional notes from various sources, covering <strong>quantum Monte Carlo sign problems, higher‑form symmetries, quantum anomalies, and conformal field theory</strong>.
</div>

<!-- 2.1 QMC & Sign Problem -->
<h3 class="subsection-header">II-A. Quantum Monte Carlo &amp; The Sign Problem</h3>

<div class="section-intro" style="border-left-color:#2ecc71;">
  A comprehensive review of the fermion sign problem in QMC, its origin, mitigation strategies, and complexity‑theoretic limits.
</div>

<div class="lecture-list">
  <div class="lecture-card card-green">
    <h3 class="lecture-title">📘 The Sign Problem in Quantum Monte Carlo</h3>
    <p class="lecture-desc">
      <strong>Author:</strong> Xiao Yan Xu<br><br>
      Pedagogical exposition covering negative/complex weights in SSE, world‑line, and DQMC; average sign and reweighting cost; structural cures (Marshall rotations, determinant pairing, fermion bags, merons, Majorana positivity); and fundamental limits (NP‑hardness, stoquastic‑basis search, topological obstructions).
    </p>
    <a href="https://maggiexheuw.github.io/PPTX/sign_problem.pdf" target="_blank" class="download-btn btn-green">📄 Download PDF</a>
  </div>
</div>

<!-- 2.2 Higher-Form Symmetry & Anomalies -->
<h3 class="subsection-header">II-B. Higher‑Form Symmetry &amp; Quantum Anomalies</h3>

<div class="section-intro" style="border-left-color:#1abc9c;">
  Modern generalized symmetries and 't Hooft anomalies, illustrated via the compact boson, Maxwell theory, and a (0+1)‑dimensional Dirac fermion.
</div>

<div class="lecture-list">
  <div class="lecture-card card-teal">
    <h3 class="lecture-title">📘 Compact Boson, Duality, Gauging, and Anomalies</h3>
    <p class="lecture-desc">
      <strong>Source:</strong> Summer school handout<br><br>
      Self‑contained tour of 1+1D compact boson: momentum/winding \(\mathrm{U}(1)\) symmetries, mixed anomaly, T‑duality, discrete gauging, emergent dual symmetries, BF pairing, Kramers‑Wannier duality, and generalization to higher‑form symmetries in Maxwell theory (electric/magnetic 1‑form symmetries, Coulomb/Higgs/confined phases).
    </p>
    <a href="https://maggiexheuw.github.io/PPTX/higher_form_symmetry.pdf" target="_blank" class="download-btn btn-teal">📄 Download PDF</a>
  </div>

  <div class="lecture-card card-gold">
    <h3 class="lecture-title">📘 量子反常简介（'t Hooft Anomaly）</h3>
    <p class="lecture-desc">
      <strong>Author:</strong> 姚元 (Yuan Yao)<br><br>
      Elementary introduction using a (0+1)‑dimensional Dirac fermion on a thermal circle. Shows the conflict between local \(\mathrm{U}(1)\) gauge invariance and global charge conjugation, leading to a mixed anomaly; demonstrates anomaly inflow and proves that anomalous symmetries imply ingappability.
    </p>
    <a href="https://maggiexheuw.github.io/PPTX/anomaly.pdf" target="_blank" class="download-btn btn-gold">📄 Download PDF</a>
  </div>
</div>

<!-- 2.3 Conformal Field Theory -->
<h3 class="subsection-header">II-C. Conformal Field Theory</h3>

<div class="section-intro" style="border-left-color:#9b59b6;">
  A pedagogical introduction to 2D CFT, from conformal transformations and the Virasoro algebra to unitary minimal models and the critical Ising model.
</div>

<div class="lecture-list">
  <div class="lecture-card card-violet">
    <h3 class="lecture-title">📘 Conformal Field Theory – From Basics to Minimal Models</h3>
    <p class="lecture-desc">
      <strong>Author:</strong> 朱伟 (Wei Zhu)<br><br>
      Self‑contained notes covering: conformal transformations in \(d\) dimensions, primary fields and correlation functions, holomorphic nature in 2D, Witt and Virasoro algebras, central charge, radial quantization, state‑operator correspondence, OPE, highest‑weight representations, Kac determinant, unitary minimal models, and the critical Ising model as a worked example (fusion rules, critical exponents).
    </p>
    <a href="https://maggiexheuw.github.io/PPTX/cft.pdf" target="_blank" class="download-btn btn-violet">📄 Download PDF</a>
  </div>
</div>

<!-- ============================================================ -->
<!-- Section 3: K.T. Law 系列课程                                 -->
<!-- ============================================================ -->

<h2 class="section-header">III. K.T. Law Lecture Series – Superconductivity, Quantum Geometry, and Disorder</h2>

<div class="section-intro">
  Graduate lecture series by <strong>Prof. K.T. Law</strong> covering the microscopic theory of superconductivity (Cooper problem &amp; BCS), quantum geometry in flat‑band superconductivity, Fermi‑velocity‑controlled length scales, Friedel oscillations, and Anderson localization in one dimension. Includes both Beamer slides and detailed lecture notes.
</div>

<!-- 3.1 Cooper Problem & BCS -->
<h3 class="subsection-header">III-A. Cooper Problem &amp; BCS Theory of Superconductivity</h3>

<div class="section-intro" style="border-left-color:#2980b9;">
  The two pillars of microscopic superconductivity: Cooper's instability and the BCS variational ground state.
</div>

<div class="lecture-list">

  <div class="lecture-card card-cooper">
    <h3 class="lecture-title">📊 Beamer – The Cooper Problem</h3>
    <p class="lecture-desc">
      <strong>Type:</strong> Presentation slides<br><br>
      Introduction to the Cooper problem: physical setup, trial wavefunction, two‑body Schrödinger equation, model interaction, self‑consistency, and bound‑state energy \( \Delta = 2\hbar\omega_D e^{-2/N(0)V} \). Emphasizes the essential singularity and instability of the normal state.
    </p>
    <a href="https://maggiexheuw.github.io/PPTX/KT-Law-beamer1.pdf" target="_blank" class="download-btn btn-cooper">📄 Download PDF</a>
  </div>

  <div class="lecture-card card-cooper">
    <h3 class="lecture-title">📊 Beamer – From the Cooper Pair Problem to the BCS Ground State</h3>
    <p class="lecture-desc">
      <strong>Type:</strong> Presentation slides<br><br>
      Bridge from single‑pair Cooper problem to full BCS many‑body state: reduced BCS Hamiltonian, variational wavefunction \( |\Psi_{\mathrm{BCS}}\rangle = \prod_k (u_k + v_k c_{k\uparrow}^\dagger c_{-k\downarrow}^\dagger)|0\rangle \), energy minimization, gap parameter, coherence factors, pair amplitude \( g_k \), gap equation, and condensation energy.
    </p>
    <a href="https://maggiexheuw.github.io/PPTX/KT-Law-beamer2.pdf" target="_blank" class="download-btn btn-cooper">📄 Download PDF</a>
  </div>

  <div class="lecture-card card-cooper">
    <h3 class="lecture-title">📘 Lecture Note – The Cooper Problem (Concise)</h3>
    <p class="lecture-desc">
      <strong>Type:</strong> Lecture notes<br><br>
      Concise written derivation of the Cooper problem: motivation, trial wavefunction, amplitude equation, model interaction, bound‑state energy, and physical significance.
    </p>
    <a href="https://maggiexheuw.github.io/PPTX/KT-Law1.pdf" target="_blank" class="download-btn btn-cooper">📄 Download PDF</a>
  </div>

  <div class="lecture-card card-cooper">
    <h3 class="lecture-title">📘 Lecture Note – Cooper Problem &amp; BCS Theory (Full)</h3>
    <p class="lecture-desc">
      <strong>Type:</strong> Comprehensive lecture notes<br><br>
      Complete written version covering both Cooper problem and full BCS variational treatment: two‑body Schrödinger equation, model interaction, bound‑state energy, reduced BCS Hamiltonian, BCS trial wavefunction, energy minimization, coherence factors, pair amplitude, gap equation, and condensation energy.
    </p>
    <a href="https://maggiexheuw.github.io/PPTX/KT-Law2.pdf" target="_blank" class="download-btn btn-cooper">📄 Download PDF</a>
  </div>

</div>

<!-- 3.2 Quantum Geometry & Flat-Band SC -->
<h3 class="subsection-header">III-B. Quantum Geometry &amp; Flat‑Band Superconductivity</h3>

<div class="section-intro" style="border-left-color:#8e44ad;">
  How quantum geometry — the quantum metric of Bloch wavefunctions — sets the Cooper pair size in flat bands where \( v_F \to 0 \).
</div>

<div class="lecture-list">

  <div class="lecture-card card-geometry">
    <h3 class="lecture-title">📊 Beamer – The Size of a Cooper Pair in a Flat Band</h3>
    <p class="lecture-desc">
      <strong>Type:</strong> Presentation slides<br><br>
      Resolves the paradox of flat‑band superconductivity: introduces the quantum metric \( g_{\mu\nu}(\mathbf{k}) \), derives \( \xi_{\mathrm{pair}}^2 = 2\langle \operatorname{tr} g \rangle_{\mathrm{BZ}} \), and proves the topological lower bound \( \xi_{\mathrm{pair}}^2 \gtrsim |C| a^2/\pi \) for Chern bands.
    </p>
    <a href="https://maggiexheuw.github.io/PPTX/KT-Law-beamer3.pdf" target="_blank" class="download-btn btn-geometry">📄 Download PDF</a>
  </div>

  <div class="lecture-card card-geometry">
    <h3 class="lecture-title">📘 Lecture Note – The Size of a Cooper Pair in a Flat Band</h3>
    <p class="lecture-desc">
      <strong>Type:</strong> Comprehensive lecture notes<br><br>
      Self‑contained pedagogical introduction: conventional BCS coherence length, flat‑band Cooper problem with band projection, binding energy \( E_b = -U/N_{\mathrm{orb}} \), quantum geometric tensor \( Q_{\mu\nu} = g_{\mu\nu} - \frac{i}{2}\Omega_{\mu\nu} \), pair size \( \xi_{\mathrm{pair}}^2 = 2\langle \operatorname{tr} g \rangle_{\mathrm{BZ}} \), topological lower bound, and equivalence to Wannier spread.
    </p>
    <a href="https://maggiexheuw.github.io/PPTX/KT-Lwa3.pdf" target="_blank" class="download-btn btn-geometry">📄 Download PDF</a>
  </div>

</div>

<!-- 3.3 Fermi-Velocity Length Scales -->
<h3 class="subsection-header">III-C. Fermi‑Velocity‑Controlled Length Scales in Solids</h3>

<div class="section-intro" style="border-left-color:#d35400;">
  The master formula \( L_E \sim \hbar v_F / E \) unifying thermal length, BCS coherence, Kondo cloud, mean free path, and localization.
</div>

<div class="lecture-list">

  <div class="lecture-card card-length">
    <h3 class="lecture-title">📊 Beamer – Fermi-Velocity-Controlled Length Scales in Solids</h3>
    <p class="lecture-desc">
      <strong>Type:</strong> Presentation slides<br><br>
      Unified framework: \( L_T = \hbar v_F/2\pi k_B T \), \( \xi_0 = \hbar v_F/\pi\Delta \), \( \xi_K = \hbar v_F/k_B T_K \), \( \ell = v_F\tau \), diffusive generalization \( L_E \to \sqrt{\hbar D/E} \), and Thouless energy \( E_{\mathrm{Th}} = \hbar D/L^2 \).
    </p>
    <a href="https://maggiexheuw.github.io/PPTX/KT-Law-beamer4.pdf" target="_blank" class="download-btn btn-length">📄 Download PDF</a>
  </div>

  <div class="lecture-card card-length">
    <h3 class="lecture-title">📘 Lecture Note – Fermi-Velocity-Controlled Length Scales in Solids</h3>
    <p class="lecture-desc">
      <strong>Type:</strong> Comprehensive lecture notes<br><br>
      Full derivations: linearization at Fermi surface, three readings of \( L_E = \hbar v_F/E \), thermal length via Matsubara poles, SNS Josephson decay, BCS coherence from pair wavefunction, Kondo cloud from Yosida variational state, mean free path and 1D localization, diffusive generalization, and Thouless energy.
    </p>
    <a href="https://maggiexheuw.github.io/PPTX/KT-Law-beamer4%20Lecture%20Note.pdf" target="_blank" class="download-btn btn-length">📄 Download PDF</a>
  </div>

</div>

<!-- 3.4 Friedel Oscillations -->
<h3 class="subsection-header">III-D. Friedel Oscillations in One Dimension</h3>

<div class="section-intro" style="border-left-color:#27ae60;">
  Green's function / \( T \)-matrix derivation of the \( 2k_F \) density oscillation and its thermal cutoff.
</div>

<div class="lecture-list">

  <div class="lecture-card card-friedel">
    <h3 class="lecture-title">📊 Beamer – Friedel Oscillations in One Dimension</h3>
    <p class="lecture-desc">
      <strong>Type:</strong> Presentation slides<br><br>
      Green's function derivation: Lehmann representation, free Green's function \( G_0^R(x,x';E) = -\dfrac{i}{\hbar v_E} e^{ik_E|x-x'|} \), Dyson equation, separable \( T \)-matrix \( t(E) = u/(1 + iu/\hbar v_E) \), zero‑temperature result \( \delta n(x) \sim -\dfrac{u\nu_F}{2}\dfrac{\cos(2k_F x)}{|x|} \), and finite‑temperature thermal length \( \xi_T = \hbar v_F/2\pi k_B T \).
    </p>
    <a href="https://maggiexheuw.github.io/PPTX/Friedel.pdf" target="_blank" class="download-btn btn-friedel">📄 Download PDF</a>
  </div>

  <div class="lecture-card card-friedel">
    <h3 class="lecture-title">📘 Lecture Note – Friedel Oscillations in One Dimension</h3>
    <p class="lecture-desc">
      <strong>Type:</strong> Detailed lecture notes<br><br>
      Six‑section written note: physical setup, Lehmann representation and spectral function, free Green's function, Dyson equation and separable \( T \)-matrix, charge density from Green's function, zero‑temperature result, and finite‑temperature thermal envelope \( (x/\xi_T)/\sinh(x/\xi_T) \).
    </p>
    <a href="https://maggiexheuw.github.io/PPTX/KT-Law-beamer5.pdf" target="_blank" class="download-btn btn-friedel">📄 Download PDF</a>
  </div>

</div>

<!-- 3.5 Anderson Localization -->
<h3 class="subsection-header">III-E. Anderson Localization in One Dimension</h3>

<div class="section-intro" style="border-left-color:#c0392b;">
  Transfer‑matrix derivation of exponential localization, Lyapunov exponent, and the Thouless formula.
</div>

<div class="lecture-list">

  <div class="lecture-card card-local">
    <h3 class="lecture-title">📊 Beamer – Anderson Localization in One Dimension</h3>
    <p class="lecture-desc">
      <strong>Type:</strong> Presentation slides<br><br>
      1D Anderson model: tight‑binding Hamiltonian, three‑term recurrence, transfer matrix \( T_n \in \mathrm{SL}(2,\mathbb{R}) \), Lyapunov exponent \( \gamma(E) = \lim_{N\to\infty} \frac{1}{N}\ln\|\mathcal{M}_N\| \), localization length \( \xi(E) = 1/\gamma(E) \), clean‑chain limit, weak‑disorder Thouless formula \( \xi(E) = 2(4t^2 - E^2)/\sigma^2 \), Prüfer variables, and numerical QR re‑orthogonalization.
    </p>
    <a href="https://maggiexheuw.github.io/PPTX/Anderrson-loc.pdf" target="_blank" class="download-btn btn-local">📄 Download PDF</a>
  </div>

  <div class="lecture-card card-local">
    <h3 class="lecture-title">📘 Lecture Note – Anderson Localization in One Dimension</h3>
    <p class="lecture-desc">
      <strong>Type:</strong> Detailed lecture notes<br><br>
      Complete written derivation: tight‑binding Hamiltonian, transfer matrix, Furstenberg and Oseledec theorems, clean‑chain limit, weak‑disorder expansion in Prüfer variables, Thouless formula, numerical implementation, and connection to scaling theory of localization.
    </p>
    <a href="https://maggiexheuw.github.io/PPTX/KT-Law-beamer7.pdf" target="_blank" class="download-btn btn-local">📄 Download PDF</a>
  </div>

</div>

<!-- 3.6 Solid-State Length Scale Overview -->
<h3 class="subsection-header">III-F. Solid‑State Length Scale Overview</h3>

<div class="section-intro" style="border-left-color:#16a085;">
  A master overview of the hierarchy of length scales in solids, unifying the entire course.
</div>

<div class="lecture-list">

  <div class="lecture-card card-scale">
    <h3 class="lecture-title">📊 Beamer – Solid-State Length Scale Overview</h3>
    <p class="lecture-desc">
      <strong>Type:</strong> Presentation slides<br><br>
      High‑level summary: hierarchy of length scales from \( \lambda_F \) to \( \xi_K \), diffusive generalization, Thouless energy as the inverse dictionary, and limits of the framework (flat bands, non‑Fermi liquids, Dirac materials).
    </p>
    <a href="https://maggiexheuw.github.io/PPTX/Solid-state.pdf" target="_blank" class="download-btn btn-scale">📄 Download PDF</a>
  </div>

</div>