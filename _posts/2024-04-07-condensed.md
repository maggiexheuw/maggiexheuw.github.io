---
layout: post
title: Condense matter theory
subtitle: Condense matter theory note
date: 2024-04-07
author: Maggie
header-img: img/Andersonphoto.jpg
catalog: true
tags:
  - Condense matter theory
---

<style>
  /* ===== 笔记页面专用样式 ===== */
  .post-content {
    font-family: "Times New Roman", Times, serif;
    font-size: 1.2rem;
    line-height: 1.8;
    color: #1e2a3a;
    max-width: 900px;
    margin: 0 auto;
    padding: 0 20px;
  }

  .post-content h1,
  .post-content h2,
  .post-content h3,
  .post-content h4 {
    font-weight: 700;
    color: #1e3c6f;
    margin-top: 2.2rem;
    margin-bottom: 1rem;
  }

  .post-content h1 {
    font-size: 2.8rem;
    border-bottom: 4px solid #eef3f9;
    padding-bottom: 0.3rem;
  }

  .post-content h2 {
    font-size: 2.2rem;
    border-left: 6px solid #4a6fa5;
    padding-left: 1rem;
  }

  .post-content h3 {
    font-size: 1.8rem;
    margin-top: 1.8rem;
  }

  .post-content h4 {
    font-size: 1.5rem;
    margin-top: 1.5rem;
  }

  /* ===== 笔记列表（PDF链接） ===== */
  .notes-grid {
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
    gap: 1rem;
    margin: 1.5rem 0 2rem 0;
  }

  .note-card {
    background: #f8faff;
    border-radius: 16px;
    padding: 1.2rem 1.4rem;
    text-align: center;
    border-left: 5px solid #4a6fa5;
    transition: all 0.25s ease;
    box-shadow: 0 2px 8px rgba(0, 0, 0, 0.04);
  }

  .note-card:hover {
    background: #f0f6ff;
    transform: translateY(-4px);
    box-shadow: 0 8px 24px rgba(74, 111, 165, 0.15);
  }

  .note-card a {
    font-size: 1.3rem;
    font-weight: 700;
    color: #1e3c6f;
    text-decoration: none;
    display: block;
  }

  .note-card a:hover {
    color: #4a6fa5;
    text-decoration: underline;
  }

  /* ===== 链接列表（有用链接） ===== */
  .link-section {
    margin: 2rem 0;
    padding: 1.2rem 1.8rem;
    background: #ffffff;
    border-radius: 20px;
    box-shadow: 0 4px 20px rgba(0, 0, 0, 0.06);
    border: 1px solid rgba(255, 255, 255, 0.6);
    transition: box-shadow 0.3s;
  }

  .link-section:hover {
    box-shadow: 0 12px 40px rgba(26, 58, 107, 0.10);
  }

  .link-section h4 {
    font-size: 1.6rem;
    color: #1e3c6f;
    border-bottom: 2px solid #eef3f9;
    padding-bottom: 0.3rem;
    margin-top: 0.2rem;
  }

  .link-list {
    list-style: none;
    padding-left: 0;
    display: flex;
    flex-wrap: wrap;
    gap: 0.5rem 1.2rem;
    margin: 0.8rem 0 0.2rem 0;
  }

  .link-list li {
    margin: 0;
  }

  .link-list li a {
    display: inline-block;
    padding: 0.2rem 1rem;
    background: #f0f4fc;
    border-radius: 40px;
    font-size: 1.1rem;
    color: #1e3c6f;
    text-decoration: none;
    font-weight: 500;
    transition: all 0.25s ease;
    border: 1px solid transparent;
  }

  .link-list li a:hover {
    background: #4a6fa5;
    color: #fff;
    transform: translateY(-2px);
    box-shadow: 0 4px 16px rgba(74, 111, 165, 0.25);
    border-color: #4a6fa5;
  }

  /* ===== 特殊标题装饰 ===== */
  .section-title {
    font-size: 2.0rem;
    font-weight: 700;
    color: #1e3c6f;
    margin-top: 2.5rem;
    display: flex;
    align-items: center;
    gap: 0.8rem;
  }

  .section-title::before {
    content: "📚";
    font-size: 2.2rem;
  }

  .section-title.links-title::before {
    content: "🔗";
  }

  /* 响应式 */
  @media (max-width: 768px) {
    .post-content {
      font-size: 1.1rem;
    }
    .post-content h1 {
      font-size: 2.2rem;
    }
    .post-content h2 {
      font-size: 1.8rem;
    }
    .post-content h3 {
      font-size: 1.5rem;
    }
    .notes-grid {
      grid-template-columns: 1fr 1fr;
    }
    .note-card a {
      font-size: 1.1rem;
    }
    .link-list li a {
      font-size: 0.95rem;
      padding: 0.15rem 0.8rem;
    }
  }

  @media (max-width: 480px) {
    .notes-grid {
      grid-template-columns: 1fr;
    }
    .link-list {
      flex-direction: column;
      align-items: flex-start;
    }
    .link-list li a {
      width: 100%;
      text-align: center;
    }
  }
</style>

<!-- ===== 内容开始 ===== -->
<div class="post-content">

  <h1>Condense matter physics notes</h1>

  <!-- PDF 笔记卡片网格 -->
  <div class="notes-grid">
    <div class="note-card">
      <a href="https://maggiexheuw.github.io/pdf/condense-19-26.pdf">Fermi liquid</a>
    </div>
    <div class="note-card">
      <a href="https://maggiexheuw.github.io/pdf/condense-27-41.pdf">Superconductivity</a>
    </div>
    <div class="note-card">
      <a href="https://maggiexheuw.github.io/pdf/condense-43-53.pdf">Green function</a>
    </div>
    <div class="note-card">
      <a href="https://maggiexheuw.github.io/pdf/condense-71-77.pdf">KT transition</a>
    </div>
  </div>

  <!-- ===== 有用链接 ===== -->
  <h2 class="section-title links-title">Some useful links</h2>

  <!-- Condense matter physics -->
  <div class="link-section">
    <h4>Condense matter physics</h4>
    <ul class="link-list">
      <li><a href="https://mareknarozniak.com/tagged/#qutip">Marek narozniak's homepage</a></li>
      <li><a href="https://www.openmx-square.org/">Open MX</a></li>
      <li><a href="https://verse-and-dimensions.fandom.com/wiki/Quaternionic_projective_line">Quaternionic projective line</a></li>
      <li><a href="https://bingweb.binghamton.edu/~suzuki/SolidStatePhysicsLN.html">Magnetism</a></li>
      <li><a href="http://txiang.iphy.ac.cn/">Taoxiang</a></li>
    </ul>
  </div>

  <!-- Group theory -->
  <div class="link-section">
    <h4>Group theory</h4>
    <ul class="link-list">
      <li><a href="https://esackinger.wordpress.com/blog/lie-groups-and-their-representations/">Group theory in physics</a></li>
      <li><a href="https://math.mit.edu/classes/18.745/classnotes.html">Lie group lectures</a></li>
      <li><a href="https://edu.itp.phys.ethz.ch/hs13/Symmetries/">Symmetry in physics</a></li>
    </ul>
  </div>

  <!-- Topological phases -->
  <div class="link-section">
    <h4>Topological phases</h4>
    <ul class="link-list">
      <li><a href="https://journals.aps.org/rmp/abstract/10.1103/RevModPhys.88.035001">Witten (2016)</a></li>
    </ul>
  </div>

  <!-- Field theory -->
  <div class="link-section">
    <h4>Field theory</h4>
    <ul class="link-list">
      <li><a href="https://archive.int.washington.edu/users/dbkaplan/">Kaplan homepage</a></li>
      <li><a href="https://archive.int.washington.edu/users/dbkaplan/572_16/">Kaplan field theory</a></li>
      <li><a href="https://eduardo.physics.illinois.edu/phys583/ch20.pdf">Fradkin chiral anomaly</a></li>
      <li><a href="https://www.icts.res.in/sites/default/files/comprty2018-2018-02-02-David-B-Kaplan-3.pdf">The mystery of chiral gauge theories (David Kaplan)</a></li>
      <li><a href="http://fma.if.usp.br/~burdman/QFT2/qft2index.html">Gustavo Burdman</a></li>
      <li><a href="https://indico.cern.ch/event/615296/contributions/2606797/attachments/1470593/2275403/Chiral_Anomalies.pdf">Chiral anomaly lecture</a></li>
      <li><a href="https://www.math.umd.edu/~tadmor/ki_net/activities/presentations/">Field theory talks</a></li>
    </ul>
  </div>

  <!-- Atiyah Singer index theorem -->
  <div class="link-section">
    <h4>Atiyah Singer index theorem</h4>
    <ul class="link-list">
      <li><a href="https://www3.nd.edu/~lnicolae/ind-thm.pdf">Nicolaescu note</a></li>
      <li><a href="https://www3.nd.edu/~lnicolae/">Nicolaescu homepage</a></li>
    </ul>
  </div>

  <!-- Teaching -->
  <div class="link-section">
    <h4>Teaching</h4>
    <ul class="link-list">
      <li><a href="https://phas.ubc.ca/~stamp/TEACHING/">Quantum mechanics</a></li>
      <li><a href="https://www.math.uci.edu/~brusso/">Bernard Russo (Lie algebra)</a></li>
      <li><a href="http://www.math.columbia.edu/~woit/">Peter Woit (Lie algebra)</a></li>
      <li><a href="http://staff.ustc.edu.cn/~wangzuoq/Courses/">Zuoqin Wang (Lie algebra)</a></li>
      <li><a href="https://bingweb.binghamton.edu/~suzuki/SolidStatePhysicsLN.html">Suzuki's homepage (Solid physics)</a></li>
    </ul>
  </div>

</div>
<!-- ===== 内容结束 ===== -->