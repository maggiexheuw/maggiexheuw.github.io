---
layout: post
title: "CCMP Reports — Condensed Matter Physics Conference 2026"
subtitle: "Collected lecture notes from Liyang CCMP 2026"
date: 2026-07-22 00:00:00 +0800
author: Maggie
header-img: img/EdWitten.jpg
catalog: true
mathjax: true
tags:
  - CCMP 2026
  - condensed matter physics
  - high-Tc superconductivity
  - topological phases
  - quantum many-body
---
<script>
(function () {
  window.MathJax = window.MathJax || {
    tex: {
      inlineMath: [['\\(', '\\)']],
      displayMath: [['\\[', '\\]']],
      processEscapes: true
    },
    options: {
      skipHtmlTags: ['script', 'noscript', 'style', 'textarea', 'pre', 'code']
    }
  };

  if (!document.querySelector('script[src*="mathjax"]')) {
    var mathJaxScript = document.createElement('script');
    mathJaxScript.defer = true;
    mathJaxScript.src = 'https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js';
    document.head.appendChild(mathJaxScript);
  }
}());
</script>


<style>
/*
  CCMP 2026 editorial page
  All selectors are namespaced under .ccmp-page to avoid conflicts with the blog theme.
*/
.ccmp-page,
.ccmp-page * {
  box-sizing: border-box;
}

.ccmp-page {
  --ccmp-paper: #f3f0e8;
  --ccmp-paper-deep: #e8e2d5;
  --ccmp-card: rgba(255, 253, 248, 0.94);
  --ccmp-ink: #16211f;
  --ccmp-muted: #68716d;
  --ccmp-line: rgba(22, 33, 31, 0.14);
  --ccmp-accent: #a53a2c;
  --ccmp-accent-deep: #76251b;
  --ccmp-green: #1f5b4b;
  --ccmp-navy: #102c37;
  --ccmp-gold: #b58c45;
  --ccmp-shadow: 0 20px 60px rgba(31, 42, 38, 0.12);
  --ccmp-serif: "Iowan Old Style", "Palatino Linotype", "Book Antiqua", "Noto Serif", Georgia, serif;
  --ccmp-sans: Inter, ui-sans-serif, -apple-system, BlinkMacSystemFont, "Segoe UI", "Noto Sans", Arial, sans-serif;

  position: relative;
  left: 50%;
  width: 100vw;
  margin-left: -50vw;
  overflow: hidden;
  color: var(--ccmp-ink);
  background:
    radial-gradient(circle at 10% 4%, rgba(165, 58, 44, 0.08), transparent 25rem),
    radial-gradient(circle at 92% 28%, rgba(31, 91, 75, 0.08), transparent 28rem),
    var(--ccmp-paper);
  font-family: var(--ccmp-sans);
  font-size: 17px;
  line-height: 1.72;
  isolation: isolate;
}

.ccmp-page[data-ccmp-theme="dark"] {
  --ccmp-paper: #111816;
  --ccmp-paper-deep: #19221f;
  --ccmp-card: rgba(24, 33, 30, 0.94);
  --ccmp-ink: #f0ede4;
  --ccmp-muted: #b2bbb6;
  --ccmp-line: rgba(240, 237, 228, 0.14);
  --ccmp-accent: #ef8577;
  --ccmp-accent-deep: #f2a196;
  --ccmp-green: #73c4a7;
  --ccmp-navy: #0c2028;
  --ccmp-gold: #d4b475;
  --ccmp-shadow: 0 22px 70px rgba(0, 0, 0, 0.32);
}

.ccmp-page button,
.ccmp-page input {
  font: inherit;
}

.ccmp-page a {
  color: inherit;
}

.ccmp-page h1,
.ccmp-page h2,
.ccmp-page h3,
.ccmp-page p {
  max-width: none;
}

.ccmp-page [hidden] {
  display: none !important;
}

.ccmp-skip-link {
  position: absolute;
  top: 12px;
  left: 12px;
  z-index: 100;
  padding: 10px 14px;
  border-radius: 8px;
  background: #fff;
  color: #111;
  transform: translateY(-160%);
}

.ccmp-skip-link:focus {
  transform: translateY(0);
}

.ccmp-hero {
  position: relative;
  min-height: 650px;
  display: grid;
  align-items: end;
  color: #f7f3e9;
  background:
    linear-gradient(115deg, rgba(8, 31, 37, 0.98), rgba(15, 50, 54, 0.90) 48%, rgba(76, 40, 29, 0.90)),
    var(--ccmp-navy);
  border-bottom: 1px solid rgba(255, 255, 255, 0.14);
}

.ccmp-hero::before,
.ccmp-hero::after {
  content: "";
  position: absolute;
  pointer-events: none;
}

.ccmp-hero::before {
  inset: 0;
  opacity: 0.35;
  background-image:
    linear-gradient(rgba(255, 255, 255, 0.045) 1px, transparent 1px),
    linear-gradient(90deg, rgba(255, 255, 255, 0.045) 1px, transparent 1px);
  background-size: 46px 46px;
  mask-image: linear-gradient(to bottom, black, transparent 85%);
}

.ccmp-hero::after {
  width: min(62vw, 820px);
  aspect-ratio: 1;
  right: -17vw;
  top: -54%;
  border: 1px solid rgba(255, 255, 255, 0.18);
  border-radius: 50%;
  box-shadow:
    0 0 0 74px rgba(255, 255, 255, 0.025),
    0 0 0 148px rgba(255, 255, 255, 0.018),
    0 0 0 222px rgba(255, 255, 255, 0.012);
}

.ccmp-hero-inner,
.ccmp-main,
.ccmp-footer-inner {
  width: min(1320px, calc(100% - 48px));
  margin-inline: auto;
}

.ccmp-hero-inner {
  position: relative;
  z-index: 2;
  padding: 30px 0 66px;
}

.ccmp-masthead {
  display: flex;
  justify-content: space-between;
  align-items: center;
  gap: 24px;
  padding: 18px 0;
  border-bottom: 1px solid rgba(255, 255, 255, 0.18);
  font-size: 0.72rem;
  font-weight: 750;
  letter-spacing: 0.16em;
  text-transform: uppercase;
}

.ccmp-brand {
  display: flex;
  align-items: center;
  gap: 12px;
}

.ccmp-brand-mark {
  display: grid;
  place-items: center;
  width: 34px;
  height: 34px;
  border: 1px solid rgba(255, 255, 255, 0.36);
  border-radius: 50%;
  font-family: var(--ccmp-serif);
  font-size: 1rem;
  letter-spacing: 0;
}

.ccmp-theme-toggle {
  display: inline-flex;
  align-items: center;
  gap: 8px;
  padding: 8px 12px;
  border: 1px solid rgba(255, 255, 255, 0.28);
  border-radius: 999px;
  color: inherit;
  background: rgba(255, 255, 255, 0.08);
  cursor: pointer;
}

.ccmp-theme-toggle:hover,
.ccmp-theme-toggle:focus-visible {
  background: rgba(255, 255, 255, 0.16);
}

.ccmp-hero-copy {
  display: grid;
  grid-template-columns: minmax(0, 1.55fr) minmax(280px, 0.55fr);
  gap: clamp(36px, 8vw, 110px);
  align-items: end;
  padding-top: clamp(78px, 13vw, 145px);
}

.ccmp-kicker {
  margin: 0 0 18px;
  color: #d6c6a8;
  font-size: 0.76rem;
  font-weight: 800;
  letter-spacing: 0.17em;
  text-transform: uppercase;
}

.ccmp-title {
  margin: 0;
  font-family: var(--ccmp-serif);
  font-size: clamp(4.2rem, 10.5vw, 9.2rem);
  font-weight: 600;
  line-height: 0.79;
  letter-spacing: -0.065em;
}

.ccmp-title span {
  display: block;
  margin-left: 0.52em;
  color: #e9dcc4;
  font-style: italic;
  font-weight: 400;
}

.ccmp-hero-summary {
  margin: 0;
  color: rgba(247, 243, 233, 0.78);
  font-family: var(--ccmp-serif);
  font-size: clamp(1.02rem, 1.6vw, 1.27rem);
  line-height: 1.72;
}

.ccmp-hero-meta {
  display: grid;
  gap: 12px;
  margin-top: 28px;
  padding-top: 24px;
  border-top: 1px solid rgba(255, 255, 255, 0.18);
  font-size: 0.78rem;
  letter-spacing: 0.04em;
}

.ccmp-hero-meta div {
  display: flex;
  justify-content: space-between;
  gap: 20px;
}

.ccmp-hero-meta dt {
  color: rgba(247, 243, 233, 0.58);
}

.ccmp-hero-meta dd {
  margin: 0;
  text-align: right;
}

.ccmp-main {
  position: relative;
  z-index: 3;
  padding: 64px 0 96px;
}

.ccmp-intro {
  display: grid;
  grid-template-columns: minmax(0, 1.4fr) minmax(290px, 0.6fr);
  gap: 28px;
  margin-bottom: 28px;
}

.ccmp-editorial,
.ccmp-index,
.ccmp-tools,
.ccmp-card,
.ccmp-modal-panel {
  border: 1px solid var(--ccmp-line);
  background: var(--ccmp-card);
  box-shadow: var(--ccmp-shadow);
}

.ccmp-editorial,
.ccmp-index {
  border-radius: 24px;
}

.ccmp-editorial {
  position: relative;
  padding: clamp(30px, 5vw, 58px);
  overflow: hidden;
}

.ccmp-editorial::after {
  content: "∮";
  position: absolute;
  right: 22px;
  bottom: -52px;
  color: var(--ccmp-accent);
  opacity: 0.08;
  font-family: var(--ccmp-serif);
  font-size: 11rem;
  line-height: 1;
}

.ccmp-section-label {
  display: flex;
  align-items: center;
  gap: 10px;
  margin-bottom: 16px;
  color: var(--ccmp-accent);
  font-size: 0.72rem;
  font-weight: 850;
  letter-spacing: 0.16em;
  text-transform: uppercase;
}

.ccmp-section-label::before {
  content: "";
  width: 28px;
  height: 1px;
  background: currentColor;
}

.ccmp-editorial h2 {
  position: relative;
  z-index: 1;
  margin: 0 0 18px;
  font-family: var(--ccmp-serif);
  font-size: clamp(2rem, 4vw, 3.25rem);
  font-weight: 600;
  line-height: 1.08;
  letter-spacing: -0.035em;
}

.ccmp-editorial p {
  position: relative;
  z-index: 1;
  margin: 0;
  color: var(--ccmp-muted);
  font-family: var(--ccmp-serif);
  font-size: 1.08rem;
}

.ccmp-stats {
  position: relative;
  z-index: 1;
  display: grid;
  grid-template-columns: repeat(3, minmax(0, 1fr));
  gap: 14px;
  margin-top: 34px;
}

.ccmp-stat {
  padding-top: 16px;
  border-top: 1px solid var(--ccmp-line);
}

.ccmp-stat strong {
  display: block;
  font-family: var(--ccmp-serif);
  font-size: 2rem;
  line-height: 1;
}

.ccmp-stat span {
  color: var(--ccmp-muted);
  font-size: 0.72rem;
  letter-spacing: 0.08em;
  text-transform: uppercase;
}

.ccmp-index {
  padding: 30px;
}

.ccmp-index h2 {
  margin: 0 0 18px;
  font-family: var(--ccmp-serif);
  font-size: 1.55rem;
  font-weight: 600;
}

.ccmp-index-list {
  display: grid;
  gap: 3px;
  max-height: 385px;
  overflow: auto;
  padding-right: 8px;
  scrollbar-width: thin;
}

.ccmp-index-list a {
  display: grid;
  grid-template-columns: 34px 1fr;
  gap: 10px;
  padding: 9px 8px;
  border-radius: 10px;
  color: var(--ccmp-muted);
  text-decoration: none;
  font-size: 0.76rem;
  line-height: 1.42;
}

.ccmp-index-list a:hover,
.ccmp-index-list a:focus-visible {
  color: var(--ccmp-ink);
  background: rgba(165, 58, 44, 0.08);
}

.ccmp-index-list b {
  color: var(--ccmp-accent);
  font-variant-numeric: tabular-nums;
}

.ccmp-tools {
  position: sticky;
  top: 14px;
  z-index: 20;
  display: grid;
  grid-template-columns: minmax(240px, 1fr) auto;
  gap: 14px 20px;
  align-items: center;
  margin-bottom: 34px;
  padding: 16px;
  border-radius: 18px;
  backdrop-filter: blur(18px);
}

.ccmp-search {
  position: relative;
}

.ccmp-search::before {
  content: "⌕";
  position: absolute;
  left: 17px;
  top: 50%;
  transform: translateY(-52%);
  color: var(--ccmp-muted);
  font-family: var(--ccmp-serif);
  font-size: 1.35rem;
}

.ccmp-search input {
  width: 100%;
  height: 50px;
  padding: 0 18px 0 46px;
  border: 1px solid var(--ccmp-line);
  border-radius: 999px;
  outline: none;
  color: var(--ccmp-ink);
  background: var(--ccmp-paper);
}

.ccmp-search input:focus {
  border-color: var(--ccmp-accent);
  box-shadow: 0 0 0 4px rgba(165, 58, 44, 0.11);
}

.ccmp-filters {
  display: flex;
  flex-wrap: wrap;
  gap: 8px;
  justify-content: flex-end;
}

.ccmp-filter,
.ccmp-tag {
  border: 1px solid var(--ccmp-line);
  border-radius: 999px;
  color: var(--ccmp-muted);
  background: transparent;
  cursor: pointer;
}

.ccmp-filter {
  padding: 9px 13px;
  font-size: 0.72rem;
  font-weight: 750;
}

.ccmp-filter:hover,
.ccmp-filter:focus-visible,
.ccmp-filter.is-active {
  border-color: var(--ccmp-ink);
  color: var(--ccmp-paper);
  background: var(--ccmp-ink);
}

.ccmp-result-count {
  grid-column: 1 / -1;
  margin: -3px 8px 0;
  color: var(--ccmp-muted);
  font-size: 0.72rem;
}

.ccmp-grid {
  display: grid;
  grid-template-columns: repeat(2, minmax(0, 1fr));
  gap: 30px;
}

.ccmp-card {
  --card-accent: var(--ccmp-accent);
  display: grid;
  grid-template-rows: 205px 1fr;
  min-width: 0;
  border-radius: 24px;
  overflow: hidden;
  scroll-margin-top: 110px;
  transition: transform 220ms ease, box-shadow 220ms ease, border-color 220ms ease;
}

.ccmp-card:nth-child(3n + 2) { --card-accent: var(--ccmp-green); }
.ccmp-card:nth-child(3n) { --card-accent: var(--ccmp-gold); }

.ccmp-card:hover {
  transform: translateY(-5px);
  border-color: color-mix(in srgb, var(--card-accent), transparent 58%);
  box-shadow: 0 28px 76px rgba(31, 42, 38, 0.17);
}

.ccmp-card.is-hidden {
  display: none;
}

.ccmp-cover {
  position: relative;
  display: flex;
  flex-direction: column;
  justify-content: space-between;
  min-height: 0;
  padding: 26px;
  overflow: hidden;
  color: #f8f3e8;
  background: linear-gradient(135deg, #173c3e, #10292f 68%, #1b2421);
  background:
    linear-gradient(135deg, color-mix(in srgb, var(--card-accent), #0c2428 72%), #10292f 68%, #1b2421);
}

.ccmp-cover::before {
  content: "";
  position: absolute;
  inset: 0;
  opacity: 0.28;
  background-image:
    linear-gradient(rgba(255,255,255,0.08) 1px, transparent 1px),
    linear-gradient(90deg, rgba(255,255,255,0.08) 1px, transparent 1px);
  background-size: 28px 28px;
  mask-image: linear-gradient(135deg, black, transparent 80%);
}

.ccmp-cover::after {
  content: "ψ";
  position: absolute;
  right: 20px;
  bottom: -42px;
  color: rgba(255, 255, 255, 0.09);
  font-family: var(--ccmp-serif);
  font-size: 10rem;
  line-height: 1;
}

.ccmp-cover-top,
.ccmp-cover-bottom {
  position: relative;
  z-index: 1;
}

.ccmp-cover-top {
  display: flex;
  justify-content: space-between;
  gap: 20px;
  align-items: start;
}

.ccmp-report-number {
  font-family: var(--ccmp-serif);
  font-size: 3.7rem;
  line-height: 0.9;
  letter-spacing: -0.06em;
}

.ccmp-report-series {
  max-width: 130px;
  color: rgba(248, 243, 232, 0.68);
  font-size: 0.63rem;
  font-weight: 800;
  letter-spacing: 0.16em;
  text-align: right;
  text-transform: uppercase;
}

.ccmp-cover-bottom {
  display: flex;
  justify-content: space-between;
  align-items: end;
  gap: 20px;
}

.ccmp-primary-topic {
  font-family: var(--ccmp-serif);
  font-size: 1.1rem;
}

.ccmp-cover-symbol {
  color: rgba(248, 243, 232, 0.62);
  font-family: var(--ccmp-serif);
  font-size: 1.3rem;
  letter-spacing: 0.08em;
}

.ccmp-card-body {
  display: flex;
  flex-direction: column;
  padding: 30px 30px 28px;
}

.ccmp-card h2 {
  margin: 0;
  font-family: var(--ccmp-serif);
  font-size: clamp(1.58rem, 2.7vw, 2.15rem);
  font-weight: 600;
  line-height: 1.14;
  letter-spacing: -0.026em;
}

.ccmp-byline {
  display: flex;
  align-items: center;
  gap: 10px;
  margin: 13px 0 0;
  color: var(--ccmp-muted);
  font-size: 0.8rem;
  font-weight: 750;
}

.ccmp-byline::before {
  content: "";
  width: 24px;
  height: 1px;
  background: var(--card-accent);
}

.ccmp-tags {
  display: flex;
  flex-wrap: wrap;
  gap: 7px;
  margin: 19px 0;
}

.ccmp-tag {
  padding: 5px 9px;
  font-size: 0.67rem;
  font-weight: 720;
}

.ccmp-tag:hover,
.ccmp-tag:focus-visible {
  border-color: var(--card-accent);
  color: var(--ccmp-ink);
  background: color-mix(in srgb, var(--card-accent), transparent 90%);
}

.ccmp-abstract {
  flex: 1;
  margin: 0;
  color: var(--ccmp-muted);
  font-family: var(--ccmp-serif);
  font-size: 0.98rem;
  line-height: 1.72;
}

.ccmp-actions {
  display: flex;
  flex-wrap: wrap;
  gap: 9px;
  margin-top: 25px;
}

.ccmp-button {
  display: inline-flex;
  align-items: center;
  justify-content: center;
  min-height: 42px;
  padding: 9px 14px;
  border: 1px solid var(--ccmp-line);
  border-radius: 999px;
  color: var(--ccmp-muted);
  background: transparent;
  text-decoration: none;
  font-size: 0.72rem;
  font-weight: 800;
  cursor: pointer;
}

.ccmp-button:hover,
.ccmp-button:focus-visible {
  transform: translateY(-1px);
  color: var(--ccmp-ink);
  border-color: var(--card-accent);
}

.ccmp-button-primary {
  border-color: var(--ccmp-ink);
  color: var(--ccmp-paper);
  background: var(--ccmp-ink);
}

.ccmp-button-primary:hover,
.ccmp-button-primary:focus-visible {
  color: var(--ccmp-paper);
  opacity: 0.88;
}

.ccmp-empty {
  padding: 86px 20px;
  color: var(--ccmp-muted);
  text-align: center;
}

.ccmp-empty strong {
  display: block;
  margin-bottom: 8px;
  color: var(--ccmp-ink);
  font-family: var(--ccmp-serif);
  font-size: 1.6rem;
}

.ccmp-footer {
  color: rgba(247, 243, 233, 0.78);
  background: #0e252b;
}

.ccmp-footer-inner {
  display: grid;
  grid-template-columns: 1fr auto;
  gap: 30px;
  align-items: end;
  padding: 46px 0 54px;
}

.ccmp-footer-title {
  margin: 0;
  color: #f7f3e9;
  font-family: var(--ccmp-serif);
  font-size: 1.8rem;
  font-weight: 500;
}

.ccmp-footer p {
  margin: 8px 0 0;
  font-size: 0.76rem;
}

.ccmp-back-top {
  color: inherit;
  font-size: 0.72rem;
  font-weight: 800;
  letter-spacing: 0.1em;
  text-decoration: none;
  text-transform: uppercase;
}

.ccmp-modal {
  position: fixed;
  inset: 0;
  z-index: 9999;
  display: grid;
  place-items: center;
  padding: 22px;
  background: rgba(4, 12, 14, 0.78);
  backdrop-filter: blur(10px);
}

.ccmp-modal-panel {
  width: min(1120px, 100%);
  height: min(86vh, 860px);
  display: grid;
  grid-template-rows: auto 1fr;
  border-radius: 20px;
  overflow: hidden;
}

.ccmp-modal-header {
  display: flex;
  justify-content: space-between;
  align-items: center;
  gap: 20px;
  padding: 14px 18px;
  border-bottom: 1px solid var(--ccmp-line);
}

.ccmp-modal-title {
  margin: 0;
  overflow: hidden;
  font-family: var(--ccmp-serif);
  font-size: 1rem;
  font-weight: 600;
  text-overflow: ellipsis;
  white-space: nowrap;
}

.ccmp-modal-close {
  flex: 0 0 auto;
  width: 38px;
  height: 38px;
  border: 1px solid var(--ccmp-line);
  border-radius: 50%;
  color: var(--ccmp-ink);
  background: transparent;
  cursor: pointer;
}

.ccmp-modal iframe {
  width: 100%;
  height: 100%;
  border: 0;
  background: #fff;
}

.ccmp-page mjx-container {
  overflow-x: auto;
  overflow-y: hidden;
  max-width: 100%;
}

@media (min-width: 1050px) {
  .ccmp-card:first-child {
    grid-column: 1 / -1;
    grid-template-columns: minmax(340px, 0.78fr) minmax(0, 1.22fr);
    grid-template-rows: 1fr;
  }

  .ccmp-card:first-child .ccmp-cover {
    min-height: 440px;
  }

  .ccmp-card:first-child .ccmp-card-body {
    padding: 42px;
  }

  .ccmp-card:first-child h2 {
    font-size: clamp(2.2rem, 4vw, 3.45rem);
  }

  .ccmp-card:first-child .ccmp-abstract {
    font-size: 1.05rem;
  }
}

@media (max-width: 960px) {
  .ccmp-hero-copy,
  .ccmp-intro {
    grid-template-columns: 1fr;
  }

  .ccmp-hero-copy {
    gap: 38px;
  }

  .ccmp-hero-summary {
    max-width: 680px;
  }

  .ccmp-index-list {
    max-height: 290px;
  }

  .ccmp-tools {
    position: relative;
    top: auto;
    grid-template-columns: 1fr;
  }

  .ccmp-filters {
    justify-content: flex-start;
  }

  .ccmp-grid {
    grid-template-columns: 1fr;
  }
}

@media (max-width: 620px) {
  .ccmp-page {
    font-size: 16px;
  }

  .ccmp-hero {
    min-height: 610px;
  }

  .ccmp-hero-inner,
  .ccmp-main,
  .ccmp-footer-inner {
    width: min(100% - 28px, 1320px);
  }

  .ccmp-masthead {
    letter-spacing: 0.1em;
  }

  .ccmp-theme-toggle span {
    display: none;
  }

  .ccmp-title {
    font-size: clamp(3.6rem, 21vw, 6rem);
  }

  .ccmp-title span {
    margin-left: 0.28em;
  }

  .ccmp-hero-copy {
    padding-top: 70px;
  }

  .ccmp-stats {
    gap: 9px;
  }

  .ccmp-editorial,
  .ccmp-index,
  .ccmp-card {
    border-radius: 18px;
  }

  .ccmp-tools {
    margin-inline: -2px;
    padding: 12px;
  }

  .ccmp-filters {
    flex-wrap: nowrap;
    justify-content: flex-start;
    overflow-x: auto;
    padding-bottom: 4px;
  }

  .ccmp-filter {
    flex: 0 0 auto;
  }

  .ccmp-cover {
    min-height: 185px;
  }

  .ccmp-card-body {
    padding: 25px 22px 24px;
  }

  .ccmp-actions {
    display: grid;
    grid-template-columns: 1fr 1fr;
  }

  .ccmp-button {
    width: 100%;
  }

  .ccmp-footer-inner {
    grid-template-columns: 1fr;
  }

  .ccmp-modal {
    padding: 8px;
  }

  .ccmp-modal-panel {
    height: 92vh;
    border-radius: 14px;
  }
}

@media (prefers-reduced-motion: reduce) {
  .ccmp-page *,
  .ccmp-page *::before,
  .ccmp-page *::after {
    scroll-behavior: auto !important;
    transition: none !important;
  }
}

/* ==========================================================
   SO(5) reduced CG calculator
   All selectors remain inside .ccmp-page / .ccmp-cg-lab.
   ========================================================== */
.ccmp-cg-lab {
  margin: 0 0 34px;
  padding: clamp(24px, 4vw, 42px);
  border: 1px solid var(--ccmp-line);
  border-radius: 26px;
  background:
    linear-gradient(145deg, color-mix(in srgb, var(--ccmp-card), transparent 2%), color-mix(in srgb, var(--ccmp-paper-deep), transparent 18%));
  box-shadow: var(--ccmp-shadow);
  scroll-margin-top: 110px;
}

.ccmp-cg-head {
  display: grid;
  grid-template-columns: minmax(0, 1fr) auto;
  gap: 24px;
  align-items: start;
  margin-bottom: 26px;
}

.ccmp-cg-head h2 {
  margin: 4px 0 10px;
  font-family: var(--ccmp-serif);
  font-size: clamp(2rem, 4vw, 3.1rem);
  font-weight: 600;
  line-height: 1.08;
  letter-spacing: -0.035em;
}

.ccmp-cg-head p {
  margin: 0;
  max-width: 850px;
  color: var(--ccmp-muted);
  font-family: var(--ccmp-serif);
}

.ccmp-cg-equation {
  min-width: 230px;
  padding: 15px 18px;
  border: 1px solid var(--ccmp-line);
  border-radius: 16px;
  color: var(--ccmp-ink);
  background: var(--ccmp-paper);
  font-family: var(--ccmp-serif);
  text-align: center;
}

.ccmp-cg-form {
  display: grid;
  gap: 22px;
}

.ccmp-cg-fieldset {
  min-width: 0;
  margin: 0;
  padding: 20px;
  border: 1px solid var(--ccmp-line);
  border-radius: 18px;
  background: color-mix(in srgb, var(--ccmp-card), transparent 10%);
}

.ccmp-cg-fieldset legend {
  padding: 0 8px;
  color: var(--ccmp-accent);
  font-size: 0.76rem;
  font-weight: 850;
  letter-spacing: 0.12em;
  text-transform: uppercase;
}

.ccmp-cg-fields {
  display: grid;
  grid-template-columns: repeat(4, minmax(120px, 1fr));
  gap: 14px;
}

.ccmp-cg-fields-3 {
  grid-template-columns: repeat(3, minmax(160px, 1fr));
}

.ccmp-cg-field {
  display: grid;
  gap: 6px;
  min-width: 0;
}

.ccmp-cg-field label,
.ccmp-cg-check-label {
  color: var(--ccmp-muted);
  font-size: 0.75rem;
  font-weight: 780;
}

.ccmp-cg-field input,
.ccmp-cg-field select {
  width: 100%;
  min-width: 0;
  height: 45px;
  padding: 0 12px;
  border: 1px solid var(--ccmp-line);
  border-radius: 11px;
  outline: none;
  color: var(--ccmp-ink);
  background: var(--ccmp-paper);
  font: inherit;
  font-variant-numeric: tabular-nums;
}

.ccmp-cg-field input:focus,
.ccmp-cg-field select:focus {
  border-color: var(--ccmp-accent);
  box-shadow: 0 0 0 4px color-mix(in srgb, var(--ccmp-accent), transparent 88%);
}

.ccmp-cg-help {
  color: var(--ccmp-muted);
  font-size: 0.68rem;
  line-height: 1.35;
}

.ccmp-cg-options {
  display: flex;
  flex-wrap: wrap;
  gap: 12px 22px;
  align-items: center;
  margin-top: 15px;
}

.ccmp-cg-check-label {
  display: inline-flex;
  align-items: center;
  gap: 8px;
  cursor: pointer;
}

.ccmp-cg-check-label input {
  width: 17px;
  height: 17px;
  accent-color: var(--ccmp-accent);
}

.ccmp-cg-advanced {
  border: 1px solid var(--ccmp-line);
  border-radius: 18px;
  background: color-mix(in srgb, var(--ccmp-card), transparent 10%);
}

.ccmp-cg-advanced summary {
  padding: 16px 20px;
  cursor: pointer;
  color: var(--ccmp-ink);
  font-weight: 800;
}

.ccmp-cg-advanced[open] summary {
  border-bottom: 1px solid var(--ccmp-line);
}

.ccmp-cg-advanced-body {
  display: grid;
  gap: 18px;
  padding: 20px;
}

.ccmp-cg-actions {
  display: flex;
  flex-wrap: wrap;
  gap: 10px;
  align-items: center;
}

.ccmp-cg-submit,
.ccmp-cg-secondary {
  min-height: 46px;
  padding: 10px 18px;
  border: 1px solid var(--ccmp-line);
  border-radius: 999px;
  cursor: pointer;
  font-weight: 820;
}

.ccmp-cg-submit {
  border-color: var(--ccmp-ink);
  color: var(--ccmp-paper);
  background: var(--ccmp-ink);
}

.ccmp-cg-submit:hover,
.ccmp-cg-submit:focus-visible {
  opacity: 0.88;
  transform: translateY(-1px);
}

.ccmp-cg-submit[disabled] {
  opacity: 0.58;
  cursor: wait;
  transform: none;
}

.ccmp-cg-secondary {
  color: var(--ccmp-muted);
  background: transparent;
}

.ccmp-cg-secondary:hover,
.ccmp-cg-secondary:focus-visible {
  color: var(--ccmp-ink);
  border-color: var(--ccmp-accent);
}

.ccmp-cg-api-note {
  margin-left: auto;
  color: var(--ccmp-muted);
  font-size: 0.7rem;
}

.ccmp-cg-status {
  margin-top: 20px;
  padding: 13px 15px;
  border: 1px solid var(--ccmp-line);
  border-radius: 13px;
  color: var(--ccmp-muted);
  background: var(--ccmp-paper);
  font-size: 0.8rem;
}

.ccmp-cg-status[data-kind="running"] {
  color: var(--ccmp-green);
  border-color: color-mix(in srgb, var(--ccmp-green), transparent 55%);
}

.ccmp-cg-status[data-kind="success"] {
  color: var(--ccmp-green);
  background: color-mix(in srgb, var(--ccmp-green), transparent 93%);
}

.ccmp-cg-status[data-kind="error"] {
  color: var(--ccmp-accent-deep);
  background: color-mix(in srgb, var(--ccmp-accent), transparent 92%);
}

.ccmp-cg-result {
  display: grid;
  gap: 22px;
  margin-top: 24px;
}

.ccmp-cg-summary {
  display: grid;
  grid-template-columns: repeat(6, minmax(115px, 1fr));
  gap: 10px;
}

.ccmp-cg-metric {
  min-width: 0;
  padding: 14px;
  border-top: 2px solid var(--ccmp-accent);
  border-radius: 10px;
  background: var(--ccmp-paper);
}

.ccmp-cg-metric span {
  display: block;
  color: var(--ccmp-muted);
  font-size: 0.65rem;
  font-weight: 760;
  letter-spacing: 0.04em;
  text-transform: uppercase;
}

.ccmp-cg-metric strong {
  display: block;
  margin-top: 6px;
  overflow: hidden;
  color: var(--ccmp-ink);
  font-family: var(--ccmp-serif);
  font-size: 1.3rem;
  line-height: 1.2;
  text-overflow: ellipsis;
  white-space: nowrap;
}

.ccmp-cg-panel {
  min-width: 0;
  padding: 20px;
  border: 1px solid var(--ccmp-line);
  border-radius: 18px;
  background: var(--ccmp-card);
}

.ccmp-cg-panel-head {
  display: flex;
  flex-wrap: wrap;
  justify-content: space-between;
  gap: 12px;
  align-items: center;
  margin-bottom: 13px;
}

.ccmp-cg-panel h3 {
  margin: 0;
  font-family: var(--ccmp-serif);
  font-size: 1.35rem;
  font-weight: 600;
}

.ccmp-cg-table-wrap {
  max-height: 560px;
  overflow: auto;
  border: 1px solid var(--ccmp-line);
  border-radius: 12px;
}

.ccmp-cg-table {
  width: 100%;
  border-collapse: collapse;
  font-size: 0.75rem;
  font-variant-numeric: tabular-nums;
}

.ccmp-cg-table th,
.ccmp-cg-table td {
  padding: 9px 10px;
  border-bottom: 1px solid var(--ccmp-line);
  text-align: left;
  white-space: nowrap;
}

.ccmp-cg-table th {
  position: sticky;
  top: 0;
  z-index: 1;
  color: var(--ccmp-paper);
  background: var(--ccmp-ink);
  font-size: 0.67rem;
  letter-spacing: 0.03em;
}

.ccmp-cg-table tbody tr:nth-child(even) {
  background: color-mix(in srgb, var(--ccmp-paper-deep), transparent 35%);
}

.ccmp-cg-coefficient {
  font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
  font-weight: 780;
}

.ccmp-cg-exact {
  margin-left: 7px;
  color: var(--ccmp-accent);
  font-family: var(--ccmp-serif);
}

.ccmp-cg-empty {
  padding: 34px 12px;
  color: var(--ccmp-muted);
  text-align: center;
}

.ccmp-cg-json-download {
  display: none;
}

@media (max-width: 1050px) {
  .ccmp-cg-fields {
    grid-template-columns: repeat(2, minmax(150px, 1fr));
  }

  .ccmp-cg-summary {
    grid-template-columns: repeat(3, minmax(120px, 1fr));
  }
}

@media (max-width: 700px) {
  .ccmp-cg-head {
    grid-template-columns: 1fr;
  }

  .ccmp-cg-equation {
    min-width: 0;
  }

  .ccmp-cg-fields,
  .ccmp-cg-fields-3 {
    grid-template-columns: 1fr 1fr;
  }

  .ccmp-cg-summary {
    grid-template-columns: 1fr 1fr;
  }

  .ccmp-cg-api-note {
    width: 100%;
    margin-left: 0;
  }
}

@media (max-width: 440px) {
  .ccmp-cg-fields,
  .ccmp-cg-fields-3,
  .ccmp-cg-summary {
    grid-template-columns: 1fr;
  }
}

</style>


<div class="ccmp-page" id="ccmp-top">
  <a class="ccmp-skip-link" href="#ccmp-reports">Skip to reports</a>

  <header class="ccmp-hero" aria-labelledby="ccmp-page-title">
    <div class="ccmp-hero-inner">
      <div class="ccmp-masthead">
        <div class="ccmp-brand"><span class="ccmp-brand-mark">C</span><span>Liyang · June 2026</span></div>
        <button class="ccmp-theme-toggle" id="ccmp-theme-toggle" type="button" aria-label="Switch color theme"><span>Theme</span><b aria-hidden="true">◐</b></button>
      </div>

      <div class="ccmp-hero-copy">
        <div>
          <p class="ccmp-kicker">Conference field notes · Volume I</p>
          <h1 class="ccmp-title" id="ccmp-page-title">CCMP <span>Reports</span></h1>
        </div>
        <div>
          <p class="ccmp-hero-summary">Ten editorial lecture notes tracing current ideas in high-temperature superconductivity, topology, quantum geometry, magnetic symmetry, fractonic matter and boundary conformal field theory.</p>
          <dl class="ccmp-hero-meta">
            <div><dt>Editor</dt><dd>Maggie</dd></div>
            <div><dt>Conference</dt><dd>CCMP 2026</dd></div>
            <div><dt>Location</dt><dd>Liyang, China</dd></div>
            <div><dt>Published</dt><dd>27 June 2026</dd></div>
          </dl>
        </div>
      </div>
    </div>
  </header>

  <main class="ccmp-main" id="ccmp-reports">
    <section class="ccmp-intro" aria-label="Collection overview">
      <div class="ccmp-editorial">
        <div class="ccmp-section-label">Editor’s note</div>
        <h2>A research notebook designed as a small proceedings volume.</h2>
        <p>This collection records selected talks from CCMP 2026. Each report pairs a concise scientific abstract with a direct link to the complete PDF note. PDF files are loaded only when requested, keeping the page fast and reliable even on mobile connections.</p>
        <div class="ccmp-stats" aria-label="Collection statistics">
          <div class="ccmp-stat"><strong>10</strong><span>Reports</span></div>
          <div class="ccmp-stat"><strong>10</strong><span>Speakers</span></div>
          <div class="ccmp-stat"><strong>1</strong><span>Conference</span></div>
        </div>
      </div>

      <aside class="ccmp-index" aria-label="Table of contents">
        <div class="ccmp-section-label">Index</div>
        <h2>Contents</h2>
        <nav class="ccmp-index-list">
          <a href="#so5-calculator"><b>CG</b><span>Interactive SO(5) Reduced Clebsch–Gordon Calculator</span></a>
          <a href="#report-01"><b>01</b><span><i>σ</i>-bonding High-<i>T</i><sub>c</sub> Superconductivity and Neural Tensor Network State</span></a>
          <a href="#report-02"><b>02</b><span>Phase Stiffness, Soft Modes, and Competing Orders in High-<i>T</i><sub>c</sub> Superconductors</span></a>
          <a href="#report-03"><b>03</b><span>Electric Current-Induced Superconductivity and Perfect Superconducting Diode</span></a>
          <a href="#report-04"><b>04</b><span>Stable Topology in Exactly Flat Bands (CTFB)</span></a>
          <a href="#report-05"><b>05</b><span>Spin Space Group, Altermagnetism, and Symmetry Classification of Magnetic Orders</span></a>
          <a href="#report-06"><b>06</b><span>Driven Fermions, Baths, Floquet Steady States, and Dissipation-Shaped Quantum Geometry</span></a>
          <a href="#report-07"><b>07</b><span>Fractonic Fractional Quantum Hall Effect – from Kane wires to lineons</span></a>
          <a href="#report-08"><b>08</b><span>Chirality Trinity – Chirality, Time Chirality, and Chirality Prime</span></a>
          <a href="#report-09"><b>09</b><span>Kagome CDW, Loop Current, and Superconductivity in CsV<sub>3</sub>Sb<sub>5</sub></span></a>
          <a href="#report-10"><b>10</b><span>Boundary Conformal Field Theory – Cardy condition and anomaly</span></a>
        </nav>
      </aside>
    </section>


    <section class="ccmp-cg-lab" id="so5-calculator" data-api-url="{{ site.so5_cg_api_url }}" aria-labelledby="so5-cg-title">
      <div class="ccmp-cg-head">
        <div>
          <div class="ccmp-section-label">Interactive research tool</div>
          <h2 id="so5-cg-title">SO(5) Reduced Clebsch–Gordon Calculator</h2>
          <p>
            输入表示与 SO(4) 通道量子数，服务器将调用 Julia 的高性能
            <span style="white-space:nowrap;">C<sub>2</sub> + C<sub>4</sub></span>
            计算核心并返回约化 CG 系数。首次请求包含 JIT 与矩阵缓存构造，随后相同磁扇区会更快。
          </p>
        </div>
        <div class="ccmp-cg-equation" aria-label="Coupling formula">
          \[(p,0)\otimes(p,0)\longrightarrow(P,Q)\]
          \[SO(5)\supset SO(4)\simeq SU(2)_L\times SU(2)_R\]
        </div>
      </div>

      <form class="ccmp-cg-form" id="so5-cg-form" novalidate>
        <fieldset class="ccmp-cg-fieldset">
          <legend>Representations</legend>
          <div class="ccmp-cg-fields">
            <div class="ccmp-cg-field">
              <label for="so5-p">Single-particle p</label>
              <input id="so5-p" name="p" type="number" min="0" max="10" step="1" value="10" required>
              <span class="ccmp-cg-help">当前公开接口默认限制 p ≤ 10。</span>
            </div>
            <div class="ccmp-cg-field">
              <label for="so5-q">Single-particle q</label>
              <input id="so5-q" name="q" type="number" value="0" readonly aria-readonly="true">
              <span class="ccmp-cg-help">核心仅支持 (p,0)。</span>
            </div>
            <div class="ccmp-cg-field">
              <label for="so5-P">Target P</label>
              <input id="so5-P" name="P" type="number" min="0" step="1" value="10" required>
              <span class="ccmp-cg-help">要求 0 ≤ Q ≤ P ≤ 2p。</span>
            </div>
            <div class="ccmp-cg-field">
              <label for="so5-Q">Target Q</label>
              <input id="so5-Q" name="Q" type="number" min="0" step="1" value="10" required>
              <span class="ccmp-cg-help">目标 SO(5) irrep (P,Q)。</span>
            </div>
          </div>
        </fieldset>

        <fieldset class="ccmp-cg-fieldset">
          <legend>SO(4) channel</legend>
          <div class="ccmp-cg-fields">
            <div class="ccmp-cg-field">
              <label for="so5-J">J</label>
              <input id="so5-J" name="J" type="number" min="0" step="0.5" value="0" required>
              <span class="ccmp-cg-help">整数或半整数。</span>
            </div>
            <div class="ccmp-cg-field">
              <label for="so5-M1">M₁</label>
              <input id="so5-M1" name="M1" type="number" step="0.5" value="0" required>
              <span class="ccmp-cg-help">必须满足 |M₁| ≤ J。</span>
            </div>
            <div class="ccmp-cg-field">
              <label for="so5-K">K</label>
              <input id="so5-K" name="K" type="number" min="0" step="0.5" value="0" required>
              <span class="ccmp-cg-help">整数或半整数。</span>
            </div>
            <div class="ccmp-cg-field">
              <label for="so5-M2">M₂</label>
              <input id="so5-M2" name="M2" type="number" step="0.5" value="0" required>
              <span class="ccmp-cg-help">必须满足 |M₂| ≤ K。</span>
            </div>
          </div>
          <div class="ccmp-cg-options">
            <label class="ccmp-cg-check-label">
              <input id="so5-highest-weight" type="checkbox">
              使用最高权分量 M₁ = J、M₂ = K（单目标通常最快）
            </label>
          </div>
        </fieldset>

        <details class="ccmp-cg-advanced">
          <summary>Advanced numerical parameters</summary>
          <div class="ccmp-cg-advanced-body">
            <div class="ccmp-cg-fields ccmp-cg-fields-3">
              <div class="ccmp-cg-field">
                <label for="so5-epsC4">εC4</label>
                <input id="so5-epsC4" name="epsC4" type="number" min="0" step="any" value="1e-5">
              </div>
              <div class="ccmp-cg-field">
                <label for="so5-epsJ">εJ</label>
                <input id="so5-epsJ" name="epsJ" type="number" min="0" step="any" value="1e-4">
              </div>
              <div class="ccmp-cg-field">
                <label for="so5-epsK">εK</label>
                <input id="so5-epsK" name="epsK" type="number" min="0" step="any" value="1e-6">
              </div>
            </div>

            <div class="ccmp-cg-fields">
              <div class="ccmp-cg-field">
                <label for="so5-tolC">tol C₂</label>
                <input id="so5-tolC" name="tolC" type="number" min="0" step="any" value="1e-6">
              </div>
              <div class="ccmp-cg-field">
                <label for="so5-tolC4">tol C₄</label>
                <input id="so5-tolC4" name="tolC4" type="number" min="0" step="any" value="1e-6">
              </div>
              <div class="ccmp-cg-field">
                <label for="so5-tolJ">tol J</label>
                <input id="so5-tolJ" name="tolJ" type="number" min="0" step="any" value="1e-6">
              </div>
              <div class="ccmp-cg-field">
                <label for="so5-tolK">tol K</label>
                <input id="so5-tolK" name="tolK" type="number" min="0" step="any" value="1e-6">
              </div>
            </div>

            <div class="ccmp-cg-fields ccmp-cg-fields-3">
              <div class="ccmp-cg-field">
                <label for="so5-cg-tol">SO(4) CG cutoff</label>
                <input id="so5-cg-tol" name="cg_tol" type="number" min="0" step="any" value="1e-12">
              </div>
              <div class="ccmp-cg-field">
                <label for="so5-max-rows">Maximum returned rows</label>
                <input id="so5-max-rows" name="max_rows" type="number" min="1" max="2000" step="1" value="200">
              </div>
              <div class="ccmp-cg-field">
                <label for="so5-digits">Displayed significant digits</label>
                <input id="so5-digits" name="digits" type="number" min="4" max="16" step="1" value="10">
              </div>
            </div>

            <div class="ccmp-cg-options">
              <label class="ccmp-cg-check-label"><input id="so5-calibrate-c4" type="checkbox" checked> Calibrate C₄</label>
              <label class="ccmp-cg-check-label"><input id="so5-show-diagnostics" type="checkbox" checked> Return diagnostics</label>
              <label class="ccmp-cg-check-label"><input id="so5-recognize-radicals" type="checkbox" checked> Guess simple radicals</label>
            </div>
          </div>
        </details>

        <div class="ccmp-cg-actions">
          <button class="ccmp-cg-submit" id="so5-cg-submit" type="submit">Compute CG coefficients</button>
          <button class="ccmp-cg-secondary" id="so5-cg-health" type="button">Check API</button>
          <button class="ccmp-cg-secondary" id="so5-cg-reset" type="reset">Reset</button>
          <span class="ccmp-cg-api-note" id="so5-cg-api-note">Julia backend required</span>
        </div>
      </form>

      <div class="ccmp-cg-status" id="so5-cg-status" data-kind="idle" role="status" aria-live="polite">
        Configure <code>so5_cg_api_url</code> in your Jekyll <code>_config.yml</code>, then submit a calculation.
      </div>

      <div class="ccmp-cg-result" id="so5-cg-result" hidden>
        <div class="ccmp-cg-summary" id="so5-cg-summary"></div>

        <section class="ccmp-cg-panel" aria-labelledby="so5-coeff-title">
          <div class="ccmp-cg-panel-head">
            <div>
              <h3 id="so5-coeff-title">Reduced coefficients</h3>
              <div class="ccmp-cg-help" id="so5-coeff-caption"></div>
            </div>
            <div class="ccmp-cg-actions">
              <button class="ccmp-cg-secondary" id="so5-copy-csv" type="button">Copy CSV</button>
              <button class="ccmp-cg-secondary" id="so5-download-json" type="button">Download JSON</button>
            </div>
          </div>
          <div class="ccmp-cg-table-wrap" id="so5-coeff-table-wrap"></div>
        </section>

        <section class="ccmp-cg-panel" id="so5-diagnostics-panel" aria-labelledby="so5-diagnostics-title" hidden>
          <div class="ccmp-cg-panel-head">
            <h3 id="so5-diagnostics-title">Selected-state diagnostics</h3>
          </div>
          <div class="ccmp-cg-table-wrap" id="so5-diagnostics-table-wrap"></div>
        </section>
      </div>
    </section>


    <section class="ccmp-tools" aria-label="Search and filter reports">
      <label class="ccmp-search">
        <span class="sr-only" style="position:absolute;width:1px;height:1px;padding:0;margin:-1px;overflow:hidden;clip:rect(0,0,0,0);white-space:nowrap;border:0;">Search reports</span>
        <input id="ccmp-search-input" type="search" placeholder="Search title, speaker, topic or abstract…" autocomplete="off">
      </label>
      <div class="ccmp-filters" id="ccmp-filters">
        <button class="ccmp-filter is-active" type="button" data-filter="all">All</button>
        <button class="ccmp-filter" type="button" data-filter="Superconductivity">Superconductivity</button>
        <button class="ccmp-filter" type="button" data-filter="Topology">Topology</button>
        <button class="ccmp-filter" type="button" data-filter="Symmetry">Symmetry</button>
        <button class="ccmp-filter" type="button" data-filter="Floquet">Floquet</button>
        <button class="ccmp-filter" type="button" data-filter="BCFT">BCFT</button>
      </div>
      <div class="ccmp-result-count" id="ccmp-result-count" aria-live="polite">Showing 10 reports</div>
    </section>

    <section class="ccmp-grid" aria-label="CCMP lecture reports">

      <article class="ccmp-card" id="report-01" data-tags="Superconductivity|Tensor Networks|Quantum Spin Liquid">
        <div class="ccmp-cover" aria-hidden="true">
          <div class="ccmp-cover-top">
            <span class="ccmp-report-number">01</span>
            <span class="ccmp-report-series">Condensed Matter Physics · CCMP 2026</span>
          </div>
          <div class="ccmp-cover-bottom">
            <span class="ccmp-primary-topic">Superconductivity</span>
            <span class="ccmp-cover-symbol">σ</span>
          </div>
        </div>
        <div class="ccmp-card-body">
          <h2><i>σ</i>-bonding High-<i>T</i><sub>c</sub> Superconductivity and Neural Tensor Network State</h2>
          <p class="ccmp-byline">Lecture by Tao Xiang</p>
          <div class="ccmp-tags" aria-label="Topics"><button class="ccmp-tag" type="button" data-filter="Superconductivity">Superconductivity</button><button class="ccmp-tag" type="button" data-filter="Tensor Networks">Tensor Networks</button><button class="ccmp-tag" type="button" data-filter="Quantum Spin Liquid">Quantum Spin Liquid</button></div>
          <p class="ccmp-abstract">This talk connects two frontiers: material design via metallizing \( \sigma \)-bonding electrons (exemplified by MgB\(_2\)), and variational many-body solvers using Neural Tensor Network States (\( \nu \)TNS). Xiang emphasizes that hole doping in cuprates creates Zhang-Rice singlets rather than simple rigid-band carriers. He presents recent \( \nu \)TNS benchmarks on the frustrated \( J_1 \)-\( J_2 \) square-lattice Heisenberg model at \( J_2/J_1 \sim 0.5 \), where the numerical data favor a gapless quantum spin liquid over competing VBS or Néel orders.</p>
          <div class="ccmp-actions">
            <button class="ccmp-button ccmp-button-primary ccmp-preview" type="button" data-pdf="https://maggiexheuw.github.io/Report/Classical_paper_reading%20(11).pdf" data-title="σ-bonding High-Tc Superconductivity and Neural Tensor Network State">Preview PDF</button>
            <a class="ccmp-button" href="https://maggiexheuw.github.io/Report/Classical_paper_reading%20(11).pdf" target="_blank" rel="noopener noreferrer">Open full report ↗</a>
            <a class="ccmp-button" href="#report-01" aria-label="Permalink to report 01">Permalink</a>
          </div>
        </div>
      </article>

      <article class="ccmp-card" id="report-02" data-tags="Superconductivity|Phase Fluctuations|Cuprates">
        <div class="ccmp-cover" aria-hidden="true">
          <div class="ccmp-cover-top">
            <span class="ccmp-report-number">02</span>
            <span class="ccmp-report-series">Condensed Matter Physics · CCMP 2026</span>
          </div>
          <div class="ccmp-cover-bottom">
            <span class="ccmp-primary-topic">Superconductivity</span>
            <span class="ccmp-cover-symbol">ρₛ</span>
          </div>
        </div>
        <div class="ccmp-card-body">
          <h2>Phase Stiffness, Soft Modes, and Competing Orders in High-<i>T</i><sub>c</sub> Superconductors</h2>
          <p class="ccmp-byline">Lecture by Yasutomo Uemura</p>
          <div class="ccmp-tags" aria-label="Topics"><button class="ccmp-tag" type="button" data-filter="Superconductivity">Superconductivity</button><button class="ccmp-tag" type="button" data-filter="Phase Fluctuations">Phase Fluctuations</button><button class="ccmp-tag" type="button" data-filter="Cuprates">Cuprates</button></div>
          <p class="ccmp-abstract">Uemura argues that \( T_{\rm c} \) in underdoped cuprates is controlled primarily by superfluid density (phase stiffness) rather than by the pairing gap. Using \( \mu \)SR measurements of \( \lambda^{-2} \propto n_s/m^* \), he presents the empirical Uemura relation \( T_{\rm c} \propto n_s/m^* \). The talk covers the Emery-Kivelson phase-fluctuation picture, Nernst-effect evidence for vortex-like excitations above \( T_{\rm c} \), and the role of stripe/CDW competing orders in suppressing global phase coherence.</p>
          <div class="ccmp-actions">
            <button class="ccmp-button ccmp-button-primary ccmp-preview" type="button" data-pdf="https://maggiexheuw.github.io/Report/Classical_paper_reading%20(10).pdf" data-title="Phase Stiffness, Soft Modes, and Competing Orders in High-Tc Superconductors">Preview PDF</button>
            <a class="ccmp-button" href="https://maggiexheuw.github.io/Report/Classical_paper_reading%20(10).pdf" target="_blank" rel="noopener noreferrer">Open full report ↗</a>
            <a class="ccmp-button" href="#report-02" aria-label="Permalink to report 02">Permalink</a>
          </div>
        </div>
      </article>

      <article class="ccmp-card" id="report-03" data-tags="Superconductivity|Nonreciprocity|FFLO">
        <div class="ccmp-cover" aria-hidden="true">
          <div class="ccmp-cover-top">
            <span class="ccmp-report-number">03</span>
            <span class="ccmp-report-series">Condensed Matter Physics · CCMP 2026</span>
          </div>
          <div class="ccmp-cover-bottom">
            <span class="ccmp-primary-topic">Superconductivity</span>
            <span class="ccmp-cover-symbol">Jc</span>
          </div>
        </div>
        <div class="ccmp-card-body">
          <h2>Electric Current-Induced Superconductivity and Perfect Superconducting Diode</h2>
          <p class="ccmp-byline">Lecture by Youichi Yanase</p>
          <div class="ccmp-tags" aria-label="Topics"><button class="ccmp-tag" type="button" data-filter="Superconductivity">Superconductivity</button><button class="ccmp-tag" type="button" data-filter="Nonreciprocity">Nonreciprocity</button><button class="ccmp-tag" type="button" data-filter="FFLO">FFLO</button></div>
          <p class="ccmp-abstract">Yanase systematically introduces the superconducting diode effect (SDE) arising from finite-momentum Cooper pairs, including Rashba helical states and FFLO-type pairing. He highlights the Daido-Yanase bilayer dissipation model that achieves perfect SDE (vanishing critical current in one direction), and the trigonal Fulde-Ferrell scenario where a finite current selects a single FF domain to produce a zero-resistance state. Light-driven and correlation-induced nonreciprocal mechanisms are also discussed.</p>
          <div class="ccmp-actions">
            <button class="ccmp-button ccmp-button-primary ccmp-preview" type="button" data-pdf="https://maggiexheuw.github.io/Report/Classical_paper_reading%20(9).pdf" data-title="Electric Current-Induced Superconductivity and Perfect Superconducting Diode">Preview PDF</button>
            <a class="ccmp-button" href="https://maggiexheuw.github.io/Report/Classical_paper_reading%20(9).pdf" target="_blank" rel="noopener noreferrer">Open full report ↗</a>
            <a class="ccmp-button" href="#report-03" aria-label="Permalink to report 03">Permalink</a>
          </div>
        </div>
      </article>

      <article class="ccmp-card" id="report-04" data-tags="Topology|Flat Bands|Chern Insulators">
        <div class="ccmp-cover" aria-hidden="true">
          <div class="ccmp-cover-top">
            <span class="ccmp-report-number">04</span>
            <span class="ccmp-report-series">Condensed Matter Physics · CCMP 2026</span>
          </div>
          <div class="ccmp-cover-bottom">
            <span class="ccmp-primary-topic">Topology</span>
            <span class="ccmp-cover-symbol">C</span>
          </div>
        </div>
        <div class="ccmp-card-body">
          <h2>Stable Topology in Exactly Flat Bands (CTFB)</h2>
          <p class="ccmp-byline">Lecture by Zhida Song</p>
          <div class="ccmp-tags" aria-label="Topics"><button class="ccmp-tag" type="button" data-filter="Topology">Topology</button><button class="ccmp-tag" type="button" data-filter="Flat Bands">Flat Bands</button><button class="ccmp-tag" type="button" data-filter="Chern Insulators">Chern Insulators</button></div>
          <p class="ccmp-abstract">Song addresses the no-go theorem that forbids strictly flat, finite-range, gapped bands with stable topology. By relaxing the gapped condition and allowing critical band touching, he constructs Critical Topological Flat Bands (CTFB) where the projector remains continuous (though non-analytic) at the touching point. Using bipartite kernel constructions and symmetry indicators (\( \Delta B \)), he shows how stable Chern/\( Z_2 \) topology can coexist with exact flatness, providing a rigorous starting point for fractional Chern insulators.</p>
          <div class="ccmp-actions">
            <button class="ccmp-button ccmp-button-primary ccmp-preview" type="button" data-pdf="https://maggiexheuw.github.io/Report/Classical_paper_reading%20(8).pdf" data-title="Stable Topology in Exactly Flat Bands (CTFB)">Preview PDF</button>
            <a class="ccmp-button" href="https://maggiexheuw.github.io/Report/Classical_paper_reading%20(8).pdf" target="_blank" rel="noopener noreferrer">Open full report ↗</a>
            <a class="ccmp-button" href="#report-04" aria-label="Permalink to report 04">Permalink</a>
          </div>
        </div>
      </article>

      <article class="ccmp-card" id="report-05" data-tags="Altermagnetism|Symmetry|Magnetism">
        <div class="ccmp-cover" aria-hidden="true">
          <div class="ccmp-cover-top">
            <span class="ccmp-report-number">05</span>
            <span class="ccmp-report-series">Condensed Matter Physics · CCMP 2026</span>
          </div>
          <div class="ccmp-cover-bottom">
            <span class="ccmp-primary-topic">Altermagnetism</span>
            <span class="ccmp-cover-symbol">𝒢</span>
          </div>
        </div>
        <div class="ccmp-card-body">
          <h2>Spin Space Group, Altermagnetism, and Symmetry Classification of Magnetic Orders</h2>
          <p class="ccmp-byline">Lecture by Qihang Liu</p>
          <div class="ccmp-tags" aria-label="Topics"><button class="ccmp-tag" type="button" data-filter="Altermagnetism">Altermagnetism</button><button class="ccmp-tag" type="button" data-filter="Symmetry">Symmetry</button><button class="ccmp-tag" type="button" data-filter="Magnetism">Magnetism</button></div>
          <p class="ccmp-abstract">Liu introduces spin space groups to redefine ferromagnets and antiferromagnets based on spin-space operations rather than net magnetization. He focuses on altermagnetism — a collinear compensated phase with zero net moment but momentum-dependent spin splitting, arising from nonrelativistic crystal symmetry. The talk also covers spin translational groups for classifying complex AFM geometries (spiral, multi-\( q \)), and how oriented spin space groups diagnose allowed responses such as AHE, magnetoelectricity, and topology.</p>
          <div class="ccmp-actions">
            <button class="ccmp-button ccmp-button-primary ccmp-preview" type="button" data-pdf="https://maggiexheuw.github.io/Report/Classical_paper_reading%20(7).pdf" data-title="Spin Space Group, Altermagnetism, and Symmetry Classification of Magnetic Orders">Preview PDF</button>
            <a class="ccmp-button" href="https://maggiexheuw.github.io/Report/Classical_paper_reading%20(7).pdf" target="_blank" rel="noopener noreferrer">Open full report ↗</a>
            <a class="ccmp-button" href="#report-05" aria-label="Permalink to report 05">Permalink</a>
          </div>
        </div>
      </article>

      <article class="ccmp-card" id="report-06" data-tags="Floquet|Open Quantum Systems|Quantum Geometry">
        <div class="ccmp-cover" aria-hidden="true">
          <div class="ccmp-cover-top">
            <span class="ccmp-report-number">06</span>
            <span class="ccmp-report-series">Condensed Matter Physics · CCMP 2026</span>
          </div>
          <div class="ccmp-cover-bottom">
            <span class="ccmp-primary-topic">Floquet</span>
            <span class="ccmp-cover-symbol">Γ</span>
          </div>
        </div>
        <div class="ccmp-card-body">
          <h2>Driven Fermions, Baths, Floquet Steady States, and Dissipation-Shaped Quantum Geometry</h2>
          <p class="ccmp-byline">Lecture by Likun Shi</p>
          <div class="ccmp-tags" aria-label="Topics"><button class="ccmp-tag" type="button" data-filter="Floquet">Floquet</button><button class="ccmp-tag" type="button" data-filter="Open Quantum Systems">Open Quantum Systems</button><button class="ccmp-tag" type="button" data-filter="Quantum Geometry">Quantum Geometry</button></div>
          <p class="ccmp-abstract">Shi shows that Floquet band structure alone does not determine occupation; the nature of the bath (fermionic vs. bosonic) selects distinct non-equilibrium steady states. Fermionic baths yield Floquet Fermi liquids with nested Fermi surfaces, while bosonic baths produce ultracritical Floquet non-Fermi liquids with sharp but non-Fermi singularities. In the second half, he warns that the \( \Gamma^0 \) term in second-order nonlinear transport — often assumed to be purely quantum-geometric — is in fact shaped by dissipation, making bath identification essential for extracting quantum metrics.</p>
          <div class="ccmp-actions">
            <button class="ccmp-button ccmp-button-primary ccmp-preview" type="button" data-pdf="https://maggiexheuw.github.io/Report/Classical_paper_reading%20(6).pdf" data-title="Driven Fermions, Baths, Floquet Steady States, and Dissipation-Shaped Quantum Geometry">Preview PDF</button>
            <a class="ccmp-button" href="https://maggiexheuw.github.io/Report/Classical_paper_reading%20(6).pdf" target="_blank" rel="noopener noreferrer">Open full report ↗</a>
            <a class="ccmp-button" href="#report-06" aria-label="Permalink to report 06">Permalink</a>
          </div>
        </div>
      </article>

      <article class="ccmp-card" id="report-07" data-tags="Fractons|FQH|Topological Order">
        <div class="ccmp-cover" aria-hidden="true">
          <div class="ccmp-cover-top">
            <span class="ccmp-report-number">07</span>
            <span class="ccmp-report-series">Condensed Matter Physics · CCMP 2026</span>
          </div>
          <div class="ccmp-cover-bottom">
            <span class="ccmp-primary-topic">Fractons</span>
            <span class="ccmp-cover-symbol">ν</span>
          </div>
        </div>
        <div class="ccmp-card-body">
          <h2>Fractonic Fractional Quantum Hall Effect – from Kane wires to lineons</h2>
          <p class="ccmp-byline">Lecture by Adolfo Grushin</p>
          <div class="ccmp-tags" aria-label="Topics"><button class="ccmp-tag" type="button" data-filter="Fractons">Fractons</button><button class="ccmp-tag" type="button" data-filter="FQH">FQH</button><button class="ccmp-tag" type="button" data-filter="Topological Order">Topological Order</button></div>
          <p class="ccmp-abstract">This talk generalizes Kane's coupled-wire construction by gluing adjacent Laughlin strips with different filling fractions \( 1/m_j \), using Haldane's null-vector criterion. The non-uniform gluing generates position-dependent anyon condensation rules, which split the ordinary 2D FQH anyons into mobility-restricted excitations: lineons (1D motion), spread lineons (local interface crossing), and fully mobile C-anyons. This provides a natural route to fractonic topological order without requiring crystalline translation symmetry.</p>
          <div class="ccmp-actions">
            <button class="ccmp-button ccmp-button-primary ccmp-preview" type="button" data-pdf="https://maggiexheuw.github.io/Report/Classical_paper_reading%20(5).pdf" data-title="Fractonic Fractional Quantum Hall Effect – from Kane wires to lineons">Preview PDF</button>
            <a class="ccmp-button" href="https://maggiexheuw.github.io/Report/Classical_paper_reading%20(5).pdf" target="_blank" rel="noopener noreferrer">Open full report ↗</a>
            <a class="ccmp-button" href="#report-07" aria-label="Permalink to report 07">Permalink</a>
          </div>
        </div>
      </article>

      <article class="ccmp-card" id="report-08" data-tags="Chirality|Magnetoelectricity|Symmetry">
        <div class="ccmp-cover" aria-hidden="true">
          <div class="ccmp-cover-top">
            <span class="ccmp-report-number">08</span>
            <span class="ccmp-report-series">Condensed Matter Physics · CCMP 2026</span>
          </div>
          <div class="ccmp-cover-bottom">
            <span class="ccmp-primary-topic">Chirality</span>
            <span class="ccmp-cover-symbol">PT</span>
          </div>
        </div>
        <div class="ccmp-card-body">
          <h2>Chirality Trinity – Chirality, Time Chirality, and Chirality Prime</h2>
          <p class="ccmp-byline">Lecture by Sang-Wook Cheong</p>
          <div class="ccmp-tags" aria-label="Topics"><button class="ccmp-tag" type="button" data-filter="Chirality">Chirality</button><button class="ccmp-tag" type="button" data-filter="Magnetoelectricity">Magnetoelectricity</button><button class="ccmp-tag" type="button" data-filter="Symmetry">Symmetry</button></div>
          <p class="ccmp-abstract">Cheong classifies handedness using \( P \) (parity), \( T \) (time-reversal), and \( PT \) symmetries. Ordinary chirality (\( \mathbf{H}\cdot\mathbf{k} \)) is \( P \)-odd/\( T \)-even; time chirality (\( \mathbf{k}\cdot\mathbf{E} \)) is \( P \)-even/\( T \)-odd; chirality prime (\( \mathbf{E}\cdot\mathbf{H} \)) is \( P \)-odd/\( T \)-odd. Superchirality occurs when all three are broken. He applies this framework to ferro-rotational (electric toroidal) orders, circular dichroism, chiral vs. time-chiral phonons, and shows how these selection rules dictate Hall, Edelstein, and magnetoelectric responses without detailed microscopic calculations.</p>
          <div class="ccmp-actions">
            <button class="ccmp-button ccmp-button-primary ccmp-preview" type="button" data-pdf="https://maggiexheuw.github.io/Report/Classical_paper_reading%20(4).pdf" data-title="Chirality Trinity – Chirality, Time Chirality, and Chirality Prime">Preview PDF</button>
            <a class="ccmp-button" href="https://maggiexheuw.github.io/Report/Classical_paper_reading%20(4).pdf" target="_blank" rel="noopener noreferrer">Open full report ↗</a>
            <a class="ccmp-button" href="#report-08" aria-label="Permalink to report 08">Permalink</a>
          </div>
        </div>
      </article>

      <article class="ccmp-card" id="report-09" data-tags="Kagome|CDW|Superconductivity">
        <div class="ccmp-cover" aria-hidden="true">
          <div class="ccmp-cover-top">
            <span class="ccmp-report-number">09</span>
            <span class="ccmp-report-series">Condensed Matter Physics · CCMP 2026</span>
          </div>
          <div class="ccmp-cover-bottom">
            <span class="ccmp-primary-topic">Kagome</span>
            <span class="ccmp-cover-symbol">Δ</span>
          </div>
        </div>
        <div class="ccmp-card-body">
          <h2>Kagome CDW, Loop Current, and Superconductivity in CsV<sub>3</sub>Sb<sub>5</sub></h2>
          <p class="ccmp-byline">Lecture by Xianhui Chen</p>
          <div class="ccmp-tags" aria-label="Topics"><button class="ccmp-tag" type="button" data-filter="Kagome">Kagome</button><button class="ccmp-tag" type="button" data-filter="CDW">CDW</button><button class="ccmp-tag" type="button" data-filter="Superconductivity">Superconductivity</button></div>
          <p class="ccmp-abstract">Chen presents experimental evidence that the CDW in CsV\(_3\)Sb\(_5\) is not a simple Peierls nesting instability. The transition is accompanied by a giant anomalous Hall effect, NMR-detected orbital ordering, electronic nematicity, and a double-dome superconducting phase under pressure. He interprets the CDW as a multi-component order parameter involving charge, bond, orbital, and possibly loop-current modulations, drawing on the Capponi-Wu-Zhang bilayer model where \( \operatorname{Im}\langle c_i^\dagger c_j \rangle \) serves as the order parameter for spontaneous circulating currents.</p>
          <div class="ccmp-actions">
            <button class="ccmp-button ccmp-button-primary ccmp-preview" type="button" data-pdf="https://maggiexheuw.github.io/Report/Classical_paper_reading%20(3).pdf" data-title="Kagome CDW, Loop Current, and Superconductivity in CsV3Sb5">Preview PDF</button>
            <a class="ccmp-button" href="https://maggiexheuw.github.io/Report/Classical_paper_reading%20(3).pdf" target="_blank" rel="noopener noreferrer">Open full report ↗</a>
            <a class="ccmp-button" href="#report-09" aria-label="Permalink to report 09">Permalink</a>
          </div>
        </div>
      </article>

      <article class="ccmp-card" id="report-10" data-tags="BCFT|Anomaly|Conformal Field Theory">
        <div class="ccmp-cover" aria-hidden="true">
          <div class="ccmp-cover-top">
            <span class="ccmp-report-number">10</span>
            <span class="ccmp-report-series">Condensed Matter Physics · CCMP 2026</span>
          </div>
          <div class="ccmp-cover-bottom">
            <span class="ccmp-primary-topic">BCFT</span>
            <span class="ccmp-cover-symbol">Lₙ</span>
          </div>
        </div>
        <div class="ccmp-card-body">
          <h2>Boundary Conformal Field Theory – Cardy condition and anomaly</h2>
          <p class="ccmp-byline">Lecture by Masaki Oshikawa</p>
          <div class="ccmp-tags" aria-label="Topics"><button class="ccmp-tag" type="button" data-filter="BCFT">BCFT</button><button class="ccmp-tag" type="button" data-filter="Anomaly">Anomaly</button><button class="ccmp-tag" type="button" data-filter="Conformal Field Theory">Conformal Field Theory</button></div>
          <p class="ccmp-abstract">This is a pedagogical note on BCFT, starting from the gluing condition \( T(z) = \bar{T}(\bar{z}) \) on the boundary, leading to Ishibashi states and Cardy states. The annulus partition function is computed in both open and closed channels, and their equality yields the Cardy condition, which is then solved via the Verlinde formula. The final section explains Oshikawa's criterion: existence of a \( G \)-invariant Cardy boundary state diagnoses the absence of 't Hooft anomaly, illustrated with \( \mathrm{SU}(N)_k \) WZW models where center symmetry dictates \( k \in N\mathbb{Z} \) for anomaly-free gapping.</p>
          <div class="ccmp-actions">
            <button class="ccmp-button ccmp-button-primary ccmp-preview" type="button" data-pdf="https://maggiexheuw.github.io/Report/Classical_paper_reading%20(2).pdf" data-title="Boundary Conformal Field Theory – Cardy condition and anomaly">Preview PDF</button>
            <a class="ccmp-button" href="https://maggiexheuw.github.io/Report/Classical_paper_reading%20(2).pdf" target="_blank" rel="noopener noreferrer">Open full report ↗</a>
            <a class="ccmp-button" href="#report-10" aria-label="Permalink to report 10">Permalink</a>
          </div>
        </div>
      </article>
    </section>

    <div class="ccmp-empty" id="ccmp-empty" hidden>
      <strong>No matching reports</strong>
      <span>Try a broader keyword or select “All”.</span>
    </div>
  </main>

  <footer class="ccmp-footer">
    <div class="ccmp-footer-inner">
      <div>
        <h2 class="ccmp-footer-title">CCMP Reports · 2026</h2>
        <p>Curated and written by Maggie · Last updated 27 June 2026</p>
      </div>
      <a class="ccmp-back-top" href="#ccmp-top">Back to top ↑</a>
    </div>
  </footer>

  <div class="ccmp-modal" id="ccmp-modal" hidden aria-hidden="true">
    <div class="ccmp-modal-panel" role="dialog" aria-modal="true" aria-labelledby="ccmp-modal-title">
      <div class="ccmp-modal-header">
        <h2 class="ccmp-modal-title" id="ccmp-modal-title">PDF preview</h2>
        <button class="ccmp-modal-close" id="ccmp-modal-close" type="button" aria-label="Close PDF preview">×</button>
      </div>
      <iframe id="ccmp-pdf-frame" title="PDF preview" loading="lazy"></iframe>
    </div>
  </div>
</div>

<script>
(function () {
  'use strict';

  function initCCMPPage() {
    var root = document.querySelector('.ccmp-page');
    if (!root) return;

    var cards = Array.prototype.slice.call(root.querySelectorAll('.ccmp-card'));
    var searchInput = root.querySelector('#ccmp-search-input');
    var filters = Array.prototype.slice.call(root.querySelectorAll('[data-filter]'));
    var filterButtons = Array.prototype.slice.call(root.querySelectorAll('.ccmp-filter'));
    var resultCount = root.querySelector('#ccmp-result-count');
    var emptyState = root.querySelector('#ccmp-empty');
    var activeFilter = 'all';

    cards.forEach(function (card) {
      card.setAttribute('data-search-text', card.textContent.toLowerCase().replace(/\s+/g, ' ').trim());
    });

    function applyFilters() {
      var query = searchInput.value.toLowerCase().trim();
      var visible = 0;

      cards.forEach(function (card) {
        var textMatch = !query || card.getAttribute('data-search-text').indexOf(query) !== -1;
        var tags = (card.getAttribute('data-tags') || '').split('|');
        var tagMatch = activeFilter === 'all' || tags.indexOf(activeFilter) !== -1;
        var show = textMatch && tagMatch;
        card.classList.toggle('is-hidden', !show);
        if (show) visible += 1;
      });

      resultCount.textContent = 'Showing ' + visible + ' report' + (visible === 1 ? '' : 's');
      emptyState.hidden = visible !== 0;
    }

    searchInput.addEventListener('input', applyFilters);

    filters.forEach(function (button) {
      button.addEventListener('click', function () {
        activeFilter = button.getAttribute('data-filter') || 'all';
        filterButtons.forEach(function (item) {
          item.classList.toggle('is-active', item.getAttribute('data-filter') === activeFilter);
        });
        applyFilters();
        if (button.classList.contains('ccmp-tag')) {
          root.querySelector('.ccmp-tools').scrollIntoView({ behavior: 'smooth', block: 'start' });
        }
      });
    });

    var themeButton = root.querySelector('#ccmp-theme-toggle');
    var savedTheme = null;
    try { savedTheme = localStorage.getItem('ccmp-theme'); } catch (error) { savedTheme = null; }
    if (savedTheme === 'dark' || savedTheme === 'light') {
      root.setAttribute('data-ccmp-theme', savedTheme);
    } else if (window.matchMedia && window.matchMedia('(prefers-color-scheme: dark)').matches) {
      root.setAttribute('data-ccmp-theme', 'dark');
    }

    themeButton.addEventListener('click', function () {
      var current = root.getAttribute('data-ccmp-theme');
      var darkBySystem = window.matchMedia && window.matchMedia('(prefers-color-scheme: dark)').matches;
      var next = current === 'dark' || (!current && darkBySystem) ? 'light' : 'dark';
      root.setAttribute('data-ccmp-theme', next);
      try { localStorage.setItem('ccmp-theme', next); } catch (error) { /* storage may be disabled */ }
    });

    var modal = root.querySelector('#ccmp-modal');
    var modalTitle = root.querySelector('#ccmp-modal-title');
    var modalClose = root.querySelector('#ccmp-modal-close');
    var pdfFrame = root.querySelector('#ccmp-pdf-frame');
    var lastFocused = null;

    function closeModal() {
      modal.hidden = true;
      modal.setAttribute('aria-hidden', 'true');
      pdfFrame.removeAttribute('src');
      document.documentElement.style.overflow = '';
      if (lastFocused) lastFocused.focus();
    }

    Array.prototype.slice.call(root.querySelectorAll('.ccmp-preview')).forEach(function (button) {
      button.addEventListener('click', function () {
        lastFocused = button;
        modalTitle.textContent = button.getAttribute('data-title') || 'PDF preview';
        pdfFrame.setAttribute('src', button.getAttribute('data-pdf') + '#page=1&view=FitH&toolbar=1&navpanes=0');
        modal.hidden = false;
        modal.setAttribute('aria-hidden', 'false');
        document.documentElement.style.overflow = 'hidden';
        modalClose.focus();
      });
    });

    modalClose.addEventListener('click', closeModal);
    modal.addEventListener('click', function (event) {
      if (event.target === modal) closeModal();
    });
    document.addEventListener('keydown', function (event) {
      if (event.key === 'Escape' && !modal.hidden) closeModal();
    });
  }

  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', initCCMPPage);
  } else {
    initCCMPPage();
  }
}());
</script>

<script>
(function () {
  'use strict';

  function initSO5Calculator() {
    var root = document.querySelector('#so5-calculator');
    if (!root) return;

    var form = root.querySelector('#so5-cg-form');
    var submitButton = root.querySelector('#so5-cg-submit');
    var healthButton = root.querySelector('#so5-cg-health');
    var resetButton = root.querySelector('#so5-cg-reset');
    var status = root.querySelector('#so5-cg-status');
    var resultBox = root.querySelector('#so5-cg-result');
    var summary = root.querySelector('#so5-cg-summary');
    var coeffWrap = root.querySelector('#so5-coeff-table-wrap');
    var coeffCaption = root.querySelector('#so5-coeff-caption');
    var diagnosticsPanel = root.querySelector('#so5-diagnostics-panel');
    var diagnosticsWrap = root.querySelector('#so5-diagnostics-table-wrap');
    var apiNote = root.querySelector('#so5-cg-api-note');
    var highestWeight = root.querySelector('#so5-highest-weight');
    var inputJ = root.querySelector('#so5-J');
    var inputK = root.querySelector('#so5-K');
    var inputM1 = root.querySelector('#so5-M1');
    var inputM2 = root.querySelector('#so5-M2');
    var latestPayload = null;

    var configuredURL = (window.SO5_CG_API_URL || root.getAttribute('data-api-url') || '').trim();
    var apiURL = configuredURL.replace(/\/$/, '');
    apiNote.textContent = apiURL ? 'API: ' + apiURL : 'API URL is not configured';

    function setStatus(message, kind) {
      status.textContent = message;
      status.setAttribute('data-kind', kind || 'idle');
    }

    function readNumber(id) {
      var element = root.querySelector(id);
      return Number(element.value);
    }

    function isHalfInteger(value) {
      return Number.isFinite(value) && Math.abs(value * 2 - Math.round(value * 2)) < 1e-9;
    }

    function syncHighestWeight() {
      if (!highestWeight.checked) return;
      inputM1.value = inputJ.value;
      inputM2.value = inputK.value;
    }

    highestWeight.addEventListener('change', syncHighestWeight);
    inputJ.addEventListener('input', syncHighestWeight);
    inputK.addEventListener('input', syncHighestWeight);

    function requestPayload() {
      return {
        p: readNumber('#so5-p'),
        P: readNumber('#so5-P'),
        Q: readNumber('#so5-Q'),
        J: readNumber('#so5-J'),
        M1: readNumber('#so5-M1'),
        K: readNumber('#so5-K'),
        M2: readNumber('#so5-M2'),
        epsC4: readNumber('#so5-epsC4'),
        epsJ: readNumber('#so5-epsJ'),
        epsK: readNumber('#so5-epsK'),
        tolC: readNumber('#so5-tolC'),
        tolC4: readNumber('#so5-tolC4'),
        tolJ: readNumber('#so5-tolJ'),
        tolK: readNumber('#so5-tolK'),
        cg_tol: readNumber('#so5-cg-tol'),
        max_rows: readNumber('#so5-max-rows'),
        digits: readNumber('#so5-digits'),
        calibrate_C4: root.querySelector('#so5-calibrate-c4').checked,
        show_diagnostics: root.querySelector('#so5-show-diagnostics').checked,
        recognize_radicals: root.querySelector('#so5-recognize-radicals').checked
      };
    }

    function validatePayload(data) {
      if (!Number.isInteger(data.p) || data.p < 0 || data.p > 10) return 'p must be an integer between 0 and 10.';
      if (!Number.isInteger(data.P) || !Number.isInteger(data.Q)) return 'P and Q must be integers.';
      if (data.Q < 0 || data.Q > data.P || data.P > 2 * data.p) return 'Need 0 ≤ Q ≤ P ≤ 2p.';
      if (![data.J, data.M1, data.K, data.M2].every(isHalfInteger)) return 'J, M₁, K and M₂ must be integers or half-integers.';
      if (data.J < 0 || data.K < 0 || data.J > data.p || data.K > data.p) return 'Need 0 ≤ J,K ≤ p.';
      if (Math.abs(data.M1) > data.J || Math.abs(data.M2) > data.K) return 'Need |M₁| ≤ J and |M₂| ≤ K.';
      if (![data.epsC4, data.epsJ, data.epsK, data.tolC, data.tolC4, data.tolJ, data.tolK, data.cg_tol].every(function (x) { return Number.isFinite(x) && x > 0; })) {
        return 'All epsilon and tolerance values must be positive finite numbers.';
      }
      return '';
    }

    function formatNumber(value, digits) {
      if (!Number.isFinite(Number(value))) return String(value);
      var number = Number(value);
      if (number === 0) return '0';
      return number.toPrecision(digits).replace(/(?:\.0+|(?:(\.\d*?)0+))(?=e|$)/, '$1');
    }

    function metric(label, value) {
      var box = document.createElement('div');
      box.className = 'ccmp-cg-metric';
      var name = document.createElement('span');
      name.textContent = label;
      var strong = document.createElement('strong');
      strong.textContent = String(value);
      box.appendChild(name);
      box.appendChild(strong);
      return box;
    }

    function buildTable(headers, rows) {
      var table = document.createElement('table');
      table.className = 'ccmp-cg-table';
      var thead = document.createElement('thead');
      var headRow = document.createElement('tr');
      headers.forEach(function (header) {
        var th = document.createElement('th');
        th.textContent = header;
        headRow.appendChild(th);
      });
      thead.appendChild(headRow);
      table.appendChild(thead);

      var tbody = document.createElement('tbody');
      rows.forEach(function (cells) {
        var tr = document.createElement('tr');
        cells.forEach(function (cell) {
          var td = document.createElement('td');
          if (cell && cell.nodeType) td.appendChild(cell);
          else td.textContent = cell === null || cell === undefined ? '' : String(cell);
          tr.appendChild(td);
        });
        tbody.appendChild(tr);
      });
      table.appendChild(tbody);
      return table;
    }

    function renderResponse(data) {
      latestPayload = data;
      var digits = Number(root.querySelector('#so5-digits').value) || 10;
      var dims = data.dimensions || {};
      var targets = data.targets || {};

      summary.textContent = '';
      summary.appendChild(metric('Runtime', formatNumber(data.elapsed_seconds, 6) + ' s'));
      summary.appendChild(metric('Single-particle dim', dims.single_particle));
      summary.appendChild(metric('Sector dim', dims.product_sector));
      summary.appendChild(metric('Selected vectors', dims.selected_vectors));
      summary.appendChild(metric('Target C₂', formatNumber(targets.C2, digits)));
      summary.appendChild(metric('Target C₄', formatNumber(targets.C4, digits)));

      var coefficients = Array.isArray(data.coefficients) ? data.coefficients : [];
      coeffWrap.textContent = '';
      coeffCaption.textContent = 'Showing ' + (dims.returned_coefficients || coefficients.length) + ' of ' + (dims.coefficient_count || coefficients.length) + ' coefficient rows.';

      if (!coefficients.length) {
        var empty = document.createElement('div');
        empty.className = 'ccmp-cg-empty';
        empty.textContent = 'No reduced coefficients were extracted for this channel.';
        coeffWrap.appendChild(empty);
      } else {
        var coeffRows = coefficients.map(function (row) {
          var coefficient = document.createElement('span');
          coefficient.className = 'ccmp-cg-coefficient';
          coefficient.textContent = formatNumber(row.coefficient, digits);
          if (row.exact) {
            var exact = document.createElement('span');
            exact.className = 'ccmp-cg-exact';
            exact.textContent = '(' + row.exact + ')';
            coefficient.appendChild(exact);
          }
          return [
            row.channel,
            '(' + row.P + ',' + row.Q + ')',
            row.J_label,
            row.K_label,
            '(' + row.j1_label + ',' + row.k1_label + ')',
            '(' + row.j2_label + ',' + row.k2_label + ')',
            coefficient
          ];
        });
        coeffWrap.appendChild(buildTable(['channel', 'SO(5)', 'J', 'K', '(j₁,k₁)', '(j₂,k₂)', 'coefficient'], coeffRows));
      }

      var diagnostics = Array.isArray(data.diagnostics) ? data.diagnostics : [];
      diagnosticsWrap.textContent = '';
      diagnosticsPanel.hidden = diagnostics.length === 0;
      if (diagnostics.length) {
        var diagnosticRows = diagnostics.map(function (row, index) {
          return [
            index + 1,
            formatNumber(row.O_eigenvalue, digits),
            formatNumber(row.C2_SO5, digits),
            formatNumber(row.C4_SO5, digits),
            formatNumber(row.J2, digits),
            formatNumber(row.K2, digits),
            formatNumber(row.resC, 4),
            formatNumber(row.resC4, 4),
            formatNumber(row.resJ, 4),
            formatNumber(row.resK, 4)
          ];
        });
        diagnosticsWrap.appendChild(buildTable(['state', 'O', 'C₂', 'C₄', 'J²', 'K²', 'res C₂', 'res C₄', 'res J', 'res K'], diagnosticRows));
      }

      resultBox.hidden = false;
      if (window.MathJax && window.MathJax.typesetPromise) {
        window.MathJax.typesetPromise([root]).catch(function () {});
      }
    }

    async function readJSONResponse(response) {
      var text = await response.text();
      try {
        return JSON.parse(text);
      } catch (error) {
        throw new Error('The API returned a non-JSON response (HTTP ' + response.status + ').');
      }
    }

    form.addEventListener('submit', async function (event) {
      event.preventDefault();
      if (!apiURL) {
        setStatus('API URL is not configured. Set so5_cg_api_url in _config.yml.', 'error');
        return;
      }

      syncHighestWeight();
      var payload = requestPayload();
      var validationError = validatePayload(payload);
      if (validationError) {
        setStatus(validationError, 'error');
        return;
      }

      submitButton.disabled = true;
      resultBox.hidden = true;
      setStatus('Julia is computing the fixed magnetic sector. The first request may be slower because of JIT compilation and cache construction…', 'running');

      try {
        var response = await fetch(apiURL, {
          method: 'POST',
          mode: 'cors',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify(payload)
        });
        var data = await readJSONResponse(response);
        if (!response.ok || !data.ok) {
          throw new Error(data.error || data.detail || ('HTTP ' + response.status));
        }
        renderResponse(data);
        setStatus('Calculation finished in ' + formatNumber(data.elapsed_seconds, 6) + ' s. Repeating the same magnetic sector can reuse the Julia caches.', 'success');
      } catch (error) {
        setStatus('Calculation failed: ' + error.message, 'error');
      } finally {
        submitButton.disabled = false;
      }
    });

    healthButton.addEventListener('click', async function () {
      if (!apiURL) {
        setStatus('API URL is not configured. Set so5_cg_api_url in _config.yml.', 'error');
        return;
      }
      var healthURL = apiURL.replace(/\/api\/cg$/, '/health');
      setStatus('Checking Julia API…', 'running');
      try {
        var response = await fetch(healthURL, { mode: 'cors' });
        var data = await readJSONResponse(response);
        if (!response.ok || !data.ok) throw new Error(data.error || ('HTTP ' + response.status));
        setStatus('API online · Julia ' + data.julia_version + ' · Julia threads ' + data.julia_threads + ' · BLAS threads ' + data.blas_threads + ' · max p = ' + data.max_p, 'success');
      } catch (error) {
        setStatus('API check failed: ' + error.message, 'error');
      }
    });

    resetButton.addEventListener('click', function () {
      window.setTimeout(function () {
        latestPayload = null;
        resultBox.hidden = true;
        highestWeight.checked = false;
        setStatus(apiURL ? 'Ready.' : 'API URL is not configured. Set so5_cg_api_url in _config.yml.', apiURL ? 'idle' : 'error');
      }, 0);
    });

    root.querySelector('#so5-copy-csv').addEventListener('click', async function () {
      if (!latestPayload || !Array.isArray(latestPayload.coefficients)) return;
      var lines = ['channel,P,Q,J,K,j1,k1,j2,k2,coefficient,exact'];
      latestPayload.coefficients.forEach(function (row) {
        lines.push([
          row.channel, row.P, row.Q, row.J, row.K,
          row.j1, row.k1, row.j2, row.k2,
          row.coefficient, row.exact || ''
        ].map(function (value) {
          var text = String(value).replace(/"/g, '""');
          return '"' + text + '"';
        }).join(','));
      });
      try {
        await navigator.clipboard.writeText(lines.join('\n'));
        setStatus('Coefficient table copied as CSV.', 'success');
      } catch (error) {
        setStatus('Could not copy CSV: ' + error.message, 'error');
      }
    });

    root.querySelector('#so5-download-json').addEventListener('click', function () {
      if (!latestPayload) return;
      var blob = new Blob([JSON.stringify(latestPayload, null, 2)], { type: 'application/json' });
      var url = URL.createObjectURL(blob);
      var link = document.createElement('a');
      link.href = url;
      link.download = 'so5-cg-result.json';
      document.body.appendChild(link);
      link.click();
      link.remove();
      URL.revokeObjectURL(url);
    });

    if (apiURL) setStatus('Ready. The first calculation may be slower; cached sectors are faster.', 'idle');
  }

  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', initSO5Calculator);
  } else {
    initSO5Calculator();
  }
}());
</script>

