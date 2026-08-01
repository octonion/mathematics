# Algebraic Exponential Integrals: Transcendence and Linear Independence

**Christopher D. Long** · Headlamp Software, Lexington, Kentucky, USA · [galizur@gmail.com](mailto:galizur@gmail.com)

August 2026

📄 **[Paper (PDF)](algebraic_exponential_integrals_simplified.pdf)** · 📝 **[LaTeX source](algebraic_exponential_integrals_simplified.tex)**

*2020 MSC:* Primary 11J81; Secondary 14F10, 34M35.
*Keywords:* $E$-functions, Siegel–Shidlovskii theorem, exponential integrals, polynomial moments, transcendence, linear independence.

---

## Abstract

Let $q\in\overline{\mathbb{Q}}[x]$ be nonconstant, let $\xi\in\overline{\mathbb{Q}}^{\times}$, and let

$$\Delta=\sum_{\alpha} d_{\alpha}[\alpha]$$

be a finite degree-zero zero-chain with $d_{\alpha},\alpha\in\overline{\mathbb{Q}}$. The entire function

$$I_{q,\Delta}(z)=\int_{\Delta} e^{zq(x)}\,dx$$

is defined by integration over any one-chain with boundary $\Delta$. We prove the algebraic-value dichotomy

$$I_{q,\Delta}(\xi)\in\overline{\mathbb{Q}}\quad\Longleftrightarrow\quad I_{q,\Delta}(z)\equiv 0 .$$

Thus every nonzero algebraic-parameter value of an algebraic exponential integral is transcendental. More generally, all $\overline{\mathbb{Q}}$-linear relations among finitely many such values, together with $1$, are exactly the specializations of their functional relations. This yields linear independence, for example, for exponential integrals from a common endpoint to algebraic points whose $q$-values are pairwise distinct regular values.

For $\deg q\ge 2$, the proof uses an explicit inhomogeneous differential system for the moment functions

$$U_{j}(z)=\int_{\Delta} x^{j} e^{zq(x)}\,dx .$$

A complete-reducibility theorem for its homogeneous part, combined with an elementary local-spectrum argument, implies that every stable rational relation space is obtained from a constant vector space. Beukers' refinement of the Siegel–Shidlovskii theorem then converts an algebraic value into an exponential-polynomial identity. A constant-coefficient annihilator and an edgewise Cauchy-transform argument force the original exponential period to vanish identically.

---

## Main results

**Theorem 1.2 (Algebraic-value rigidity).** Let $q\in\overline{\mathbb{Q}}[x]$ be nonconstant, let $\Delta$ be a finite algebraic zero-chain of degree zero, and let $\xi\in\overline{\mathbb{Q}}^{\times}$. Then

$$I_{q,\Delta}(\xi)\in\overline{\mathbb{Q}}\quad\Longleftrightarrow\quad I_{q,\Delta}(z)\equiv 0,$$

and in the equivalent cases $I_{q,\Delta}(\xi)=0$.

**Corollary 1.3 (Algebraic exponential-integral theorem).** For nonconstant $q\in\overline{\mathbb{Q}}[x]$, distinct $A,B\in\overline{\mathbb{Q}}$, and $\xi\in\overline{\mathbb{Q}}^{\times}$, the number

$$\int_{A}^{B} e^{\xi q(x)}\,dx$$

is transcendental; in particular it is nonzero.

**Theorem 1.4 (Transfer of linear relations).** With $I_{j}=I_{q,\Delta_{j}}$ and $a_{0},\dots,a_{r}\in\overline{\mathbb{Q}}$,

$$a_{0}+\sum_{j=1}^{r} a_{j}I_{j}(\xi)=0\quad\Longleftrightarrow\quad a_{0}=0\ \text{ and }\ \sum_{j=1}^{r} a_{j}I_{j}(z)\equiv 0 .$$

Consequently

$$\dim_{\overline{\mathbb{Q}}}\mathrm{span}_{\overline{\mathbb{Q}}}\lbrace 1,I_{1}(\xi),\dots,I_{r}(\xi)\rbrace=1+\dim_{\overline{\mathbb{Q}}}\mathrm{span}_{\overline{\mathbb{Q}}}\lbrace I_{1},\dots,I_{r}\rbrace,$$

the right-hand span being taken in the space of entire functions.

**Corollary 1.5 (Linear independence over separated regular values).** If $q(a_{0}),q(a_{1}),\dots,q(a_{r})$ are pairwise distinct regular values of $q$, then for every $\xi\in\overline{\mathbb{Q}}^{\times}$ the numbers

$$1,\quad \int_{a_{0}}^{a_{1}} e^{\xi q(x)}\,dx,\quad\dots,\quad \int_{a_{0}}^{a_{r}} e^{\xi q(x)}\,dx$$

are linearly independent over $\overline{\mathbb{Q}}$.

Functional relations do occur and cannot be dropped from Theorem 1.4: if $q$ is even, then

$$\int_{a}^{b} e^{zq(x)}\,dx+\int_{-a}^{-b} e^{zq(x)}\,dx=0$$

identically in $z$.

---

## Structure of the proof

Write $n=\deg q\ge 2$ and, for $0\le j\le n-2$, perform the Euclidean divisions $q(x)x^{j}=h_{j}(x)q'(x)+r_{j}(x)$ with $\deg r_{j}\le n-2$. Let $\mathsf{C},\mathsf{D}$ be the matrices whose $j$-th columns are the coefficient vectors of $r_{j}$ and $-h_{j}'$, and put $C=\mathsf{C}^{\mathsf{T}}$, $D=\mathsf{D}^{\mathsf{T}}$.

1. **Moment system (§3).** The vector $\mathbf{U}=(U_{0},\dots,U_{n-2})^{\mathsf{T}}$ satisfies

   $$\mathbf{U}'=\Bigl(C+\frac{D}{z}\Bigr)\mathbf{U}+\frac{1}{nz}\sum_{\lambda\in\Lambda_{\Delta}} \mathbf{b}_{\lambda}e^{\lambda z},\qquad \mathrm{Spec}(D)=\Bigl\lbrace -\tfrac{1}{n},-\tfrac{2}{n},\dots,-\tfrac{n-1}{n}\Bigr\rbrace .$$

2. **Arithmetic (§4).** Each $U_{j}$ is an $E$-function, with the explicit joint-height bound $h(u_{j,0},\dots,u_{j,m})=O(m)$.

3. **Structure (§5).** The row system $\boldsymbol{\phi}\mapsto\boldsymbol{\phi}'+\boldsymbol{\phi}A(z)$, $A=C+D/z$, is completely reducible — via the decomposition theorem, Riemann–Hilbert, and the generic Fourier transform of $q_{+}\mathcal{O}_{\mathbb{A}^{1}_{x}}$. Combined with the local spectrum of the compressed operator at a critical point $\tau$ of multiplicity $m_{\tau}=\mathrm{ord}_{\tau} q'$,

   $$-\frac{m_{\tau}}{m_{\tau}+1},\ -\frac{m_{\tau}-1}{m_{\tau}+1},\ \dots,\ -\frac{1}{m_{\tau}+1},$$

   every $\mathcal{T}$-stable subspace $\mathcal{R}\subseteq\mathbb{C}(z)^{1\times(n-1)}$ has the form $\mathcal{R}=W\otimes_{\mathbb{C}}\mathbb{C}(z)$ for a constant $(C,D)$-invariant $W$.

4. **Value to exponential polynomial (§6).** Beukers' lifting theorem plus constant splitting forces

   $$U_{0}(z)=\sum_{\lambda\in\Lambda} p_{\lambda}(z)e^{\lambda z},\qquad p_{\lambda}\in\mathbb{C}[z].$$

5. **Annihilation and rigidity (§7–§8).** A constant-coefficient operator $L=\prod_{\lambda\in\Lambda}(\partial_{z}-\lambda)^{N}$ yields a nonzero $H\in\mathbb{C}[t]$ with vanishing phase moments $\int_{\Gamma} q(x)^{m}H(q(x))\,dx=0$ for all $m\ge 0$. An edgewise Cauchy-transform argument on a polynomial constellation then kills every edge density, giving $I_{q,\Delta}\equiv 0$ — and, as a by-product, the exact edge-density criterion for functional vanishing.

---

## Building

```sh
pdflatex algebraic_exponential_integrals_simplified.tex
pdflatex algebraic_exponential_integrals_simplified.tex
pdflatex algebraic_exponential_integrals_simplified.tex
```

Requires `amsmath`, `amssymb`, `mathtools`, `amsthm`, `aliascnt`, `mathrsfs`, `microtype`, `enumitem`, `hyperref`, `cleveref`, `lmodern`. Bibliography is inline (`thebibliography`); no BibTeX pass is needed. Output is 20 pages.

---

## Status

Preprint. **Not yet peer reviewed and not formally verified.**

This manuscript was developed through interactive work between the author and ChatGPT 5.6 Pro, which assisted in proof discovery, organization, symbolic checking, reference verification, adversarial auditing, and drafting. Claude Opus 5 contributed the complete-reducibility/projector strategy, the reduction to phase-polynomial moments, and extensive adversarial audits. The AI systems are not authors. The human author bears full responsibility for all statements, proofs, citations, and any remaining errors.

---

## References

- B. Adamczewski and T. Rivoal, *Exceptional values of $E$-functions at algebraic points*, Bull. Lond. Math. Soc. **50** (2018), no. 4, 697–708. [doi:10.1112/blms.12168](https://doi.org/10.1112/blms.12168)
- A. Baker, *Transcendental Number Theory*, Cambridge Mathematical Library, Cambridge Univ. Press, 1990.
- F. Beukers, *A refined version of the Siegel–Shidlovskii theorem*, Ann. of Math. (2) **163** (2006), no. 1, 369–379. [doi:10.4007/annals.2006.163.369](https://doi.org/10.4007/annals.2006.163.369)
- A. Bostan, T. Rivoal, and B. Salvy, *Minimization of differential equations and algebraic values of $E$-functions*, Math. Comp. **93** (2024), no. 347, 1427–1472. [doi:10.1090/mcom/3912](https://doi.org/10.1090/mcom/3912)
- M. A. A. de Cataldo and L. Migliorini, *The decomposition theorem, perverse sheaves and the topology of algebraic maps*, Bull. Amer. Math. Soc. (N.S.) **46** (2009), no. 4, 535–633. [doi:10.1090/S0273-0979-09-01260-9](https://doi.org/10.1090/S0273-0979-09-01260-9)
- É. Delaygue, *A Lindemann–Weierstrass theorem for $E$-functions*, J. Reine Angew. Math. **820** (2025), 75–85. [doi:10.1515/crelle-2024-0090](https://doi.org/10.1515/crelle-2024-0090)
- R. Hotta, K. Takeuchi, and T. Tanisaki, *$D$-Modules, Perverse Sheaves, and Representation Theory*, Progr. Math. **236**, Birkhäuser, 2008. [doi:10.1007/978-0-8176-4523-6](https://doi.org/10.1007/978-0-8176-4523-6)
- Yu. V. Nesterenko and A. B. Shidlovskii, *Linear independence of values of $E$-functions*, Mat. Sb. **187** (1996), no. 8, 93–108; Sb. Math. **187** (1996), no. 8, 1197–1211. [doi:10.1070/SM1996v187n08ABEH000152](https://doi.org/10.1070/SM1996v187n08ABEH000152)
- F. Pakovich and M. Muzychuk, *Solution of the polynomial moment problem*, Proc. Lond. Math. Soc. (3) **99** (2009), no. 3, 633–657. [doi:10.1112/plms/pdp010](https://doi.org/10.1112/plms/pdp010)
- F. Pakovich, N. Roytvarf, and Y. Yomdin, *Cauchy-type integrals of algebraic functions*, Israel J. Math. **144** (2004), 221–291. [doi:10.1007/BF02916714](https://doi.org/10.1007/BF02916714)
