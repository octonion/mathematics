# The Two-Variable Factorial Conjecture via Mittag–Leffler Constellations and Fourier–Laplace Rigidity

**Christopher D. Long** · Headlamp Software, Lexington, Kentucky, USA · [galizur@gmail.com](mailto:galizur@gmail.com)

Revised research draft, August 1, 2026

📄 **[Paper (PDF)](factorial_conjecture_two_variables_unified_revised_updated.pdf)** · 📝 **[LaTeX source](factorial_conjecture_two_variables_unified_revised_updated.tex)**

*Keywords:* Factorial Conjecture, polynomial moments, Mittag–Leffler function, constellation monodromy, $E$-function, Fourier–Laplace transform.

---

## Abstract

Let $\mathcal{L}(x^{a}y^{b})=a!\,b!$ on $\mathbb{C}[x,y]$. We prove the two-variable Factorial Conjecture: if $\mathcal{L}(f^{m})=0$ for every $m\ge 1$, then $f=0$. The proof has three parts. First, algebraic specialization reduces a hypothetical counterexample to coefficients in $\overline{\mathbb{Q}}$. Second, an exact Mittag–Leffler generating period rules out every nonflat projective leading form. Its parameterized radial reversion has a computable logarithmic symbol; the positive-interval constellation of Pakovich and Muzychuk forces a primitive Fourier relation, and inverse-branch asymptotics at infinity contradict that relation. In the remaining flat-top case $f_{D}=(x+y)^{D}$, the leading singular coefficient of the normalized generating function is

$$\int_{0}^{1}\exp\!\left(\frac{f_{D-1}(t,1-t)}{D}\right)dt,$$

so vanishing of all factorial moments forces this integral to vanish. Finally, we prove that for every nonconstant $q\in\overline{\mathbb{Q}}[x]$ and distinct $A,B\in\overline{\mathbb{Q}}$, the integral $\int_{A}^{B}e^{q(x)}\,dx$ is transcendental. The arithmetic proof uses an inhomogeneous $E$-function moment system, semisimplicity of its Fourier–Laplace connection, horizontal-endomorphism rigidity in the Jacobian algebra of $q$, Beukers' refinement of the Siegel–Shidlovskii theorem, and a Cauchy-transform moment argument. Applied to $q(t)=f_{D-1}(t,1-t)/D$, this gives the final contradiction.

---

## Main results

**Theorem (Two-variable Factorial Conjecture, §1).** For every $f\in\mathbb{C}[x,y]$,

$$\mathcal{L}(f^{m})=0\ \text{ for all } m\ge 1\qquad\Longrightarrow\qquad f=0 .$$

The functional is the one introduced by van den Essen, Wright, and Zhao in their study of the Image Conjecture; the homogeneous two-variable case is due to Liu and Sun. The inhomogeneous problem does not formally reduce to the homogeneous one, because radial integration couples the homogeneous layers through Gamma factors.

**Theorem (Nonflat-top exclusion, §7).** Let $0\ne f\in\mathbb{C}[x,y]$ have exact degree $D$. If $A_{0}(t)=f_{D}(t,1-t)$ is nonconstant, then $\mathcal{L}(f^{n})\ne 0$ for infinitely many $n\ge 1$.

**Corollary (Flat-top reduction, §7).** Every hypothetical nonzero counterexample has $f_{D}(x,y)=a(x+y)^{D}$ for some $a\in\mathbb{C}^{\times}$.

**Proposition (First flat singular coefficient, §8).** In the flat-top case, the normalized generating function has leading singular coefficient

$$c_{0}(f)=\int_{0}^{1}\exp\!\left(\frac{A_{1}(t)}{D}\right)dt,\qquad\text{and}\qquad \mathcal{L}(f^{n})=0\ (n\ge 1)\ \Longrightarrow\ c_{0}(f)=0 .$$

**Theorem (Algebraic exponential-integral theorem, §9).** Let $q\in\overline{\mathbb{Q}}[x]$ be nonconstant and let $A,B\in\overline{\mathbb{Q}}$ with $A\ne B$. Then

$$\int_{A}^{B}e^{q(x)}\,dx$$

is transcendental; in particular it is nonzero.

---

## Structure of the proof

Write a nonzero $f$ of exact degree $D$ as $f=f_{D}+f_{D-1}+\cdots+f_{0}$ with $f_{j}$ homogeneous of degree $j$, and set the projective layers

$$A_{j}(t)=f_{D-j}(t,1-t),\qquad 0\le j\le D .$$

The argument follows the chain

$$\begin{aligned}
\text{counterexample over }\mathbb{C}
&\Longrightarrow \text{counterexample over }\overline{\mathbb{Q}},\\
&\Longrightarrow A_{0}(t)\ \text{is constant},\\
&\Longrightarrow f_{D}=(x+y)^{D}\ \text{after scaling},\\
&\Longrightarrow \int_{0}^{1}e^{A_{1}(t)/D}\,dt=0,\\
&\Longrightarrow \text{contradiction}.
\end{aligned}$$

1. **Specialization (§2).** A Rabinowitsch-style Nullstellensatz argument moves a hypothetical counterexample from $\mathbb{C}[x,y]$ to $\overline{\mathbb{Q}}[x,y]$ with the same monomial support.

2. **Radial–projective coordinates (§3).** The Gamma integral gives $\mathcal{L}(g)=\int_{0}^{\infty}\!\int_{0}^{\infty}g(x,y)e^{-x-y}\,dx\,dy$, and $x=rt$, $y=r(1-t)$ turns the normalized moments $I_{n}(f)=\mathcal{L}(f^{n})/\Gamma(Dn+2)$ into an exact radial–projective integral.

3. **Mittag–Leffler generating periods (§4).** With $P(t,r)=\sum_{j=0}^{D}A_{j}(t)r^{D-j}$,

   $$\mathcal{G}(s)=\sum_{n\ge 0}I_{n}s^{n}=\int_{0}^{1}\!\int_{0}^{\infty}e^{-r}r\,E_{D,2}\bigl(sP(t,r)\bigr)\,dr\,dt,$$

   $$\mathcal{J}(s)=\sum_{n\ge 0}\frac{I_{n}}{Dn+2}s^{n}=\int_{0}^{1}\!\int_{0}^{\infty}e^{-r}r\,E_{D,3}\bigl(sP(t,r)\bigr)\,dr\,dt .$$

4. **Local logarithmic symbol (§5).** Parameterized radial reversion plus a root-of-unity filter yields, near a point with $A_{0}(t_{0})\ne 0$,

   $$\mathcal{B}(w,t)=-\frac{1}{DA_{0}(t)}\exp\!\left(\frac{A_{1}(t)}{DA_{0}(t)}\right)\log\bigl(w-A_{0}(t)\bigr)+\bigl(\text{holomorphic}\bigr),$$

   so the lowest logarithmic symbol sees only $A_{0}$ and $A_{1}$.

5. **Constellations and infinity (§6–§7).** Rationality of $\mathcal{J}$ would force an incidence-orthogonality relation on the positive-interval constellation of Pakovich–Muzychuk, hence a primitive Fourier relation; inverse-branch asymptotics at infinity contradict it. Therefore $A_{0}$ is constant, i.e. $f_{D}=a(x+y)^{D}$.

6. **The flat case (§8).** In that case the singular expansion has a simple pole whose coefficient is $c_{0}(f)$, so all moments vanishing forces $\int_{0}^{1}e^{A_{1}(t)/D}\,dt=0$.

7. **Transcendence (§9–§17).** With $n=\deg q$ and the Euclidean divisions $q(x)x^{j}=h_{j}(x)q'(x)+r_{j}(x)$, the moments $U_{j}(z)=\int_{A}^{B}x^{j}e^{zq(x)}\,dx$ are $E$-functions satisfying

   $$\mathbf{U}'=\Bigl(C+\frac{D}{z}\Bigr)\mathbf{U}+\frac{1}{nz}\sum_{\lambda\in\Lambda_{\partial}}\mathbf{b}_{\lambda}e^{\lambda z},\qquad \mathrm{Spec}(D)=\Bigl\lbrace-\tfrac{1}{n},-\tfrac{2}{n},\dots,-\tfrac{n-1}{n}\Bigr\rbrace .$$

   The homogeneous part is a semisimple Fourier–Laplace connection; horizontal-endomorphism rigidity shows its endomorphisms are exactly $Z(C)\cap Z(D)$ and constant, hence the connection is multiplicity-free. Beukers' lifting theorem then forces $U_{0}(z)=\sum_{\lambda\in\Lambda}p_{\lambda}(z)e^{\lambda z}$ with $p_{\lambda}\in\mathbb{C}[z]$; a constant-coefficient annihilator $L=\prod_{\lambda\in\Lambda}(\partial_{z}-\lambda)^{N}$ produces a nonzero $H\in\mathbb{C}[t]$ with $\int_{A}^{B}q(x)^{m}H(q(x))\,dx=0$ for all $m\ge 0$, and a Cauchy-transform argument rules this out for every nonconstant phase.

8. **Completion (§18).** Taking $q(t)=f_{D-1}(t,1-t)/D$ contradicts $\int_{0}^{1}e^{q(t)}\,dt=0$: the integral is $e^{q}\ne 0$ if $q$ is constant, and transcendental otherwise.

---

## Building

```sh
pdflatex factorial_conjecture_two_variables_unified_revised_updated.tex
pdflatex factorial_conjecture_two_variables_unified_revised_updated.tex
pdflatex factorial_conjecture_two_variables_unified_revised_updated.tex
```

Requires `amsmath`, `amssymb`, `mathtools`, `amsthm`, `aliascnt`, `mathrsfs`, `microtype`, `enumitem`, `xcolor`, `url`, `hyperref`, `cleveref`, `lmodern`. Bibliography is inline (`thebibliography`); no BibTeX pass is needed. Output is 37 pages.

---

## Status

Research draft. **Not yet peer reviewed and not formally verified.**

This manuscript was developed through interactive work between the author and ChatGPT 5.6 Sol, which assisted in proof discovery, organization, symbolic checking, reference verification, adversarial auditing, and drafting. Claude Opus 5 contributed the semisimple-projector strategy, the reduction to phase-polynomial moments, and extensive adversarial audits. The AI systems are not authors. The human author bears full responsibility for all statements, proofs, citations, and any remaining errors.

---

## References

- A. Baker, *Transcendental Number Theory*, Cambridge Mathematical Library, Cambridge Univ. Press, 1990.
- F. Beukers, *A refined version of the Siegel–Shidlovskii theorem*, Ann. of Math. (2) **163** (2006), no. 1, 369–379. [doi:10.4007/annals.2006.163.369](https://doi.org/10.4007/annals.2006.163.369)
- M. A. A. de Cataldo and L. Migliorini, *The decomposition theorem, perverse sheaves and the topology of algebraic maps*, Bull. Amer. Math. Soc. (N.S.) **46** (2009), no. 4, 535–633. [doi:10.1090/S0273-0979-09-01260-9](https://doi.org/10.1090/S0273-0979-09-01260-9)
- R. Hotta, K. Takeuchi, and T. Tanisaki, *$D$-Modules, Perverse Sheaves, and Representation Theory*, Progr. Math. **236**, Birkhäuser, 2008. [doi:10.1007/978-0-8176-4523-6](https://doi.org/10.1007/978-0-8176-4523-6)
- D. Liu and X. Sun, *The Factorial Conjecture and images of locally nilpotent derivations*, Bull. Aust. Math. Soc. **101** (2020), no. 1, 71–79. [doi:10.1017/S0004972719000546](https://doi.org/10.1017/S0004972719000546)
- F. Pakovich and M. Muzychuk, *Solution of the polynomial moment problem*, Proc. Lond. Math. Soc. (3) **99** (2009), no. 3, 633–657. [doi:10.1112/plms/pdp010](https://doi.org/10.1112/plms/pdp010)
- A. van den Essen, D. Wright, and W. Zhao, *On the Image Conjecture*, J. Algebra **340** (2011), 211–224. [doi:10.1016/j.jalgebra.2011.04.036](https://doi.org/10.1016/j.jalgebra.2011.04.036)
