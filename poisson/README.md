# An Explicit Counterexample to the Rank-Two Poisson Conjecture

## Abstract

Let

$$
\mathcal{P}_2=\mathbb{C}[x,q,p,z]
$$

carry the canonical Poisson bracket determined by

$$
\{p,x\}=\{z,q\}=1
$$

and by the vanishing of the other brackets between distinct generators. Here and throughout, “rank two” means two canonical pairs in the standard indexing of the canonical Poisson algebras; thus there are four polynomial generators and the Poisson tensor has geometric rank four. We give explicit polynomials

$$
R,T,D,S\in\mathbb{Q}[x,q,p,z]
$$

satisfying

$$
\{D,R\}=1,\qquad \{S,T\}=1,
$$

and

$$
\{R,S\}=\{R,T\}=\{D,S\}=\{D,T\}=0,
$$

while

$$
R=x(2-3xq).
$$

Consequently, the assignment

$$
(x,q,p,z)\mapsto(R,T,D,S)
$$

defines a Poisson endomorphism of $\mathcal{P}_2$ that is not an automorphism. This disproves the Poisson Conjecture for two canonical pairs, and hence for every number of canonical pairs at least two. The associated polynomial map of $\mathbb{A}^4$ preserves the canonical symplectic form, has Jacobian determinant one, and has an explicit fiber consisting of exactly three points. The proof uses a polynomial source coordinate system in which the symplectic identity reduces to three displayed coefficient identities. A separate appendix uses the same four output polynomials and their Hamiltonian duals to construct an explicit nonautomorphic endomorphism of the fourth Weyl algebra.

## Repository contents

### Source

- [`explicit_rank_two_poisson_counterexample_updated.pdf`](explicit_rank_two_poisson_counterexample_updated.pdf)
- [`explicit_rank_two_poisson_counterexample_updated.tex`](explicit_rank_two_poisson_counterexample_updated.tex)

### Verification scripts

- [`verify_rank2_poisson_claude_audit.py`](verify_rank2_poisson_claude_audit.py)
- [`verify_rank2_poisson_sparse.py`](verify_rank2_poisson_sparse.py)
- [`verify_rank2_poisson_counterexample.py`](verify_rank2_poisson_counterexample.py)
- [`verify_rank2_poisson_sympy.py`](verify_rank2_poisson_sympy.py)

### Output text

- [`verify_rank2_poisson_claude_audit.txt`](verify_rank2_poisson_claude_audit.txt)
- [`verify_rank2_poisson_sparse.txt`](verify_rank2_poisson_sparse.txt)
- [`verify_rank2_poisson_counterexample.txt`](verify_rank2_poisson_counterexample.txt)
- [`verify_rank2_poisson_sympy.txt`](verify_rank2_poisson_sympy.txt)
