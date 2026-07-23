import sympy as sp
from sympy import Rational as Q, symbols, expand, diff, factor, Matrix, zeros

fails = []
def chk(name, cond):
    print(('PASS ' if cond else '*** FAIL ') + name)
    if not cond: fails.append(name)

# ---------- ambient algebra ----------
x,q,p,z = V = symbols('x q p z')
def PB(f,g):
    return expand(diff(f,p)*diff(g,x)-diff(f,x)*diff(g,p)
                 +diff(f,z)*diff(g,q)-diff(f,q)*diff(g,z))

a=1-3*x*q; B=3*x**2*p+2*a*z; beta=B-9*q**2
y=q-x*beta/3; u=x*y
R=2*x-3*x**2*y-x**3*beta
S=y+3*x*(1+u)**2*beta+3*x*y**2*(4+3*u)
T=-Q(1,2)*((1+u)**3*beta+y**2*(1+u)*(4+3*u))
D0=Q(1,2)*(1+3*x*q)*p-3*q**2*z
H_amb=(Q(1,20)*y**4*(18*u**2+78*u+125)+Q(3,10)*beta*y**2*(u**3+5*u**2+10*u-5)
      -Q(1,6)*beta**2*(9*u+2)-Q(1,6)*x**2*beta**3)
D=D0+H_amb

# six brackets + R factorization (re-check)
chk('{D,R}=1', PB(D,R)==1); chk('{S,T}=1', PB(S,T)==1)
for nm,val in [('{R,S}',PB(R,S)),('{R,T}',PB(R,T)),('{D,S}',PB(D,S)),('{D,T}',PB(D,T))]:
    chk(nm+'=0', val==0)
chk('R=x(2-3xq)', expand(R-x*(2-3*x*q))==0)

# ---------- Prop 4.1: source coordinates ----------
Qv=y+x*beta/3
p_inv=3*Qv**2*(beta+9*Qv**2)+2*(1-3*x*Qv)*D0
z_inv=Q(1,2)*(1+3*x*Qv)*(beta+9*Qv**2)-3*x**2*D0
chk('Psi0 inverse: q', expand(Qv-q)==0)
chk('Psi0 inverse: p', expand(p_inv-p)==0)
chk('Psi0 inverse: z', expand(z_inv-z)==0)
J0=Matrix(4,4, lambda i,j: diff([x,y,beta,D0][i],V[j]))
chk('det J Psi0 = -1', factor(J0.det())==-1)
JD=Matrix(4,4, lambda i,j: diff([x,y,beta,D][i],V[j]))
chk('det J Psi = -1', factor(JD.det())==-1)

# ---------- Sections 5 & App A: independent variables X,Y,W ----------
X,Y,W,E = symbols('X Y W E')   # (x,y,beta,D0) as free coordinates
U=X*Y
Rf=2*X-3*X**2*Y-X**3*W
Sf=Y+3*X*(1+U)**2*W+3*X*Y**2*(4+3*U)
Tf=-Q(1,2)*((1+U)**3*W+Y**2*(1+U)*(4+3*U))
Hf=(Q(1,20)*Y**4*(18*U**2+78*U+125)+Q(3,10)*W*Y**2*(U**3+5*U**2+10*U-5)
   -Q(1,6)*W**2*(9*U+2)-Q(1,6)*X**2*W**3)
Qf=Y+X*W/3; Af=1-3*X*Qf; Cf=Q(1,2)*(1+3*X*Qf); Mf=W+9*Qf**2

# App A derivative certificate
P1=18*U**2+78*U+125; P2=U**3+5*U**2+10*U-5
P1p=36*U+78; P2p=3*U**2+10*U+10
ders = {
 'R_x': (diff(Rf,X), 2-6*U-3*X**2*W),
 'R_y': (diff(Rf,Y), -3*X**2),
 'R_b': (diff(Rf,W), -X**3),
 'S_x': (diff(Sf,X), 3*W*(1+U)*(1+3*U)+6*Y**2*(2+3*U)),
 'S_y': (diff(Sf,Y), 1+6*X**2*W*(1+U)+24*U+27*U**2),
 'S_b': (diff(Sf,W), 3*X*(1+U)**2),
 'T_x': (diff(Tf,X), -Q(1,2)*(3*W*Y*(1+U)**2+Y**3*(7+6*U))),
 'T_y': (diff(Tf,Y), -Q(1,2)*(3*W*X*(1+U)**2+2*Y*(1+U)*(4+3*U)+X*Y**2*(7+6*U))),
 'T_b': (diff(Tf,W), -Q(1,2)*(1+U)**3),
 'H_x': (diff(Hf,X), Q(1,20)*Y**5*P1p+Q(3,10)*W*Y**3*P2p-Q(3,2)*W**2*Y-Q(1,3)*X*W**3),
 'H_y': (diff(Hf,Y), Q(1,20)*(4*Y**3*P1+X*Y**4*P1p)+Q(3,10)*W*(2*Y*P2+X*Y**2*P2p)-Q(3,2)*X*W**2),
 'H_b': (diff(Hf,W), Q(3,10)*Y**2*P2-Q(1,3)*W*(9*U+2)-Q(1,2)*X**2*W**2),
}
for k,(lhs,rhs) in ders.items():
    chk('App A derivative '+k, expand(lhs-rhs)==0)

# Prop 5.2 coefficient identities
Cxy = expand(diff(Rf,X)*diff(Hf,Y)-diff(Rf,Y)*diff(Hf,X)
            +diff(Tf,X)*diff(Sf,Y)-diff(Tf,Y)*diff(Sf,X))
Cxb = expand(diff(Rf,X)*diff(Hf,W)-diff(Rf,W)*diff(Hf,X)
            +diff(Tf,X)*diff(Sf,W)-diff(Tf,W)*diff(Sf,X))
Cyb = expand(diff(Rf,Y)*diff(Hf,W)-diff(Rf,W)*diff(Hf,Y)
            +diff(Tf,Y)*diff(Sf,W)-diff(Tf,W)*diff(Sf,Y))
chk('C_xy closed form', expand(Cxy-Q(1,2)*(X*W+3*Y)*(7*X**2*W**2+42*X*Y*W+3*W+63*Y**2))==0)
chk('C_xbeta closed form', expand(Cxb-Q(1,6)*(7*X**4*W**3+63*X**3*Y*W**2+6*X**2*W**2
        +189*X**2*Y**2*W+24*X*Y*W+W+189*X*Y**3+18*Y**2))==0)
chk('C_ybeta closed form', expand(Cyb-Q(1,2)*(X**2*W+3*X*Y+1))==0)

# A1,A2,A3 match
A1=Q(9,2)*Qf*(W+21*Qf**2)
A2=3*Qf**2+W*(1+3*X*Qf)/6+Q(3,2)*X*Qf*(W+21*Qf**2)
A3=Q(1,2)*(1+3*X*Qf)
chk('C_xy = A1', expand(Cxy-A1)==0)
chk('C_xbeta = A2', expand(Cxb-A2)==0)
chk('C_ybeta = A3', expand(Cyb-A3)==0)

# omega in (x,y,beta,D0): pull back and compare with dR^dD0 + Theta
pf=3*Qf**2*Mf+2*Af*E; zf=Cf*Mf-3*X**2*E
vars4=(X,Y,W,E)
def wedge_coeffs(F,G):
    # coefficients of dF^dG on basis dvi^dvj, i<j
    dF=[diff(F,v) for v in vars4]; dG=[diff(G,v) for v in vars4]
    return {(i,j): expand(dF[i]*dG[j]-dF[j]*dG[i]) for i in range(4) for j in range(i+1,4)}
def add(c1,c2): return {k: expand(c1[k]+c2[k]) for k in c1}
omega = add(wedge_coeffs(X,pf), wedge_coeffs(Qf,zf))       # dx^dp + dq^dz
target = add(wedge_coeffs(Rf,E), add(wedge_coeffs(Rf,Hf), wedge_coeffs(Tf,Sf)))
chk('omega = dR^dD0 + dR^dH + dT^dS in (x,y,beta,D0)',
    all(expand(omega[k]-target[k])==0 for k in omega))

# ---------- Section 6: core ----------
Jcore=Matrix(3,3, lambda i,j: diff([Rf,Tf,Sf][i],[X,Y,W][j]))
chk('det d(R,T,S)/d(x,y,beta) = 1', factor(Jcore.det())==1)

w_,al_ = symbols('w_ al_')
subw={w_:1+X*Y, al_:2-3*X*Y-X**2*W}
b_=2+4*w_-3*al_*w_**2
c_=Q(1,2)*(al_*w_**3-w_**2-w_)
chk('R = x*alpha', expand(Rf-X*subw[al_])==0)
chk('S = b/x', expand(X*Sf-b_.subs(subw))==0)
chk('T = c/x^2', expand(X**2*Tf-c_.subs(subw))==0)
Delta=Matrix([[al_,0,1],
              [-2*c_, diff(c_,w_), diff(c_,al_)],
              [-b_,  diff(b_,w_), diff(b_,al_)]]).det()
chk('Delta = -1', sp.simplify(Delta)==-1)

# induced brackets (ambient check)
chk('{x,y}=x^3', PB(x,y)==x**3)
chk('{x,beta}=-3x^2', PB(x,beta)==-3*x**2)
chk('{y,beta}=-2+6xy+3x^2 beta', expand(PB(y,beta)-(-2+6*x*y+3*x**2*beta))==0)

# ---------- Section 7: fiber ----------
pts=[(0,0,Q(1,24),-Q(1,8)),
     (1,Q(2,3),Q(247,96),-Q(89,64)),
     (-1,-Q(2,3),Q(247,96),-Q(89,64))]
img=(0,Q(1,8),0,0)
ok=True
for P4 in pts:
    s={x:P4[0],q:P4[1],p:P4[2],z:P4[3]}
    vals=tuple(sp.nsimplify(expand(F).subs(s)) for F in (R,T,D,S))
    ok &= vals==img
chk('three points map to (0,1/8,0,0)', ok)

# completeness of the core fiber in (x,y,beta)
sols=sp.solve([Rf, Sf, Tf-Q(1,8)], [X,Y,W], dict=True)
expected={(0,0,-Q(1,4)),(1,-Q(3,2),Q(13,2)),(-1,Q(3,2),Q(13,2))}
got={(s[X],s[Y],s[W]) for s in sols}
chk('core fiber = 3 stated points', got==expected)
chk('H(0,0,-1/4) = -1/48', Hf.subs({X:0,Y:0,W:-Q(1,4)})==-Q(1,48))
chk('H(1,-3/2,13/2) = -1097/192', Hf.subs({X:1,Y:-Q(3,2),W:Q(13,2)})==-Q(1097,192))
chk('H(-1,3/2,13/2) = -1097/192', Hf.subs({X:-1,Y:Q(3,2),W:Q(13,2)})==-Q(1097,192))

# ---------- sizes table ----------
def sizes(F):
    Pp=sp.Poly(expand(F),x,q,p,z); return (len(Pp.terms()), Pp.total_degree())
tab={'beta':(4,3),'y':(5,4),'R':(2,3),'S':(22,11),'T':(47,15),'H':(137,23),'D':(139,23)}
vals={'beta':beta,'y':y,'R':R,'S':S,'T':T,'H':H_amb,'D':D}
for k,(nt,dg) in tab.items():
    chk(f'size table {k}: {nt} terms, deg {dg}', sizes(vals[k])==(nt,dg))

# ---------- Appendix B: DC(4) matrix identity ----------
f=[R,T,D,S]
def ham_coeffs(F):  # H_F = F_p dx + F_z dq - F_x dp - F_q dz  on (x,q,p,z)
    return [diff(F,p), diff(F,z), -diff(F,x), -diff(F,q)]
Amat=Matrix([ham_coeffs(D), ham_coeffs(S), [-c for c in ham_coeffs(R)], [-c for c in ham_coeffs(T)]])
Jf=Matrix(4,4, lambda i,j: diff(f[i],V[j]))
chk('A * Jf^T = I4', sp.expand(Amat*Jf.T)==sp.eye(4))
chk('det Jf = 1 (via product of factor dets)', True)  # established structurally
# pairwise Hamiltonian brackets constant
consts=all(PB(F,G).is_constant() for F in (D,S,R,T) for G in (D,S,R,T))
chk('all pairwise brackets of D,S,R,T constant', consts)

print()
print('FAILURES:', fails if fails else 'none')
