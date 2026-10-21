"""Audit exact bulk/component identities in the R1 continuation derivation.

These are the same twelve SymPy identities executed during derivation review.
They check tensor/component algebra, not a new physical field solution.
Use --output PATH to save a separate receipt; no original artifact is replaced
when the option is omitted. Any nonzero identity raises and exits nonzero.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import sympy as s


def _zero(checks, name, expression):
    residual = s.simplify(expression)
    if residual != 0:
        raise AssertionError(f'{name}: identity residual {residual}')
    checks.append({'name': name, 'identity_residual': '0'})


def stress_current_and_ricci(checks):
    E,P,h00,h01,h11,ho,g=s.symbols('E P H00 H01 H11 HO g3', real=True)
    eta=s.diag(-1,1,1,1)
    U=s.Matrix([E,P,0,0])
    V=eta*U
    H=s.Matrix([[h00,h01,0,0],[h01,h11,0,0],[0,0,ho,0],[0,0,0,ho]])
    S=s.trace(eta*H)
    W=H*V
    K=(V.T*H*V)[0]
    Y=(U.T*V)[0]
    T=g*(-S*(U*U.T)+U*W.T+W*U.T-eta*K)
    expected=[g*(h11*Y-2*ho*E**2),g*(h00*Y-2*ho*P**2),
              -g*K,g*(h01*Y-2*ho*E*P)]
    actual=[T[0,0],T[1,1],T[2,2],T[0,1]]
    names=['Cubic stress energy-density projection',
           'Cubic stress radial-pressure projection',
           'Cubic stress tangential-pressure projection',
           'Cubic stress mixed temporal-radial projection']
    for name,a,b in zip(names,actual,expected):
        _zero(checks,name,a-b)
    rho,pr,pt,tau,kappa=s.symbols('rho pr pt tau kappa')
    TT=s.Matrix([[rho,tau,0,0],[tau,pr,0,0],[0,0,pt,0],[0,0,0,pt]])
    RR=kappa*(TT-eta*s.trace(eta*TT)/2)
    ricci=(V.T*RR*V)[0]
    expectedR=kappa*((rho+pr+2*pt)*E**2-4*tau*E*P
                     +(rho+pr-2*pt)*P**2)/2
    _zero(checks,'Trace-reversed Ricci contraction with scalar gradient',
          ricci-expectedR)
    J=V+g*(eta*H*eta*U-S*V)
    _zero(checks,'Scalar-current temporal component',
          J[0]-(-E+g*((h11+2*ho)*E-h01*P)))
    _zero(checks,'Scalar-current radial component',
          J[1]-(P+g*(-h01*E+(h00-2*ho)*P)))


def einstein_components(checks):
    coords=s.symbols('t r theta phi', real=True)
    t,r,theta,phi=coords
    nu=s.Function('nu')(t,r)
    ga=s.Function('gamma')(t,r)
    N=s.exp(nu)
    A=s.exp(ga)
    metric=s.diag(-N**2,A**2,r**2,r**2*s.sin(theta)**2)
    inv=metric.inv()
    connection={}
    for k in range(4):
        for i in range(4):
            for j in range(4):
                v=sum(inv[k,l]*(s.diff(metric[l,j],coords[i])
                     +s.diff(metric[l,i],coords[j])
                     -s.diff(metric[i,j],coords[l])) for l in range(4))/2
                connection[k,i,j]=s.simplify(v)
    ricci={}
    for i,j in [(0,0),(1,1),(2,2),(3,3),(0,1)]:
        v=sum(s.diff(connection[k,i,j],coords[k])
             -s.diff(connection[k,i,k],coords[j])
             +sum(connection[k,i,j]*connection[l,k,l]
                  -connection[l,i,k]*connection[k,j,l] for l in range(4))
             for k in range(4))
        ricci[i,j]=s.simplify(v)
    scalar=s.simplify(sum(inv[i,i]*ricci[i,i] for i in range(4)))
    einstein={ij:s.simplify(v-metric[ij[0],ij[1]]*scalar/2)
              for ij,v in ricci.items()}
    f=A**-2
    wanted={
        (0,0):N**2*((1-f)/r**2+2*f*s.diff(ga,r)/r),
        (1,1):A**2*(-(1-f)/r**2+2*f*s.diff(nu,r)/r),
        (0,1):2*s.diff(ga,t)/r,
        (2,2):r**2*(f*(s.diff(nu,r,2)+s.diff(nu,r)**2
                       -s.diff(nu,r)*s.diff(ga,r)
                       +(s.diff(nu,r)-s.diff(ga,r))/r)
                     -N**-2*(s.diff(ga,t,2)+s.diff(ga,t)**2
                             -s.diff(nu,t)*s.diff(ga,t)))
    }
    names={(0,0):'tt',(1,1):'rr',(0,1):'tr',(2,2):'angular'}
    for ij,v in wanted.items():
        _zero(checks,'Einstein '+names[ij]+' row from exact Christoffel/Ricci construction',
              einstein[ij]-v)


def fluid_euler_projection(checks):
    t,r=s.symbols('t r', real=True)
    nu=s.Function('nu')(t,r)
    ga=s.Function('gamma')(t,r)
    b=s.Function('beta')(t,r)
    N=s.exp(nu)
    A=s.exp(ga)
    gamma=(1-b*b)**s.Rational(-1,2)
    coords=[t,r]
    metric=s.diag(-N*N,A*A)
    inv=metric.inv()
    u=s.Matrix([gamma/N,gamma*b/A])
    n=s.Matrix([gamma*b/N,gamma/A])
    connection={}
    for k in range(2):
        for i in range(2):
            for j in range(2):
                connection[k,i,j]=s.simplify(sum(inv[k,l]*(
                    s.diff(metric[l,j],coords[i])+s.diff(metric[l,i],coords[j])
                    -s.diff(metric[i,j],coords[l])) for l in range(2))/2)
    acceleration=s.Matrix([
        sum(u[i]*s.diff(u[k],coords[i]) for i in range(2))
        +sum(connection[k,i,j]*u[i]*u[j] for i in range(2) for j in range(2))
        for k in range(2)])
    projection=(n.T*metric*acceleration)[0]/gamma
    wanted=gamma**2*(s.diff(b,t)/N+b*s.diff(b,r)/A)\
           +s.diff(nu,r)/A+b*s.diff(ga,t)/N
    _zero(checks,'Perfect-fluid radial Euler acceleration projection',projection-wanted)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,help='Write a separate JSON algebra receipt.')
    args=parser.parse_args()
    checks=[]
    stress_current_and_ricci(checks)
    einstein_components(checks)
    fluid_euler_projection(checks)
    if len(checks)!=12:
        raise AssertionError('Expected exactly twelve component identities')
    result={
        'status':'PASS',
        'checked_at_utc':datetime.now(timezone.utc).isoformat(),
        'sympy_version':s.__version__,
        'scope':'Exact component algebra; no evolved physical solution or repeat of numerical lens tests.',
        'checks':checks,
        'verification_script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'physical_admission':False
    }
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        with args.output.open('w',encoding='utf-8',newline='\n') as f:
            json.dump(result,f,indent=2)
            f.write('\n')
    print(json.dumps({'status':result['status'],'identities_passed':len(checks),
                      'output':str(args.output) if args.output else None}))


if __name__=='__main__':
    main()
