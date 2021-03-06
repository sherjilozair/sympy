from sympy.physics.qho_1d import psi_n, E_n
from sympy.physics.quantum.constants import hbar
from sympy import var, simplify, integrate, oo
from sympy.core import S, pi, Rational
from sympy.functions import sqrt, exp

var("x m omega")
nu = m * omega / hbar

def test_wavefunction():
  Psi = {
    0: (nu/pi)**(S(1)/4) * exp(-nu * x**2 /2),
    1: (nu/pi)**(S(1)/4) * sqrt(2*nu) * x * exp(-nu * x**2 /2),
    2: (nu/pi)**(S(1)/4) * (2 * nu * x**2 - 1)/sqrt(2) * exp(-nu * x**2 /2),
    3: (nu/pi)**(S(1)/4) * sqrt(nu/3) * (2 * nu * x**3 - 3 * x) * exp(-nu * x**2 /2)
  }
  for n in Psi:
    assert simplify(psi_n(n, x, m, omega) - Psi[n]) == 0

def test_norm(n=1):
  # Maximum "n" which is tested:
  for i in range(n+1):
    assert integrate(psi_n(i, x, 1, 1)**2, (x,-oo,oo)) == 1

def test_orthogonality(n=1):
  # Maximum "n" which is tested:
  for i in range(n+1):
    for j in range(i+1,n+1):
      assert integrate(psi_n(i, x, 1, 1)*psi_n(j, x, 1, 1), (x,-oo,oo)) == 0

def test_energies(n=1):
  # Maximum "n" which is tested:
  for i in range(n+1):
    assert E_n(i,omega) == hbar * omega * (i + Rational(1,2))
