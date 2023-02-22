kk

K := frac(QQ[a0, a1, a2, a3, a4])
R := K[x0, x1]


f := a0 * x0^4 + 4 * a1 * x0^3 * x1 + 6 * a2 * x0^2 * x1^2 + 4 * a3 * x0 * x1^3 + a4 * x1^4


fpolar = inverseSystem f


m = matrix{{a0, a1, a2}, {a1, a2, a3}, {a2, a3, a4}}
det2 = det m



ex = ((((a0+2)*(a0+1) + 2) * (a0 + 2)*(a0 + 1)) / 4 - (2 * a0 +1)*(2 * a0 +2)) / 2
exfac = factor(ex)