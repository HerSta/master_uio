R = ZZ[x,y,z, MonomialOrder => Lex]

--f_1 = x^2 + 2*y^2 - y - 2*z
--f_2 = x^2 - 8*y^2 + 10*z - 1
--f_3 = x^2 - 7 * y * z
--
--I = ideal(f_1, f_2, f_3)
--
--G = gens gb I

-- USERGUIDE:
-- CLICK shift + f12 to open buffer
-- shift + f11 to evaluate expression


g_1 = x^2 + y^2 + z^2 - 4
g_2 = x^2 + 2* y^2 - 5
g_3 = x*z - 1 

I = ideal (g_1, g_2, g_3)

G = gens gb I





R = ZZ/101[y0,y1]
S = ZZ/101[x0,x1]
F = x0^2 + 2 * x0 * x1 +  x1^2
Fperp = inverseSystem F

S/Fperp


-- Computes the resolution
C = res Fperp

-- dd for Differentials
C.dd

hilbertFunction(0,Fperp)
hilbertFunction(1,Fperp)
hilbertFunction(2,Fperp)
hilbertFunction(3,Fperp)
hilbertFunction(4,Fperp)
hilbertFunction(5,Fperp)
hilbertSeries(Fperp)
hilbertPolynomial(Fperp)

hilbertFunction(0, S/Fperp)
hilbertFunction(1, S/Fperp)
hilbertFunction(2, S/Fperp)
hilbertFunction(3, S/Fperp)



S = ZZ/101[x1,x2,x3, x4]
f = x1^3 * x2 + x3^3 *x1 + x4^2 * x2^2
fp = inverseSystem f

h = hilbertFunction(0,S/fp)
h = hilbertFunction(1,S/fp)
h = hilbertFunction(2,S/fp)
h = hilbertFunction(3,S/fp)
h = hilbertFunction(4,S/fp)



-- Finding cactusrank 21


R = ZZ/101[y0,y1]
S = ZZ/101[x0,x1]
F = x0^1*x1^9 + x1^5
Fperp = inverseSystem F
S/Fperp

hilbertPolynomial(S/Fperp)
apply(10, k -> hilbertFunction(k,Fperp))

hilbertFunction(0,Fperp)
hilbertFunction(1,Fperp)
hilbertFunction(2,Fperp)
hilbertFunction(3,Fperp)
hilbertFunction(4,Fperp)
hilbertFunction(5,Fperp)
hilbertFunction(6,Fperp)
hilbertFunction(7,Fperp)
hilbertFunction(8,Fperp)
hilbertFunction(9,Fperp)
hilbertFunction(10,Fperp)
hilbertFunction(11,Fperp)
hilbertSeries(Fperp)

-- CATALECTICANT
S = Proj(QQ[x0,x1,x2])
F = x0^5*x1^5 + x0^5 *x2^5 + x1^5*x2^5 + x0^4*x1^3*x2^3 --+ x0^3*x1^3*x2^4
Fperp = inverseSystem F
h = hilbertPolynomial(Fperp)
apply(11, k -> hilbertFunction(k,Fperp))
sum oo

-- Specializing at x_0
S = Proj(QQ[x1,x2])
 F = x1^5 + x2^5 + x1^5*x2^5 + x1^3*x2^3
Fperp = inverseSystem F
h = hilbertPolynomial(Fperp)
apply(11, k -> hilbertFunction(k,Fperp))
sum oo
-- Specializing at x_1
S = Proj(QQ[x0,x2])
F = x0^5 + x0^5 *x2^5 + x2^5 + x0^4*x2^3
Fperp = inverseSystem F
h = hilbertPolynomial(Fperp)
apply(11, k -> hilbertFunction(k,Fperp))
sum oo

-- Specializing at x_1 and x_2
S = Proj(QQ[x0])
F = x0^5
Fperp = inverseSystem F
h = hilbertPolynomial(Fperp)
apply(11, k -> hilbertFunction(k,Fperp))
