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


-- SPECIALIZED
R = QQ[y0,y1]
S = QQ[x0,x1]
F = x0^10+ x1^10 + x0^3 *x1^2
Fperp = inverseSystem F
S/Fperp
H = apply(12, k -> hilbertFunction(k,S/Fperp))
sum H

-- HOMOGNEOUES
S = QQ[x0,x1,x2]
F = x0^10+ x1^10 + x0^3 *x1^2* x2^5
Fperp = inverseSystem F
H = apply(12, k -> hilbertFunction(k,S/Fperp))
sum H



hilbertSum = (f) ->(
    Fperp := inverseSystem f;
    H := apply(15, k -> hilbertFunction(k,Fperp));
    sum H);


--Finds all monomials with a hilbert sum of 22
findMonoms = () -> (
    for i from 0 to 10  list
    (for j from 0 to 10 when j <= 30 list
        if (hilbertSum( x0^i * x1^j ) == 22)
        then print (i,j)
        else "");
    "");

-- Finds all polynomials with a hilbert sum of 22
findPolys = () -> (
    for i from 0 to 10 list
    (for j from 0 to 10 list
    (for k from 0 to 10 list
    (for l from 0 to 10 list
        if (hilbertSum( x0^i + x1^j + x0^k * x1^l ) == 22)
        then print (i,j,k,l)
        else "")));
    "");

getPoly = (a,b,c,d) -> (
        x0^a + x1^b + x0^c*x1^d
    )

-- expects numbers in order x0, x1, x0*x1
-- returns numbers in the order x0, x1, x0x1, x0x2, x1x2, x0x1x2
specialize = (a,b,c,d) -> (
        z1 := 10-a;
        z2 := 10-b;
        z3 := 10 - c - d;
        if a == 10 then l0:=10 else l0:=0;
        if b == 10 then l1:=10 else l1:=0;
        l2 := 3;
        (l0, l1,(l0, l1),(l0, l2), (l1, l2),(l0, l1, l2))
    )

getHomPoly = (a0,a1,a2,a3) -> (
        a1_0 + a1_1
    )

-- expects (1, 1, 2, 2, 2, 3) tuples
getHomPoly = ( a0, a1, a2, a3,a4,a5 ) -> (
    x0^a0 + x1^a1 + x0^(a2_0) * x1^(a2_1) + x0^(a3_0) * x2^(a3_1) + x1^(a4_0) * x2^(a4_1) + x0^(a5_0) * x1^(a5_1) *  x2^(a5_2)
    );

hilbertSumHom = (f) ->(
    S := QQ[x0,x1,x2];
    Fperp := inverseSystem f;
    H := apply(15, k -> hilbertFunction(k,Fperp));
    sum H);



-- CATALECTICANT
S = Proj(QQ[x0,x1,x2])
F = x0^5*x1^5 + x0^5 *x2^5 + x1^5*x2^5 + x0^4*x1^3*x2^3 --+ x0^3*x1^3*x2^4
Fperp = inverseSystem F
h = hilbertPolynomial(Fperp)
apply(11, k -> hilbertFunction(k,Fperp))
--sum oo

-- Specializing at x_0
S = Proj(QQ[x1,x2])
 F = x1^5 + x2^5 + x1^5*x2^5 + x1^3*x2^3
Fperp = inverseSystem F
h = hilbertPolynomial(Fperp)
apply(11, k -> hilbertFunction(k,Fperp))
--sum oo
-- Specializing at x_1
S = Proj(QQ[x0,x2])
F = x0^5 + x0^5 *x2^5 + x2^5 + x0^4*x2^3
Fperp = inverseSystem F
h = hilbertPolynomial(Fperp)
apply(11, k -> hilbertFunction(k,Fperp))
--sum oo

-- Specializing at x_1 and x_2
S = Proj(QQ[x0])
F = x0^5
Fperp = inverseSystem F
h = hilbertPolynomial(Fperp)
apply(11, k -> hilbertFunction(k,Fperp))



-- BUCHSBAUM EISENBUD
S = Proj(QQ[x0,x1,x2])
F = x0^5*x1^5 + x0^5 *x2^5 + x1^5*x2^5 + x0^4*x1^3*x2^3 --+ x0^3*x1^3*x2^4
Fperp = inverseSystem F
betti res Fperp
M = res Fperp
M.dd
M.dd_2
rank M.dd_2
--resBE M.dd_2
-- betti res I
-- resBe
--


kk=QQ[x,y,z]
J = inverseSystem(x^6+y^6+z^6)
mat = matrix{{hilbertFunction(0,J),hilbertFunction(1,J),hilbertFunction(2,J),hilbertFunction(3,J),hilbertFunction(4,J),hilbertFunction(5,J)}}
betti res J
M=res J

M.dd
K = M.dd_2
--resBE(K)

-- TO LOAD THESE METHODS:
-- load "apolar.m2"
testModulesForDeg17CY = (N,k,p) -> (
    x:=symbol x;R:=(ZZ/p)[x_0..x_6];
    numberOfGoodModules:=0;i:=0;
    usedTime:=timing while (i<N) do (
        b:=random(R^3,R^{16:-1});
        --we put SyzygyLimit=>60 because we expect
        --k<16 syzygies, so 16+28+k<=60
        fb:=res(coker b,
            DegreeLimit =>0,SyzygyLimit=>60,LengthLimit =>3);
        if rank fb_3>=k and (dim coker b)==0 then (
            fb=res(coker b, DegreeLimit =>0,LengthLimit =>4);
            if rank fb_4==0
            then numberOfGoodModules=numberOfGoodModules+1;);
        i=i+1;);
    collectGarbage();
    timeForNModules:=usedTime#0;
    {timeForNModules,numberOfGoodModules});

randomModuleForDeg17CY = (k,R) -> (
    isGoodModule:=false;i:=0;
    while not isGoodModule do (
        b:=random(R^3,R^{16:-1});
        --we put SyzygyLimit=>60 because we expect
        --k<16 syzygies, so 16+28+k<=60
        fb:=res(coker b,
            DegreeLimit =>0,SyzygyLimit=>60,LengthLimit =>3);
        if rank fb_3>=k and (dim coker b)==0 then (
            fb=res(cokEr b, DegreeLimit =>0,LengthLimit =>4);
            if rank fb_4==0 then isGoodModule=true;);
        i=i+1;);
    <<" -- Trial n. " << i <<", k="<< rank fb_3 <<endl;
    b);

skewSymMorphismsForDeg17CY = (b) -> (
    --we create a parameter ring for the morphisms:
    K:=coefficientRing ring b;
    u:=symbol u;
    U:=K[u_0..u_(binomial(6,2)-1)];
    --now we compute the equations for the u_i’s:
    UU:=U**ring b;
    equationsInUU:=flatten (substitute(b,UU)*
        substitute(genericSkewMatrix(U,u_0,6),UU));
    uu:=substitute(vars U,UU);
    equations:=substitute(
        diff(uu,transpose equationsInUU),ring b);
    syz(equations,DegreeLimit =>0));
