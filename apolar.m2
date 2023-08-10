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


--g_1 = x^2 + y^2 + z^2 - 4
--g_2 = x^2 + 2* y^2 - 5
--g_3 = x*z - 1

--I = ideal (g_1, g_2, g_3)

--G = gens gb I





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
--F = x0^9+ x1^7 + x0^8 *x1^1
F = x0^8*x1^2+x0^9+x1^9+x0^3*x1^3
F = x0^10+x0^5*x1^5+x0^8+x1^8+x0^3*x1^3
F = x0^5*x1^5+x0^9+x1^9+x0^6+x0^3*x1^3
F = x0^9 + x0^3*x1^4 + x1^7 + x0^4*x1^2
F = x1^10 + x0^3*x1^4 + x0^6 + x0^4*x1
--F = x0^3*x2^3 + x0^6*x2^4 + x0^4*x2^5
--F = x1^10 + x1^4*x2^3 + x2^4 +x1*x2^5
--F = x0^8* x1 +x0^3 * x1^4+x0^6+x1^5
--F = x0^8 *x2+x0^3 * x2^3+x0^6* x2^4+x2^5
--F =  x1 *x2+ x1^4* x2^3+x2^4+x1^5 *x2^5
Fperp = inverseSystem F
S/Fperp
H = apply(11, k -> hilbertFunction(k,S/Fperp))
sum H

-- HOMOGNEOUES
S = QQ[x0,x1,x2]
--F = x0^9*x2+ x1^7*x2^3 + x0^8 *x1^1*x2
F = x0^5*x1^5+x0^9*x2+x1^9*x2+x0^3*x1^3*x2^4+x2^10
F = x0^10+x0^5*x1^5+x0^8*x2^2+x1^8*x2^2+x0^3*x1^3*x2^4+x2^10
F = x0^5*x1^5+x0^9*x2+x1^9*x2+x0^6*x2^4+x0^3*x1^3*x2^4+x2^10+3
F =  x0^9*x2 + x0^3*x1^4*x2^3 + x1^7*x2^3 + x0^4*x1^2*x2^4 + x2^10
F = (x1^10 + x0^3*x1^4*x2^3 + x0^6*x2^4 + x0^4*x1*x2^5 )*x2^(17)
F = x0^8* x1 *x2+x0^3 * x1^4* x2^3+x0^6* x2^4+x1^5 *x2^5
Fperp = inverseSystem F
betti Fperp
H = apply(12, k -> hilbertFunction(k,S/Fperp))
sum H



hilbertSum = (f) ->(
    Fperp := inverseSystem f;
    H := apply(11, k -> hilbertFunction(k,Fperp));
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
findPolysPrint = () -> (
    for i from 0 to 10 list -- 0 to 10
    (for j from 0 to 10 list -- 0 to 10
    (for k from 1 to 9 list
    (for l from 1 to 9 list
        if (hilbertSum( x0^i + x1^j + x0^k * x1^l ) == 22)
        then print (i,j,k,l)
        else "")));
    "");



findPolyLists = (d) -> (
        numMonoms := d+1;
        polys := ();
        for i from 0 to numMonoms list
        (
            l := ();
            for j from 0 to numMonoms list
            (
                    if (i == j) then l = append(l,1) else l=append(l,0);

            );
            polys = append(polys, l);
        );
    polys);


getPolyFromList = l -> (
        apply(l, p -> getSinglePolyFromList toList(p))
    );


getSinglePolyFromList = (l) -> (
        var := l_0 * x0^(10) +
            l_1 * x0^9 *x1^1 +
            l_2 * x0^8 *x1^2 +
            l_3 * x0^7 *x1^3 +
            l_4 * x0^6 *x1^4 +
            l_5 * x0^5 *x1^5 +
            l_6 * x0^4 *x1^6 +
            l_7 * x0^3 *x1^7 +
            l_8 * x0^2 *x1^8 +
            l_9 * x0^1 *x1^9 +
            l_10 * x0^0 *x1^(10);
        --for i from 0 to 10 list(
        --        mon := poly(l_i * x0^(10-i) * x1^(i));
        --        var = append(var, mon);
        --    );
        var);

-- input deg d
findPolys = (d, g) -> (
    var := ();
    numPartDiff := ceiling (binomial(2+d, 2)/3);
    for i from 0 to d list -- 0 to 10
    (for j from 0 to d list -- 0 to 10
    (for k from 0 to d list -- 1 to 9
    (for l from 0 to (d-k) list -- 1 to 9
        if (hilbertSum(x0^i + x1^j + x0^k * x1^l + g_0 + g_1) ==  numPartDiff)
        then var = append(var, (i,j,k,l))
        else "")));
    var);

getPoly = (a,b,c,d) -> (
        x0^a + x1^b + x0^c*x1^d
    )

-- expects numbers in order x0, x1, x0*x1
-- returns numbers in the order x0, x1, x0x1, x0x2, x1x2, x0x1x2
specialize = A -> (
        a:= A_0;
        b:= A_1;
        c:= A_2;
        d:= A_3;
        l0:=0;
        l1:=0;
        l2:=0;
        l3:=0;
        l4:=0;
        l5:=0;
        l6:=0;
        l7:=0;
        l8:=0;
        l9:=0;
        l10:=0;
        z1 := 10-a;
        z2 := 10-b;
        z3 := 10 - c - d;
        if a == 10 then l0=10;
        if b == 10 then l1=10;
        if c+d == 10 then (l2=c; l3=d);
        if a < 10 then (l4=a; l5=z1);
        if b < 10 then (l6=b; l7=z2);
        if c+d < 10 then(
            l8=c;
            l9=d;
            l10=z3;);
        (l0, l1, (l2, l3), (l4, l5), (l6, l7), (l8, l9, l10))
    )


-- INPUT
 --(0, 0, (0, 0), (3, 7), (4, 6), (2, 3, 5))
getHomPoly = (a, g) -> (
    x0^(a_0) +
    x1^(a_1) +
    x0^((a_2)_0) * x1^((a_2)_1) +
    x0^((a_3)_0) * x2^((a_3)_1) +
    x1^((a_4)_0) * x2^((a_4)_1) +
    x0^((a_5)_0) * x1^((a_5)_1) * x2^((a_5)_2) +
     g_0 * x2^(10 - (degree g_0)_0) + g_1 * x2^(10 - (degree g_1)_0) - 3 -x2^(10)

    );
-- expects (1, 1, 2, 2, 2, 3) tuples
getHomPoly2 = ( a0, a1, a2, a3,a4,a5 ) -> (
    x0^a0 + x1^a1 + x0^(a2_0) * x1^(a2_1) + x0^(a3_0) * x2^(a3_1) + x1^(a4_0) * x2^(a4_1) + x0^(a5_0) * x1^(a5_1) *  x2^(a5_2)
    );

hilbertSumHom = (f) ->(
    S := QQ[x0,x1,x2];
    Fperp := inverseSystem f;
    H := apply(15, k -> hilbertFunction(k,Fperp));
    sum H);

maxHilbert = (f) ->(
    Fperp := inverseSystem f;
    H := apply(15, k -> hilbertFunction(k,Fperp));
    --print H;
    max H);

sys = () -> (
        d := 10;
        S = QQ[x0,x1,x2];
        polsL := findPolyLists(d);
        pols := getPolyFromList polsL;
        print("number of pols:" | toString(length pols));
        homs := apply(pols, s -> homogenize( s, x2));
        print("number of homs:" | toString(length pols));
        --print special;
        --print homs;
        num := length pols - 1;
        for i from 0 to num list (
                --print( "polynomial:" );
                --print homs_i;
                ----print homs_i;
                --print( "biggest hilbert element:" | toString(maxHilbert(homs_i)));
                print maxHilbert(homs_i);
                --
                if (maxHilbert(homs_i) == 21) then print "ELEMENT FOUND!!";
                --if (maxHilbert(homs_i) == 21) then print homs_i;

                --if (maxHilbert(homs_i) == 20) then print homs_i;
                if (maxHilbert(homs_i) == 19) then print homs_i;
            )
    );

combo = () -> (
        d := 10;
        S = QQ[x0,x1];
        pols := findPolys(d);
        print("number of pols:" | toString(length pols));
        S = QQ[x0,x1,x2];
        special := apply(pols, p -> specialize p);
        homs := apply(special, s -> getHomPoly s);
        print("number of homs:" | toString(length pols));
        --print special;
        --print homs;
        num := length pols - 1;
        for i from 0 to num list (
                --print( "polynomial:" );
                --print homs_i;
                ----print homs_i;
                --print( "biggest hilbert element:" | toString(maxHilbert(homs_i)));
                --print maxHilbert(homs_i);
                if (maxHilbert(homs_i) == 21) then print toString(homs_i);

            )
    );


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


S = Proj(QQ[x0,x1,x2])
F = x1^(10)+ x1^4 * x2^3+x2^4+x1* x2^5 -- YAY
Fperp = inverseSystem F
FperpHom = homogenize(Fperp,x0)
print toString(FperpHom)
apply(11, k -> hilbertFunction(k,Fperp))
apply(11, k -> hilbertFunction(k,FperpHom))

S = Proj(QQ[x0,x2])
F = x0^3 * x2^3+x0^6 *x2^4+x0^4* x2^5 -- YAY
Fperp = inverseSystem F
--h = hilbertPolynomial(Fperp)
apply(11, k -> hilbertFunction(k,Fperp))

S = Proj(QQ[x0,x1])
F = x1^(10)+x0^3 * x1^4 +x0^6 +x0^4 *x1
Fperp = inverseSystem F
--h = hilbertPolynomial(Fperp)
apply(11, k -> hilbertFunction(k,Fperp))


-- BUCHSBAUM EISENBUD
S = QQ[x0,x1,x2]
--F = x0^5*x1^5 + x0^5 *x2^5 + x1^5*x2^5 + x0^4*x1^3*x2^3 --+ x0^3*x1^3*x2^4
F = x0^8* x1* x2+x0^3* x1^4 * x2^3+x0^6* x2^4+x1^5* x2^5
F = x0^8* x2^2+x1^8 *x2^2+x0^3* x1^4* x2^3+x0^4* x1* x2^5 -- VERY CLOSE
F= x0^9* x2+x0^3* x1^4* x2^3+x1^6* x2^4+x0^4* x1* x2^5
F = x1^(10)+x0^3 * x1^4 * x2^3+x0^6 *x2^4+x0^4 *x1* x2^5 -- YAY
Fperp = inverseSystem F
betti res Fperp
M = res Fperp
M.dd
B = M.dd_2
rank M.dd_2
--loadPackage "ResLengthThree"
A=resLengthThreeAlg M -- define the algebra
netList multTableOneTwo A -- this presents the multiplication table in readable form
-- we need to cut out the first row and column to get a matrix corresponding to the multiplication
-- {1..13} removes first row and column
some = ((matrix((multTableOneTwo(A))_{1..13}))_{1..13})
-- g_1=>1 replaces placeholder values with 1
H=sub(some, g_1=>1)
X=transpose(H)*B -- this is I think the skew-symmetric matrix corresponding to B after the change of basis described above
--print X
--print (X*50*50*9*2*500)
--print (X*50 * (1/27))

-- Extracting matrix to the left of zero-block
subM = X^{7,8,9,10,11,12}_{0,1,2,3,4,5,6}
myIdeal = minors(6,subM)
use S;
myIdeal = substitute(myIdeal, S)
v = variety myIdeal
dim v
degree v


--resBE(K)

-- TO LOAD THESE METHODS:
-- load "apolar.m2"




--R=ZZ/5[x1,x2,z3,x4]
--T=random(R^5, R^{5:-1})
--N=T-transpose T
--I=pfaffians(4,N)
--J=resolution I -- I start with a pfaffian resolution
--B=J.dd_2 -- here I lost skew symmetry of the middle matrix

-- ANNES EKSEMPEL!
kk=QQ[x,y,z]
Fperp = inverseSystem(x^6+y^6+z^6)
betti res Fperp
J=res Fperp
J.dd
B = J.dd_2
A=resLengthThreeAlg J -- define the algebra
netList multTableOneTwo A -- this presents the multiplication table in readable form
-- we need to cut out the first row and column to get a matrix corresponding to the multiplication
H=sub(((matrix((multTableOneTwo(A))_{1..5}))_{1..5}), g_1=>1)
X=transpose(H)*B -- this is I think the skew-symmetric matrix corresponding to B after the change of basis described above

-- NOW I NEED TO MAKE A BASIS CHANGE BACK to the original ring kk
subM = X^{0,1}_{2,3,4}
--myIdeal = pfaffians(2,subM)
myIdeal = minors(2,subM)
use kk
myIdeal= substitute(myIdeal, kk)
v = variety myIdeal
dim v
degree v






-- Full Procedure
--
--
F = x1^(10)+x0^3 * x1^4 * x2^3+x0^6 *x2^4+x0^4 *x1* x2^5 -- YAY
polyHasDim22 = (F) -> (
    Fperp := inverseSystem F;
    M := res Fperp;
    B := M.dd_2;
    A := resLengthThreeAlg M;
    netList multTableOneTwo A;
    some := ((matrix((multTableOneTwo(A))_{1..13}))_{1..13});
    H:=sub(some, g_1=>1);
    X:=transpose(H)*B;
    zeroMatrix := matrix{{0,0,0,0,0,0,0},{0,0,0,0,0,0,0},{0,0,0,0,0,0,0},{0,0,0,0,0,0,0},{0,0,0,0,0,0,0},{0,0,0,0,0,0,0}};
    subM := X^{7..12}_{0..6};
    if (subM ==  zeroMatrix) then
    (
        print "Zero block found!";
        print F;
        myIdeal := minors(6,subM);
        use S;
        myIdeal := substitute(myIdeal, S);
        v := variety myIdeal;
        if (degree v == 22) then return true;
    )
    else false;
    false);

proc = () -> (
        d := 10;
        S = QQ[x0,x1];
        for j from 3 to (d-1) list
        (
            for k from 4 to (d-1) list
            (
                spice := {x0^(j)*x1^(k) , x0^(k)*x1^5};
                pols := findPolys(d, spice);
                print("number of polynomials with natural rank 22 to analyze:" | toString(length pols));
                S = QQ[x0,x1,x2];
                g := {x0^(j)*x1^(k) , x0^k*x1^5};
                special := apply(pols, p -> specialize p);
                homs := apply(special, s -> getHomPoly(s, g));
                num := length pols - 1;
                for i from 0 to num list (
                    if (maxHilbert(homs_i) == 21) then
                    (
                        print "This polynomial has cactus rank ?, natural rank 22, and catalecticant rank 21:";
                        print toString(homs_i);
                        -- Any polynomial in here has natural rank 22 and catalecticant rank 21
                        -- Now we check if the cactus rank is 22 or 21
                        if (polyHasDim22 homs_i) then
                        (
                            print "This polynomial has cactus rank 22, natural rank 22, and catalecticant rank 21:";
                            print homs_i;
                            );
                        );
                    );
              );
        );
    );
