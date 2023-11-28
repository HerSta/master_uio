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
     g_0 * x2^(10 - (degree g_0)_0) + g_1 * x2^(10 - (degree g_1)_0) - 3

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
        num := length pols - 1;
        for i from 0 to num list (
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
loadPackage "ResLengthThree"
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


-- GRASSMANIAN APPROACH
S = QQ[x0,x1,x2]
F = x1^(10)+x0^3 * x1^4 * x2^3+x0^6 *x2^4+x0^4 *x1* x2^5 -- YAY
F = x0^8* x1* x2+x0^3* x1^4 * x2^3+x0^6* x2^4+x1^5* x2^5
Fperp = inverseSystem F
betti res Fperp
M = res Fperp
M.dd
B = M.dd_2

-- simple
S = QQ[x0,x1,x2]
m = 9; -- number of rows
n = 9; --number of columns
T=random(S^m, S^{n:-1}) -- {a:-b} means a cols, degree b
B=T-transpose T

numVars = 4*5
R = QQ[c,a_0..a_numVars, MonomialOrder=>Lex]
eqM = matrix{{1,0,0,0,a_0,  a_1,  a_2,  a_3, a_4},
             {0,1,0,0,a_5,  a_6,  a_7,  a_8, a_9},
             {0,0,1,0,a_10, a_11, a_12, a_13, a_14},
             {0,0,0,1,a_15, a_16, a_17, a_18, a_19}}
B0 = sub(sub(sub(sub(B, x1=>0), x2=>0), x0=>1), QQ)
B1 = sub(sub(sub(sub(B, x0=>0), x2=>0), x1=>1), QQ)
B2 = sub(sub(sub(sub(B, x1=>0), x0=>0), x2=>1), QQ)

B0 + B1 + B2  -- Checking that decomp works
eqs = {}
for i from 0 to 3 list
(
    for j from (i+1) to 3 list
    (
        eqs = append(eqs,eqM^{i} * B0 * transpose(eqM^{j}));
        eqs = append(eqs,eqM^{i} * B1 * transpose(eqM^{j}));
        eqs = append(eqs,eqM^{i} * B2 * transpose(eqM^{j}));
    )
)
eqs

I = ideal(eqs)
dim I
v = variety(I)
degree v
dim v

-- Starting with matrices
S = ZZ/2[x0,x1,x2]
T=random(S^13, S^{13:-1}) -- {a:-b} means a cols, degree b
S = QQ[x0,x1,x2]

B = matrix({
    {0,x0,0,0,0,0,0,0,0,0,0,0,x1},
    {-x0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,x2},
    {-x1,0,0,0,0,0,0,0,0,0,0,-x2,0}
    })

assert(B + transpose(B) == 0)

B0 = sub(sub(sub(sub(B, x1=>0), x2=>0), x0=>1), QQ)
B1 = sub(sub(sub(sub(B, x0=>0), x2=>0), x1=>1), QQ)
B2 = sub(sub(sub(sub(B, x1=>0), x0=>0), x2=>1), QQ)

R = QQ[b, a_0..a_41]
eqM = matrix{
    {1,0,0,0,0,0,a_0, a_1, a_2, a_3, a_4, a_5, a_6 },
    {0,1,0,0,0,0,a_7, a_8, a_9, a_10,a_11,a_12,a_13},
    {0,0,1,0,0,0,a_14,a_15,a_16,a_17,a_18,a_19,a_20},
    {0,0,0,1,0,0,a_21,a_22,a_23,a_24,a_25,a_26,a_27},
    {0,0,0,0,1,0,a_28,a_29,a_30,a_31,a_32,a_33,a_34},
    {0,0,0,0,0,1,a_35,a_36,a_37,a_38,a_39,a_40,a_41}
    }
eqs = {};
for i from 0 to 5 list
(
    for j from (i+1) to 5 list
    (
        eqs = append(eqs,eqM^{i} * B0 * transpose(eqM^{j}));
        eqs = append(eqs,eqM^{i} * B1 * transpose(eqM^{j}));
        eqs = append(eqs,eqM^{i} * B2 * transpose(eqM^{j}));
    )
)

eqs
I = ideal(eqs)
dim I



rank M.dd_2

A=resLengthThreeAlg M
netList multTableOneTwo A
some = ((matrix((multTableOneTwo(A))_{1..13}))_{1..13})
H=sub(some, g_1=>1)
X=transpose(H)*B



Y = substitute(X, R);

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

loadPackage "ResLengthThree"
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
    potZero := X^{7..12}_{7..12};
    if (potZero ==  0) then
    (
        subM := X^{7..12}_{0..6};
        print "Zero block found!";
        print F;
        myIdeal := minors(6,subM);
        use S;
        myIdeal = substitute(myIdeal, S);
        v := variety myIdeal;
        if (degree v == 22) then return true;
    )
    else false;
    false);

proc = () -> (
        d := 10;
        S = QQ[x0,x1];
        for j from 2 to (d-1) list
        (
            for k from 2 to (d-1) list
            (
                spice := {x0^(j)*x1^(k) ,x0^3*x1^4};
                pols := findPolys(d, spice);
                print("number of polynomials with natural rank 22 to analyze:" | toString(length pols));
                S = QQ[x0,x1,x2];
                g := {x0^(j)*x1^(k) , x0^3*x1^4};
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



K = QQ[a0,a1,a2,a3,a4]
R = K[x0,x1]
F = a0*x0^4 + 4*a1*x0^3*x1 + 6*a2*x0^2*x1^2 + 4*a3*x0*x1^3 + a4*x1^4
F = x0^4 + 4*x0^3*x1 + 6*x0^2*x1^2 + 4*x0*x1^3 + x1^4
fperp = inverseSystem(F)
A = matrix{
        {a0,a1,a2},
        {a1,a2,a3},
        {a2,a3,a4}
        }
b = adjoint(A, R^3, R^3)


R = QQ[x0,x1,x2, MonomialOrder => Lex];
F = x0^2*x1^2 + x0^2*x2^2 + x1^2 * x2^2 + x0^2 * x1^1 * x2^1
F = x0^2*x1^2 + x0^2*x2^2
Fx0 = x1^2 + x2^2
fperp = inverseSystem(F)
fperpx0 = inverseSystem(Fx0)



-- F binary deg 4

R = QQ[a0,a1,a2,a3,a4]
S = frac R
A = matrix{
    {a0,a1,a2},
    {a1,a2,a3},
    {a2,a3,a4}
    }

co = inverse A * det A

use R
I = co_0_2 - co_1_1
J = ideal substitute(I,R)
    variety J
dim J
degree J
degree variety J
dim variety J


-- deg 4
R = QQ[a0,a1,a2,a3,a4]
S = frac R
A = matrix{
    {a0,a1,a2},
    {a1,a2,a3},
    {a2,a3,a4}
    }

co = inverse A * det A

I = ideal (co_0_2 - co_1_1)
J = substitute(I,R)
pdI = primaryDecomposition J
L = length pdI
dim variety pdI_0
degree variety pdI_0
dim variety pdI_1
degree variety pdI_1

pd0 = substitute(pdI_0, R)
pd1 = substitute(pdI_1, R)
h = substitute(ideal det A, R)
dim variety h
degree variety h

-- F binary deg 6

R = QQ[a0,a1,a2,a3,a4,a5,a6]
S = frac R
A = matrix{
    {a0,a1,a2,a3},
    {a1,a2,a3,a4},
    {a2,a3,a4,a5},
    {a3,a4,a5,a6}
    }


co = inverse A * det A

A = matrix{
    {a0,a1,a2,a3},
    {a1,a2,a3,a0},
    {a2,a3,a0,a1},
    {a3,a0,a1,a2}
    }

co = inverse A * det A


I = ideal (co_0_2 - co_1_1, co_0_3 - co_1_2,
        co_1_3 - co_2_2,
    co_2_3 - co_3_2)
J = substitute(I,R)
pdI = primaryDecomposition J
L = length pdI
m = minors(3,A)
dim variety pdI_0
degree variety pdI_0
dim variety pdI_1
degree variety pdI_1

pd0 = substitute(pdI_0, R)
pd1 = substitute(pdI_1, R)
h = substitute(ideal det A, R)
dim variety h
degree variety h

v = variety (pd0 + h)
dim v
degree v

v = variety (pd1 + h)
dim v
degree v

inter0 = intersect(h, pd0)
inter1 = intersect(h, pd1)
dim variety ideal inter_0
degree variety ideal inter_0
dim variety ideal inter_1
degree variety ideal inter_1

dim J
variety J
degree variety J
dim variety J

-- deg 8
R = QQ[a0,a1,a2,a3,a4,a5,a6,a7,a8]
S = frac R
A = matrix{
    {a0,a1,a2,a3,a4},
    {a1,a2,a3,a4,a5},
    {a2,a3,a4,a5,a6},
    {a3,a4,a5,a6,a7},
    {a4,a5,a6,a7,a8}
    }

co = inverse A * det A

I = ideal (co_0_2 - co_1_1, co_0_3 - co_1_2, co_0_4 - co_1_3,
            co_1_3 - co_2_2, co_1_4 - co_2_3,
            co_2_4 - co_3_3)
J = substitute(I,R)
pdI = primaryDecomposition J
L = length pdI
dim variety pdI_0
degree variety pdI_0
dim variety pdI_1
degree variety pdI_1

dim J
variety J
degree variety J
dim variety J

-- deg 10
R = QQ[a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10]
S = frac R
A = matrix{
    {a0,a1,a2,a3,a4,a5},
    {a1,a2,a3,a4,a5,a6},
    {a2,a3,a4,a5,a6,a7},
    {a3,a4,a5,a6,a7,a8},
    {a4,a5,a6,a7,a8,a9},
    {a5,a6,a7,a8,a9,a10}
    }

co = inverse A * det A;


I = ideal (co_0_2 - co_1_1, co_0_3 - co_1_2, co_0_4 - co_1_3, co_0_5 - co_1_4,
            co_1_3 - co_2_2, co_1_4 - co_2_3, co_1_5 - co_2_4,
            co_2_4 - co_3_3, co_2_5 - co_3_4,
            co_3_5 - co_4_4)
J = substitute(I,R)
pdI = primaryDecomposition J
L = length pdI
dim variety pdI_0
degree variety pdI_0
dim variety pdI_1
degree variety pdI_1
dim variety pdI_2
degree variety pdI_2
dim variety pdI_3
degree variety pdI_3
dim variety pdI_4
degree variety pdI_4
dim variety pdI_5
degree variety pdI_5
dim J
variety J
degree variety J
dim variety J


    -- deg 10
R = QQ[a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12]
S = frac R
A = matrix{
    {a0,a1,a2,a3,a4,a5,a6},
    {a1,a2,a3,a4,a5,a6,a7},
    {a2,a3,a4,a5,a6,a7,a8},
    {a3,a4,a5,a6,a7,a8,a9},
    {a4,a5,a6,a7,a8,a9,a10},
    {a5,a6,a7,a8,a9,a10,a11},
    {a6,a7,a8,a9,a10,a11,a12}
    }

co = inverse A * det A

I = ideal (co_0_2 - co_1_1, co_0_3 - co_1_2, co_0_4 - co_1_3, co_0_5 - co_1_4,
            co_1_3 - co_2_2, co_1_4 - co_2_3, co_1_5 - co_2_4,
            co_2_4 - co_3_3, co_2_5 - co_3_4,
            co_3_5 - co_4_4)
J = substitute(I,R)
pdI = primaryDecomposition J
L = length pdI
dim variety pdI_0
degree variety pdI_0
dim variety pdI_1
degree variety pdI_1
dim variety pdI_2
degree variety pdI_2
dim variety pdI_3
degree variety pdI_3
dim variety pdI_4
degree variety pdI_4
dim variety pdI_5
degree variety pdI_5


-- pfaffians(2,A)
-- Cofac = minors(2,A)
