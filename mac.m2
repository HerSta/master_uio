-- Here are some sample commands:
R = CC[x_0,x_1,x_2,MonomialOrder => Lex]


I = ideal(a^2 + b^2 + c^2, a^2 - b^2)
X = variety (I)

l_1 = (a + b + c)^4
l_2 = (a + sqrt(-1) *b)^4
l_3 = (a + b + sqrt(-1) * sqrt(2) * c)^4
l_4 = (a - b + sqrt(-1) * sqrt(2) * c)^4
l_5 = (a + b - sqrt(-1) * sqrt(2) * c)^4
l_6 = (a - b - sqrt(-1) * sqrt(2) * c)^4


Y = l_1 + l_2 + l_3 + l_4 + l_5+ l_6



K = ideal (a^2, c^2)
KX = variety K

I = ideal(x_0^2 + x_1^2 + x_2^2 + x_0*x_1 + x_0*x_2 + x_1*x_2)
J = ideal(x_0 + x_1 + x_2)
int = intersect(I, J)


A = ideal(x_0^2 + x_1^2 + x_2^2 + x_0*x_1 + x_0*x_2 + x_1*x_2, x_0 + x_1 + x_2)
g = gens gb A
V = variety A
degree V
dim V
V = variety ideal(a^2 + b^2 + c^2 + a*b + a*c + b*c, a+b+c)
