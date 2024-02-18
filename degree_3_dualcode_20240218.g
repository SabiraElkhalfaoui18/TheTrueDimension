LoadPackage("Hermitian");

if not IsBound(q) then q:=11; fi;

Y:=HermitianIndeterminates(GF(q^2),"Y1","Y2");
Hq:=Hermitian_Curve(Y[1]);
###
pt:=RandomPlaceOfGivenDegreeOfHermitian_Curve(Hq,3);
fr:=FrobeniusAutomorphismOfHermitian_Curve(Hq);
P:=Sum(AC_FrobeniusAutomorphismOrbit(fr,pt));
Q_infty:=Hermitian_Place(Hq,[infinity]);
###
n:=q^3;  # length of the code
g:=q*(q-1)/2;   # genus
s:=2*g;

D:=Sum(AllRationalAffinePlacesOfHermitian_Curve(Hq));
hcode:=Hermitian_FunctionalCode(s*P,D);

R:=Y[1]*Product(Filtered(GF(q^2),c->c<>-c^q),c->Y[2]-c);
dR:=Derivative(R,Y[1])+Y[1]^q*Derivative(R,Y[2]);
dR_values:=List(AllRationalAffinePlacesOfHermitian_Curve(Hq),x->Value(dR,x));;
ForAny(dR_values,IsZero);

hcode_predual:=Hermitian_FunctionalCode((q^3+2*g-2)*Q_infty-s*P,D);
a1:=Random(hcode);; a2:=Random(hcode_predual);; Sum([1..q^3],i->a1[i]*a2[i]/dR_values[i]);

###

hcode_plus:=Hermitian_FunctionalCode(s*P,D+Q_infty);
rcode:=RestrictVectorSpace(hcode,GF(q));
[Length(Zero(rcode)),Dimension(rcode)];
rcode_plus:=RestrictVectorSpace(hcode_plus,GF(q));
[Length(Zero(rcode_plus)),Dimension(rcode_plus)];