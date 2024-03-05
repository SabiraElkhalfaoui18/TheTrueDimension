LoadPackage("Hermitian");

if not IsBound(q) then q:=5; fi;

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
R:=Y[1]*Product(Filtered(GF(q^2),c->c<>-c^q),c->Y[2]-c);

sup_pt:=Support(pt)[1];
ell:=sup_pt[1]^q*Y[1]-Y[2]-sup_pt[2]^q;
Value(ell,pt);
Valuation(PrincipalHermitian_Divisor(Hq,ell*ell^fr*ell^(fr^2)),pt);

pls:=AllRationalPlacesOfHermitian_Curve(Hq);;

T0:=ell^(q^2)/(ell^fr);
T:=(T0+T0^fr+T0^(fr^2))/3;
T1:=(ell*ell^fr*ell^(fr^2))^((q^2-1)/3)/T;;
T1_vals:=List(pls,x->Value(T1,x));
T1_mat:=DiagonalMat(T1_vals);;

rr:=Hermitian_RiemannRochSpaceBasis(s*P);;
genmat:=List(rr,x->List(pls,p->Value(x,p)));;
drr:=Hermitian_RiemannRochSpaceBasis(((q^2-1)*(q+1)/3-s)*P);;
dgenmat:=List(drr,x->List(pls,p->Value(x,p)));;
Size(rr)+Size(drr);
IsZero(genmat*T1_mat*TransposedMat(dgenmat));

