LoadPackage("Hermitian");

if not IsBound(q) then q:=7; fi;

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

sup_pt:=Support(pt)[1];
ell:=sup_pt[1]^q*Y[1]-Y[2]-sup_pt[2]^q;
Value(ell,pt);
Valuation(PrincipalHermitian_Divisor(Hq,ell*ell^fr*ell^(fr^2)),pt);

rr:=Hermitian_RiemannRochSpaceBasis(s*P);;
pls:=AllRationalPlacesOfHermitian_Curve(Hq);;
#genmat:=List(rr,x->List(pls,p->Value(x,p)));;


f1:=(ell^fr/ell)^(q-1);;
f:=f1+f1^fr+f1^(fr^2);;
List(AllRationalPlacesOfHermitian_Curve(Hq),x->Order(Value(f,x)));

f1:=ell*(ell^(fr))^q;
cs:=[];
for c in GF(q^3) do
    f:=c*f1+(c*f1)^fr+(c*f1)^(fr^2);;
    if not IsZero(c) and ForAll(AllRationalAffinePlacesOfHermitian_Curve(Hq),x->IsZero(Value(f,x))) then 
        Print("hit! ",c,"\n"); 
        Add(cs,c); 
    fi;
od;
Product(cs,c->Y[1]-c);

f1/f1^(fr^2);
