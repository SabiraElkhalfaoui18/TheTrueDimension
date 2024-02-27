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

###

R:=Y[1]*Product(Filtered(GF(q^2),c->c<>-c^q),c->Y[2]-c);
dR:=Derivative(R,Y[1])+Y[1]^q*Derivative(R,Y[2]);
dR_values:=List(AllRationalAffinePlacesOfHermitian_Curve(Hq),x->Value(dR,x));;

B:=function(s,i)
    local u,v,f;
    if 0=s mod (q+1) then 
        u:=s/(q+1);
        v:=0;
    else
        u:=Int(s/(q+1))+1;
        v:=u*(q+1)-s;
    fi;
    f:=ell^(2*u-v)*(ell^fr)^(v-u)/(ell^(fr^2))^u;
    if i>0 then f:=f^(fr^i); fi;
    return f;
end;

for s in [q-1..q^2-q+1] do 
    Print(s,"\t",Value(B(s,0),Q_infty),"\n"); 
od;

###

hform:=function(u,v) return u[1]*v[1]^q-u[2]*v[3]^q-u[3]*v[2]^q; end;
xi:=hform([Y[1],Y[2],1],sup_pt)/hform(sup_pt^(fr^2),sup_pt);
Sum([0,1,2],i->xi^(fr^(i+1))*sup_pt^(fr^i));
Sum([0,1,2],i->(ell*(ell^fr)^q/hform(sup_pt,sup_pt^fr)^q)^(fr^i));