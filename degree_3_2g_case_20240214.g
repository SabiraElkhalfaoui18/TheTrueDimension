LoadPackage("Hermitian");

type:="degree_3";
res:=[];

for q in [3,5,7,9] do
	Print("##################################\n");
	Print("### q = ",q," , "," q^2 "," = ",StringPP(q^2),"\n");
	Print("##################################\n");
	Y:=HermitianIndeterminates(GF(q^2),"Y1","Y2");
	Hq:=Hermitian_Curve(Y[1]);
	###
	pt:=RandomPlaceOfGivenDegreeOfHermitian_Curve(Hq,3);
	fr:=FrobeniusAutomorphismOfHermitian_Curve(Hq);
	P:=Sum(AC_FrobeniusAutomorphismOrbit(fr,pt));
	###
	n:=q^3;  # length of the code
	g:=q*(q-1)/2;   # genus
	r:=q;   # size of the subfield
	for s in [2*g-1,2*g,2*g+1] do
		degG:=s;
		# if degG+1-g>=n then break; fi;
		if IsBound(hcode) and Dimension(hcode)=n then break; fi;
		hcode:=Hermitian_FunctionalCode(s*P);
		rcode:=RestrictVectorSpace(hcode,GF(r));
		Add(res,[q,r,s,n,degG,Dimension(hcode),Dimension(rcode)]);
		Print("# degG = ", 3*degG,", ",Dimension(rcode), "\n");
		s:=s+1;
	od;
	Print("Done for q = ",q,". \n");
od;

####

sup_pt:=Support(pt)[1];
ell:=sup_pt[1]^q*Y[1]-Y[2]-sup_pt[2]^q;
Value(ell,pt);
Valuation(PrincipalHermitian_Divisor(Hq,ell*ell^fr*ell^(fr^2)),pt);

beta:=Z(q^6)^((q^3-1)*(q+1));
mm:=[sup_pt,sup_pt^fr,sup_pt^(fr^2)];;
mm:=mm^(-1)*DiagonalMat([beta,beta^q,1])*mm;;
mm:=mm/mm[1][1];
mm^(q^2-q+1);
aa:=Hermitian_CurveAutomorphism(Hq,mm);

hcurveaut_on_fnct:=function(f,a)
    local y,b;
    y:=IndeterminatesOfHermitianRatFunc(f);
    if y=[] then return f; fi;
    b:=[y[1],y[2],1]*(a!.mat^(-1));
    return Value(f,y,[b[1]/b[3],b[2]/b[3]]);
end;

NumeratorOfRationalFunction(hcurveaut_on_fnct(Y[1]^(q+1)-Y[2]-Y[2]^q,aa))/(Y[1]^(q+1)-Y[2]-Y[2]^q);

hcurveaut_on_fnct(ell/ell^fr,aa)=beta*(ell/ell^fr);
hcurveaut_on_fnct(ell^fr/ell^(fr^2),aa)=beta^(q^2)*(ell^fr/ell^(fr^2));
hcurveaut_on_fnct(ell^(fr^2)/ell,aa)=beta^(q^4)*(ell^(fr^2)/ell);

hcurveaut_on_fnct(ell/ell^fr,aa)=beta*(ell/ell^fr);
hcurveaut_on_fnct(ell^fr/ell^(fr^2),aa)=beta^(q-1)*(ell^fr/ell^(fr^2));
hcurveaut_on_fnct(ell^(fr^2)/ell,aa)=beta^(-q)*(ell^(fr^2)/ell);

f:=ell^(1-q)*(ell^fr)^q/ell^(fr^2);;
hcurveaut_on_fnct(f,aa)/f;
ff:=f+f^fr+f^(fr^2);;
#hcurveaut_on_fnct(ff,aa)/ff;
List(AllRationalPlacesOfHermitian_Curve(Hq),x->Value(ff,x));;
List(Set(last),x->Number(last,y->x=y));
