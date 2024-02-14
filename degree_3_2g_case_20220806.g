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
	s:=2*g-1;
	while true do
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

LoadPackage("Hermitian");
q:=5;

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
s:=2*g;
s=(q-1)*(q+1)-(q-1);

rr:=Hermitian_RiemannRochSpaceBasis(s*P);;
pts:=AllRationalAffinePlacesOfHermitian_Curve(Hq);;
cmat:=List(rr,r->List(pts,p->Value(r,p)));;
hcode:=VectorSpace(GF(q^2),cmat);
rcode:=RestrictVectorSpace(hcode,GF(r));

sol:=List(Basis(rcode),x->SolutionMat(cmat,x));
bas:=sol*rr;
ForAll(bas,r->ForAll(pts,p->Value(r,p) in GF(q)));

ell1:=Support(pt)[1];
ell1:=[ell1[2]*ell1[3]^(q^2)-ell1[3]*ell1[2]^(q^2),-ell1[1]*ell1[3]^(q^2)+ell1[3]*ell1[1]^(q^2),ell1[1]*ell1[2]^(q^2)-ell1[2]*ell1[1]^(q^2)];
ell1:=ell1*[Y[1],Y[2],1];
ell1*ell1^fr*ell1^(fr^2);


