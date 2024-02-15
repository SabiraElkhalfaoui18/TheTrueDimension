LoadPackage("Hermitian");

type:="degree_3";
res:=[];

for q in [3,5,7,9,11] do
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
m_pt:=[sup_pt,sup_pt^fr,sup_pt^(fr^2)];;
mm_1:=DiagonalMat([beta,beta^q,1]);;
mm_2:=m_pt^(-1)*mm_1*m_pt;;
mm_2*[[1,0,0],[0,0,-1],[0,-1,0]]*List(TransposedMat(mm_2),x->List(x,y->y^q)) = beta^q*[[1,0,0],[0,0,-1],[0,-1,0]];
mm:=mm_2/mm_2[1][1];
mm^(q^2-q+1);
Display(mm*[[1,0,0],[0,0,-1],[0,-1,0]]*List(TransposedMat(mm),x->List(x,y->y^q)));

aa:=Hermitian_CurveAutomorphism(Hq,mm);
NumeratorOfRationalFunction((Y[1]^(q+1)-Y[2]-Y[2]^q)^aa)/(Y[1]^(q+1)-Y[2]-Y[2]^q);

(ell/ell^fr)^aa=beta*(ell/ell^fr);
(ell^fr/ell^(fr^2))^aa=beta^(q^2)*(ell^fr/ell^(fr^2));
(ell^(fr^2)/ell)^aa=beta^(q^4)*(ell^(fr^2)/ell);

# (ell/ell^fr)^aa=beta*(ell/ell^fr);
# (ell^fr/ell^(fr^2))^aa=beta^(q-1)*(ell^fr/ell^(fr^2));
# (ell^(fr^2)/ell)^aa=beta^(-q)*(ell^(fr^2)/ell);

f:=ell^(1-q)*(ell^fr)^q/ell^(fr^2);;
(f)^aa/f;
ff:=f+f^fr+f^(fr^2);;
# (ff)^aa/ff;
Collected(List(AllRationalPlacesOfHermitian_Curve(Hq),x->Value(ff,x)));

gamma:=sup_pt[1]^(q^3+1)-sup_pt[2]-sup_pt[2]^(q^3);
u:=Z(q^6)^((1-q^4)*LogFFE(gamma,Z(q^6))/(q^3+1));
mc_1:=[[0,0,u],[u^(q^2),0,0],[0,u^(q^4),0]];;
mc:=m_pt^(-1)*mc_1*m_pt;;
mc:=mc/mc[1][1];
mm^mc = mm^(q^2);
Display(mc);
Display(mc^3); 

ac:=Hermitian_CurveAutomorphism(Hq,mc);
NumeratorOfRationalFunction((Y[1]^(q+1)-Y[2]-Y[2]^q)^ac)/(Y[1]^(q+1)-Y[2]-Y[2]^q);
Filtered(AllRationalPlacesOfHermitian_Curve(Hq),x->x=x^ac);
