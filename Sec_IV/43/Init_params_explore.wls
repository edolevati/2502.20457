#!/usr/bin/env wolframscript
(* ::Package:: *)

(*metric in Boyer-Lindquist coordinates*)
ClearAll["Global`*"]
xc = {t,r,\[Theta],\[CurlyPhi]};
u = {Derivative[1][t],Derivative[1][r],Derivative[1][\[Theta]],Derivative[1][\[CurlyPhi]]};
acc = {Derivative[2][t],Derivative[2][r],Derivative[2][\[Theta]],Derivative[2][\[CurlyPhi]]};
n = Length[xc];
rPlus = M + (M^2 - a^2 M^2)^(1/2);
rMinus = M - (M^2 - a^2 M^2)^(1/2);
\[CapitalSigma][r_,\[Theta]_] = r^2 + a^2 M^2 Cos[\[Theta]]^2; (*a = dimensionless BH spin*)
\[CapitalDelta][r_] = r^2 + a^2 M^2 - 2M r;
gtt[r_,\[Theta]_] = -(1-(2 M r)/\[CapitalSigma][r,\[Theta]]);
grr[r_,\[Theta]_] = \[CapitalSigma][r,\[Theta]]/\[CapitalDelta][r];
g\[Theta]\[Theta][r_,\[Theta]_] = \[CapitalSigma][r,\[Theta]];
g\[CurlyPhi]\[CurlyPhi][r_,\[Theta]_] = (r^2 + a^2 M^2 + (2 M^3 r a^2 Sin[\[Theta]]^2)/\[CapitalSigma][r,\[Theta]])Sin[\[Theta]]^2;
gt\[CurlyPhi][r_,\[Theta]_] = -((2 M^2 r a Sin[\[Theta]]^2)/\[CapitalSigma][r,\[Theta]]);
metric = {{gtt[r,\[Theta]], 0, 0, gt\[CurlyPhi][r,\[Theta]]},{0, grr[r,\[Theta]], 0, 0},{0, 0, g\[Theta]\[Theta][r,\[Theta]], 0},{gt\[CurlyPhi][r,\[Theta]], 0, 0, g\[CurlyPhi]\[CurlyPhi][r,\[Theta]]}};
inversemetric = Simplify[Inverse[metric]];


(*eqns for orbital elements*)
rOrbitMinEq[pOrbit_,ecc_] = pOrbit / (1+ecc);
rOrbitMaxEq[pOrbit_,ecc_] = pOrbit / (1-ecc);


(*initial parameters*)
k = -1; 
M = 1.0;
\[Mu] = 1.0*10^-5;
initPar = Import["/THEORY/USERS/edoardo.levati/DATA/kick/paper/part_1/43/Init_params_explore.csv"]; (*load initial condition for a, i, ecc, p *)
(*initPar = Import["/Users/edolevati/Documents/Magistrale/Tesi/kick/paper/part_1/Init_params_explore.csv"];*) (*load initial condition for a, i, ecc, p, rp *)
nIter = Dimensions[initPar][[1]];


(*coefficients and elliptic integrals from Schmidt 0202090, e\[NotEqual]0*)
ff[r_] = r^4 + a^2 (r(r+2) + zMinus^2*(\[CapitalDelta][r*M]/M^2));
gg[r_] = 2*a*r;
hh[r_] = r(r-2) + zMinus^2*(\[CapitalDelta][r*M]/M^2)/(1 - zMinus^2);
dd[r_] = (r^2 + a^2*zMinus^2)*(\[CapitalDelta][r*M]/M^2);
KK[kk_] := NIntegrate[1/(1 - kk^2*Sin[thetaAus]^2)^(1/2), {thetaAus, 0, Pi/2}]; 
EE[kk_] := NIntegrate[(1 - kk^2*Sin[thetaAus]^2)^(1/2), {thetaAus, 0, Pi/2}];
\[CapitalPi][hhh_, kk_] := NIntegrate[1 / ((1 + hhh*Sin[thetaAus]^2)*(1 - kk^2*Sin[thetaAus]^2)^(1/2)), {thetaAus, 0, Pi/2}];


(* ***Trigger function - Ref. 2103.06306*** *)
(* ***We start the inspiral from the 4:3 resonance - Ref. 1208.3906*** *)
Csi[omegaTheta_, omegaRadial_] = Abs[4/3 - omegaTheta/omegaRadial]; 


For[j=2, j<=IntegerPart[nIter], j++, 
{
	Print["\t", j-1, " / ", nIter-1];
	a = initPar[[j, 1]]; 
	iOrbitInit = initPar[[j, 2]]; 
	eccInit = initPar[[j, 3]]; 
	pOrbitInit = initPar[[j, 4]]; 
	
	(* ***E/Lz/Q for given orbital parameters p, e, i*** *)
	rPeriapsis = rOrbitMinEq[pOrbitInit,eccInit];
	rApoapsis = rOrbitMaxEq[pOrbitInit,eccInit];
	thetaMinus = Pi/2 - iOrbitInit;
	zMinus = Cos[thetaMinus];
	(*coefficients - use dimensionless quantities*)
	ff1 = ff[rPeriapsis/M];
	gg1 = gg[rPeriapsis/M];
	hh1 = hh[rPeriapsis/M];
	dd1 = dd[rPeriapsis/M];
	ff2 = ff[rApoapsis/M];
	gg2 = gg[rApoapsis/M];
	hh2 = hh[rApoapsis/M];
	dd2 = dd[rApoapsis/M];
	kDet = dd1*hh2 - hh1*dd2;
	epsDet = dd1*gg2 - gg1*dd2;
	roDet = ff1*hh2 - hh1*ff2;
	etaDet = ff1*gg2 - gg1*ff2;
	sigmaDet = gg1*hh2 - hh1*gg2;
	(*E/Lz*)
	eInitPrograde = +((kDet*roDet + 2*epsDet*sigmaDet - 2*(sigmaDet*(sigmaDet*epsDet^2 + roDet*epsDet*kDet - etaDet*kDet^2))^(1/2)) / (roDet^2 + 4*etaDet*sigmaDet))^(1/2);
	(*eInitRetrograde = +((kDet*roDet + 2*epsDet*sigmaDet + 2*(sigmaDet*(sigmaDet*epsDet^2 + roDet*epsDet*kDet - etaDet*kDet^2))^(1/2)) / (roDet^2 + 4*etaDet*sigmaDet))^(1/2);*)
	LzInitEq[LzInit_, eInit_] = ff1*eInit^2 - 2*gg1*eInit*LzInit - hh1*LzInit^2 - dd1;
	Off[Reduce::ratnz];
	LzInitPrograde = Reduce[LzInitEq[LzInit, eInitPrograde]==0 && LzInit>=0, {LzInit}, Reals][[2]];
	(*LzInitRetrograde = Reduce[LzInitEq[LzInit, eInitRetrograde]==0 && LzInit<=0, {LzInit}, Reals][[2]];*)
	Lz = LzInitPrograde \[Mu] M;
	e = eInitPrograde \[Mu];
	Qdimless = zMinus^2*(a^2*(1 - (e/\[Mu])^2) + (Lz/(\[Mu]*M))^2/(1-zMinus^2));
	Q = Qdimless*M^2*\[Mu]^2;
	
	(* ***fundamental frequencies in Kerr - Ref. 0906.1420*** *)
	r1 = rApoapsis; (*turning point*)
	r2 = rPeriapsis; (*turning point*)
	r3 = ((2*M/(1 - (e/\[Mu])^2) - (r1 + r2)) + ((2*M/(1 - (e/\[Mu])^2) - (r1 + r2))^2 - 4*(a^2*M^2*(Q/\[Mu]^2)/((1 - (e/\[Mu])^2)*r1*r2)))^(1/2)) / 2;
	r4 = (a^2*M^2*(Q/\[Mu]^2)/((1 - (e/\[Mu])^2)*r1*r2)) / r3;
	eps0 = a^2*M^2*(1 - (e/\[Mu])^2) / (Lz/\[Mu])^2;
	zPlus = (Q/\[Mu]^2) / ((Lz/\[Mu])^2 * eps0 * zMinus^2);
	kkr = (((r1 -r2)/(r1 -r3))*((r3 - r4)/(r2 - r4)))^(1/2);
	kktheta = zMinus^2*(Lz/\[Mu])*(eps0/(Q/\[Mu]^2))^(1/2);
	(*frequencies in proper time*)
	omegaRadialProper = (Pi * ((1 - (e/\[Mu])^2)*(r1 - r3)*(r2 - r4))^(1/2)) / (2*KK[kkr]);
	omegaThetaProper = (Pi*(Lz/\[Mu])*(eps0*zPlus)^(1/2)) / (2*KK[kktheta]);
	(*frequencies in coord time*)
	hhhr = (r1 - r2)/(r1 - r3); 
	hhhPlus = ((r1 - r2)*(r3 - rPlus)) / ((r1 - r3)*(r2 - rPlus)); 
	hhhMinus = ((r1 - r2)*(r3 - rMinus)) / ((r1 - r3)*(r2 - rMinus)); 
	\[CapitalGamma] = 4*M^2*(e/\[Mu]) + (2*a^2*M^2*(e/\[Mu])*zPlus*omegaThetaProper)/(Pi*(Lz/\[Mu])*(eps0*zPlus)^(1/2))*(KK[kktheta] - EE[kktheta]) + (2*omegaRadialProper)/(Pi * ((1 - (e/\[Mu])^2)*(r1 - r3)*(r2 - r4))^(1/2)) * ((e/\[Mu])/2*((r3*(r1 + r2 + r3) - r1*r2)*KK[kkr] + (r2 - r3)*(r1 + r2 + r3 + r4)*\[CapitalPi][-hhhr, kkr] + (r1 - r3)*(r2 - r4)*EE[kkr]) + 2*M*(e/\[Mu])*(r3*KK[kkr] + (r2 - r3)*\[CapitalPi][-hhhr, kkr]) + 2*M/(rPlus - rMinus)*( (((4*M^2*(e/\[Mu]) - a*M*(Lz/\[Mu]))*rPlus - 2*M^3*a^2*(e/\[Mu]))/(r3 - rPlus)*(KK[kkr] - \[CapitalPi][-hhhPlus, kkr]*(r2 - r3)/(r2 - rPlus))) - (((4*M^2*(e/\[Mu]) - a*M*(Lz/\[Mu]))*rMinus - 2*M^3*a^2*(e/\[Mu]))/(r3 - rMinus)*(KK[kkr] - \[CapitalPi][-hhhMinus, kkr]*(r2 - r3)/(r2 - rMinus))) ));
	omegaRadialCoord = omegaRadialProper / \[CapitalGamma];
	omegaThetaCoord = omegaThetaProper / \[CapitalGamma];
	
	(* ***save data*** *)
	aData[j-1] = a;
	iOrbitData[j-1] = iOrbitInit;
	eccData[j-1] = eccInit; 
	pOrbitData[j-1] = pOrbitInit; 
	energyData[j-1] = e; 
	angularData[j-1] = Lz; 
	omegaRadialData[j-1] = omegaRadialCoord; 
	omegaThetaData[j-1] = omegaThetaCoord; 
	triggerData[j-1] = Csi[omegaThetaCoord, omegaRadialCoord];
}];


(* ***export data*** *)
data = Table[{aData[i], iOrbitData[i], eccData[i], pOrbitData[i], energyData[i], angularData[i], omegaRadialData[i], omegaThetaData[i], triggerData[i]}, {i, 1, IntegerPart[nIter]-1}];
fname = StringJoin["/THEORY/USERS/edoardo.levati/DATA/kick/paper/part_1/43/Init_params", (*ToString[j],*) ".h5"];
(*fname = StringJoin["/Users/edolevati/Documents/Magistrale/Tesi/kick/paper/part_1/Init_params", (*ToString[j],*) ".h5"];*)
Export[fname, data, OverwriteTarget -> True];
Exit[];
