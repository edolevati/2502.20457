#!/usr/bin/env wolframscript
(* ::Package:: *)

(*metric in Boyer-Lindquist coordinates*)
ClearAll["Global`*"]
Print["-> Geodesics"];
xc = {t,r,\[Theta],\[CurlyPhi]};
u = {Derivative[1][t],Derivative[1][r],Derivative[1][\[Theta]],Derivative[1][\[CurlyPhi]]};
acc = {Derivative[2][t],Derivative[2][r],Derivative[2][\[Theta]],Derivative[2][\[CurlyPhi]]};
nn = Length[xc];
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
affine = Simplify[Table[(1/2)*Sum[(inversemetric[[i,s]])*(D[metric[[s,j]],xc[[l]]]+D[metric[[s,l]], xc[[j]]]-D[metric[[j,l]],xc[[s]]]),{s,1,nn}], {i,1,nn}, {j,1,nn}, {l,1,nn}]];


(*geodesic equations*)
geodesic = Simplify[-affine . u . u];


(*momentum eqns for t,\[CurlyPhi] from stationarity and axisimmetry*)
killing1 = {1,0,0,0};
killing2 = {0,0,0,1};
metric2 = {{gtt, 0, 0, gt\[CurlyPhi]},{0, grr, 0, 0},{0, 0, g\[Theta]\[Theta], 0},{gt\[CurlyPhi], 0, 0, g\[CurlyPhi]\[CurlyPhi]}};
isoTime = -metric2 . killing1 . u;
isoRot = metric2 . killing2 . u;
sol = Solve[{isoTime==e/\[Mu], isoRot==Lz/\[Mu]}, {Derivative[1][t], Derivative[1][\[CurlyPhi]]}];
momentumEqt = Derivative[1][t]/.sol[[1]][[1]];
momentumEq\[CurlyPhi] = Derivative[1][\[CurlyPhi]]/.sol[[1]][[2]];
momentumEqtCalc[r_,\[Theta]_] = ReplaceAll[momentumEqt, {g\[CurlyPhi]\[CurlyPhi]-> g\[CurlyPhi]\[CurlyPhi][r,\[Theta]], gtt-> gtt[r,\[Theta]], gt\[CurlyPhi]-> gt\[CurlyPhi][r,\[Theta]]}];
momentumEq\[CurlyPhi]Calc[r_,\[Theta]_] = ReplaceAll[momentumEq\[CurlyPhi], {g\[CurlyPhi]\[CurlyPhi]-> g\[CurlyPhi]\[CurlyPhi][r,\[Theta]], gtt-> gtt[r,\[Theta]], gt\[CurlyPhi]-> gt\[CurlyPhi][r,\[Theta]]}];


(*transformation to coordinate time*)
coefft = Simplify[ReplaceAll[geodesic[[1]]/(Derivative[1][t])^2, {Derivative[1][r]->Derivative[1][r]*momentumEqtCalc[r,\[Theta]], Derivative[1][\[Theta]]->Derivative[1][\[Theta]]*momentumEqtCalc[r,\[Theta]], Derivative[1][\[CurlyPhi]]->Derivative[1][\[CurlyPhi]]*momentumEqtCalc[r,\[Theta]], Derivative[1][t]->momentumEqtCalc[r,\[Theta]]}]];
coeffr = Simplify[ReplaceAll[-geodesic[[2]]/(Derivative[1][t])^2, {Derivative[1][r]->Derivative[1][r]*momentumEqtCalc[r,\[Theta]], Derivative[1][\[Theta]]->Derivative[1][\[Theta]]*momentumEqtCalc[r,\[Theta]], Derivative[1][\[CurlyPhi]]->Derivative[1][\[CurlyPhi]]*momentumEqtCalc[r,\[Theta]], Derivative[1][t]->momentumEqtCalc[r,\[Theta]]}]];
coeff\[Theta] = Simplify[ReplaceAll[-geodesic[[3]]/(Derivative[1][t])^2, {Derivative[1][r]->Derivative[1][r]*momentumEqtCalc[r,\[Theta]], Derivative[1][\[Theta]]->Derivative[1][\[Theta]]*momentumEqtCalc[r,\[Theta]], Derivative[1][\[CurlyPhi]]->Derivative[1][\[CurlyPhi]]*momentumEqtCalc[r,\[Theta]], Derivative[1][t]->momentumEqtCalc[r,\[Theta]]}]];


(*eqns of motion for r,\[Theta],\[CurlyPhi] to be integrated*)
d\[CurlyPhi]dt = Simplify[ReplaceAll[u[[4]] - momentumEq\[CurlyPhi]Calc[r,\[Theta]]/momentumEqtCalc[r,\[Theta]], {Derivative[1][\[CurlyPhi]]->Derivative[1][\[CurlyPhi]][tt], r->r[tt], \[Theta]->\[Theta][tt]}]];
geortt = Simplify[ReplaceAll[acc[[2]]+u[[2]]*coefft+coeffr, {Derivative[2][r]->Derivative[2][r][tt], Derivative[1][r]->Derivative[1][r][tt], r->r[tt], Derivative[1][\[Theta]]->Derivative[1][\[Theta]][tt], \[Theta]->\[Theta][tt], Derivative[1][\[CurlyPhi]]-> Simplify[ReplaceAll[momentumEq\[CurlyPhi]Calc[r,\[Theta]]/momentumEqtCalc[r,\[Theta]], {r->r[tt], \[Theta]->\[Theta][tt]}]]}]];
geo\[Theta]tt = Simplify[ReplaceAll[acc[[3]]+u[[3]]*coefft+coeff\[Theta], {Derivative[2][\[Theta]]->Derivative[2][\[Theta]][tt], Derivative[1][\[Theta]]->Derivative[1][\[Theta]][tt], \[Theta]->\[Theta][tt], Derivative[1][r]->Derivative[1][r][tt], r->r[tt], Derivative[1][\[CurlyPhi]]-> Simplify[ReplaceAll[momentumEq\[CurlyPhi]Calc[r,\[Theta]]/momentumEqtCalc[r,\[Theta]], {r->r[tt], \[Theta]->\[Theta][tt]}]]}]];


(*4-velocity conservation*)
constraintEqCoord[r_, \[Theta]_, vrt_, v\[Theta]t_] = Simplify[ReplaceAll[metric . u . u-k, {Derivative[1][t]->momentumEqtCalc[r, \[Theta]], Derivative[1][r]->vrt*momentumEqtCalc[r, \[Theta]], Derivative[1][\[Theta]]->v\[Theta]t*momentumEqtCalc[r, \[Theta]], Derivative[1][\[CurlyPhi]]->momentumEq\[CurlyPhi]Calc[r, \[Theta]]}]];


(*carter constant*)
carterEqCoord[r_, \[Theta]_, v\[Theta]t_] = Simplify[g\[Theta]\[Theta][r,\[Theta]]^2*v\[Theta]t^2*momentumEqtCalc[r, \[Theta]]^2 + Cos[\[Theta]]^2*(a^2*M^2*(1 - (e/\[Mu])^2) + (Lz/\[Mu])^2 / Sin[\[Theta]]^2)];


(*fluxes (2PN) - from Gair 0510129 eqns 39-44-45-46*)
g1[ecc_] = 1+73/24*ecc^2+37/96*ecc^4;
g2[ecc_] = 73/12+823/24*ecc^2+949/32*ecc^4+491/192*ecc^6;
g3[ecc_] = 1247/336+9181/672*ecc^2;
g4[ecc_] = 4+1375/48*ecc^2;
g5[ecc_] = 44711/9072+172157/2592*ecc^2;
g6[ecc_] = 33/16+359/32*ecc^2;
g9[ecc_] = 1+7/8*ecc^2;
g10a[ecc_] = 61/24+63/8*ecc^2+95/64*ecc^4;
g10b[ecc_] = 61/8+91/4*ecc^2+461/64*ecc^4;
g11[ecc_] = 1247/336+425/336*ecc^2;
g12[ecc_] = 4+97/8*ecc^2;
g13[ecc_] = 44711/9072+302893/6048*ecc^2;
g14[ecc_] = 33/16+95/16*ecc^2;
eDot[pOrbit_, ecc_, iOrbit_] =  -32.0/5.0 * (\[Mu]/M)^2 * (M/pOrbit)^5 * (1-ecc^2)^(3/2) * (g1[ecc] - a*(M/pOrbit)^(3/2)*g2[ecc]*Cos[iOrbit] - (M/pOrbit)*g3[ecc] + Pi*(M/pOrbit)^(3/2)*g4[ecc] - (M/pOrbit)^2*g5[ecc] + a^2*(M/pOrbit)^2*g6[ecc] -527/96*a^2*(M/pOrbit)^2*Sin[iOrbit]^2);
LzDot[pOrbit_, ecc_, iOrbit_] = -32.0/5.0 * \[Mu]^2/M * (M/pOrbit)^(7/2) * (1-ecc^2)^(3/2) * (g9[ecc]*Cos[iOrbit] + a*(M/pOrbit)^(3/2)*(g10a[ecc] - Cos[iOrbit]^2*g10b[ecc]) - (M/pOrbit)*g11[ecc]*Cos[iOrbit] + Pi*(M/pOrbit)^(3/2)*g12[ecc]*Cos[iOrbit] - (M/pOrbit)^2*g13[ecc]*Cos[iOrbit] + a^2*(M/pOrbit)^2*Cos[iOrbit]*(g14[ecc] - 45/8*Sin[iOrbit]^2));


(*Q fluxes (2PN) - from Gair 0510129*)
QDot[pOrbit_, ecc_, iOrbit_] = -64.0/5.0 * \[Mu]^2/M * (M/pOrbit)^(7/2) * (1-ecc^2)^(3/2) * Sin[iOrbit] * Q0^(1/2) * (g9[ecc] - a*(M/pOrbit)^(3/2)*Cos[iOrbit]*g10b[ecc] - (M/pOrbit)*g11[ecc] + Pi*(M/pOrbit)^(3/2)*g12[ecc] - (M/pOrbit)^2*g13[ecc] + a^2*(M/pOrbit)^2*(g14[ecc] - 45/8*Sin[iOrbit]^2));
d1a = \[Minus]10.7420;
d1b = 28.5942;
d1c = \[Minus]9.07738;
d2a = \[Minus]1.42836;
d2b = 10.7003;
d2c = \[Minus]33.7090;
c1a = \[Minus]28.1517;
c1b = 60.9607;
c1c = 40.9998;
c2a = \[Minus]0.348161;
c2b = 2.37258;
c2c = \[Minus]66.6584;
c3a = \[Minus]0.715392;
c3b = 3.21593;
c3c = 5.28888;
c4a = \[Minus]7.61034;
c4b = 128.878;
c4c = \[Minus]475.465;
c5a = 12.2908;
c5b = \[Minus]113.125;
c5c = 306.119;
c6a = 40.9259;
c6b = \[Minus]347.271;
c6c = 886.503;
c7a = \[Minus]25.4831;
c7b = 224.227;
c7c = \[Minus]490.982;
c8a = \[Minus]9.00634;
c8b = 91.1767;
c8c = \[Minus]297.002;
c9a = \[Minus]0.645000;
c9b = \[Minus]5.13592;
c9c = 47.1982;
c10a = \[Minus]0.0309341;
c10b = \[Minus]22.2416;
c10c = 7.55265;
c11a = \[Minus]3.33476;
c11b = 22.7013;
c11c = \[Minus]12.4700;
f1a = \[Minus]283.955;
f1b = 736.209;
f2a = 483.266;
f2b = \[Minus]1325.19;
f3a = \[Minus]219.224;
f3b = 634.499;
f4a = \[Minus]25.8203;
f4b = 82.0780;
f5a = 301.478;
f5b = \[Minus]904.161;
f6a = \[Minus]271.966;
f6b = 827.319;
f7a = \[Minus]162.268;
f7b = 247.168;
f8a = 152.125;
f8b = \[Minus]182.165;
f9a = 184.465;
f9b = \[Minus]267.553;
f10a = \[Minus]188.132;
f10b = 254.067;
LzDotFit[pOrbit_, iOrbit_] = -32.0/5.0 * \[Mu]^2/M * (M/pOrbit)^(7/2) * (Cos[iOrbit] + a*(M/pOrbit)^(3/2)*(61/24 - 61/8*Cos[iOrbit]^2) - 1247/336*(M/pOrbit)*Cos[iOrbit] +4*Pi*(M/pOrbit)^(3/2)*Cos[iOrbit] - 44711/9072*(M/pOrbit)^(2)*Cos[iOrbit] + a^2*(M/pOrbit)^(2)*Cos[iOrbit]*(33/16 - 45/8*Sin[iOrbit]^2) + (M/pOrbit)^(5/2)*( a*(d1a + d1b*(M/pOrbit)^(1/2) + d1c*(M/pOrbit)) + a^3*(d2a + d2b*(M/pOrbit)^(1/2) + d2c*(M/pOrbit)) + Cos[iOrbit]*(c1a + c1b*(M/pOrbit)^(1/2) + c1c*(M/pOrbit)) + a^2*Cos[iOrbit]*(c2a + c2b*(M/pOrbit)^(1/2) + c2c*(M/pOrbit)) + a^4*Cos[iOrbit]*(c3a + c3b*(M/pOrbit)^(1/2) + c3c*(M/pOrbit)) + a*Cos[iOrbit]^2*(c4a + c4b*(M/pOrbit)^(1/2) + c4c*(M/pOrbit)) + a^3*Cos[iOrbit]^2*(c5a + c5b*(M/pOrbit)^(1/2) + c5c*(M/pOrbit)) + a^2*Cos[iOrbit]^3*(c6a + c6b*(M/pOrbit)^(1/2) + c6c*(M/pOrbit)) + a^4*Cos[iOrbit]^3*(c7a + c7b*(M/pOrbit)^(1/2) + c7c*(M/pOrbit)) + a^3*Cos[iOrbit]^4*(c8a + c8b*(M/pOrbit)^(1/2) + c8c*(M/pOrbit)) + a^4*Cos[iOrbit]^5*(c9a + c9b*(M/pOrbit)^(1/2) + c9c*(M/pOrbit))) + (M/pOrbit)^(7/2)*a*Cos[iOrbit]*(f1a + f1b*(M/pOrbit)^(1/2) + a*(f2a + f2b*(M/pOrbit)^(1/2)) + a^2*(f3a + f3b*(M/pOrbit)^(1/2)) + Cos[iOrbit]^2*(f4a + f4b*(M/pOrbit)^(1/2)) + a*Cos[iOrbit]^2*(f5a + f5b*(M/pOrbit)^(1/2)) + a^2*Cos[iOrbit]^2*(f6a + f6b*(M/pOrbit)^(1/2))));
iOrbitFit[pOrbit_, iOrbit_] = 32.0/5.0 * \[Mu]^2/M * a * Sin[iOrbit]^2 / Q0^(1/2) * (M/pOrbit)^(5) * (61/24 + (M/pOrbit)*(d1a + d1b*(M/pOrbit)^(1/2) + d1c*(M/pOrbit)) + a^2*(M/pOrbit)*(d2a + d2b*(M/pOrbit)^(1/2) + d2c*(M/pOrbit)) + a*Cos[iOrbit]*(M/pOrbit)^(1/2)*(c10a + c10b*(M/pOrbit) + c10c*(M/pOrbit)^(3/2)) + a^2*Cos[iOrbit]^2*(M/pOrbit)*(c11a + c11b*(M/pOrbit)^(1/2) + c11c*(M/pOrbit)) + a^3*Cos[iOrbit]*(M/pOrbit)^(5/2)*(f7a + f7b*(M/pOrbit)^(1/2) + a*(f8a + f8b*(M/pOrbit)^(1/2)) + Cos[iOrbit]^2*(f9a + f9b*(M/pOrbit)^(1/2)) + a*Cos[iOrbit]^2*(f10a + f10b*(M/pOrbit)^(1/2))));
QDotMod[pOrbit_, ecc_, iOrbit_] = (1-ecc^2)^(3/2) * Q0^(1/2) * ((1-ecc^2)^(-3/2) * QDot[pOrbit, ecc, iOrbit]/Q0^(1/2) - QDot[pOrbit, 0, iOrbit]/Q0^(1/2) + 2*Tan[iOrbit]*(LzDotFit[pOrbit, iOrbit] + Q0^(1/2)/(Sin[iOrbit]^2)*iOrbitFit[pOrbit, iOrbit]));
LzDotMod[pOrbit_, ecc_, iOrbit_] = (1-ecc^2)^(3/2) * ((1-ecc^2)^(-3/2) * LzDot[pOrbit, ecc, iOrbit] - LzDot[pOrbit, 0, iOrbit] + LzDotFit[pOrbit, iOrbit]);


(*eqns for orbital elements*)
rOrbitMinEq[pOrbit_,ecc_] = pOrbit / (1+ecc);
rOrbitMaxEq[pOrbit_,ecc_] = pOrbit / (1-ecc);


(*conversione to physical units*)
c = N[Entity["PhysicalConstant", "SpeedOfLight"]["Value"][[1]]]; (*m/s*)
G = N[Entity["PhysicalConstant", "GravitationalConstant"]["Value"][[1]]]; (*N m^2/Kg^2*)
DL = 2*10^8*3.08567758*10^16 ; (*luminosity distance = 200Mpc, z = 0.04487*)(*m*)
mSun = 1.9891 10^30; (*Kg*)
convMass = 1.0*10^6 * mSun;(*Kg*)
convLenght = convMass * G / c^2; (*m*)
convTime = convMass * G / c^3; (*s*)
day = 60*60*24;


(*TT gauge projection*)
n = {0, 0, 1};
deltaKron =  {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
p = Table[deltaKron[[j, k]] - n[[j]] * n[[k]], {j, 1, 3}, {k, 1, 3}];
pTT = Table[p[[j, m]] * p[[k, l]] - 1/2 * p[[j, k]] * p[[m, l]], {j, 1, 3}, {k, 1, 3}, {m, 1, 3}, {l, 1, 3}];


(*LISA antennae response patterns, some assumptions have been made to simplify the patterns (see Appendix in 2009.00028)*)
\[Theta]S = 0; 
\[Phi]S = Pi/2;
(*\[Theta]K = Pi/2; 
\[Phi]K = Pi/2;*)
T = 60*60*24*365;(*1yr*)
\[Theta]LISA[t_] = ArcCos[1/2*Cos[\[Theta]S] - Sqrt[3]/2*Sin[\[Theta]S]*Cos[2*Pi*t/T-\[Phi]S]];
\[CapitalPhi]LISA = Pi/2;
\[Phi]LISA[t_] = 2*Pi*t/T + \[CapitalPhi]LISA;
\[Psi]LISA[t_] = -2*Pi*t/T;
FIp[t_] = 1/2 * (1 + Cos[\[Theta]LISA[t]]^2) * Cos[2*\[Phi]LISA[t]] * Cos[2*\[Psi]LISA[t]] - Cos[\[Theta]LISA[t]] * Sin[2 * \[Phi]LISA[t]] * Sin[2*\[Psi]LISA[t]];
FIc[t_] = 1/2 * (1 + Cos[\[Theta]LISA[t]]^2) * Cos[2*\[Phi]LISA[t]] * Cos[2*\[Psi]LISA[t]] + Cos[\[Theta]LISA[t]] * Sin[2 * \[Phi]LISA[t]] * Sin[2*\[Psi]LISA[t]];
FIIp[t_] = 1/2 * (1 + Cos[\[Theta]LISA[t]]^2) * Sin[2*\[Phi]LISA[t]] * Cos[2*\[Psi]LISA[t]] + Cos[\[Theta]LISA[t]] * Cos[2 * \[Phi]LISA[t]] * Sin[2*\[Psi]LISA[t]];
FIIc[t_] = 1/2 * (1 + Cos[\[Theta]LISA[t]]^2) * Sin[2*\[Phi]LISA[t]] * Sin[2*\[Psi]LISA[t]] - Cos[\[Theta]LISA[t]] * Cos[2 * \[Phi]LISA[t]] * Cos[2*\[Psi]LISA[t]];


(*coefficients and elliptic integrals from Schmidt 0202090, e\[NotEqual]0*)
ff[r_] = r^4 + a^2 (r(r+2) + zMinus^2*(\[CapitalDelta][r]/M^2));
gg[r_] = 2*a*r;
hh[r_] = r(r-2) + zMinus^2*(\[CapitalDelta][r]/M^2)/(1 - zMinus^2);
dd[r_] = (r^2 + a^2*zMinus^2)*(\[CapitalDelta][r]/M^2);
KK[kk_] := NIntegrate[1/(1 - kk^2*Sin[thetaAus]^2)^(1/2), {thetaAus, 0, Pi/2}]; 
EE[kk_] := NIntegrate[(1 - kk^2*Sin[thetaAus]^2)^(1/2), {thetaAus, 0, Pi/2}];
\[CapitalPi][hhh_, kk_] := NIntegrate[1 / ((1 + hhh*Sin[thetaAus]^2)*(1 - kk^2*Sin[thetaAus]^2)^(1/2)), {thetaAus, 0, Pi/2}];


(* ***Effective Resonance Model from Speri, Gair 2103.06306*** *)
Csi[omegaTheta_, omegaRadial_] = Abs[1.5 - omegaTheta/omegaRadial]; (*trigger function*)
normFactor = NIntegrate[1 + Cos[4*Pi*xAus^2], {xAus, -0.5, 0.5}];
\[Omega][tt_] = (1 + Cos[4*Pi*((tt - tRes0) / tRes)^2]) / normFactor; (*impulse function*)


(*initial parameters*)
k = -1; 
M = 1.0;
\[Mu] = 1.0*10^-5;
Mtot = (M+\[Mu])*convMass; (*Kg*)
Mreduced = (M*\[Mu]*convMass^2)/Mtot;(*reduced mass [Kg]*)
chirpM = Mreduced^(3/5) * Mtot^(2/5); (*Kg*)
(*initPar = Import["/THEORY/USERS/edoardo.levati/DATA/kick/kerr/parameter_space_final.csv"]; (* load initial orbital parameters (a, i, ecc, p, tRes) *)
nIter = Dimensions[initPar][[1]];*)


(*inspiral*)
Needs["Developer`"];
LaunchKernels[5];
For[ii=1, ii<=(*IntegerPart[nIter]*) 1, ii++, 
{
	Print["-> EMRI = ", 1, " / ", 1];
	
	(* ***initial orbital parameters (a, i, ecc, p, tResStart, tRes) *** *)
	a = 0.95; 
	iOrbitInit = 1.3; 
	eccInit = 0.1; 
	pOrbitInit = 8.905099999992064 M; 
	tResStart = 96200 M; 
	tRes = 5057.87145012503 M; 
	
	(* ***E/Lz from orbital elements p, e, i*** *)
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
	
	(* ***initial conditions for NDSolve*** *)
	r0tt = rPeriapsis; (*rApoapsis;*)
	\[Theta]0tt = Pi/2; 
	\[CurlyPhi]0tt = 0.0;
	vr0tt = 0.0; 
	v\[Theta]0tt = (v\[Theta]t /. Solve[constraintEqCoord[r0tt, \[Theta]0tt, vr0tt, v\[Theta]t] == 0, v\[Theta]t, Reals][[2]]);
	(*v\[Theta]0tt = Sqrt[-1 / g\[Theta]\[Theta][r0tt,\[Theta]0tt] * (1 - ((e/\[Mu])^2 * g\[CurlyPhi]\[CurlyPhi][r0tt,\[Theta]0tt] + (Lz/\[Mu])^2 * gtt[r0tt,\[Theta]0tt] + 2*gt\[CurlyPhi][r0tt,\[Theta]0tt]*(e/\[Mu])*(Lz/\[Mu])) / (gt\[CurlyPhi][r0tt,\[Theta]0tt]^2 - gtt[r0tt,\[Theta]0tt]*g\[CurlyPhi]\[CurlyPhi][r0tt,\[Theta]0tt]))] * (gt\[CurlyPhi][r0tt,\[Theta]0tt]^2 - gtt[r0tt,\[Theta]0tt]*g\[CurlyPhi]\[CurlyPhi][r0tt,\[Theta]0tt]) / ((e/\[Mu]) * g\[CurlyPhi]\[CurlyPhi][r0tt,\[Theta]0tt] + (Lz/\[Mu]) * gt\[CurlyPhi][r0tt,\[Theta]0tt])*)
	tCoalesc = 5/256 * c^5/G^3 * (r0tt*convLenght)^4 / (Mtot^2*Mreduced); (*coalescence time quadrupole formula [s]*) 
	(*LuminosityGW = (32/5 * G^4 * Mreduced^2 * Mtot^3) / (c^5*(r0tt*convLenght)^5);  (*Luminosity quadrupole formula [W]*) *)
	tMax = 3.2*10^6 M;
	nDiscr = 10^6; (*on Vera -> Nmax = 10^9*)
	dttArray1 = 100 M; (*approach resonance*)
	dttArray2 = 200 M; (*in resonance*)
	dttArray3 = 200 M; (*from tResEnd to tMax*)
	dttSolve = 10^3 M;	
	dttDiscr = dttSolve / nDiscr; (*discretizing the interpolating functions coming from NDSolve*) 
	dttDiscrGW1 = 0.1 M;
	dttDiscrGW2 = 0.2 M;
	dttDiscrGW3 = 0.2 M;
	h0 = 4/DL * (M*convLenght) * (\[Mu]*convLenght) * 1/(r0tt*convLenght); (*GW strain from quadrupole approximation*)
	
	(* ***print initial conditions and parameters*** *)
	initHeader = {"M = ", "\[Mu] = ", "a = ", "p = ", "ecc = ", "i = ", "e / \[Mu] = ", "Lz / \[Mu] = ", "\[Tau]res (hours) = ", "t_max (days) = ", "t_coalescence (days) = ", "DL (Mpc) = ", "chirp mass (Msun) = ", "h0 = "};
	initData = {M, \[Mu], Round[a, 0.001], Round[pOrbitInit, 0.001], Round[eccInit, 0.001], Round[iOrbitInit, 0.001], Round[e/\[Mu], 0.001], Round[Lz/(\[Mu]), 0.001], Round[tRes*convTime/(60*60), 0.01], Round[tMax*convTime/day], Round[tCoalesc/day], Round[DL/(3.08567758*10^16*10^6)], Round[chirpM/mSun], Round[h0, N[10^-22]]}; 
	Print["\t", "-> Initial conditions"];
	Do[Print["\t\t", initHeader[[iii]], initData[[iii]]], {iii, Length[initData]}];
	
	(* ***dissipation loop t < tResStart*** *)
	Print["\t","-> Dissipation loop | t < tResStart"];
	tt0 = 0.0;
	ttmax = tResStart; 
	ttArray = Range[tt0, ttmax, dttArray1];
	nUpdates = IntegerPart[Round[(ttmax-tt0)/dttArray1]]; (*number of iterations used to update integrals of motion (E, Lz)*)
	For[j=1, j<=nUpdates, j++, 
	{
		Print["\t\t","update = ", j, " / ", nUpdates];
	
		(* ***set initial values for E, Lz, Q*** *)
		e0 = e;
		Lz0 = Lz;
		Q0 = Q;
		
		If[j>1,
		{
			(* ***initial positions/velocities*** *)
			r0tt = rtt[ttArray[[j]]]; 
			\[Theta]0tt = \[Theta]tt[ttArray[[j]]]; 
			\[CurlyPhi]0tt = \[CurlyPhi]tt[ttArray[[j]]]; 
			vr0ttTemp = vrtt[ttArray[[j]]]; 
			v\[Theta]0ttTemp = v\[Theta]tt[ttArray[[j]]];
			eq1 = Simplify[carterEqCoord[r0tt, \[Theta]0tt, v\[Theta]0Sol] - Q0/\[Mu]^2];
			If[Sign[v\[Theta]0ttTemp] < 0, v\[Theta]SignCond = (v\[Theta]0Sol<0), v\[Theta]SignCond = (v\[Theta]0Sol>0)]; 
			v\[Theta]0tt = Solve[eq1==0 && v\[Theta]SignCond, {v\[Theta]0Sol}, Reals][[1, 1, 2]];
			eq2 = Simplify[constraintEqCoord[r0tt, \[Theta]0tt, vr0Sol, v\[Theta]0tt]];
			If[Sign[vr0ttTemp] < 0, vrSignCond = (vr0Sol<0), vrSignCond = (vr0Sol>0)]; 
			solInit = Solve[eq2==0 && vrSignCond, {vr0Sol}, Reals];
			If[Length[solInit]!=0, vr0tt = solInit[[1, 1, 2]], 
			{
				(* ***if not possible to satisfy constraintEq == 0 introduce arbitrary error*** *) 
				errPower = 20;
				errConstr = 1.0*10^-(errPower);
				jj = 1;
				count = 0;
				While[Length[Solve[eq2==errConstr && vrSignCond, {vr0Sol}, Reals]]==0,
				{
					jj += 1;
					If[Mod[jj, 10] == 0, {jj=1; count+=1; errConstr += 10^-(errPower-(count-1));}, errConstr += 10^-(errPower-count);];
				}];
				vr0tt = Solve[eq2==errConstr && vrSignCond, {vr0Sol}, Reals][[1, 1, 2]];
			}];
		}];
		
		(* ***integration*** *)
		Print["\t\t\t carter err = ", carterEqCoord[r0tt, \[Theta]0tt, v\[Theta]0tt] - Q0/\[Mu]^2];
		Print["\t\t\t 4-velocity err = ", constraintEqCoord[r0tt, \[Theta]0tt, vr0tt, v\[Theta]0tt]];
		invariant = {carterEqCoord[r[tt], \[Theta][tt], Derivative[1][\[Theta]][tt]], constraintEqCoord[r[tt], \[Theta][tt], Derivative[1][r][tt], Derivative[1][\[Theta]][tt]]};
		ProjImpRK = {"Projection", Method->"ImplicitRungeKutta", "Invariants"->invariant};
		eqnstt = {geortt==0, geo\[Theta]tt==0, d\[CurlyPhi]dt==0, r[ttArray[[j]]]==r0tt, Derivative[1][r][ttArray[[j]]]==vr0tt, \[Theta][ttArray[[j]]]==\[Theta]0tt, Derivative[1][\[Theta]][ttArray[[j]]]==v\[Theta]0tt, \[CurlyPhi][ttArray[[j]]]==\[CurlyPhi]0tt};
		varstt = {r, \[Theta], \[CurlyPhi], Derivative[1][r], Derivative[1][\[Theta]]};
		orbitTimett = NDSolve[eqnstt, varstt, {tt, ttArray[[j]], ttArray[[j]] + dttSolve}, MaxSteps->\[Infinity], PrecisionGoal->MachinePrecision, Method->ProjImpRK];
		rtt = r/.orbitTimett[[1, 1]];
		\[Theta]tt = \[Theta]/.orbitTimett[[1, 2]];
		\[CurlyPhi]tt = \[CurlyPhi]/.orbitTimett[[1, 3]];
		vrtt = Derivative[1][r]/.orbitTimett[[1, 4]];
		v\[Theta]tt = Derivative[1][\[Theta]]/.orbitTimett[[1, 5]];
		
		(* ***discretizing NDSolve solution*** *)
		f = varstt/.orbitTimett[[1]];
		data = Transpose[ToPackedArray[Parallelize[Table[f[[kindex]][Range[ttArray[[j]], ttArray[[j]] + dttSolve, dttDiscr]],{kindex, 1, 5, 1}], ProgressReporting->False]]];
		
		(* ***orbital elements from discretized data*** *)
		rOrbitMin = Min[data[[1;;-1,1]]];
		rOrbitMax = Max[data[[1;;-1,1]]];
		pOrbit = (2*rOrbitMin*rOrbitMax) / (rOrbitMax + rOrbitMin);
		ecc = (rOrbitMax - rOrbitMin) / (rOrbitMax + rOrbitMin);
		\[Theta]OrbitMin = Min[data[[1;;-1,2]]];
		\[Theta]OrbitMax = Max[data[[1;;-1,2]]];
		iOrbit = (\[Theta]OrbitMax - \[Theta]OrbitMin)/2.;
		
		(* ***fundamental frequencies in Kerr - from Fujita, Hikida 0906.1420*** *)
		zMinus = Cos[\[Theta]OrbitMin];
		r1 = rOrbitMax; (*turning point*)
		r2 = rOrbitMin; (*turning point*)
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
		csiVal = Csi[omegaThetaCoord, omegaRadialCoord];
		
		(* ***save energy and orbital frequencies data*** *)
		tData[j] = ttArray[[j]]; 
		energyData[j] = e; 
		angularData[j] = Lz; 
		carterData[j] = Q; 
		omegaRadialData[j] = omegaRadialCoord; 
		omegaThetaData[j] = omegaThetaCoord; 
		csiData[j] = csiVal;
		pData[j] = pOrbit;
		eccData[j] = ecc;
		iData[j] = iOrbit;
		
		(* ***fluxes*** *)
		eFlux = eDot[pOrbit, ecc, iOrbit] * dttArray1;
		LzFlux = LzDot[pOrbit, ecc, iOrbit] * dttArray1;
		QFlux = QDotMod[pOrbit, ecc, iOrbit] * dttArray1;
		e = e0+eFlux; 
		Lz = Lz0+LzFlux;
		Q = Q0+QFlux;
	}];
	
	(* ***dissipation loop tResStart <= t <= tResEnd*** *)
	Print["\t","-> Dissipation loop | tResStart <= t <= tResEnd"];
	nUpdatesTot = nUpdates;
	tRes0 = tResStart + tRes/2;
	tResEnd = tResStart + tRes;
	tt0 = ttArray[[-1]];
	ttmax = Round[tResEnd, 200]; 
	ttArray = Range[tt0, ttmax, dttArray2];
	nUpdates = IntegerPart[Round[(ttmax-tt0)/dttArray2]]; (*number of iterations used to update integrals of motion (E, Lz)*)
	For[j=1, j<=nUpdates, j++, 
	{
		Print["\t\t","update = ", j, " / ", nUpdates];
		
		(* ***set initial values for E, Lz*** *)
		e0 = e;
		Lz0 = Lz;
		Q0 = Q;
		
		(* ***initial positions/velocities*** *)
		r0tt = rtt[ttArray[[j]]]; 
		\[Theta]0tt = \[Theta]tt[ttArray[[j]]]; 
		\[CurlyPhi]0tt = \[CurlyPhi]tt[ttArray[[j]]]; 
		vr0ttTemp = vrtt[ttArray[[j]]]; 
		v\[Theta]0ttTemp = v\[Theta]tt[ttArray[[j]]];
		eq1 = Simplify[carterEqCoord[r0tt, \[Theta]0tt, v\[Theta]0Sol] - Q0/\[Mu]^2];
		If[Sign[v\[Theta]0ttTemp] < 0, v\[Theta]SignCond = (v\[Theta]0Sol<0), v\[Theta]SignCond = (v\[Theta]0Sol>0)]; 
		v\[Theta]0tt = Solve[eq1==0 && v\[Theta]SignCond, {v\[Theta]0Sol}, Reals][[1, 1, 2]];
		eq2 = Simplify[constraintEqCoord[r0tt, \[Theta]0tt, vr0Sol, v\[Theta]0tt]];
		If[Sign[vr0ttTemp] < 0, vrSignCond = (vr0Sol<0), vrSignCond = (vr0Sol>0)]; 
		solInit = Solve[eq2==0 && vrSignCond, {vr0Sol}, Reals];
		If[Length[solInit]!=0, vr0tt = solInit[[1, 1, 2]], 
		{
			(* ***if not possible to satisfy constraintEq == 0 introduce arbitrary error*** *) 
			errPower = 20;
			errConstr = 1.0*10^-(errPower);
			jj = 1;
			count = 0;
			While[Length[Solve[eq2==errConstr && vrSignCond, {vr0Sol}, Reals]]==0,
			{
				jj += 1;
				If[Mod[jj, 10] == 0, {jj=1; count+=1; errConstr += 10^-(errPower-(count-1));}, errConstr += 10^-(errPower-count);];
			}];
			vr0tt = Solve[eq2==errConstr && vrSignCond, {vr0Sol}, Reals][[1, 1, 2]];
		}];
		
		(* ***integration*** *)
		Print["\t\t\t carter err = ", carterEqCoord[r0tt, \[Theta]0tt, v\[Theta]0tt] - Q0/\[Mu]^2];
		Print["\t\t\t 4-velocity err = ", constraintEqCoord[r0tt, \[Theta]0tt, vr0tt, v\[Theta]0tt]];
		invariant = {carterEqCoord[r[tt], \[Theta][tt], Derivative[1][\[Theta]][tt]], constraintEqCoord[r[tt], \[Theta][tt], Derivative[1][r][tt], Derivative[1][\[Theta]][tt]]};
		ProjImpRK = {"Projection", Method->"ImplicitRungeKutta", "Invariants"->invariant};
		eqnstt = {geortt==0, geo\[Theta]tt==0, d\[CurlyPhi]dt==0, r[ttArray[[j]]]==r0tt, Derivative[1][r][ttArray[[j]]]==vr0tt, \[Theta][ttArray[[j]]]==\[Theta]0tt, Derivative[1][\[Theta]][ttArray[[j]]]==v\[Theta]0tt, \[CurlyPhi][ttArray[[j]]]==\[CurlyPhi]0tt};
		varstt = {r, \[Theta], \[CurlyPhi], Derivative[1][r], Derivative[1][\[Theta]]};
		orbitTimett = NDSolve[eqnstt, varstt, {tt, ttArray[[j]], ttArray[[j]] + dttSolve}, MaxSteps->\[Infinity], PrecisionGoal->MachinePrecision, Method->ProjImpRK];
		rtt = r/.orbitTimett[[1, 1]];
		\[Theta]tt = \[Theta]/.orbitTimett[[1, 2]];
		\[CurlyPhi]tt = \[CurlyPhi]/.orbitTimett[[1, 3]];
		vrtt = Derivative[1][r]/.orbitTimett[[1, 4]];
		v\[Theta]tt = Derivative[1][\[Theta]]/.orbitTimett[[1, 5]];
		
		(* ***discretizing NDSolve solution*** *)
		f = varstt/.orbitTimett[[1]];
		data = Transpose[ToPackedArray[Parallelize[Table[f[[kindex]][Range[ttArray[[j]], ttArray[[j]] + dttSolve, dttDiscr]],{kindex, 1, 5, 1}], ProgressReporting->False]]];
		
		(* ***orbital elements from discretized data*** *)
		rOrbitMin = Min[data[[1;;-1,1]]];
		rOrbitMax = Max[data[[1;;-1,1]]];
		pOrbit = (2*rOrbitMin*rOrbitMax) / (rOrbitMax + rOrbitMin);
		ecc = (rOrbitMax - rOrbitMin) / (rOrbitMax + rOrbitMin);
		\[Theta]OrbitMin = Min[data[[1;;-1,2]]];
		\[Theta]OrbitMax = Max[data[[1;;-1,2]]];
		iOrbit = (\[Theta]OrbitMax - \[Theta]OrbitMin)/2.;
		
		(* ***fundamental frequencies in Kerr - from Fujita, Hikida 0906.1420*** *)
		zMinus = Cos[\[Theta]OrbitMin];
		r1 = rOrbitMax; (*turning point*)
		r2 = rOrbitMin; (*turning point*)
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
		csiVal = Csi[omegaThetaCoord, omegaRadialCoord];
		
		(* ***save energy and frequencies data*** *)
		tData[j+nUpdatesTot] = ttArray[[j]]; 
		energyData[j+nUpdatesTot] = e; 
		angularData[j+nUpdatesTot] = Lz; 
		carterData[j+nUpdatesTot] = Q;
		omegaRadialData[j+nUpdatesTot] = omegaRadialCoord; 
		omegaThetaData[j+nUpdatesTot] = omegaThetaCoord; 
		csiData[j+nUpdatesTot] = csiVal;
		pData[j+nUpdatesTot] = pOrbit;
		eccData[j+nUpdatesTot] = ecc;
		iData[j+nUpdatesTot] = iOrbit;
		
		(* ***fluxes*** *)
		eFlux = (eDot[pOrbit, ecc, iOrbit]) * dttArray2;
		LzFlux = (LzDot[pOrbit, ecc, iOrbit]) * dttArray2;
		QFlux = QDotMod[pOrbit, ecc, iOrbit] * dttArray2;
		e = e0+eFlux; 
		Lz = Lz0+LzFlux;
		Q = Q0+QFlux;
		
		(* ***GWs*** *)
		x[t_] = rtt[t] * Sin[\[Theta]tt[t]] * Cos[\[CurlyPhi]tt[t]] * convLenght;
		y[t_] = rtt[t] * Sin[\[Theta]tt[t]] * Sin[\[CurlyPhi]tt[t]] * convLenght;
		z[t_] = rtt[t] * Cos[\[Theta]tt[t]] * convLenght;
		q[t_] = {{x[t]^2, x[t]*y[t], x[t]*z[t]}, {x[t]*y[t], y[t]^2, y[t]*z[t]}, {x[t]*z[t], y[t]*z[t], z[t]^2}};
		qReduced[t_] = Table[q[t][[j, k]] - 1/3*deltaKron[[j, k]]*Tr[q[t]], {j, 1, 3}, {k, 1, 3}];
		(*QTT[t_] = Table[Sum[pTT[[j, k, m, l]] * q[t][[m, l]], {m, 1, 3}, {l, 1, 3}], {j, 1, 3}, {k, 1, 3}];*)
		QTT[t_] = Table[Sum[pTT[[j, k, m, l]] * qReduced[t][[m, l]], {m, 1, 3}, {l, 1, 3}], {j, 1, 3}, {k, 1, 3}];
		hTT[t_] = 2*G*(\[Mu]*convMass)/(c^4 DL) * Derivative[2][QTT][t];
		hTTplus[t_] = hTT[t][[1]][[1]];
		hTTcross[t_] = hTT[t][[1]][[2]];
		hI[t_] = (\[Sqrt]3/2)*(FIp[t]*hTTplus[t] + FIc[t]*hTTcross[t]); 
		hII[t_] = (\[Sqrt]3/2)*(FIIp[t]*hTTplus[t] + FIIc[t]*hTTcross[t]); 
		hdata = ToPackedArray[Transpose[{Range[ttArray[[j]], ttArray[[j+1]], dttDiscrGW2], hI[Range[ttArray[[j]], ttArray[[j+1]], dttDiscrGW2]], hII[Range[ttArray[[j]], ttArray[[j+1]], dttDiscrGW2]]}]];
		
		(* ***export GWs data*** *)
		fname = StringJoin["/THEORY/USERS/edoardo.levati/DATA/kick/kerr/hdata_test/Method1_", ToString[ii], "/hdata_", ToString[j+nUpdatesTot], ".h5"];
		Export[fname, hdata, OverwriteTarget -> True];
	}];
	
	(* ***dissipation loop t > tResEnd*** *)
	Print["\t","-> Dissipation loop | t > tResEnd"];
	nUpdatesTot = nUpdatesTot + nUpdates;
	tt0 = ttArray[[-1]];
	nUpdates = IntegerPart[Round[(tMax-tt0)/dttArray3]]; (*number of iterations used to update integrals of motion (E, Lz)*)
	ttmax = tt0 + nUpdates*dttArray3;
	ttArray = Range[tt0, ttmax, dttArray3];
	For[j=1, j<=nUpdates, j++, 
	{
		Print["\t\t","update = ", j, " / ", nUpdates];
	
		(* ***set initial values for E, Lz*** *)
		e0 = e;
		Lz0 = Lz;
		Q0 = Q;
		
		(* ***initial positions/velocities*** *)
		r0tt = rtt[ttArray[[j]]]; 
		\[Theta]0tt = \[Theta]tt[ttArray[[j]]]; 
		\[CurlyPhi]0tt = \[CurlyPhi]tt[ttArray[[j]]]; 
		vr0ttTemp = vrtt[ttArray[[j]]]; 
		v\[Theta]0ttTemp = v\[Theta]tt[ttArray[[j]]];
		eq1 = Simplify[carterEqCoord[r0tt, \[Theta]0tt, v\[Theta]0Sol] - Q0/\[Mu]^2];
		If[Sign[v\[Theta]0ttTemp] < 0, v\[Theta]SignCond = (v\[Theta]0Sol<0), v\[Theta]SignCond = (v\[Theta]0Sol>0)]; 
		v\[Theta]0tt = Solve[eq1==0 && v\[Theta]SignCond, {v\[Theta]0Sol}, Reals][[1, 1, 2]];
		eq2 = Simplify[constraintEqCoord[r0tt, \[Theta]0tt, vr0Sol, v\[Theta]0tt]];
		If[Sign[vr0ttTemp] < 0, vrSignCond = (vr0Sol<0), vrSignCond = (vr0Sol>0)]; 
		solInit = Solve[eq2==0 && vrSignCond, {vr0Sol}, Reals];
		If[Length[solInit]!=0, vr0tt = solInit[[1, 1, 2]], 
		{
			(* ***if not possible to satisfy constraintEq == 0 introduce arbitrary error*** *) 
			errPower = 20;
			errConstr = 1.0*10^-(errPower);
			jj = 1;
			count = 0;
			While[Length[Solve[eq2==errConstr && vrSignCond, {vr0Sol}, Reals]]==0,
			{
				jj += 1;
				If[Mod[jj, 10] == 0, {jj=1; count+=1; errConstr += 10^-(errPower-(count-1));}, errConstr += 10^-(errPower-count);];
			}];
			vr0tt = Solve[eq2==errConstr && vrSignCond, {vr0Sol}, Reals][[1, 1, 2]];
		}];

		(* ***integration*** *)
		Print["\t\t\t carter err = ", carterEqCoord[r0tt, \[Theta]0tt, v\[Theta]0tt] - Q0/\[Mu]^2];
		Print["\t\t\t 4-velocity err = ", constraintEqCoord[r0tt, \[Theta]0tt, vr0tt, v\[Theta]0tt]];
		invariant = {carterEqCoord[r[tt], \[Theta][tt], Derivative[1][\[Theta]][tt]], constraintEqCoord[r[tt], \[Theta][tt], Derivative[1][r][tt], Derivative[1][\[Theta]][tt]]};
		ProjImpRK = {"Projection", Method->"ImplicitRungeKutta", "Invariants"->invariant};
		eqnstt = {geortt==0, geo\[Theta]tt==0, d\[CurlyPhi]dt==0, r[ttArray[[j]]]==r0tt, Derivative[1][r][ttArray[[j]]]==vr0tt, \[Theta][ttArray[[j]]]==\[Theta]0tt, Derivative[1][\[Theta]][ttArray[[j]]]==v\[Theta]0tt, \[CurlyPhi][ttArray[[j]]]==\[CurlyPhi]0tt};
		varstt = {r, \[Theta], \[CurlyPhi], Derivative[1][r], Derivative[1][\[Theta]]};
		orbitTimett = Quiet[Check[NDSolve[eqnstt, varstt, {tt, ttArray[[j]], ttArray[[j]] + dttSolve}, MaxSteps->\[Infinity], PrecisionGoal->MachinePrecision, Method->ProjImpRK], err, {NDSolve::impsitsf, NIntegrate::slwcon, NIntegrate::ncvb, General::stop, InterpolatingFunction::dmval}]];
		If[orbitTimett == err,
		{
			Print["Plunge -> exiting dissipation loop"];
			Break[]; (* exit loops if system plunges *)
		}];
		rtt = r/.orbitTimett[[1, 1]];
		\[Theta]tt = \[Theta]/.orbitTimett[[1, 2]];
		\[CurlyPhi]tt = \[CurlyPhi]/.orbitTimett[[1, 3]];
		vrtt = Derivative[1][r]/.orbitTimett[[1, 4]];
		v\[Theta]tt = Derivative[1][\[Theta]]/.orbitTimett[[1, 5]];
		
		(* ***discretizing NDSolve solution*** *)
		f = varstt/.orbitTimett[[1]];
		data = Transpose[ToPackedArray[Parallelize[Table[f[[kindex]][Range[ttArray[[j]], ttArray[[j]] + dttSolve, dttDiscr]],{kindex, 1, 5, 1}], ProgressReporting->False]]];
		
		(* ***orbital elements from discretized data*** *)
		rOrbitMin = Min[data[[1;;-1,1]]];
		rOrbitMax = Max[data[[1;;-1,1]]];
		pOrbit = (2*rOrbitMin*rOrbitMax) / (rOrbitMax + rOrbitMin);
		ecc = (rOrbitMax - rOrbitMin) / (rOrbitMax + rOrbitMin);
		\[Theta]OrbitMin = Min[data[[1;;-1,2]]];
		\[Theta]OrbitMax = Max[data[[1;;-1,2]]];
		iOrbit = (\[Theta]OrbitMax - \[Theta]OrbitMin)/2.;
		
		(* ***fundamental frequencies in Kerr - from Fujita, Hikida 0906.1420*** *)
		zMinus = Cos[\[Theta]OrbitMin];
		r1 = rOrbitMax; (*turning point*)
		r2 = rOrbitMin; (*turning point*)
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
		csiVal = Csi[omegaThetaCoord, omegaRadialCoord];
		
		(* ***save energy and frequencies data*** *)
		tData[j+nUpdatesTot] = ttArray[[j]]; 
		energyData[j+nUpdatesTot] = e;
		angularData[j+nUpdatesTot] = Lz;
		carterData[j+nUpdatesTot] = Q;
		omegaRadialData[j+nUpdatesTot] = omegaRadialCoord; 
		omegaThetaData[j+nUpdatesTot] = omegaThetaCoord; 
		csiData[j+nUpdatesTot] = csiVal;
		pData[j+nUpdatesTot] = pOrbit;
		eccData[j+nUpdatesTot] = ecc;
		iData[j+nUpdatesTot] = iOrbit;
		
		(* ***fluxes*** *)
		eFlux = eDot[pOrbit, ecc, iOrbit] * dttArray3;
		LzFlux = LzDot[pOrbit, ecc, iOrbit] * dttArray3;
		QFlux = QDotMod[pOrbit, ecc, iOrbit] * dttArray3;
		e = e0+eFlux; 
		Lz = Lz0+LzFlux;
		Q = Q0+QFlux;
		
		(* ***GWs*** *)
		x[t_] = rtt[t] * Sin[\[Theta]tt[t]] * Cos[\[CurlyPhi]tt[t]] * convLenght;
		y[t_] = rtt[t] * Sin[\[Theta]tt[t]] * Sin[\[CurlyPhi]tt[t]] * convLenght;
		z[t_] = rtt[t] * Cos[\[Theta]tt[t]] * convLenght;
		q[t_] = {{x[t]^2, x[t]*y[t], x[t]*z[t]}, {x[t]*y[t], y[t]^2, y[t]*z[t]}, {x[t]*z[t], y[t]*z[t], z[t]^2}};
		qReduced[t_] = Table[q[t][[j, k]] - 1/3*deltaKron[[j, k]]*Tr[q[t]], {j, 1, 3}, {k, 1, 3}];
		(*QTT[t_] = Table[Sum[pTT[[j, k, m, l]] * q[t][[m, l]], {m, 1, 3}, {l, 1, 3}], {j, 1, 3}, {k, 1, 3}];*)
		QTT[t_] = Table[Sum[pTT[[j, k, m, l]] * qReduced[t][[m, l]], {m, 1, 3}, {l, 1, 3}], {j, 1, 3}, {k, 1, 3}];
		hTT[t_] = 2*G*(\[Mu]*convMass)/(c^4 DL) * Derivative[2][QTT][t];
		hTTplus[t_] = hTT[t][[1]][[1]];
		hTTcross[t_] = hTT[t][[1]][[2]];
		hI[t_] = (\[Sqrt]3/2)*(FIp[t]*hTTplus[t] + FIc[t]*hTTcross[t]); 
		hII[t_] = (\[Sqrt]3/2)*(FIIp[t]*hTTplus[t] + FIIc[t]*hTTcross[t]); 
		hdata = ToPackedArray[Transpose[{Range[ttArray[[j]], ttArray[[j+1]], dttDiscrGW3], hI[Range[ttArray[[j]], ttArray[[j+1]], dttDiscrGW3]], hII[Range[ttArray[[j]], ttArray[[j+1]], dttDiscrGW3]]}]];
		
		(* ***export GWs data*** *)
		fname = StringJoin["/THEORY/USERS/edoardo.levati/DATA/kick/kerr/hdata_test/Method1_", ToString[ii], "/hdata_", ToString[j+nUpdatesTot], ".h5"];
		Export[fname, hdata, OverwriteTarget -> True];
	}];
	
	(*(* ***export energy data*** *)
	nUpdatesTot = nUpdatesTot + (j-1);
	data2 = Table[{tData[i], energyData[i], angularData[i], carterData[i], omegaRadialData[i], omegaThetaData[i], csiData[i], pData[i], eccData[i], iData[i]}, {i, 1, nUpdatesTot}];
	fname2 = StringJoin["/THEORY/USERS/edoardo.levati/DATA/kick/kerr/energyData_adiabatic_TEST_EMRI2_Coupled_Method_1/energy_", ToString[1], ".h5"];
	Export[fname2, data2, OverwriteTarget -> True];*)
}];
CloseKernels[];
Exit[];
