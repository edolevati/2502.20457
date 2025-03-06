#!/usr/bin/env wolframscript
(* ::Package:: *)

(*Keplerian parametrization to simplify the treatment of turning points in the r and \[Theta] motion*)
ClearAll["Global`*"]
r[\[Psi]_] = pOrbit0 / (1 + ecc0*Cos[\[Psi]]);
cos\[Theta][\[Chi]_] = zMinus * Cos[\[Chi]];


(*Kerr metric*)
rPlus = M + (M^2 - a^2*M^2)^(1/2);
rMinus = M - (M^2 - a^2*M^2)^(1/2);
\[CapitalSigma][\[Psi]_,\[Chi]_] = r[\[Psi]]^2 + a^2*M^2*cos\[Theta][\[Chi]]^2; (*a = dimensionless BH spin*)
\[CapitalDelta][\[Psi]_] = r[\[Psi]]^2 + a^2*M^2 - 2*M*r[\[Psi]];
\[CapitalDelta]\[CapitalDelta][rr_] = rr^2 + a^2*M^2 - 2*M*rr;
gtt[\[Psi]_,\[Chi]_] = -(1-(2*M*r[\[Psi]])/\[CapitalSigma][\[Psi],\[Chi]]);
grr[\[Psi]_,\[Chi]_] = \[CapitalSigma][\[Psi],\[Chi]]/\[CapitalDelta][\[Psi]];
g\[Theta]\[Theta][\[Psi]_,\[Chi]_] = \[CapitalSigma][\[Psi],\[Chi]];
g\[CurlyPhi]\[CurlyPhi][\[Psi]_,\[Chi]_] = (r[\[Psi]]^2 + a^2*M^2 + (2*M^3*r[\[Psi]]*a^2*(1 - cos\[Theta][\[Chi]]^2))/\[CapitalSigma][\[Psi],\[Chi]])*(1 - cos\[Theta][\[Chi]]^2);
gt\[CurlyPhi][\[Psi]_,\[Chi]_] = -((2*M^2*r[\[Psi]]*a*(1 - cos\[Theta][\[Chi]]^2))/\[CapitalSigma][\[Psi],\[Chi]]);


(*Momentum eqns for \[Psi], \[Chi], \[CurlyPhi], t in Mino's time - from Fujita, Hikida 0906.1420 eqn 46*)
P[\[Psi]_] = (e0/\[Mu]) *(r[\[Psi]]^2 + a^2*M^2) - a*M*(Lz0/\[Mu]);
\[CapitalPhi]r[\[Psi]_] = a*M / \[CapitalDelta][\[Psi]] * P[\[Psi]];
\[CapitalPhi]\[Theta][\[Chi]_] = (Lz0/\[Mu]) / (1 - cos\[Theta][\[Chi]]^2);
TTr[\[Psi]_] = (r[\[Psi]]^2 + a^2*M^2) / \[CapitalDelta][\[Psi]] * P[\[Psi]];
TT\[Theta][\[Chi]_] = -a^2*M^2*(e0/\[Mu])*(1 - cos\[Theta][\[Chi]]^2);
d\[Psi]d\[Lambda] = (M * Sqrt[(1 - (e0/\[Mu])^2)*((pOrbit0 / M - p3) - ecc0*(pOrbit0 / M + p3*Cos[\[Psi]]))*((pOrbit0 / M - p4) + ecc0*(pOrbit0 / M - p4*Cos[\[Psi]]))])/(1 - ecc0^2);
d\[Chi]d\[Lambda] = Sqrt[a^2*M^2 * (1 - (e0/\[Mu])^2) * (zPlus - zMinus^2*Cos[\[Chi]]^2)];
d\[CurlyPhi]d\[Lambda] = \[CapitalPhi]r[\[Psi]] + \[CapitalPhi]\[Theta][\[Chi]] - a*M*(e0/\[Mu]);
dtd\[Lambda] = TTr[\[Psi]] + TT\[Theta][\[Chi]] + a*M*(Lz0/\[Mu]);


(*Transformation to coordinate time*)
Print["-> Geodesics"];
d\[Psi]dt = Simplify[ReplaceAll[Derivative[1][\[Psi]] - d\[Psi]d\[Lambda]/dtd\[Lambda], {Derivative[1][\[Psi]]->Derivative[1][\[Psi]][tt], \[Psi]->\[Psi][tt], \[Chi]->\[Chi][tt]}]];
d\[Chi]dt = Simplify[ReplaceAll[Derivative[1][\[Chi]] - d\[Chi]d\[Lambda]/dtd\[Lambda], {Derivative[1][\[Chi]]->Derivative[1][\[Chi]][tt], \[Psi]->\[Psi][tt], \[Chi]->\[Chi][tt]}]];
d\[CurlyPhi]dt = Simplify[ReplaceAll[Derivative[1][\[CurlyPhi]] - d\[CurlyPhi]d\[Lambda]/dtd\[Lambda], {Derivative[1][\[CurlyPhi]]->Derivative[1][\[CurlyPhi]][tt], \[Psi]->\[Psi][tt], \[Chi]->\[Chi][tt]}]];


(*Momentum eqns for t, \[CurlyPhi] in proper time from stationarity and axisimmetry*)
u = {Derivative[1][t],Derivative[1][rr],Derivative[1][\[Theta]],Derivative[1][\[CurlyPhi]]};
killing1 = {1,0,0,0};
killing2 = {0,0,0,1};
metric = {{gtt, 0, 0, gt\[CurlyPhi]},{0, grr, 0, 0},{0, 0, g\[Theta]\[Theta], 0},{gt\[CurlyPhi], 0, 0, g\[CurlyPhi]\[CurlyPhi]}};
isoTime = -metric . killing1 . u;
isoRot = metric . killing2 . u;
sol = Solve[{isoTime==e/\[Mu], isoRot==Lz/\[Mu]}, {Derivative[1][t], Derivative[1][\[CurlyPhi]]}];
momentumEqt = Derivative[1][t]/.sol[[1]][[1]];
momentumEq\[CurlyPhi] = Derivative[1][\[CurlyPhi]]/.sol[[1]][[2]];
momentumEqtCalc[\[Psi]_,\[Chi]_] = Simplify[ReplaceAll[momentumEqt, {g\[CurlyPhi]\[CurlyPhi]-> g\[CurlyPhi]\[CurlyPhi][\[Psi],\[Chi]], gtt-> gtt[\[Psi],\[Chi]], gt\[CurlyPhi]-> gt\[CurlyPhi][\[Psi],\[Chi]]}]];
momentumEq\[CurlyPhi]Calc[\[Psi]_,\[Chi]_] = Simplify[ReplaceAll[momentumEq\[CurlyPhi], {g\[CurlyPhi]\[CurlyPhi]-> g\[CurlyPhi]\[CurlyPhi][\[Psi],\[Chi]], gtt-> gtt[\[Psi],\[Chi]], gt\[CurlyPhi]-> gt\[CurlyPhi][\[Psi],\[Chi]]}]];


(*Momentum eqns for r, \[Theta] in proper time from HJ equation*)
R[\[Psi]_] = \[CapitalDelta][\[Psi]]*(k*r[\[Psi]]^2 - (Lz/\[Mu] - a*M*(e/\[Mu]))^2 - Q/\[Mu]^2)+(e/\[Mu]*(r[\[Psi]]^2 + a^2*M^2) - a*M*Lz/\[Mu])^2; 
\[CapitalTheta][\[Chi]_] = Q/\[Mu]^2 + cos\[Theta][\[Chi]]^2 * ((k + (e/\[Mu])^2)*a^2*M^2 - (Lz/\[Mu])^2/(1 - cos\[Theta][\[Chi]]^2));
momentumEqrCalc[\[Psi]_, \[Chi]_] = Simplify[Sqrt[R[\[Psi]]] / \[CapitalSigma][\[Psi],\[Chi]]];
momentumEq\[Theta]Calc[\[Psi]_, \[Chi]_] = Simplify[Sqrt[\[CapitalTheta][\[Chi]]] / \[CapitalSigma][\[Psi],\[Chi]]];


(*Constraint eqn from 4-velocity conservation*)
constraintEqCalc[\[Psi]_, \[Chi]_] = Simplify[ReplaceAll[metric . u . u-k, {Derivative[1][t]->momentumEqtCalc[\[Psi],\[Chi]], Derivative[1][rr]->momentumEqrCalc[\[Psi],\[Chi]], Derivative[1][\[Theta]]->momentumEq\[Theta]Calc[\[Psi],\[Chi]], Derivative[1][\[CurlyPhi]]->momentumEq\[CurlyPhi]Calc[\[Psi],\[Chi]], grr-> grr[\[Psi],\[Chi]], g\[CurlyPhi]\[CurlyPhi]-> g\[CurlyPhi]\[CurlyPhi][\[Psi],\[Chi]], gtt-> gtt[\[Psi],\[Chi]], gt\[CurlyPhi]-> gt\[CurlyPhi][\[Psi],\[Chi]], g\[Theta]\[Theta]-> g\[Theta]\[Theta][\[Psi],\[Chi]]}]];
If[constraintEqCalc[\[Psi], \[Chi]] == 0, Print["-> Constraint equation check = TRUE"], Print["\t","-> Constraint equation check = FALSE"]];


(*E/Lz fluxes (2PN) - from Gair 0510129 eqns 39-44-45-46*)
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


(*p/ecc/i fluxes - from 2102.02713 Appendix B*)
xIotaEQ[iOrbit_] = Cos[iOrbit];
DD[rr_] = 2*M*(Q0/\[Mu]^2 + (Lz0/\[Mu] - a*M*e0/\[Mu])^2) - 2*rr*((Lz0/\[Mu])^2 + Q0/\[Mu]^2 + a^2*M^2*(1 - (e0/\[Mu])^2)) + 6*M*rr^2 - 4*rr^3*(1 - (e0/\[Mu])^2);
JE[rr_] = (4*a*M^2*(Lz0/\[Mu] - a*M*e0/\[Mu])*rr - 2*e0/\[Mu]*rr^2*(a^2*M^2 + rr^2)) / DD[rr];
JLz[rr_] = (4*M*(a*M*e0/\[Mu] - Lz0/\[Mu])*rr + 2*Lz0/\[Mu]*rr^2) / DD[rr];
JQ[rr_] = (rr^2 - 2*M*rr + a^2*M^2) / DD[rr];
xIotaDot1[pOrbit_, ecc_, iOrbit_] = (1 - xIotaEQ[iOrbit]^2) * (2*(1 - xIotaEQ[iOrbit]^2)^(1/2)*(Q0/\[Mu]^2 - a^2*M^2*(1 - (e0/\[Mu])^2)*(1 - xIotaEQ[iOrbit]^2))^(1/2)*(LzDot[pOrbit, ecc, iOrbit]/\[Mu]) - xIotaEQ[iOrbit]*(QDotMod[pOrbit, ecc, iOrbit]/\[Mu]^2) - 2*xIotaEQ[iOrbit]*(1 - xIotaEQ[iOrbit]^2)*a^2*M^2*(e0/\[Mu])*(eDot[pOrbit, ecc, iOrbit]/\[Mu])) / (2*(Q0/\[Mu]^2 - a^2*M^2*(1 - (e0/\[Mu])^2)*(1 - xIotaEQ[iOrbit]^2)^2));
xIotaDot2[pOrbit_, ecc_, iOrbit_] = (2*xIotaEQ[iOrbit]*(1 - xIotaEQ[iOrbit]^2)*(Lz0/\[Mu])*(LzDot[pOrbit, ecc, iOrbit]/\[Mu]) - xIotaEQ[iOrbit]^3*((QDotMod[pOrbit, ecc, iOrbit]/\[Mu]^2) + 2*(1 - xIotaEQ[iOrbit]^2)*a^2*M^2*(e0/\[Mu])*(eDot[pOrbit, ecc, iOrbit]/\[Mu]))) / (2*((Lz0/\[Mu])^2 + xIotaEQ[iOrbit]^4*a^2*M^2*(1 - (e0/\[Mu])^2)));
rMaxDot[rOrbitMax_, pOrbit_, ecc_, iOrbit_] = JE[rOrbitMax]*eDot[pOrbit, ecc, iOrbit]/\[Mu] + JLz[rOrbitMax]*LzDot[pOrbit, ecc, iOrbit]/\[Mu] + JQ[rOrbitMax]*QDotMod[pOrbit, ecc, iOrbit]/(\[Mu]^2);
rMinDot[rOrbitMin_, pOrbit_, ecc_, iOrbit_] = JE[rOrbitMin]*eDot[pOrbit, ecc, iOrbit]/\[Mu] + JLz[rOrbitMin]*LzDot[pOrbit, ecc, iOrbit]/\[Mu] + JQ[rOrbitMin]*QDotMod[pOrbit, ecc, iOrbit]/(\[Mu]^2);
pDot[rOrbitMax_, rOrbitMin_, pOrbit_, ecc_, iOrbit_] = ((1 - ecc)^2 / 2) * rMaxDot[rOrbitMax, pOrbit, ecc, iOrbit] + ((1 + ecc)^2 / 2) * rMinDot[rOrbitMin, pOrbit, ecc, iOrbit];
eccDot[rOrbitMax_, rOrbitMin_, pOrbit_, ecc_, iOrbit_] = ((1 - ecc^2) / (2 * pOrbit)) * ((1 - ecc) * rMaxDot[rOrbitMax, pOrbit, ecc, iOrbit] - (1 + ecc) * rMinDot[rOrbitMin, pOrbit, ecc, iOrbit]);


(*p/ecc/i fluxes in resonance*)
xIotaDotRes1[tKick_, ResFlag_, pOrbit_, ecc_, iOrbit_] = (1 - xIotaEQ[iOrbit]^2) * (2*(1 - xIotaEQ[iOrbit]^2)^(1/2)*(Q0/\[Mu]^2 - a^2*M^2*(1 - (e0/\[Mu])^2)*(1 - xIotaEQ[iOrbit]^2))^(1/2)*((LzDot[pOrbit, ecc, iOrbit] * (1 + resCoeffLz*\[Omega][tKick]*ResFlag))/\[Mu]) - xIotaEQ[iOrbit]*((QDotMod[pOrbit, ecc, iOrbit] * (1 + resCoeffQ*\[Omega][tKick]*ResFlag))/\[Mu]^2) - 2*xIotaEQ[iOrbit]*(1 - xIotaEQ[iOrbit]^2)*a^2*M^2*(e0/\[Mu])*((eDot[pOrbit, ecc, iOrbit]* (1 + resCoeffE*\[Omega][tKick]*ResFlag))/\[Mu])) / (2*(Q0/\[Mu]^2 - a^2*M^2*(1 - (e0/\[Mu])^2)*(1 - xIotaEQ[iOrbit]^2)^2));
xIotaDotRes2[tKick_, ResFlag_, pOrbit_, ecc_, iOrbit_] = (2*xIotaEQ[iOrbit]*(1 - xIotaEQ[iOrbit]^2)*(Lz0/\[Mu])*((LzDot[pOrbit, ecc, iOrbit] * (1 + resCoeffLz*\[Omega][tKick]*ResFlag))/\[Mu]) - xIotaEQ[iOrbit]^3*(((QDotMod[pOrbit, ecc, iOrbit] * (1 + resCoeffQ*\[Omega][tKick]*ResFlag))/\[Mu]^2) + 2*(1 - xIotaEQ[iOrbit]^2)*a^2*M^2*(e0/\[Mu])*((eDot[pOrbit, ecc, iOrbit]* (1 + resCoeffE*\[Omega][tKick]*ResFlag))/\[Mu]))) / (2*((Lz0/\[Mu])^2 + xIotaEQ[iOrbit]^4*a^2*M^2*(1 - (e0/\[Mu])^2)));
rMaxDotRes[tKick_, ResFlag_, rOrbitMax_, pOrbit_, ecc_, iOrbit_] = JE[rOrbitMax]*(eDot[pOrbit, ecc, iOrbit]* (1 + resCoeffE*\[Omega][tKick]*ResFlag))/\[Mu] + JLz[rOrbitMax]*(LzDot[pOrbit, ecc, iOrbit] * (1 + resCoeffLz*\[Omega][tKick]*ResFlag))/\[Mu] + JQ[rOrbitMax]*(QDotMod[pOrbit, ecc, iOrbit] * (1 + resCoeffQ*\[Omega][tKick]*ResFlag))/(\[Mu]^2);
rMinDotRes[tKick_, ResFlag_, rOrbitMin_, pOrbit_, ecc_, iOrbit_] = JE[rOrbitMin]*(eDot[pOrbit, ecc, iOrbit]* (1 + resCoeffE*\[Omega][tKick]*ResFlag))/\[Mu] + JLz[rOrbitMin]*(LzDot[pOrbit, ecc, iOrbit] * (1 + resCoeffLz*\[Omega][tKick]*ResFlag))/\[Mu] + JQ[rOrbitMin]*(QDotMod[pOrbit, ecc, iOrbit] * (1 + resCoeffQ*\[Omega][tKick]*ResFlag))/(\[Mu]^2);
pDotRes[tKick_, ResFlag_, rOrbitMax_, rOrbitMin_, pOrbit_, ecc_, iOrbit_] = ((1 - ecc)^2 / 2) * rMaxDotRes[tKick, ResFlag, rOrbitMax, pOrbit, ecc, iOrbit] + ((1 + ecc)^2 / 2) * rMinDotRes[tKick, ResFlag, rOrbitMin, pOrbit, ecc, iOrbit];
eccDotRes[tKick_, ResFlag_, rOrbitMax_, rOrbitMin_, pOrbit_, ecc_, iOrbit_] = ((1 - ecc^2) / (2 * pOrbit)) * ((1 - ecc) * rMaxDotRes[tKick, ResFlag, rOrbitMax, pOrbit, ecc, iOrbit] - (1 + ecc) * rMinDotRes[tKick, ResFlag, rOrbitMin, pOrbit, ecc, iOrbit]);


(*Eqns for rp, ra*)
rOrbitMinEq[pOrbit_,ecc_] = pOrbit / (1+ecc);
rOrbitMaxEq[pOrbit_,ecc_] = pOrbit / (1-ecc);


(*Conversion to physical units*)
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


(*Coefficients and elliptic integrals from Schmidt 0202090, e\[NotEqual]0*)
ff[\[Psi]_] = (r[\[Psi]]/M)^4 + a^2 ((r[\[Psi]]/M)((r[\[Psi]]/M)+2) + zMinus^2*(\[CapitalDelta][\[Psi]]/M^2));
gg[\[Psi]_] = 2*a*(r[\[Psi]]/M);
hh[\[Psi]_] = (r[\[Psi]]/M)((r[\[Psi]]/M)-2) + zMinus^2*(\[CapitalDelta][\[Psi]]/M^2)/(1 - zMinus^2);
dd[\[Psi]_] = ((r[\[Psi]]/M)^2 + a^2*zMinus^2)*(\[CapitalDelta][\[Psi]]/M^2);
KK[kk_] := NIntegrate[1/(1 - kk^2*Sin[thetaAus]^2)^(1/2), {thetaAus, 0, Pi/2}]; 
EE[kk_] := NIntegrate[(1 - kk^2*Sin[thetaAus]^2)^(1/2), {thetaAus, 0, Pi/2}];
\[CapitalPi][hhh_, kk_] := NIntegrate[1 / ((1 + hhh*Sin[thetaAus]^2)*(1 - kk^2*Sin[thetaAus]^2)^(1/2)), {thetaAus, 0, Pi/2}];


(* Effective Resonance Model from Speri, Gair 2103.06306 *)
Csi[omegaTheta_, omegaRadial_] = Abs[3/2 - omegaTheta/omegaRadial]; (*trigger function*)
normFactor = NIntegrate[1 + Cos[4*Pi*xAus^2], {xAus, -0.5, 0.5}];
\[Omega][tt_] = (1 + Cos[4*Pi*((tt - tRes0) / tRes)^2]) / normFactor; (*impulse function*)


(*Initial parameters*)
k = -1; 
M = 1.0;
\[Mu] = 1.0*10^-5;
Mtot = (M+\[Mu])*convMass; (*Kg*)
Mreduced = (M*\[Mu]*convMass^2)/Mtot;(*reduced mass [Kg]*)
chirpM = Mreduced^(3/5) * Mtot^(2/5); (*Kg*)
initPar = Import["/THEORY/USERS/edoardo.levati/DATA/kick/paper/part_1/32/Parameter_space_32.csv"]; (* load initial orbital parameters (a, i, ecc, p, tRes) *)
nIter = Dimensions[initPar][[1]];


(*Inspiral*)
Needs["Developer`"];
For[ii=5, ii<=IntegerPart[nIter], ii++, 
{
	Print["-> EMRI = ", ii-1, " / ", nIter-1];
	
	(* ***Initial orbital parameters (a, i, ecc, p, tResStart, tRes) *** *)
	a = initPar[[ii, 1]]; 
	iOrbit0 = initPar[[ii, 2]]; 
	ecc0 = initPar[[ii, 3]]; 
	pOrbit0 = initPar[[ii, 4]]*M; 
	tResStart = initPar[[ii, 5]]*M; 
	tRes = Round[initPar[[ii, 6]], 200]*M; (*round to 100/200 depending on the value of dttArray2 at line 302;*) 
	(*resCoeffE = -0.003; (*resonance coefficient E*)
	resCoeffLz = -0.003; (*resonance coefficient Lz*)
	resCoeffQ = -0.003; (*resonance coefficient Q*)*)
	resCoeffE = -initPar[[ii, 7]]/100; (*resonance coefficient E - Ref. 1208.3906*)
	resCoeffLz = -initPar[[ii, 8]]/100; (*resonance coefficient Lz - Ref. 1208.3906*)
	resCoeffQ = -initPar[[ii, 9]]/100; (*resonance coefficient Q - Ref. 1208.3906*)
	rOrbitMin0 = rOrbitMinEq[pOrbit0,ecc0];
	rOrbitMax0 = rOrbitMaxEq[pOrbit0,ecc0]; 
	\[Theta]OrbitMin0 = Pi/2 - iOrbit0;
	zMinus = Cos[\[Theta]OrbitMin0];
	
	(* ***E/Lz/Q from orbital elements p, e, i*** *)
	\[Psi]Periapsis = 0;
	\[Psi]Apoapsis = Pi;
	(*Coefficients - use dimensionless quantities*)
	ff1 = ff[\[Psi]Periapsis];
	gg1 = gg[\[Psi]Periapsis];
	hh1 = hh[\[Psi]Periapsis];
	dd1 = dd[\[Psi]Periapsis];
	ff2 = ff[\[Psi]Apoapsis];
	gg2 = gg[\[Psi]Apoapsis];
	hh2 = hh[\[Psi]Apoapsis];
	dd2 = dd[\[Psi]Apoapsis];
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
	AuxQ = Reduce[LzInitEq[LzInit, eInitPrograde]==0 && LzInit>=0, {LzInit}, Reals];
	LzInitPrograde = If[Head[AuxQ] === Or, AuxQ[[1]][[2]], AuxQ[[2]]];
	(*LzInitPrograde = Reduce[LzInitEq[LzInit, eInitPrograde]==0 && LzInit>=0, {LzInit}, Reals][[2]];*)
	(*LzInitRetrograde = Reduce[LzInitEq[LzInit, eInitRetrograde]==0 && LzInit<=0, {LzInit}, Reals][[2]];*)
	Lz0 = LzInitPrograde \[Mu] M;
	e0 = eInitPrograde \[Mu];
	(*Q*)
	Qdimless = zMinus^2*(a^2*(1 - (e0/\[Mu])^2) + (Lz0/(\[Mu]*M))^2/(1-zMinus^2));
	Q0 = Qdimless*M^2*\[Mu]^2;
	(*Q0 = Lz0^2 * Tan[iOrbit0]^2;*)
	
	(* ***Initial conditions for NDSolve*** *)
	\[Psi]0tt = 0; 
	\[Chi]0tt = Pi/2;
	\[CurlyPhi]0tt = 0;
	tMax = 6.4*10^6 M; (*1 year*)
	nDiscr = 10^6; (*on Vera -> Nmax = 10^9*)
	dttArray1 = 100 M; (*t0 - tResStart*)
	dttArray2 = 200 M; (*tResStart - tResEnd*)
	dttArray3 = 200 M; (*tResEnd - tMax*)
	dttSolve = 10^3 M;	
	dttDiscr = dttSolve / nDiscr; (*discretizing the interpolating functions coming from NDSolve*) 
	dttDiscrGW1 = 0.1 M;
	dttDiscrGW2 = 0.2 M;
	dttDiscrGW3 = 0.2 M;
	h0 = 4/DL * (M*convLenght) * (\[Mu]*convLenght) * 1/(rOrbitMin0*convLenght); (*GW strain from quadrupole approximation*)
	
	(* ***Print initial conditions and parameters*** *)
	initHeader = {"M = ", "\[Mu] = ", "a = ", "p = ", "ecc = ", "i = ", "e / \[Mu] = ", "Lz / (\[Mu]M) = ", "Q / (M^2 \[Mu]^2) = ", "tRes (days) = ", "tMax (days) = ", "DL (Mpc) = ", "chirp mass (Msun) = ", "h0 = "};
	initData = {M, \[Mu], Round[a, 0.001], Round[pOrbit0, 0.001], Round[ecc0, 0.001], Round[iOrbit0, 0.001], Round[e0/\[Mu], 0.001], Round[Lz0/(M \[Mu]), 0.001], Round[Q0/(M^2 \[Mu]^2), 0.001], Round[tRes*convTime/day, 0.01], Round[tMax*convTime/day], Round[DL/(3.08567758*10^16*10^6)], Round[chirpM/mSun], Round[h0, N[10^-22]]}; 
	Print["\t", "-> Initial conditions"];
	Do[Print["\t\t", initHeader[[iii]], initData[[iii]]], {iii, Length[initData]}];
	
	(* ***Dissipation loop t < tResStart*** *)
	Print["\t","-> Dissipation loop | t < tResStart"];
	tt0 = 0.0;
	ttmax = tResStart; 
	ttArray = Range[tt0, ttmax, dttArray1];
	nUpdates = IntegerPart[Round[(ttmax-tt0)/dttArray1]]; (*number of iterations used to update integrals of motion (E, Lz)*)
	For[j=1, j<=nUpdates, j++, 
	{
		Print["\t\t","update = ", j, " / ", nUpdates];
	
		(* ***Set initial values for E, Lz, Q*** *)
		If[j>1,
		{
			e0 = e;
			Lz0 = Lz;
			Q0 = Q;
			rOrbitMin0 = rOrbitMin;
			rOrbitMax0 = rOrbitMax;
			\[Theta]OrbitMin0 = \[Theta]OrbitMin;
			zMinus = Cos[\[Theta]OrbitMin0];
			pOrbit0 = pOrbit;
			ecc0 = ecc;
			iOrbit0 = iOrbit;
		}];
		
		(* ***Inizialize auxiliary functions in eqs of motion*** *)
		r1 = rOrbitMax0;
		r2 = rOrbitMin0;
		r3 = ((2*M/(1 - (e0/\[Mu])^2) - (r1 + r2)) + ((2*M/(1 - (e0/\[Mu])^2) - (r1 + r2))^2 - 4*(a^2*M^2*(Q0/\[Mu]^2)/((1 - (e0/\[Mu])^2)*r1*r2)))^(1/2)) / 2;
		r4 = (a^2*M^2*(Q0/\[Mu]^2)/((1 - (e0/\[Mu])^2)*r1*r2)) / r3;
		p3 = (r3/M)*(1 - ecc0);
		p4 = (r4/M)*(1 + ecc0);
		eps0 = a^2*M^2*(1 - (e0/\[Mu])^2) / (Lz0/\[Mu])^2;
		zPlus = (Q0/\[Mu]^2) / ((Lz0/\[Mu])^2 * eps0 * zMinus^2);
		
		(* ***Integration*** *)
		ImpRK = {"ImplicitRungeKutta", "Coefficients"->"ImplicitRungeKuttaGaussCoefficients"};
		If[j>1,
		{
			\[Psi]0tt = \[Psi]tt[ttArray[[j]]]; 
			\[Chi]0tt = \[Chi]tt[ttArray[[j]]]; 
			\[CurlyPhi]0tt = \[CurlyPhi]tt[ttArray[[j]]]; 
		}];
		eqnstt = {d\[Psi]dt==0, d\[Chi]dt==0, d\[CurlyPhi]dt==0, \[Psi][ttArray[[j]]]==\[Psi]0tt, \[Chi][ttArray[[j]]]==\[Chi]0tt, \[CurlyPhi][ttArray[[j]]]==\[CurlyPhi]0tt};
		varstt = {\[Psi], \[Chi], \[CurlyPhi]};
		orbitTimett = NDSolve[eqnstt, varstt, {tt, ttArray[[j]], ttArray[[j+1]]}, MaxSteps->\[Infinity], PrecisionGoal->MachinePrecision, Method->ImpRK];
		\[Psi]tt = \[Psi]/.orbitTimett[[1, 1]];
		\[Chi]tt = \[Chi]/.orbitTimett[[1, 2]];
		\[CurlyPhi]tt = \[CurlyPhi]/.orbitTimett[[1, 3]];
		
		(*
		(* ***Discretizing NDSolve solution*** *)
		(*f = varstt/.orbitTimett[[1]];
		data = Transpose[ToPackedArray[Parallelize[Table[f[[kindex]][Range[ttArray[[j]],  ttArray[[j+1]], dttDiscr]],{kindex, 1, 3, 1}], ProgressReporting->False]]];*)
		data = ToPackedArray[Transpose[{\[Psi]tt[Range[ttArray[[j]], ttArray[[j+1]], dttDiscr]], \[Chi]tt[Range[ttArray[[j]], ttArray[[j+1]], dttDiscr]], \[CurlyPhi]tt[Range[ttArray[[j]], ttArray[[j+1]], dttDiscr]]}]];
		data1 = r[data[[1;;-1,1]]];
		data2 = ArcCos[zMinus * Cos[data[[1;;-1,2]]]];
		data3 = data[[1;;-1,3]]; 
		
		(* ***Export discretized data*** *)
		fname = StringJoin["/THEORY/USERS/edoardo.levati/DATA/kick/kerr/data/adiabatic_", ToString[ii-1], "/data1_", ToString[j], ".h5"];
		Export[fname, data1, OverwriteTarget -> True];
		fname = StringJoin["/THEORY/USERS/edoardo.levati/DATA/kick/kerr/data/adiabatic_", ToString[ii-1], "/data2_", ToString[j], ".h5"];
		Export[fname, data2, OverwriteTarget -> True];
		fname = StringJoin["/THEORY/USERS/edoardo.levati/DATA/kick/kerr/data/adiabatic_", ToString[ii-1], "/data3_", ToString[j], ".h5"];
		Export[fname, data3, OverwriteTarget -> True];
		*)
		
		(* ***Fundamental frequencies in Kerr - from Fujita, Hikida 0906.1420*** *)
		kkr = (((r1 -r2)/(r1 -r3))*((r3 - r4)/(r2 - r4)))^(1/2);
		kktheta = zMinus^2*(Lz0/\[Mu])*(eps0/(Q0/\[Mu]^2))^(1/2);
		(*Frequencies in proper time*)
		omegaRadialProper = (Pi * ((1 - (e0/\[Mu])^2)*(r1 - r3)*(r2 - r4))^(1/2)) / (2*KK[kkr]);
		omegaThetaProper = (Pi*(Lz0/\[Mu])*(eps0*zPlus)^(1/2)) / (2*KK[kktheta]);
		(*Frequencies in coord time*)
		hhhr = (r1 - r2)/(r1 - r3); 
		hhhPlus = ((r1 - r2)*(r3 - rPlus)) / ((r1 - r3)*(r2 - rPlus)); 
		hhhMinus = ((r1 - r2)*(r3 - rMinus)) / ((r1 - r3)*(r2 - rMinus)); 
		\[CapitalGamma] = 4*M^2*(e0/\[Mu]) + (2*a^2*M^2*(e0/\[Mu])*zPlus*omegaThetaProper)/(Pi*(Lz0/\[Mu])*(eps0*zPlus)^(1/2))*(KK[kktheta] - EE[kktheta]) + (2*omegaRadialProper)/(Pi * ((1 - (e0/\[Mu])^2)*(r1 - r3)*(r2 - r4))^(1/2)) * ((e0/\[Mu])/2*((r3*(r1 + r2 + r3) - r1*r2)*KK[kkr] + (r2 - r3)*(r1 + r2 + r3 + r4)*\[CapitalPi][-hhhr, kkr] + (r1 - r3)*(r2 - r4)*EE[kkr]) + 2*M*(e0/\[Mu])*(r3*KK[kkr] + (r2 - r3)*\[CapitalPi][-hhhr, kkr]) + 2*M/(rPlus - rMinus)*( (((4*M^2*(e0/\[Mu]) - a*M*(Lz0/\[Mu]))*rPlus - 2*M^3*a^2*(e0/\[Mu]))/(r3 - rPlus)*(KK[kkr] - \[CapitalPi][-hhhPlus, kkr]*(r2 - r3)/(r2 - rPlus))) - (((4*M^2*(e0/\[Mu]) - a*M*(Lz0/\[Mu]))*rMinus - 2*M^3*a^2*(e0/\[Mu]))/(r3 - rMinus)*(KK[kkr] - \[CapitalPi][-hhhMinus, kkr]*(r2 - r3)/(r2 - rMinus)))));
		omegaRadialCoord = omegaRadialProper / \[CapitalGamma];
		omegaThetaCoord = omegaThetaProper / \[CapitalGamma];
		csiVal = Csi[omegaThetaCoord, omegaRadialCoord];
		
		(* ***Save orbital data*** *)
		tData[j] = ttArray[[j]]; 
		energyData[j] = e0; 
		angularData[j] = Lz0; 
		carterData[j] = Q0; 
		omegaRadialData[j] = omegaRadialCoord; 
		omegaThetaData[j] = omegaThetaCoord; 
		csiData[j] = csiVal;
		pData[j] = pOrbit0;
		eccData[j] = ecc0;
		iData[j] = iOrbit0;
		psiData[j] = \[Psi]tt[ttArray[[j]]];
		chiData[j] = \[Chi]tt[ttArray[[j]]];
		phiData[j] = \[CurlyPhi]tt[ttArray[[j]]];
		
		(* ***Fluxes*** *)
		eFlux = eDot[pOrbit0, ecc0, iOrbit0] * dttArray1;
		LzFlux = LzDot[pOrbit0, ecc0, iOrbit0] * dttArray1;
		QFlux = QDotMod[pOrbit0, ecc0, iOrbit0] * dttArray1;
		rOrbitMaxFlux = rMaxDot[rOrbitMax0, pOrbit0, ecc0, iOrbit0] * dttArray1;
		rOrbitMinFlux = rMinDot[rOrbitMin0, pOrbit0, ecc0, iOrbit0] * dttArray1;
		pOrbitFlux = pDot[rOrbitMax0, rOrbitMin0, pOrbit0, ecc0, iOrbit0] * dttArray1;
		eccFlux = eccDot[rOrbitMax0, rOrbitMin0, pOrbit0, ecc0, iOrbit0] * dttArray1;
		If[Abs[xIotaEQ[iOrbit0]]<=0.5, xIotaFlux = xIotaDot1[pOrbit0, ecc0, iOrbit0] * dttArray1, xIotaFlux = xIotaDot2[pOrbit0, ecc0, iOrbit0] * dttArray1];
		e = e0+eFlux; 
		Lz = Lz0+LzFlux;
		Q = Q0+QFlux;
		rOrbitMax = rOrbitMax0 + rOrbitMaxFlux;
		rOrbitMin = rOrbitMin0 + rOrbitMinFlux;
		pOrbit = pOrbit0+pOrbitFlux; 
		ecc = ecc0+eccFlux;
		xIota = xIotaEQ[iOrbit0] + xIotaFlux; 
		iOrbit = ArcCos[xIota];
		\[Theta]OrbitMin = Pi/2 - iOrbit;
	}];
	
	(* ***Dissipation loop tResStart <= t <= tResEnd*** *)
	Print["\t","-> Dissipation loop | tResStart <= t <= tResEnd"];
	nUpdatesTot = nUpdates;
	tRes0 = tResStart + tRes/2;
	tResEnd = tResStart + tRes;
	tt0 = ttArray[[-1]];
	ttmax = tResEnd;
	ttArray = Range[tt0, ttmax, dttArray2];
	nUpdates = IntegerPart[Round[(ttmax-tt0)/dttArray2]]; (*number of iterations used to update integrals of motion (E, Lz)*)
	For[j=1, j<=nUpdates, j++, 
	{
		Print["\t\t","update = ", j, " / ", nUpdates];
		
		(* ***Set initial values for E, Lz, Q*** *)
		e0 = e;
		Lz0 = Lz;
		Q0 = Q;
		rOrbitMin0 = rOrbitMin;
		rOrbitMax0 = rOrbitMax;
		\[Theta]OrbitMin0 = \[Theta]OrbitMin;
		zMinus = Cos[\[Theta]OrbitMin0];
		pOrbit0 = pOrbit;
		ecc0 = ecc;
		iOrbit0 = iOrbit;
		
		(* ***Inizialize auxiliary functions in eqs of motion*** *)
		r1 = rOrbitMax0;
		r2 = rOrbitMin0;
		r3 = ((2*M/(1 - (e0/\[Mu])^2) - (r1 + r2)) + ((2*M/(1 - (e0/\[Mu])^2) - (r1 + r2))^2 - 4*(a^2*M^2*(Q0/\[Mu]^2)/((1 - (e0/\[Mu])^2)*r1*r2)))^(1/2)) / 2;
		r4 = (a^2*M^2*(Q0/\[Mu]^2)/((1 - (e0/\[Mu])^2)*r1*r2)) / r3;
		p3 = (r3/M)*(1 - ecc0);
		p4 = (r4/M)*(1 + ecc0);
		eps0 = a^2*M^2*(1 - (e0/\[Mu])^2) / (Lz0/\[Mu])^2;
		zPlus = (Q0/\[Mu]^2) / ((Lz0/\[Mu])^2 * eps0 * zMinus^2);
		
		(* ***Integration*** *)
		ImpRK = {"ImplicitRungeKutta", "Coefficients"->"ImplicitRungeKuttaGaussCoefficients"};
		\[Psi]0tt = \[Psi]tt[ttArray[[j]]]; 
		\[Chi]0tt = \[Chi]tt[ttArray[[j]]]; 
		\[CurlyPhi]0tt = \[CurlyPhi]tt[ttArray[[j]]]; 
		eqnstt = {d\[Psi]dt==0, d\[Chi]dt==0, d\[CurlyPhi]dt==0, \[Psi][ttArray[[j]]]==\[Psi]0tt, \[Chi][ttArray[[j]]]==\[Chi]0tt, \[CurlyPhi][ttArray[[j]]]==\[CurlyPhi]0tt};
		varstt = {\[Psi], \[Chi], \[CurlyPhi]};
		orbitTimett = NDSolve[eqnstt, varstt, {tt, ttArray[[j]], ttArray[[j+1]]}, MaxSteps->\[Infinity], PrecisionGoal->MachinePrecision, Method->ImpRK];
		\[Psi]tt = \[Psi]/.orbitTimett[[1, 1]];
		\[Chi]tt = \[Chi]/.orbitTimett[[1, 2]];
		\[CurlyPhi]tt = \[CurlyPhi]/.orbitTimett[[1, 3]];
		
		(*
		(* ***Discretizing NDSolve solution*** *)
		(*f = varstt/.orbitTimett[[1]];
		data = Transpose[ToPackedArray[Parallelize[Table[f[[kindex]][Range[ttArray[[j]],  ttArray[[j+1]], dttDiscr]],{kindex, 1, 3, 1}], ProgressReporting->False]]];*)
		data = ToPackedArray[Transpose[{\[Psi]tt[Range[ttArray[[j]], ttArray[[j+1]], dttDiscr]], \[Chi]tt[Range[ttArray[[j]], ttArray[[j+1]], dttDiscr]], \[CurlyPhi]tt[Range[ttArray[[j]], ttArray[[j+1]], dttDiscr]]}]];
		data1 = r[data[[1;;-1,1]]];
		data2 = ArcCos[zMinus * Cos[data[[1;;-1,2]]]];
		data3 = data[[1;;-1,3]]; 
		
		(* ***Export discretized data*** *)
		fname = StringJoin["/THEORY/USERS/edoardo.levati/DATA/kick/kerr/data/adiabatic_", ToString[ii-1], "/data1_", ToString[j+nUpdatesTot], ".h5"];
		Export[fname, data1, OverwriteTarget -> True];
		fname = StringJoin["/THEORY/USERS/edoardo.levati/DATA/kick/kerr/data/adiabatic_", ToString[ii-1], "/data2_", ToString[j+nUpdatesTot], ".h5"];
		Export[fname, data2, OverwriteTarget -> True];
		fname = StringJoin["/THEORY/USERS/edoardo.levati/DATA/kick/kerr/data/adiabatic_", ToString[ii-1], "/data3_", ToString[j+nUpdatesTot], ".h5"];
		Export[fname, data3, OverwriteTarget -> True];
		*)
		
		(* ***Fundamental frequencies in Kerr - from Fujita, Hikida 0906.1420*** *)
		kkr = (((r1 -r2)/(r1 -r3))*((r3 - r4)/(r2 - r4)))^(1/2);
		kktheta = zMinus^2*(Lz0/\[Mu])*(eps0/(Q0/\[Mu]^2))^(1/2);
		(*Frequencies in proper time*)
		omegaRadialProper = (Pi * ((1 - (e0/\[Mu])^2)*(r1 - r3)*(r2 - r4))^(1/2)) / (2*KK[kkr]);
		omegaThetaProper = (Pi*(Lz0/\[Mu])*(eps0*zPlus)^(1/2)) / (2*KK[kktheta]);
		(*Frequencies in coord time*)
		hhhr = (r1 - r2)/(r1 - r3); 
		hhhPlus = ((r1 - r2)*(r3 - rPlus)) / ((r1 - r3)*(r2 - rPlus)); 
		hhhMinus = ((r1 - r2)*(r3 - rMinus)) / ((r1 - r3)*(r2 - rMinus)); 
		\[CapitalGamma] = 4*M^2*(e0/\[Mu]) + (2*a^2*M^2*(e0/\[Mu])*zPlus*omegaThetaProper)/(Pi*(Lz0/\[Mu])*(eps0*zPlus)^(1/2))*(KK[kktheta] - EE[kktheta]) + (2*omegaRadialProper)/(Pi * ((1 - (e0/\[Mu])^2)*(r1 - r3)*(r2 - r4))^(1/2)) * ((e0/\[Mu])/2*((r3*(r1 + r2 + r3) - r1*r2)*KK[kkr] + (r2 - r3)*(r1 + r2 + r3 + r4)*\[CapitalPi][-hhhr, kkr] + (r1 - r3)*(r2 - r4)*EE[kkr]) + 2*M*(e0/\[Mu])*(r3*KK[kkr] + (r2 - r3)*\[CapitalPi][-hhhr, kkr]) + 2*M/(rPlus - rMinus)*( (((4*M^2*(e0/\[Mu]) - a*M*(Lz0/\[Mu]))*rPlus - 2*M^3*a^2*(e0/\[Mu]))/(r3 - rPlus)*(KK[kkr] - \[CapitalPi][-hhhPlus, kkr]*(r2 - r3)/(r2 - rPlus))) - (((4*M^2*(e0/\[Mu]) - a*M*(Lz0/\[Mu]))*rMinus - 2*M^3*a^2*(e0/\[Mu]))/(r3 - rMinus)*(KK[kkr] - \[CapitalPi][-hhhMinus, kkr]*(r2 - r3)/(r2 - rMinus)))));
		omegaRadialCoord = omegaRadialProper / \[CapitalGamma];
		omegaThetaCoord = omegaThetaProper / \[CapitalGamma];
		csiVal = Csi[omegaThetaCoord, omegaRadialCoord];
		
		(* ***Save orbital data*** *)
		tData[j+nUpdatesTot] = ttArray[[j]]; 
		energyData[j+nUpdatesTot] = e0; 
		angularData[j+nUpdatesTot] = Lz0; 
		carterData[j+nUpdatesTot] = Q0; 
		omegaRadialData[j+nUpdatesTot] = omegaRadialCoord; 
		omegaThetaData[j+nUpdatesTot] = omegaThetaCoord; 
		csiData[j+nUpdatesTot] = csiVal;
		pData[j+nUpdatesTot] = pOrbit0;
		eccData[j+nUpdatesTot] = ecc0;
		iData[j+nUpdatesTot] = iOrbit0;
		psiData[j+nUpdatesTot] = \[Psi]tt[ttArray[[j]]];
		chiData[j+nUpdatesTot] = \[Chi]tt[ttArray[[j]]];
		phiData[j+nUpdatesTot] = \[CurlyPhi]tt[ttArray[[j]]];
		
		(* ***Fluxes*** *)
		If[ttArray[[j]] >= tResStart && ttArray[[j]] <= tResEnd, ResFlag = 1, ResFlag = 0];
		tKick = Min[ttArray[[j+1]], tResEnd];
		eFlux = (eDot[pOrbit0, ecc0, iOrbit0] * (1 + resCoeffE*\[Omega][tKick]*ResFlag)) * dttArray2;
		LzFlux = (LzDot[pOrbit0, ecc0, iOrbit0]* (1 + resCoeffLz*\[Omega][tKick]*ResFlag)) * dttArray2;
		QFlux = (QDotMod[pOrbit0, ecc0, iOrbit0] * (1 + resCoeffQ*\[Omega][tKick]*ResFlag)) * dttArray2;
		rOrbitMaxFlux = rMaxDotRes[tKick, ResFlag, rOrbitMax0, pOrbit0, ecc0, iOrbit0] * dttArray2;
		rOrbitMinFlux = rMinDotRes[tKick, ResFlag, rOrbitMin0, pOrbit0, ecc0, iOrbit0] * dttArray2;
		pOrbitFlux = pDotRes[tKick, ResFlag, rOrbitMax0, rOrbitMin0, pOrbit0, ecc0, iOrbit0] * dttArray2;
		eccFlux = eccDotRes[tKick, ResFlag, rOrbitMax0, rOrbitMin0, pOrbit0, ecc0, iOrbit0] * dttArray2;
		If[Abs[xIotaEQ[iOrbit0]]<=0.5, xIotaFlux = xIotaDotRes1[tKick, ResFlag, pOrbit0, ecc0, iOrbit0] * dttArray2, xIotaFlux = xIotaDotRes2[tKick, ResFlag, pOrbit0, ecc0, iOrbit0] * dttArray2];
		e = e0+eFlux; 
		Lz = Lz0+LzFlux;
		Q = Q0+QFlux;
		rOrbitMax = rOrbitMax0 + rOrbitMaxFlux;
		rOrbitMin = rOrbitMin0 + rOrbitMinFlux;
		pOrbit = pOrbit0+pOrbitFlux; 
		ecc = ecc0+eccFlux;
		xIota = xIotaEQ[iOrbit0] + xIotaFlux; 
		iOrbit = ArcCos[xIota];
		\[Theta]OrbitMin = Pi/2 - iOrbit;
		
		(* ***GWs*** *)
		rtt[t_] = pOrbit0 / (1 + ecc0*Cos[\[Psi]tt[t]]);
		\[Theta]tt[t_] = ArcCos[cos\[Theta][\[Chi]tt[t]]];
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
		
		(* ***Export GWs data*** *)
		fname = StringJoin["/THEORY/USERS/edoardo.levati/DATA/kick/paper/part_1/hdata/32/resonance_hughes_coeff_", ToString[ii-1], "/hdata_", ToString[j+nUpdatesTot], ".h5"];
		Export[fname, hdata, OverwriteTarget -> True];
	}];
	
	(* ***Dissipation loop t > tResEnd*** *)
	Print["\t","-> Dissipation loop | t > tResEnd"];
	nUpdatesTot = nUpdatesTot + nUpdates;
	tt0 = ttArray[[-1]];
	nUpdates = IntegerPart[Round[(tMax-tt0)/dttArray3]]; (*number of iterations used to update integrals of motion (E, Lz)*)
	ttmax = tt0 + nUpdates*dttArray3;
	ttArray = Range[tt0, ttmax, dttArray3];	
	For[j=1, j<=nUpdates, j++, 
	{
		Print["\t\t","update = ", j, " / ", nUpdates];
		
		(* ***Set initial values for E, Lz, Q*** *)
		e0 = e;
		Lz0 = Lz;
		Q0 = Q;
		rOrbitMin0 = rOrbitMin;
		rOrbitMax0 = rOrbitMax;
		\[Theta]OrbitMin0 = \[Theta]OrbitMin;
		zMinus = Cos[\[Theta]OrbitMin0];
		pOrbit0 = pOrbit;
		ecc0 = ecc;
		iOrbit0 = iOrbit;
		
		(* ***Inizialize auxiliary functions in eqs of motion*** *)
		r1 = rOrbitMax0;
		r2 = rOrbitMin0;
		r3 = ((2*M/(1 - (e0/\[Mu])^2) - (r1 + r2)) + ((2*M/(1 - (e0/\[Mu])^2) - (r1 + r2))^2 - 4*(a^2*M^2*(Q0/\[Mu]^2)/((1 - (e0/\[Mu])^2)*r1*r2)))^(1/2)) / 2;
		r4 = (a^2*M^2*(Q0/\[Mu]^2)/((1 - (e0/\[Mu])^2)*r1*r2)) / r3;
		p3 = (r3/M)*(1 - ecc0);
		p4 = (r4/M)*(1 + ecc0);
		eps0 = a^2*M^2*(1 - (e0/\[Mu])^2) / (Lz0/\[Mu])^2;
		zPlus = (Q0/\[Mu]^2) / ((Lz0/\[Mu])^2 * eps0 * zMinus^2);
		
		(* ***Integration*** *)
		ImpRK = {"ImplicitRungeKutta", "Coefficients"->"ImplicitRungeKuttaGaussCoefficients"};
		\[Psi]0tt = \[Psi]tt[ttArray[[j]]]; 
		\[Chi]0tt = \[Chi]tt[ttArray[[j]]]; 
		\[CurlyPhi]0tt = \[CurlyPhi]tt[ttArray[[j]]]; 
		eqnstt = {d\[Psi]dt==0, d\[Chi]dt==0, d\[CurlyPhi]dt==0, \[Psi][ttArray[[j]]]==\[Psi]0tt, \[Chi][ttArray[[j]]]==\[Chi]0tt, \[CurlyPhi][ttArray[[j]]]==\[CurlyPhi]0tt};
		varstt = {\[Psi], \[Chi], \[CurlyPhi]};
		(*orbitTimett = NDSolve[eqnstt, varstt, {tt, ttArray[[j]], ttArray[[j+1]]}, MaxSteps->\[Infinity], PrecisionGoal->MachinePrecision, Method->ImpRK];*)
		orbitTimett = Quiet[Check[NDSolve[eqnstt, varstt, {tt, ttArray[[j]], ttArray[[j+1]]}, MaxSteps->\[Infinity], PrecisionGoal->MachinePrecision, Method->ImpRK], err, {NDSolve::impsitsf, NIntegrate::slwcon, NIntegrate::ncvb, General::stop, InterpolatingFunction::dmval}]];
		If[orbitTimett == err,
		{
			Print["\t\t\t NDSolve error (integration) -> exiting dissipation loop"];
			Break[]; 
		}];
		\[Psi]tt = \[Psi]/.orbitTimett[[1, 1]];
		\[Chi]tt = \[Chi]/.orbitTimett[[1, 2]];
		\[CurlyPhi]tt = \[CurlyPhi]/.orbitTimett[[1, 3]];
		
		(*
		(* ***Discretizing NDSolve solution*** *)
		(*f = varstt/.orbitTimett[[1]];
		data = Transpose[ToPackedArray[Parallelize[Table[f[[kindex]][Range[ttArray[[j]],  ttArray[[j+1]], dttDiscr]],{kindex, 1, 3, 1}], ProgressReporting->False]]];*)
		data = ToPackedArray[Transpose[{\[Psi]tt[Range[ttArray[[j]], ttArray[[j+1]], dttDiscr]], \[Chi]tt[Range[ttArray[[j]], ttArray[[j+1]], dttDiscr]], \[CurlyPhi]tt[Range[ttArray[[j]], ttArray[[j+1]], dttDiscr]]}]];
		data1 = r[data[[1;;-1,1]]];
		data2 = ArcCos[zMinus * Cos[data[[1;;-1,2]]]];
		data3 = data[[1;;-1,3]]; 
		
		(* ***Export discretized data*** *)
		fname = StringJoin["/THEORY/USERS/edoardo.levati/DATA/kick/kerr/data/adiabatic_", ToString[ii-1], "/data1_", ToString[j+nUpdatesTot], ".h5"];
		Export[fname, data1, OverwriteTarget -> True];
		fname = StringJoin["/THEORY/USERS/edoardo.levati/DATA/kick/kerr/data/adiabatic_", ToString[ii-1], "/data2_", ToString[j+nUpdatesTot], ".h5"];
		Export[fname, data2, OverwriteTarget -> True];
		fname = StringJoin["/THEORY/USERS/edoardo.levati/DATA/kick/kerr/data/adiabatic_", ToString[ii-1], "/data3_", ToString[j+nUpdatesTot], ".h5"];
		Export[fname, data3, OverwriteTarget -> True];
		*)
		
		(* ***Fundamental frequencies in Kerr - from Fujita, Hikida 0906.1420*** *)
		kkr = (((r1 -r2)/(r1 -r3))*((r3 - r4)/(r2 - r4)))^(1/2);
		kktheta = zMinus^2*(Lz0/\[Mu])*(eps0/(Q0/\[Mu]^2))^(1/2);
		(*Frequencies in proper time*)
		omegaRadialProper = Quiet[Check[(Pi * ((1 - (e0/\[Mu])^2)*(r1 - r3)*(r2 - r4))^(1/2)) / (2*KK[kkr]), err, {NIntegrate::slwcon, NIntegrate::ncvb, General::stop}]];
		omegaThetaProper = Quiet[Check[(Pi*(Lz0/\[Mu])*(eps0*zPlus)^(1/2)) / (2*KK[kktheta]), err, {NIntegrate::slwcon, NIntegrate::ncvb, General::stop}]];
		(*Frequencies in coord time*)
		hhhr = (r1 - r2)/(r1 - r3); 
		hhhPlus = ((r1 - r2)*(r3 - rPlus)) / ((r1 - r3)*(r2 - rPlus)); 
		hhhMinus = ((r1 - r2)*(r3 - rMinus)) / ((r1 - r3)*(r2 - rMinus)); 
		\[CapitalGamma] = Quiet[Check[4*M^2*(e0/\[Mu]) + (2*a^2*M^2*(e0/\[Mu])*zPlus*omegaThetaProper)/(Pi*(Lz0/\[Mu])*(eps0*zPlus)^(1/2))*(KK[kktheta] - EE[kktheta]) + (2*omegaRadialProper)/(Pi * ((1 - (e0/\[Mu])^2)*(r1 - r3)*(r2 - r4))^(1/2)) * ((e0/\[Mu])/2*((r3*(r1 + r2 + r3) - r1*r2)*KK[kkr] + (r2 - r3)*(r1 + r2 + r3 + r4)*\[CapitalPi][-hhhr, kkr] + (r1 - r3)*(r2 - r4)*EE[kkr]) + 2*M*(e0/\[Mu])*(r3*KK[kkr] + (r2 - r3)*\[CapitalPi][-hhhr, kkr]) + 2*M/(rPlus - rMinus)*( (((4*M^2*(e0/\[Mu]) - a*M*(Lz0/\[Mu]))*rPlus - 2*M^3*a^2*(e0/\[Mu]))/(r3 - rPlus)*(KK[kkr] - \[CapitalPi][-hhhPlus, kkr]*(r2 - r3)/(r2 - rPlus))) - (((4*M^2*(e0/\[Mu]) - a*M*(Lz0/\[Mu]))*rMinus - 2*M^3*a^2*(e0/\[Mu]))/(r3 - rMinus)*(KK[kkr] - \[CapitalPi][-hhhMinus, kkr]*(r2 - r3)/(r2 - rMinus))))), err, {NIntegrate::slwcon, NIntegrate::ncvb, General::stop}]];
		If[omegaRadialProper == err || omegaThetaProper == err || \[CapitalGamma] == err,
		{
			Print["\t\t\t NIntegrate error (orbital frequencies) -> separatrix limit -> exiting dissipation loop"];
			Break[]; 
		}];
		omegaRadialCoord = omegaRadialProper / \[CapitalGamma];
		omegaThetaCoord = omegaThetaProper / \[CapitalGamma];
		csiVal = Csi[omegaThetaCoord, omegaRadialCoord];
		
		(* ***Save orbital data*** *)
		tData[j+nUpdatesTot] = ttArray[[j]]; 
		energyData[j+nUpdatesTot] = e0; 
		angularData[j+nUpdatesTot] = Lz0; 
		carterData[j+nUpdatesTot] = Q0; 
		omegaRadialData[j+nUpdatesTot] = omegaRadialCoord; 
		omegaThetaData[j+nUpdatesTot] = omegaThetaCoord; 
		csiData[j+nUpdatesTot] = csiVal;
		pData[j+nUpdatesTot] = pOrbit0;
		eccData[j+nUpdatesTot] = ecc0;
		iData[j+nUpdatesTot] = iOrbit0;
		psiData[j+nUpdatesTot] = \[Psi]tt[ttArray[[j]]];
		chiData[j+nUpdatesTot] = \[Chi]tt[ttArray[[j]]];
		phiData[j+nUpdatesTot] = \[CurlyPhi]tt[ttArray[[j]]];
		
		(* ***Fluxes*** *)
		eFlux = eDot[pOrbit0, ecc0, iOrbit0] * dttArray3;
		LzFlux = LzDot[pOrbit0, ecc0, iOrbit0] * dttArray3;
		QFlux = QDotMod[pOrbit0, ecc0, iOrbit0] * dttArray3;
		rOrbitMaxFlux = rMaxDot[rOrbitMax0, pOrbit0, ecc0, iOrbit0] * dttArray3;
		rOrbitMinFlux = rMinDot[rOrbitMin0, pOrbit0, ecc0, iOrbit0] * dttArray3;
		pOrbitFlux = pDot[rOrbitMax0, rOrbitMin0, pOrbit0, ecc0, iOrbit0] * dttArray3;
		eccFlux = eccDot[rOrbitMax0, rOrbitMin0, pOrbit0, ecc0, iOrbit0] * dttArray3;
		If[Abs[xIotaEQ[iOrbit0]]<=0.5, xIotaFlux = xIotaDot1[pOrbit0, ecc0, iOrbit0] * dttArray3, xIotaFlux = xIotaDot2[pOrbit0, ecc0, iOrbit0] * dttArray3];
		e = e0+eFlux; 
		Lz = Lz0+LzFlux;
		Q = Q0+QFlux;
		rOrbitMax = rOrbitMax0 + rOrbitMaxFlux;
		rOrbitMin = rOrbitMin0 + rOrbitMinFlux;
		pOrbit = pOrbit0 + pOrbitFlux; 
		ecc = ecc0 + eccFlux;
		xIota = xIotaEQ[iOrbit0] + xIotaFlux; 
		iOrbit = ArcCos[xIota];
		\[Theta]OrbitMin = Pi/2 - iOrbit;
		
		(* ***GWs*** *)
		rtt[t_] = pOrbit0 / (1 + ecc0*Cos[\[Psi]tt[t]]);
		\[Theta]tt[t_] = ArcCos[cos\[Theta][\[Chi]tt[t]]];
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
		
		(* ***Export GWs data*** *)
		fname = StringJoin["/THEORY/USERS/edoardo.levati/DATA/kick/paper/part_1/hdata/32/resonance_hughes_coeff_", ToString[ii-1], "/hdata_", ToString[j+nUpdatesTot], ".h5"];
		Export[fname, hdata, OverwriteTarget -> True];
		
		(* ***Check circular orbit limit*** *)
		eccLim = 0.01;
		If[(ecc < eccLim),
		{
			Print["\t\t\t Circular orbit limit -> exiting dissipation loop"];
			Break[]; 
		}];
		
		(* ***Check NK limit (rp < 5M) *** *)
		(*rOrbitLim = 5*M;
		If[rOrbitMin < rOrbitLim,
		{
			Print["\t\t\t NK limit (rp < 5M) -> exiting dissipation loop"];
			Break[]; 
		}];*)
	}];
	
	(* ***Export orbital data*** *)
	nUpdatesTot = nUpdatesTot + (j-1);
	data4 = Table[{tData[i], energyData[i], angularData[i], carterData[i], omegaRadialData[i], omegaThetaData[i], csiData[i], pData[i], eccData[i], iData[i], psiData[i], chiData[i], phiData[i]}, {i, 1, nUpdatesTot}];
	fname = StringJoin["/THEORY/USERS/edoardo.levati/DATA/kick/paper/part_1/data/32/resonance_hughes_coeff_", ToString[ii-1], ".h5"];
	Export[fname, data4, OverwriteTarget -> True];
}];
Exit[];
