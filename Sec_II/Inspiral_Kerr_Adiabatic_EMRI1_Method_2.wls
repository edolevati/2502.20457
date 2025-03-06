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


(* ***geodesic equations*** *)
geodesic = Simplify[-affine . u . u];


(*geodesics proper time*)
geot = ReplaceAll[acc[[1]]-geodesic[[1]], {Derivative[2][t]->Derivative[2][t][\[Tau]], Derivative[1][t]->Derivative[1][t][\[Tau]], Derivative[1][r]->Derivative[1][r][\[Tau]], r->r[\[Tau]], Derivative[1][\[Theta]]->Derivative[1][\[Theta]][\[Tau]], \[Theta]->\[Theta][\[Tau]], Derivative[1][\[CurlyPhi]]->Derivative[1][\[CurlyPhi]][\[Tau]]}];
geor = ReplaceAll[acc[[2]]-geodesic[[2]], {Derivative[2][r]->Derivative[2][r][\[Tau]], Derivative[1][r]->Derivative[1][r][\[Tau]], r->r[\[Tau]], Derivative[1][t]->Derivative[1][t][\[Tau]], Derivative[1][\[Theta]]->Derivative[1][\[Theta]][\[Tau]], \[Theta]->\[Theta][\[Tau]], Derivative[1][\[CurlyPhi]]->Derivative[1][\[CurlyPhi]][\[Tau]]}];
geo\[Theta] = ReplaceAll[acc[[3]]-geodesic[[3]], {Derivative[2][\[Theta]]->Derivative[2][\[Theta]][\[Tau]], Derivative[1][\[Theta]]->Derivative[1][\[Theta]][\[Tau]], \[Theta]->\[Theta][\[Tau]], Derivative[1][r]->Derivative[1][r][\[Tau]], r->r[\[Tau]], Derivative[1][t]->Derivative[1][t][\[Tau]], Derivative[1][\[CurlyPhi]]->Derivative[1][\[CurlyPhi]][\[Tau]]}];
geo\[CurlyPhi] = ReplaceAll[acc[[4]]-geodesic[[4]], {Derivative[2][\[CurlyPhi]]->Derivative[2][\[CurlyPhi]][\[Tau]], Derivative[1][\[CurlyPhi]]->Derivative[1][\[CurlyPhi]][\[Tau]], Derivative[1][t]->Derivative[1][t][\[Tau]], Derivative[1][r]->Derivative[1][r][\[Tau]], r->r[\[Tau]], Derivative[1][\[Theta]]->Derivative[1][\[Theta]][\[Tau]], \[Theta]->\[Theta][\[Tau]]}];


(*1st order GSF - PhysRevD.100.084055*) (*!!!we find different expressions for Fr/F\[Theta]!!!*)
FtEq = ((eDot2PN+(2 M r[\[Tau]] (a M (-LzDot2PN+a eDot2PN M)+eDot2PN r[\[Tau]]^2))/((a^2 M^2 Cos[\[Theta][\[Tau]]]^2+r[\[Tau]]^2) (a^2 M^2-2 M r[\[Tau]]+r[\[Tau]]^2))) Derivative[1][t][\[Tau]])/\[Mu];
F\[CurlyPhi]Eq = ((a^2 LzDot2PN M^2 Cot[\[Theta][\[Tau]]]^2+r[\[Tau]] (2 a eDot2PN M^2+LzDot2PN Csc[\[Theta][\[Tau]]]^2 (-2 M+r[\[Tau]]))) Derivative[1][t][\[Tau]])/(\[Mu] (a^2 M^2 Cos[\[Theta][\[Tau]]]^2+r[\[Tau]]^2) (a^2 M^2-2 M r[\[Tau]]+r[\[Tau]]^2));
FrEq = -((-4 \[Mu] (a^2 M^2-2 M r[\[Tau]]+r[\[Tau]]^2) Derivative[1][t][\[Tau]]^2 (r[\[Tau]] (2 a M^2 (-LzDot2PN+a eDot2PN M) Cos[\[Theta][\[Tau]]]^2+a^2 eDot2PN M^2 Cos[\[Theta][\[Tau]]]^2 r[\[Tau]]+eDot2PN r[\[Tau]]^3)+a^2 M^2 \[Mu] Cos[\[Theta][\[Tau]]] (a^2 M^2 Cos[\[Theta][\[Tau]]]^2-2 M r[\[Tau]]+r[\[Tau]]^2) Sin[\[Theta][\[Tau]]] Derivative[1][\[Theta]][\[Tau]])+(a^2 M^2-2 M r[\[Tau]]+r[\[Tau]]^2) Derivative[1][t][\[Tau]] (2 QDot2PN (a^2 M^2 Cos[\[Theta][\[Tau]]]^2+r[\[Tau]]^2)+\[Mu] r[\[Tau]] (4 LzDot2PN r[\[Tau]]^3 Sin[\[Theta][\[Tau]]]^2+a^2 LzDot2PN M^2 r[\[Tau]] Sin[2 \[Theta][\[Tau]]]^2+8 a^2 M^3 Cos[\[Theta][\[Tau]]] Sin[\[Theta][\[Tau]]]^2 ((-LzDot2PN+a eDot2PN M) Cos[\[Theta][\[Tau]]]-2 a M \[Mu] Sin[\[Theta][\[Tau]]] Derivative[1][\[Theta]][\[Tau]])) Derivative[1][\[CurlyPhi]][\[Tau]])+a^2 M^2 \[Mu]^2 Sin[2 \[Theta][\[Tau]]] Derivative[1][\[Theta]][\[Tau]] (-4 M r[\[Tau]]^5 (Derivative[1][\[Theta]][\[Tau]]^2+Sin[\[Theta][\[Tau]]]^2 Derivative[1][\[CurlyPhi]][\[Tau]]^2)+2 r[\[Tau]]^6 (Derivative[1][\[Theta]][\[Tau]]^2+Sin[\[Theta][\[Tau]]]^2 Derivative[1][\[CurlyPhi]][\[Tau]]^2)+2 a^4 M^4 Cos[\[Theta][\[Tau]]]^2 (1+Cos[\[Theta][\[Tau]]]^2 (Derivative[1][r][\[Tau]]^2+a^2 M^2 Derivative[1][\[Theta]][\[Tau]]^2)+a^2 M^2 Sin[\[Theta][\[Tau]]]^2 Derivative[1][\[CurlyPhi]][\[Tau]]^2)+a^2 M^2 r[\[Tau]]^2 (3+Cos[2 \[Theta][\[Tau]]]+Cos[\[Theta][\[Tau]]]^2 (4 Derivative[1][r][\[Tau]]^2+a^2 M^2 (5+Cos[2 \[Theta][\[Tau]]]) Derivative[1][\[Theta]][\[Tau]]^2)+2 M^2 (2 (-1+a^2)+(2+a^2) Cos[2 \[Theta][\[Tau]]]) Sin[\[Theta][\[Tau]]]^2 Derivative[1][\[CurlyPhi]][\[Tau]]^2)-4 a^2 M^3 r[\[Tau]] (Cos[\[Theta][\[Tau]]]^2+a^2 M^2 (Cos[\[Theta][\[Tau]]]^4 Derivative[1][\[Theta]][\[Tau]]^2+Cos[2 \[Theta][\[Tau]]] Sin[\[Theta][\[Tau]]]^2 Derivative[1][\[CurlyPhi]][\[Tau]]^2))+r[\[Tau]]^4 (2+2 Derivative[1][r][\[Tau]]^2+a^2 M^2 (2 (2+Cos[2 \[Theta][\[Tau]]]) Derivative[1][\[Theta]][\[Tau]]^2+(5+Cos[2 \[Theta][\[Tau]]]) Sin[\[Theta][\[Tau]]]^2 Derivative[1][\[CurlyPhi]][\[Tau]]^2))-2 M r[\[Tau]]^3 (2+a^2 M^2 (4 Cos[\[Theta][\[Tau]]]^2 Derivative[1][\[Theta]][\[Tau]]^2+Sin[2 \[Theta][\[Tau]]]^2 Derivative[1][\[CurlyPhi]][\[Tau]]^2))))/(4 \[Mu]^2 (a^2 M^2 Cos[\[Theta][\[Tau]]]^2+r[\[Tau]]^2)^3 Derivative[1][r][\[Tau]]));
F\[Theta]Eq = (4 a M^2 \[Mu] Cos[\[Theta][\[Tau]]] (a^2 M^2-2 M r[\[Tau]]+r[\[Tau]]^2) Derivative[1][t][\[Tau]]^2 (a^3 eDot2PN M^2 Cos[\[Theta][\[Tau]]]^3+Cos[\[Theta][\[Tau]]] r[\[Tau]] (2 LzDot2PN-2 a eDot2PN M+a eDot2PN r[\[Tau]])-a \[Mu] (a^2 M^2 Cos[\[Theta][\[Tau]]]^2-2 M r[\[Tau]]+r[\[Tau]]^2) Sin[\[Theta][\[Tau]]] Derivative[1][\[Theta]][\[Tau]])-2 (a^2 M^2-2 M r[\[Tau]]+r[\[Tau]]^2) Derivative[1][t][\[Tau]] (-QDot2PN (a^2 M^2 Cos[\[Theta][\[Tau]]]^2+r[\[Tau]]^2)+2 \[Mu] Cos[\[Theta][\[Tau]]] (Cos[\[Theta][\[Tau]]] (LzDot2PN (a^2 M^2+r[\[Tau]]^2) (a^2 M^2 Cos[\[Theta][\[Tau]]]^2+r[\[Tau]]^2)+2 a^2 M^3 (LzDot2PN-a eDot2PN M) r[\[Tau]] Sin[\[Theta][\[Tau]]]^2)+4 a^3 M^4 \[Mu] r[\[Tau]] Sin[\[Theta][\[Tau]]]^3 Derivative[1][\[Theta]][\[Tau]]) Derivative[1][\[CurlyPhi]][\[Tau]])+a^2 M^2 \[Mu]^2 Sin[2 \[Theta][\[Tau]]] Derivative[1][\[Theta]][\[Tau]] (-4 M r[\[Tau]]^5 (Derivative[1][\[Theta]][\[Tau]]^2+Sin[\[Theta][\[Tau]]]^2 Derivative[1][\[CurlyPhi]][\[Tau]]^2)+2 r[\[Tau]]^6 (Derivative[1][\[Theta]][\[Tau]]^2+Sin[\[Theta][\[Tau]]]^2 Derivative[1][\[CurlyPhi]][\[Tau]]^2)+2 a^4 M^4 Cos[\[Theta][\[Tau]]]^2 (1+Cos[\[Theta][\[Tau]]]^2 (Derivative[1][r][\[Tau]]^2+a^2 M^2 Derivative[1][\[Theta]][\[Tau]]^2)+a^2 M^2 Sin[\[Theta][\[Tau]]]^2 Derivative[1][\[CurlyPhi]][\[Tau]]^2)+a^2 M^2 r[\[Tau]]^2 (3+Cos[2 \[Theta][\[Tau]]]+Cos[\[Theta][\[Tau]]]^2 (4 Derivative[1][r][\[Tau]]^2+a^2 M^2 (5+Cos[2 \[Theta][\[Tau]]]) Derivative[1][\[Theta]][\[Tau]]^2)+2 M^2 (2 (-1+a^2)+(2+a^2) Cos[2 \[Theta][\[Tau]]]) Sin[\[Theta][\[Tau]]]^2 Derivative[1][\[CurlyPhi]][\[Tau]]^2)-4 a^2 M^3 r[\[Tau]] (Cos[\[Theta][\[Tau]]]^2+a^2 M^2 (Cos[\[Theta][\[Tau]]]^4 Derivative[1][\[Theta]][\[Tau]]^2+Cos[2 \[Theta][\[Tau]]] Sin[\[Theta][\[Tau]]]^2 Derivative[1][\[CurlyPhi]][\[Tau]]^2))+r[\[Tau]]^4 (2+2 Derivative[1][r][\[Tau]]^2+a^2 M^2 (2 (2+Cos[2 \[Theta][\[Tau]]]) Derivative[1][\[Theta]][\[Tau]]^2+(5+Cos[2 \[Theta][\[Tau]]]) Sin[\[Theta][\[Tau]]]^2 Derivative[1][\[CurlyPhi]][\[Tau]]^2))-2 M r[\[Tau]]^3 (2+a^2 M^2 (4 Cos[\[Theta][\[Tau]]]^2 Derivative[1][\[Theta]][\[Tau]]^2+Sin[2 \[Theta][\[Tau]]]^2 Derivative[1][\[CurlyPhi]][\[Tau]]^2))))/(4 \[Mu]^2 (a^2 M^2 Cos[\[Theta][\[Tau]]]^2+r[\[Tau]]^2)^3 (a^2 M^2-2 M r[\[Tau]]+r[\[Tau]]^2) Derivative[1][\[Theta]][\[Tau]]);


(*geodesics + GSF*)
geotGSF = Simplify[geot - FtEq];
georGSF = Simplify[geor - FrEq];
geo\[Theta]GSF = Simplify[geo\[Theta] - F\[Theta]Eq];
geo\[CurlyPhi]GSF = Simplify[geo\[CurlyPhi] - F\[CurlyPhi]Eq];


(*momentum equations for t,\[CurlyPhi] from stationarity and axisimmetry*)
killing1 = {1, 0, 0, 0};
killing2 = {0, 0, 0, 1};
metric2 = {{gtt, 0, 0, gt\[CurlyPhi]},{0, grr, 0, 0},{0, 0, g\[Theta]\[Theta], 0},{gt\[CurlyPhi], 0, 0, g\[CurlyPhi]\[CurlyPhi]}};
isoTime = -metric2 . killing1 . u;
isoRot = metric2 . killing2 . u;
sol = Solve[{isoTime==e/\[Mu], isoRot==Lz/\[Mu]}, {Derivative[1][t], Derivative[1][\[CurlyPhi]]}];
momentumEqt = Derivative[1][t]/.sol[[1]][[1]];
momentumEq\[CurlyPhi] = Derivative[1][\[CurlyPhi]]/.sol[[1]][[2]];
momentumEqtCalc[r_,\[Theta]_] = ReplaceAll[momentumEqt, {g\[CurlyPhi]\[CurlyPhi]-> g\[CurlyPhi]\[CurlyPhi][r,\[Theta]], gtt-> gtt[r,\[Theta]], gt\[CurlyPhi]-> gt\[CurlyPhi][r,\[Theta]]}];
momentumEq\[CurlyPhi]Calc[r_,\[Theta]_] = ReplaceAll[momentumEq\[CurlyPhi], {g\[CurlyPhi]\[CurlyPhi]-> g\[CurlyPhi]\[CurlyPhi][r,\[Theta]], gtt-> gtt[r,\[Theta]], gt\[CurlyPhi]-> gt\[CurlyPhi][r,\[Theta]]}];


(*integrals of motion*)
energyEq[r_, \[Theta]_, vt_, v\[CurlyPhi]_] = Simplify[ReplaceAll[-metric . killing1 . u, {Derivative[1][t]->vt, Derivative[1][\[CurlyPhi]]->v\[CurlyPhi]}]];
angularEq[r_, \[Theta]_, vt_, v\[CurlyPhi]_] = Simplify[ReplaceAll[metric . killing2 . u, {Derivative[1][t]->vt, Derivative[1][\[CurlyPhi]]->v\[CurlyPhi]}]];
carterEq[r_, \[Theta]_, v\[Theta]_, vt_, v\[CurlyPhi]_] = Simplify[g\[Theta]\[Theta][r,\[Theta]]^2*v\[Theta]^2 + Cos[\[Theta]]^2*(a^2*M^2*(1 - (energyEq[r, \[Theta], vt, v\[CurlyPhi]])^2) + (angularEq[r, \[Theta], vt, v\[CurlyPhi]])^2 / Sin[\[Theta]]^2)];
constraintEq[r_, \[Theta]_, vt_, vr_, v\[Theta]_, v\[CurlyPhi]_] = Simplify[ReplaceAll[metric . u . u-k, {Derivative[1][t]->vt, Derivative[1][r]->vr, Derivative[1][\[Theta]]->v\[Theta], Derivative[1][\[CurlyPhi]]->v\[CurlyPhi]}]];


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
(*normFactor = NIntegrate[1 + Cos[4*Pi*xAus^2], {xAus, -0.5, 0.5}];
\[Omega][tt_] = (1 + Cos[4*Pi*((tt - tRes0) / tRes)^2]) / normFactor; (*impulse function*)*)


(*initial parameters*)
k = -1; 
M = 1.0;
\[Mu] = 3.0*10^-6;
Mtot = (M+\[Mu])*convMass; (*Kg*)
Mreduced = (M*\[Mu]*convMass^2)/Mtot;(*reduced mass [Kg]*)
chirpM = Mreduced^(3/5) * Mtot^(2/5); (*Kg*)
(*initPar = Import["/THEORY/USERS/edoardo.levati/DATA/kick/kerr/parameter_space_final.csv"]; (* load initial orbital parameters (a, i, ecc, p, tRes) *)
nIter = Dimensions[initPar][[1]];*)


(*inspiral*)
Needs["Developer`"];
(*LaunchKernels[8];*)
For[ii=1, ii<=(*IntegerPart[nIter]*) 1, ii++, 
{
	Print["-> EMRI = ", 1, " / ", 1];
	
	(* ***initial orbital parameters (a, i, ecc, p, tResStart, tRes) *** *)
	a = 0.95; 
	iOrbitInit = 1.2; 
	eccInit = 0.7; 
	pOrbitInit = 8.493749999996451 M; 
	tResStart = 509800 M; 
	tRes = 10932.793661964544 M;  
	
	(* ***E/Lz from orbital elements p, e, i*** *)
	pOrbit = pOrbitInit;
	ecc = eccInit;
	iOrbit = iOrbitInit;
	rOrbitMin = rOrbitMinEq[pOrbit,ecc];
	rOrbitMax = rOrbitMaxEq[pOrbit,ecc];
	\[Theta]OrbitMin = Pi/2 - iOrbit;
	zMinus = Cos[\[Theta]OrbitMin];
	(*coefficients - use dimensionless quantities*)
	ff1 = ff[rOrbitMin/M];
	gg1 = gg[rOrbitMin/M];
	hh1 = hh[rOrbitMin/M];
	dd1 = dd[rOrbitMin/M];
	ff2 = ff[rOrbitMax/M];
	gg2 = gg[rOrbitMax/M];
	hh2 = hh[rOrbitMax/M];
	dd2 = dd[rOrbitMax/M];
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
	t0 = 0.0;
	r0 = rOrbitMin; (*(rOrbitMin + rOrbitMax)/2;*)
	\[Theta]0 = Pi/2; 
	\[CurlyPhi]0 = 0.0;
	vt0 = momentumEqtCalc[r0,\[Theta]0];
	v\[CurlyPhi]0 = momentumEq\[CurlyPhi]Calc[r0,\[Theta]0];
	v\[Theta]SignCond = (v\[Theta]0Sol>0);
	v\[Theta]0 = Solve[carterEq[r0, \[Theta]0, v\[Theta]0Sol, vt0, v\[CurlyPhi]0] - Q/\[Mu]^2==0 && v\[Theta]SignCond, {v\[Theta]0Sol}, Reals][[1, 1, 2]];
	vrSignCond = (vr0Sol>0);
	(*vr0 = Solve[constraintEq[r0, \[Theta]0, vt0, vr0Sol, v\[Theta]0, v\[CurlyPhi]0] == 0 && vrSignCond, {vr0Sol}, Reals][[1, 1, 2]];*)
	vr0 = Reduce[constraintEq[r0, \[Theta]0, vt0, vr0Sol, v\[Theta]0, v\[CurlyPhi]0] <= 10^{-14} && vrSignCond, {vr0Sol}, Reals][[5]];
	\[Tau]Max = 0.6*10^6 M;
	nDiscr = 10^6; (*on Vera -> Nmax = 10^9*)
	d\[Tau]Array1 = 100 M; (*approach resonance*)
	d\[Tau]Array2 = 100 M; (*in resonance*)
	d\[Tau]Array3 = 200 M; (*from tResEnd to tMax*)
	d\[Tau]Solve = 10^3 M;	
	d\[Tau]Discr = d\[Tau]Solve / nDiscr; (*discretizing the interpolating functions coming from NDSolve*) 
	
	(* ***integration methods for geodesics + GSF*** *)
	SWITCH1 = {"StiffnessSwitching", Method->{"ExplicitRungeKutta", "LinearlyImplicitMidpoint"}};
	SWITCH2 = {"StiffnessSwitching", Method->{"ExplicitRungeKutta", Automatic}};
	SWITCH3 = {"StiffnessSwitching", Method->{Automatic, "LinearlyImplicitMidpoint"}};
	SWITCH4 = {"StiffnessSwitching", Method->{Automatic, Automatic}};
	SWITCH = {SWITCH1, SWITCH2, SWITCH3, SWITCH4};
	
	(* ***print initial conditions and parameters*** *)
	initHeader = {"M = ", "\[Mu] = ", "a = ", "p = ", "ecc = ", "i = ", "e / \[Mu] = ", "Lz / \[Mu] = ", "Q / \[Mu]^2 = ", "tRes (hours) = ", "\[Tau]_max (days) = "};
	initData = {M, \[Mu], Round[a, 0.001], Round[pOrbitInit, 0.001], Round[eccInit, 0.001], Round[iOrbitInit, 0.001], Round[e/\[Mu], 0.001], Round[Lz/(\[Mu]), 0.001], Round[Q/(\[Mu]^2), 0.001], Round[tRes*convTime/(60*60), 0.01], Round[\[Tau]Max*convTime/day]}; 
	Print["\t", "-> Initial conditions"];
	Do[Print["\t\t", initHeader[[iii]], initData[[iii]]], {iii, Length[initData]}];
	
	(* ***dissipation loop t < tResStart*** *)
	Print["\t","-> Dissipation loop | t < tResStart"];
	\[Tau]0 = 0.0;
	\[Tau]max = tResStart; 
	\[Tau]Array = Range[\[Tau]0, \[Tau]max, d\[Tau]Array1];
	nUpdates = IntegerPart[Round[(\[Tau]max-\[Tau]0)/d\[Tau]Array1]]; (*number of iterations used to update integrals of motion (E, Lz)*)
	For[j=1, j<=nUpdates, j++, 
	{
		Print["\t\t","update = ", j, " / ", nUpdates];
	
		(* ***initial E, Lz, Q*** *)
		e0 = e;
		Lz0 = Lz;
		Q0 = Q;
		
		(* ***integrate geodesics to compute orbital elements*** *)
		If[j>1,
		{
			Print["\t\t\t Compute orbital elements"];
			
			(* ***initial positions/velocities*** *)
			t0GEO = t\[Tau][\[Tau]Array[[j]]];
			r0GEO = r\[Tau][\[Tau]Array[[j]]];
			\[Theta]0GEO = \[Theta]\[Tau][\[Tau]Array[[j]]];
			\[CurlyPhi]0GEO = \[CurlyPhi]\[Tau][\[Tau]Array[[j]]];	
			vt0GEO = vt\[Tau][\[Tau]Array[[j]]];
			vr0GEO = vr\[Tau][\[Tau]Array[[j]]];
			v\[Theta]0GEO = v\[Theta]\[Tau][\[Tau]Array[[j]]];
			v\[CurlyPhi]0GEO = v\[CurlyPhi]\[Tau][\[Tau]Array[[j]]];
			
			(* ***integration*** *)
			invariant = {energyEq[r[\[Tau]], \[Theta][\[Tau]], Derivative[1][t][\[Tau]], Derivative[1][\[CurlyPhi]][\[Tau]]], angularEq[r[\[Tau]], \[Theta][\[Tau]], Derivative[1][t][\[Tau]], Derivative[1][\[CurlyPhi]][\[Tau]]], carterEq[r[\[Tau]], \[Theta][\[Tau]], Derivative[1][\[Theta]][\[Tau]], Derivative[1][t][\[Tau]], Derivative[1][\[CurlyPhi]][\[Tau]]], constraintEq[r[\[Tau]], \[Theta][\[Tau]], Derivative[1][t][\[Tau]], Derivative[1][r][\[Tau]], Derivative[1][\[Theta]][\[Tau]], Derivative[1][\[CurlyPhi]][\[Tau]]]};
			ProjImpRK = {"Projection", Method->"ImplicitRungeKutta", "Invariants"->invariant};
			eqns\[Tau]GEO = {geot==0, geor==0, geo\[Theta]==0, geo\[CurlyPhi]==0, t[\[Tau]Array[[j]]]==t0GEO, Derivative[1][t][\[Tau]Array[[j]]]==vt0GEO, r[\[Tau]Array[[j]]]==r0GEO, Derivative[1][r][\[Tau]Array[[j]]]==vr0GEO, \[Theta][\[Tau]Array[[j]]]==\[Theta]0GEO, Derivative[1][\[Theta]][\[Tau]Array[[j]]]==v\[Theta]0GEO, \[CurlyPhi][\[Tau]Array[[j]]]==\[CurlyPhi]0GEO, Derivative[1][\[CurlyPhi]][\[Tau]Array[[j]]]==v\[CurlyPhi]0GEO};
			vars\[Tau]GEO = {t, r, \[Theta], \[CurlyPhi], Derivative[1][t], Derivative[1][r], Derivative[1][\[Theta]], Derivative[1][\[CurlyPhi]]};
			orbitTime\[Tau]GEO = NDSolve[eqns\[Tau]GEO, vars\[Tau]GEO, {\[Tau], \[Tau]Array[[j]], \[Tau]Array[[j]] + d\[Tau]Solve}, MaxSteps->\[Infinity], PrecisionGoal->MachinePrecision, Method->ProjImpRK];
			t\[Tau]GEO = t/.orbitTime\[Tau]GEO[[1, 1]];
			r\[Tau]GEO = r/.orbitTime\[Tau]GEO[[1, 2]];
			\[Theta]\[Tau]GEO = \[Theta]/.orbitTime\[Tau]GEO[[1, 3]];
			\[CurlyPhi]\[Tau]GEO = \[CurlyPhi]/.orbitTime\[Tau]GEO[[1, 4]];
			vt\[Tau]GEO = Derivative[1][t]/.orbitTime\[Tau]GEO[[1, 5]];
			vr\[Tau]GEO = Derivative[1][r]/.orbitTime\[Tau]GEO[[1, 6]];
			v\[Theta]\[Tau]GEO = Derivative[1][\[Theta]]/.orbitTime\[Tau]GEO[[1, 7]];
			v\[CurlyPhi]\[Tau]GEO = Derivative[1][\[CurlyPhi]]/.orbitTime\[Tau]GEO[[1, 8]];
		
			(* ***discretizing NDSolve solution*** *)
			(*f = vars\[Tau]GEO/.orbitTime\[Tau]GEO[[1]];
			data = Transpose[ToPackedArray[Parallelize[Table[f[[kindex]][Range[\[Tau]Array[[j]], \[Tau]Array[[j]] + d\[Tau]Solve, d\[Tau]Discr]],{kindex, 1, 8, 1}], ProgressReporting->False]]];*)
			data = ToPackedArray[Transpose[{r\[Tau]GEO[Range[\[Tau]Array[[j]], \[Tau]Array[[j]] + d\[Tau]Solve, d\[Tau]Discr]], \[Theta]\[Tau]GEO[Range[\[Tau]Array[[j]], \[Tau]Array[[j]] + d\[Tau]Solve, d\[Tau]Discr]]}]];
		
			(* ***orbital elements from discretized data*** *)
			rOrbitMin = Min[data[[1;;-1,1]]];
			rOrbitMax = Max[data[[1;;-1,1]]];
			pOrbit = (2*rOrbitMin*rOrbitMax) / (rOrbitMax + rOrbitMin);
			ecc = (rOrbitMax - rOrbitMin) / (rOrbitMax + rOrbitMin);
			\[Theta]OrbitMin = Min[data[[1;;-1,2]]];
			\[Theta]OrbitMax = Max[data[[1;;-1,2]]];
			iOrbit = (\[Theta]OrbitMax - \[Theta]OrbitMin)/2.;
		}];
		
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
		
		(* ***integrate geodesics + GSF to compute trajectory*** *)
		Print["\t\t\t Compute trajectory"];
		eDot2PN = eDot[pOrbit, ecc, iOrbit];
		LzDot2PN = LzDot[pOrbit, ecc, iOrbit];
		QDot2PN = QDotMod[pOrbit, ecc, iOrbit];
		If[j>1,
		{
			(* ***initial positions/velocities*** *)
			t0 = t\[Tau][\[Tau]Array[[j]]];
			r0 = r\[Tau][\[Tau]Array[[j]]];
			\[Theta]0 = \[Theta]\[Tau][\[Tau]Array[[j]]];
			\[CurlyPhi]0 = \[CurlyPhi]\[Tau][\[Tau]Array[[j]]];	
			vt0 = vt\[Tau][\[Tau]Array[[j]]];
			vr0 = vr\[Tau][\[Tau]Array[[j]]];
			v\[Theta]0 = v\[Theta]\[Tau][\[Tau]Array[[j]]];
			v\[CurlyPhi]0 = v\[CurlyPhi]\[Tau][\[Tau]Array[[j]]];
		}];		
		Print["\t\t\t\t energy err = ", energyEq[r0, \[Theta]0, vt0, v\[CurlyPhi]0] - e0/\[Mu]];
		Print["\t\t\t\t angular err = ", angularEq[r0, \[Theta]0, vt0, v\[CurlyPhi]0] - Lz0/\[Mu]];
		Print["\t\t\t\t carter err = ", carterEq[r0, \[Theta]0, v\[Theta]0, vt0, v\[CurlyPhi]0] - Q0/\[Mu]^2];
		Print["\t\t\t\t 4-velocity err = ", constraintEq[r0, \[Theta]0, vt0, vr0, v\[Theta]0, v\[CurlyPhi]0]];
		eqns\[Tau] = {geotGSF==0, georGSF==0, geo\[Theta]GSF==0, geo\[CurlyPhi]GSF==0, t[\[Tau]Array[[j]]]==t0, Derivative[1][t][\[Tau]Array[[j]]]==vt0, r[\[Tau]Array[[j]]]==r0, Derivative[1][r][\[Tau]Array[[j]]]==vr0, \[Theta][\[Tau]Array[[j]]]==\[Theta]0, Derivative[1][\[Theta]][\[Tau]Array[[j]]]==v\[Theta]0, \[CurlyPhi][\[Tau]Array[[j]]]==\[CurlyPhi]0, Derivative[1][\[CurlyPhi]][\[Tau]Array[[j]]]==v\[CurlyPhi]0};
		vars\[Tau] = {t, r, \[Theta], \[CurlyPhi], Derivative[1][t], Derivative[1][r], Derivative[1][\[Theta]], Derivative[1][\[CurlyPhi]]};
		jMethod = 1;
		orbitTime\[Tau] = Quiet[Check[NDSolve[eqns\[Tau], vars\[Tau], {\[Tau], \[Tau]Array[[j]], \[Tau]Array[[j+1]]}, MaxSteps->\[Infinity], PrecisionGoal->MachinePrecision, Method->SWITCH[[jMethod]]], err, {NDSolve::ndsz, General::stop, InterpolatingFunction::dmval}]];
		While[orbitTime\[Tau] == err && jMethod <= 4,
		{
			jMethod += 1; 
			orbitTime\[Tau] = Quiet[Check[NDSolve[eqns\[Tau], vars\[Tau], {\[Tau], \[Tau]Array[[j]], \[Tau]Array[[j+1]]}, MaxSteps->\[Infinity], PrecisionGoal->MachinePrecision, Method->SWITCH[[jMethod]]], err, {NDSolve::ndsz, General::stop, InterpolatingFunction::dmval}]];
		}];
		If[orbitTime\[Tau] == err,
		{
			Print["\t\t\t", "NDSolve error -> try different integration method"];
			Break[]; 
		}];
		t\[Tau] = t/.orbitTime\[Tau][[1, 1]];
		r\[Tau] = r/.orbitTime\[Tau][[1, 2]];
		\[Theta]\[Tau] = \[Theta]/.orbitTime\[Tau][[1, 3]];
		\[CurlyPhi]\[Tau] = \[CurlyPhi]/.orbitTime\[Tau][[1, 4]];
		vt\[Tau] = Derivative[1][t]/.orbitTime\[Tau][[1, 5]];
		vr\[Tau] = Derivative[1][r]/.orbitTime\[Tau][[1, 6]];
		v\[Theta]\[Tau] = Derivative[1][\[Theta]]/.orbitTime\[Tau][[1, 7]];
		v\[CurlyPhi]\[Tau] = Derivative[1][\[CurlyPhi]]/.orbitTime\[Tau][[1, 8]];
		
		(* ***save energy and orbital data*** *)
		tData[j] = t\[Tau][\[Tau]Array[[j]]];
		energyData[j] = energyEq[r0, \[Theta]0, vt0, v\[CurlyPhi]0]*\[Mu]; 
		angularData[j] = angularEq[r0, \[Theta]0, vt0, v\[CurlyPhi]0]*\[Mu]; 
		carterData[j] = carterEq[r0, \[Theta]0, v\[Theta]0, vt0, v\[CurlyPhi]0]*\[Mu]^2; 
		omegaRadialData[j] = omegaRadialCoord; 
		omegaThetaData[j] = omegaThetaCoord; 
		csiData[j] = csiVal;
		pData[j] = pOrbit;
		eccData[j] = ecc;
		iData[j] = iOrbit;
		
		(* ***fluxes*** *)
		dtFlux = t\[Tau][\[Tau]Array[[j+1]]] - t\[Tau][\[Tau]Array[[j]]];
		eFlux = eDot[pOrbit, ecc, iOrbit] * dtFlux;
		LzFlux = LzDot[pOrbit, ecc, iOrbit] * dtFlux;
		QFlux = QDotMod[pOrbit, ecc, iOrbit] * dtFlux;
		e = e0+eFlux; 
		Lz = Lz0+LzFlux;
		Q = Q0+QFlux;
	}];
	
	(* ***dissipation loop tResStart <= t <= tResEnd*** *)
	Print["\t","-> Dissipation loop | tResStart <= t <= tResEnd"];
	nUpdatesTot = nUpdates;
	tRes0 = tResStart + tRes/2;
	tResEnd = tResStart + tRes;
	\[Tau]0 = \[Tau]Array[[-1]];
	\[Tau]max = Round[tResEnd, 100]; 
	\[Tau]Array = Range[\[Tau]0, \[Tau]max, d\[Tau]Array2];
	nUpdates = IntegerPart[Round[(\[Tau]max-\[Tau]0)/d\[Tau]Array2]]; (*number of iterations used to update integrals of motion (E, Lz)*)
	For[j=1, j<=nUpdates, j++, 
	{
		Print["\t\t","update = ", j, " / ", nUpdates];
		
		(* ***initial E, Lz, Q*** *)
		e0 = e;
		Lz0 = Lz;
		Q0 = Q;
		
		(* ***integrate geodesics to compute orbital elements*** *)
		Print["\t\t\t Compute orbital elements"];
		
		(* ***initial positions/velocities*** *)
		t0GEO = t\[Tau][\[Tau]Array[[j]]];
		r0GEO = r\[Tau][\[Tau]Array[[j]]];
		\[Theta]0GEO = \[Theta]\[Tau][\[Tau]Array[[j]]];
		\[CurlyPhi]0GEO = \[CurlyPhi]\[Tau][\[Tau]Array[[j]]];	
		vt0GEO = vt\[Tau][\[Tau]Array[[j]]];
		vr0GEO = vr\[Tau][\[Tau]Array[[j]]];
		v\[Theta]0GEO = v\[Theta]\[Tau][\[Tau]Array[[j]]];
		v\[CurlyPhi]0GEO = v\[CurlyPhi]\[Tau][\[Tau]Array[[j]]];
		
		(* ***integration*** *)
		invariant = {energyEq[r[\[Tau]], \[Theta][\[Tau]], Derivative[1][t][\[Tau]], Derivative[1][\[CurlyPhi]][\[Tau]]], angularEq[r[\[Tau]], \[Theta][\[Tau]], Derivative[1][t][\[Tau]], Derivative[1][\[CurlyPhi]][\[Tau]]], carterEq[r[\[Tau]], \[Theta][\[Tau]], Derivative[1][\[Theta]][\[Tau]], Derivative[1][t][\[Tau]], Derivative[1][\[CurlyPhi]][\[Tau]]], constraintEq[r[\[Tau]], \[Theta][\[Tau]], Derivative[1][t][\[Tau]], Derivative[1][r][\[Tau]], Derivative[1][\[Theta]][\[Tau]], Derivative[1][\[CurlyPhi]][\[Tau]]]};
		ProjImpRK = {"Projection", Method->"ImplicitRungeKutta", "Invariants"->invariant};
		eqns\[Tau]GEO = {geot==0, geor==0, geo\[Theta]==0, geo\[CurlyPhi]==0, t[\[Tau]Array[[j]]]==t0GEO, Derivative[1][t][\[Tau]Array[[j]]]==vt0GEO, r[\[Tau]Array[[j]]]==r0GEO, Derivative[1][r][\[Tau]Array[[j]]]==vr0GEO, \[Theta][\[Tau]Array[[j]]]==\[Theta]0GEO, Derivative[1][\[Theta]][\[Tau]Array[[j]]]==v\[Theta]0GEO, \[CurlyPhi][\[Tau]Array[[j]]]==\[CurlyPhi]0GEO, Derivative[1][\[CurlyPhi]][\[Tau]Array[[j]]]==v\[CurlyPhi]0GEO};
		vars\[Tau]GEO = {t, r, \[Theta], \[CurlyPhi], Derivative[1][t], Derivative[1][r], Derivative[1][\[Theta]], Derivative[1][\[CurlyPhi]]};
		orbitTime\[Tau]GEO = NDSolve[eqns\[Tau]GEO, vars\[Tau]GEO, {\[Tau], \[Tau]Array[[j]], \[Tau]Array[[j]] + d\[Tau]Solve}, MaxSteps->\[Infinity], PrecisionGoal->MachinePrecision, Method->ProjImpRK];
		t\[Tau]GEO = t/.orbitTime\[Tau]GEO[[1, 1]];
		r\[Tau]GEO = r/.orbitTime\[Tau]GEO[[1, 2]];
		\[Theta]\[Tau]GEO = \[Theta]/.orbitTime\[Tau]GEO[[1, 3]];
		\[CurlyPhi]\[Tau]GEO = \[CurlyPhi]/.orbitTime\[Tau]GEO[[1, 4]];
		vt\[Tau]GEO = Derivative[1][t]/.orbitTime\[Tau]GEO[[1, 5]];
		vr\[Tau]GEO = Derivative[1][r]/.orbitTime\[Tau]GEO[[1, 6]];
		v\[Theta]\[Tau]GEO = Derivative[1][\[Theta]]/.orbitTime\[Tau]GEO[[1, 7]];
		v\[CurlyPhi]\[Tau]GEO = Derivative[1][\[CurlyPhi]]/.orbitTime\[Tau]GEO[[1, 8]];
	
		(* ***discretizing NDSolve solution*** *)
		(*f = vars\[Tau]GEO/.orbitTime\[Tau]GEO[[1]];
		data = Transpose[ToPackedArray[Parallelize[Table[f[[kindex]][Range[\[Tau]Array[[j]], \[Tau]Array[[j]] + d\[Tau]Solve, d\[Tau]Discr]],{kindex, 1, 8, 1}], ProgressReporting->False]]];*)
		data = ToPackedArray[Transpose[{r\[Tau]GEO[Range[\[Tau]Array[[j]], \[Tau]Array[[j]] + d\[Tau]Solve, d\[Tau]Discr]], \[Theta]\[Tau]GEO[Range[\[Tau]Array[[j]], \[Tau]Array[[j]] + d\[Tau]Solve, d\[Tau]Discr]]}]];
	
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
		
		(* ***integrate geodesics + GSF to compute trajectory*** *)
		Print["\t\t\t Compute trajectory"];
		eDot2PN = eDot[pOrbit, ecc, iOrbit];
		LzDot2PN = LzDot[pOrbit, ecc, iOrbit];
		QDot2PN = QDotMod[pOrbit, ecc, iOrbit];
		
		(* ***initial positions/velocities*** *)
		t0 = t\[Tau][\[Tau]Array[[j]]];
		r0 = r\[Tau][\[Tau]Array[[j]]];
		\[Theta]0 = \[Theta]\[Tau][\[Tau]Array[[j]]];
		\[CurlyPhi]0 = \[CurlyPhi]\[Tau][\[Tau]Array[[j]]];	
		vt0 = vt\[Tau][\[Tau]Array[[j]]];
		vr0 = vr\[Tau][\[Tau]Array[[j]]];
		v\[Theta]0 = v\[Theta]\[Tau][\[Tau]Array[[j]]];
		v\[CurlyPhi]0 = v\[CurlyPhi]\[Tau][\[Tau]Array[[j]]];
		
		(* ***integration*** *)
		Print["\t\t\t\t energy err = ", energyEq[r0, \[Theta]0, vt0, v\[CurlyPhi]0] - e0/\[Mu]];
		Print["\t\t\t\t angular err = ", angularEq[r0, \[Theta]0, vt0, v\[CurlyPhi]0] - Lz0/\[Mu]];
		Print["\t\t\t\t carter err = ", carterEq[r0, \[Theta]0, v\[Theta]0, vt0, v\[CurlyPhi]0] - Q0/\[Mu]^2];
		Print["\t\t\t\t 4-velocity err = ", constraintEq[r0, \[Theta]0, vt0, vr0, v\[Theta]0, v\[CurlyPhi]0]];
		eqns\[Tau] = {geotGSF==0, georGSF==0, geo\[Theta]GSF==0, geo\[CurlyPhi]GSF==0, t[\[Tau]Array[[j]]]==t0, Derivative[1][t][\[Tau]Array[[j]]]==vt0, r[\[Tau]Array[[j]]]==r0, Derivative[1][r][\[Tau]Array[[j]]]==vr0, \[Theta][\[Tau]Array[[j]]]==\[Theta]0, Derivative[1][\[Theta]][\[Tau]Array[[j]]]==v\[Theta]0, \[CurlyPhi][\[Tau]Array[[j]]]==\[CurlyPhi]0, Derivative[1][\[CurlyPhi]][\[Tau]Array[[j]]]==v\[CurlyPhi]0};
		vars\[Tau] = {t, r, \[Theta], \[CurlyPhi], Derivative[1][t], Derivative[1][r], Derivative[1][\[Theta]], Derivative[1][\[CurlyPhi]]};
		jMethod = 1;
		orbitTime\[Tau] = Quiet[Check[NDSolve[eqns\[Tau], vars\[Tau], {\[Tau], \[Tau]Array[[j]], \[Tau]Array[[j+1]]}, MaxSteps->\[Infinity], PrecisionGoal->MachinePrecision, Method->SWITCH[[jMethod]]], err, {NDSolve::ndsz, General::stop, InterpolatingFunction::dmval}]];
		While[orbitTime\[Tau] == err && jMethod <= 4,
		{
			jMethod += 1; 
			orbitTime\[Tau] = Quiet[Check[NDSolve[eqns\[Tau], vars\[Tau], {\[Tau], \[Tau]Array[[j]], \[Tau]Array[[j+1]]}, MaxSteps->\[Infinity], PrecisionGoal->MachinePrecision, Method->SWITCH[[jMethod]]], err, {NDSolve::ndsz, General::stop, InterpolatingFunction::dmval}]];
		}];
		If[orbitTime\[Tau] == err,
		{
			Print["\t\t\t", "NDSolve error -> try different integration method"];
			Break[]; 
		}];
		t\[Tau] = t/.orbitTime\[Tau][[1, 1]];
		r\[Tau] = r/.orbitTime\[Tau][[1, 2]];
		\[Theta]\[Tau] = \[Theta]/.orbitTime\[Tau][[1, 3]];
		\[CurlyPhi]\[Tau] = \[CurlyPhi]/.orbitTime\[Tau][[1, 4]];
		vt\[Tau] = Derivative[1][t]/.orbitTime\[Tau][[1, 5]];
		vr\[Tau] = Derivative[1][r]/.orbitTime\[Tau][[1, 6]];
		v\[Theta]\[Tau] = Derivative[1][\[Theta]]/.orbitTime\[Tau][[1, 7]];
		v\[CurlyPhi]\[Tau] = Derivative[1][\[CurlyPhi]]/.orbitTime\[Tau][[1, 8]];
		
		(* ***save energy and orbital data*** *)
		tData[j+nUpdatesTot] = t\[Tau][\[Tau]Array[[j]]];
		energyData[j+nUpdatesTot] = energyEq[r0, \[Theta]0, vt0, v\[CurlyPhi]0]*\[Mu]; 
		angularData[j+nUpdatesTot] = angularEq[r0, \[Theta]0, vt0, v\[CurlyPhi]0]*\[Mu]; 
		carterData[j+nUpdatesTot] = carterEq[r0, \[Theta]0, v\[Theta]0, vt0, v\[CurlyPhi]0]*\[Mu]^2;
		omegaRadialData[j+nUpdatesTot] = omegaRadialCoord; 
		omegaThetaData[j+nUpdatesTot] = omegaThetaCoord; 
		csiData[j+nUpdatesTot] = csiVal;
		pData[j+nUpdatesTot] = pOrbit;
		eccData[j+nUpdatesTot] = ecc;
		iData[j+nUpdatesTot] = iOrbit;
		
		(* ***fluxes*** *)
		dtFlux = t\[Tau][\[Tau]Array[[j+1]]] - t\[Tau][\[Tau]Array[[j]]];
		eFlux = eDot[pOrbit, ecc, iOrbit] * dtFlux;
		LzFlux = LzDot[pOrbit, ecc, iOrbit] * dtFlux;
		QFlux = QDotMod[pOrbit, ecc, iOrbit] * dtFlux;
		e = e0+eFlux; 
		Lz = Lz0+LzFlux;
		Q = Q0+QFlux;
	}];
	
	(* ***dissipation loop t > tResEnd*** *)
	Print["\t","-> Dissipation loop | t > tResEnd"];
	nUpdatesTot = nUpdatesTot + nUpdates;
	\[Tau]0 = \[Tau]Array[[-1]];
	nUpdates = IntegerPart[Round[(\[Tau]Max-\[Tau]0)/d\[Tau]Array3]]; (*number of iterations used to update integrals of motion (E, Lz)*)
	\[Tau]max = \[Tau]0 + nUpdates*d\[Tau]Array3;
	\[Tau]Array = Range[\[Tau]0, \[Tau]max, d\[Tau]Array3];
	For[j=1, j<=nUpdates, j++, 
	{
		Print["\t\t","update = ", j, " / ", nUpdates];
		
		(* ***initial E, Lz, Q*** *)
		e0 = e;
		Lz0 = Lz;
		Q0 = Q;
		
		(* ***integrate geodesics to compute orbital elements*** *)
		Print["\t\t\t Compute orbital elements"];
		
		(* ***initial positions/velocities*** *)
		t0GEO = t\[Tau][\[Tau]Array[[j]]];
		r0GEO = r\[Tau][\[Tau]Array[[j]]];
		\[Theta]0GEO = \[Theta]\[Tau][\[Tau]Array[[j]]];
		\[CurlyPhi]0GEO = \[CurlyPhi]\[Tau][\[Tau]Array[[j]]];	
		vt0GEO = vt\[Tau][\[Tau]Array[[j]]];
		vr0GEO = vr\[Tau][\[Tau]Array[[j]]];
		v\[Theta]0GEO = v\[Theta]\[Tau][\[Tau]Array[[j]]];
		v\[CurlyPhi]0GEO = v\[CurlyPhi]\[Tau][\[Tau]Array[[j]]];
		
		(* ***integration*** *)
		invariant = {energyEq[r[\[Tau]], \[Theta][\[Tau]], Derivative[1][t][\[Tau]], Derivative[1][\[CurlyPhi]][\[Tau]]], angularEq[r[\[Tau]], \[Theta][\[Tau]], Derivative[1][t][\[Tau]], Derivative[1][\[CurlyPhi]][\[Tau]]], carterEq[r[\[Tau]], \[Theta][\[Tau]], Derivative[1][\[Theta]][\[Tau]], Derivative[1][t][\[Tau]], Derivative[1][\[CurlyPhi]][\[Tau]]], constraintEq[r[\[Tau]], \[Theta][\[Tau]], Derivative[1][t][\[Tau]], Derivative[1][r][\[Tau]], Derivative[1][\[Theta]][\[Tau]], Derivative[1][\[CurlyPhi]][\[Tau]]]};
		ProjImpRK = {"Projection", Method->"ImplicitRungeKutta", "Invariants"->invariant};
		eqns\[Tau]GEO = {geot==0, geor==0, geo\[Theta]==0, geo\[CurlyPhi]==0, t[\[Tau]Array[[j]]]==t0GEO, Derivative[1][t][\[Tau]Array[[j]]]==vt0GEO, r[\[Tau]Array[[j]]]==r0GEO, Derivative[1][r][\[Tau]Array[[j]]]==vr0GEO, \[Theta][\[Tau]Array[[j]]]==\[Theta]0GEO, Derivative[1][\[Theta]][\[Tau]Array[[j]]]==v\[Theta]0GEO, \[CurlyPhi][\[Tau]Array[[j]]]==\[CurlyPhi]0GEO, Derivative[1][\[CurlyPhi]][\[Tau]Array[[j]]]==v\[CurlyPhi]0GEO};
		vars\[Tau]GEO = {t, r, \[Theta], \[CurlyPhi], Derivative[1][t], Derivative[1][r], Derivative[1][\[Theta]], Derivative[1][\[CurlyPhi]]};
		orbitTime\[Tau]GEO = NDSolve[eqns\[Tau]GEO, vars\[Tau]GEO, {\[Tau], \[Tau]Array[[j]], \[Tau]Array[[j]] + d\[Tau]Solve}, MaxSteps->\[Infinity], PrecisionGoal->MachinePrecision, Method->ProjImpRK];
		t\[Tau]GEO = t/.orbitTime\[Tau]GEO[[1, 1]];
		r\[Tau]GEO = r/.orbitTime\[Tau]GEO[[1, 2]];
		\[Theta]\[Tau]GEO = \[Theta]/.orbitTime\[Tau]GEO[[1, 3]];
		\[CurlyPhi]\[Tau]GEO = \[CurlyPhi]/.orbitTime\[Tau]GEO[[1, 4]];
		vt\[Tau]GEO = Derivative[1][t]/.orbitTime\[Tau]GEO[[1, 5]];
		vr\[Tau]GEO = Derivative[1][r]/.orbitTime\[Tau]GEO[[1, 6]];
		v\[Theta]\[Tau]GEO = Derivative[1][\[Theta]]/.orbitTime\[Tau]GEO[[1, 7]];
		v\[CurlyPhi]\[Tau]GEO = Derivative[1][\[CurlyPhi]]/.orbitTime\[Tau]GEO[[1, 8]];
	
		(* ***discretizing NDSolve solution*** *)
		(*f = vars\[Tau]GEO/.orbitTime\[Tau]GEO[[1]];
		data = Transpose[ToPackedArray[Parallelize[Table[f[[kindex]][Range[\[Tau]Array[[j]], \[Tau]Array[[j]] + d\[Tau]Solve, d\[Tau]Discr]],{kindex, 1, 8, 1}], ProgressReporting->False]]];*)
		data = ToPackedArray[Transpose[{r\[Tau]GEO[Range[\[Tau]Array[[j]], \[Tau]Array[[j]] + d\[Tau]Solve, d\[Tau]Discr]], \[Theta]\[Tau]GEO[Range[\[Tau]Array[[j]], \[Tau]Array[[j]] + d\[Tau]Solve, d\[Tau]Discr]]}]];
	
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
		
		(* ***integrate geodesics + GSF to compute trajectory*** *)
		Print["\t\t\t Compute trajectory"];
		eDot2PN = eDot[pOrbit, ecc, iOrbit];
		LzDot2PN = LzDot[pOrbit, ecc, iOrbit];
		QDot2PN = QDotMod[pOrbit, ecc, iOrbit];
		
		(* ***initial positions/velocities*** *)
		t0 = t\[Tau][\[Tau]Array[[j]]];
		r0 = r\[Tau][\[Tau]Array[[j]]];
		\[Theta]0 = \[Theta]\[Tau][\[Tau]Array[[j]]];
		\[CurlyPhi]0 = \[CurlyPhi]\[Tau][\[Tau]Array[[j]]];	
		vt0 = vt\[Tau][\[Tau]Array[[j]]];
		vr0 = vr\[Tau][\[Tau]Array[[j]]];
		v\[Theta]0 = v\[Theta]\[Tau][\[Tau]Array[[j]]];
		v\[CurlyPhi]0 = v\[CurlyPhi]\[Tau][\[Tau]Array[[j]]];
		
		(* ***integration*** *)
		Print["\t\t\t\t energy err = ", energyEq[r0, \[Theta]0, vt0, v\[CurlyPhi]0] - e0/\[Mu]];
		Print["\t\t\t\t angular err = ", angularEq[r0, \[Theta]0, vt0, v\[CurlyPhi]0] - Lz0/\[Mu]];
		Print["\t\t\t\t carter err = ", carterEq[r0, \[Theta]0, v\[Theta]0, vt0, v\[CurlyPhi]0] - Q0/\[Mu]^2];
		Print["\t\t\t\t 4-velocity err = ", constraintEq[r0, \[Theta]0, vt0, vr0, v\[Theta]0, v\[CurlyPhi]0]];
		eqns\[Tau] = {geotGSF==0, georGSF==0, geo\[Theta]GSF==0, geo\[CurlyPhi]GSF==0, t[\[Tau]Array[[j]]]==t0, Derivative[1][t][\[Tau]Array[[j]]]==vt0, r[\[Tau]Array[[j]]]==r0, Derivative[1][r][\[Tau]Array[[j]]]==vr0, \[Theta][\[Tau]Array[[j]]]==\[Theta]0, Derivative[1][\[Theta]][\[Tau]Array[[j]]]==v\[Theta]0, \[CurlyPhi][\[Tau]Array[[j]]]==\[CurlyPhi]0, Derivative[1][\[CurlyPhi]][\[Tau]Array[[j]]]==v\[CurlyPhi]0};
		vars\[Tau] = {t, r, \[Theta], \[CurlyPhi], Derivative[1][t], Derivative[1][r], Derivative[1][\[Theta]], Derivative[1][\[CurlyPhi]]};
		jMethod = 1;
		orbitTime\[Tau] = Quiet[Check[NDSolve[eqns\[Tau], vars\[Tau], {\[Tau], \[Tau]Array[[j]], \[Tau]Array[[j+1]]}, MaxSteps->\[Infinity], PrecisionGoal->MachinePrecision, Method->SWITCH[[jMethod]]], err, {NDSolve::ndsz, General::stop, InterpolatingFunction::dmval}]];
		While[orbitTime\[Tau] == err && jMethod <= 4,
		{
			jMethod += 1; 
			orbitTime\[Tau] = Quiet[Check[NDSolve[eqns\[Tau], vars\[Tau], {\[Tau], \[Tau]Array[[j]], \[Tau]Array[[j+1]]}, MaxSteps->\[Infinity], PrecisionGoal->MachinePrecision, Method->SWITCH[[jMethod]]], err, {NDSolve::ndsz, General::stop, InterpolatingFunction::dmval}]];
		}];
		If[orbitTime\[Tau] == err,
		{
			Print["\t\t\t", "NDSolve error -> try different integration method"];
			Break[]; 
		}];
		t\[Tau] = t/.orbitTime\[Tau][[1, 1]];
		r\[Tau] = r/.orbitTime\[Tau][[1, 2]];
		\[Theta]\[Tau] = \[Theta]/.orbitTime\[Tau][[1, 3]];
		\[CurlyPhi]\[Tau] = \[CurlyPhi]/.orbitTime\[Tau][[1, 4]];
		vt\[Tau] = Derivative[1][t]/.orbitTime\[Tau][[1, 5]];
		vr\[Tau] = Derivative[1][r]/.orbitTime\[Tau][[1, 6]];
		v\[Theta]\[Tau] = Derivative[1][\[Theta]]/.orbitTime\[Tau][[1, 7]];
		v\[CurlyPhi]\[Tau] = Derivative[1][\[CurlyPhi]]/.orbitTime\[Tau][[1, 8]];
		
		(* ***save energy and orbital data*** *)
		tData[j+nUpdatesTot] = t\[Tau][\[Tau]Array[[j]]];
		energyData[j+nUpdatesTot] = energyEq[r0, \[Theta]0, vt0, v\[CurlyPhi]0]*\[Mu]; 
		angularData[j+nUpdatesTot] = angularEq[r0, \[Theta]0, vt0, v\[CurlyPhi]0]*\[Mu]; 
		carterData[j+nUpdatesTot] = carterEq[r0, \[Theta]0, v\[Theta]0, vt0, v\[CurlyPhi]0]*\[Mu]^2;
		omegaRadialData[j+nUpdatesTot] = omegaRadialCoord; 
		omegaThetaData[j+nUpdatesTot] = omegaThetaCoord; 
		csiData[j+nUpdatesTot] = csiVal;
		pData[j+nUpdatesTot] = pOrbit;
		eccData[j+nUpdatesTot] = ecc;
		iData[j+nUpdatesTot] = iOrbit;
		
		(* ***fluxes*** *)
		dtFlux = t\[Tau][\[Tau]Array[[j+1]]] - t\[Tau][\[Tau]Array[[j]]];
		eFlux = eDot[pOrbit, ecc, iOrbit] * dtFlux;
		LzFlux = LzDot[pOrbit, ecc, iOrbit] * dtFlux;
		QFlux = QDotMod[pOrbit, ecc, iOrbit] * dtFlux;
		e = e0+eFlux; 
		Lz = Lz0+LzFlux;
		Q = Q0+QFlux;
	}];
	
	(* ***export data*** *)
	nUpdatesTot = nUpdatesTot + (j-1);
	data2 = Table[{tData[i], energyData[i], angularData[i], carterData[i], omegaRadialData[i], omegaThetaData[i], csiData[i], pData[i], eccData[i], iData[i]}, {i, 1, nUpdatesTot}];
	fname2 = StringJoin["/THEORY/USERS/edoardo.levati/DATA/kick/kerr/energyData_adiabatic_TEST_EMRI1_Coupled_Method_2/energy_", ToString[1], ".h5"];
	Export[fname2, data2, OverwriteTarget -> True];
}];
(*CloseKernels[];*)
Exit[];
