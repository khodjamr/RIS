(* ::Package:: *)

(* Mathematica 9 compatible numerical package for the manuscript
   "Null-Space-Assisted Inverse Synthesis of Dual-Mode Conformal
   Metasurfaces With Shared Impedance States". *)

BeginPackage["DualModeRIS`"];

SelfTest::usage =
  "SelfTest[] checks the integral-equation convention against the analytic PEC-cylinder series, verifies a Sherman-Morrison update, and tests null-space invariance.";
RunStudy::usage =
  "RunStudy[] reproduces the reported configuration, inverse-source, scattering, telemetry, mesh, frequency, and uncertainty studies; exports descriptive CSV and PDF files; and displays the five manuscript figures.";
StudyDirectory::usage = "StudyDirectory is an option for RunStudy specifying the output root.";
MonteCarloTrials::usage = "MonteCarloTrials is an option for RunStudy specifying the number of impedance-error trials.";
RefineCodebooks::usage = "RefineCodebooks is an option for RunStudy; True performs deterministic coordinate polishing before evaluation.";
DisplayFigures::usage = "DisplayFigures is an option for RunStudy; True prints the five manuscript figures in the launcher notebook.";

Options[RunStudy] = {StudyDirectory -> Automatic, MonteCarloTrials -> 200,
  RefineCodebooks -> False, DisplayFigures -> True};

Begin["`Private`"];

$PackageDirectory = DirectoryName[$InputFileName];

value[m_, key_String] := key /. m;
wrapAngle[x_] := Mod[x + Pi, 2 Pi] - Pi;
db10[x_] := 10 Log[10, Max[N[x], 10^-300]];

buildModel[n_Integer, a_?NumericQ, fscale_?NumericQ, nobs_Integer] :=
 Module[{k, eta, phi, r, ds, theta, self, s, amat},
  k = 2 Pi fscale; eta = 1.;
  phi = N[Range[0, n - 1] 2 Pi/n];
  r = Transpose[{a Cos[phi], a Sin[phi]}];
  ds = 2 Pi a/n;
  self = ds (1 + (2 I/Pi) (Log[k ds/4] + EulerGamma - 1));
  s = Table[
    If[p == q, (k eta/4) self,
      (k eta/4) HankelH1[0, k Norm[r[[p]] - r[[q]]]] ds],
    {p, n}, {q, n}];
  theta = N[Range[0, nobs - 1] 2 Pi/nobs - Pi];
  amat = Table[
    ds Exp[-I k ({Cos[theta[[l]]], Sin[theta[[l]]]}.r[[q]])],
    {l, nobs}, {q, n}];
  {"N" -> n, "a" -> a, "k" -> k, "eta" -> eta,
   "phi" -> phi, "r" -> r, "ds" -> ds, "theta" -> theta,
   "S" -> s, "A" -> amat}
  ];

planeWave[m_, alpha_?NumericQ] := Module[{d},
  d = {Cos[alpha], Sin[alpha]};
  Exp[I value[m, "k"] (value[m, "r"].d)]
  ];

lineFeed[m_, feed_List] := Module[{k, eta, r},
  k = value[m, "k"]; eta = value[m, "eta"]; r = value[m, "r"];
  Table[-k eta/4 HankelH1[0, k Norm[r[[p]] - feed]],
    {p, value[m, "N"]}]
  ];

solveCurrent[m_, z_List, b_List] :=
  LinearSolve[DiagonalMatrix[z] + value[m, "S"], b];

farField[m_, j_List, feed_: None] := Module[{f, k, eta, theta, khat},
  k = value[m, "k"]; eta = value[m, "eta"];
  f = -k eta/4 (value[m, "A"].j);
  If[feed =!= None,
   theta = value[m, "theta"];
   khat = Transpose[{Cos[theta], Sin[theta]}];
   f = f - k eta/4 Exp[-I k (khat.feed)]];
  f
  ];

farMatrix[m_, angles_List] := Module[{k, eta, ds, r, khat},
  k = value[m, "k"]; eta = value[m, "eta"];
  ds = value[m, "ds"]; r = value[m, "r"];
  khat = Transpose[{Cos[angles], Sin[angles]}];
  -k eta/4 ds Table[
    Exp[-I k (khat[[p]].r[[q]])],
    {p, Length[angles]}, {q, value[m, "N"]}]
  ];

expandCode[q_List, mult_Integer] := Flatten[(Table[#, {mult}] &) /@ q];
stateProfile[states_List, q_List, mult_Integer] := states[[expandCode[q, mult]]];

protectedIndices[m_, alpha_?NumericQ, half_?NumericQ] := Module[{theta},
  theta = value[m, "theta"];
  Select[Range[Length[theta]],
    Abs[wrapAngle[theta[[#]] - (alpha + Pi)]] <= half + 100 $MachineEpsilon &]
  ];

protectedRatios[m_, states_List, q_List, mult_Integer,
   incdeg_List : Range[-30, 30, 5], halfdeg_: 10] :=
 Module[{z, alpha, b, f, f0, ids},
  z = stateProfile[states, q, mult];
  Table[
   alpha = N[deg Degree]; b = planeWave[m, alpha];
   f = farField[m, solveCurrent[m, z, b]];
   f0 = farField[m, solveCurrent[m,
       ConstantArray[0. + 0. I, value[m, "N"]], b]];
   ids = protectedIndices[m, alpha, N[halfdeg Degree]];
   Max[Abs[f[[ids]]]^2]/Max[Abs[f0[[ids]]]^2],
   {deg, incdeg}]
  ];

telemetryMetrics[m_, states_List, q_List, mult_Integer] :=
 Module[{feed, z, f, p, theta, main, side, eff, sll, direct, imax,
   edgeCorrection},
  feed = {value[m, "a"] - .55, 0.};
  z = stateProfile[states, q, mult];
  f = farField[m, solveCurrent[m, z, lineFeed[m, feed]], feed];
  p = Abs[f]^2; theta = value[m, "theta"];
  main = Select[Range[Length[theta]],
    Abs[wrapAngle[theta[[#]]]] <= 8 Degree + 100 $MachineEpsilon &];
  side = Select[Range[Length[theta]],
    Abs[wrapAngle[theta[[#]]]] >= 15 Degree - 100 $MachineEpsilon &];
  (* Uniform angular quadrature: the full circle is periodic, while the
     two main-sector end samples receive trapezoidal half weights. *)
  edgeCorrection = (p[[First[main]]] + p[[Last[main]]])/2.;
  eff = (Total[p[[main]]] - edgeCorrection)/Total[p];
  sll = Max[p[[side]]]/Max[p[[main]]];
  direct = Max[p]/Mean[p]; imax = First[First[Position[p, Max[p]]]];
  {"Efficiency" -> eff, "SLL" -> sll, "Directivity" -> direct,
   "BeamAngle" -> theta[[imax]]/Degree, "Field" -> f}
  ];

publishedStates[] := .03 + I {-6.75602121, -1.61813346, -.56199291,
    -.22403657, .18768728, .60721582, 1.34454867, 4.93228419};

publishedQS[] := {7, 6, 7, 7, 8, 5, 4, 1, 8, 4, 7, 4, 6, 5, 6, 7,
  1, 5, 5, 5, 1, 4, 1, 4, 8, 3, 5, 6, 4, 1, 6, 5,
  6, 3, 7, 4, 1, 4, 1, 2, 1, 2, 1, 3, 1, 5, 6, 5,
  5, 4, 7, 4, 5, 4, 5, 3, 2, 5, 6, 6, 5, 8, 4, 7};

publishedQT[] := {4, 4, 1, 4, 5, 4, 6, 4, 5, 7, 3, 5, 4, 4, 4, 4,
  3, 5, 5, 4, 4, 5, 7, 5, 4, 4, 3, 3, 5, 7, 1, 1,
  3, 1, 6, 5, 3, 4, 4, 2, 5, 3, 6, 5, 4, 4, 5, 4,
  4, 5, 4, 4, 6, 6, 3, 6, 6, 2, 7, 6, 5, 4, 8, 4};

stealthObjective[m_, states_List, q_List, mult_Integer] := Module[{r},
  r = protectedRatios[m, states, q, mult];
  {.45 Mean[r] + .55 Max[r], r}
  ];

telemetryObjective[m_, states_List, q_List, mult_Integer] :=
 Module[{t, v},
  t = telemetryMetrics[m, states, q, mult];
  v = (1 - ("Efficiency" /. t)) + .15 ("SLL" /. t) +
    .002 ("BeamAngle" /. t)^2;
  {v, t}
  ];

coordinatePolish[m_, states_List, q0_List, mult_Integer, which_String,
   sweeps_Integer : 3] := Module[
  {q, fn, current, metric, order, old, best, bestk, trial, test,
   k, sw, changed},
  q = q0;
  fn = If[which == "S", stealthObjective, telemetryObjective];
  {current, metric} = fn[m, states, q, mult];
  Do[
   changed = 0; SeedRandom[1000 + sw]; order = RandomSample[Range[Length[q]]];
   Do[
    old = q[[idx]]; best = current; bestk = old;
    Do[If[k =!= old,
      trial = ReplacePart[q, idx -> k]; test = First[fn[m, states, trial, mult]];
      If[test < best - 10^-11, best = test; bestk = k]],
     {k, Length[states]}];
    If[bestk =!= old, q[[idx]] = bestk; current = best; changed++],
    {idx, order}];
   If[changed == 0, Break[]], {sw, sweeps}];
  {q, fn[m, states, q, mult]}
  ];

conditionedNullProjector[c_?MatrixQ] := Module[
  {rows, columns, work, gram, solved, projector, residual,
   idempotence},
  {rows, columns} = Dimensions[c];
  If[rows >= columns,
   Print["Null-projector construction failed: operator dimensions are ",
    {rows, columns}, "; fewer rows than columns are required."];
   Return[$Failed]];
  (* For a full-row-rank restricted radiation matrix C, the orthogonal
     projector onto Null[C] is I-C^H (C C^H)^(-1) C.  The sampled operators
     are very ill-conditioned, so every operation that forms this projector
     is carried out at 80-digit precision.  This avoids Mathematica 9 rank
     decisions in NullSpace and QRDecomposition while preserving every
     prescribed far-field constraint. *)
  work = SetPrecision[c, 80];
  gram = work.ConjugateTranspose[work];
  solved = Quiet[LinearSolve[gram, work]];
  If[! MatrixQ[solved, NumericQ],
   Print["Null-projector construction failed: the 80-digit Gram solve did ",
    "not return a numerical matrix for operator dimensions ",
    {rows, columns}, "."];
   Return[$Failed]];
  projector = IdentityMatrix[columns] -
    ConjugateTranspose[work].solved;
  projector = N[(projector + ConjugateTranspose[projector])/2];
  If[Dimensions[projector] =!= {columns, columns},
   Print["Null-projector construction failed: unexpected projector ",
    "dimensions ", Dimensions[projector], "."];
   Return[$Failed]];
  residual = Norm[N[c].projector]/
    Max[Norm[N[c]] Norm[projector], 10^-30];
  idempotence = Norm[projector.projector - projector]/
    Max[Norm[projector], 10^-30];
  If[! NumericQ[residual] || ! NumericQ[idempotence],
   Print["Null-projector construction failed: nonnumerical verification ",
    "metric for operator dimensions ", {rows, columns}, "."];
   Return[$Failed]];
  If[residual >= 10^-8 || idempotence >= 10^-8,
   Print["Null-projector construction failed verification for operator ",
    "dimensions ", {rows, columns}, ": normalized residual = ",
    ScientificForm[residual, 4], ", idempotence error = ",
    ScientificForm[idempotence, 4], "."];
   Return[$Failed]];
  projector
  ];

projectorNullBasis[projector_?MatrixQ] := Module[
  {hermitian, eigenvalues, eigenvectors, ids, nmat, basisResidual,
   orthogonalityError, columns},
  columns = First[Dimensions[projector]];
  hermitian = SetPrecision[
    N[(projector + ConjugateTranspose[projector])/2], MachinePrecision];
  {eigenvalues, eigenvectors} = Quiet[Eigensystem[hermitian]];
  If[! VectorQ[eigenvalues, NumericQ] ||
    ! MatrixQ[eigenvectors, NumericQ],
   Print["Null-basis extraction failed: Eigensystem did not return ",
    "numerical values and vectors."];
   Return[$Failed]];
  ids = Select[Range[Length[eigenvalues]],
    Re[eigenvalues[[#]]] > .5 &];
  If[Length[ids] == 0,
   Print["Null-basis extraction failed: no unit eigenvalues were found."];
   Return[$Failed]];
  nmat = Transpose[eigenvectors[[ids]]];
  If[Dimensions[nmat][[1]] =!= columns,
   Print["Null-basis extraction failed: unexpected basis dimensions ",
    Dimensions[nmat], "."];
   Return[$Failed]];
  basisResidual = Norm[(IdentityMatrix[columns] - projector).nmat]/
    Max[Norm[nmat], 10^-30];
  orthogonalityError = Norm[
     ConjugateTranspose[nmat].nmat -
      IdentityMatrix[Dimensions[nmat][[2]]]]/
    Max[1., Sqrt[Dimensions[nmat][[2]]]];
  If[! NumericQ[basisResidual] || ! NumericQ[orthogonalityError] ||
    basisResidual >= 10^-8 || orthogonalityError >= 10^-8,
   Print["Null-basis extraction failed verification: projector residual = ",
    ScientificForm[basisResidual, 4], ", orthogonality error = ",
    ScientificForm[orthogonalityError, 4], "."];
   Return[$Failed]];
  nmat
  ];

particularAndNull[c_, rhs_, jref_, relreg_] := Module[
  {sv, lambda, jp, projector},
  sv = SingularValueList[c]; lambda = relreg First[sv]^2;
  jp = jref + ConjugateTranspose[c].LinearSolve[
      c.ConjugateTranspose[c] + lambda IdentityMatrix[Length[rhs]],
      rhs - c.jref];
  projector = conditionedNullProjector[c];
  {jp, projector, sv, lambda}
  ];

retrieveImpedance[m_, b_, j_] := (b - value[m, "S"].j)/j;

assignShared[zS_List, zT_List, states_List] := Module[{qS, qT, costs},
  qS = Table[
    costs = Table[Total[Abs[zS[[All, n]] - states[[k]]]^2],
      {k, Length[states]}]; First[Ordering[costs, 1]],
    {n, Length[zT]}];
  qT = Table[
    First[Ordering[Abs[zT[[n]] - states]^2, 1]], {n, Length[zT]}];
  {qS, qT}
  ];

updateStates[zS_List, zT_List, qS_List, qT_List, states_List] :=
 Module[{out, vals, x, k, r},
  out = states;
  Do[
   vals = Flatten[Join[
       Table[Pick[zS[[r]], qS, k], {r, Length[zS]}],
       {Pick[zT, qT, k]}]];
   If[Length[vals] > 0,
    x = Min[6., Max[-6., Mean[Im[vals]]]]; out[[k]] = .03 + I x],
   {k, Length[states]}];
  SortBy[out, Im]
  ];

libraryMismatch[zS_List, zT_List, qS_List, qT_List, states_List] :=
 (Total[Table[Total[Abs[zS[[r]] - states[[qS]]]^2],
      {r, Length[zS]}]] + Total[Abs[zT - states[[qT]]]^2])/
   ((Length[zS] + 1) Length[zT]);

nullStep[m_, b_, j_, nmat_, target_, ridge_: 10^-4] := Module[
  {z, linearizedMap, bmat, nullity, scale, normalMatrix, normalRHS,
   cc, delta, e0, candidate, et, bestEt, relativeDelta, accepted,
   acceptedCurrent, acceptedEt, acceptedRelativeDelta, t},
  If[nmat === $Failed, Return[{j, False, Infinity, Infinity, 0.}]];
  z = retrieveImpedance[m, b, j];
  linearizedMap = -DiagonalMatrix[1/j].
    (value[m, "S"] + DiagonalMatrix[z]);
  bmat = linearizedMap.nmat;
  nullity = Dimensions[nmat][[2]];
  (* Total[matrix] is a vector in Mathematica.  Flatten is required for the
     scalar Frobenius-norm scale used in the Tikhonov term. *)
  scale = Total[Flatten[Abs[bmat]^2]]/nullity;
  normalMatrix = ConjugateTranspose[bmat].bmat +
    ridge scale IdentityMatrix[nullity];
  normalRHS = ConjugateTranspose[bmat].(target - z);
  cc = Quiet[LinearSolve[normalMatrix, normalRHS]];
  e0 = Mean[Abs[z - target]^2];
  If[! VectorQ[cc, NumericQ],
   Return[{j, False, e0, Infinity, 0.}]];
  delta = nmat.cc;
  relativeDelta = Norm[delta]/Max[Norm[j], 10^-30];
  bestEt = Infinity; accepted = False; acceptedCurrent = j;
  acceptedEt = Infinity; acceptedRelativeDelta = 0.;
  (* In Mathematica 9, Return inside Do exits the loop expression rather
     than reliably returning from this Module.  Store the accepted trial
     explicitly and leave the loop with Break. *)
  Do[
   candidate = j + t delta;
   If[Min[Abs[candidate]] > 10^-5 Max[Abs[candidate]],
    et = Mean[Abs[retrieveImpedance[m, b, candidate] - target]^2];
    bestEt = Min[bestEt, et];
    If[et < e0,
     accepted = True; acceptedCurrent = candidate; acceptedEt = et;
     acceptedRelativeDelta = t relativeDelta; Break[]]],
   {t, {1., .5, .25, .125, .0625}}];
  If[TrueQ[accepted],
   {acceptedCurrent, True, e0, acceptedEt, acceptedRelativeDelta},
   {j, False, e0, bestEt, relativeDelta}]
  ];

(* The pre-audit inverse implementation was removed after verification. *)
(* Corrected inverse diagnostic.  The original record-replacement expression
   did not update the stored current vectors in Mathematica 9.  This version
   stores the currents explicitly and measures invariance of C.j itself. *)
inverseAblation[niter_Integer : 12, verbose_ : True] := Module[
  {m, inc, stealth, alpha, b, jref, angles, c, p, feed, bT, zref,
   angT, cT, khat, ffeed, fref, g, a0, desired, bS, jS, nS, cS,
   jT, nT, states, zS, zT, qS, qT, mismatch, hist, r, data0S,
   data0T, dataChange, allProjectors, projectorNullities,
   projectorResidual, allBases, basisS, basisT, basisNullities,
   stepS, stepT, acceptedSteps, inverseFailure},
  m = buildModel[64, 2., 1., 360]; inc = {-30, 0, 30} Degree;
  stealth = Table[
    b = planeWave[m, alpha];
    jref = solveCurrent[m, ConstantArray[0. + 0. I, 64], b];
    angles = alpha + Pi + Range[-12, 12, 1.5] Degree;
    c = farMatrix[m, angles];
    p = particularAndNull[c, ConstantArray[0. + 0. I, Length[angles]],
      jref, .003];
    {b, p[[1]], p[[2]], c}, {alpha, inc}];
  bS = stealth[[All, 1]]; jS = stealth[[All, 2]];
  nS = stealth[[All, 3]]; cS = stealth[[All, 4]];

  feed = {value[m, "a"] - .55, 0.}; bT = lineFeed[m, feed];
  zref = ConstantArray[.03 + .4 I, 64]; jref = solveCurrent[m, zref, bT];
  angT = Join[Range[-170, -20, 10], Range[-12, 12, 1.5],
     Range[20, 170, 10]] Degree;
  cT = farMatrix[m, angT];
  khat = Transpose[{Cos[angT], Sin[angT]}];
  ffeed = -value[m, "k"] value[m, "eta"]/4 Exp[
     -I value[m, "k"] (khat.feed)];
  fref = cT.jref + ffeed; g = Exp[-.5 (wrapAngle[angT]/(9 Degree))^2];
  a0 = 2 Max[Abs[fref]] Exp[I Arg[fref[[First[Ordering[Abs[angT], 1]]]]]];
  desired = a0 g;
  p = particularAndNull[cT, desired - ffeed, jref, 10^-5];
  jT = p[[1]]; nT = p[[2]];

  allProjectors = Join[nS, {nT}];
  If[MemberQ[allProjectors, $Failed],
   Print["ERROR: high-precision null-projector construction failed."];
   Return[$Failed]];
  projectorNullities = Round[Re[Tr[#]]] & /@ allProjectors;
  projectorResidual = Max[Join[
     Table[Norm[cS[[r]].nS[[r]]]/Max[Norm[cS[[r]]], 10^-30],
      {r, Length[cS]}],
     {Norm[cT.nT]/Max[Norm[cT], 10^-30]}]];
  If[TrueQ[verbose],
   Print["RunStudy: inverse null-space dimensions = ",
    projectorNullities, "; maximum projector residual = ",
    ScientificForm[projectorResidual, 4], "."]];
  If[projectorNullities =!= {47, 47, 47, 15} ||
    projectorResidual >= 10^-8,
   Print["ERROR: inverse null-projector verification failed."];
   Return[$Failed]];
  allBases = projectorNullBasis /@ allProjectors;
  If[MemberQ[allBases, $Failed],
   Print["ERROR: reduced null-basis construction failed."];
   Return[$Failed]];
  basisNullities = Dimensions[#][[2]] & /@ allBases;
  If[basisNullities =!= projectorNullities,
   Print["ERROR: reduced null-basis dimensions ", basisNullities,
    " do not match the projector dimensions ", projectorNullities, "."];
   Return[$Failed]];
  basisS = Take[allBases, Length[nS]]; basisT = Last[allBases];
  If[TrueQ[verbose],
   Print["RunStudy: reduced null-space bases verified with dimensions = ",
    basisNullities, "."]];

  data0S = Table[cS[[r]].jS[[r]], {r, Length[jS]}];
  data0T = cT.jT;
  states = SortBy[.03 + I Cot[(2 Range[0, 3] + 1) Pi/(2 4)], Im];
  zS = Table[retrieveImpedance[m, bS[[r]], jS[[r]]], {r, Length[jS]}];
  zT = retrieveImpedance[m, bT, jT];
  {qS, qT} = assignShared[zS, zT, states];
  states = updateStates[zS, zT, qS, qT, states];
  {qS, qT} = assignShared[zS, zT, states];
  mismatch = libraryMismatch[zS, zT, qS, qT, states];
  hist = {{0., mismatch, 0.}}; inverseFailure = False;

  Do[
   stepS = Table[
     nullStep[m, bS[[r]], jS[[r]], basisS[[r]], states[[qS]]],
     {r, Length[jS]}];
   stepT = nullStep[m, bT, jT, basisT, states[[qT]]];
   acceptedSteps = Count[Join[stepS[[All, 2]], {stepT[[2]]}], True];
   If[acceptedSteps == 0,
    Print["ERROR: reduced null-space line search accepted no update at ",
     "iteration ", it, ". Initial nonlinear errors = ",
     Join[stepS[[All, 3]], {stepT[[3]]}],
     "; best trial errors = ", Join[stepS[[All, 4]], {stepT[[4]]}],
     "; proposed relative step norms = ",
     Join[stepS[[All, 5]], {stepT[[5]]}], "."];
    inverseFailure = True; Break[]];
   jS = stepS[[All, 1]]; jT = stepT[[1]];
   zS = Table[retrieveImpedance[m, bS[[r]], jS[[r]]], {r, Length[jS]}];
   zT = retrieveImpedance[m, bT, jT];
   {qS, qT} = assignShared[zS, zT, states];
   states = updateStates[zS, zT, qS, qT, states];
   {qS, qT} = assignShared[zS, zT, states];
   mismatch = libraryMismatch[zS, zT, qS, qT, states];
   dataChange = Max[Join[
      Table[Norm[cS[[r]].jS[[r]] - data0S[[r]]]/
        Max[Norm[data0S[[r]]], 10^-30], {r, Length[jS]}],
      {Norm[cT.jT - data0T]/Max[Norm[data0T], 10^-30]}]];
   AppendTo[hist, {N[it], mismatch, dataChange}], {it, 1, niter}];
  If[TrueQ[inverseFailure], $Failed, hist]
  ];

analyticPEC[m_, rho_?NumericQ, alpha_?NumericQ, nmax_Integer : 80] :=
 Module[{k, a, theta},
  k = value[m, "k"]; a = value[m, "a"]; theta = value[m, "theta"];
  Table[Sum[-I^n BesselJ[n, k a]/HankelH1[n, k a]
       HankelH1[n, k rho] Exp[I n (th - alpha)],
     {n, -nmax, nmax}], {th, theta}]
  ];

numericNearField[m_, j_List, rho_?NumericQ] := Module[{k, eta, ds, r, theta, obs},
  k = value[m, "k"]; eta = value[m, "eta"];
  ds = value[m, "ds"]; r = value[m, "r"]; theta = value[m, "theta"];
  obs = Transpose[{rho Cos[theta], rho Sin[theta]}];
  Table[-k eta/4 ds Sum[
      HankelH1[0, k Norm[obs[[p]] - r[[q]]]] j[[q]],
      {q, value[m, "N"]}], {p, Length[theta]}]
  ];

SelfTest[] := Module[
  {m, b, j, rho, numeric, exact, pecerr, z, mat, inv, idx, delta,
   direct, updated, smerr, inverseTestModel, operators, projectors,
   probe, nullDelta, nserr, projectorError, projectorRanks,
   nullitypass, dimensionpass, bases, basisRanks, basispass,
   inverseSmoke, inversepass, inverseInitial, inverseFinal,
   inverseDataChange, pass, k, p},
  Print["DualModeRIS SelfTest: starting..."];
  m = buildModel[256, 2., 1., 72]; b = planeWave[m, 0.];
  j = solveCurrent[m, ConstantArray[0. + 0. I, 256], b];
  rho = 5 value[m, "a"]; numeric = numericNearField[m, j, rho];
  exact = analyticPEC[m, rho, 0., 80];
  pecerr = Norm[numeric - exact]/Norm[exact];
  z = ConstantArray[.03 + .4 I, 256]; mat = DiagonalMatrix[z] + value[m, "S"];
  inv = Inverse[mat]; idx = 37; delta = .1 - .07 I;
  direct = LinearSolve[mat + delta SparseArray[{{idx, idx} -> 1.}, {256, 256}], b];
  updated = inv.b - delta inv[[All, idx]] (inv.b)[[idx]]/
      (1 + delta inv[[idx, idx]]);
  smerr = Norm[direct - updated]/Norm[direct];
  inverseTestModel = buildModel[64, 2., 1., 72];
  operators = {
    farMatrix[inverseTestModel, Range[-12, 12, 1.5] Degree],
    farMatrix[inverseTestModel,
     Join[Range[-170, -20, 10], Range[-12, 12, 1.5],
       Range[20, 170, 10]] Degree]};
  projectors = conditionedNullProjector /@ operators;
  dimensionpass = And @@ Table[
     MatrixQ[projectors[[k]]] &&
      Dimensions[projectors[[k]]] === {64, 64},
     {k, Length[projectors]}];
  If[dimensionpass,
   probe = Table[Exp[I p], {p, 64}];
   nserr = Max[Table[
      nullDelta = projectors[[k]].probe;
      Norm[operators[[k]].nullDelta]/
       Max[Norm[operators[[k]]] Norm[nullDelta], 10^-30],
      {k, Length[operators]}]];
   projectorError = Max[Table[
      Norm[projectors[[k]].projectors[[k]] - projectors[[k]]]/
       Max[Norm[projectors[[k]]], 10^-30],
      {k, Length[projectors]}]];
   projectorRanks = Round[Re[Tr[#]]] & /@ projectors;
   nullitypass = projectorRanks === {47, 15},
   nserr = Infinity; projectorError = Infinity;
   projectorRanks = {}; nullitypass = False];
  If[dimensionpass && nullitypass,
   bases = projectorNullBasis /@ projectors;
   basispass = ! MemberQ[bases, $Failed] &&
     And @@ (MatrixQ[#, NumericQ] & /@ bases);
   If[basispass,
    basisRanks = Dimensions[#][[2]] & /@ bases;
    basispass = basisRanks === {47, 15},
    basisRanks = {}],
   basispass = False; basisRanks = {}];
  If[basispass,
   inverseSmoke = inverseAblation[1, False];
   inversepass = inverseSmoke =!= $Failed && Length[inverseSmoke] >= 2;
   If[inversepass,
    inverseInitial = inverseSmoke[[1, 2]];
    inverseFinal = inverseSmoke[[-1, 2]];
    inverseDataChange = Max[inverseSmoke[[All, 3]]];
    inversepass = inverseFinal < .9 inverseInitial &&
      inverseDataChange < 10^-8,
    inverseInitial = Infinity; inverseFinal = Infinity;
    inverseDataChange = Infinity],
   inversepass = False; inverseInitial = Infinity;
   inverseFinal = Infinity; inverseDataChange = Infinity];
  pass = pecerr < .012 && smerr < 10^-10 && dimensionpass &&
    nullitypass && nserr < 10^-9 && projectorError < 10^-9 &&
    basispass && inversepass;
  Print["  PEC-series relative L2 error = ", ScientificForm[pecerr, 4]];
  Print["  Sherman-Morrison relative error = ", ScientificForm[smerr, 4]];
  Print["  Null-space invariance error = ", ScientificForm[nserr, 4]];
  Print["  Null-projector idempotence error = ",
    ScientificForm[projectorError, 4]];
  Print["  Null-projector dimensions = ", projectorRanks,
    "; dimension/rank checks = ", dimensionpass && nullitypass];
  Print["  Reduced null-basis dimensions = ", basisRanks,
    "; basis check = ", basispass];
  Print["  One-step inverse mismatch = ", inverseInitial, " -> ",
    inverseFinal, "; maximum relative data change = ",
    ScientificForm[inverseDataChange, 4]];
  Print[If[pass, "SelfTest PASSED.", "SelfTest FAILED."]];
  {"Passed" -> pass, "PECError" -> pecerr,
   "ShermanMorrisonError" -> smerr, "NullSpaceError" -> nserr,
   "ProjectorError" -> projectorError,
   "NullProjectorDimensionRankPassed" ->
    (dimensionpass && nullitypass), "NullBasisPassed" -> basispass,
   "InverseStepPassed" -> inversepass}
  ];

exportCSV[path_, header_List, rows_List] :=
  Export[path, Prepend[N[rows], header], "CSV"];

ensureDirectory[path_] := If[! DirectoryQ[path], CreateDirectory[path,
   CreateIntermediateDirectories -> True]];

(* The pre-audit RunStudy implementation was removed after verification. *)
(* ====================================================================== *)
(* PLOTTING AND REPRODUCIBILITY HELPERS -- MATHEMATICA 9 COMPATIBLE       *)
(* ====================================================================== *)

stepCoordinates[v_List] := Join[
  Flatten[Table[{{k, v[[k]]}, {k + 1, v[[k]]}},
    {k, 1, Length[v] - 1}], 1], {{Length[v], Last[v]}}];

normalizedSpectrum[c_] := Module[{s},
  s = SingularValueList[c]; s/First[s]];

(* ListPolarPlot in Mathematica 9 can pass an infinite radial range to its
   internal polar-axis routine.  These helpers draw the same polar data in
   Cartesian coordinates and avoid that version-specific failure. *)
polarCartesian[data_List] := Module[{xy},
  xy = ({#[[2]] Cos[#[[1]]], #[[2]] Sin[#[[1]]]} &) /@ data;
  Append[xy, First[xy]]
  ];

polarGridLines[maxRadius_?NumericQ, radialTicks_List] := Join[
  {Directive[GrayLevel[.82], Thin]},
  (Circle[{0, 0}, #[[1]]] &) /@
    Select[radialTicks, #[[1]] > 0 &],
  (Line[{{0, 0}, maxRadius {Cos[#], Sin[#]}}] &) /@
    (Range[0, 7] Pi/4)
  ];

polarGridLabels[maxRadius_?NumericQ, radialTicks_List] := Join[
  (Text[Style[#[[2]], 8, FontFamily -> "Times"],
       {#[[1]], 0}, {-1.15, -1.15}] &) /@ radialTicks,
  {Text[Style["0 deg", 9, FontFamily -> "Times"],
    {1.12 maxRadius, 0}],
   Text[Style["90 deg", 9, FontFamily -> "Times"],
    {0, 1.12 maxRadius}],
   Text[Style["180 deg", 9, FontFamily -> "Times"],
    {-1.12 maxRadius, 0}],
   Text[Style["-90 deg", 9, FontFamily -> "Times"],
    {0, -1.12 maxRadius}]}
  ];

codebookPlot[coords_List, q_List, colors_List, label_String,
  annotations_List] := Show[
  ListPlot[Table[Pick[coords, q, k], {k, Length[colors]}],
   PlotStyle -> Table[Directive[colors[[k]], PointSize[.018]],
     {k, Length[colors]}], AspectRatio -> 1, Frame -> True,
   PlotRange -> {{-2.35, 2.35}, {-2.35, 2.35}},
   FrameLabel -> {"x/lambda0", "y/lambda0"},
   PlotLabel -> Style[label, 12], GridLines -> Automatic,
   GridLinesStyle -> Directive[GrayLevel[.88]],
   BaseStyle -> {FontFamily -> "Times", 10}],
  Graphics[Join[{{GrayLevel[.7], Thin, Circle[{0, 0}, 2.]}}, annotations]]];

legacyDataFiles[] := {"frequency_sweep.csv", "incidence_sweep.csv",
  "inverse_convergence.csv", "mesh_convergence.csv", "monte_carlo.csv",
  "monte_carlo_percentiles.csv", "shared_states.csv",
  "telemetry_pattern.csv"};

removeLegacyData[dir_String] := Scan[
  If[FileExistsQ[#], DeleteFile[#]] &,
  (FileNameJoin[{dir, #}] &) /@ legacyDataFiles[]];

(* Park--Miller plus Box--Muller gives a documented random sequence that is
   identical across Mathematica versions and independent audit codes. *)
portableNormalSamples[count_Integer, state_Integer] := Module[
  {s, u, pairs, z},
  s = state;
  u = Table[s = Mod[16807 s, 2147483647]; N[s/2147483647],
    {2 Ceiling[count/2]}];
  pairs = Partition[u, 2];
  z = Flatten[Table[{Sqrt[-2 Log[pairs[[k, 1]]]]
       Cos[2 Pi pairs[[k, 2]]],
      Sqrt[-2 Log[pairs[[k, 1]]]] Sin[2 Pi pairs[[k, 2]]]},
     {k, Length[pairs]}]];
  {Take[z, count], s}
  ];

(* ====================================================================== *)
(* COMPLETE NUMERICAL STUDY AND ALL MANUSCRIPT OUTPUTS                    *)
(* ====================================================================== *)

RunStudy[OptionsPattern[]] := Module[
  {root, figdir, datadir, trials, states, qS, qT, design, valid,
   rDesign, rValid, tDesign, tValid, polished, meshRows, n, mult, m,
   rr, tt, freqRows, fs, invhist, p1, p2, p3, p4, fig1, fig2, fig3,
   fig4, fig5, figures, theta, thetaRad, alpha, b, f, f0, norm, feed,
   power, freePowerDB, mc, dr, dx, zs, zt, cases, deg, main, side,
   edgeCorrection, percentiles, percentileRows, metricNames, separateS,
   separateT, sepR, sepTM, summary, colors, coords, phi, geometryRows,
   inverseModel, specS, specT, specRows, maxSpec, plotHist, patternRows,
   controlRows, stateRows, scatterRows, solver, k, fieldPeak, table2Rows,
   row128, row256, row384, trialRows, randomState, normals},
  root = OptionValue[StudyDirectory];
  If[root === Automatic, root = $PackageDirectory];
  figdir = FileNameJoin[{root, "figures"}];
  datadir = FileNameJoin[{root, "data"}];
  ensureDirectory[figdir]; ensureDirectory[datadir]; removeLegacyData[datadir];
  trials = OptionValue[MonteCarloTrials];
  states = publishedStates[]; qS = publishedQS[]; qT = publishedQT[];

  Print["DualModeRIS revision: accepted-step loop repair (2026-10-24e)."];
  Print["RunStudy: building 128-unknown synthesis model..."];
  design = buildModel[128, 2., 1., 360];
  If[TrueQ[OptionValue[RefineCodebooks]],
   Print["RunStudy: coordinate polishing enabled."];
   polished = coordinatePolish[design, states, qS, 2, "S", 3];
   qS = First[polished];
   polished = coordinatePolish[design, states, qT, 2, "T", 3];
   qT = First[polished]];
  rDesign = protectedRatios[design, states, qS, 2];
  tDesign = telemetryMetrics[design, states, qT, 2];
  Print["RunStudy: building independent 256-unknown validation model..."];
  valid = buildModel[256, 2., 1., 720];
  rValid = protectedRatios[valid, states, qS, 4];
  tValid = telemetryMetrics[valid, states, qT, 4];

  (* ==================================================================== *)
  (* FIGURE 1 AND CSVs: GEOMETRY, SHARED STATES, AND TWO CONTROL WORDS     *)
  (* ==================================================================== *)
  colors = Table[ColorData["DarkRainbow"][N[(k - 1)/7]], {k, 8}];
  phi = Range[0, 63] 2 Pi/64; coords = Transpose[{2 Cos[phi], 2 Sin[phi]}];
  stateRows = Table[{k, Re[states[[k]]], Im[states[[k]]]},
    {k, Length[states]}];
  geometryRows = Table[{k, phi[[k]]/Degree, coords[[k, 1]], coords[[k, 2]],
     qS[[k]], Re[states[[qS[[k]]]]], Im[states[[qS[[k]]]]],
     qT[[k]], Re[states[[qT[[k]]]]], Im[states[[qT[[k]]]]]}, {k, 64}];
  exportCSV[FileNameJoin[{datadir, "fig1_shared_impedance_states.csv"}],
    {"state_index", "normalized_resistance_re_z_over_eta",
     "normalized_reactance_im_z_over_eta"}, stateRows];
  exportCSV[FileNameJoin[{datadir, "fig1_configuration_codebooks.csv"}],
    {"physical_cell_index", "cell_center_angle_deg", "x_over_lambda0",
     "y_over_lambda0", "low_scattering_state_index",
     "low_scattering_re_z_over_eta", "low_scattering_im_z_over_eta",
     "telemetry_state_index", "telemetry_re_z_over_eta",
     "telemetry_im_z_over_eta"}, geometryRows];
  p1 = codebookPlot[coords, qS, colors, "(a) Low-scattering codebook",
    {{Black, Thick, Arrow[{{-1.78, 1.72}, {-1.92, 1.08}}]},
      Text[Style["plane wave", 10], {-1.55, 1.92}],
      Arrow[{{1.25, -.95}, {2.02, 0.}}],
      Text[Style["protected monostatic sector", 9], {.55, -1.35}]}];
  p2 = codebookPlot[coords, qT, colors, "(b) Telemetry codebook",
    {{Red, PointSize[.025], Point[{1.45, 0.}]},
      {Black, Thick, Arrow[{{1.45, 0.}, {2.05, 0.}}]},
      Text[Style["theta0 = 0 deg", 10], {1.45, 1.45}],
      Text[Style["line feed, xf = 1.45 lambda0", 9], {.15, -1.55}]}];
  fig1 = GraphicsGrid[{{p1, p2}}, ImageSize -> 900, Spacings -> .35];
  Export[FileNameJoin[{figdir, "fig1_configuration.pdf"}], fig1];

  (* ==================================================================== *)
  (* FIGURE 2 AND CSVs: INVERSE SPECTRA AND CORRECTED ITERATION HISTORY    *)
  (* ==================================================================== *)
  Print["RunStudy: corrected inverse-source ablation..."];
  inverseModel = buildModel[64, 2., 1., 360];
  specS = normalizedSpectrum[farMatrix[inverseModel,
      Range[-12, 12, 1.5] Degree]];
  specT = normalizedSpectrum[farMatrix[inverseModel,
      Join[Range[-170, -20, 10], Range[-12, 12, 1.5],
       Range[20, 170, 10]] Degree]];
  maxSpec = Max[Length[specS], Length[specT]];
  specRows = Table[{k - 1,
     If[k <= Length[specS], specS[[k]], ""],
     If[k <= Length[specT], specT[[k]], ""]}, {k, maxSpec}];
  exportCSV[FileNameJoin[{datadir,
     "fig2a_inverse_normalized_singular_spectra.csv"}],
    {"singular_value_index_zero_based",
     "low_scattering_operator_sigma_i_over_sigma_1",
     "telemetry_operator_sigma_i_over_sigma_1"}, specRows];
  invhist = inverseAblation[12];
  If[invhist === $Failed, Return[$Failed]];
  If[Last[invhist][[2]] >= .05 || Max[invhist[[All, 3]]] >= 10^-8,
   Print["ERROR: inverse-source null-space diagnostic failed. ",
    "Initial mismatch = ", invhist[[1, 2]],
    "; final mismatch = ", invhist[[-1, 2]],
    "; maximum relative data change = ", Max[invhist[[All, 3]]], "."];
   Return[$Failed]];
  exportCSV[FileNameJoin[{datadir,
     "fig2b_inverse_iteration_history.csv"}],
    {"alternating_iteration", "normalized_shared_library_mismatch",
     "maximum_relative_change_in_restricted_data"}, invhist];
  p1 = ListLogPlot[{Transpose[{Range[0, Length[specS] - 1], specS}],
      Transpose[{Range[0, Length[specT] - 1], specT}]}, Joined -> True,
    PlotMarkers -> Automatic, Frame -> True,
    FrameLabel -> {"singular-value index", "sigma_i / sigma_1"},
    PlotLabel -> Style["(a) Restricted radiation operators", 12],
    PlotStyle -> {{Blue, Thick}, {Orange, Thick}},
    PlotLegends -> Placed[{"scattering operator", "telemetry operator"},
      {0.34, .18}], GridLines -> Automatic,
    BaseStyle -> {FontFamily -> "Times", 10}];
  plotHist = invhist;
  plotHist[[All, 3]] = Max[#, 10^-16] & /@ plotHist[[All, 3]];
  p2 = ListLogPlot[{plotHist[[All, {1, 2}]], plotHist[[All, {1, 3}]]},
    Joined -> True, PlotMarkers -> Automatic, Frame -> True,
    FrameLabel -> {"alternating iteration", "normalized error"},
    PlotLabel -> Style["(b) Null-space iterations", 12],
    PlotStyle -> {{Blue, Thick}, {Orange, Dashed, Thick}},
    PlotLegends -> Placed[{"shared-library mismatch",
       "relative data change"}, {0.67, .52}], GridLines -> Automatic,
    BaseStyle -> {FontFamily -> "Times", 10}];
  fig2 = GraphicsGrid[{{p1, p2}}, ImageSize -> 900, Spacings -> .35];
  Export[FileNameJoin[{figdir, "fig2_inverse.pdf"}], fig2];

  (* ==================================================================== *)
  (* FIGURE 3 AND CSVs: INCIDENCE SWEEP AND NORMAL-INCIDENCE PATTERN       *)
  (* ==================================================================== *)
  Print["RunStudy: low-scattering validation patterns..."];
  scatterRows = Table[{Range[-30, 30, 5][[k]], rValid[[k]],
      db10[rValid[[k]]]}, {k, Length[rValid]}];
  exportCSV[FileNameJoin[{datadir, "fig3a_incidence_sweep.csv"}],
    {"incidence_angle_deg", "worst_protected_sector_power_ratio_linear",
     "worst_protected_sector_power_ratio_db_relative_to_pec"}, scatterRows];
  p1 = ListLinePlot[{scatterRows[[All, {1, 3}]],
      Transpose[{Range[-30, 30, 5], ConstantArray[0., Length[rValid]]}]},
    Joined -> True, PlotMarkers -> {Automatic, None}, Frame -> True,
    PlotRange -> {{-32, 32}, {-17, 1}},
    FrameLabel -> {"incidence angle alpha (deg)",
      "worst protected-sector ratio (dB)"},
    PlotLabel -> Style["(a) +/-10 deg monostatic sectors", 12],
    PlotStyle -> {{Blue, Thick}, {GrayLevel[.3], Dashed, Thick}},
    PlotLegends -> Placed[{"shared eight-state design", "PEC reference"},
      {0.68, .55}], GridLines -> Automatic,
    BaseStyle -> {FontFamily -> "Times", 10}];
  alpha = 0.; b = planeWave[valid, alpha];
  f = farField[valid, solveCurrent[valid, stateProfile[states, qS, 4], b]];
  f0 = farField[valid, solveCurrent[valid,
      ConstantArray[0. + 0. I, 256], b]];
  norm = Max[Abs[f0]^2]; theta = value[valid, "theta"]/Degree;
  patternRows = Table[{theta[[k]], db10[Abs[f0[[k]]]^2/norm],
      db10[Abs[f[[k]]]^2/norm]}, {k, Length[theta]}];
  exportCSV[FileNameJoin[{datadir,
     "fig3b_normal_incidence_scattering_pattern.csv"}],
    {"observation_angle_deg", "pec_scattered_power_db_relative_to_pec_peak",
     "shared_design_scattered_power_db_relative_to_pec_peak"}, patternRows];
  p2 = ListLinePlot[{patternRows[[All, {1, 2}]],
      patternRows[[All, {1, 3}]]}, Joined -> True,
    PlotRange -> {{-180, 180}, {-50, 3}}, Frame -> True,
    FrameLabel -> {"observation angle theta (deg)",
      "scattered power / PEC peak (dB)"},
    PlotLabel -> Style["(b) Pattern at alpha = 0 deg", 12],
    PlotStyle -> {{GrayLevel[.35], Dashed, Thick}, {Blue, Thick}},
    PlotLegends -> Placed[{"PEC", "shared eight-state design"},
      {0.55, .16}], GridLines -> Automatic,
    Prolog -> {Directive[Red, Opacity[.08]],
      Rectangle[{-180, -50}, {-170, 3}], Rectangle[{170, -50}, {180, 3}]},
    BaseStyle -> {FontFamily -> "Times", 10}];
  fig3 = GraphicsGrid[{{p1, p2}}, ImageSize -> 900, Spacings -> .35];
  Export[FileNameJoin[{figdir, "fig3_scattering.pdf"}], fig3];

  (* ==================================================================== *)
  (* FIGURE 4 AND CSVs: POLAR TELEMETRY PATTERN AND 64-CELL CONTROL WORD  *)
  (* ==================================================================== *)
  Print["RunStudy: telemetry polar pattern and control word..."];
  power = Abs["Field" /. tValid]^2; fieldPeak = Max[power];
  power = power/fieldPeak; thetaRad = value[valid, "theta"];
  freePowerDB = db10[(value[valid, "k"] value[valid, "eta"]/4)^2/
     fieldPeak];
  patternRows = Table[{theta[[k]], power[[k]], db10[power[[k]]],
      freePowerDB}, {k, Length[theta]}];
  exportCSV[FileNameJoin[{datadir, "fig4a_telemetry_polar_pattern.csv"}],
    {"observation_angle_deg", "shared_design_normalized_power_linear",
     "shared_design_normalized_power_db",
     "free_space_line_source_power_db_relative_to_design_peak"}, patternRows];
  controlRows = Table[{k, phi[[k]]/Degree, qT[[k]],
      Re[states[[qT[[k]]]]], Im[states[[qT[[k]]]]]}, {k, 64}];
  exportCSV[FileNameJoin[{datadir, "fig4b_telemetry_control_word.csv"}],
    {"physical_cell_index", "cell_center_angle_deg", "selected_state_index",
     "selected_normalized_resistance_re_z_over_eta",
     "selected_normalized_reactance_im_z_over_eta"}, controlRows];
  p3 = ListLinePlot[{
      polarCartesian[Transpose[{thetaRad,
         35 + (Max[-35., #] & /@ (db10 /@ power))}]],
      polarCartesian[Transpose[{thetaRad, ConstantArray[
         35 + freePowerDB, Length[thetaRad]]}]]},
    Joined -> True, PlotRange -> {{-41, 41}, {-41, 41}},
    AspectRatio -> 1, Axes -> False, Frame -> False,
    Prolog -> polarGridLines[35,
      {{5, "-30"}, {15, "-20"}, {25, "-10"}, {35, "0"}}],
    Epilog -> polarGridLabels[35,
      {{5, "-30"}, {15, "-20"}, {25, "-10"}, {35, "0"}}],
    PlotStyle -> {{Blue, Thick}, {GrayLevel[.35], Dashed, Thick}},
    PlotLabel -> Style["(a) Normalized telemetry pattern", 12],
    PlotLegends -> Placed[{"shared eight-state design",
       "free-space line source"}, Below],
    BaseStyle -> {FontFamily -> "Times", 10}];
  p4 = ListLinePlot[stepCoordinates[qT], Frame -> True,
    PlotRange -> {{.5, 64.5}, {.7, 8.3}},
    FrameLabel -> {"physical cell index", "selected state"},
    PlotLabel -> Style["(b) Telemetry control word", 12],
    PlotStyle -> {Blue, Thick}, GridLines -> Automatic,
    BaseStyle -> {FontFamily -> "Times", 10}];
  fig4 = GraphicsGrid[{{p3, p4}}, ImageSize -> 900, Spacings -> .35];
  Export[FileNameJoin[{figdir, "fig4_telemetry.pdf"}], fig4];

  (* ==================================================================== *)
  (* FIGURE 5 AND CSVs: MESH CONVERGENCE AND FREQUENCY SENSITIVITY        *)
  (* ==================================================================== *)
  Print["RunStudy: mesh and frozen-impedance frequency checks..."];
  meshRows = Table[
    mult = n/64; m = buildModel[n, 2., 1., 720];
    rr = protectedRatios[m, states, qS, mult];
    tt = telemetryMetrics[m, states, qT, mult];
    {n, db10[Max[rr]], db10[Mean[rr]], "Efficiency" /. tt,
     db10["Directivity" /. tt], db10["SLL" /. tt], "BeamAngle" /. tt},
    {n, {128, 192, 256, 384}}];
  exportCSV[FileNameJoin[{datadir, "fig5a_mesh_convergence.csv"}],
    {"mom_unknowns", "worst_protected_sector_power_db_relative_to_pec",
     "mean_protected_sector_power_db_relative_to_pec",
     "telemetry_main_sector_efficiency", "telemetry_directivity_db",
     "telemetry_peak_sidelobe_level_db", "telemetry_beam_angle_deg"},
    meshRows];
  freqRows = Table[
    m = buildModel[256, 2., fs, 720];
    rr = protectedRatios[m, states, qS, 4];
    tt = telemetryMetrics[m, states, qT, 4];
    {fs, db10[Max[rr]], db10[Mean[rr]],
     db10["Directivity" /. tt], db10["SLL" /. tt]},
    {fs, N[Range[.94, 1.06, .01]]}];
  exportCSV[FileNameJoin[{datadir,
     "fig5b_frozen_impedance_frequency_sweep.csv"}],
    {"frequency_ratio_f_over_f0",
     "worst_protected_sector_power_db_relative_to_pec",
     "mean_protected_sector_power_db_relative_to_pec",
     "telemetry_directivity_db", "telemetry_peak_sidelobe_level_db"},
    freqRows];
  p1 = ListLinePlot[{meshRows[[All, {1, 2}]], meshRows[[All, {1, 3}]],
      meshRows[[All, {1, 5}]]}, Joined -> True, PlotMarkers -> Automatic,
    Frame -> True, FrameLabel -> {"MoM unknowns", "metric (dB)"},
    PlotLabel -> Style["(a) Mesh refinement", 12],
    PlotStyle -> {{Blue, Thick}, {Orange, Thick}, {Darker[Green], Thick}},
    PlotLegends -> Placed[{"worst sector ratio", "mean sector ratio",
       "telemetry directivity"}, {0.69, .47}], GridLines -> Automatic,
    BaseStyle -> {FontFamily -> "Times", 10}];
  p2 = ListLinePlot[{freqRows[[All, {1, 2}]], freqRows[[All, {1, 4}]]},
    Joined -> True, PlotMarkers -> Automatic, Frame -> True,
    FrameLabel -> {"frequency f/f0", "metric (dB)"},
    PlotLabel -> Style["(b) Frozen-impedance sensitivity", 12],
    PlotStyle -> {{Blue, Thick}, {Orange, Thick}},
    PlotLegends -> Placed[{"worst sector ratio", "telemetry directivity"},
      {0.32, .2}], GridLines -> Automatic,
    Epilog -> {GrayLevel[.35], Dashed, Thick, Line[{{1., -12}, {1., 14}}]},
    BaseStyle -> {FontFamily -> "Times", 10}];
  fig5 = GraphicsGrid[{{p1, p2}}, ImageSize -> 900, Spacings -> .35];
  Export[FileNameJoin[{figdir, "fig5_validation.pdf"}], fig5];

  (* ==================================================================== *)
  (* TABLE 3 CSVs: FIVE-PERCENT CELL-IMPEDANCE MONTE CARLO STUDY           *)
  (* ==================================================================== *)
  Print["RunStudy: ", trials, " five-percent impedance-error trials..."];
  m = buildModel[256, 2., 1., 360];
  cases = Table[
    alpha = deg Degree; b = planeWave[m, alpha];
    f0 = farField[m, solveCurrent[m, ConstantArray[0. + 0. I, 256], b]];
    {b, protectedIndices[m, alpha, 10 Degree],
     Max[Abs[f0[[protectedIndices[m, alpha, 10 Degree]]]]^2]},
    {deg, Range[-30, 30, 5]}];
  feed = {value[m, "a"] - .55, 0.}; thetaRad = value[m, "theta"];
  main = Select[Range[Length[thetaRad]],
    Abs[wrapAngle[thetaRad[[#]]]] <= 8 Degree + 100 $MachineEpsilon &];
  side = Select[Range[Length[thetaRad]],
    Abs[wrapAngle[thetaRad[[#]]]] >= 15 Degree - 100 $MachineEpsilon &];
  randomState = 271828; mc = Table[
    {normals, randomState} = portableNormalSamples[64, randomState];
    dr = .05 normals;
    {normals, randomState} = portableNormalSamples[64, randomState];
    dx = .05 normals;
    zs = expandCode[.03 (1 + dr) + I Im[states[[qS]]] (1 + dx), 4];
    solver = LinearSolve[DiagonalMatrix[zs] + value[m, "S"]];
    rr = Table[
      f = farField[m, solver[cases[[k, 1]]]];
      Max[Abs[f[[cases[[k, 2]]]]]^2]/cases[[k, 3]], {k, Length[cases]}];
    {normals, randomState} = portableNormalSamples[64, randomState];
    dr = .05 normals;
    {normals, randomState} = portableNormalSamples[64, randomState];
    dx = .05 normals;
    zt = expandCode[.03 (1 + dr) + I Im[states[[qT]]] (1 + dx), 4];
    f = farField[m, solveCurrent[m, zt, lineFeed[m, feed]], feed];
    power = Abs[f]^2;
    edgeCorrection = (power[[First[main]]] + power[[Last[main]]])/2.;
    {db10[Max[rr]], db10[Mean[rr]],
     (Total[power[[main]]] - edgeCorrection)/Total[power],
     db10[Max[power[[side]]]/Max[power[[main]]]],
     db10[Max[power]/Mean[power]],
     thetaRad[[First[First[Position[power, Max[power]]]]]]/Degree},
    {trials}];
  trialRows = MapIndexed[Prepend[#1, First[#2]] &, mc];
  exportCSV[FileNameJoin[{datadir, "table3_monte_carlo_trials.csv"}],
    {"trial_index", "worst_protected_sector_power_db_relative_to_pec",
     "mean_protected_sector_power_db_relative_to_pec",
     "telemetry_main_sector_efficiency",
     "telemetry_peak_sidelobe_level_db", "telemetry_directivity_db",
     "telemetry_beam_angle_deg"}, trialRows];
  metricNames = {"worst_protected_sector_power_db_relative_to_pec",
    "mean_protected_sector_power_db_relative_to_pec",
    "telemetry_main_sector_efficiency", "telemetry_peak_sidelobe_level_db",
    "telemetry_directivity_db", "telemetry_beam_angle_deg"};
  percentiles = Table[Quantile[mc[[All, k]], {.05, .5, .95}],
    {k, Dimensions[mc][[2]]}];
  percentileRows = Table[Join[{metricNames[[k]]}, percentiles[[k]]],
    {k, Length[metricNames]}];
  exportCSV[FileNameJoin[{datadir,
     "table3_monte_carlo_percentiles.csv"}],
    {"metric", "fifth_percentile", "median_fiftieth_percentile",
     "ninety_fifth_percentile"}, percentileRows];

  (* ==================================================================== *)
  (* TABLE 2 CSV: NOMINAL FORWARD-VERIFIED PERFORMANCE                    *)
  (* ==================================================================== *)
  separateS = .03 + I {-4.79869572, -1.78377682, -.49220321, -.25632662,
      .18545799, .58122016, 1.28946036, 4.19551841};
  separateT = .03 + I {-7.99999977, -1.01923389, -.32473110, -.11836519,
      .16402620, .52228471, 1.08825453, 1.45341219};
  sepR = protectedRatios[valid, separateS, qS, 4];
  sepTM = telemetryMetrics[valid, separateT, qT, 4];
  row128 = First[Select[meshRows, First[#] == 128 &]];
  row256 = First[Select[meshRows, First[#] == 256 &]];
  row384 = First[Select[meshRows, First[#] == 384 &]];
  table2Rows = {
    {"shared", "synthesis", 128, db10[Max[rDesign]],
     db10[Mean[rDesign]], "Efficiency" /. tDesign,
     db10["Directivity" /. tDesign], db10["SLL" /. tDesign],
     "BeamAngle" /. tDesign},
    Join[{"shared", "validation"}, row256],
    Join[{"shared", "validation"}, row384],
    {"separate_8_plus_8", "validation", 256, db10[Max[sepR]],
     db10[Mean[sepR]], "Efficiency" /. sepTM,
     db10["Directivity" /. sepTM], db10["SLL" /. sepTM],
     "BeamAngle" /. sepTM}};
  exportCSV[FileNameJoin[{datadir, "table2_forward_verified_performance.csv"}],
    {"impedance_library", "model_role", "mom_unknowns",
     "worst_protected_sector_power_db_relative_to_pec",
     "mean_protected_sector_power_db_relative_to_pec",
     "telemetry_main_sector_efficiency", "telemetry_directivity_db",
     "telemetry_peak_sidelobe_level_db", "telemetry_beam_angle_deg"},
    table2Rows];

  summary = {"DesignWorstSectorDB" -> db10[Max[rDesign]],
    "ValidationWorstSectorDB" -> db10[Max[rValid]],
    "ValidationMeanSectorDB" -> db10[Mean[rValid]],
    "ValidationDirectivityDB" -> db10["Directivity" /. tValid],
    "ValidationSLLDB" -> db10["SLL" /. tValid],
    "ValidationEfficiency" -> ("Efficiency" /. tValid),
    "InverseInitialMismatch" -> invhist[[1, 2]],
    "InverseFinalMismatch" -> invhist[[-1, 2]],
    "InverseMaximumRelativeDataChange" -> Max[invhist[[All, 3]]],
    "SeparateLibraryWorstSectorDB" -> db10[Max[sepR]],
    "SeparateLibraryDirectivityDB" -> db10["Directivity" /. sepTM]};
  Print["RunStudy complete. Summary: ", summary];
  figures = {fig1, fig2, fig3, fig4, fig5};
  If[TrueQ[OptionValue[DisplayFigures]],
   Print[Column[figures, Spacings -> 2]]];
  summary
  ];

End[];
EndPackage[];
