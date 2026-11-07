(* ::Package:: *)

(* Mathematica 9 compatible independent Fourier--Galerkin validation for
   the manuscript "Null-Space-Assisted Inverse Synthesis of Dual-Mode
   Conformal Metasurfaces With Shared Impedance States".

   This package is intentionally self-contained.  It does not load, call,
   or reuse DualModeRIS.m.  The frozen physical parameters, impedance
   states, and codebooks are repeated below so that this calculation can
   serve as an independently implemented cross-formulation check. *)

BeginPackage["FourierGalerkinValidation`"];

FGSelfTest::usage =
  "FGSelfTest[] runs normalization, Fourier-coefficient, multiplication-matrix, closed-form constant-impedance, PEC quadrature, and projected-residual checks.";

FGRunValidation::usage =
  "FGRunValidation[] independently evaluates the frozen low-scattering and telemetry configurations with a Fourier--Galerkin formulation, exports traceable CSV and PDF files, and displays the validation figures.";

FGValidationDirectory::usage =
  "FGValidationDirectory is an option for FGRunValidation specifying the result directory. Automatic creates fourier_galerkin_results beside this package.";

FGFourierOrders::usage =
  "FGFourierOrders is an option for FGRunValidation specifying the retained half-orders M. The number of complex Fourier unknowns is 2 M + 1.";

FGAngularSamples::usage =
  "FGAngularSamples is an option for FGRunValidation specifying the number of uniformly spaced far-field observation angles.";

FGWorkingPrecision::usage =
  "FGWorkingPrecision specifies the decimal precision used to form stable high-order Bessel/Hankel products before the dense systems are converted to machine precision.";

FGDisplayFigures::usage =
  "FGDisplayFigures is an option for FGRunValidation. True prints the exported validation figures in the launcher notebook.";

FGRunSelfTest::usage =
  "FGRunSelfTest is an option for FGRunValidation. True runs FGSelfTest before the production calculation and stops if a check fails.";

Options[FGSelfTest] = {FGWorkingPrecision -> 50};

Options[FGRunValidation] = {
  FGValidationDirectory -> Automatic,
  FGFourierOrders -> {64, 128, 192, 256},
  FGAngularSamples -> 720,
  FGWorkingPrecision -> 50,
  FGDisplayFigures -> True,
  FGRunSelfTest -> True
  };

Begin["`Private`"];

$FGPackageDirectory = DirectoryName[$InputFileName];

(* ====================================================================== *)
(* FROZEN INPUT DATA                                                      *)
(* ====================================================================== *)

fgStates[] := {
  0.03 - 6.75602121 I,
  0.03 - 1.61813346 I,
  0.03 - 0.56199291 I,
  0.03 - 0.22403657 I,
  0.03 + 0.18768728 I,
  0.03 + 0.60721582 I,
  0.03 + 1.34454867 I,
  0.03 + 4.93228419 I
  };

fgLowCodebook[] := {
  7, 6, 7, 7, 8, 5, 4, 1, 8, 4, 7, 4, 6, 5, 6, 7,
  1, 5, 5, 5, 1, 4, 1, 4, 8, 3, 5, 6, 4, 1, 6, 5,
  6, 3, 7, 4, 1, 4, 1, 2, 1, 2, 1, 3, 1, 5, 6, 5,
  5, 4, 7, 4, 5, 4, 5, 3, 2, 5, 6, 6, 5, 8, 4, 7
  };

fgTelemetryCodebook[] := {
  4, 4, 1, 4, 5, 4, 6, 4, 5, 7, 3, 5, 4, 4, 4, 4,
  3, 5, 5, 4, 4, 5, 7, 5, 4, 4, 3, 3, 5, 7, 1, 1,
  3, 1, 6, 5, 3, 4, 4, 2, 5, 3, 6, 5, 4, 4, 5, 4,
  4, 5, 4, 4, 6, 6, 3, 6, 6, 2, 7, 6, 5, 4, 8, 4
  };

fgMoM256[] := {
  "WorstSectorDB" -> -11.575114266751022,
  "MeanSectorDB" -> -13.08853693643261,
  "Efficiency" -> 0.6490546670574209,
  "DirectivityDB" -> 13.098468196464502,
  "SLLDB" -> -13.60126936499137,
  "BeamAngleDeg" -> 0.9999999999999887
  };

fgMoM384[] := {
  "WorstSectorDB" -> -11.384771355642515,
  "MeanSectorDB" -> -12.878771183356072,
  "Efficiency" -> 0.643548549120784,
  "DirectivityDB" -> 13.0818961413309,
  "SLLDB" -> -13.472225800508063,
  "BeamAngleDeg" -> 1.5000000000000084
  };

(* ====================================================================== *)
(* GENERAL UTILITIES                                                      *)
(* ====================================================================== *)

fgValue[rules_List, key_String] := key /. rules;

fgWrapAngle[x_] := Mod[x + Pi, 2 Pi] - Pi;

fgDB10[x_?NumericQ] := 10 Log[10, Max[N[x], 10.^-300]];

fgEnsureDirectory[dir_String] :=
  If[! DirectoryQ[dir],
   CreateDirectory[dir, CreateIntermediateDirectories -> True]];

fgExportCSV[path_String, header_List, rows_List] :=
  Export[path, Prepend[rows, header], "CSV"];

fgStatus[condition_] := If[TrueQ[condition], "PASS", "FAIL"];

fgFiniteNumberQ[x_] := NumberQ[N[x]] &&
  FreeQ[N[x], Indeterminate | ComplexInfinity | DirectedInfinity[_]];

fgRelativeError[x_, y_] :=
  Norm[N[x - y]]/Max[Norm[N[y]], 10.^-300];

fgPolarCartesian[data_List] := Module[{xy},
  xy = ({#[[2]] Cos[#[[1]]], #[[2]] Sin[#[[1]]]} &) /@ data;
  Append[xy, First[xy]]
  ];

fgPolarGridLines[maxRadius_?NumericQ, radialTicks_List] := Join[
  {Directive[GrayLevel[.82], Thin]},
  (Circle[{0, 0}, #[[1]]] &) /@
    Select[radialTicks, #[[1]] > 0 &],
  (Line[{{0, 0}, maxRadius {Cos[#], Sin[#]}}] &) /@
    (Range[0, 7] Pi/4)
  ];

fgPolarGridLabels[maxRadius_?NumericQ, radialTicks_List] := Join[
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

(* ====================================================================== *)
(* FOURIER COEFFICIENTS OF THE PIECEWISE-CONSTANT PROFILE                 *)
(* ====================================================================== *)

(* Cell l occupies [l Delta,(l+1) Delta), where Delta=2 Pi/L.  This is
   the continuum interpretation of the repeated-cell midpoint MoM data.
   The exact Fourier coefficient is

     v_p = sinc(p Pi/L)/L Sum_l v_l exp[-I p (l+1/2) Delta].

   Nonzero integer multiples of L vanish exactly. *)

fgPiecewiseCoefficient[values_List, p_Integer] := Module[
  {l, delta, count, factor},
  count = Length[values];
  If[p == 0, Return[N[Mean[values]]]];
  If[Mod[Abs[p], count] == 0, Return[0. + 0. I]];
  delta = 2 Pi/count;
  factor = Sin[p delta/2]/(p delta/2);
  N[factor/count Sum[
      values[[l + 1]] Exp[-I p (l + 1/2) delta],
      {l, 0, count - 1}]]
  ];

fgMultiplicationMatrix[values_List, m_Integer] := Module[
  {coefficients},
  coefficients = fgPiecewiseCoefficient[values, #] & /@
    Range[-2 m, 2 m];
  Table[coefficients[[n - q + 2 m + 1]],
    {n, -m, m}, {q, -m, m}]
  ];

fgPiecewiseSamples[values_List, angles_List] := Module[
  {count},
  count = Length[values];
  Table[values[[1 + Mod[Floor[count Mod[angles[[q]], 2 Pi]/(2 Pi)],
        count]]], {q, Length[angles]}]
  ];

(* ====================================================================== *)
(* CIRCULAR FOURIER--GALERKIN OPERATORS                                   *)
(* ====================================================================== *)

(* The time convention is exp[-I omega t].  For radius a and k=2 Pi,

   S_n = (Pi k eta a/2) J_n(k a) H_n^(1)(k a).

   Products are formed before conversion to machine precision.  This
   avoids spurious overflow/underflow at the larger retained orders. *)

fgModalData[m_Integer, wp_Integer] := Module[
  {indices, kHP, etaHP, aHP, rfHP, xHP, sPositive, linePositive,
   jPositive, invHPositive, parityValue, s, line, jn, invH,
   farWeights, k, eta, a, rf},
  indices = Range[-m, m];
  kHP = N[2 Pi, wp]; etaHP = N[1, wp];
  aHP = N[2, wp]; rfHP = N[1.45, wp]; xHP = kHP aHP;
  sPositive = Table[
    N[(Pi kHP etaHP aHP/2) BesselJ[q, xHP] HankelH1[q, xHP],
      MachinePrecision], {q, 0, m}];
  linePositive = Table[
    N[-kHP etaHP/4 HankelH1[q, xHP] BesselJ[q, kHP rfHP],
      MachinePrecision], {q, 0, m}];
  jPositive = Table[
    N[BesselJ[q, xHP], MachinePrecision], {q, 0, m}];
  invHPositive = Table[
    N[1/HankelH1[q, xHP], MachinePrecision], {q, 0, m}];
  parityValue[list_List, q_Integer] :=
    If[q >= 0, list[[q + 1]], (-1)^Abs[q] list[[Abs[q] + 1]]];
  s = Table[sPositive[[Abs[q] + 1]], {q, indices}];
  line = Table[linePositive[[Abs[q] + 1]], {q, indices}];
  jn = Table[parityValue[jPositive, q], {q, indices}];
  invH = Table[parityValue[invHPositive, q], {q, indices}];
  k = N[kHP, MachinePrecision]; eta = 1.;
  a = 2.; rf = 1.45;
  farWeights = N[-k eta/4 2 Pi a (-I)^indices jn];
  {
   "M" -> m, "Indices" -> indices, "k" -> k, "eta" -> eta,
   "a" -> a, "FeedRadius" -> rf, "SingleLayer" -> s,
   "LineExcitation" -> line, "BesselJka" -> jn,
   "InverseHankel" -> invH, "FarWeights" -> farWeights
   }
  ];

fgPlaneExcitation[modal_List, alpha_?NumericQ] := Module[{n},
  n = fgValue[modal, "Indices"];
  N[I^n fgValue[modal, "BesselJka"] Exp[-I n alpha]]
  ];

fgPECcurrent[modal_List, alpha_?NumericQ] := Module[
  {n, k, eta, a},
  n = fgValue[modal, "Indices"];
  k = fgValue[modal, "k"]; eta = fgValue[modal, "eta"];
  a = fgValue[modal, "a"];
  N[2 I^n Exp[-I n alpha] fgValue[modal, "InverseHankel"]/
    (Pi k eta a)]
  ];

(* The primary unknown is u=Z J, the total surface field.  With Y=1/Z,
   the continuous equation becomes u+S(Y u)=b.  Galerkin projection gives

      [delta_nm + S_n y_(n-m)] u_m = b_n.

   Expanding the smoother u and applying the inverse admittance product
   rule is important for the discontinuous 64-cell profile. *)

fgBuildSystem[states_List, codebook_List, modal_List] := Module[
  {m, impedance, admittance, yMatrix, matrix},
  m = fgValue[modal, "M"];
  impedance = states[[codebook]];
  admittance = 1/impedance;
  yMatrix = fgMultiplicationMatrix[admittance, m];
  matrix = IdentityMatrix[2 m + 1] +
    DiagonalMatrix[fgValue[modal, "SingleLayer"]].yMatrix;
  {"Impedance" -> impedance, "Admittance" -> admittance,
   "YMatrix" -> yMatrix, "Matrix" -> matrix}
  ];

fgFarBasis[modal_List, angles_List] :=
  Exp[I Outer[Times, angles, fgValue[modal, "Indices"]]];

fgFarFieldWithAngles[modal_List, basis_, angles_List, current_List,
   includeFeed_] := Module[{field, k, eta, rf},
  field = basis.(fgValue[modal, "FarWeights"] current);
  If[TrueQ[includeFeed],
   k = fgValue[modal, "k"]; eta = fgValue[modal, "eta"];
   rf = fgValue[modal, "FeedRadius"];
   field = field - k eta/4 Exp[-I k rf Cos[angles]]];
  field
  ];

fgProtectedIndices[angles_List, alpha_?NumericQ,
   halfWidth_?NumericQ] :=
  Select[Range[Length[angles]],
   Abs[fgWrapAngle[angles[[#]] - (alpha + Pi)]] <=
      halfWidth + 100 $MachineEpsilon &];

(* ====================================================================== *)
(* ONE FOURIER ORDER: BOTH FROZEN CONFIGURATIONS                          *)
(* ====================================================================== *)

fgOneOrderStudy[m_Integer, angularSamples_Integer, wp_Integer] := Module[
  {timing, result, states, qLow, qTelemetry, modal, angles, basis,
   lowSystem, lowMatrix, lowY, lowSolver, incidenceDegrees,
   lowRatios, lowResiduals, lowFields, pecFields, alpha, b, u, j,
   field, jPEC, fieldPEC, protected, denominator, telemetrySystem,
   telemetryMatrix, telemetryY, telemetrySolver, bTelemetry,
   uTelemetry, jTelemetry, telemetryField, telemetryResidual,
   power, main, side, mainIndices, sideIndices, edgeCorrection,
   efficiency, directivity, sll, maximumIndex, worstIndex, degree},
  timing = AbsoluteTiming[
    states = fgStates[]; qLow = fgLowCodebook[];
    qTelemetry = fgTelemetryCodebook[];
    modal = fgModalData[m, wp];
    angles = N[Range[0, angularSamples - 1] 2 Pi/angularSamples - Pi];
    basis = fgFarBasis[modal, angles];

    lowSystem = fgBuildSystem[states, qLow, modal];
    lowMatrix = fgValue[lowSystem, "Matrix"];
    lowY = fgValue[lowSystem, "YMatrix"];
    lowSolver = LinearSolve[lowMatrix];
    incidenceDegrees = Range[-30, 30, 5];
    lowRatios = {}; lowResiduals = {};
    lowFields = {}; pecFields = {};
    Do[
      alpha = N[degree Degree];
      b = fgPlaneExcitation[modal, alpha];
      u = lowSolver[b]; j = lowY.u;
      field = fgFarFieldWithAngles[modal, basis, angles, j, False];
      jPEC = fgPECcurrent[modal, alpha];
      fieldPEC = fgFarFieldWithAngles[modal, basis, angles, jPEC, False];
      protected = fgProtectedIndices[angles, alpha, 10 Degree];
      denominator = Max[Abs[fieldPEC[[protected]]]^2];
      AppendTo[lowRatios,
       Max[Abs[field[[protected]]]^2]/denominator];
      AppendTo[lowResiduals,
       Norm[lowMatrix.u - b]/Max[Norm[b], 10.^-300]];
      AppendTo[lowFields, field]; AppendTo[pecFields, fieldPEC],
      {degree, incidenceDegrees}];
    worstIndex = First[First[Position[lowRatios, Max[lowRatios]]]];

    telemetrySystem = fgBuildSystem[states, qTelemetry, modal];
    telemetryMatrix = fgValue[telemetrySystem, "Matrix"];
    telemetryY = fgValue[telemetrySystem, "YMatrix"];
    telemetrySolver = LinearSolve[telemetryMatrix];
    bTelemetry = fgValue[modal, "LineExcitation"];
    uTelemetry = telemetrySolver[bTelemetry];
    jTelemetry = telemetryY.uTelemetry;
    telemetryField = fgFarFieldWithAngles[modal, basis, angles,
      jTelemetry, True];
    telemetryResidual = Norm[telemetryMatrix.uTelemetry - bTelemetry]/
      Max[Norm[bTelemetry], 10.^-300];
    power = Abs[telemetryField]^2;
    mainIndices = Select[Range[Length[angles]],
      Abs[fgWrapAngle[angles[[#]]]] <=
        8 Degree + 100 $MachineEpsilon &];
    sideIndices = Select[Range[Length[angles]],
      Abs[fgWrapAngle[angles[[#]]]] >=
        15 Degree - 100 $MachineEpsilon &];
    main = power[[mainIndices]]; side = power[[sideIndices]];
    edgeCorrection = (First[main] + Last[main])/2.;
    efficiency = (Total[main] - edgeCorrection)/Total[power];
    directivity = Max[power]/Mean[power];
    sll = Max[side]/Max[main];
    maximumIndex = First[First[Position[power, Max[power]]]];

    result = {
      "M" -> m, "NumberOfModes" -> 2 m + 1,
      "Angles" -> angles,
      "IncidenceDegrees" -> incidenceDegrees,
      "LowRatios" -> lowRatios,
      "LowResiduals" -> lowResiduals,
      "LowFields" -> lowFields, "PECFields" -> pecFields,
      "WorstSectorDB" -> fgDB10[Max[lowRatios]],
      "MeanSectorDB" -> fgDB10[Mean[lowRatios]],
      "WorstIncidenceDeg" -> incidenceDegrees[[worstIndex]],
      "WorstIncidenceIndex" -> worstIndex,
      "MaximumLowProjectedResidual" -> Max[lowResiduals],
      "TelemetryField" -> telemetryField,
      "Efficiency" -> efficiency,
      "DirectivityDB" -> fgDB10[directivity],
      "SLLDB" -> fgDB10[sll],
      "BeamAngleDeg" -> angles[[maximumIndex]]/Degree,
      "TelemetryProjectedResidual" -> telemetryResidual
      };
    ];
  Append[result, "WallTimeSeconds" -> First[timing]]
  ];

(* ====================================================================== *)
(* SELF-TESTS                                                            *)
(* ====================================================================== *)

FGSelfTest[OptionsPattern[]] := Module[
  {wp, states, qLow, qTelemetry, rows, passed, inputPassed,
   pRange, exactCoefficients, quadratureCoefficients, quadratureCount,
   quadratureAngles, yValues, ySamples, coefficientError, mSmall,
   modesSmall, trialU, yMatrixSmall, projectedProduct,
   quadratureProduct, uSamples, multiplicationError, modal, zConstant,
   constantY, constantMatrix, alpha, b, u, current, exactCurrent,
   constantError, pecM, pecModal, pecCurrent, phi, sourcePoints,
   currentSamples, rho, observationAngles, modalNear, directNear,
   nearError, k, eta, a, ds, obs, distances, actualModal,
   actualSystem, actualMatrix, actualY, actualSolver, actualB,
   actualU, lowResidual, telemetrySystem, telemetryMatrix,
   telemetrySolver, telemetryB, telemetryU, telemetryResidual,
   finitePassed},

  Print["FGSelfTest: starting independent Fourier--Galerkin checks..."];
  wp = OptionValue[FGWorkingPrecision];
  states = fgStates[]; qLow = fgLowCodebook[];
  qTelemetry = fgTelemetryCodebook[];
  rows = {};

  inputPassed = Length[states] == 8 && Length[qLow] == 64 &&
    Length[qTelemetry] == 64 && Min[Re /@ states] >= 0 &&
    Min[qLow] >= 1 && Max[qLow] <= 8 &&
    Min[qTelemetry] >= 1 && Max[qTelemetry] <= 8;
  AppendTo[rows, {"frozen input dimensions and passivity",
    Min[Re /@ states], "minimum Re[Z] >= 0 and correct lengths",
    fgStatus[inputPassed]}];

  (* Exact coefficients versus an independent midpoint quadrature. *)
  pRange = Range[-12, 12];
  yValues = 1/states[[qLow]];
  exactCoefficients = fgPiecewiseCoefficient[yValues, #] & /@ pRange;
  quadratureCount = 8192;
  quadratureAngles = N[(Range[0, quadratureCount - 1] + 1/2)
      2 Pi/quadratureCount];
  ySamples = fgPiecewiseSamples[yValues, quadratureAngles];
  quadratureCoefficients = Table[
    Mean[ySamples Exp[-I p quadratureAngles]], {p, pRange}];
  coefficientError = fgRelativeError[quadratureCoefficients,
    exactCoefficients];
  AppendTo[rows, {"piecewise Fourier coefficients",
    coefficientError, "relative error <= 2*10^-5",
    fgStatus[coefficientError <= 2 10^-5]}];

  (* Toeplitz orientation versus pointwise multiplication and projection. *)
  mSmall = 8; modesSmall = Range[-mSmall, mSmall];
  trialU = N[Exp[-Abs[modesSmall]/3.] (1 + 0.07 I modesSmall)];
  yMatrixSmall = fgMultiplicationMatrix[yValues, mSmall];
  projectedProduct = yMatrixSmall.trialU;
  uSamples = Exp[I Outer[Times, quadratureAngles, modesSmall]].trialU;
  quadratureProduct = Table[
    Mean[ySamples uSamples Exp[-I n quadratureAngles]],
    {n, modesSmall}];
  multiplicationError = fgRelativeError[quadratureProduct,
    projectedProduct];
  AppendTo[rows, {"admittance multiplication matrix orientation",
    multiplicationError, "relative error <= 3*10^-5",
    fgStatus[multiplicationError <= 3 10^-5]}];

  (* Constant impedance has an independent closed-form modal solution. *)
  modal = fgModalData[18, wp];
  zConstant = 0.7 + 0.3 I;
  constantY = (1/zConstant) IdentityMatrix[37];
  constantMatrix = IdentityMatrix[37] +
    DiagonalMatrix[fgValue[modal, "SingleLayer"]].constantY;
  alpha = 0.37; b = fgPlaneExcitation[modal, alpha];
  u = LinearSolve[constantMatrix, b]; current = constantY.u;
  exactCurrent = b/(zConstant + fgValue[modal, "SingleLayer"]);
  constantError = fgRelativeError[current, exactCurrent];
  AppendTo[rows, {"constant-impedance closed-form solution",
    constantError, "relative error <= 10^-11",
    fgStatus[constantError <= 10^-11]}];

  (* Modal single-layer field versus direct physical-space quadrature. *)
  pecM = 28; pecModal = fgModalData[pecM, wp];
  alpha = -0.29; pecCurrent = fgPECcurrent[pecModal, alpha];
  k = fgValue[pecModal, "k"]; eta = fgValue[pecModal, "eta"];
  a = fgValue[pecModal, "a"]; rho = 5.;
  observationAngles = {-0.9, -0.2, 0.4, 1.1};
  modalNear = Table[
    -k eta/4 2 Pi a Total[Table[
       HankelH1[n, k rho] BesselJ[n, k a]
        pecCurrent[[n + pecM + 1]] Exp[I n theta],
       {n, -pecM, pecM}]], {theta, observationAngles}];
  quadratureCount = 4096;
  phi = N[(Range[0, quadratureCount - 1] + 1/2)
      2 Pi/quadratureCount];
  sourcePoints = Transpose[{a Cos[phi], a Sin[phi]}];
  currentSamples = Exp[I Outer[Times, phi,
       fgValue[pecModal, "Indices"]]].pecCurrent;
  ds = 2 Pi a/quadratureCount;
  directNear = Table[
    obs = {rho Cos[theta], rho Sin[theta]};
    distances = Sqrt[(#[[1]]^2 + #[[2]]^2) & /@
       (sourcePoints - ConstantArray[obs, quadratureCount])];
    -k eta/4 ds Total[HankelH1[0, k distances] currentSamples],
    {theta, observationAngles}];
  nearError = fgRelativeError[directNear, modalNear];
  AppendTo[rows, {"modal field versus direct boundary quadrature",
    nearError, "relative error <= 10^-6",
    fgStatus[nearError <= 10^-6]}];

  (* Actual frozen systems: projected residuals must reach roundoff. *)
  actualModal = fgModalData[32, wp];
  actualSystem = fgBuildSystem[states, qLow, actualModal];
  actualMatrix = fgValue[actualSystem, "Matrix"];
  actualY = fgValue[actualSystem, "YMatrix"];
  actualSolver = LinearSolve[actualMatrix];
  actualB = fgPlaneExcitation[actualModal, 15 Degree];
  actualU = actualSolver[actualB];
  lowResidual = Norm[actualMatrix.actualU - actualB]/Norm[actualB];
  AppendTo[rows, {"frozen low-scattering projected residual",
    lowResidual, "relative residual <= 10^-11",
    fgStatus[lowResidual <= 10^-11]}];

  telemetrySystem = fgBuildSystem[states, qTelemetry, actualModal];
  telemetryMatrix = fgValue[telemetrySystem, "Matrix"];
  telemetrySolver = LinearSolve[telemetryMatrix];
  telemetryB = fgValue[actualModal, "LineExcitation"];
  telemetryU = telemetrySolver[telemetryB];
  telemetryResidual = Norm[telemetryMatrix.telemetryU - telemetryB]/
    Norm[telemetryB];
  AppendTo[rows, {"frozen telemetry projected residual",
    telemetryResidual, "relative residual <= 10^-11",
    fgStatus[telemetryResidual <= 10^-11]}];

  finitePassed = And @@ (fgFiniteNumberQ[#[[2]]] & /@ rows);
  AppendTo[rows, {"all reported test values finite", Boole[finitePassed],
    "value = 1", fgStatus[finitePassed]}];
  passed = And @@ (Last[#] == "PASS" & /@ rows);

  Print[TableForm[rows,
    TableHeadings -> {None, {"check", "measured value", "criterion",
       "status"}}]];
  Print[If[passed, "FGSelfTest PASSED.", "FGSelfTest FAILED."]];
  {"Passed" -> passed, "Checks" -> rows}
  ];

(* ====================================================================== *)
(* OUTPUT TABLES, FIGURES, AND ACCEPTANCE CHECKS                          *)
(* ====================================================================== *)

fgConvergenceRows[results_List] := Table[
  {fgValue[r, "M"], fgValue[r, "NumberOfModes"],
   fgValue[r, "WorstSectorDB"], fgValue[r, "MeanSectorDB"],
   fgValue[r, "WorstIncidenceDeg"], fgValue[r, "Efficiency"],
   fgValue[r, "DirectivityDB"], fgValue[r, "SLLDB"],
   fgValue[r, "BeamAngleDeg"],
   fgValue[r, "MaximumLowProjectedResidual"],
   fgValue[r, "TelemetryProjectedResidual"],
   fgValue[r, "WallTimeSeconds"]}, {r, results}];

fgAcceptanceRows[results_List] := Module[
  {previous, final, rows, difference, test},
  previous = results[[-2]]; final = Last[results]; rows = {};
  test[name_, key_, tolerance_, unit_] := (
    difference = Abs[fgValue[final, key] - fgValue[previous, key]];
    AppendTo[rows, {name, difference, tolerance, unit,
      fgStatus[difference <= tolerance]}]);
  test["final-two-order worst-sector change", "WorstSectorDB", 0.5,
    "dB"];
  test["final-two-order mean-sector change", "MeanSectorDB", 0.5,
    "dB"];
  test["final-two-order telemetry directivity change", "DirectivityDB",
    0.2, "dB"];
  test["final-two-order telemetry sidelobe change", "SLLDB", 0.5,
    "dB"];
  test["final-two-order beam-angle change", "BeamAngleDeg", 0.5,
    "deg"];
  test["final-two-order main-sector-efficiency change", "Efficiency",
    0.01, "absolute"];
  AppendTo[rows, {"maximum final low-scattering projected residual",
    fgValue[final, "MaximumLowProjectedResidual"], 10^-10, "relative",
    fgStatus[fgValue[final, "MaximumLowProjectedResidual"] <= 10^-10]}];
  AppendTo[rows, {"final telemetry projected residual",
    fgValue[final, "TelemetryProjectedResidual"], 10^-10, "relative",
    fgStatus[fgValue[final, "TelemetryProjectedResidual"] <= 10^-10]}];
  rows
  ];

fgCrossSolverRows[final_List] := Module[
  {mom256, mom384, specification},
  mom256 = fgMoM256[]; mom384 = fgMoM384[];
  specification = {
    {"worst protected-sector power ratio", "WorstSectorDB", "dB", 1.0},
    {"mean protected-sector power ratio", "MeanSectorDB", "dB", 1.0},
    {"telemetry main-sector efficiency", "Efficiency", "absolute", 0.04},
    {"telemetry two-dimensional directivity", "DirectivityDB", "dB", .25},
    {"telemetry peak sidelobe level", "SLLDB", "dB", .75},
    {"telemetry beam angle", "BeamAngleDeg", "deg", 1.0}
    };
  Table[
   {row[[1]], row[[3]], fgValue[mom256, row[[2]]],
    fgValue[mom384, row[[2]]], fgValue[final, row[[2]]],
    fgValue[final, row[[2]]] - fgValue[mom256, row[[2]]],
    row[[4]], fgStatus[
      Abs[fgValue[final, row[[2]]] - fgValue[mom256, row[[2]]]] <=
       row[[4]] + 10^-12]},
   {row, specification}]
  ];

fgStateRows[] := Module[{states, qLow, qTelemetry, delta},
  states = fgStates[]; qLow = fgLowCodebook[];
  qTelemetry = fgTelemetryCodebook[]; delta = 360./64;
  Table[
   {cell, (cell - 1) delta, (cell - 1/2) delta, cell delta,
    qLow[[cell]], Re[states[[qLow[[cell]]]]],
    Im[states[[qLow[[cell]]]]], qTelemetry[[cell]],
    Re[states[[qTelemetry[[cell]]]]],
    Im[states[[qTelemetry[[cell]]]]]},
   {cell, 1, 64}]
  ];

fgMakeFigures[results_List, final_List] := Module[
  {orders, convergence, p1, p2, fig1, angles, angleDegrees,
   worstIndex, incidenceDegrees, alpha, lowField, pecField, protected,
   normalization, lowDB, pecDB, pLow, telemetryField, power,
   normalizedPower, freePowerDB, polarData, freeData, pTelemetry,
   k, eta, figures},
  orders = fgValue[#, "M"] & /@ results;
  convergence = fgConvergenceRows[results];
  p1 = ListLinePlot[
    {convergence[[All, {1, 3}]], convergence[[All, {1, 4}]]},
    Joined -> True, PlotMarkers -> Automatic, Frame -> True,
    FrameLabel -> {"Fourier half-order M", "power ratio (dB)"},
    PlotLabel -> Style["(a) Low-scattering convergence", 12],
    PlotStyle -> {{Blue, Thick}, {Darker[Green], Thick}},
    PlotLegends -> Placed[{"worst protected sector", "mean protected sector"},
      Below], GridLines -> Automatic,
    BaseStyle -> {FontFamily -> "Times", 10}];
  p2 = ListLinePlot[
    {convergence[[All, {1, 7}]], convergence[[All, {1, 8}]]},
    Joined -> True, PlotMarkers -> Automatic, Frame -> True,
    FrameLabel -> {"Fourier half-order M", "metric (dB)"},
    PlotLabel -> Style["(b) Telemetry convergence", 12],
    PlotStyle -> {{Blue, Thick}, {Orange, Thick}},
    PlotLegends -> Placed[{"2-D directivity", "peak sidelobe level"},
      Below], GridLines -> Automatic,
    BaseStyle -> {FontFamily -> "Times", 10}];
  fig1 = GraphicsGrid[{{p1, p2}}, ImageSize -> 900, Spacings -> .35];

  angles = fgValue[final, "Angles"];
  angleDegrees = angles/Degree;
  worstIndex = fgValue[final, "WorstIncidenceIndex"];
  incidenceDegrees = fgValue[final, "IncidenceDegrees"];
  alpha = incidenceDegrees[[worstIndex]] Degree;
  lowField = fgValue[final, "LowFields"][[worstIndex]];
  pecField = fgValue[final, "PECFields"][[worstIndex]];
  protected = fgProtectedIndices[angles, alpha, 10 Degree];
  normalization = Max[Abs[pecField[[protected]]]^2];
  lowDB = (fgDB10[#/normalization] &) /@ (Abs[lowField]^2);
  pecDB = (fgDB10[#/normalization] &) /@ (Abs[pecField]^2);
  pLow = ListLinePlot[
    {Transpose[{angleDegrees, Max[-50., #] & /@ pecDB}],
     Transpose[{angleDegrees, Max[-50., #] & /@ lowDB}]},
    Joined -> True, PlotRange -> {{-180, 180}, {-50, 15}},
    Frame -> True,
    FrameLabel -> {"observation angle theta (deg)",
      "power / protected-sector PEC peak (dB)"},
    PlotLabel -> Style[StringJoin[
       "Independent low-scattering pattern, alpha = ",
       ToString[incidenceDegrees[[worstIndex]]], " deg"], 12],
    PlotStyle -> {{GrayLevel[.35], Dashed, Thick}, {Blue, Thick}},
    PlotLegends -> Placed[{"PEC reference", "Fourier--Galerkin surface"},
      Below], GridLines -> Automatic,
    BaseStyle -> {FontFamily -> "Times", 10}, ImageSize -> 760];

  telemetryField = fgValue[final, "TelemetryField"];
  power = Abs[telemetryField]^2;
  normalizedPower = power/Max[power];
  k = 2 Pi; eta = 1.;
  freePowerDB = fgDB10[(k eta/4)^2/Max[power]];
  polarData = Transpose[{angles,
      40 + (Max[-40., #] & /@ (fgDB10 /@ normalizedPower))}];
  freeData = Transpose[{angles,
      ConstantArray[40 + Max[-40., freePowerDB], Length[angles]]}];
  pTelemetry = ListLinePlot[
    {fgPolarCartesian[polarData], fgPolarCartesian[freeData]},
    Joined -> True, PlotRange -> {{-47, 47}, {-47, 47}},
    AspectRatio -> 1, Axes -> False, Frame -> False,
    Prolog -> fgPolarGridLines[40,
      {{0, "-40"}, {10, "-30"}, {20, "-20"},
       {30, "-10"}, {40, "0"}}],
    Epilog -> fgPolarGridLabels[40,
      {{0, "-40"}, {10, "-30"}, {20, "-20"},
       {30, "-10"}, {40, "0"}}],
    PlotStyle -> {{Blue, Thick}, {GrayLevel[.35], Dashed, Thick}},
    PlotLabel -> Style["Independent normalized telemetry pattern", 12],
    PlotLegends -> Placed[{"Fourier--Galerkin surface",
       "free-space line source"}, Below],
    BaseStyle -> {FontFamily -> "Times", 10}, ImageSize -> 620];
  figures = {fig1, pLow, pTelemetry};
  figures
  ];

(* ====================================================================== *)
(* COMPLETE INDEPENDENT VALIDATION                                        *)
(* ====================================================================== *)

FGRunValidation[OptionsPattern[]] := Module[
  {root, dataDirectory, figureDirectory, orders, angularSamples, wp,
   selfTest, results, result, final, convergenceRows, acceptanceRows,
   crossRows, incidenceRows, telemetryRows, lowPatternRows, angles,
   angleDegrees, lowRatios, lowResiduals, worstIndex, alpha,
   protected, lowField, pecField, normalization, telemetryField,
   telemetryPower, telemetryPeak, mainFlags, sideFlags, stateRows,
   figures, figureFiles, summaryLines, allAcceptancePassed,
   crossPassed},

  root = OptionValue[FGValidationDirectory];
  If[root === Automatic,
   root = FileNameJoin[{$FGPackageDirectory,
       "fourier_galerkin_results"}]];
  dataDirectory = FileNameJoin[{root, "data"}];
  figureDirectory = FileNameJoin[{root, "figures"}];
  fgEnsureDirectory[dataDirectory]; fgEnsureDirectory[figureDirectory];

  orders = Sort[DeleteDuplicates[OptionValue[FGFourierOrders]]];
  angularSamples = OptionValue[FGAngularSamples];
  wp = OptionValue[FGWorkingPrecision];
  If[Length[orders] < 2 || ! And @@ (IntegerQ[#] && # >= 16 & /@ orders),
   Print["FGRunValidation stopped: provide at least two integer Fourier orders M >= 16."];
   Return[$Failed]];
  If[! IntegerQ[angularSamples] || angularSamples < 360,
   Print["FGRunValidation stopped: FGAngularSamples must be an integer >= 360."];
   Return[$Failed]];

  If[TrueQ[OptionValue[FGRunSelfTest]],
   selfTest = FGSelfTest[FGWorkingPrecision -> wp];
   If[! TrueQ[fgValue[selfTest, "Passed"]],
    Print["FGRunValidation stopped because FGSelfTest did not pass."];
    Return[$Failed]]];

  Print["FGRunValidation: the original DualModeRIS.m package is not loaded."];
  Print["FGRunValidation: retained half-orders M = ", orders, "."];
  results = {};
  Do[
   Print["  Solving M = ", order, " (", 2 order + 1,
    " complex Fourier modes)..."];
   result = fgOneOrderStudy[order, angularSamples, wp];
   AppendTo[results, result];
   Print["    worst sector = ", NumberForm[fgValue[result,
       "WorstSectorDB"], {7, 3}], " dB; telemetry D_2D = ",
    NumberForm[fgValue[result, "DirectivityDB"], {7, 3}], " dB; ",
    "wall time = ", NumberForm[fgValue[result, "WallTimeSeconds"],
      {8, 2}], " s."],
   {order, orders}];
  final = Last[results];

  convergenceRows = fgConvergenceRows[results];
  fgExportCSV[FileNameJoin[{dataDirectory,
      "fg_validation_convergence.csv"}],
    {"fourier_half_order_M", "number_of_complex_modes_2M_plus_1",
     "worst_protected_sector_power_db_relative_to_pec",
     "mean_protected_sector_power_db_relative_to_pec",
     "incidence_angle_of_worst_sector_deg",
     "telemetry_main_sector_efficiency",
     "telemetry_two_dimensional_directivity_db",
     "telemetry_peak_sidelobe_level_db",
     "telemetry_beam_angle_deg",
     "maximum_low_scattering_projected_relative_residual",
     "telemetry_projected_relative_residual", "wall_time_seconds"},
    convergenceRows];

  acceptanceRows = fgAcceptanceRows[results];
  fgExportCSV[FileNameJoin[{dataDirectory,
      "fg_numerical_acceptance_checks.csv"}],
    {"check", "measured_absolute_change_or_residual", "tolerance",
     "unit", "status"}, acceptanceRows];

  crossRows = fgCrossSolverRows[final];
  fgExportCSV[FileNameJoin[{dataDirectory,
      "fg_cross_solver_comparison.csv"}],
    {"metric", "unit", "published_midpoint_mom_256",
     "published_midpoint_mom_384", "fourier_galerkin_final",
     "signed_difference_fg_minus_mom256", "declared_audit_tolerance",
     "within_declared_tolerance"}, crossRows];

  angles = fgValue[final, "Angles"];
  angleDegrees = angles/Degree;
  lowRatios = fgValue[final, "LowRatios"];
  lowResiduals = fgValue[final, "LowResiduals"];
  incidenceRows = Table[
    {fgValue[final, "IncidenceDegrees"][[q]], lowRatios[[q]],
     fgDB10[lowRatios[[q]]], lowResiduals[[q]]},
    {q, Length[lowRatios]}];
  fgExportCSV[FileNameJoin[{dataDirectory,
      "fg_final_low_scattering_incidence_sweep.csv"}],
    {"incidence_angle_deg", "worst_protected_sector_power_ratio_linear",
     "worst_protected_sector_power_ratio_db_relative_to_pec",
     "projected_relative_residual"}, incidenceRows];

  worstIndex = fgValue[final, "WorstIncidenceIndex"];
  alpha = fgValue[final, "IncidenceDegrees"][[worstIndex]] Degree;
  protected = fgProtectedIndices[angles, alpha, 10 Degree];
  lowField = fgValue[final, "LowFields"][[worstIndex]];
  pecField = fgValue[final, "PECFields"][[worstIndex]];
  normalization = Max[Abs[pecField[[protected]]]^2];
  lowPatternRows = Table[
    {angleDegrees[[q]], Re[lowField[[q]]], Im[lowField[[q]]],
     Abs[lowField[[q]]]^2,
     fgDB10[Abs[lowField[[q]]]^2/normalization],
     Abs[pecField[[q]]]^2,
     fgDB10[Abs[pecField[[q]]]^2/normalization],
     Boole[MemberQ[protected, q]]}, {q, Length[angles]}];
  fgExportCSV[FileNameJoin[{dataDirectory,
      "fg_final_worst_incidence_scattering_pattern.csv"}],
    {"observation_angle_deg", "fg_scattered_field_real",
     "fg_scattered_field_imag", "fg_scattered_power_linear",
     "fg_scattered_power_db_relative_to_protected_sector_pec_peak",
     "pec_scattered_power_linear",
     "pec_scattered_power_db_relative_to_protected_sector_pec_peak",
     "inside_protected_sector_1_yes_0_no"}, lowPatternRows];

  telemetryField = fgValue[final, "TelemetryField"];
  telemetryPower = Abs[telemetryField]^2;
  telemetryPeak = Max[telemetryPower];
  mainFlags = Boole[Abs[fgWrapAngle[#]] <=
       8 Degree + 100 $MachineEpsilon] & /@ angles;
  sideFlags = Boole[Abs[fgWrapAngle[#]] >=
       15 Degree - 100 $MachineEpsilon] & /@ angles;
  telemetryRows = Table[
    {angleDegrees[[q]], Re[telemetryField[[q]]],
     Im[telemetryField[[q]]], telemetryPower[[q]],
     telemetryPower[[q]]/telemetryPeak,
     fgDB10[telemetryPower[[q]]/telemetryPeak], mainFlags[[q]],
     sideFlags[[q]]}, {q, Length[angles]}];
  fgExportCSV[FileNameJoin[{dataDirectory,
      "fg_final_telemetry_pattern.csv"}],
    {"observation_angle_deg", "total_field_real_including_direct_feed",
     "total_field_imag_including_direct_feed", "total_power_linear",
     "normalized_power_linear", "normalized_power_db",
     "inside_main_sector_1_yes_0_no",
     "inside_sidelobe_region_1_yes_0_no"}, telemetryRows];

  stateRows = fgStateRows[];
  fgExportCSV[FileNameJoin[{dataDirectory,
      "fg_frozen_inputs_and_codebooks.csv"}],
    {"physical_cell_index", "cell_start_angle_deg",
     "cell_center_angle_deg", "cell_end_angle_deg",
     "low_scattering_state_index", "low_scattering_re_z_over_eta",
     "low_scattering_im_z_over_eta", "telemetry_state_index",
     "telemetry_re_z_over_eta", "telemetry_im_z_over_eta"}, stateRows];

  figures = fgMakeFigures[results, final];
  figureFiles = {
    FileNameJoin[{figureDirectory, "fg_fig1_modal_convergence.pdf"}],
    FileNameJoin[{figureDirectory,
      "fg_fig2_worst_incidence_scattering_pattern.pdf"}],
    FileNameJoin[{figureDirectory, "fg_fig3_telemetry_pattern.pdf"}]
    };
  Do[Export[figureFiles[[q]], figures[[q]]], {q, Length[figures]}];

  allAcceptancePassed = And @@ (Last[#] == "PASS" & /@ acceptanceRows);
  crossPassed = And @@ (Last[#] == "PASS" & /@ crossRows);
  summaryLines = {
    "INDEPENDENT FOURIER--GALERKIN VALIDATION SUMMARY",
    "",
    StringJoin["Final half-order M: ", ToString[fgValue[final, "M"]]],
    StringJoin["Worst protected-sector ratio (dB): ",
      ToString[fgValue[final, "WorstSectorDB"], InputForm]],
    StringJoin["Mean protected-sector ratio (dB): ",
      ToString[fgValue[final, "MeanSectorDB"], InputForm]],
    StringJoin["Worst sampled incidence (deg): ",
      ToString[fgValue[final, "WorstIncidenceDeg"]]],
    StringJoin["Telemetry main-sector efficiency: ",
      ToString[fgValue[final, "Efficiency"], InputForm]],
    StringJoin["Telemetry 2-D directivity (dB): ",
      ToString[fgValue[final, "DirectivityDB"], InputForm]],
    StringJoin["Telemetry peak sidelobe level (dB): ",
      ToString[fgValue[final, "SLLDB"], InputForm]],
    StringJoin["Telemetry beam angle (deg): ",
      ToString[fgValue[final, "BeamAngleDeg"], InputForm]],
    StringJoin["Order-convergence and residual checks: ",
      If[allAcceptancePassed, "PASS", "FAIL"]],
    StringJoin["Cross-solver comparison within declared audit tolerances: ",
      If[crossPassed, "PASS", "FAIL"]],
    "",
    "The Fourier--Galerkin result is an independent formulation. Exact equality with the midpoint MoM values is neither imposed nor expected."
    };
  Export[FileNameJoin[{root, "fg_validation_summary.txt"}],
    StringJoin[Riffle[summaryLines, "\n"]], "Text"];

  Print["FGRunValidation: numerical acceptance checks"];
  Print[TableForm[acceptanceRows,
    TableHeadings -> {None, {"check", "measured", "tolerance",
       "unit", "status"}}]];
  Print["FGRunValidation: cross-solver comparison"];
  Print[TableForm[crossRows,
    TableHeadings -> {None, {"metric", "unit", "MoM 256", "MoM 384",
       "FG final", "FG-MoM256", "tolerance", "status"}}]];
  Print["FGRunValidation: results written to ", root];
  If[TrueQ[OptionValue[FGDisplayFigures]], Scan[Print, figures]];

  {"PassedNumericalChecks" -> allAcceptancePassed,
   "WithinCrossSolverAuditTolerances" -> crossPassed,
   "ResultDirectory" -> root, "FinalResult" -> final,
   "AcceptanceChecks" -> acceptanceRows,
   "CrossSolverComparison" -> crossRows}
  ];

End[];
EndPackage[];

