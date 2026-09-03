#!/usr/bin/env python3
"""
verify_files.py - transfer-integrity check for the Paper 1 Code & Data repo.

Run from anywhere:  python verification/verify_files.py
Confirms that (1) every expected script/function/model/data/R input is present, and
(2) every model/data path referenced by the MATLAB scripts and R notebooks resolves
inside the repo. Exits non-zero if anything is missing.
"""
import os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)                      # repo root = parent of verification/
ok, problems = [], []

def check(relpath, note=""):
    p = os.path.join(ROOT, relpath)
    if os.path.exists(p):
        ok.append(relpath)
    else:
        problems.append(f"MISSING: {relpath}" + (f"  ({note})" if note else ""))

# ---- 1. Entry scripts ----
ENTRY = [
    "matlab/scripts/01_HSA_target_identification/HSA_model_Generation_Verification.m",
    "matlab/scripts/02_common_target_identification/common_target_identification_validation.m",
    "matlab/scripts/03_transcription_factor_targets/tf_target_analysis.m",
    "matlab/scripts/04_nutrient_cost/metCostDetermination.m",
    "matlab/scripts/04_nutrient_cost/metCostKOBenchmark.m",
]
for f in ENTRY: check(f)

# ---- 2. Colocated MATLAB inputs ----
for f in [
    "matlab/scripts/01_HSA_target_identification/HSA_carbonScreen_metadata.xlsx",
    "matlab/scripts/02_common_target_identification/carbonScreen_metadata.xlsx",
    "matlab/scripts/02_common_target_identification/commonTargets.csv",
    "matlab/scripts/03_transcription_factor_targets/all_tf_targets.xlsx",
]: check(f)
for p in ["HSA","RBDJ","nb3B2","trastuzumab"]:
    check(f"matlab/scripts/02_common_target_identification/FSEOF_results_{p}.csv")

# ---- 3. Function dependency closure (expect 37 .m files) ----
EXPECTED_FUNCS = 17
fdir = os.path.join(ROOT, "matlab", "functions")
nfun = sum(1 for _,_,fs in os.walk(fdir) for x in fs if x.endswith(".m"))
if nfun == EXPECTED_FUNCS:
    ok.append(f"matlab/functions: {nfun} .m files")
else:
    problems.append(f"function count = {nfun}, expected {EXPECTED_FUNCS}")
for f in ["modelBuilding/generateSynthesisReactions4.m",
          "modelBuilding/generatemAbSynthesisReactions.m",
          "modelBuilding/addTranslationReaction.m",
          "modelBuilding/addDSBReactions.m",
          "modelBuilding/addNglyGolgiReactions.m",
          "modelBuilding/addFoldingReaction.m",
          "modelBuilding/addInterchainDSBReaction.m"]:
    check("matlab/functions/"+f)

# ---- 4. Models & data ----
for f in ["models/MODEL1612130000.xml","models/mbt212871sup0003AppendixS3.xml","models/iMT1026v3.mat",
          "data/baselineModifications.xlsx","data/ppPSIMv3_RP.csv","data/geneConversion.csv","data/essentiality.csv",
          "setup_solver.m","README.md","LICENSE","CITATION.cff",".gitignore"]:
    check(f)

# ---- 5. R notebooks + figure_data ----
for f in ["R/Figs.Rmd","R/geneTargetAnalysis.Rmd","R/03_TF_Graphs.Rmd",
          "R/figure_data/goSLIMmod.gmt","R/figure_data/StrainCultivationResultCompilation.xlsx",
          "R/figure_data/commonTargets.csv","R/figure_data/all_tf_targets.xlsx",
          "R/figure_data/01 Round 1 HSA Upregulations/FSEOF_results.xlsx",
          "R/figure_data/02 Round 2 - Multimolecule FSEOF/Multimolecule_FSEOF_Results.xlsx",
          "R/figure_data/03 Transcription Factors/TranscriptionFactorResults.xlsx",
          "R/figure_data/04 High Cost KO/ppPSIMv3_modified_withCost.csv"]:
    check(f)

# ---- 6. Referenced-path resolution: MATLAB fullfile(repoRoot,'models'|'data','X') ----
ref_re = re.compile(r"fullfile\(repoRoot,'(models|data)','([^']+)'\)")
for s in ENTRY:
    sp = os.path.join(ROOT, s)
    if not os.path.exists(sp): continue
    for sub, name in ref_re.findall(open(sp, encoding="utf-8", errors="replace").read()):
        rel = f"{sub}/{name}"
        if not os.path.exists(os.path.join(ROOT, rel)):
            problems.append(f"{s} references {rel} but it is missing")

# ---- 7. Referenced-path resolution: R figure_data/... ----
rpath_re = re.compile(r'"(figure_data/[^"]+)"')
for s in ["R/Figs.Rmd","R/geneTargetAnalysis.Rmd","R/03_TF_Graphs.Rmd"]:
    sp = os.path.join(ROOT, s)
    if not os.path.exists(sp): continue
    for rel in set(rpath_re.findall(open(sp, encoding="utf-8", errors="replace").read())):
        if not os.path.exists(os.path.join(ROOT, "R", rel)):
            problems.append(f"{s} references R/{rel} but it is missing")

# ---- 8. Leftover non-portable references (should be none) ----
import glob
for s in ENTRY:
    txt = open(os.path.join(ROOT,s), encoding="utf-8", errors="replace").read()
    if re.search(r"[A-Za-z]:\\\\|gurobi_setup", txt) and "setup_solver" not in txt.split("gurobi",1)[0][-200:]:
        pass
    if re.search(r"[A-Za-z]:\\\\", txt):
        problems.append(f"{s} still contains an absolute Windows path")

# ---- report ----
print(f"Repo root: {ROOT}")
print(f"Checks passed: {len(ok)}")
for o in ok: print("  OK  ", o)
print()
if problems:
    print(f"PROBLEMS ({len(problems)}):")
    for p in problems: print("  X  ", p)
    sys.exit(1)
print("ALL CHECKS PASSED - transfer looks complete.")
