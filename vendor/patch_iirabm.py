#!/usr/bin/env python3
"""Reproducible port of the upstream IIRABM (cockrell/iirabm_public V1) to a
portable *serial* single-run build.

The upstream code targets Cray/GCC + the deprecated MPI C++ bindings. This
patch:
  (1) adds <random> where Apple clang needs it and drops the mpi.h includes,
  (2) converts a small set of run-mode / timing CONSTANTS (compile-time in
      Parameters.h) into runtime globals, so one binary covers baseline,
      intervention-replay and variable schedules,
  (3) replaces the MPI orchestration main()/postSimulationData() with a serial
      single-run driver that reads params/seed/intervention from argv and
      writes the full trajectory to a file.

Physics (mainSimulation loop, cell/sim/agent functions) is byte-for-byte
upstream. Parallelism is done at the Python bridge level (one process per
replicate), which matches how the paper's MPI build parallelised: each rank was
an independent simulation.

Usage: patch_iirabm.py <pristine_V1_dir> <out_dir>
"""
import re, sys, shutil, pathlib

src = pathlib.Path(sys.argv[1])
out = pathlib.Path(sys.argv[2])
if out.exists():
    shutil.rmtree(out)
out.mkdir(parents=True)
for f in src.glob("*.cpp"):
    shutil.copy(f, out / f.name)
for f in src.glob("*.h"):
    shutil.copy(f, out / f.name)

def read(p):  return (out / p).read_text()
def write(p, s): (out / p).write_text(s)

# ---------------------------------------------------------------- Parameters.h
# Convert run-mode / timing constants to extern runtime globals (defined in the
# serial driver). Keep array-sizing constants compile-time.
params = read("Parameters.h")
RUNTIME_INT = [
    "numTimeSteps", "numReseeds", "reseedTime", "GA", "reseedPostInt",
    "interventionTimeStart", "interventionTimeSeperation", "interventionLength",
    "numIndependentInterventions", "testIntervention", "testStart", "numTests",
]
for name in RUNTIME_INT:
    # drop the `const int NAME = ...;` definition; replace with extern decl
    params, n = re.subn(
        r'const\s+int\s+' + name + r'\s*=\s*[^;]+;',
        f'extern int {name};',
        params, count=1)
    assert n == 1, f"failed to externalize {name}"
# intervention file path -> runtime extern string
params, n = re.subn(r'const\s+string\s+intFile\s*=\s*[^;]+;',
                    'extern string intFile;', params, count=1)
assert n == 1, "failed to externalize intFile"
# Make the output trajectory arrays big enough for long sims regardless of the
# compile-time numWriteSteps (we size by numWriteSteps, used as an array bound).
params, n = re.subn(r'const\s+int\s+numWriteSteps\s*=\s*[^;]+;',
                    'const int numWriteSteps=15000;', params, count=1)
assert n == 1
# genomeLength / testLength sized for up to 8 intervention stages (8*12=96)
params = re.sub(r'const\s+int\s+genomeLength\s*=\s*[^;]+;',
                'const int genomeLength=96;', params, count=1)
params = re.sub(r'const\s+int\s+testLength\s*=\s*[^;]+;',
                'const int testLength=96;', params, count=1)
# numIndependentInterventions assertion const removed above (now runtime)
write("Parameters.h", params)

# ---------------------------------------------- add <random> / drop mpi.h
def ensure_random(fn):
    s = read(fn)
    if "#include <random>" not in s:
        s = s.replace('#include "agents.h"',
                      '#include <random>\n#include "agents.h"', 1)
    write(fn, s)

for fn in ["cellFunctions.cpp", "simulationFunctions.cpp", "gaFunctions.cpp"]:
    ensure_random(fn)

# gaFunctions: drop mpi.h, add extern gaGenerator
s = read("gaFunctions.cpp")
s = s.replace('#include "mpi.h"\n', '')
if "extern mt19937 gaGenerator;" not in s:
    s = s.replace('using namespace std;',
                  'using namespace std;\nextern mt19937 gaGenerator;', 1)
write("gaFunctions.cpp", s)

# ---------------------------------------------------------- IIR_Sepsis.cpp
s = read("IIR_Sepsis.cpp")
# drop mpi include, add headers we need
s = s.replace('#include "mpi.h"\n', '#include <cstring>\n#include <string>\n#include <sstream>\n')

# Define the runtime globals (formerly Parameters.h consts) once, here.
runtime_defs = """
// ---- runtime globals (ported from compile-time Parameters.h consts) ----
int numTimeSteps=10000;
int numReseeds=0;
int reseedTime=505;
int GA=0;
int reseedPostInt=0;
int interventionTimeStart=505;
int interventionTimeSeperation=100;
int interventionLength=99;
int numIndependentInterventions=0;
int testIntervention=0;
int testStart=205;
int numTests=0;
string intFile="";
// ------------------------------------------------------------------------
"""
# insert the runtime defs right after the `using namespace std;` line
s = s.replace("using namespace std;\n",
              "using namespace std;\n" + runtime_defs, 1)

# The upstream `extern const int ... numTimeSteps,numWriteSteps, ...` line now
# conflicts (numTimeSteps is no longer const). Fix that extern declaration.
s = s.replace(
    "extern const int cellCapacity,xDim,yDim,injuryStep,parameterInput,numTimeSteps,numWriteSteps,\n\tga_popSize,maxIter;",
    "extern const int cellCapacity,xDim,yDim,injuryStep,parameterInput,numWriteSteps,\n\tga_popSize,maxIter;\nextern int numTimeSteps;")

# ---- excise the MPI main() and postSimulationData(), keep everything else ----
# main() spans from 'int main(){' to its closing. We delete main() and
# postSimulationData() by regex between known anchors.
# Delete main(): from 'int main(){' up to (but not including) the line
# 'void mainSimulation('
m_start = s.index("int main(){")
# anchor on the DEFINITION of mainSimulation (2nd occurrence), not its
# forward declaration near the top of the file.
first = s.index("void mainSimulation(")
ms_anchor = s.index("void mainSimulation(", first + 1)
assert ms_anchor > m_start, "mainSimulation definition must follow main()"
s = s[:m_start] + s[ms_anchor:]

# Delete postSimulationData(): from 'void postSimulationData(int index){' up to
# 'void printIntervention(){'
pd_start = s.index("void postSimulationData(int index){")
pi_anchor = s.index("void printIntervention(){")
s = s[:pd_start] + s[pi_anchor:]

# mainSimulation still *declares* `void postSimulationData(int index);` at top &
# calls it? No — main called it. mainSimulation doesn't. Remove dangling proto.
s = s.replace("void postSimulationData(int index);\n", "")

# Append the serial driver (new main + writer) at the end of the file.
serial_driver = r'''

// ======================= SERIAL SINGLE-RUN DRIVER =======================
// Replaces the upstream MPI orchestration. Runs ONE simulation with params /
// seed / intervention taken from argv and writes the full trajectory to a file.
// Output format (whitespace-separated), matching upstream allOutputWriter:
//   line 1 : inj oxyHeal infectSpread numRecurInj numInfectRepeat seed nSteps
//   then 21 series lines, each nSteps values:
//     pmn mono TH1 TH2 oxy infect TNF sTNFr IL10 GCSF pTH1 pTH2
//     PAF IFNg IL4 IL8 IL12 sIL1r IL1ra IL1
// A dead patient's oxy series is -1 from the death step onward (upstream rule).

extern void initialize();

static void writeSeries(ofstream& f, float* a, int n){
    for(int j=0;j<n;j++){ f<<a[j]<<" "; } f<<"\n";
}
static void writeSeriesI(ofstream& f, int* a, int n){
    for(int j=0;j<n;j++){ f<<a[j]<<" "; } f<<"\n";
}

int main(int argc, char** argv){
    procID=0; numProcs=1;
    // defaults = paper's trained parameter set (I2 H0.1 T5 E2 inj33)
    int   inj      = (argc>1)? atoi(argv[1]) : 33;
    float oxyHeal  = (argc>2)? atof(argv[2]) : 0.1f;
    int   toxi     = (argc>3)? atoi(argv[3]) : 5;    // infectSpread
    int   envTox   = (argc>4)? atoi(argv[4]) : 2;    // numRecurInj
    int   invas    = (argc>5)? atoi(argv[5]) : 2;    // numInfectRepeat
    int   seed     = (argc>6)? atoi(argv[6]) : 1;
    int   nSteps   = (argc>7)? atoi(argv[7]) : 10000;
    int   doInt    = (argc>8)? atoi(argv[8]) : 0;    // testIntervention 0/1
    int   nStages  = (argc>9)? atoi(argv[9]) : 0;    // numIndependentInterventions
    int   intStart = (argc>10)? atoi(argv[10]) : 505;
    int   intSep   = (argc>11)? atoi(argv[11]) : 100;
    int   intLen   = (argc>12)? atoi(argv[12]) : 99;
    string intPath = (argc>13)? string(argv[13]) : string("");
    string outPath = (argc>14)? string(argv[14]) : string("iirabm_out.dat");

    if(nSteps > numWriteSteps) nSteps = numWriteSteps;

    numTimeSteps = nSteps;
    GA = 0;
    numReseeds = 0;
    reseedPostInt = 0;
    testIntervention = doInt;
    numIndependentInterventions = nStages;
    interventionTimeStart = intStart;
    interventionTimeSeperation = intSep;
    interventionLength = intLen;
    numTests = 0;
    intFile = intPath;

    for(int i=0;i<genomeLength;i++){ myIntervention[i]=1; testInt[i]=1; }
    if(doInt==1 && !intPath.empty()){
        readIntervention();                 // fills myIntervention (GA==0 path)
        for(int i=0;i<genomeLength;i++){ testInt[i]=myIntervention[i]; }
    }

    mainSimulation(inj, oxyHeal, toxi, envTox, invas, seed);

    ofstream f(outPath.c_str());
    f<<inj<<" "<<oxyHeal<<" "<<toxi<<" "<<envTox<<" "<<invas<<" "<<seed<<" "<<nSteps<<"\n";
    writeSeriesI(f, d_pmn, nSteps);
    writeSeriesI(f, d_mono, nSteps);
    writeSeriesI(f, d_TH1, nSteps);
    writeSeriesI(f, d_TH2, nSteps);
    writeSeries(f, d_oxy, nSteps);
    writeSeries(f, d_infect, nSteps);
    writeSeries(f, d_TNF, nSteps);
    writeSeries(f, d_sTNFr, nSteps);
    writeSeries(f, d_IL10, nSteps);
    writeSeries(f, d_GCSF, nSteps);
    writeSeries(f, d_pTH1, nSteps);
    writeSeries(f, d_pTH2, nSteps);
    writeSeries(f, d_PAF, nSteps);
    writeSeries(f, d_IFNg, nSteps);
    writeSeries(f, d_IL4, nSteps);
    writeSeries(f, d_IL8, nSteps);
    writeSeries(f, d_IL12, nSteps);
    writeSeries(f, d_sIL1r, nSteps);
    writeSeries(f, d_IL1ra, nSteps);
    writeSeries(f, d_IL1, nSteps);
    f.close();
    return 0;
}
'''
s = s + serial_driver
write("IIR_Sepsis.cpp", s)

print("patched ->", out)
