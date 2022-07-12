#include "FullFold.h"
#include "FullEnsemble.h"
#include "EmscriptenUtils.h"

#include "src/thermo/utils/pfuncUtilsConstants.h"
#include "src/thermo/utils/pfuncUtilsHeader.h"
#include "src/shared/utilsHeader.h"
#include "src/thermo/utils/DNAExternals.h"
#include <vector>
#include <utility>

FullFoldResult* FullFoldDefault (const std::string& seqString, bool const pseudoknotted = false) {
    auto autoSeqString = MakeCString(seqString);
    char* string = autoSeqString.get();

    int seqNum[MAXSEQLENGTH+1];
    int tmpLength = strlen(string);
    dnaStructures mfeStructs = {NULL, 0, 0, 0, 0};

    convertSeq(string, seqNum, tmpLength);

    if ( pseudoknotted ) {
        mfeFullWithSym(seqNum, tmpLength, &mfeStructs, 5, RNA,
                        1 /*DANGLETYPE*/, 37, TRUE, 1, SODIUM_CONC, MAGNESIUM_CONC,
                        USE_LONG_HELIX_FOR_SALT_CORRECTION);
    } else {
        mfeFullWithSym(seqNum, tmpLength, &mfeStructs, 3, RNA,
                        1 /*DANGLETYPE*/, 37, TRUE, 1, SODIUM_CONC, MAGNESIUM_CONC,
                        USE_LONG_HELIX_FOR_SALT_CORRECTION);
    }

    
    //initialize the result to return
    FullFoldResult* result = new FullFoldResult();
   
    
    //get dot bracket notation from data
    
    oneDnaStruct currentStruct = mfeStructs.validStructs[0];
    std::string singlestructure = GenerateDotBracketPairsList(string, currentStruct.theStruct);

    double energyError = currentStruct.error;
    double correctedEnergy = currentStruct.correctedEnergy;  

    result->structure = singlestructure;
    result->mfe = mfeStructs.validStructs[0].correctedEnergy;

    clearDnaStructures(&mfeStructs);


    return result;
}

FullFoldResult* FullFoldTemperature(double temperature_in, const std::string& seqString, bool const pseudoknotted = false) {
    auto autoSeqString = MakeCString(seqString);
    char* string = autoSeqString.get();

    int seqNum[MAXSEQLENGTH+1];
    int tmpLength = strlen(string);
    dnaStructures mfeStructs = {NULL, 0, 0, 0, 0};

    convertSeq(string, seqNum, tmpLength);

    // Note: perversely, though most conditionals suggest that >=5 is pseudoknotted,
    // in fact 6 is totally broken. The magic number is five. I take no pleasure in
    // reporting this.
    if ( pseudoknotted ) {
        mfeFullWithSym(seqNum, tmpLength, &mfeStructs, 5, RNA,
                        1 /*DANGLETYPE*/, temperature_in, TRUE, 1, SODIUM_CONC, MAGNESIUM_CONC,
                        USE_LONG_HELIX_FOR_SALT_CORRECTION);
    } else {
        mfeFullWithSym(seqNum, tmpLength, &mfeStructs, 3, RNA,
                        1 /*DANGLETYPE*/, temperature_in, TRUE, 1, SODIUM_CONC, MAGNESIUM_CONC,
                        USE_LONG_HELIX_FOR_SALT_CORRECTION);
    }

    FullFoldResult* result = new FullFoldResult();
    
    //get dot bracket notation from data
    
    oneDnaStruct currentStruct = mfeStructs.validStructs[0];
    std::string singlestructure = GenerateDotBracketPairsList(string, currentStruct.theStruct);

    double energyError = currentStruct.error;
    double correctedEnergy = currentStruct.correctedEnergy;  

    result->structure = singlestructure;
    result->mfe = mfeStructs.validStructs[0].correctedEnergy;

    clearDnaStructures(&mfeStructs);

    return result;
}

DotPlotResult* GetDotPlot (double temperature_in, const std::string& seqString) {
    auto autoSeqString = MakeCString(seqString);
    char* string = autoSeqString.get();

    double energy = 0;

    int seqNum[MAXSEQLENGTH+1];
    int tmpLength;

    DotPlotResult* result = new DotPlotResult();

    tmpLength = strlen(string);
    convertSeq(string, seqNum, tmpLength);

    pairPr = (DBL_TYPE*) calloc( (tmpLength+1)*(tmpLength+1), sizeof(DBL_TYPE));

    energy = pfuncFull( seqNum, 3, RNA, 1, temperature_in, 1, 1.0, 0.0, 0);

    for (int i = 0; i < tmpLength; i++) {
        for (int j = i+1; j < tmpLength; j++) {
            int k = (tmpLength+1)*i + j;
            if (pairPr[k] < 1e-5) continue;

            result->plot.push_back(i + 1);
            result->plot.push_back(j + 1);
            result->plot.push_back(pairPr[k]);
        }
    }

    if (pairPr) {
        free(pairPr);
        pairPr = NULL;
    }

    result->energy = energy;
    return result;
}

// binding site data
int g_site_i, g_site_j, g_site_p, g_site_q, g_site_bonus;
int _binding_cb(int i, int j, int* d, int* e) {
    if ((i == g_site_i) && (j == g_site_j)) {
        // fprintf(stderr, "query i:%d j:%d\n", i, j);
        (*d) = g_site_p;
        (*e) = g_site_q;
        return 1;
    }
    return 0;
}

DBL_TYPE _binding_site_cb(int i, int j, int p, int q) {
    // fprintf(stderr, "cb i:%d j:%d p:%d q:%d\n", i, j, p, q);
    if ((i == g_site_i) && (j == g_site_j) && (p == g_site_p) && (q == g_site_q)) {
        // fprintf(stderr, "cb match\n");
        return g_site_bonus * -.01;
    }

    return 0.;
}

extern int (*binding_cb)(int i, int j, int* d, int* e);
extern DBL_TYPE (*binding_site_cb)(int i, int j, int p, int q);

FullFoldResult* FullFoldWithBindingSite (const std::string& seqString, int site_i, int site_p, int site_j, int site_q, int site_bonus) {
    auto autoSeqString = MakeCString(seqString);
    char* string = autoSeqString.get();

    int seqNum[MAXSEQLENGTH+1];
    int tmpLength;
    dnaStructures mfeStructs = {NULL, 0, 0, 0, 0};
    int j;

    tmpLength = strlen(string);
    convertSeq(string, seqNum, tmpLength);

    // activate binding site callbacks
    binding_cb = _binding_cb;
    binding_site_cb = _binding_site_cb;
    g_site_i = site_i;
    g_site_p = site_p;
    g_site_j = site_j;
    g_site_q = site_q;
    g_site_bonus = site_bonus;
    mfeFullWithSym(seqNum, tmpLength, &mfeStructs, 3, RNA,
                   1 /*DANGLETYPE*/, 37, TRUE, 1, SODIUM_CONC, MAGNESIUM_CONC,
                   USE_LONG_HELIX_FOR_SALT_CORRECTION);
    // clean up
    binding_site_cb = NULL;
    binding_cb = NULL;

    FullFoldResult* result = new FullFoldResult();

    
    //get dot bracket notation from data
    
    oneDnaStruct currentStruct = mfeStructs.validStructs[0];
    std::string singlestructure = GenerateDotBracketPairsList(string, currentStruct.theStruct);

    double energyError = currentStruct.error;
    double correctedEnergy = currentStruct.correctedEnergy;  

    result->structure = singlestructure;
    result->mfe = mfeStructs.validStructs[0].correctedEnergy;

    clearDnaStructures(&mfeStructs);

    return result;
}

FullFoldResult* CoFoldSequence (const std::string& seqString) {
    auto autoSeqString = MakeCString(seqString);
    char* string = autoSeqString.get();

    int seqNum[MAXSEQLENGTH+1];
    dnaStructures mfeStructs = {NULL, 0, 0, 0, 0};
    int i, j;
    char* pc;
    int oligonumber=0;
    do {
        pc = strchr(string, '&');
        if (pc) 
        {
            (*pc) = '+';
            oligonumber++;
        }
    } while(pc);

    int seqLength = strlen(string);
    convertSeq(string, seqNum, seqLength);

    mfeFullWithSym(seqNum, seqLength, &mfeStructs, 3, RNA,
                   1 /*DANGLETYPE*/, 37, TRUE, 1, SODIUM_CONC, MAGNESIUM_CONC,
                   USE_LONG_HELIX_FOR_SALT_CORRECTION);

    
    //get dot bracket notation from data
    FullFoldResult* result = new FullFoldResult();
    oneDnaStruct currentStruct = mfeStructs.validStructs[0];
    std::string singlestructure = GenerateDotBracketPairsList(string , currentStruct.theStruct);

    double energyError = currentStruct.error;
    double correctedEnergy = currentStruct.correctedEnergy;  

    result->structure = singlestructure;
    result->mfe = mfeStructs.validStructs[0].correctedEnergy;

    clearDnaStructures(&mfeStructs);


    return result;
}

FullFoldResult* CoFoldSequenceWithBindingSite (const std::string& seqString, int site_i, int site_p, int site_j, int site_q, int site_bonus) {
    auto autoSeqString = MakeCString(seqString);
    char* string = autoSeqString.get();

    int seqNum[MAXSEQLENGTH+1];
    dnaStructures mfeStructs = {NULL, 0, 0, 0, 0};
    int i, j;
    char* pc;
    int oligonumber=0;
    do {
        pc = strchr(string, '&');
        if (pc) 
        {
            (*pc) = '+';
            oligonumber++;
        }
    } while(pc);

    int seqLength = strlen(string);
    convertSeq(string, seqNum, seqLength);

    // activate binding site callbacks
    binding_cb = _binding_cb;
    binding_site_cb = _binding_site_cb;
    g_site_i = site_i;
    g_site_p = site_p;
    g_site_j = site_j;
    g_site_q = site_q;
    g_site_bonus = site_bonus;
    mfeFullWithSym(seqNum, seqLength, &mfeStructs, 3, RNA,
                   1 /*DANGLETYPE*/, 37, TRUE, 1, SODIUM_CONC, MAGNESIUM_CONC,
                   USE_LONG_HELIX_FOR_SALT_CORRECTION);
    // clean up
    binding_site_cb = NULL;
    binding_cb = NULL;

    //get dot bracket notation from data
    FullFoldResult* result = new FullFoldResult();
    oneDnaStruct currentStruct = mfeStructs.validStructs[0];
    std::string singlestructure = GenerateDotBracketPairsList(string, currentStruct.theStruct);

    double energyError = currentStruct.error;
    double correctedEnergy = currentStruct.correctedEnergy;  

    result->structure = singlestructure;
    result->mfe = mfeStructs.validStructs[0].correctedEnergy;

    clearDnaStructures(&mfeStructs);

    return result;
}


