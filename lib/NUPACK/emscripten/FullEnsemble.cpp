#include "FullEnsemble.h"
#include "EmscriptenUtils.h"

#include "src/thermo/utils/pfuncUtilsConstants.h"
#include "src/thermo/utils/pfuncUtilsHeader.h"
#include "src/shared/utilsHeader.h"
#include "src/thermo/utils/DNAExternals.h"
#include <vector>
#include <utility>
#include <string> 


FullAdvancedResult* FullEnsembleWithOligos(const std::string& seqString, int temperature, float kcalDeltaRange, bool const pseudoknotted = false)
{

    //this chuck of code sets up the variables used to represent adn configure the sequence for the C code to use
    auto autoSeqString = MakeCString(seqString);
    char* string = autoSeqString.get();
    int seqNum[MAXSEQLENGTH+1];
    


    //runtime_constants.h defines NAD_INFINITY
    //#define NAD_INFINITY 100000 //artificial value for positive infinity
    //this is how teh dnastructure call looks for subopt.c which is waht iran when you ron normal nupack compiled code
    //here are comments from pfuncutilheaders.h on hwat the attributes of dnastructures consist of
    
    //oneDnaStruct and dnaStructures are used for enumerating sequences
       // typedef struct {
       // int *theStruct; //describes what is paired to what
       // DBL_TYPE error; //accumulated error (from the mfe) for a structure
       // DBL_TYPE correctedEnergy; //actual energy of a structure
       // int slength;
       // //(accounting for symmetry).
       //  } oneDnaStruct;
       //
       // typedef struct {
       // oneDnaStruct *validStructs;
       //int nStructs; //# of structures stored
       //int nAlloc; //# of structures allocated
       // int seqlength;
       // DBL_TYPE minError; //minimum deviation from mfe for all seqs
       //in validStructs
       //
       //} dnaStructures;
    
    //this struct will store
    //all the structures within the given range
    dnaStructures suboptStructs = {NULL, 0, 0, 0, NAD_INFINITY}; 
 
    //convert from how it comes from eterna to how nuapck needs it for joined oligos '+'
   
    char* pc;
    do {
        pc = strchr(string, '&');
        if (pc)(*pc) = '+';
    } while(pc);
  

    int seqStringLength = strlen(string);    
    convertSeq(string, seqNum, seqStringLength);
  

    //first get the ensemble through subopt
    if ( pseudoknotted ) {
        mfeFullWithSym_SubOpt(seqNum, seqStringLength, &suboptStructs, 5, RNA, 1 /*DANGLETYPE*/, 
                                temperature, TRUE, (DBL_TYPE) kcalDeltaRange, 
                                0, SODIUM_CONC, MAGNESIUM_CONC, USE_LONG_HELIX_FOR_SALT_CORRECTION);
    } else {
        mfeFullWithSym_SubOpt(seqNum, seqStringLength, &suboptStructs, 3, RNA, 1 /*DANGLETYPE*/, 
                                temperature, TRUE, (DBL_TYPE) kcalDeltaRange, 
                                0, SODIUM_CONC, MAGNESIUM_CONC, USE_LONG_HELIX_FOR_SALT_CORRECTION);
    }

    //initialize the result to return
    FullAdvancedResult* result = new FullAdvancedResult();
    
    int i, j;
    //get dot bracket notation from data
    



   for (i = 0; i < suboptStructs.nStructs; i++ ) {
        oneDnaStruct currentStruct = suboptStructs.validStructs[i];
        std::string singlestructure = GenerateDotBracketPairsList(string, currentStruct.theStruct);

        double energyError = currentStruct.error;
        double correctedEnergy = currentStruct.correctedEnergy;  

        result->suboptStructures.push_back(singlestructure);
        result->suboptEnergyError.push_back(energyError);
        result->suboptFreeEnergy.push_back(correctedEnergy);
    }

    clearDnaStructures(&suboptStructs);

    //get defect here later for now popuolate with 0
    result->ensembleDefect=0;


    return result;
}


FullAdvancedResult* FullEnsembleNoBindingSite(const std::string& seqString, int temperature, float kcalDeltaRange, bool const pseudoknotted = false)
{

    //this chuck of code sets up the variables used to represent adn configure the sequence for the C code to use
    auto autoSeqString = MakeCString(seqString);
    char* string = autoSeqString.get();
    int seqNum[MAXSEQLENGTH+1];
    int seqStringLength = strlen(string);
    //runtime_constants.h defines NAD_INFINITY
    //#define NAD_INFINITY 100000 //artificial value for positive infinity
    //this is how teh dnastructure call looks for subopt.c which is waht iran when you ron normal nupack compiled code
    //here are comments from pfuncutilheaders.h on hwat the attributes of dnastructures consist of
    
    //oneDnaStruct and dnaStructures are used for enumerating sequences
       // typedef struct {
       // int *theStruct; //describes what is paired to what
       // DBL_TYPE error; //accumulated error (from the mfe) for a structure
       // DBL_TYPE correctedEnergy; //actual energy of a structure
       // int slength;
       // //(accounting for symmetry).
       //  } oneDnaStruct;
       //
       // typedef struct {
       // oneDnaStruct *validStructs;
       //int nStructs; //# of structures stored
       //int nAlloc; //# of structures allocated
       // int seqlength;
       // DBL_TYPE minError; //minimum deviation from mfe for all seqs
       //in validStructs
       //
       //} dnaStructures;
    
    //this struct will store
    //all the structures within the given range
    dnaStructures suboptStructs = {NULL, 0, 0, 0, NAD_INFINITY};     
    convertSeq(string, seqNum, seqStringLength);



    //first get the ensemble through subopt
    if ( pseudoknotted ) {
        mfeFullWithSym_SubOpt(seqNum, seqStringLength, &suboptStructs, 5, RNA, 1 /*DANGLETYPE*/, 
                                temperature, TRUE, (DBL_TYPE) kcalDeltaRange, 
                                0, SODIUM_CONC, MAGNESIUM_CONC, USE_LONG_HELIX_FOR_SALT_CORRECTION);
    } else {
        mfeFullWithSym_SubOpt(seqNum, seqStringLength, &suboptStructs, 3, RNA, 1 /*DANGLETYPE*/, 
                                temperature, TRUE, (DBL_TYPE) kcalDeltaRange, 
                                0, SODIUM_CONC, MAGNESIUM_CONC, USE_LONG_HELIX_FOR_SALT_CORRECTION);
    }

    //initialize the result to return
    FullAdvancedResult* result = new FullAdvancedResult();
    
    int i, j;
    //get dot bracket notation from data
    



    for (i = 0; i < suboptStructs.nStructs; i++ ) {
        oneDnaStruct currentStruct = suboptStructs.validStructs[i];
        std::string singlestructure = GenerateDotBracketPairsList(string, currentStruct.theStruct);

        double energyError = currentStruct.error;
        double correctedEnergy = currentStruct.correctedEnergy;  

        result->suboptStructures.push_back(singlestructure);
        result->suboptEnergyError.push_back(energyError);
        result->suboptFreeEnergy.push_back(correctedEnergy);
    }

    clearDnaStructures(&suboptStructs);

    //get defect here later for now popuolate with 0
    result->ensembleDefect=0;


    return result;
}

std::string GenerateDotBracketPairsList(char* seq, const int *thepairs) {
    /*
    This prints the structure of the fold using a '.' for 
    unpaired bases, and ( ), { }, [ ], < > for pairs. 

    The file specified by file name is appended.

    Initially, thefold[i] = '.' for all i

    If this ever becomes the slow step, it can be optimized to run faster
    */
   //count oligos
    char* pc;
    do {
        pc = strchr(seq, '&');
        if (pc)(*pc) = '+';
    } while(pc);

    
  
    int seqNum[ MAXSEQLENGTH+1]; 
    int isNicked[ MAXSEQLENGTH];
    int nNicks = 0;

    int nicks[MAXSTRANDS];
    int nickIndex;
    int **etaN;
    int length, tmpLength;
    int seqlength;

   
    //the rest is for printing purposes
    tmpLength = length = strlen( seq);
    int i,j, pos;

    for( i = 0; i < tmpLength; i++) {
        isNicked[i] = 0;
        if( seq[i] == '+') {        
        length--;
        isNicked[ i - nNicks++ -1] = 1;
        } 
    }

    

    //initialize nicks
    for( i = 0; i < MAXSTRANDS; i++) {
        nicks[i] = -1;
    }

    nickIndex = 0;
    for( i = 0; i < length; i++) {
        if( isNicked[i])
        nicks[ nickIndex++] = i;
    }

    seqlength = length;
    //overkill, but convenient
    etaN = (int**) malloc( (length*(length+1)/2 + (length+1))*sizeof( int*));
    InitEtaN( etaN, nicks, length);

    int nStrands = etaN[ EtaNIndex( 0.5, seqlength-0.5, seqlength)][0]+1;
    char *thefold = (char*) malloc( (seqlength + nStrands) * sizeof(char));
    
    char pairSymbols[] = { '(', ')', '{','}', '[', ']', '<', '>' };
    int type = 0;
    int nTypes = 4;
   
    int **pairlist; // Each row is i,j pair
    int npairs; // number of pairs in structure
   

    char *parensString;
    parensString = ( char*) malloc( (seqlength+1)*sizeof( char) );
    int lastL, lastR;

    
    // Allocate memory for pairlist (this is more than we need, but be safe)
    pairlist = (int **) malloc(seqlength * sizeof(int *));
    for (i = 0; i < seqlength; i++) {
        pairlist[i] = (int *) malloc(2 * sizeof(int));
    }

    // Create pairlist from thepairs
    npairs = 0;
    for( j = 0; j < seqlength; j++) {
        if(thepairs[j] > j) {
        pairlist[npairs][0] = j;
        pairlist[npairs++][1] = thepairs[j];
        }
    }

    // Creat dot-paren structure
    for( i = 0; i < seqlength+1; i++) {
        parensString[i] = '.';
    }

    //offSet = 0;
    lastL = -1; 
    lastR = seqlength;
    for( i = 0; i < seqlength; i++) {
        if( thepairs[i] != -1 && thepairs[i] > i) {
        if( i > lastR || thepairs[i] > lastR) {
            for( j = 0; j < i; j++) {
            if( thepairs[j] > i && thepairs[j] < thepairs[i]) {
                type = (type + 1) % nTypes;
                break;
            }
            }
        }

        parensString[i] = pairSymbols[ 2*type];
        parensString[ thepairs[i]] = pairSymbols[2*type + 1];
        lastL = i;
        lastR = thepairs[i];
        }
    }

    for( i = 0; i < seqlength+nStrands-1; i++) {
        thefold[i] = '.';
    }

    pos = 0;
    for( i = 0; i < length; i++) {
        thefold[ pos++] = parensString[ i];
        if( etaN[ EtaNIndex_same(i+0.5, seqlength)][0] == 1) {
        thefold[ pos++] = '+';
        }
    }

    //now create a string from the characer array to make the string version of the structure
    
    
    std::string dotBracketStructure;
    for (int k = 0; k < length; k++ ) {
          
        dotBracketStructure.push_back(thefold[k]);
          
    }  

    return dotBracketStructure;
    
}

