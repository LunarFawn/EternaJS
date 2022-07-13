#include <vector>
#include <utility>
#include <string> 


#include "FullEnsemble.h"
#include "EmscriptenUtils.h"


#include "src/thermo/utils/pfuncUtilsConstants.h"
#include "src/thermo/utils/pfuncUtilsHeader.h"
#include "src/shared/utilsHeader.h"
#include "src/thermo/utils/DNAExternals.h"

extern DBL_TYPE *pairPrPbg;  //for pseudoknots
extern DBL_TYPE *pairPrPb;  //for pseudoknots

extern double CUTOFF;
extern int Multistranded;

bool forEternaDigest = true;
bool notforEternaDigest= false;

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
        if (pc) (*pc) = '+';
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
    
   std::string mfeStructure;


   for (i = 0; i < suboptStructs.nStructs; i++ ) {
        oneDnaStruct currentStruct = suboptStructs.validStructs[i];
        //for eternadigest sets true 
        //notforeternadigest sets false
        
        std::string singleStructure = GetDotParensStructureFromFoldStructure(string, currentStruct.theStruct, forEternaDigest);
        if(i == 0)
        {
          //this is the first one so it is the mfe
          
          mfeStructure =  GetDotParensStructureFromFoldStructure(string, currentStruct.theStruct, notforEternaDigest);
        }

        double energyError = currentStruct.error;
        double correctedEnergy = currentStruct.correctedEnergy;  

        result->suboptStructures.push_back(singleStructure);
        result->suboptEnergyError.push_back(energyError);
        result->suboptFreeEnergy.push_back(correctedEnergy);
    }

    clearDnaStructures(&suboptStructs);
    
    //get defect here later for now popuolate with 0
    auto autoStructString = MakeCString(mfeStructure);
    char* newSecondStruct = autoStructString.get();

    //get ensemble defect for mfe structure
     double* ensembleDefectArray = new double[2];
     double* mfeEnsembleDefectArray = new double[2];
     double ensembleDefect = -2;
     double ensembleDefectNormalized = -3;
     double mfeDefect = -4;
     double mfeDefectNormalized = -5;
     bool getDefectMFE = false;
    GetEnsembleDefectFromDotParens(string, newSecondStruct, temperature, pseudoknotted, false, ensembleDefectArray, getDefectMFE);
    ensembleDefect = ensembleDefectArray[0];
    ensembleDefectNormalized = ensembleDefectArray[1];
    
    getDefectMFE=true;
    GetEnsembleDefectFromDotParens(string, newSecondStruct, temperature, pseudoknotted, false, mfeEnsembleDefectArray, getDefectMFE);
    mfeDefect = mfeEnsembleDefectArray[0];
    mfeDefectNormalized = mfeEnsembleDefectArray[1];

    result->ensembleDefect=ensembleDefect;
    result->ensembleDefectNormalized=ensembleDefectNormalized;
    result->mfeDefect=mfeDefect;
    result->mfeDefectNormalized=mfeDefectNormalized;


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
    
    std::string mfeStructure;


    for (i = 0; i < suboptStructs.nStructs; i++ ) {
        oneDnaStruct currentStruct = suboptStructs.validStructs[i];
        bool forEternaDigest = true;
        std::string singleStructure = GetDotParensStructureFromFoldStructure(string, currentStruct.theStruct, forEternaDigest);
        if(i == 0)
        {
          //this is the first one so it is the mfe
          mfeStructure = singleStructure;
        }

        double energyError = currentStruct.error;
        double correctedEnergy = currentStruct.correctedEnergy;  

        result->suboptStructures.push_back(singleStructure);
        result->suboptEnergyError.push_back(energyError);
        result->suboptFreeEnergy.push_back(correctedEnergy);
    }

    clearDnaStructures(&suboptStructs);

    //get defect here later for now popuolate with 0

    //get ensemble defect for mfe structure
   // double ensembleDefect;
    //double ensembleDefectNormalized;
   // double mfeDefect;
    //double mfeDefectNormalized;
    //bool getDefectMFE = false;
    //ensembleDefect, ensembleDefectNormalized = GetEnsembleDefectFromDotParens(seqString, mfeStructure, temperature, pseudoknotted, false, getDefectMFE);
    
    //getDefectMFE=true;
    //mfeDefect, mfeDefectNormalized = GetEnsembleDefectFromDotParens(seqString, mfeStructure, temperature, pseudoknotted, false, getDefectMFE);
    
    //result->ensembleDefect=ensembleDefect;
    //result->ensembleDefectNormalized=ensembleDefectNormalized;
    //result->mfeDefect=mfeDefect;
    //result->mfeDefectNormalized=mfeDefectNormalized;

    result->ensembleDefect=0;
    result->ensembleDefectNormalized=1;
    result->mfeDefect=2;
    result->mfeDefectNormalized=3;


    return result;
}


std::string GetDotParensStructureFromFoldStructure(char* RNAseq, const int *thepairs, bool makeForEterna) {
    /*
    This prints the structure of the fold using a '.' for 
    unpaired bases, and ( ), { }, [ ], < > for pairs. 

    The file specified by file name is appended.

    Initially, thefold[i] = '.' for all i

    If this ever becomes the slow step, it can be optimized to run faster
    */
   //count oligos
    char* newSeq;
    do {
        newSeq = strchr(RNAseq, '&');
        if (newSeq) (*newSeq) = '+';
    } while(newSeq);
  
    int seqNum[ MAXSEQLENGTH+1]; 
    int isNicked[ MAXSEQLENGTH];
    int nNicks = 0;

    int nicks[MAXSTRANDS];
    int nickIndex;
    int **etaN;
    int length, tmpLength, seqlength;
    

   
    //the rest is for printing purposes
    tmpLength = length = strlen(RNAseq);
    int i,j, pos;

    for( i = 0; i < tmpLength; i++) {
        isNicked[i] = 0;
        if( RNAseq[i] == '+') {        
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
    for( i = 0; i < seqlength; i++) {
        thefold[ pos++] = parensString[ i];
        if( etaN[ EtaNIndex_same(i+0.5, seqlength)][0] == 1) {
            if (makeForEterna==true)
            {
                thefold[ pos++] = '&';
            }
            else
            {
                thefold[ pos++] = '+';
            }   
        }
    }

    //now create a string from the characer array to make the string version of the structure
    
    
    std::string dotBracketStructure = "";
    for (int k = 0; k < seqlength; k++ ) {
          
        dotBracketStructure.push_back(thefold[k]);
          
    }  

    return dotBracketStructure;
    
}

void GetEnsembleDefectFromDotParens(char* seqChar, char* dotParensStructure, int temperature, bool pseudoknot, bool multiFold, double *returnEnsembleDefect, bool mfeDefect)
{   

    
    
  DANGLETYPE=1;  

  if (multiFold==true)
  {
    Multistranded = 1;
  }
  else
  {
    Multistranded = 0;
  }

  if (pseudoknot == true)
  {
    DO_PSEUDOKNOTS = 1;
  }
  else
  {
    DO_PSEUDOKNOTS = 0;
  }

  if (mfeDefect==true)
  {
    USE_MFE = 1;
    // -degenerate flag not used for defect calcs, force ONLY_ONE_MFE
    ONLY_ONE_MFE = 1;
  }
  else
  {
    USE_MFE = 0;
    ONLY_ONE_MFE = 0;
  }

  int i, j;
  int trySymetry=TRUE;
  
  
  DBL_TYPE nsStar;
  DBL_TYPE mfe; // Minimum free energy (not really used)
  DBL_TYPE ene; // Free energy (used to check if the structure is legal)
  int vs;
  int complexity;
  int tmpLength;
  int seqlength;
  int nbases; // Number of bases as read from ppairs or mfe file
 
  int nNicks;  // Number of strands
  int doCalc; // Whether we need to compute the pair probability matrix/mfe or not
  
   int seqStringLength;
  char *tok; // Token
  char tokseps[] = " \t\n"; // Token separators

  //convert structure and sequence to char

 
  
  int seqNum[MAXSEQLENGTH+1];
  seqlength = seqStringLength = tmpLength = strlen(seqChar);
  

  convertSeq(seqChar, seqNum, seqStringLength);

 
  
  // Get the number of strand breaks
  nNicks = 0;
  for (i = 0; i < tmpLength; i++) {
    if (seqChar[i] == '+') {
      nNicks++;
      seqlength--;
    }
  }

  int thepairs[MAXSEQLENGTH+1];
  getStructureFromParens( dotParensStructure, thepairs, seqlength);
  
  ene = naEnergyPairsOrParensFullWithSym( thepairs, NULL, seqNum, RNA, 1 /*DANGLETYPE*/,
          temperature, trySymetry,
          SODIUM_CONC, MAGNESIUM_CONC,
          USE_LONG_HELIX_FOR_SALT_CORRECTION);

  // Allocate memory for storing pair probabilities
   pairPr = (DBL_TYPE*) calloc( (seqlength+1)*(seqlength+1), sizeof(DBL_TYPE));
  
  //stuff for folding and getting probs
    dnaStructures mfeStructs = {NULL, 0, 0, 0, NAD_INFINITY};
  
 //printInputs( argc, argv, seqChar, vs, NULL, parens, "screen");

    if (USE_MFE) {
      if( !DO_PSEUDOKNOTS ) {
        complexity = 3;
      }
      else {
        complexity = 5;
      }

      // Compute MFE and MFE structure
      mfe = mfeFullWithSym( seqNum, tmpLength, &mfeStructs, complexity, RNA, DANGLETYPE, temperature, trySymetry,
          ONLY_ONE_MFE, SODIUM_CONC, MAGNESIUM_CONC,
          USE_LONG_HELIX_FOR_SALT_CORRECTION);

      // Compute nsStar from output
      nsStar = 0.0;
      for (i = 0; i < seqlength; i++) {
        if (thepairs[i] != (mfeStructs.validStructs)[0].theStruct[i]) {
          nsStar += 1.0;
        }
      }
    }
    else
    {
      // Allocate memory for storing pair probabilities
      pairPrPbg = (DBL_TYPE*) calloc( (seqlength+1)*(seqlength+1), sizeof(DBL_TYPE));
      pairPrPb = (DBL_TYPE*) calloc( (seqlength+1)*(seqlength+1), sizeof(DBL_TYPE));

      if( !DO_PSEUDOKNOTS ) {
        complexity = 3;
      }
      else {
        complexity = 5;
      }

      nsStar = nsStarPairsOrParensFull(seqlength, seqNum, thepairs, NULL,
               complexity, RNA, 1 /*DANGLETYPE*/,
               temperature, SODIUM_CONC,
               MAGNESIUM_CONC, USE_LONG_HELIX_FOR_SALT_CORRECTION);


      free( pairPrPbg);
      free( pairPrPb);
    }

  

  double EnsembleDefect;
  double EnsembleDefectNormalized;
  
  
  if (USE_MFE) {
    //Fraction of correct nucleotides vs. MFE
    EnsembleDefect = (long double) nsStar;
    EnsembleDefectNormalized = (long double) nsStar/seqlength;
  }
  else {
    //Ensemble defect n(s,phi) and normalized ensemble defect n(s,phi)/N;
    EnsembleDefect = (long double) nsStar;
    EnsembleDefectNormalized = (long double) nsStar/seqlength;
  }
  returnEnsembleDefect[0]=EnsembleDefect;
  returnEnsembleDefect[1]=EnsembleDefectNormalized;
  
  

}
