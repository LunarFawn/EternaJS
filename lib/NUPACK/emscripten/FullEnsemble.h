#ifndef NUPACK_FULLENSEMBLE_H
#define NUPACK_FULLENSEMBLE_H

#include <string>
#include <vector>


struct FullAdvancedResult {
    double ensembleDefect;
    double ensembleDefectNormalized;
    double mfeDefect;
    double mfeDefectNormalized;
    std::vector<std::string> suboptStructures;
    std::vector<double> suboptEnergyError;
    std::vector<double> suboptFreeEnergy;
};

FullAdvancedResult* FullEnsembleNoBindingSite (const std::string& seqString, int temperature, float kcalDeltaRange, bool const pseudoknotted);
FullAdvancedResult* FullEnsembleWithOligos (const std::string& seqString, int temperature, float kcalDeltaRange,  bool const pseudoknotted);
std::string GetDotParensStructureFromFoldStructure(char* seq, const int *thepairs, bool makeForEterna);
void GetEnsembleDefectFromDotParens(char* seqString,char* dotParensStructure,
                                             int temperature, bool pseudoknot, bool multiFold, double *returnEnsembleDefect,  bool mfeDefect);
#endif //NUPACK_FULLENSEMBLE_H
