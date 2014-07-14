
#include "allvars.h"

#ifdef SUBFIND
#include "subfind.h"
#include "fof.h"


int Ncollective;
int MaxNsubgroups;
int Nsubgroups;
int TotNsubgroups;

#ifdef SUBFIND_COUNT_BIG_HALOS
int Nbiggroups;
int TotNbiggroups;
#endif

#ifdef LT_ADD_GAL_TO_SUB
int bc_use_index[12] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
char bc_name[12][2] = { "u ", "V ", "g ", "r ", "i ", "z ", "Y ", "J ", "H ", "K ", "L ", "M " };
#endif

struct subgroup_properties *SubGroup;


struct nearest_r2_data *R2Loc;

struct nearest_ngb_data *NgbLoc;


struct r2data *R2list;

//struct nearest_ngb_data *NgbLoc;

//struct nearest_r2_data *R2Loc;

double *Dist2list;

#endif
