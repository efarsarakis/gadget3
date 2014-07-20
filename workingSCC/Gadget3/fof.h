#ifndef FOF_H
#define FOF_H

#include "allvars.h"

void fof_exchange_group_data(void);
int fof_compare_FOF_PList_MinID(const void *a, const void *b);
int fof_compare_FOF_GList_MinID(const void *a, const void *b);
int fof_compare_FOF_GList_MinIDTask(const void *a, const void *b);
int fof_compare_FOF_GList_MinIDTask_MinID(const void *a, const void *b);
int fof_compare_FOF_GList_LocCountTaskDiffMinID(const void *a, const void *b);
int fof_compare_FOF_GList_ExtCountMinID(const void *a, const void *b);
int fof_compare_Group_GrNr(const void *a, const void *b);
int fof_compare_Group_MinIDTask(const void *a, const void *b);
int fof_compare_Group_MinID(const void *a, const void *b);
int fof_compare_ID_list_GrNrID(const void *a, const void *b);
int fof_compare_Group_MinIDTask_MinID(const void *a, const void *b);
int fof_compare_Group_Len(const void *a, const void *b);

void fof_finish_group_properties(void);

int compare_group_mass_ID(const void *a, const void *b);
void fof_assign_HostHaloMass(void);

extern int Ngroups, TotNgroups;
extern long long TotNids;

typedef struct
{
  MyIDType MinID;
  MyIDType MinIDTask;
  int LocCount;
  int ExtCount;
#ifdef DENSITY_SPLIT_BY_TYPE
  int LocDMCount;
  int ExtDMCount;
#endif
  int GrNr;
} fof_group_list; 

typedef struct
{
  int Len;
  unsigned int Offset;
  MyIDType MinID;
  MyIDType MinIDTask;
  int GrNr;
  int LenType[6];
  MyOutputFloat MassType[6];
  MyOutputFloat Mass;
  MyOutputFloat CM[3];
  MyOutputFloat Vel[3];
  MyDouble FirstPos[3];
#ifdef SFR
  double Sfr;
#endif
#ifdef BLACK_HOLES
  MyOutputFloat BH_Mass;
  MyOutputFloat BH_Mdot;
  MyOutputFloat MaxDens;
  int index_maxdens, task_maxdens;
#endif

#ifdef SUBFIND
  int Nsubs;
  int FirstSub;
  MyDouble Pos[3];
  MyOutputFloat M_Mean200, R_Mean200;
  MyOutputFloat M_Crit200, R_Crit200;
  MyOutputFloat M_TopHat200, R_TopHat200;
  int ContaminationLen;
  MyOutputFloat ContaminationMass;
#ifdef SO_VEL_DISPERSIONS  
  MyOutputFloat VelDisp_Mean200, VelDisp_Crit200, VelDisp_TopHat200;
#endif
#ifdef SO_BAR_INFO
  MyOutputFloat gas_mass[3], star_mass[3], temp[3], xlum[3]; 
#endif
#endif

} group_properties;

typedef struct
{
  MyIDType ID;
  unsigned int GrNr;
} fof_id_list;

extern group_properties *Group;

#ifdef ALTERNATIVE_PSORT
void fof_sort_FOF_GList_LocCountTaskDiffMinID (fof_group_list *data, int ndata);

void fof_sort_FOF_GList_ExtCountMinID (fof_group_list *data, int ndata);

void fof_sort_Group_GrNr (group_properties *data, int ndata);

void fof_sort_ID_list_GrNrID (fof_id_list *data, int ndata);
#endif



#endif
