#ifndef ADJ_BOX_POWERSPEC_PROTO_H
#define ADJ_BOX_POWERSPEC_PROTO_H

int adj_box_powerspec_find_nearest_evaluate(int target, int mode, int *nexport, int *nsend_local);
double adj_box_powerspec_obtain_fields(void);
void adj_box_powerspec_save(const char *fname);
void adj_box_powerspec_gather_results(void);
void adj_box_powerspec(void);
int adj_box_powerspec_treefind(MyDouble searchcenter[3], MyFloat hsml, int target, int *startnode, int mode, int *nexport, int *nsend_local);
void adj_box_powerspec_calc_dispersion(void);

#endif

