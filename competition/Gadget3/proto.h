#ifndef ALLVARS_H
#include "allvars.h"
#endif
#include "forcetree.h"
#include "domain.h"
#include "cooling.h"
#ifdef COSMIC_RAYS
#include "cosmic_rays.h"
#endif
#ifdef LT_STELLAREVOLUTION
#include "lt_sfr/lt.h"
#endif
#ifdef CS_MODEL
#include "cs_model/cs_metals.h"
#endif
#ifdef CHEMCOOL
#include "f2c.h"
#endif
#ifdef VORONOI
#include "voronoi.h"
#endif

#ifdef ADJ_BOX_POWERSPEC
#include "power_spec/adj_box_powerspec_proto.h"
#endif


#ifdef HAVE_HDF5
#include <hdf5.h>
void write_header_attributes_in_hdf5(hid_t handle);
void read_header_attributes_in_hdf5(char *fname);
void write_parameters_attributes_in_hdf5(hid_t handle);
void write_units_attributes_in_hdf5(hid_t handle);
void write_constants_attributes_in_hdf5(hid_t handle);
#endif
void report_VmRSS(void);
void output_compile_time_options(void);

void report_pinning(void);
void pin_to_core_set(void);
void get_core_set(void);

void init_turb(void);
void set_turb_ampl(void);
void add_turb_accel(void);
void reset_turb_temp(void);
void log_turb_temp(void);
void init_static_nfw();

void create_snapshot_if_desired(void);

void mpi_report_comittable_memory(long long BaseMem);
long long report_comittable_memory(long long *MemTotal,
                                   long long *Committed_AS,
                                   long long *SwapTotal,
                                   long long *SwapFree);

void do_first_halfstep_kick(void);
void do_second_halfstep_kick(void);
void find_timesteps(void);
void compute_hydro_accelerations(void);
void compute_densities(void);
void compute_grav_accelerations(void);
void calc_memory_checksum(void *base, size_t bytes);

void get_disk_forces(double RR, double zz, double *f_R, double *f_z);
double get_disk_mass(double time);
void growing_disk_init(void);

double get_turb_pot(double x, double y, double z);

void   sub_turb_move_perturbers(double t0, double t1);
void   sub_turb_add_forces(void);
void   sub_turb_read_table(void);
void   sub_turb_parent_halo_accel(double dx, double dy, double dz, double *acc);
double sub_turb_enclosed_mass(double r, double msub, double vmax, double radvmax, double c);


int powerspec_turb_treefind(MyDouble searchcenter[3], MyFloat hsml, int target, int *startnode, int mode, int *nexport, int *nsend_local);
int powerspec_turb_find_nearest_evaluate(int target, int mode, int *nexport, int *nsend_local);
void powerspec_turb_calc_dispersion(void);
double powerspec_turb_obtain_fields(void);
void powerspec_turb_save(char *fname, double *disp);
void powerspec_turb_collect(void);
void powerspec_turb(int filenr);


void set_cosmo_factors_for_current_time(void);
void drift_sph_extra_physics(int i, integertime tstart, integertime tend, double dt_entr);

void set_non_standard_physics_for_current_time(void);
void calculate_non_standard_physics(void);
void compute_statistics(void);
void execute_resubmit_command(void);
void make_list_of_active_particles(void);
void output_extra_log_messages(void);

static inline double DMAX(double a, double b) { return (a > b) ? a : b; } 
static inline double DMIN(double a, double b) { return (a < b) ? a : b; }
static inline int IMAX(int a, int b) { return (a > b) ? a : b; } 
static inline int IMIN(int a, int b) { return (a < b) ? a : b; } 


void do_distortion_tensor_kick(int i, double dt_gravkick);
void set_predicted_sph_quantities_for_extra_physics(int i);
void do_sph_kick_for_extra_physics(int i, integertime tstart, integertime tend, double dt_entr);

void check_particle_for_temperature_minimum(int i);

double get_pressure(int i);

void read_fof(int num);
int fof_compare_ID_list_ID(const void *a, const void *b);

void myfree_msg(void *p, char *msg);
void kspace_neutrinos_init(void);

void twopoint(void);
void twopoint_save(void);
int twopoint_ngb_treefind_variable(MyDouble searchcenter[3], MyFloat rsearch, int target, int *startnode,
				   int mode, int *nexport, int *nsend_local);
int twopoint_count_local(int target, int mode, int *nexport, int *nsend_local);

void powerspec(int flag, int *typeflag);
double PowerSpec_Efstathiou(double k);
void powerspec_save(void);
void foldonitself(int *typelist);
void dump_potential(void);

int snIaheating_evaluate(int target, int mode, int *nexport, int *nSend_local);
void snIa_heating(void);
void voronoi_setup_exchange(void);

double get_neutrino_powerspec(double k, double ascale);
double get_powerspec(double k, double ascale);
void init_transfer_functions(void);


int MPI_Check_Sendrecv(void *sendbuf, int sendcount, MPI_Datatype sendtype,
                       int dest, int sendtag, void *recvbufreal, int recvcount,
                       MPI_Datatype recvtype, int source, int recvtag, MPI_Comm comm, MPI_Status * status);

int MPI_Sizelimited_Sendrecv(void *sendbuf, int sendcount, MPI_Datatype sendtype,
			     int dest, int sendtag, void *recvbuf, int recvcount,
			     MPI_Datatype recvtype, int source, int recvtag, MPI_Comm comm, MPI_Status * status);

int mpi_calculate_offsets(int *send_count, int *send_offset, int *recv_count, int *recv_offset, int send_identical);
void sort_based_on_field(void *data, int field_offset, int n_items, int item_size, void **data2ptr);
void mpi_distribute_items_to_tasks(void *data, int task_offset, int *n_items, int *max_n, int item_size);

void parallel_sort_special_P_GrNr_ID(void);
void calculate_power_spectra(int num, long long *ntot_type_all);

int pmforce_is_particle_high_res(int type, MyDouble *pos);

void compare_partitions(void);
void assign_unique_ids(void);
int permut_data_compare(const void *a, const void *b);
void  generate_permutation_in_active_list(void);
void get_particle_numbers(char *fname, int num_files);

void conduction(void);
void conduction_matrix_multiply(double *in, double *out);
double conduction_vector_multiply(double *a, double *b);
int conduction_evaluate(int target, int mode, double *in, double *out, double *sum,
			int *nexport, int *nsend_local);


void fof_get_group_center(double *cm, int gr);
void fof_get_group_velocity(double *cmvel, int gr);
int fof_find_dmparticles_evaluate(int target, int mode, int *nexport, int *nsend_local);
void fof_compute_group_properties(int gr, int start, int len);

void parallel_sort(void *base, size_t nmemb, size_t size, int (*compar) (const void *, const void *));
void parallel_sort_comm(void *base, size_t nmemb, size_t size, int (*compar) (const void *, const void *), MPI_Comm comm);
int compare_IDs(const void *a, const void *b);
void test_id_uniqueness(void);

void check_particles_info(const char *func, const char *file, int linenr);

int io_compare_P_ID(const void *a, const void *b);
int io_compare_P_GrNr_SubNr(const void *a, const void *b);


void drift_particle(int i, integertime time1);
int ShouldWeDoDynamicUpdate(void);

void put_symbol(double t0, double t1, char c);
void write_cpu_log(void);

int get_timestep_bin(integertime ti_step);

const char* svn_version(void);

void find_particles_and_save_them(int num);
void lineofsight_output(void);
void sum_over_processors_and_normalize(void);
void absorb_along_lines_of_sight(void);
void output_lines_of_sight(int num);
integertime find_next_lineofsighttime(integertime time0);
integertime find_next_gridoutputtime(integertime ti_curr);
void add_along_lines_of_sight(void);
double los_periodic(double x);
void do_the_kick(int i, integertime tstart, integertime tend, integertime tcurrent);


void x86_fix(void) ;

void *mymalloc_fullinfo(const char *varname, size_t n, const char *func, const char *file, int linenr);
void *mymalloc_movable_fullinfo(void *ptr, const char *varname, size_t n, const char *func, const char *file, int line);

void *myrealloc_fullinfo(void *p, size_t n, const char *func, const char *file, int line);
void *myrealloc_movable_fullinfo(void *p, size_t n, const char *func, const char *file, int line);

void myfree_fullinfo(void *p, const char *func, const char *file, int line);
void myfree_movable_fullinfo(void *p, const char *func, const char *file, int line);

void mymalloc_init(void);
void dump_memory_table(void);
void report_detailed_memory_usage_of_largest_task(size_t *OldHighMarkBytes, const char *label, const char *func, const char *file, int line);

double get_shear_viscosity(int i);

void kinetic_feedback_mhm(void);
int kin_compare_key(const void *a, const void *b);
void kinetic_evaluate(int target, int mode);

void bubble(void);
void multi_bubbles(void);
void bh_bubble(double bh_dmass, MyFloat center[3], MyIDType BH_id);
void find_CM_of_biggest_group(void);
int compare_length_values(const void *a, const void *b);
double rho_dot(double z, void *params);
double bhgrowth(double z1, double z2);

void second_order_ics(void);
double F1_Omega(double);
double F2_Omega(double);

int smoothed_evaluate(int target, int mode, int *nexport, int *nsend_local);
void smoothed_values(void);

int fof_find_dmparticles_evaluate(int target, int mode, int *nexport, int *nsend_local);

double INLINE_FUNC hubble_function(double a);
#ifdef DARKENERGY
double DarkEnergy_a(double);
double DarkEnergy_t(double);
#ifdef TIMEDEPDE
void fwa_init(void);
double INLINE_FUNC fwa(double);
double INLINE_FUNC get_wa(double);
#ifdef TIMEDEPGRAV
double INLINE_FUNC dHfak(double a);
double INLINE_FUNC dGfak(double a);
#endif
#endif
#endif

#ifdef EXTERNALHUBBLE
double INLINE_FUNC hubble_function_external(double a);
#endif

#ifdef MODGRAV
double INLINE_FUNC hubble_function_modgrav(double a);
double INLINE_FUNC mg_R_background(double a);
double INLINE_FUNC mg_phi_background(double a);
double INLINE_FUNC mg_V(double phi);
#ifdef MODGRAV_FoR
double INLINE_FUNC mg_R(double phi);
#endif
double INLINE_FUNC mg_dV(double phi);
double INLINE_FUNC mg_ddV(double phi);
double INLINE_FUNC mg_beta0(double phi);
double INLINE_FUNC mg_beta1(double phi);
double INLINE_FUNC mg_dm0(double phi);
double INLINE_FUNC mg_dm1(double phi);

void mg_fieldsolve ();
void mg_save_timing_fieldsolve ();
void save_amr_data ();
void save_amr_data_to_file (char *filenamebase);
void save_multigrid_data (int minlevel, int maxlevel);
void map_multigrid_cells_to_new_tree (struct amr_cell_data *old_multigrid_cells, int num_old_multigrid_cells, int find_tasks);
#endif


void blackhole_accretion(void);
int blackhole_evaluate(int target, int mode, int *nexport, int *nsend_local);
int blackhole_evaluate_swallow(int target, int mode, int *nexport, int *nsend_local);

int  blackhole_compare_key(const void *a, const void *b);


int ngb_treefind_fof_nearest(MyDouble searchcenter[3], MyFloat hsml, int target, int *startnode, int mode,
			     int *nexport, int *nsend_local, int MyFOF_PRIMARY_LINK_TYPES);

void ngb_init(void);

void fof_fof(int num);
void fof_import_ghosts(void);
void fof_course_binning(void);
void fof_find_groups(void);
void fof_check_cell(int p, int i, int j, int k);
void fof_find_minids(void);
int fof_link_accross(void);
void fof_exchange_id_lists(void);
int fof_grid_compare(const void *a, const void *b);
void fof_compile_catalogue(void);
void fof_save_groups(int num);
void fof_save_local_catalogue(int num);
double fof_periodic(double x);
double fof_periodic_wrap(double x);
void fof_find_nearest_dmparticle(void);
int fof_find_nearest_dmparticle_evaluate(int target, int mode, int *nexport, int *nsend_local);

int fof_compare_key(const void *a, const void *b);
void fof_link_special(void);
void fof_link_specialpair(int p, int s);
void fof_make_black_holes(void);

int io_compare_P_GrNr_ID(const void *a, const void *b);

void write_file(char *fname, int readTask, int lastTask);

void distribute_file(int nfiles, int firstfile, int firsttask, int lasttask, int *filenr, int *master,
		     int *last);

int get_values_per_blockelement(enum iofields blocknr);

int get_datatype_in_block(enum iofields blocknr);
void get_dataset_name(enum iofields blocknr, char *buf);


int blockpresent(enum iofields blocknr);
void fill_write_buffer(enum iofields blocknr, int *pindex, int pc, int type);
void empty_read_buffer(enum iofields blocknr, int offset, int pc, int type);

int get_particles_in_block(enum iofields blocknr, int *typelist);

int get_bytes_per_blockelement(enum iofields blocknr, int mode);

void read_file(char *fname, int readTask, int lastTask);

void get_Tab_IO_Label(enum iofields blocknr, char *label);


void long_range_init_regionsize(void);

int find_files(char *fname);

int metals_compare_key(const void *a, const void *b);
void enrichment_evaluate(int target, int mode);

int hydro_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist);
void *hydro_evaluate_primary(void *p);
void *hydro_evaluate_secondary(void *p);

void pm_init_nonperiodic_allocate(void);

void  pm_init_nonperiodic_free(void);

double get_random_number(MyIDType id);
void set_random_numbers(void);

int grav_tree_compare_key(const void *a, const void *b);
int dens_compare_key(const void *a, const void *b);
int hydro_compare_key(const void *a, const void *b);

int data_index_compare(const void *a, const void *b);
int peano_compare_key(const void *a, const void *b);

#if defined(PM_FROM_L3) && defined(PMGRID)
void create_plans(void);
#endif

#ifdef SORT_FROM_L3
void quicksort(void *a, size_t n, size_t es, int (*cmp) (const void *, const void *));
void serial_sort(void *base, size_t nmemb, size_t size, int (*compar) (const void *, const void *));
#ifdef OMP_SORT
void omp_qsort(void *a, size_t n, size_t es, int (*cmp) (const void *, const void *));
void serial_sort_omp(void *base, size_t nmemb, size_t size, int (*compar) (const void *, const void *));
int floorLog2(unsigned int n);
#endif
#endif

void mysort_dataindex(void *b, size_t n, size_t s, int (*cmp) (const void *, const void *));
void mysort_domain(void *b, size_t n, size_t s);
void mysort_idlist(void *b, size_t n, size_t s, int (*cmp) (const void *, const void *));
void mysort_pmperiodic(void *b, size_t n, size_t s, int (*cmp) (const void *, const void *));
void mysort_pmnonperiodic(void *b, size_t n, size_t s, int (*cmp) (const void *, const void *));
void mysort_peano(void *b, size_t n, size_t s, int (*cmp) (const void *, const void *));

void check_wind_creation(void);
void treat_outflowing_particles(void);
void set_injection_accel(void);


int density_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist);
void *density_evaluate_primary(void *p);
void *density_evaluate_secondary(void *p);
int density_isactive(int n);

void GetMachNumberCR( struct sph_particle_data *Particle );
#ifndef CS_MODEL
void GetMachNumber( struct sph_particle_data* Particle );
#else
void GetMachNumber( struct sph_particle_data* Particle, struct particle_data* Particle_Hsml );
#endif
void GetShock_DtEnergy( struct sph_particle_data* Particle );

#ifdef MAGNETIC
#ifdef BFROMROTA
void rot_a(void);
void rot_a_evaluate(int i, int mode);
#endif
#endif
size_t sizemax(size_t a, size_t b);


void reconstruct_timebins(void);

void init_peano_map(void);
peanokey peano_hilbert_key(int x, int y, int z, int bits);
peanokey peano_and_morton_key(int x, int y, int z, int bits, peanokey *morton);
peanokey morton_key(int x, int y, int z, int bits);

void catch_abort(int sig);
void catch_fatal(int sig);
void terminate_processes(void);
void enable_core_dumps_and_fpu_exceptions(void);
void write_pid_file(void);

#ifdef PAUSE_RUN_TO_ATTACH_DEBUGGER
void pause_run_to_attach_debugger();
#endif

void pm_init_periodic_allocate(void);

void pm_init_periodic_free(void);

void move_particles(integertime time1);


void find_next_sync_point_and_drift(void);
void find_dt_displacement_constraint(double hfac);
#ifdef MAKEGLASS
void do_glass_making_step(void);
#endif
#ifdef RELAXOBJECT
void determine_relaxfac(void);
#endif
#ifdef WAKEUP
void process_wake_ups(void);
#endif

void set_units_sfr(void);

void gravity_forcetest(void);

void allocate_commbuffers(void);
void allocate_memory(void);
void begrun(void);
void check_omega(void);
void close_outputfiles(void);
void compute_accelerations(void);
void compute_global_quantities_of_system(void);
void compute_potential(void);
void construct_timetree(void);
void cooling_and_starformation(void);
void cooling_only(void);
void count_hot_phase(void);
void delete_node(int i);
void density(void);
void density_decouple(void);
void determine_interior(void);
int dissolvegas(void);
void do_box_wrapping(void);
double enclosed_mass(double R);
#ifndef LT_STELLAREVOLUTION
void endrun(int);
#else
void EndRun(int, const char *, const char *, const int);
#endif
void energy_statistics(void);
void ensure_neighbours(void);

void output_log_messages(void);
void ewald_corr(double dx, double dy, double dz, double *fper);

void ewald_force(int ii, int jj, int kk, double x[3], double force[3]);
void ewald_force_ni(int iii, int jjj, int kkk, double x[3], double force[3]);

void ewald_init(void);
double ewald_psi(double x[3]);
double ewald_pot_corr(double dx, double dy, double dz);
int find_ancestor(int i);
integertime find_next_outputtime(integertime time);
void find_next_time(void);
integertime find_next_time_walk(int node);
void free_memory(void);
void advance_and_find_timesteps(void);
integertime get_timestep(int p, double *a, int flag);

void determine_PMinterior(void);
void gravity_tree(void);
void hydro_force(void);
void init(void);
#ifndef LT_STELLAREVOLUTION
double get_starformation_rate(int i);
void init_clouds(void);
void integrate_sfr(void);
#else
double get_starformation_rate(int i, float *Temperature, float *xclouds);
void init_clouds(int, double, double, double, double, double*, double*);
void integrate_sfr(double, double, double, double, double, double);
#endif
void insert_node(int i);
int mark_targets(void);
size_t my_fwrite(void *ptr, size_t size, size_t nmemb, FILE * stream);
size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE * stream);
void mpi_printf(const char *fmt, ...);
void open_outputfiles(void);
void write_outputfiles_header(void);
void peano_hilbert_order(void);
double pot_integrand(double xx);
void predict(double time);
void predict_collisionless_only(double time);
void predict_sph_particles(double time);
void prepare_decouple(void);
void read_ic(char *fname);
void read_ic_cluster(char *fname);
void read_ic_cluster_gas(char *fname);
void read_ic_cluster_wimp(char *fname);
int read_outputlist(char *fname);
void read_parameter_file(char *fname);
void rearrange_particle_sequence(void);
void reorder_gas(void);
void reorder_particles(void);
void restart(int mod);
void run(void);
void savepositions(int num);
void savepositions_ioformat1(int num);
double second(void);
void set_softenings(void);
void set_sph_kernel(void);
void set_units(void);
void setup_smoothinglengths(void);

void minimum_large_ints(int n, long long *src, long long *res);
void sumup_large_ints(int n, int *src, long long *res);
void sumup_longs(int n, long long *src, long long *res);

void statistics(void);
double timediff(double t0, double t1);
void veldisp(void);
void veldisp_ensure_neighbours(int mode);

void gravity_tree_shortrange(void);


double get_hydrokick_factor(integertime time0, integertime time1);
double get_gravkick_factor(integertime time0, integertime time1);
double drift_integ(double a, void *param);
double gravkick_integ(double a, void *param);
double growthfactor_integ(double a, void *param);
double hydrokick_integ(double a, void *param);
void init_drift_table(void);
double get_drift_factor(integertime time0, integertime time1);
#ifdef MAGNETIC
double get_magkick_factor(integertime time0, integertime time1);
#endif
double measure_time(void);
double report_time(void);

/* on some DEC Alphas, the correct prototype for pow() is missing,
   even when math.h is included ! */

double pow(double, double);


void long_range_init(void);
void long_range_force(void);
void pm_init_periodic(void);
void pmforce_periodic(int mode, int *typelist);
void pm_init_regionsize(void);
void pm_init_nonperiodic(void);
int pmforce_nonperiodic(int grnr);

int pmpotential_nonperiodic(int grnr);
void pmpotential_periodic(void);

void readjust_timebase(double TimeMax_old, double TimeMax_new);

double enclosed_mass(double R);
void pm_setup_nonperiodic_kernel(void);

#if defined(CHEMISTRY) || defined(UM_CHEMISTRY)
double dmax(double, double);
double dmin(double, double);
#endif

#ifdef LT_STELLAREVOLUTION
void               fsolver_error_handler     (const char *, const char *, int, int);
int                get_Yset                  (int);

void               ReadYields                (int, int, int);
void               reading_thresholds_for_thermal_instability   (void);
double             GetMetalLambda            (double, double);

void               init_SN                   (void);
unsigned int       evolve_SN                 (void);
double INLINE_FUNC get_NextChemTime          (double, int, int*);
double INLINE_FUNC get_cost_SE               (int);
int    INLINE_FUNC get_chemstep              (int, int, double*, double);
int    INLINE_FUNC get_current_chem_bin      (double, int);
int                get_chemstep_bin          (double, double, int*, int);
int                compare_steps             (const void*, const void*);

int                append_chemicallyactive_particles(unsigned int *);
void               drop_chemicallyactive_particles  (void);


void               read_metals               (void);

void               init_clouds_cm            (int, double*, double*, double, double, double, int, double*);

int                load_SFs_IMFs             (void);
int                get_SF_index              (int, int*, int*);

int                write_eff_model           (int, int);
int                read_eff_model            (int, int);

void               read_SnIa_yields          (void);
void               read_SnII_yields          (void);
void               read_AGB_yields           (void);

void               read_metalcool_table      (void);

void               initialize_star_lifetimes (void);
double INLINE_FUNC get_age                   (double);

void               recalculate_stellarevolution_vars(void *, int);
void               recalc_eff_model          (void);

int                calculate_effective_yields(double, double, int, int);
double             calculate_FactorSN        (double, double, void*);

double             get_imf_params            (double);
double INLINE_FUNC normalize_imf             (double, double, void*);

double             *get_meanZ(void);
void               write_metallicity_stat    (void);
void               get_metals_sumcheck       (int mode);
double INLINE_FUNC get_metalmass             (float *);
double INLINE_FUNC get_metallicity           (int, int);
#ifdef SUBFIND
double INLINE_FUNC get_metallicity_subfind   (int);
double             get_starformation_rate_subfind(int i, float *Temperature, float *xclouds);
int                get_SF_index_subfind      (int, int*, int*);
#endif
double INLINE_FUNC get_metallicity_solarunits(MyFloat);

double INLINE_FUNC get_entropy               (int);

int    INLINE_FUNC getindex                  (double*, int, int, double*, int*);

int    INLINE_FUNC perturb                   (double*, double*);
int                is_chemically_active      (int);
#endif

#if defined (UM_METAL_COOLING)/* not necessarily UM_CHEMISTRY */

#ifdef UM_e_MET_IMPACTS
double um_metal_line_cooling(double, double, int, double);
#else
double um_metal_line_cooling(double, double, int);
#endif

#endif

#if defined(CHEMISTRY)
int compute_abundances(int mode, int ithis, double a_start, double a_end);
#endif

#if defined(CHEMISTRY) || defined(UM_CHEMISTRY)
int InitChem(void);
int init_rad(double);
#endif


#ifdef UM_CHEMISTRY

/* main call of chemistry treatment */
int compute_abundances(int mode, int ithis, double a_start, double a_end, double *um_energy);

/* main prototypes of chemistry treatment */
double INLINE_FUNC Um_Compute_MeanMolecularWeight(int);
double Um_AbundanceRatios(double, double, double*, double*, double*, double, double, int);


#if (defined(SFR) || defined(COOLING)) && !(defined(LT_METAL_COOLING) || defined(LT_METAL_COOLING_WAL))
double Um_DoCooling(double, double, double, double*, int, int);    // std SFR with um_chemistry
double Um_GetCoolingTime(double, double, double*, int); // std SFR + um_chem
#else
#ifndef LT_METAL_COOLING_WAL
double Um_DoCooling(double, double, double, double*, double, int, int);
double Um_GetCoolingTime(double, double, double*, double, int);
#else
double Um_DoCooling(double, double, double, double*, double*, double, double, int, int);
double Um_GetCoolingTime(double, double, double*, double*, double, double, int);
#endif
double Um_AbundanceRatios(double, double, double*, double*, double*, double, double, int);

#endif  // close std SFR with um_chemistry


/* additional prototypes of chemistry treatment */
double INLINE_FUNC Um_Compute_basenumden(double, int);
double INLINE_FUNC Um_Compute_n_den(double, int);
double INLINE_FUNC Um_Compute_sum_nuclei(int);
double INLINE_FUNC Um_Compute_mean_nden(int);

#if defined (UM_METAL_COOLING) && defined (LT_METAL_COOLING) && defined (UM_MET_IN_NONEQ_COOLING) && defined (UM_MET_NONEQ_CORRECTIONS)
void Um_met_corrections(int);
#else
void Std_corrections(int);
#endif

#if defined (UM_METAL_COOLING) && defined (LT_METAL_COOLING) && defined (UM_MET_IN_NONEQ_COOLING)
double INLINE_FUNC Um_Compute_invden(int i);
#endif

/* these fcts are not passed anywhere outside chemistry (see um_chemistry_noneq.c) */
//int energy_solver(int, int, double, double double, double, double *, double *, double *); 
//int energy_solver(int, int, double, double, double, double, double *, double *)
//int species_solver(int, int, double, double, double, double, double *);
//int radiative_rates(int, double, double);

#endif

#ifdef UM_CHECK
void Um_cooling_check(void);
#endif



#ifdef RADTRANSFER
void rt_get_J_nu(int ithis);

void rt_write_stats(void);

/* Eddington tensor computation */
int eddington_treeevaluate(int target, int mode, int *nexport, int *nsend_local);
void eddington(void);

int n_treeevaluate(int target, int mode, int *nexport, int *nsend_local);
void n(void);

/* radiative transfer integration */
void radtransfer(void);
double radtransfer_vector_multiply(double *a, double *b);
double radtransfer_vector_sum(double *a);
void radtransfer_matrix_multiply(double *in, double *out, double *sum);
int radtransfer_evaluate(int target, int mode, double *in, double *out, double *sum, int *nexport, int *nsend_local);
void radtransfer_set_simple_inits(void);
void radtransfer_update_chemistry(void);

void rt_get_sigma(void);
void rt_get_lum_stars(void);
void rt_get_lum_gas(int target, double *je);

double rt_GetCoolingTime(int i, double u, double rho);
double radtransfer_cooling_photoheating(int i, double dt);
double rt_DoCooling(int, double);
double rt_DoHeating(int, double);
double rt_get_cooling_rate(int i, double entropy);

/* emission approximation */
int ngb_treefind_stars(MyDouble searchcenter[3], MyFloat hsml, int target, int *startnode, int mode,
		       int *nexport, int *nsend_local);
void density_sfr(void);
int density_sfr_evaluate(int target, int mode, int *nexport, int *nsend_local);

void sfr_lum(void);
int sfr_lum_evaluate(int target, int mode, int *nexport, int *nsend_local);

void gas_lum(void);
void star_lum(void);
int star_lum_evaluate(int target, int mode, int *nexport, int *nsend_local);

void bh_lum(void);
int bh_lum_evaluate(int target, int mode, int *nexport, int *nsend_local);

#endif //end radtransfer


#ifdef AUTO_SWAP_ENDIAN_READIC
void swap_Nbyte(char *data, int n, int m);
void swap_header(void);
#endif
void find_block(char *label,FILE *fd);

#ifdef DISTORTIONTENSORPS
void get_half_kick_distortion(int pindex, MyBigFloat half_kick_add[6][6]);
void analyse_phase_space(int pindex,
                         MyBigFloat *s_1, MyBigFloat *s_2, MyBigFloat *s_3,
                         MyBigFloat *smear_x, MyBigFloat *smear_y, MyBigFloat *smear_z,
                         MyBigFloat *second_deriv, MyBigFloat *sigma);
MyBigFloat get_analytic_annihilation(MyBigFloat s_1, MyBigFloat s_2, MyBigFloat s_3,
                                   MyBigFloat s_1_prime, MyBigFloat s_2_prime, MyBigFloat s_3_prime,
                                   MyBigFloat second_deriv, MyBigFloat second_deriv_prime,
                                   MyBigFloat dt, MyBigFloat sigma);
MyBigFloat get_max_caustic_density(MyBigFloat s_1, MyBigFloat s_2,
                                 MyBigFloat s_1_prime, MyBigFloat s_2_prime,
                                 MyBigFloat second_deriv, MyBigFloat second_deriv_prime,
                                 MyBigFloat sigma,
                                 int pindex);
void get_current_ps_info(int pindex, MyBigFloat * flde, MyBigFloat * psde);
void do_phase_space_drift(int i, double dt_drift);
void do_the_phase_space_kick(int i, double dt_gravkick);
void do_long_range_phase_space_kick(int i, double dt_gravkick);
/* some math functions we need from phasespace_math.c */
void ludcmp(MyBigFloat ** a, int n, int *indx, MyBigFloat * d);
MyBigFloat **matrix(long nrl, long nrh, long ncl, long nch);
MyBigFloat *vector(long nl, long nh);
void free_matrix(MyBigFloat ** m, long nrl, long nrh, long ncl, long nch);
void free_vector(MyBigFloat * v, long nl, long nh);
void mult_matrix(MyBigFloat ** matrix_a, MyBigFloat ** matrix_b, int dimension, MyBigFloat ** matrix_result);
void mult_matrix_transpose_A(MyBigFloat ** matrix_a, MyBigFloat ** matrix_b, int dimension, MyBigFloat ** matrix_result);
void luinvert(MyBigFloat ** input_matrix, int n, MyBigFloat ** inverse_matrix);
void eigsrt(MyBigFloat d[], MyBigFloat ** v, int n);
void jacobi(MyBigFloat ** a, int n, MyBigFloat d[], MyBigFloat ** v, int *nrot);
#ifdef PMGRID
void pmtidaltensor_periodic_diff(void);
void pmtidaltensor_periodic_fourier(int component);
int pmtidaltensor_nonperiodic_diff(int grnr);
int pmtidaltensor_nonperiodic_fourier(int component, int grnr);
void check_tidaltensor_periodic(int particle_ID);
void check_tidaltensor_nonperiodic(int particle_ID);
#endif
#endif

#ifdef SINKS
void do_sinks(void);
#endif

#if defined(HEALPIX)

void mk_xy2pix(int *x2pix, int *y2pix);
void vec2pix_nest( const long nside, double *vec, long *ipix);
void healpix_halo_cond(float *res);

#endif

#ifdef CHEMCOOL
void coolinmo_(void);
void cheminmo_(void);
void init_tolerances_(void);
void load_h2_table_(void);
void init_temperature_lookup_(void);
double do_chemcool(int part_index, double dt);
double evolve_abundances_(double* dt, double* dl, double* yn, double* divv, double* energy, double* abundances, double* column_est);
void rate_eq_(int* nsp, double* t_start, double* y, double* ydot, double* rpar, int* ipar);
#endif

#ifdef JD_DPP
void compute_Dpp();
#endif


#ifdef SCFPOTENTIAL
void SCF_do_center_of_mass_correction(double fac_rad, double start_rad, double fac_part, int max_iter);
void SCF_write(int task);
void SCF_calc_from_random(long *seed);
void SCF_calc_from_particles(void);
void SCF_init(void);
void SCF_reset(void);
void SCF_free(void);
void SCF_evaluate(MyDouble x, MyDouble y, MyDouble z, MyDouble *potential, MyDouble *ax, MyDouble *ay, MyDouble *az);
void SCF_collect_update(void);

void sphere_acc(double x, double y, double z, double *xa, double *ya, double *za);
void to_unit(double x, double y, double z, double *xs, double *ys, double *zs);
double ran1(long *idum);
double gasdev(long *idum);
double factrl(int n);
int nlm_all(int num, int n, int l, int m);
int nlm(int n, int l, int m);
int nl(int n, int l);
int lm(int l, int m);
int kdelta(int a, int b);
double gnlm_var(int n, int l, int m);
double hnlm_var(int n, int l, int m);
#endif


#ifdef MOL_CLOUDS
int do_mol_clouds();
void init_mol_clouds();
#endif

#ifdef BP_REAL_CRs
void bp_cr_evol();
#endif

#ifdef SIDM
void sidm_Init_Particles(void);
void sidm_Init_CrossSection(void);
void sidm_DoScatter(void);
inline MyDouble sidm_cross_sigma(MyDouble rel_vel);
#endif

#ifdef ADAPTGRAVSOFT
void ags_setup_smoothinglengths(void);
int ags_ngb_treefind_variable_threads(MyDouble searchcenter[3], MyFloat hsml, int target, int *startnode,
				      int mode, int *exportflag, int *exportnodecount, int *exportindex,
				      int *ngblist);
void ags_density(void);
int ags_density_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist);
void *ags_density_evaluate_primary(void *p);
void *ags_density_evaluate_secondary(void *p);
int ags_density_isactive (int n);
void ags_force_update_hmax(void);
#endif

#ifdef ALTERNATIVE_PSORT
void init_sort_ID(MyIDType *data, int ndata);
#endif

#ifdef STATICBRANDT
inline double OmegaR(int i, int mode);
#endif

#ifdef DISKPOT
void gravity_tree_subfunc_diskpot(void);
#endif

#ifdef KD_ACTIVE_PARTICLE_LIST_FOR_OPENMP
void build_active_particle_list(void);
#endif

