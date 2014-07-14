#ifndef INLINE_FUNC
#ifdef INLINE
#define INLINE_FUNC inline
#else
#define INLINE_FUNC
#endif
#endif

double AbundanceRatios(double u, double rho, double *ne_guess, double *nH0_pointer, double *nHeII_pointer);
void   find_abundances_and_rates(double logT, double rho, double *ne_guess);
void   InitCool(void);
void   InitCoolMemory(void);
void   IonizeParams(void);
void   IonizeParamsFunction(void);
void   IonizeParamsTable(void);
double INLINE_FUNC LogTemp(double u, double ne);
void   MakeCoolingTable(void);
void   ReadIonizeParams(char *fname);
void   SetZeroIonization(void);
void   TestCool(void);

#if !defined(LT_METAL_COOLING) && !defined(LT_METAL_COOLING_WAL)

double convert_u_to_temp(double u, double rho, double *ne_guess);
double CoolingRate(double logT, double rho, double *nelec);
double CoolingRateFromU(double u, double rho, double *ne_guess);
double DoCooling(double u_old, double rho, double dt, double *ne_guess);
double GetCoolingTime(double u_old, double rho,  double *ne_guess);
double DoInstabilityCooling(double m_old, double u, double rho, double dt, double fac, double *ne_guess);

#else

#if defined(LT_METAL_COOLING)
double convert_u_to_temp(double u, double rho, double *ne_guess);
double GetCoolingTime(double u_old, double rho,  double *ne_guess, double Z, double *temp);
double CoolingRateFromU(double u, double rho, double *ne_guess, double Z, double *temp);
double DoCooling(double u_old, double rho, double dt, double *ne_guess, double Z, double *temp);
double CoolingRate(double logT, double rho, double *nelec, double Z);
double DoInstabilityCooling(double m_old, double u, double rho, double dt, double fac, double *ne_guess, double Z, double *temp);
#endif

#if defined(LT_METAL_COOLING_WAL)
double CoolingRate(double Temp, double Redshift, double *Metallicities);
double CoolingRateFromU(double U, double Redshift, double DZ, double *Metallicities, double *temp);
double DoCooling(double U_in, double Rho, double *Metallicities, double Redshift, double DZ, double dt, double *temp);
double convert_u_to_temp(double U, double Redshift, double DZ, double *Metallicities);
double GetCoolingTime(double U, double Rho, double Redshift, double DZ, double *Metallicities, double *temp);
int    get_cool_redshift(double, double*);
double get_max_cool_redshift();
int    get_cool_n_el();
int    Is_a_Coolant(int);
void   WalCool_set_PID(MyIDType);
void   WalCool_tables_load(double);
void   WalCool_get_collis_table();
void   set_cooltable_index(int);
double *set_metallicities(int, double*, double);
#if defined(SUBFIND)  
double *set_metallicities_subfind(int, double*, double);
#endif
#endif

#endif



