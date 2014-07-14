#include <stdio.h>
void output_compile_time_options(void)
{
printf(
"        PERIODIC\n"
"        PMGRID=512\n"
"        MULTIPLEDOMAINS=4\n"
"        PEANOHILBERT\n"
"        WALLCLOCK\n"
"        MYSORT\n"
"        NO_ISEND_IRECV_IN_DOMAIN\n"
"        FIX_PATHSCALE_MPI_STATUS_IGNORE_BUG\n"
"        LONGIDS\n"
"        DEBUG\n"
"\n");
}
