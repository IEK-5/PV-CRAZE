
void sc2fit(double xi[],double yi[],double sdi[],long nxy,double b[],long nb,double *w,long nw,double yknot[],double ypknot[],double *sigfac,long *ierr1);

#define SPLNFIT_NBINV 			1	// nb < 2
#define SPLNFIT_INSUFFDATA 		2	// nxy < nb + 2
#define SPLNFIT_KNOTSSTRICINC	3	// knot x coordinates not strictly increasing
#define SPLNFIT_INSUFFWS		4	// insufficient workspace allocated
#define SPLNFIT_SUNEQZERO		5	// S[0] must not be 0
#define SPLNFIT_STRICTINC		6	// Xi's not strictly increasing
#define SPLNFIT_XIOUTRANGE		7	// Xi's not within knot range
#define SPLNFIT_ICRFWS			8	// increase workspace
#define SPLNFIT_SBACC			9	// sbacc fail
#define SPLNFIT_SNOTGTZERO		10	// require S[1+]>0 if S[0]>0
#define SPLNFIT_SINGULAR		11	// Singular system


