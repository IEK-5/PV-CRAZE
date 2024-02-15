
void nelmin ( double fn ( double x[], void * args ), void *args, int n, double start[], double xmin[], 
	double *ynewlo, double reqmin, double step[], int konvge, int kcount, 
	int *icount, int *numres, int *ifault );
