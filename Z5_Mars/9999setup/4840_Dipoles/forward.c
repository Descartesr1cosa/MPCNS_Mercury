/**************************************************************************/
/*  dbcomp_mars.c                                                         */
/*                                                                        */
/* Calculate the Br, Bt, Bp field on a regular grid at X km from an       */
/*  equivalent source magnetisation distribution.                         */
/*                                                                        */
/* langlais@ltpmail.gsfc.nasa.gov
/**************************************************************************/

/************************/
/****    includes    ****/
/************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <malloc.h>

/***********************/
/****    Defines    ****/
/***********************/


#define R_M    3393.5
#define R R_M

#define M_PI   3.14159265358979323846
#define PI_180 M_PI/180.

/***********************/
/****    Formats    ****/
/***********************/

#define READ_raw "%lf %lf %lf %lf %lf %lf %lf\n"
#define READ_bin "%lf %lf %lf %lf %lf %lf %d %lf %lf %lf\n"

#define READ_dip "%lf %lf %lf %lf %lf %lf\n"
#define WRITE_grd "%8.4lf %8.4lf %8.4lf %12.5lf %12.5lf %12.5lf\n"
#define WRITE_bin "%7.3lf %7.3lf %7.3lf %12.5lf %12.5lf %10.3lf %12.5lf %12.5lf %10.3lf %12.5lf %12.5lf %10.3lf\n"

#define WRITE_raw "%13.9lf %7.3lf %7.3lf %7.3lf %12.5lf %12.5lf %12.5lf %12.5lf %12.5lf %12.5lf\n"
#define WRITE_raw1 "%13.9lf %8.3lf %8.3lf %8.3lf %10.3lf %10.3lf %10.3lf\n"

/****    Main routine    *****/

int main(argc, argv)
int argc;
char *argv[];
      {
	FILE   *dip_in, *obs_in, *fout;    /* Files */
	char   str[150];
	double **B_R;                     /* synthetic Br */
	double **B_T;                     /* synthetic Bt */
	double **B_P;                     /* synthetic Bp */
	double thresh;                   /* max obs-dip distance */
	double thresh2;                   /* max obs-dip distance */
	int    N_DIP;                             /* number of dipoles */
	double spacg, deep; /* Dipole parameters */
	double *phi_dip, *tet_dip, *tck_dip; /* Dip. Position */
	double *ray_dip; /* Dip. Position */
	double *M_r_dip, *M_t_dip, *M_p_dip; /* Dip. magn. */
	double *M_f_dip, *inc_dip, *dec_dip; /* Dip. magn. */
	double dip1, dip2, dip3;
	double *cos_t_dip, *sin_t_dip;           /* computation variables */
	int    A, i, ii, j, jj, k;                       /* indices */
	double lat_min, lat_max, lon_min, lon_max; /* min/max for a grid */
	int    N_PT_LAT, N_PT_LON, N_PT;      /* Grid parameters    */
	double INC;
	double lat, lon, alt;
	double lat_obs, lon_obs, alt_obs;
	double dphi;
	int    DIP_FORM;
	int    INPUT_FORMAT;
	int    N_DATA; /* Input number of data per bin */
	int    N_DATA_MIN; /* associated minimum value */
	double ALT; /* output altitude */
	double dumb;                           /* temp variable */
	double tet_obs, phi_obs, ray_obs;      /* computation variables */
	double cos_t_obs, sin_t_obs;           /* computation variables */
	double sin_dphi, cos_dphi;             /* computation variables */
	double FR, FT, FP;                     /* computation variables */
	double A1, A2, A3;                     /* computation variables */
	double B1, B2, B3;                     /* computation variables */
	double C1, C2, C3;                     /* computation variables */
	double D1, D2, D3;                     /* computation variables */
	double E1, E2, E3;                     /* computation variables */
	double F1, F2, F3;                     /* computation variables */
 	double scaling, volume, C_E;           /* computation variables */
	double DD11, DD12, DD13;               /* computation variables */
	double DD21, DD22, DD23;               /* computation variables */
	double DD31, DD32, DD33;               /* computation variables */
	int    N_OBS; /* number of observations */

	/****    Introduction    ****/

	printf("####################################################\n");
	printf("\tComputation of Br, Bt, Bp on a grid from an \n");
	printf("\tequivalent source magnetization distribution\n\n");
	printf("####################################################\n");

	if (argc != 10)
	  {
	    printf("\t10 parameters must be defined! Type\n");
	    printf("\tforward_code.exe <Dipole_in> <output> <alt> <INC> ");
	    printf("<dipole format> <lat_min> <lat_max> <lon_min> <lon_max>\n");
	    printf("\t<Dip_in> is the input dipole file\n");
	    printf("\t<alt> is the output altitude (km)\n");
	    printf("\t<INC> is the grid increment (for synthetic grid)\n");
	    printf("\t<dipole format> is 1 (Mr, Mt, Mp are given) or\n");
	    printf("\t                   2 (Mf, inc, dec are given)\n");
	    printf("\t<lat_min>, <lat_max>, <lon_min>, <lon_max> are in degrees\n");
	    printf("\tExit!\n");
	    exit(-1);
	  }


	/****    File opening    ****/

	dip_in = fopen(argv[1], "r");
	if (dip_in == NULL)
	  {
	    printf("\tCan not open file <%s>. Exit!\n", argv[1]);
	    exit(-1);
	  }

	fout = fopen(argv[2], "r");
	if (fout != NULL)
	  {
	    printf("\n\tThe output file %s exists. Continue (y/n)[y]?",argv[2]);
	    fgets(str, 150, stdin);
	    if (str[0] == 'n' || str[0] == 'N' )
	      {
		exit(-1);
	      }
	    fclose(fout);
	  }

	fout = fopen(argv[2], "w");
	if (fout == NULL)
	  {
	    printf("\tCan not open file <%s>. Exit!\n", argv[2]);
	    exit(-1);
	  }

	/* input parameters */


	ALT    = atof(argv[3]); /* output altitude */
 
	thresh = 1500.; /* max obs-dip distance */
	thresh2 = thresh*thresh;

	INC    = atof(argv[4]); /* data increment */

	DIP_FORM = atoi(argv[5]); /* dipole format */
	if (DIP_FORM != 1 && DIP_FORM != 2)
	  {
	    printf("\tBad format option for dipoles. Please check.\n");
	    exit(-1);
	  }

	lat_min = atof(argv[6]);
	lat_max = atof(argv[7]);
	lon_min = atof(argv[8]);
	lon_max = atof(argv[9]);

	scaling = 100.;

	/****    Read the dipoles    ****/

	fscanf(dip_in, "%lf %d %lf\n", &spacg, &N_DIP, &deep);

	/* test altitude */

	if (ALT < spacg*(3393.5 * 2. * M_PI / 360.))
	  {
	    printf("\n\tThe minimum altitude should be lower than the man spacing of the model\n");
	    printf("\tMean spacing is %6.3lf deg * %6.3lf km/deg = %6.3lf km\n", spacg, (3393.5 * 2. * M_PI / 360.), spacg*(3393.5 * 2. * M_PI / 360.));
	    printf("\tRequested altitude is %6.3lf\n", ALT);
	    printf("\tDo you REALLY want to continue (y/n)[y]?");
	    fgets(str, 150, stdin);
	    if (str[0] == 'n' || str[0] == 'N' )
	      {
		fclose(fout);
		fclose(dip_in);
		exit(-1);
	      }
	  }
	
	phi_dip = (double *)malloc(N_DIP*sizeof(double));
	tet_dip = (double *)malloc(N_DIP*sizeof(double));
	tck_dip = (double *)malloc(N_DIP*sizeof(double));

	if (DIP_FORM == 1)
	  {
	    M_r_dip = (double *)malloc(N_DIP*sizeof(double));
	    M_t_dip = (double *)malloc(N_DIP*sizeof(double));
	    M_p_dip = (double *)malloc(N_DIP*sizeof(double));
	  }

	else if (DIP_FORM == 2)
	  {
	    M_f_dip = (double *)malloc(N_DIP*sizeof(double));
	    inc_dip = (double *)malloc(N_DIP*sizeof(double));
	    dec_dip = (double *)malloc(N_DIP*sizeof(double));
	  }

	ray_dip   = (double *)malloc(N_DIP*sizeof(double));
	sin_t_dip = (double *)malloc(N_DIP*sizeof(double));
	cos_t_dip = (double *)malloc(N_DIP*sizeof(double));

	i = 0;

	while(!feof(dip_in))
	  {
	    fscanf(dip_in, READ_dip, &lon, &lat, &tck_dip[i], 
				&dip1, &dip2, &dip3);
	    if (DIP_FORM == 1)
	      {
		M_r_dip[i] = dip1;
		M_t_dip[i] = dip2;
		M_p_dip[i] = dip3;
	      }
	    else if (DIP_FORM == 2)
	      {
		M_f_dip[i] = dip1;
		inc_dip[i] = dip2;
		dec_dip[i] = dip3;
	      }

	    tet_dip[i] = (90.-lat)*PI_180;
	    phi_dip[i] = lon * PI_180;
	    ray_dip[i] = R - deep;
	    cos_t_dip[i] = cos(tet_dip[i]);
	    sin_t_dip[i] = sin(tet_dip[i]);

	    i ++;
	  }

	if (N_DIP != i)
	  {
	    printf("\tWrong number of dipoles. Please check!\n\n");
	    exit (-1);
	  }

	fprintf(stdout,"\t\nNumber of dipole read : %d \n\n", N_DIP);
	fflush(stdout);

	/****    Initialisation    ****/

	k = 0;

	N_PT_LAT = (int)((lat_max - lat_min)/INC);
	N_PT_LON = (int)((lon_max - lon_min)/INC);
	B_R = (double **)malloc(N_PT_LAT * sizeof(double*));
	B_T = (double **)malloc(N_PT_LAT * sizeof(double*));
	B_P = (double **)malloc(N_PT_LAT * sizeof(double*));
	for (i=0;i<N_PT_LAT;i++)
	  {
	    B_R[i] = (double *)malloc(N_PT_LON * sizeof(double));
	    B_T[i] = (double *)malloc(N_PT_LON * sizeof(double));
	    B_P[i] = (double *)malloc(N_PT_LON * sizeof(double));
	  }
	for (i=0;i<N_PT_LAT;i++)
	  {
	    for (ii=0;ii<N_PT_LON;ii++)
	      {
		B_R[i][ii] = -9999.;
		B_T[i][ii] = -9999.;
		B_P[i][ii] = -9999.;
		k++;
	      }
	  }
	fprintf(stdout,"\tMemory allocated!\n");
	fflush(stdout);

	/****    Compute the field    ****/

	N_PT = k;
	fprintf(stdout,"\tComputation on a synthetic grid (%d points) :   0%%", N_PT);
	fflush(stdout);
	k = 0;
	
	for (i=0;i<N_PT_LAT;i++)
	  {
	    lat_obs = (double)i * INC + lat_min + INC/2.;
	    for (ii=0;ii<N_PT_LON;ii++)
	      {
		k ++;
		if (k%(N_PT/100) == 0)
		  {
		    fprintf(stdout, "\b\b\b\b%3d%%", k/(N_PT/100));
		    fflush(stdout);
		  }
		lon_obs = (double)ii * INC + lon_min + INC/2.;
		alt_obs = ALT;
		
		tet_obs = (90. - lat_obs) * M_PI / 180.;
		phi_obs = lon_obs * M_PI / 180.;
		ray_obs = R + alt_obs;
		/* cosinus and sinus, observation */
		cos_t_obs = cos(tet_obs);
		sin_t_obs = sin(tet_obs);
		/* Initialization */
		FR = 0.;
		FT = 0.;
		FP = 0.;
		for (j= 0; j< N_DIP; j ++)
		  {
		    /* dipole location*/
		    dphi       = phi_obs - phi_dip[j];
		    /* cosinus and sinus, dipole */
		    cos_dphi   = cos(dphi);
		    sin_dphi   = sin(dphi);
		    /* distance dipole-obs */
		    E2 = ray_obs*ray_obs + ray_dip[j]*ray_dip[j] - 
		      2.*ray_dip[j]*ray_obs*
		      (cos_t_obs*cos_t_dip[j]+
		       sin_t_obs*sin_t_dip[j]*cos_dphi);
		    /* test on distance */
		    if (E2 <= thresh2)
		      {
			E1 = sqrt(E2); /* =C */
			A1 =  cos_t_obs*cos_t_dip[j] + sin_t_obs*sin_t_dip[j]*cos_dphi;
			A2 = -sin_t_obs*cos_t_dip[j] + cos_t_obs*sin_t_dip[j]*cos_dphi;
			A3 = -sin_t_dip[j]*sin_dphi;
			
			B1 =  cos_t_obs*sin_t_dip[j] - sin_t_obs*cos_t_dip[j]*cos_dphi;
			B2 = -sin_t_obs*sin_t_dip[j] - cos_t_obs*cos_t_dip[j]*cos_dphi;
			B3 =  cos_t_dip[j]*sin_dphi;
			
			C1 =  sin_t_obs*sin_dphi;
			C2 =  cos_t_obs*sin_dphi;
			C3 =  cos_dphi;
			
			E3 = pow(E1, 3);
			
			volume = pow(R*2.*M_PI/360.*spacg,2) * tck_dip[j];
			/* volume = (59.2*59.2)*(1.89*1.89)*40.0; */
			C_E = (scaling * volume ) / E3;
			
			D1 = ( ray_obs - ray_dip[j]*A1) * 3.; 
			D2 = (         - ray_dip[j]*A2) * 3.; 
			D3 = (         - ray_dip[j]*A3) * 3.;
			
			F1 = ( ray_obs*A1 - ray_dip[j]) / E2;
			F2 = (-ray_obs*B1             ) / E2;
			F3 = ( ray_obs*C1             ) / E2;
			
			DD11 = (D1 * F1 - A1) * C_E;
			DD21 = (D2 * F1 - A2) * C_E;
			DD31 = (D3 * F1 - A3) * C_E;
			
			DD12 = (D1 * F2 + B1) * C_E;
			DD22 = (D2 * F2 + B2) * C_E;
			DD32 = (D3 * F2 + B3) * C_E;
			
			DD13 = (D1 * F3 - C1) * C_E;
			DD23 = (D2 * F3 - C2) * C_E;
			DD33 = (D3 * F3 - C3) * C_E;
			
			if (DIP_FORM == 1)
			  {
			    FR += DD11*M_r_dip[j]+DD12*M_t_dip[j]+DD13*M_p_dip[j];
			    FT += DD21*M_r_dip[j]+DD22*M_t_dip[j]+DD23*M_p_dip[j];
			    FP += DD31*M_r_dip[j]+DD32*M_t_dip[j]+DD33*M_p_dip[j];
			  }
			else if (DIP_FORM == 2)
			  {
			    FR += (DD11*sin(inc_dip[j]*PI_180)
				   + DD12*cos(inc_dip[j]*PI_180)*cos(dec_dip[j]*PI_180)
				   + DD13*cos(inc_dip[j]*PI_180)*sin(dec_dip[j]*PI_180)) * M_f_dip[j];
			    FT += (DD21*sin(inc_dip[j]*PI_180)
				   + DD22*cos(inc_dip[j]*PI_180)*cos(dec_dip[j]*PI_180)
				   + DD23*cos(inc_dip[j]*PI_180)*sin(dec_dip[j]*PI_180)) * M_f_dip[j];
			    FP += (DD31*sin(inc_dip[j]*PI_180)
				   + DD32*cos(inc_dip[j]*PI_180)*cos(dec_dip[j]*PI_180)
				   + DD33*cos(inc_dip[j]*PI_180)*sin(dec_dip[j]*PI_180)) * M_f_dip[j];
			  }
		      }
		  }
		B_R[i][ii] = FR;
		B_T[i][ii] = FT;
		B_P[i][ii] = FP;
	      }
	  }
	fprintf(stdout, " OK!\n");fflush(stdout);
	    
	printf("\tWriting the output file (%d points)         :   0%%", N_PT);
	k = 0;
	for (i=0;i<N_PT_LAT;i++)
	  {
	    lat = (double)i * INC + lat_min + INC/2.;
	    for (ii=0;ii<N_PT_LON;ii++)
	      {
		k ++;
		if (k%(N_PT/100) == 0)
		  {
		    fprintf(stdout, "\b\b\b\b%3d%%", k/(N_PT/100));
		    fflush(stdout);
		  }
		lon = (double)ii * INC + lon_min + INC/2.;
		alt = ALT;
		fprintf(fout, WRITE_grd, lon, lat, alt, B_R[i][ii], B_T[i][ii], B_P[i][ii]);
	      }
	  }
	fprintf(stdout, " OK!\n");fflush(stdout);
	fclose(dip_in);
	fclose(fout);
      }
