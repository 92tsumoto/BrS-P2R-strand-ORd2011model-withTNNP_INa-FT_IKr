#include <string.h>
#include "syspara.h"

int mode = 1;
int P = 8;
FILE *fopen(), *fpin, *fp0, *fp1, *fp2, *fp3, *fp4, *fp5;
FILE *fp6, *fp7, *fp8, *fp9, *fp10, *fp11, *fp12, *fp13, *fp14, *fp15;
FILE *fp16, *fp17, *fp18, *fp19, *fp20, *fp21, *fp22, *fp23, *fp24, *fp25;

typedef double Number;
typedef long long Lint;

main(int argc, char **argv)
{

	int i,k,m,l,z;
	int ii=0;
	long j;
	unsigned int count=0;
	double t = 0.0;
	double time;
	double h,R_all,V_all;
	double d1,d2;
	double v_thre,v_thre2;
	double cut;
	char *tmpname;
	char cmd[BUFSIZ];
	double tend;

/* Action Potential Duration and Max. Info */
	double **v_old,**v_old2,**dvdt_new;
	double ***vmax1,***vmax2; // Max. Voltage (mV)
	double ***dvdtmax1,***dvdtmax2; // Max. dv/dt (mV/ms)
	double ***apd1,***apd2; // Action Potential Duration
	double ***toneapd1,***toneapd2; // Time of dv/dt Max.
	double ***ttwoapd1,***ttwoapd2; // Time of 90% Repolarization
	double ***rmbp; // Resting Membrane Potential
	double ***nair; // Intracellular Na At Rest
	double ***cair; // Intracellular Ca At Rest
	double ***kir; // Intracellular K At Rest
	double ***caimax; // Peak Intracellular Ca
	double ***inamax; // Peak INa
	double **total_inamax; // Peak total INa
	double *total_ina; // Peak total INa within a cell
	double v_thre_time[beats][MEDIA_SITE][MEDIA_PATCH] ; // time when AP passes throght the threshold value
	double v_thre_time2[beats][MEDIA_SITE][MEDIA_PATCH] ; // time when AP passes throght the threshold value
	double diff[beats][MEDIA_SITE][MEDIA_PATCH] ; // time when AP passes throght the threshold value
	int thre_flag[MEDIA_SITE][MEDIA_PATCH] ; // time flag when AP passes throght the threshold value
	int thre_flag2[MEDIA_SITE][MEDIA_PATCH] ; // time flag when AP passes throght the threshold value
	int fflag[MEDIA_SITE][MEDIA_PATCH] ; // time flag when AP passes throght the threshold value


	tmpname = "temp";

	sprintf(cmd, "/usr/bin/cpp -P %s > %s", argv[1],tmpname);
	if(system(cmd) == -1){
		fprintf(stderr,"cannot open %s\n",argv[1]);
		exit(1);
	}
	if((fpin=fopen(tmpname,"r"))==NULL){
		fprintf(stderr,"cannot open %s\n",argv[1]);
		exit(1);
	}

	// prepare for output data files 
	if ((fp1 = fopen("status.out","w")) == NULL){
		printf("Can't open File\n");
		exit(1);
	}
	if ((fp2 = fopen("initdat.out","w")) == NULL){
		printf("Can't open File\n");
		exit(1);
	}
	if ((fp3 = fopen("ndata.out","w")) == NULL){
		printf("Can't open File\n");
		exit(1);
	}
	if ((fp4 = fopen("ndata_final.out","w")) == NULL){
		printf("Can't open File\n");
		exit(1);
	}
	if ((fp5 = fopen("act_time_1st.out","w")) == NULL){
		printf("Can't open File\n");
		exit(1);
	}
	if ((fp6 = fopen("act_time.out","w")) == NULL){
		printf("Can't open File\n");
		exit(1);
	}
	if ((fp7 = fopen("vm_data.out","w")) == NULL){
		printf("Can't open File\n");
		exit(1);
	}

	// input initial parameters
	printf("start parameter input\n");
	input_para(fpin);
	printf("end parameter input\n");

	if (var.write){
		if ((fp0 = fopen(argv[2],"w"))==NULL){
			fprintf(stderr, "%s cannot open.\n",argv[2]);
			exit(-1);
		}
	}

	for (ii = 0; ii < var.datas; ii++){
		time = 0.0;
		tend = var.tend[ii];

		// Eular's step size
		h = 1.0 / var.m;
		// AP threshold value
		v_thre = -40.0;
		v_thre2 = -10.0;
		// a variable for data output	
		cut = var.m/10.0;

		h *= var.tsign[ii];

	// initialized memory
		initial_mem();

		v_old = (double**)malloc2d(sizeof(Number),MEDIA_SITE,MEDIA_PATCH);
		v_old2 = (double**)malloc2d(sizeof(Number),MEDIA_SITE,MEDIA_PATCH);
		dvdt_new = (double**)malloc2d(sizeof(Number),MEDIA_SITE,MEDIA_PATCH);

		vmax1 = (double***)malloc3d(sizeof(Number),beats,MEDIA_SITE,MEDIA_PATCH);
		vmax2 = (double***)malloc3d(sizeof(Number),beats,MEDIA_SITE,MEDIA_PATCH);
		dvdtmax1 = (double***)malloc3d(sizeof(Number),beats,MEDIA_SITE,MEDIA_PATCH);
		dvdtmax2 = (double***)malloc3d(sizeof(Number),beats,MEDIA_SITE,MEDIA_PATCH);
		apd1 = (double***)malloc3d(sizeof(Number),beats,MEDIA_SITE,MEDIA_PATCH);
		apd2 = (double***)malloc3d(sizeof(Number),beats,MEDIA_SITE,MEDIA_PATCH);
		toneapd1 = (double***)malloc3d(sizeof(Number),beats,MEDIA_SITE,MEDIA_PATCH);
		toneapd2 = (double***)malloc3d(sizeof(Number),beats,MEDIA_SITE,MEDIA_PATCH);
		ttwoapd1 = (double***)malloc3d(sizeof(Number),beats,MEDIA_SITE,MEDIA_PATCH);
		ttwoapd2 = (double***)malloc3d(sizeof(Number),beats,MEDIA_SITE,MEDIA_PATCH);
		rmbp = (double***)malloc3d(sizeof(Number),beats,MEDIA_SITE,MEDIA_PATCH);
		nair = (double***)malloc3d(sizeof(Number),beats,MEDIA_SITE,MEDIA_PATCH);
		cair = (double***)malloc3d(sizeof(Number),beats,MEDIA_SITE,MEDIA_PATCH);
		kir = (double***)malloc3d(sizeof(Number),beats,MEDIA_SITE,MEDIA_PATCH);
		caimax = (double***)malloc3d(sizeof(Number),beats,MEDIA_SITE,MEDIA_PATCH);
		inamax = (double***)malloc3d(sizeof(Number),beats,MEDIA_SITE,MEDIA_PATCH);

		total_inamax = (double**)malloc2d(sizeof(Number),beats,MEDIA_SITE);
		//total_inamax = (double**)calloc(beats, sizeof(double*));

		total_ina = (double*)calloc(MEDIA_SITE, sizeof(double));
		var.total_ina = (double*)calloc(MEDIA_SITE, sizeof(double));
		if (total_ina==NULL || var.total_ina==NULL ) exit(1);

		printf("finished memory initialization\n");

	// initial values input.
		for (i = 0; i < NN; i++){ 
			for (m = 0; m < MEDIA_SITE; m++){ 
				for (l = 0; l < MEDIA_PATCH; l++){ 
					var.u[i][m][l] = var.x0[ii][i];
					var.inner_v[m][l] = var.u[0][m][l];
				}
			}
		}


	// input static parameters
	static_paras(fp1);
	printf("finished input of static parameters\n");

	if(var.out_data){
		if ((fp8 = fopen("ina_data.out","w")) == NULL){
			printf("Can't open File\n");
			exit(1);
		}
		if ((fp9 = fopen("ical_data.out","w")) == NULL){
			printf("Can't open File\n");
			exit(1);
		}
		if ((fp10 = fopen("ik1_data.out","w")) == NULL){
			printf("Can't open File\n");
			exit(1);
		}
		if ((fp11 = fopen("iks_data.out","w")) == NULL){
			printf("Can't open File\n");
			exit(1);
		}
		if ((fp12 = fopen("ikr_data.out","w")) == NULL){
			printf("Can't open File\n");
			exit(1);
		}
		if ((fp13 = fopen("ito_data.out","w")) == NULL){
			printf("Can't open File\n");
			exit(1);
		}
		if ((fp14 = fopen("inak_data.out","w")) == NULL){
			printf("Can't open File\n");
			exit(1);
		}
		if ((fp15 = fopen("incx_data.out","w")) == NULL){
			printf("Can't open File\n");
			exit(1);
		}
	}
	if(var.out_data_plus){
		if ((fp16 = fopen("cleft_data.out","w")) == NULL){
			printf("Can't open File\n");
			exit(1);
		}
		if ((fp17 = fopen("intra_v_data.out","w")) == NULL){
			printf("Can't open File\n");
			exit(1);
		}
	}
	if(var.out_data_plus2){
		if ((fp18 = fopen("nai_data.out","w")) == NULL){
			printf("Can't open File\n");
			exit(1);
		}
		if ((fp19 = fopen("ki_data.out","w")) == NULL){
			printf("Can't open File\n");
			exit(1);
		}
	}

	// for solving algebra equation.
	var.par_k = 9*(MEDIA_SITE-2)+6*2;
	printf("var.par_k = %d\n",var.par_k);
	printf("R matrix size = %d X %d\n",MAT_SIZE,MAT_SIZE);
	printf("k=%d, Mem_size=%lf MB\n",var.par_k,((var.par_k)*(int)sizeof(Number))/1024/1024);
	// setting of each element of Resistance Matrix
	var.ia=(Lint *)calloc(var.par_k,sizeof(Lint));
	if(var.ia==NULL){
		printf(" ia null \n");
		exit(1);
	}
	var.ja=(Lint *)calloc(var.par_k,sizeof(Lint));
	if(var.ja==NULL){
		printf(" ja null \n");
		exit(1);
	}
	var.AT=(Number *)calloc(var.par_k,sizeof(Number));
	if(var.AT==NULL){
		printf(" AT null \n");
		exit(1);
	}
	printf("coeff input start\n");
	linear_coeff(); // setting of each element of Resistance Matrix
	printf("coeff input end\n");
	var.mtype = 11; /* Real unsymmetric matrix */
	var.nrhs = 1; /* Number of right hand sides. */
	init_pardiso();
	printf("pardiso initiation ok\n");

	var.b=(Number *)calloc(MAT_SIZE,sizeof(Number));
	if(var.b==NULL){
		printf(" B matrix is null \n");
		exit(1);
	}

	var.solv=(Number *)calloc(MAT_SIZE,sizeof(Number));
	if(var.solv==NULL){
		printf(" solution matrix is null \n");
		exit(1);
	}

	// Tablize exp functions.       
	printf("start tablization\n");
	make_ExpTable();
	printf("finished tablization\n");

	if(var.drug_mode == 0){
		var.normal_block = 1.0/(1.0+pow((var.drug_conc/var.ic50),var.hill));
		var.block_rate = var.normal_block;
	} else if(var.drug_mode == 1){ 
		var.normal_block=var.block_rate;
		var.drug_conc = var.ic50*pow(1.0/var.block_rate-1.0,1/var.hill);
		fprintf(fp1,"concentration=%e\n",var.drug_conc);
	}

	// affected fraction (with facilitation)
	if(var.drug_type==0){ // without facilitation
		var.fraction_facil = 0.0;
		printf("without facilitation model\n");
	} else if(var.drug_type==1){ // with facilitation
		var.fraction_facil = 1.0/(1.0+pow((var.drug_conc/var.ec50),var.ehill));
		printf("with facilitation model\n");
	//} else if(var.drug_type==2){  // only facilitation
	//  var.fraction_facil = 1.0/(1.0+pow((var.drug_conc/var.ec50),var.ehill));
	//  printf("with facilitation model\n");
	//  if(var.model_type == 0){
	//      var.normal_block = 1.0;
	//      var.block_rate = var.normal_block;
	//  } else if(var.model_type == 1){
	//      var.normal_block=1.0;
	//      var.drug_conc = var.ic50*pow(1.0/var.block_rate-1.0,1/var.hill);
	//  }
	}

	printf("fraction=%lf,noaffected=%lf,affected=%lf,total=%lf\n",
	var.fraction_facil,var.block_rate*var.fraction_facil,var.block_rate*(1.0-var.fraction_facil),var.block_rate);
	fprintf(fp1,"fraction=%lf,noaffected=%lf,affected=%lf,total=%lf\n",
	var.fraction_facil,var.block_rate*var.fraction_facil,var.block_rate*(1.0-var.fraction_facil),var.block_rate);

	// Initialization time
	time -= h;
	var.dt = h;
	var.beat = 0;
	
	for(var.beat=0; var.beat < beats; var.beat++){
		if(var.beat==beats-1){
			var.l = 2.0*var.BCL;
		} else {
			var.l = var.BCL;
		}

		for (j = 0; j< (var.m * var.BCL ); j++){
			t = h*j;
			time += h;
			if(j==0){
				for(k=0;k<MEDIA_SITE;k++){
					for(m=0;m<MEDIA_PATCH;m++){
						fflag[k][m]=0;
						diff[var.beat][k][m]=3.0*var.BCL;
					}
				}
			}

			for(k=0;k<MEDIA_SITE;k++){
				for(m=0;m<MEDIA_PATCH;m++){
					v_old[k][m] = var.u[0][k][m];
					//v_old2[k][m] = var.inner_v[k][m];
					v_old2[k][m] = var.u[0][k][m];
				}
			}

			if ( time-(var.BCL*var.beat+10) >= 0.0 && time-(var.BCL*var.beat+10) < h ){
				for(k=0;k<MEDIA_SITE;k++){
					for(m=0;m<MEDIA_PATCH;m++){
						rmbp[var.beat][k][m] = var.u[0][k][m]; nair[var.beat][k][m] = var.u[33][k][m];
						kir[var.beat][k][m]  = var.u[35][k][m]; cair[var.beat][k][m] = var.u[37][k][m];
						thre_flag[k][m] = 0; thre_flag2[k][m] = 0; 
						vmax1[var.beat][k][m] = -60.0; vmax2[var.beat][k][m] = -60.0;
						//diff[var.beat][k][m]=1000.0;
					}
				}
			}

			if (time-(var.BCL*(double)var.beat+10.0) >= 0.0 && time-(var.BCL*(double)var.beat+10.0) < 1.0){
				var.Istim = var.Istim_base;
				//var.Istim = var.Istim_base/(var.s2[(int)floor(MEDIA_PATCH/2.0)]/(2*M_PI*var.a*var.a+2*M_PI*var.a*var.length));
				//printf("Ist=%lf\n",var.Istim);
			} else {
				var.Istim = 0.0;
			}

			if (fabs(time) > tend &&  tend != 0.0) break;

			eular(NN,h,t);

			for (k=0;k<MEDIA_SITE;k++){
				for (m=0;m<MEDIA_PATCH;m++){

					dvdt_new[k][m] = (var.u[0][k][m]-v_old[k][m])/h;

					if(var.beat>=0){
						if (var.u[37][k][m] > caimax[var.beat][k][m] ){
							caimax[var.beat][k][m] = var.u[37][k][m];
						}
						if (ina.total[k][m] < inamax[var.beat][k][m] ){
							inamax[var.beat][k][m] = ina.total[k][m];
						}
					
						if(fflag[k][m]==0){
							if (dvdt_new[k][m] > dvdtmax1[var.beat][k][m]){
								dvdtmax1[var.beat][k][m] = dvdt_new[k][m];
								toneapd1[var.beat][k][m] = time;
							}
							if (var.u[0][k][m] > vmax1[var.beat][k][m]){
								vmax1[var.beat][k][m] = var.u[0][k][m];
							}
							if (var.u[0][k][m] >= (vmax1[var.beat][k][m] -0.9*(vmax1[var.beat][k][m] -rmbp[var.beat][k][m] )) ){
								ttwoapd1[var.beat][k][m] = time;
							}
							if (diff[var.beat][k][m]==ttwoapd1[var.beat][k][m] && ttwoapd1[var.beat][k][m]!=0.0){
								fflag[k][m]=1;
							} else {
								diff[var.beat][k][m]=ttwoapd1[var.beat][k][m];
							}
						}

						if(fflag[k][m]>0){
							if(var.dvdt[k][m]>0){
								fflag[k][m]=2;
							}
							if(fflag[k][m]==2){
								if (dvdt_new[k][m] > dvdtmax2[var.beat][k][m] ){
									dvdtmax2[var.beat][k][m] = dvdt_new[k][m];
									toneapd2[var.beat][k][m] = time;
								}
								if (var.u[0][k][m] > vmax2[var.beat][k][m] ){
									vmax2[var.beat][k][m] = var.u[0][k][m];
								}
								if (var.u[0][k][m] >= (vmax2[var.beat][k][m] -0.9*(vmax2[var.beat][k][m] -rmbp[var.beat][k][m] )) ){
									ttwoapd2[var.beat][k][m] = time;
								}
							}
						}

						// chack for activation time
						//if( (v_old2[k][m] - v_thre)<0 && (v_old2[k][m] - v_thre)*(var.inner_v[k][m] - v_thre)<0){
						if( (v_old[k][m] - v_thre)<0 && (v_old[k][m] - v_thre)*(var.u[0][k][m] - v_thre)<0){
							if(thre_flag[k][m] != 1){
								v_thre_time[var.beat][k][m] = time;
								if(m==1){
									printf("1st pass thre unit[%d][%d] = %lf\n",k,m,v_thre_time[var.beat][k][m]);
									fprintf(fp5,"%d %d %d %lf\n",var.beat,k,m,v_thre_time[var.beat][k][m]);
								}
								thre_flag[k][m] = 1;
							}
						}

						if( (v_old2[k][m] - v_thre2)<0 && (v_old2[k][m] - v_thre2)*(var.u[0][k][m] - v_thre2)<0){
							if(thre_flag2[k][m] != 1){
								v_thre_time2[var.beat][k][m] = time;
								if(m==1){
									printf("2nd pass thre unit[%d][%d] = %lf\n",k,m,v_thre_time2[var.beat][k][m]);
									fprintf(fp6,"%d %d %d %lf\n",var.beat,k,m,v_thre_time2[var.beat][k][m]);
								}
								thre_flag2[k][m] = 1;
							}
						}
						// check end for activation time
					
					var.dvdt[k][m] = dvdt_new[k][m];

					}
				}
			}

			for (k=0;k<MEDIA_SITE;k++){
				for (m=0;m<MEDIA_PATCH;m++){ total_ina[k]+=ina.total[k][m]*var.s2[m]; }
					if (var.total_ina[k] < total_inamax[var.beat][k] ){ total_inamax[var.beat][k] = var.total_ina[k]; }
					var.total_ina[k]=total_ina[k];
					total_ina[k]=0;
			}

			if(var.out_data){
				if (time>= (beats-2)*var.BCL && time <= beats*var.BCL){
					if(j%10 == 0){ // To reduce the recording data.
						ina_data(fp8,time);
					}
					if(j%(int)cut == 0){ // To reduce the recording data.
						ical_data(fp9,time);
						ik1_data(fp10,time);
						iks_data(fp11,time);
						ikr_data(fp12,time);
						ito_data(fp13,time);
						inak_data(fp14,time);
						incx_data(fp15,time);
					}
				}
			}
			if(var.out_data_plus){
				if (time>= (beats-2)*var.BCL){
					//if(j%(int)cut == 0){ // To reduce the recording data.
					if(j%10 == 0){ // To reduce the recording data.
						intra_v_data(fp17,time);
					}
				}
				if (time>= (beats-2)*var.BCL){
					if(j%10 == 0){ // To reduce the recording data.
						cleft(fp16,time);
					}
				}
			}
			if(var.out_data_plus2){
				if (time>= (beats-2)*var.BCL){
					if(j%(int)cut == 0){ // To reduce the recording data.
						nai_data(fp18,time);
						ki_data(fp19,time);
					}
				}
			}
			if (time>= (beats-5)*var.BCL){
				if(j%(int)cut == 0){ // To reduce the recording data.
					vm_data(fp7,time);
				}
			}

		} // for-loop end; j

		// reset of total_ina
		for (k=0;k<MEDIA_SITE;k++){
			var.total_ina[k]=0.0;
		}

			// for output of initial values
			fprintf(fp2,"beats=%d\n",var.beat+1);
			for(m=0;m<MEDIA_SITE;m++){
				for(k=0;k<MEDIA_PATCH;k++){
					for(i=0;i<NN;i++){ fprintf(fp2,"%16.15e\n",var.u[i][m][k]);}
				}
			} fprintf(fp2,"\n");

		} // for-loop end; var.beats

		//printf("data out\n");
		for(z=0;z<beats;z++){
			for(m=0;m<MEDIA_SITE;m++){
				for(k=0;k<MEDIA_PATCH;k++){
					apd1[z][m][k] = ttwoapd1[z][m][k] - toneapd1[z][m][k] ;
					apd2[z][m][k] = ttwoapd2[z][m][k] - toneapd2[z][m][k] ;
				}
			}

			for(m=0;m<MEDIA_SITE;m++){
				for(k=0;k<MEDIA_PATCH;k++){
					fprintf(fp3,"%i\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%e\n",
						z+1,vmax1[z][m][k],dvdtmax1[z][m][k],v_thre_time[z][m][k],v_thre_time2[z][m][k],toneapd1[z][m][k],ttwoapd1[z][m][k],apd1[z][m][k],
							vmax2[z][m][k],dvdtmax2[z][m][k],toneapd2[z][m][k],ttwoapd2[z][m][k],apd2[z][m][k],rmbp[z][m][k],inamax[z][m][k]);
				}
			}

			if(z==beats-1){
				for(m=0;m<MEDIA_SITE;m++){
					for(k=0;k<MEDIA_PATCH;k++){
					//fprintf(fp5,"%d-%d-%d\n",beats-1,m,k);
						fprintf(fp4,"%i\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%e\n",
							beats,vmax1[z][m][k],dvdtmax1[z][m][k],v_thre_time[z][m][k],v_thre_time2[z][m][k],
							toneapd1[z][m][k],ttwoapd1[z][m][k],apd1[z][m][k],vmax2[z][m][k],dvdtmax2[z][m][k],
							toneapd2[z][m][k],ttwoapd2[z][m][k],apd2[z][m][k],rmbp[z][m][k],inamax[z][m][k]);
					}
				}
			}
		}

	} // for-loop end; ii

/* -------------------------------------------------------------------- */
/* .. Termination and release of memory. */
/* -------------------------------------------------------------------- */
	release_of_memory();

	fclose(fp1); fclose(fp2); fclose(fp3); fclose(fp4); fclose(fp5);
	fclose(fp6); fclose(fp7); 
	if(var.out_data){
		fclose(fp8); fclose(fp9); fclose(fp10); fclose(fp11); fclose(fp12); fclose(fp13);
		fclose(fp14); fclose(fp15);
	}
	if(var.out_data_plus){
		fclose(fp16); fclose(fp17);
	}
	if(var.out_data_plus2){
		fclose(fp18); fclose(fp19);
	}

	free(v_old); free(v_old2); free(dvdt_new);

	free(vmax1); free(vmax2); free(dvdtmax1); free(dvdtmax2); 
	free(apd1); free(apd2); free(toneapd1); free(toneapd2); free(ttwoapd1); free(ttwoapd2);
	free(rmbp); free(nair); free(cair); free(kir); free(caimax); free(inamax);

	free(total_inamax); free(total_ina); free(var.total_ina);

	free(var.ia);free(var.ja);free(var.AT);free(var.b);free(var.solv);

	close_mem();

} // end main

