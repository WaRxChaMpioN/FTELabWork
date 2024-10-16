/*--------------------------------------*- C++ -*----------------------------------------------------------*\
| ***** ***** *****	     R adiation     | Code developed in                                             |
|   *     *     * 	     A nd     	    | Radiation and Fluid Flow Physics Laboratory		    |
|   *     *     *            F luid	    |   					                    |
|   *     *     *	     F low          | Developers:-Dr. Pradeep, Chanakya, Ankur, Naman, Roshan,      |
|   *     *     *	     P hysics       | Kamal, Shreesh, Harsh, Pankaj, Abhishek                                 |
| ***** *****   *   	     L aboratory    |                                                		    |
\*---------------------------------------------------------------------------------------------------------*/
#include<stdio.h>
#include<stdlib.h>
#include<malloc.h>
#include<math.h>
main()
{
	FILE *fptr,*fptr1,*fptr2,*fp1,*fp2;
	int ni =22, nj=22, nt = 8, np=16,Tr=0,Tl=0,Tb=500,Tt=0,i,j,l,m,d,l1,m1;
	double Int_old[ni][nj][nt][np],Int_new[ni][nj][nt][np],S[ni][nj][nt][np],B[ni][nj][nt][np],*D_x,*D_y,*D_omega,dx,dy,dtheta,dphi,theta,phi,rms,Sum_n,I_w,k=0,s=1;
	double Nume,Deno,Source,Sum_o,flux_bot[ni],ctr_temp[nj],PI,q_w, Tg=1000,sigma=5.67e-8;
	PI = 4.0*atan(1.0);
	fp1=fopen("Heat_Flux_Isotropic_Scattering_Twalls_(K=0,S=1).csv","w");
	fp2=fopen("Centerline_Temp_Isotropic_Scattering_Twalls_(K=0,S=1).csv","w");
	if(fp1==NULL)
	printf("Unable to open flux Bot file\n");
	if(fp2==NULL)
	printf("Unable to open Centerline Temp file\n");
	dx=1.0/(ni-2.0);
	dy=1.0/(nj-2.0);
	dtheta=PI/(nt);
	dphi=2.0*PI/np;

		D_x=(double *)malloc((nt*np)*sizeof(double));
		D_y=(double *)malloc((nt*np)*sizeof(double));
		D_omega=(double *)malloc((nt*np)*sizeof(double));
		for(l=0;l<nt;l++)
			{
			 theta=dtheta*(l+0.5);
			 for(m=0;m<np;m++)
			 	{
				 phi=dphi*(m+0.5);
				 *(D_x+np*l+m)=cos(phi)*sin(dphi/2.0)*(dtheta-cos(2.0*theta)*sin(dtheta));
				 *(D_y+np*l+m)=sin(phi)*sin(dphi/2.0)*(dtheta-cos(2.0*theta)*sin(dtheta));
				 *(D_omega+np*l+m)=2.0*dphi*sin(theta)*sin(dtheta/2.0);
				}
			}




		 for(l=0;l<nt;l++)
			 for(m=0;m<np;m++)
				 for(j=0;j<nj;j++)
					 for(i=0;i<ni;i++)
						 Int_new[i][j][l][m]=0;

			do												//Convergence Loop
			 {
            for(l=0;l<nt;l++)
                for(m=0;m<np;m++)
                    for(j=0;j<nj;j++)
                        for(i=0;i<ni;i++)
                            S[i][j][l][m]=0;

            for(l=0;l<nt;l++)
                for(m=0;m<np;m++)
                    for(j=0;j<nj;j++)
                        for(i=0;i<ni;i++)
                            B[i][j][l][m]=0;

			  for(l=0;l<nt;l++)
			  	for(m=0;m<np;m++)
			  		for(j=0;j<nj;j++)
						for(i=0;i<ni;i++)
							Int_old[i][j][l][m]=Int_new[i][j][l][m];

			  for(l=0;l<nt;l++)                              //Scattering Integral Calculation
			  	for(m=0;m<np;m++)
			  		for(j=0;j<nj;j++)
						for(i=0;i<ni;i++){
                            for(l1=0;l1<nt;l1++){
                                for(m1=0;m1<np;m1++)
                                {


                                if(l1==l && m1==m)
                                {
                                    B[i][j][l][m]=0;//(*(D_omega+np*l1+m1))*(s/(4*PI));
                                }
                                else
                                {
                                    S[i][j][l][m]+=0;//(*(D_omega+np*l1+m1))*(s/(4*PI))*Int_old[i][j][l1][m1];
                                }
                                }
                            }
						}
////////////////////////////////////////////////Boundary Conditions///////////////////////////////////////////////////////////////////////////////////////////////////////

			for(i=1;i<=ni-1;i++)											//For Bottom Wall
				{

				for(l=0;l<nt;l++)
					for(m=0;m<np/2;m++)
						Int_new[i][0][l][m]=2.0*(sigma*pow(Tb,4)/PI)-(Int_new[i][1][l][m]);

				}

			for(i=1;i<=ni-1;i++)											//For Top Wall
				{

				for(l=0;l<nt;l++)
					for(m=np/2;m<np;m++)
					Int_new[i][nj-1][l][m]=2.0*(sigma*pow(Tt,4)/PI)-(Int_new[i][nj-2][l][m]);
				}

			for(j=1;j<=nj-1;j++)											//For Left Wall
				{


				for(l=0;l<nt;l++)
					for(m=3*np/2;m<np/2;m--)
						(Int_new[0][j][l][m])=2.0*(sigma*pow(Tl,4)/PI)-(Int_new[1][j][l][m]);

				}


			for(j=1;j<=nj-1;j++)											//For Right Wall
				{


				for(l=0;l<nt;l++)
					for(m=np/2;m>3*np/2;m++)
						(Int_new[ni-1][j][l][m])=2.0*(sigma*pow(Tr,4)/PI)-(Int_new[ni-2][j][l][m]);

				}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			for(l=0;l<nt;l++)
				for(m=0;m<np;m++)
					{
					 if(*(D_x+np*l+m)>0&&*(D_y+np*l+m)>0)
					 	{
						 for(j=1;j<nj;j++)
						 	for(i=1;i<ni;i++)
								{
								 Nume=dy*(fabs(*(D_x+np*l+m)))*(Int_new[i-1][j][l][m])+dx*(fabs(*(D_y+np*l+m)))*(Int_new[i][j-1][l][m]);
								 Source=(k*sigma*pow(Tg,4)/PI+S[i][j][l][m])*dx*dy*(*(D_omega+np*l+m));
								 Deno=dy*(fabs(*(D_x+np*l+m)))+dx*(fabs(*(D_y+np*l+m)))+((k+s)-B[i][j][l][m])*dx*dy*(*(D_omega+np*l+m));
								 Int_new[i][j][l][m]=(Nume+Source)/Deno;
							 	if(Int_new[i][j][l][m]<0)
								 		Int_new[i][j][l][m]=0.0;
								 }
						 }
					if(*(D_x+np*l+m)>0&&*(D_y+np*l+m)<0)
						{
						 for(j=nj-2;j>=0;j--)
						 	for(i=1;i<ni;i++)
								{
								 Nume=dy*(fabs(*(D_x+np*l+m)))*(Int_new[i-1][j][l][m])+dx*(fabs(*(D_y+np*l+m)))*(Int_new[i][j+1][l][m]);
								 Source=(k*sigma*pow(Tg,4)/PI+S[i][j][l][m])*dx*dy*(*(D_omega+np*l+m));
								 Deno=dy*(fabs(*(D_x+np*l+m)))+dx*(fabs(*(D_y+np*l+m)))+((k+s)-B[i][j][l][m])*dx*dy*(*(D_omega+np*l+m));
								 Int_new[i][j][l][m]=(Nume+Source)/Deno;
								 	if(Int_new[i][j][l][m]<0)
								 		Int_new[i][j][l][m]=0.0;
								 }
						 }
					if(*(D_x+np*l+m)<0&&*(D_y+np*l+m)>0)
						{
						 for(j=1;j<nj;j++)
						 	for(i=ni-2;i>=0;i--)
								{
								 Nume=dy*(fabs(*(D_x+np*l+m)))*(Int_new[i+1][j][l][m])+dx*(fabs(*(D_y+np*l+m)))*(Int_new[i][j-1][l][m]);
								 Source=(k*sigma*pow(Tg,4)/PI+S[i][j][l][m])*dx*dy*(*(D_omega+np*l+m));
								 Deno=dy*(fabs(*(D_x+np*l+m)))+dx*(fabs(*(D_y+np*l+m)))+((k+s)-B[i][j][l][m])*dx*dy*(*(D_omega+np*l+m));
								 Int_new[i][j][l][m]=(Nume+Source)/Deno;
										if(Int_new[i][j][l][m]<0)
								 		Int_new[i][j][l][m]=0.0;
								 }
						 }
					if(*(D_x+np*l+m)<0&&*(D_y+np*l+m)<0)
					 	{
						 for(j=nj-2;j>=0;j--)
						 	for(i=ni-2;i>=0;i--)
								{
								 Nume=dy*(fabs(*(D_x+np*l+m)))*(Int_new[i+1][j][l][m])+dx*(fabs(*(D_y+np*l+m)))*(Int_new[i][j+1][l][m]);
								 Source=(k*sigma*pow(Tg,4)/PI+S[i][j][l][m])*dx*dy*(*(D_omega+np*l+m));
								 Deno=dy*(fabs(*(D_x+np*l+m)))+dx*(fabs(*(D_y+np*l+m)))+((k+s)-B[i][j][l][m])*dx*dy*(*(D_omega+np*l+m));
								 Int_new[i][j][l][m]=(Nume+Source)/Deno;
							 	 if(Int_new[i][j][l][m]<0)
								 		Int_new[i][j][l][m]=0.0;
								 }
						 }
					 }
			         Sum_n=0.0;
				 Sum_o=0.0;
			for(l=0;l<nt;l++)
				for(m=0;m<np;m++)
					for(j=1;j<(nj-1);j++)
						for(i=1;i<(ni-1);i++)
							{
							Sum_n=Sum_n+Int_new[i][j][l][m];
							Sum_o=Sum_o+(Int_old[i][j][l][m]);
							 }
            rms=fabs((Sum_n-Sum_o))/fabs(Sum_n);

			}while((rms)>1.0e-6);

///////////////////////////////////////////////////////////////////////////// Heat Flux ///////////////////////////////////////////////////////////////////////////////

		for(i=1;i<(ni-1);i++)
			{
			 q_w=0.0;
			for(l=0;l<nt;l++)
				for(m=np/2;m<np;m++)
				{
					q_w=q_w+Int_new[i][1][l][m]*(*(D_y+np*l+m));
				}
			  flux_bot[i]=(q_w+sigma*pow(Tb,4))/(sigma*pow(Tb,4)-sigma*pow(Tr,4));
			 }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////// Centerline Temperature Distribution ///////////////////////////////////////////////////////////////////////////////

		for(j=1;j<(nj-1);j++)
			{
			 q_w=0.0;
			for(l=0;l<nt;l++)
				for(m=1/4.*np;m<3/4.*np;m++)
				{
					q_w=q_w-Int_new[10][j][l][m]*(*(D_x+np*l+m));// in each iteration q is subtracted with Int_new just to match the sign convention.
				}
			  ctr_temp[j]=((q_w)/(sigma)-pow(Tr,4))/(pow(Tb,4)-pow(Tr,4));
			 }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		for(i=1;i<(ni-1);i++)
			{
			fprintf(fp1,"%e\t, %e\n",(1./(ni-2))*(2.0*i-1.0)/2.0,flux_bot[i]);
			}
	fclose(fp1);
		for(j=1;j<(nj-1);j++)
			{
			fprintf(fp2,"%e\t, %e\n",(1./(nj-2))*(2.0*j-1.0)/2.0,ctr_temp[j]);
			}
	fclose(fp2);

	return 0;
}

