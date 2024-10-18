/*--------------------------------------*- C++ -*----------------------------------------------------------*\
| ***** ***** *****	     R adiation     | Code developed in                                             |
|   *     *     * 	     A nd     	    | Radiation and Fluid Flow Physics Laboratory		    |
|   *     *     *        F luid	   		|   					                    |
|   *     *     *	     F low          | Developers:-Dr. Pradeep, Chanakya, Ankur, Naman, Roshan,      |
|   *     *     *	     P hysics       | Kamal, Shreesh, Harsh, Pankaj                                 |
| ***** *****   *   	 L aboratory    |                                                		    |
\*---------------------------------------------------------------------------------------------------------*/
#include<stdio.h>
#include<stdlib.h>
#include<malloc.h>
#include<math.h>
#include<unistd.h>
#include<time.h>

//#include<cstring>
int main()
{

    FILE *fptr,*fptr1,*fptr2,*fptr3,*fp1,*fp2,*fp3,*fp5,*fp6,*fp7,*fp8;
    int A1=1,sia1,p,ni =22, nj=22, nt =6 , np=24,Tr=0,Tl=0,Tb=0,Tt=1000,i,j,l,m,l1,m1,d,N=1000000,count=0;
    double Int_new[ni][nj][nt][np],*D_x,*D_y,*D_omega,dx,dy,dtheta,dphi,theta,phi,rms,Sum_n,*k,*k1,*kp,*w,*a,*s,I_w,coun=0;
    double Nume,Deno,Source,Sum_o,flux_bot[ni],flux_bot_t[ni],flux_bot_l[ni],flux_bot_r[ni],div_q[ni],dq, PI,q_w, Tg=0,sigma=5.6696e-8,aa=0,AA,q2=0,q3=0,q_w1=0,q_w6=0,q_w7=0,q_w8=0;
    double Int_new1[ni][nj][nt][np],Iold[ni][nj][nt][np],sum1=0,sum2=0,rms1,urf=0.1,small=1e-10;
    
    double S[ni][nj][nt][np],S_old[ni][nj][nt][np],B[ni][nj][nt][np],phase[200];
    //memset( Int_new1, 0, ni*nj*nt*np*sizeof(double) );
    //double as[ni][nj][nt][np]={0};
    double as[ni][nj][nt][np];
    double d_omg,ang_sum,theta1,phase_fun,phase_funsum=0.0,phasefun,sia;
    k = (double*)calloc(N,sizeof(double*));
    k1 = (double*)calloc(N,sizeof(double*));
    w = (double*)calloc(N,sizeof(double*));
    kp = (double*)calloc(N,sizeof(double*));
    s = (double*)calloc(N,sizeof(double*));
    fptr2= fopen("kappa1e5.dat","r");
    fptr3= fopen("sigma1e5.dat","r");
    fp5=fopen("co2_k_1000_005.dat","r");


    //fptr = fopen("K_H2O1000_particle505e_2.dat", "r");
    fptr = fopen("H2O_k_1000_010.dat", "r");
    fptr1 = fopen("s_1000.dat", "r"); // "w" defines "writing mode"
    for (d = 1; d<=N; d++) 
    {
        fscanf(fptr,"%lf", &k[d]);
        fscanf(fp5,"%lf", &k1[d]);
        fscanf(fptr1,"%lf", &w[d]);
    }
    for (d = 1; d<(N-90000); d++) // kp length
    {
        fscanf(fptr2,"%lf", &kp[d]);
        fscanf(fptr3,"%lf", &s[d]);
        s[d]=100*s[d];
        //s[d]=0;
        //k[d]=0;//100*(k[d] + k1[d]) + kp[d];    // absorption in cm^-1
    }
    
    for (d = (N-90000); d<=N; d++) 
        s[d]=0;
    
    fclose(fptr);
    fclose(fptr1);
    fclose(fptr2);
    fclose(fptr3);
    PI = 4.0*atan(1.0);
    for(d=1;d<=N;d++)
    {
        k[d]=0;
        k1[d]=0;
    }
    fp1=fopen("flux_LBL_gas_part.dat","w");
    fp2=fopen("div_LBL_gas_part.dat","w");
    fp3=fopen("spectral_flux_bottom.dat","w");
    fp6=fopen("spectral_flux_left.dat","w");
    fp7=fopen("spectral_flux_right.dat","w");
    fp8=fopen("spectral_flux_top.dat","w");
    if(fp1==NULL)
        printf("Unable to open flux Bot file\n");
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
            for(j=1;j<nj-1;j++)
                for(i=1;i<ni-1;i++)
                {
                    Int_new[i][j][l][m]=0;//1123.479;
                    S[i][j][l][m]=0;
                    B[i][j][l][m]=0;
                    as[i][j][l][m]=0;
                    Iold[i][j][l][m]=0;
                }
    
    time_t begin = time(NULL);
        for(d=1; d<=N; d++) // Start of absorption coefficient
        {   
            // (sigma*pow(Tg,4))/PI;

            // intitalisation
            //////////////////////////////////////////////////    Boundary Conditions /////////////////////////////////////////////////////////////////////////////////////
            for(i=1; i<ni-1; i++) // Bottom Wall
            {
            for(l=0; l<nt; l++)
                for(m=0; m<np/2; m++)
                Int_new[i][0][l][m] = w[d] * (sigma * pow(Tb,4)) / PI; // -Int_new[i][1][l][m];
            }

            for(i=1; i<ni-1; i++) // Top Wall
            {
            for(l=0; l<nt; l++)
                for(m=np/2; m<np; m++)
                Int_new[i][nj-1][l][m] = w[d] * (sigma * pow(Tt,4)) / PI; // -Int_new[i][nj-2][l][m];
            }

            for(j=1; j<nj-1; j++) // Left Wall
            {
            for(l=0; l<nt; l++)
            {
                for(m=0; m<np/4; m++)
                Int_new[0][j][l][m] = w[d] * (sigma * pow(Tl,4)) / PI; // -Int_new[1][j][l][m];
                for(m=3*np/4; m<np; m++)
                Int_new[0][j][l][m] = w[d] * (sigma * pow(Tl,4)) / PI; // -Int_new[1][j][l][m];
            }
            }

            for(j=1; j<nj-1; j++) // Right Wall
            {
            for(l=0; l<nt; l++)
                for(m=np/4; m<3*np/4; m++)
                Int_new[ni-1][j][l][m] = w[d] * (sigma * pow(Tr,4)) / PI; // -Int_new[ni-2][j][l][m];
            }

            do // to check convergence
            {   
            for(l=0; l<nt; l++)
                for(m=0; m<np; m++)
                for(j=1; j<nj-1; j++)
                    for(i=1; i<ni-1; i++)
                    {
                    // S[i][j][l][m]=0;// to fullfill convergence criteria frm chai and patankar
                    Iold[i][j][l][m] = Int_new[i][j][l][m];
                    }

            for(l=0; l<nt; l++) // InScattering Integral Calculation
                for(m=0; m<np; m++)
                for(j=1; j<nj-1; j++)
                    for(i=1; i<ni-1; i++)
                    {
                    S[i][j][l][m] = 0;
                    for(l1=0; l1<nt; l1++)
                    {
                        for(m1=0; m1<np; m1++)
                        {
                        if(l1 == l && m1 == m)
                        {
                            B[i][j][l][m] = (*(D_omega + np*l1 + m1)) * (s[d] / (4*PI));
                            // B[i][j][l][m]=0;
                        }
                        else
                        {
                            // S[i][j][l][m]=0;
                            S[i][j][l][m] += (*(D_omega + np*l1 + m1)) * (s[d] / (4*PI)) * Int_new[i][j][l1][m1];
                        }
                        }
                    }
                    }

            for(l=0; l<nt; l++)
                for(m=0; m<np/4; m++)
                for(j=1; j<nj-1; j++)
                    for(i=1; i<ni-1; i++)
                    {
                    Nume = dy * (fabs(*(D_x + np*l + m))) * (Int_new[i-1][j][l][m]) + dx * (fabs(*(D_y + np*l + m))) * (Int_new[i][j-1][l][m]);
                    Source = (k[d] * sigma * pow(Tg,4) * w[d] / PI + S[i][j][l][m]) * dx * dy * (*(D_omega + np*l + m));
                    Deno = dy * (fabs(*(D_x + np*l + m))) + dx * (fabs(*(D_y + np*l + m))) + (k[d] + s[d] - B[i][j][l][m]) * dx * dy * (*(D_omega + np*l + m));
                    // Iold[i][j][l][m]=Int_new[i][j][l][m];
                    // count+=1;
                    Int_new[i][j][l][m] = (Nume + Source) / Deno;
                    // rms=fabs(Int_new[i][j][l][m]-Iold[i][j][l][m])/Int_new[i][j][l][m];
                    }

            for(l=0; l<nt; l++)
                for(m=np/4; m<np/2; m++)
                for(j=1; j<nj-1; j++)
                    for(i=ni-2; i>=1; i--)
                    {
                    Nume = dy * (fabs(*(D_x + np*l + m))) * (Int_new[i+1][j][l][m]) + dx * (fabs(*(D_y + np*l + m))) * (Int_new[i][j-1][l][m]);
                    Source = (k[d] * sigma * pow(Tg,4) * w[d] / PI + S[i][j][l][m]) * dx * dy * (*(D_omega + np*l + m));
                    Deno = dy * (fabs(*(D_x + np*l + m))) + dx * (fabs(*(D_y + np*l + m))) + (k[d] + s[d] - B[i][j][l][m]) * dx * dy * (*(D_omega + np*l + m));
                    // Iold[i][j][l][m]=Int_new[i][j][l][m];
                    Int_new[i][j][l][m] = (Nume + Source) / Deno;
                    // rms=fabs(Int_new[i][j][l][m]-Iold[i][j][l][m])/Int_new[i][j][l][m]; // printf("RMS=%lf\n",rms);
                    }

            for(l=0; l<nt; l++)
                for(m=3*np/4; m<np; m++)
                for(j=nj-2; j>=1; j--)
                    for(i=1; i<ni-1; i++)
                    {
                    Nume = dy * (fabs(*(D_x + np*l + m))) * (Int_new[i-1][j][l][m]) + dx * (fabs(*(D_y + np*l + m))) * (Int_new[i][j+1][l][m]);
                    Source = (k[d] * sigma * pow(Tg,4) * w[d] / PI + S[i][j][l][m]) * dx * dy * (*(D_omega + np*l + m));
                    Deno = dy * (fabs(*(D_x + np*l + m))) + dx * (fabs(*(D_y + np*l + m))) + (k[d] + s[d] - B[i][j][l][m]) * dx * dy * (*(D_omega + np*l + m));
                    // Iold[i][j][l][m]=Int_new[i][j][l][m];
                    Int_new[i][j][l][m] = (Nume + Source) / Deno;
                    // rms=fabs(Int_new[i][j][l][m]-Iold[i][j][l][m])/Int_new[i][j][l][m];
                    }

            for(l=0; l<nt; l++)
                for(m=np/2; m<3*np/4; m++)
                for(j=nj-2; j>=1; j--)
                    for(i=ni-2; i>=1; i--)
                    {
                    Nume = dy * (fabs(*(D_x + np*l + m))) * (Int_new[i+1][j][l][m]) + dx * (fabs(*(D_y + np*l + m))) * (Int_new[i][j+1][l][m]);
                    Source = (k[d] * sigma * pow(Tg,4) * w[d] / PI + S[i][j][l][m]) * dx * dy * (*(D_omega + np*l + m));
                    Deno = dy * (fabs(*(D_x + np*l + m))) + dx * (fabs(*(D_y + np*l + m))) + (k[d] + s[d] - B[i][j][l][m]) * dx * dy * (*(D_omega + np*l + m));
                    // Iold[i][j][l][m]=Int_new[i][j][l][m];
                    Int_new[i][j][l][m] = (Nume + Source) / Deno;
                    // rms=fabs(Int_new[i][j][l][m]-Iold[i][j][l][m])/fabs(Int_new[i][j][l][m]);
                    }

            for(l=0; l<nt; l++)
                for(m=0; m<np; m++)
                for(j=1; j<nj-1; j++)
                    for(i=1; i<ni-1; i++)
                    {
                    // S[i][j][l][m]=(1-urf)*S_old[i][j][l][m]+urf*S[i][j][l][m];
                    sum1 += Int_new[i][j][l][m];
                    sum2 += Iold[i][j][l][m];
                    }
            rms1 = fabs(sum1 - sum2) / fabs(sum1 + small);

            printf("%d\t%lf      ", d, rms1);
            time_t current_time = time(NULL);
            printf("Time elapsed: %ld seconds\n", current_time - begin);

            } while(rms1 > 1.0e-4);
            // sum1=0;
            // sum2=0;
            // printf("d=%d\n",d);
            for(j=1; j<(nj-1); j++)
            for(i=1; i<(ni-1); i++)
                for(l=0; l<nt; l++)
                for(m=0; m<np; m++)
                {
                    Int_new1[i][j][l][m] += Int_new[i][j][l][m]; // I=Ig*dg
                    as[i][j][l][m] += k[d] * Int_new[i][j][l][m];
                }
            aa += 4 * (k[d] * sigma * pow(Tg,4) * (w[d]));
            q_w1 = 0.0;
            q_w6 = 0.0;
            q_w7 = 0.0;
            q_w8 = 0.0;

            for(l=0; l<nt; l++)
            for(m=np/2; m<np; m++)
            {   
                q_w1 = q_w1 + (Int_new[(ni-2)/2][1][l][m]) * fabs(*(D_y + np*l + m));
            }

            for(l=0; l<nt; l++)
            for(m=0; m<np/2; m++)
            {   
                q_w8 = q_w8 + (Int_new[(ni-2)/2][nj-2][l][m]) * fabs(*(D_y + np*l + m)); // printf("Int=%lf\n",Int_new1[i][1][l][m]);
            }

            for(l=0; l<nt; l++)
            for(m=np/4; m<3*np/4; m++)
            {
                q_w6 = q_w6 + (Int_new[1][(ni-2)/2][l][m]) * fabs(*(D_x + np*l + m)); // printf("Int=%lf\n",Int_new1[i][1][l][m]);
            }

            for(l=0; l<nt; l++)
            for(m=0; m<np/4; m++)
                q_w7 = q_w7 + (Int_new[ni-2][(ni-2)/2][l][m]) * fabs(*(D_x + np*l + m)); // -Int_new[1][j][l][m];
            for(m=3*np/4; m<np; m++)
                q_w7 = q_w7 + (Int_new[ni-2][(ni-2)/2][l][m]) * fabs(*(D_x + np*l + m));

            coun = coun + 0.01;
            fprintf(fp6, "%lf\t%lf\n", coun, q_w6);
            fprintf(fp7, "%lf\t%lf\n", coun, q_w7);
            fprintf(fp8, "%lf\t%lf\n", coun, q_w8);
            fprintf(fp3, "%lf\t%lf\n", coun, q_w1);
        } // End of absorption coefficient loop
        sleep(3);
        time_t end = time(NULL);
        printf("The time taken by while loop is %ld\n", end - begin);
//////////////////////////////////////////////////////////////////////Heat Flux///////////////////////////////////////////////////////////////////////////////////	
    
// bottom;
for(i=1;i<(ni-1);i++)
{ 
    q_w=0.0;
    for(l=0;l<nt;l++)
        for(m=np/2;m<np;m++)
        {	
            q_w=q_w+(Int_new1[i][1][l][m])*fabs(*(D_y+np*l+m)); //printf("Int=%lf\n",Int_new1[i][1][l][m]);
        }				
    q3+=q_w;
    flux_bot[i]=(sigma*pow(Tb,4)-q_w)/(sigma*pow(Tt,4));
}
             
             
/////////////top
for(i=1;i<(ni-1);i++)
{ 
    q_w=0.0;
    for(l=0;l<nt;l++)
        for(m=0;m<np/2;m++)
        {	
            q_w=q_w+(Int_new1[i][nj-2][l][m])*fabs(*(D_y+np*l+m)); //printf("Int=%lf\n",Int_new1[i][1][l][m]);
        }				
    flux_bot_t[i]=(sigma*pow(Tb,4)-q_w)/(sigma*pow(1000,4));
}
             
             
///////left
for(j=1;j<nj-1;j++)							
{	
    q_w=0.0;				 
    for(l=0;l<nt;l++)
        for(m=np/4;m<3*np/4;m++)
        {
            q_w=q_w+(Int_new1[1][j][l][m])*fabs(*(D_x+np*l+m)); //printf("Int=%lf\n",Int_new1[i][1][l][m]);
        }				
    flux_bot_l[j]=(sigma*pow(Tb,4)-q_w)/(sigma*pow(1000,4));
}
                
                
///////right
for(j=1;j<nj-1;j++)							
{	
    q_w=0.0;				 
    for(l=0;l<nt;l++)
        for(m=0;m<np/4;m++)
            q_w=q_w+(Int_new1[ni-2][j][l][m])*fabs(*(D_x+np*l+m)); //-Int_new[1][j][l][m];
        for(m=3*np/4;m<np;m++)
            q_w=q_w+(Int_new1[ni-2][j][l][m])*fabs(*(D_x+np*l+m)); 
    flux_bot_r[j]=(sigma*pow(Tb,4)-q_w)/(sigma*pow(1000,4));
}

/////////////////////////////////////////////////////////////Contour//////////////////////////////////////////////////////////////////////////////////////////////

for(i=1;i<ni-1;i++)
{ 
    q_w=0.0,dq=0.0;
    for(l=0;l<nt;l++)
        for(m=0;m<np;m++)
        {	
            dq+=as[i][10][l][m]*(*(D_omega+np*l+m));									
        }				
    div_q[i]=-aa+dq;
}

for(j=1;j<(nj-1);j++)
    for(i=1;i<ni-1;i++)
    { 
        q_w=0.0,dq=0.0;
        for(l=0;l<nt;l++)
            for(m=0;m<np;m++)
            {	
                dq+=as[i][j][l][m]*(*(D_omega+np*l+m));									
            }				
        q2+=aa-dq;
    }
    
printf("Div is %lf\t and 4*flux is %lf\nthe error is %lf\n b-ve 1e4\n",q2*dx*dy,4*q3*dx,(fabs(q2*dx*dy)-fabs(4*q3*dx))*100/(fabs(q2*dx*dy)));
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
    
for(i=1;i<(ni-1);i++)
{
    fprintf(fp1,"%e\t%e\t%e\t%e\t%e\n",(i-0.5)*dx,flux_bot[i],flux_bot_t[i],flux_bot_l[i],flux_bot_r[i]);			
    //fprintf(fp1,"%e\t%e\n",(i-0.5)*dx,flux_bot[i]);
    fprintf(fp2,"%e\t%e\n",(i-0.5)*dx,div_q[i]);
}
          
fclose(fp1);
fclose(fp2);
return 0;
}
