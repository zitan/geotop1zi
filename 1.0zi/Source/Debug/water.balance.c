    
#include "constant.h"
#include "keywords_file.h"
#include "struct.geotop.h"
#include "pedo.funct.h"
#include "sparse_matrix.h"
#include "util_math.h"
#include "water.balance.h"

extern T_INIT *UV;
extern STRINGBIN *files;
extern long Nl, Nr, Nc;
extern double NoV;
extern DOUBLEVECTOR *outdata_basin;

#define tol_max_GC 1.E+5
#define tol_min_GC 1.E-13
#define max_res_adm 1.E-2
#define min_lambda 1.E-6
#define max_Courant_sup 0.25
#define max_Courant_channel 0.25
#define iM 1
#define ni 1.E-4

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void water_balance(ALLDATA *adt){

	DOUBLETENSOR *P;
	
	short out, sux;
	long r, c, l, iter, sy,i;
	double n, Dt0, Dt, tb, te, te0, Loss, Vout=0., Eout=0., OCout=0.,Vsubout=0.;	

	FILE *f;		

	P=new_doubletensor0(Nl, Nr, Nc);

	te0=0.0;
	
	do{

		Dt0=adt->P->Dt;
		if(te0+Dt0 > adt->P->Dt) Dt0=adt->P->Dt-te0;
		te0 += Dt0;	
		n=1.;
		te=0;
	
		do{
	
			tb=te;				
		
			do{
		
				Dt=Dt0/n;

				f=fopen(files->co[ferr]+1,"a");
				fprintf(f,"\nn:%f Dt:%f Dt0:%f te0:%f tb:%f \n",n,Dt,Dt0,te0,tb);	
				fclose(f);
				
				if(Dt<adt->P->Dt) printf("Time step:%f/%f\n",Dt,te0);
				/*
				for(i=1;i<=adt->C->r->nh;i++){
                if(r>0) adt->C->Qsub->co[i] =0;
                } */
				
				sux = Richards(Dt, P, &Loss, &iter, adt);
								
				out=1;
								
				if(sux==0 && Dt>adt->P->DtminWb){
				
					n*=adt->P->nredDtWb;
					out=0;
					
					f=fopen(files->co[ferr]+1,"a");
					fprintf(f,"Reduced time step: dt:%f Loss:%f \n\n",Dt0/n,Loss);
					fclose(f);					
				
				}else{
					
					f=fopen(files->co[ferr]+1,"a");
					fprintf(f,"Time step ok: dt:%f Loss:%f \n\n",Dt0/n,Loss);
					fclose(f);	
				
				}
									
			}while( out==0 ); 
			
			if(sux==0){
				printf("ERROR:Water balance does not converge\n");
				printf("It is not possible to continue, Please check the parameters in the block 2 of the parameter file\n");
				printf("or reduce the time step or review the initial conditions\n\n");
				printf("If you think that everything is right, please report this error to the GEOtop community\n\n");
				t_error("Code execution interrupted");
			}
			
			te=tb+Dt;                                 

			for(r=1;r<=Nr;r++){
				for(c=1;c<=Nc;c++){
					if(adt->L->LC->co[r][c]!=NoV){
						for(l=0;l<=Nl;l++){
							adt->S->P->co[l][r][c] = P->co[l][r][c];
                        }
						adt->W->hrill->co[r][c]=adt->S->P->co[0][r][c]/adt->P->wfrac;	
					}
				}
			}
			
			supflow(Dt, adt->W->hrill->co, adt->W->dh_sup->co, adt->T, adt->L, adt->W, adt->S, adt->C, adt->P, &Vout, &Vsubout, &Eout, &OCout);
					
		  for(r=1;r<=Nr;r++){
				for(c=1;c<=Nc;c++){
					if(adt->L->LC->co[r][c]!=NoV){
					      adt->S->P->co[0][r][c]=adt->W->hrill->co[r][c]*adt->P->wfrac;
				    }	
			    }
			}	
			
			outdata_basin->co[oomasserror] += fabs(Loss);
		
		}while(te<Dt0);
		
	}while(te0<adt->P->Dt);
	
	adt->C->Q_out = Vout/adt->P->Dt;
	adt->C->Qsub_out = Vsubout/adt->P->Dt;
	adt->S->outero = Eout;
	adt->S->OTEOC = OCout;
	
	for(r=1;r<=Nr;r++){
		for(c=1;c<=Nc;c++){
			if(adt->L->LC->co[r][c]!=NoV){
				for(l=1;l<=Nl;l++){
					sy = adt->S->type->co[r][c];
					
					adt->S->th->co[l][r][c] = theta_from_psi( adt->S->P->co[l][r][c], l, r, c, adt->S, PsiMin );

					//total water pressure (liq + ice)
					adt->S->Ptot->co[l][r][c] = psi_teta(adt->S->th->co[l][r][c]+adt->S->thice->co[l][r][c], 0.0, adt->S->pa->co[sy][jsat][l],
														 adt->S->pa->co[sy][jres][l], adt->S->pa->co[sy][ja][l], adt->S->pa->co[sy][jns][l], 
														 1.-1./adt->S->pa->co[sy][jns][l], PsiMin, adt->S->pa->co[sy][jss][l]);
					
					//th is the water content (always <= saturation)
					adt->S->th->co[l][r][c] = Fmin( adt->S->th->co[l][r][c] , adt->S->pa->co[sy][jsat][l]-adt->S->thice->co[l][r][c] );
				}
			}
		}
	}
	
	free_doubletensor(P);
				
}
	

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/* The code solves non-linear system F(H), where H is the total head (piezometric) [mm], F is vectorial function with components equal to the number of points,
 
 H = P + z, where P is the pressure [mm] and z the elevation over a reference height (given by the average elevation of the domain)
 
 F(H) comes out from the temporal discretization of Richards' equation PDE
 
 In this code: F(H) = I*f(H) + K(H)*H,
 
 where I is the identity matrix, f(H) = volume storage + sink, is a vectorial function (not a matrix), and K(H) is the hydraulic conductivity matrix
 
 K(H) has a number of interesting properties:
 
 1) is SPD (symmetrical positive defined)
 
 2) if kij is a component kii = sum(j!=i) kji
 
 3) the sum of all the components of the matrix is 0
 
 4) if T is a vector, the vector K T has the properties that the sum of its components is 0
 
 Here it is possible to update or not K with the new values of H, In any case the derivative of K(H) is not considered in the Jacobian. So the matrix K(H) is also
 
 used in the JACOBIAN.
 
 K is described storing only its strict lower component (strict = without diagonal) with the 3 vectors Li, Lp, Lx (in the same way as UFMPACK) */
 
 
short Richards(double Dt, DOUBLETENSOR *P, double *loss, long *iter, ALLDATA *adt){
	
	//DOUBLEVECTOR *diag, *udiag;

	double res=0.0, res0[3], res_prev[iM], res_av, res00, lambda[3], epsilon, mu,kch;

	long i, l, r, c, m, sy,ch,cont, cont2, iter_tot=0, n=(Nl+1)*adt->P->total_pixel;	
	static long cum_iter;	
	short out, out2;	
	int sux;
		
	FILE *f;
			
	if(adt->I->time == 0.0) cum_iter = 0;
		
	//diag=new_doublevector(n);
	//udiag=new_doublevector(n-1);
	
	*iter = 0;	//conjugated gradient iteration number
	
	
	for(i=1;i<=n;i++){
		
		l = adt->T->lrc_cont->co[i][1];
		r = adt->T->lrc_cont->co[i][2];
		c = adt->T->lrc_cont->co[i][3];
		
			
		
		/*
		Layer 0 is water on the surface. This allows a more robust description of infiltration processes.
		However, surface water flow is described separately in the subroutine supflow
		
		Layers > 0 , represent subsurface
		*/
		
		adt->W->H0->co[i] = adt->S->P->co[l][r][c] + adt->T->Z->co[l][r][c];		
		
		adt->W->H1->co[i] = adt->W->H0->co[i];
				
		if(adt->W->H1->co[i] != adt->W->H1->co[i]) printf("no value in r:%ld c:%ld l:%ld\n",r,c,l);
		
	}
	
	sux = find_matrix_K(adt->W->Lx, adt, adt->W->H1, Dt);
	
	
	find_f(adt->W->f, adt, adt->W->H1, Dt);
	product_matrix_using_lower_part_by_vector_plus_vector(-1., adt->W->B, adt->W->f, adt->W->H1, adt->T->Li, adt->T->Lp, adt->W->Lx);	
	res = norm_2(adt->W->B, n);
	
	res00 = res; //initial norm of the residual
	epsilon = adt->P->TolVWb + adt->P->RelTolVWb * Fmin( res00 , sqrt((double)n) );

	if(res!=res) printf("res in no value\n");

	cont=0;

	do{
		
		cont++;
		
		for(i=1;i<=n;i++){	
			adt->W->H0->co[i] = adt->W->H1->co[i];
		}

		initialize_doublevector(adt->W->dH, 0.);
		
		//mu is the forcing term of Newton Raphson, defined such that ||F( x + d )|| < mu * ||F(x)|| see references
		if (cont==1) {
			mu = adt->P->TolCG;
		}else{
			mu *= Fmin( 1.0 , res/res0[0] );
			if(mu < 0.5*epsilon/res) mu = 0.5*epsilon/res;
		}
		
		//CALCOLATE AND STORE JACOBIAN AS SPARSE MATRIX
		sux = find_dfdH(adt->W->df, adt, adt->W->H1, Dt);	//it calcolates only df/dH, J = I*df/dH + K, K is calculated above
		
		//NEGLECTING THE LATERAL FLUXES (neglecting the extra-tridiagonal part)
		/*get_diag_lower_matrix(diag, udiag, adt->T->Ai, adt->T->Ap, adt->W->Ax);
		 tridiag(0, 0, 0, n, udiag, diag, udiag, adt->W->B, adt->W->dH);*/
		
		*iter = BiCGSTAB_strict_lower_matrix_plus_identity_by_vector(mu, tol_min_GC, tol_max_GC, adt->W->dH, adt->W->B, adt->W->df, adt->T->Li, adt->T->Lp, adt->W->Lx);
		iter_tot += (*iter);//The number of the cumulated GC iterations is stored
		
						
		//non-monotonic line search (it is monotonic if M==1)
		for(m=Fminlong(cont,iM);m>1;m--){
			res_prev[m-1]=res_prev[m-2];
		}
		res_prev[0]=res;
		
		res_av=0.0;
		for(m=1;m<=Fminlong(cont,iM);m++){
			res_av=Fmax(res_prev[m-1],res_av);
		}				
		cont2 = 0.0;
		res0[0] = res;
		
		do{
			
			cont2++;
			
			/*The damping factor (or path length) lambda, defined as H1 = H0 + lambda*dH is chosen by minimizing the merit function :
			 0.5*norm(F_water(H0+lambda*dH)) interpolated with a quadratic polynome.
			 This method could not always be suited, because it has the disadvantage that it can make the solution stagnate to a point.
			 A relatively low number of MaxiterCorr (around 3-5) can prevent stagnation*/

			if(cont2 == 1){
				lambda[0] = 1.0;
				
			}else if(cont2 == 2){
				lambda[1] = lambda[0];
				res0[1] = res;
				lambda[0] = thmax;
				
			}else{
				lambda[2] = lambda[1];
				res0[2] = res0[1];
				lambda[1] = lambda[0];
				res0[1] = res;
				lambda[0] = minimize_merit_function(res0[0], lambda[1], res0[1], lambda[2], res0[2]);
				
			}
			
			for(i=1;i<=n;i++){
				
				adt->W->H1->co[i] = adt->W->H0->co[i] + lambda[0] * adt->W->dH->co[i];
				
				if(adt->W->H1->co[i] != adt->W->H1->co[i]) {
					printf("no value psi Richards3D l:%ld r:%ld c:%ld\n",l,r,c);
					stop_execution();
				}
				
			}
			
			if(adt->P->UpdateK == 1) sux = find_matrix_K(adt->W->Lx, adt, adt->W->H1, Dt);		
			//sux = find_matrix_K(adt->W->Lx, adt, adt->W->H1, Dt);
									
			find_f(adt->W->f, adt, adt->W->H1, Dt);
			product_matrix_using_lower_part_by_vector_plus_vector(-1., adt->W->B, adt->W->f, adt->W->H1, adt->T->Li, adt->T->Lp, adt->W->Lx);	
			res = norm_2(adt->W->B, n);		
			
			out2=0;
			
			if(res <= (1.0 - ni*lambda[0]*(1.-mu))*res_av) out2=1;
			if(lambda[0] <= min_lambda) out2=1;
			if(cont==1) out2=1;
					
		}while(out2==0);	
		
		out=0;
		//The condition to go out from the Newton cycle is made on the norm of the residual < Absolute tolerance + Relative tolerance * res00
		if( res <= Fmin( epsilon , max_res_adm ) ) out=1;
		//Max iteration number
		if( cont >= adt->P->MaxiterTol ) out=1;	
		
		f=fopen(files->co[ferr]+1,"a");
		fprintf(f,"iter:%ld/%ld cont:%ld cont2:%ld res:%e res_av:%e epsilon:%e mu:%e lambda:%e\n",*iter,iter_tot,cont,cont2,res,res_av,epsilon,mu,lambda[0]);
		fclose(f);

	}while(out==0);
					
	cum_iter += iter_tot;
	
	//it can be shown that massloss per unit pixel [mm] is the linear norm of B * Dt / total_pixel
	*loss = sum(adt->W->B, n)*Dt/adt->P->total_pixel;
	f=fopen(files->co[ferr]+1,"a");
	fprintf(f,"cumiter:%ld massloss:%e\n",cum_iter,*loss);
	fclose(f);
	
	//assign updated state variables
	for(i=1;i<=n;i++){
		l=adt->T->lrc_cont->co[i][1];
		r=adt->T->lrc_cont->co[i][2];
		c=adt->T->lrc_cont->co[i][3];
		sy=adt->S->type->co[r][c];
		ch=adt->C->ch->co[r][c];
		
		P->co[l][r][c] = adt->W->H1->co[i] - adt->T->Z->co[l][r][c];
		
		if(ch!=0 && l>0){
		if( P->co[l][r][c] > adt->C->h_sup->co[ch]-adt->P->depr_channel){
		kch=k_from_psi(jKh, P->co[l][r][c], l, r, c, adt->S, adt->P->imp);
		adt->C->Qsub->co[ch]+=kch*(P->co[l][r][c] - adt->C->h_sup->co[ch]+adt->P->depr_channel)
		                     *1E-9*UV->U->co[1]*adt->S->pa->co[sy][jdz][l]/UV->U->co[2];}
		}else{adt->C->Qsub->co[ch]=0;}
		
	}	
	
	//calculate Qsub
	/*
	for(i=1;i<=adt->C->r->nh;i++){
		r = adt->C->r->co[i];
		c = adt->C->c->co[i];
		//there are channel pixels, if there are not channels r vector has 1 component and is equal to 0
		if(r>0) adt->C->Qsub->co[i] = adt->P->Kch_b * (adt->P->w_dx*UV->U->co[1]*adt->C->length->co[i]) 
		* 1.E-3*Fmax( P->co[0][r][c] - adt->C->h_sup->co[i] + adt->P->depr_channel, 0.0 );//m3/s
	
	}
	*/
		
	//end subroutine
	if( res <= epsilon ){
		return 1;//converges
	}else{
		return 0;//does not converge
	}
		
	
}


/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double cm_h(double cm0, double h, double h_thres){
	
	double cm;
	
	if(h > h_thres){
		cm = cm0;
	}else if(h > 0){
		cm = ( 0.0 * (h_thres - h) + cm0*(h) ) / h_thres;
	}else{
		cm = 0.0;
	}
	
	return(cm);
	
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void supflow(double Dt, double **h, double **dh, TOPO *top, LAND *land, WATER *wat, SOIL *sl, CHANNEL *cnet, PAR *par, double *Vout, double *Vsubout, double *Eout, double *OCout)

{
	long r,c,R,C,ch;                                    
	static short r_DD[11]={0,0,-1,-1,-1,0,1,1,1,-9999,0}; // differential of number-pixel for rows and
	static short c_DD[11]={0,1,1,0,-1,-1,-1,0,1,-9999,0}; // columns, depending on Drainage Directions
	double dx,dy;                                         // the two dimensions of a pixel
	double Ks;											  // the Strickler's coefficent calculated with a standard deviation
	double b[10],lth[10];                                 // area perpendicular to the superficial flow divided by hsup
	double i;											  // hydraulic gradient
	double tb,te,dt,qch;
	short lu;
	DOUBLEMATRIX *cl;
    DOUBLEMATRIX *dl;	
	DOUBLEMATRIX *TC;
	DOUBLEMATRIX *Conc;
	DOUBLEMATRIX *D;
	DOUBLEMATRIX *q;
	DOUBLEMATRIX *v;


	if(par->point_sim!=1){	//distributed simulations
		
		dx=UV->U->co[1];                                    
		dy=UV->U->co[2]; 
		
        lth[1]=0.0;  
		lth[2]=dy;             
		lth[3]=sqrt(dx*dx+dy*dy); 
		lth[4]=dx;            
		lth[5]=sqrt(dx*dx+dy*dy); 
		lth[6]=dy;   
		lth[7]=sqrt(dx*dx+dy*dy); 
		lth[8]=dx;             
		lth[9]=sqrt(dx*dx+dy*dy);	
	
	
		b[1]=0.0;  
		b[2]=dy*par->wfrac;             
		b[3]=sqrt(dx*dx+dy*dy)*par->wfrac; 
		b[4]=dx*par->wfrac;            
		b[5]=sqrt(dx*dx+dy*dy)*par->wfrac; 
		b[6]=dy*par->wfrac;   
		b[7]=sqrt(dx*dx+dy*dy)*par->wfrac; 
		b[8]=dx*par->wfrac;             
		b[9]=sqrt(dx*dx+dy*dy)*par->wfrac; 
		
	
	/*
		b[1]=0.0;  
		b[2]=dy;             
		b[3]=sqrt(dx*dx+dy*dy); 
		b[4]=dx;            
		b[5]=sqrt(dx*dx+dy*dy); 
		b[6]=dy;   
		b[7]=sqrt(dx*dx+dy*dy); 
		b[8]=dx;             
		b[9]=sqrt(dx*dx+dy*dy); 	
	*/	
		
	cl=new_doublematrix(Nr,Nc); 
    initialize_doublematrix(cl,0);	
	dl=new_doublematrix(Nr,Nc);
	initialize_doublematrix(dl,0);
	TC=new_doublematrix(Nr,Nc);
	initialize_doublematrix(TC,0);
	v=new_doublematrix(Nr,Nc);
	initialize_doublematrix(v,0);
	q=new_doublematrix(Nr,Nc);
	initialize_doublematrix(q,0);
	D=new_doublematrix(Nr,Nc);
	initialize_doublematrix(D,0);
	Conc=new_doublematrix(Nr,Nc);
	initialize_doublematrix(Conc,0);
 			
	  EroParaInit(cl, dl, land, sl);    //  initialize the erosion calculation parameters
	
		te=0.0;
		
		do{
			
			tb=te;
			dt=Dt;
			
			
			find_dt_max(max_Courant_sup, h, land, top, cnet, par, &dt);
			
			te=tb+dt;
			if(te>Dt){
				te=Dt;
				dt=te-tb;
			}			
						
			for(r=1;r<=Nr;r++){
				for(c=1;c<=Nc;c++){
					if(land->LC->co[r][c]!=NoV){
							//wat->Min->co[r][c]=0;					
						if(h[r][c] > 0){
							lu=(short)land->LC->co[r][c];
							Ks=cm_h(land->ty->co[lu][jcm], h[r][c], par->thres_hsup);
							
							R=r+r_DD[top->DD->co[r][c]];
							C=c+c_DD[top->DD->co[r][c]];
							if(land->LC->co[R][C]!=NoV){
						    
								i=0.001*( h[r][c] - h[R][C])/lth[top->DD->co[r][c]+1] + top->i_DD->co[r][c];				
								if(i<0) i=0.0;							
								q->co[r][c] = b[top->DD->co[r][c]+1]*Ks*pow(h[r][c]/1000.0,1.0+par->gamma_m)*sqrt(i)*1000.0/dx/dy;  //mm/s
								
							}else{//pixel draining out of the basin
								
								i=top->i_DD->co[r][c];				
								if(i<0) i=0.0;							
								q->co[r][c] = b[top->DD->co[r][c]+1]*Ks*pow(h[r][c]/1000.0,1.0+par->gamma_m)*sqrt(i)*1000.0/dx/dy;
								
								*Vout = *Vout + q->co[r][c]*dx*dy*dt/1000.;	//m3
								//*Eout+=q->co[r][c]*dx*dy*dt/1000*Conc->co[r][c];
							}
								
						        v->co[r][c]=q->co[r][c]*dx*dy/h[r][c]/b[top->DD->co[r][c]+1];   //Calculate mean flow velocity  m/s
				    
						}else{
							q->co[r][c] = 0.0;
							v->co[r][c]=0;
						}
						
						dh[r][c] = Fmin( q->co[r][c]*dx*dy*dt/b[top->DD->co[r][c]+1]/lth[top->DD->co[r][c]+1], Fmax( h[r][c] , 0.0 ) );
										
					}
				}
			}
				
	    TranCap(TC,v, h, q, cl, dl, top, land, sl);        // calculate transport capacity
        D_Cal(h, D, TC, q, Conc,top, land, wat, sl, par);  // Calculate the flux btw soil and water
		
			for(r=1;r<=Nr;r++){
				for(c=1;c<=Nc;c++){
					if(land->LC->co[r][c]!=NoV){
				
						h[r][c] -= dh[r][c];
					   /* if(h[r][c]<0){
						printf("haha\n");
						}
                         */						
						R=r+r_DD[top->DD->co[r][c]];
						C=c+c_DD[top->DD->co[r][c]];
						
						if(land->LC->co[R][C]!=NoV){
			
							if(h[R][C]>0){
								h[R][C] += dh[r][c];
							}else{
								if( dh[r][c] > 0) h[R][C] = dh[r][c];
							}
						
						}				
					}
				}
			}
			
          Erosion(dt, h, Conc, q, D, top, land, wat, sl, par, cnet, OCout);   // calculate erosion based on mass balance of sediments
 
			//Superficial flow to the channels
			for(ch=1;ch<=cnet->r->nh;ch++){
				r = cnet->r->co[ch];
				c = cnet->c->co[ch];
				
				if(r>0){
				
					if(cnet->h_sup->co[ch] <= par->depr_channel){	//free flow
									
						if(h[r][c] > 0){
							qch = Cd*(2./3.)*sqrt(2.*g*1.E-3*h[r][c])*((2.*cnet->length->co[ch])/(dx*dy))*h[r][c];//mm/s
							if((qch*dt*dx*dy/b[top->DD->co[r][c]+1]/lth[top->DD->co[r][c]+1])> h[r][c]) qch = h[r][c]*b[top->DD->co[r][c]+1]*lth[top->DD->co[r][c]+1]/dx/dy/dt;
						}else{
							qch = 0.0;
						}
						
					}else if(h[r][c] > cnet->h_sup->co[ch] - par->depr_channel){//submerged flow towards channel
					
						qch = Cd*(2./3.)*sqrt(2.*g*1.E-3*(h[r][c]-(cnet->h_sup->co[ch] - par->depr_channel)))*((2.*cnet->length->co[ch])/(dx*dy))*h[r][c];//mm/s
						if((qch*dt*dx*dy/b[top->DD->co[r][c]+1]/lth[top->DD->co[r][c]+1])> h[r][c]) qch = h[r][c]*b[top->DD->co[r][c]+1]*lth[top->DD->co[r][c]+1]/dx/dy/dt;			
					
					}else{	//submerged flow towards land
						
						qch = 0.0;
						
					}
				
					h[r][c] -= qch*dx*dy*dt/b[top->DD->co[r][c]+1]/lth[top->DD->co[r][c]+1];
										
					cnet->Qsup->co[ch] = 1.E-3*qch*dx*dy;	//m3/s
				//	cnet->Inmass->co[ch] = Conc->co[r][c]*cnet->Qsup->co[ch];    //sediments inflow into the channel  kg/s
				
				}
				
			}
		 
	
			channel_flow(dt, cnet->h_sup, cnet->dh_sup, top, cnet, par, Vout, Vsubout, Eout);
	 
		}while(te<Dt);
	
	    free_doublematrix(cl);
	    free_doublematrix(dl);
	    free_doublematrix(TC);
	    free_doublematrix(Conc);
	    free_doublematrix(v);
		free_doublematrix(q);
		free_doublematrix(D);

		
	}else{	//point simulation  
		
		for(r=1;r<=Nr;r++){
			for(c=1;c<=Nc;c++){
				if(land->LC->co[r][c]!=NoV){
					h[r][c]=Fmin(h[r][c], 0.0);
				}
			}
		}
	}
}	


/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void find_dt_max(double Courant, double **h, LAND *land, TOPO *top, CHANNEL *cnet, PAR *par, double *dt){
	
	double b[10],lth[10];    //area perpendicular to the superficial flow divided by hsup
	double i, q, dx, dy, Ks;											  
	short lu;
	long r, c, R, C, ch;
	
	short r_DD[11]={0,0,-1,-1,-1,0,1,1,1,-9999,0}; // differential of number-pixel for rows and
	short c_DD[11]={0,1,1,0,-1,-1,-1,0,1,-9999,0}; // columns, depending on Drainage Directions
	
	dx=UV->U->co[1];                                    
	dy=UV->U->co[2];   
	
		lth[1]=0.0;
	    lth[2]=dx;
	    lth[3]=sqrt(dx*dx+dy*dy);
	    lth[4]=dy;
	    lth[5]=sqrt(dx*dx+dy*dy);
	    lth[6]=dx;
	    lth[7]=sqrt(dx*dx+dy*dy);
	    lth[8]=dy;
	    lth[9]=sqrt(dx*dx+dy*dy);	
		
		
	    b[1]=0.0;  
		b[2]=dy*par->wfrac;             
		b[3]=sqrt(dx*dx+dy*dy)*par->wfrac; 
		b[4]=dx*par->wfrac;            
		b[5]=sqrt(dx*dx+dy*dy)*par->wfrac; 
		b[6]=dy*par->wfrac;   
		b[7]=sqrt(dx*dx+dy*dy)*par->wfrac; 
		b[8]=dx*par->wfrac;             
		b[9]=sqrt(dx*dx+dy*dy)*par->wfrac; 	
		
		/*
		b[1]=0.0;  
		b[2]=dy;             
		b[3]=sqrt(dx*dx+dy*dy); 
		b[4]=dx;            
		b[5]=sqrt(dx*dx+dy*dy); 
		b[6]=dy;   
		b[7]=sqrt(dx*dx+dy*dy); 
		b[8]=dx;             
		b[9]=sqrt(dx*dx+dy*dy); 	
		*/

	for(r=1;r<=Nr;r++){
		for(c=1;c<=Nc;c++){
			if(h[r][c]>0 && land->LC->co[r][c]!=NoV){
				
				lu=(short)land->LC->co[r][c];
				Ks=cm_h(land->ty->co[lu][jcm], h[r][c], par->thres_hsup);
				
				R=r+r_DD[top->DD->co[r][c]];
				C=c+c_DD[top->DD->co[r][c]];
				
				if(land->LC->co[R][C]!=NoV){
					i=0.001*( h[r][c] - h[R][C] )/lth[top->DD->co[r][c]+1] + top->i_DD->co[r][c];				
				}else{
					i=top->i_DD->co[r][c];
				}
				if(i<0) i=0;  //correct false slopes
				
				q=b[top->DD->co[r][c]+1]*Ks*pow(h[r][c]/1000.0,1.0+par->gamma_m)*sqrt(i)*1000.0/dx/dy;	//mm/s
				
				if(q>0){
					if(Courant*h[r][c]*lth[top->DD->co[r][c]+1]*b[top->DD->co[r][c]+1]/q/dx/dy<(*dt)) *dt=Courant*h[r][c]
					*lth[top->DD->co[r][c]+1]*b[top->DD->co[r][c]+1]/q/dx/dy; 
				}
				
				if(top->pixel_type->co[r][c]>=10){	//channel pixel
					
					ch = cnet->ch->co[r][c];
					q = 0.;
					
					if(cnet->h_sup->co[ch] <= par->depr_channel){	//free flow
						
						q = Cd*(2./3.)*sqrt(2.*g*1.E-3*h[r][c])*((2.*cnet->length->co[ch])/(dx*dy))*h[r][c];//mm/s
						
					}else if(h[r][c] > cnet->h_sup->co[ch] - par->depr_channel){//submerged flow towards channel
						
						q = Cd*(2./3.)*sqrt(2.*g*1.E-3*(h[r][c]-(cnet->h_sup->co[ch] - par->depr_channel)))*((2.*cnet->length->co[ch])/(dx*dy))*h[r][c];//mm/s
					
					}
						
					if(q>0){
						if(Courant*h[r][c]*lth[top->DD->co[r][c]+1]*b[top->DD->co[r][c]+1]/q/dx/dy<(*dt)) *dt=Courant*h[r][c]
						*lth[top->DD->co[r][c]+1]*b[top->DD->co[r][c]+1]/q/dx/dy; 
					}
				}
				
					
			}
		}
	}
}
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
//cnt is the counter of Li Lp Lx (lower diagonal without diagonal)
int find_matrix_K(DOUBLEVECTOR *Lx, ALLDATA *adt, DOUBLEVECTOR *H, double Dt){
	
	long i, l, r, c, I, R, C, sy, cnt=0;
	double dz=0.0, dzn=0.0, k=0.0, kn=0.0, klim=0.0, ds;
	
	for(i=1;i<=H->nh;i++){
		
		l=adt->T->lrc_cont->co[i][1];
		r=adt->T->lrc_cont->co[i][2];
		c=adt->T->lrc_cont->co[i][3];
		sy=adt->S->type->co[r][c];
		
		
		//vertical hydraulic conductivity
		if(l>0){
			if(adt->P->UpdateK==1){
				k = k_from_psi(jKv, H->co[i] - adt->T->Z->co[l][r][c], l, r, c, adt->S, adt->P->imp );
			}else{
				k = k_from_psi(jKv, adt->S->P->co[l][r][c], l, r, c, adt->S, adt->P->imp );
			}
			dz = adt->S->pa->co[sy][jdz][l];
		}
		
				
		
		//flux from cell below
		if (l<Nl) {
			cnt++;
			I = i+1;
			
			if(l==0){	//overland flow
				
				if( H->co[i]< H->co[I] ){	//upward flux
					kn = k_from_psi(jKv,  H->co[I] - adt->T->Z->co[l+1][r][c], l+1, r, c, adt->S, adt->P->imp );
				}else{	//downward flow
					kn = k_from_psi(jKv,  psisat_from(l+1, r, c, adt->S), l+1, r, c, adt->S, adt->P->imp );
				}
				
				dzn = 0.5*adt->S->pa->co[sy][jdz][l+1];

			}else{	//subsurface flow
				
				if(adt->P->UpdateK==1){
					kn = k_from_psi(jKv,  H->co[I] - adt->T->Z->co[l+1][r][c], l+1, r, c, adt->S, adt->P->imp );
				}else{
					kn = k_from_psi(jKv,  adt->S->P->co[l+1][r][c], l+1, r, c, adt->S, adt->P->imp );
				}
				
				dzn = adt->S->pa->co[sy][jdz][l+1];
				kn = Mean(adt->P->harm_or_arit_mean, dz, dzn, k, kn);
				klim = Harmonic_Mean( dz, dzn, k_from_psi(jKv,  psisat_from(l, r, c, adt->S), l, r, c, adt->S, adt->P->imp ),
									 k_from_psi(jKv,  psisat_from(l+1, r, c, adt->S), l+1, r, c, adt->S, adt->P->imp ) );
				if(kn > klim) kn=klim;
				dzn = 0.5*dz + 0.5*dzn;
								
			}
						
			Lx->co[cnt] = -kn/dzn;
			
		}		
		
		
		//lateral hydraulic conductivity
		if(l>0){
			if(adt->P->UpdateK==1){
				k = k_from_psi(jKh, H->co[i] - adt->T->Z->co[l][r][c], l, r, c, adt->S, adt->P->imp );
			}else{
				k = k_from_psi(jKh, adt->S->P->co[l][r][c], l, r, c, adt->S, adt->P->imp );
			}
		}		
		
		
		//lateral fluxes
		//4 neighbouring cells
		R = r-1;
		C = c;
		ds = 1.E3*UV->U->co[2];	//[mm]
		//1.
		if(R>=1 && R<=Nr && C>=1 && C<=Nc && l>0 && adt->P->point_sim!=1){
			if(adt->L->LC->co[R][C]!=NoV && adt->T->i_cont[l][R][C]>i){ 
				
				I=adt->T->i_cont[l][R][C];	
				
				if(adt->P->UpdateK==1){
					kn = k_from_psi(jKh, H->co[I] - adt->T->Z->co[l][R][C], l, R, C, adt->S, adt->P->imp );
				}else{
					kn = k_from_psi(jKh, adt->S->P->co[l][R][C], l, R, C, adt->S, adt->P->imp );
				}
				
				kn = Mean(adt->P->harm_or_arit_mean, ds, ds, k, kn);
				klim = Harmonic_Mean( dz, dzn, k_from_psi(jKv,  psisat_from(l, r, c, adt->S), l, r, c, adt->S, adt->P->imp ),
									 k_from_psi(jKv,  psisat_from(l, R, C, adt->S), l, R, C, adt->S, adt->P->imp ) );
				if(kn > klim) kn=klim;
				
				cnt++;
				Lx->co[cnt] = -(dz/ds)*kn/ds;
				
			}
		}
		
		R = r+1;
		C = c;
		ds = 1.E3*UV->U->co[2];	//[mm]
		//2.
		if(R>=1 && R<=Nr && C>=1 && C<=Nc && l>0 && adt->P->point_sim!=1){
			if(adt->L->LC->co[R][C]!=NoV && adt->T->i_cont[l][R][C]>i){ 
				
				I=adt->T->i_cont[l][R][C];	
				
				if(adt->P->UpdateK==1){
					kn = k_from_psi(jKh, H->co[I] - adt->T->Z->co[l][R][C], l, R, C, adt->S, adt->P->imp );
				}else{
					kn = k_from_psi(jKh, adt->S->P->co[l][R][C], l, R, C, adt->S, adt->P->imp );
				}
				
				kn = Mean(adt->P->harm_or_arit_mean, ds, ds, k, kn);
				klim = Harmonic_Mean( dz, dzn, k_from_psi(jKv,  psisat_from(l, r, c, adt->S), l, r, c, adt->S, adt->P->imp ),
									 k_from_psi(jKv,  psisat_from(l, R, C, adt->S), l, R, C, adt->S, adt->P->imp ) );
				if(kn > klim) kn=klim;
				
				cnt++;
				Lx->co[cnt] = -(dz/ds)*kn/ds;
				
				
			}
		}
		
		R = r;
		C = c-1;
		ds = 1.E3*UV->U->co[1];	//[mm]
		//3.
		if(R>=1 && R<=Nr && C>=1 && C<=Nc && l>0 && adt->P->point_sim!=1){
			if(adt->L->LC->co[R][C]!=NoV && adt->T->i_cont[l][R][C]>i){ 
				
				I=adt->T->i_cont[l][R][C];	
				
				if(adt->P->UpdateK==1){
					kn = k_from_psi(jKh, H->co[I] - adt->T->Z->co[l][R][C], l, R, C, adt->S, adt->P->imp );
				}else{
					kn = k_from_psi(jKh, adt->S->P->co[l][R][C], l, R, C, adt->S, adt->P->imp );
				}
				
				kn = Mean(adt->P->harm_or_arit_mean, ds, ds, k, kn);
				klim = Harmonic_Mean( dz, dzn, k_from_psi(jKv,  psisat_from(l, r, c, adt->S), l, r, c, adt->S, adt->P->imp ),
									 k_from_psi(jKv,  psisat_from(l, R, C, adt->S), l, R, C, adt->S, adt->P->imp ) );
				if(kn > klim) kn=klim;
				
				cnt++;
				Lx->co[cnt] = -(dz/ds)*kn/ds;
				
				
			}
		}
		
		
		R = r;
		C = c+1;
		ds = 1.E3*UV->U->co[1];	//[mm]
		//4.
		if(R>=1 && R<=Nr && C>=1 && C<=Nc && l>0 && adt->P->point_sim!=1){
			if(adt->L->LC->co[R][C]!=NoV && adt->T->i_cont[l][R][C]>i){ 
				
				I=adt->T->i_cont[l][R][C];	
				
				if(adt->P->UpdateK==1){
					kn = k_from_psi(jKh, H->co[I] - adt->T->Z->co[l][R][C], l, R, C, adt->S, adt->P->imp );
				}else{
					kn = k_from_psi(jKh, adt->S->P->co[l][R][C], l, R, C, adt->S, adt->P->imp );
				}
				
				kn = Mean(adt->P->harm_or_arit_mean, ds, ds, k, kn);
				klim = Harmonic_Mean( dz, dzn, k_from_psi(jKv,  psisat_from(l, r, c, adt->S), l, r, c, adt->S, adt->P->imp ),
									 k_from_psi(jKv,  psisat_from(l, R, C, adt->S), l, R, C, adt->S, adt->P->imp ) );
				if(kn > klim) kn=klim;
				
				cnt++;
				Lx->co[cnt] = -(dz/ds)*kn/ds;
				
			}
		}
	}
	
	return 0;
	
}


/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

int find_dfdH(DOUBLEVECTOR *df, ALLDATA *adt, DOUBLEVECTOR *H, double Dt){
	
	long i,r,l,c,sy,ch;
	double h,dz,kch;
		
	for(i=1;i<=H->nh;i++){
		
		l=adt->T->lrc_cont->co[i][1];
		r=adt->T->lrc_cont->co[i][2];
		c=adt->T->lrc_cont->co[i][3];		
		sy=adt->S->type->co[r][c];
		ch=adt->C->ch->co[r][c];
		
		//hydraulic capacity (diagonal term) = (dV/dH)/(Ah*Dt)
		if(l==0){
			h = H->co[i] - 1.E3*adt->T->Z0dp->co[r][c];
			if(h>0){
				df->co[i] = 1./Dt;
			}else{
				df->co[i] = 0.;
			}
		}else{
			dz = adt->S->pa->co[sy][jdz][l];
			df->co[i] = dtheta_dpsi_from_psi(H->co[i] - adt->T->Z->co[l][r][c], l, r, c, adt->S, PsiMin)*dz/Dt;
		}
	
		
		//add subsuperficial flow to the channels (sink depending on H)
		
		if(ch!=0 && l>0){ 
			h = H->co[i] - adt->T->Z->co[l][r][c];			
			if( h > adt->C->h_sup->co[ch]-adt->P->depr_channel) {
			kch=k_from_psi(jKh, h, l, r, c, adt->S, adt->P->imp );
			df->co[i] +=kch*1E-6*adt->S->pa->co[sy][jdz][l]/UV->U->co[1]/UV->U->co[1];
			}				 
		}
		
		
		
	/*	if(ch!=0 && l==0){
			h = H->co[i] - adt->T->Z->co[l][r][c];
			if( h > adt->C->h_sup->co[ch]-adt->P->depr_channel) df->co[i] +=  (adt->P->w_dx*adt->C->length->co[ch]/UV->U->co[1]) 
			* adt->P->Kch_b;			
		}  */
	}
	return 0;
}
		
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

int find_f(DOUBLEVECTOR *f, ALLDATA *adt, DOUBLEVECTOR *H, double Dt){

	long i, l, r, c, sy, ch;
	double dz, h, V0, V1,kch;
	
	for(i=1;i<=H->nh;i++){
		
		l=adt->T->lrc_cont->co[i][1];
		r=adt->T->lrc_cont->co[i][2];
		c=adt->T->lrc_cont->co[i][3];
		sy=adt->S->type->co[r][c];
		ch=adt->C->ch->co[r][c];
		
		/*if(l==0 && r==45 && c==55){
		printf("hehe");
		}*/				
	    //hydraulic capacity (diagonal term)
		if(l==0){
			V1 = Fmax(0.0, H->co[i] - adt->T->Z->co[l][r][c]);
			V0 = Fmax(0.0, adt->S->P->co[l][r][c]);
		}else{
			dz = adt->S->pa->co[sy][jdz][l];		
			V1 = dz * theta_from_psi( H->co[i] - adt->T->Z->co[l][r][c], l, r, c, adt->S, PsiMin );
			V0 = dz * theta_from_psi( adt->S->P->co[l][r][c], l, r, c, adt->S, PsiMin );
		}
		f->co[i] = (V1-V0)/Dt;
			
		//subsuperficial flow to the channels
	/*	if(ch!=0 && l==0){ 
			h = H->co[i] - adt->T->Z->co[l][r][c];
			if( h > adt->C->h_sup->co[ch]-adt->P->depr_channel) f->co[i] += (adt->P->w_dx*adt->C->length->co[ch]/UV->U->co[1]) * 
			adt->P->Kch_b * ( h - adt->C->h_sup->co[ch]+ adt->P->depr_channel );
		}  */
	
		if(ch!=0 && l>0){ 
			h = H->co[i] - adt->T->Z->co[l][r][c];
			if( h > adt->C->h_sup->co[ch]-adt->P->depr_channel) {
			kch=k_from_psi(jKh, h, l, r, c, adt->S, adt->P->imp );
			f->co[i] +=kch*(h - adt->C->h_sup->co[ch]+adt->P->depr_channel)
		                     *1E-6*adt->S->pa->co[sy][jdz][l]/UV->U->co[1]/UV->U->co[1];
			
			}				 
		}
		
				//subsuperficial flow to the channels
		/*		
		if(ch!=0 && l==0){ 
			h = H->co[i] - adt->T->Z->co[l][r][c];
			if( h > adt->C->h_sup->co[ch]-adt->P->depr_channel) f->co[i] += (adt->P->w_dx*adt->C->length->co[ch]/UV->U->co[1]) * adt->P->Kch_b * ( h - adt->C->h_sup->co[ch] );
		}
		
		*/
	
		//evaporation and precipitation
		if(l>0){
			f->co[i] += adt->S->ET->co[l][r][c];
		}else{
			f->co[i] -= adt->W->Pnet->co[r][c];	
		}
				
						
	}
	return 0;
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void find_dt_max_channel(double Courant, DOUBLEVECTOR *h, TOPO *top, CHANNEL *cnet, PAR *par, double *dt){

	short r_DD[11]={0,0,-1,-1,-1,0,1,1,1,-9999,0}; // differential of number-pixel for rows and
	short c_DD[11]={0,1,1,0,-1,-1,-1,0,1,-9999,0}; // columns, depending on Drainage Directions
	long r, c, ch, r_down, c_down, ch_down;		
	double b[10];    //channel length
	double dn, Ks, q, Vmax, i;

	dn = par->w_dx * UV->U->co[1];		//transversal [m]
	
	b[1]=UV->U->co[1];					//distance between channel centers
	b[2]=UV->U->co[1];             
	b[3]=UV->U->co[1]*sqrt(2.); 
	b[4]=UV->U->co[1];            
	b[5]=UV->U->co[1]*sqrt(2.); 
	b[6]=UV->U->co[1];   
	b[7]=UV->U->co[1]*sqrt(2.); 
	b[8]=UV->U->co[1];             
	b[9]=UV->U->co[1]*sqrt(2.); 
	
	for(ch=1;ch<=cnet->r->nh;ch++){
		
		if(h->co[ch] > 0){
			Ks=cm_h(par->Ks_channel, h->co[ch], par->thres_hchannel);
			
			r = cnet->r->co[ch];
			c = cnet->c->co[ch];
			
			r_down = r+r_DD[top->DD->co[r][c]];
			c_down = c+c_DD[top->DD->co[r][c]];
			ch_down = cnet->ch->co[r_down][c_down];
			
			if(ch_down==0) ch_down=ch;
						
			i = Fmax( 0.0 , 1.E-3 * ( h->co[ch] - h->co[ch_down] ) / b[top->DD->co[r][c]+1] + top->i_DD->co[r][c] );	
			
			q = dn * Ks * pow( 1.E-3 * h->co[ch] , 1.0+par->gamma_m ) * sqrt(i);	//m3/s
			Vmax = 1.E-3*h->co[ch]*dn*cnet->length->co[ch];	//m3
			
			if(q>0){
				if(Courant*Vmax/q<(*dt)) *dt=Courant*Vmax/q; 
			}
		}
    }
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void channel_flow(double Dt, DOUBLEVECTOR *h, DOUBLEVECTOR *dV, TOPO *top, CHANNEL *cnet, PAR *par, double *Vout, double *Vsubout, double *Eout)

{
	static short r_DD[11]={0,0,-1,-1,-1,0,1,1,1,-9999,0}; // differential of number-pixel for rows and
	static short c_DD[11]={0,1,1,0,-1,-1,-1,0,1,-9999,0}; // columns, depending on Drainage Directions
	long r,c,ch,r_down,c_down,ch_down;                                    
	double dn;                                         // the two dimensions of a channel pixel
	double Ks;											  // the Strickler's coefficent
	double b[10];                                         // distance between channels centers
	double i;											  // hydraulic gradient
	double tb,te,dt;
	double dVsub[cnet->r->nh];
	
	DOUBLEVECTOR *Mcin; 
	DOUBLEVECTOR *qch;
	
	Mcin=new_doublevector(cnet->r->nh); 
    initialize_doublevector(Mcin,0);
	qch=new_doublevector(cnet->r->nh); 
    initialize_doublevector(qch,0);
	
	for(ch=1;ch<=cnet->r->nh;ch++){
	dVsub[ch]=0;
	}
	
	if( par->point_sim!=1 && cnet->r->co[1]!=0 ){	//if it is not point simulation and there are channels
		
		dn = par->w_dx * UV->U->co[1];		//transversal length [m]
		
		b[1]=UV->U->co[1];					//distance between channel centers
		b[2]=UV->U->co[1];             
		b[3]=UV->U->co[1]*sqrt(2.); 
		b[4]=UV->U->co[1];            
		b[5]=UV->U->co[1]*sqrt(2.); 
		b[6]=UV->U->co[1];   
		b[7]=UV->U->co[1]*sqrt(2.); 
		b[8]=UV->U->co[1];             
		b[9]=UV->U->co[1]*sqrt(2.); 
		
		te=0.0;
		
		do{
			
			tb=te;
			dt=Dt;
			
			find_dt_max_channel(max_Courant_channel, h, top, cnet, par, &dt);
			
			te=tb+dt;
			if(te>Dt){
				te=Dt;
				dt=te-tb;
			}		
			
			for(ch=1;ch<=cnet->r->nh;ch++){
				
				if(h->co[ch] > 0){
					
					r = cnet->r->co[ch];
					c = cnet->c->co[ch];
					r_down = r+r_DD[top->DD->co[r][c]];
					c_down = c+c_DD[top->DD->co[r][c]];
					ch_down = cnet->ch->co[r_down][c_down];
					
					Ks=cm_h(par->Ks_channel, h->co[ch], par->thres_hchannel);
					
					if(ch_down==0){//outlet section
						
						qch->co[ch] = Cd*(2./3.)*sqrt(2.*g*1.E-3*h->co[ch])*(1.E-3*h->co[ch])*dn;	//[m3/s]
						
					}else{

						i = Fmax( 0.0 , 1.E-3 * ( h->co[ch] - h->co[ch_down] ) / b[top->DD->co[r][c]+1] + top->i_DD->co[r][c] );	
						qch->co[ch] = dn * Ks * pow( 1.E-3 * h->co[ch] , 1.0+par->gamma_m ) * sqrt(i);	//m3/s
						
					}
					
					if(qch->co[ch]*dt > 1.E-3*h->co[ch]*dn*cnet->length->co[ch]){
                    qch->co[ch] = 1.E-3*h->co[ch]*dn*cnet->length->co[ch]/dt;
                    }	
					
					dV->co[ch] = Fmin( qch->co[ch]*dt , 1.E-3*h->co[ch]*dn*cnet->length->co[ch] );	//m3
				   
                    /*inflow mass calculation kg*/				   
					Mcin->co[ch]=1.E-3*h->co[ch]*dn*cnet->length->co[ch]*cnet->Conc->co[ch]+cnet->Inmass->co[ch]*dt;
					
				
				}else{
					
					dV->co[ch] = 0.0;
			    	Mcin->co[ch] = 0.0;
				}
			}
			
			for(ch=1;ch<=cnet->r->nh;ch++){
			   if(h->co[ch] > 0){					
					r = cnet->r->co[ch];
					c = cnet->c->co[ch];
					r_down = r+r_DD[top->DD->co[r][c]];
					c_down = c+c_DD[top->DD->co[r][c]];
					ch_down = cnet->ch->co[r_down][c_down];
				    Mcin->co[ch_down]+=qch->co[ch]*cnet->Conc->co[ch]*dt;	
                }
            }	
		   
			for(ch=1;ch<=cnet->r->nh;ch++){
				
				r = cnet->r->co[ch];
				c = cnet->c->co[ch];
				r_down = r+r_DD[top->DD->co[r][c]];
				c_down = c+c_DD[top->DD->co[r][c]];
				ch_down = cnet->ch->co[r_down][c_down];
				if(ch_down==0) {ch_down=ch;}

				h->co[ch] -= 1.E3*dV->co[ch]/(dn*cnet->length->co[ch]);
												
				if (ch_down != ch){	//not the outlet
					h->co[ch_down] += 1.E3*dV->co[ch]/(dn*cnet->length->co[ch_down]);	//mm
				}else{	//outlet				
					*Vout = *Vout + dV->co[ch];	//m3
												    
				}
			}
			
			for(ch=1;ch<=cnet->r->nh;ch++){
				
				dVsub[ch]=cnet->Qsub->co[ch]*dt; //m3
				dV->co[ch] = (cnet->Qsup->co[ch]+cnet->Qsub->co[ch])*dt; //m3				
				r = cnet->r->co[ch];
				c = cnet->c->co[ch];
				h->co[ch] += 1.E3*dV->co[ch]/(dn*cnet->length->co[ch]);	//mm
				*Vsubout = *Vsubout + dVsub[ch];	//m3

			}
			
			//  sedi transport calculation
		channel_sedimport(dt, h, Mcin, qch, cnet, par, top, Eout);
			  
		}while(te<Dt);
	}
}
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void EroParaInit(DOUBLEMATRIX *cl, DOUBLEMATRIX *dl, LAND *land, SOIL *sl)
{
    long r,c,sy;
	short lu;
	double ssg;                                             /* submerged specific gravity */
	double kinv;                                            /*kinematic viscosity [m2/s] */
	double para1,para2;                                     /*parameters used to calculate settling velocity*/
	double Kev,Ker;                                      //kinetic engery of rainfall and leaffall
	double Cohr;                                          // root additional cohesion  kPa

	kinv=0.000001;    // kinitic viscosity m2/s
	para1=24;        // parameter for settling velocity
	para2=1.2;	     // parameter for settling velocity
		
	for(r=1;r<=Nr;r++){
		for(c=1;c<=Nc;c++){
		    if(land->LC->co[r][c]!=NoV){
			    sy=sl->type->co[r][c];
			    lu=(short)land->LC->co[r][c];
			    Cohr=1.04*land->ty->co[lu][jvts]*land->ty->co[lu][jrar]*1000;     //1000   Mpa to Kpa
				ssg=sl->rho->co[1][r][c]/1000;
				sl->cohad->co[r][c]=sl->coh->co[1][r][c]*Fmin(pow((sl->th->co[1][r][c]/sl->pa->co[sy][jsat][1]),2.0),1.0)+Cohr;
				sl->ys->co[r][c]=0.79*exp(-0.85*(sl->cohad->co[r][c]));
			    //sl->ys->co[r][c]=1/(0.89+0.56*sl->cohad->co[r][c]);
               // printf("jsat=%d,ratio=%g\n",jsat,sl->th->co[1][r][c]/sl->pa->co[sy][jsat][1]);			  
			   sl->vs->co[r][c]=ssg*9.81*pow(sl->D50->co[1][r][c]*0.000001,2)
				/(para1*kinv+pow(0.75*para2*ssg*9.81*pow(sl->D50->co[1][r][c]*0.000001,3),0.5));
		
		   		if(sl->D50->co[1][r][c]>0){
			    cl->co[r][c]=pow(((sl->D50->co[1][r][c]+5)/0.32),-0.6);
		        dl->co[r][c]=pow(((sl->D50->co[1][r][c]+5)/300),0.25);				
				}else{cl->co[r][c]=0.0;
		              dl->co[r][c]=0.0;
				}
			}else{sl->ys->co[r][c]=0.0;
				  sl->vs->co[r][c]=0.0;}	 
		}
    }
	
}	

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void TranCap(DOUBLEMATRIX *TC,DOUBLEMATRIX *v, double **h, DOUBLEMATRIX *q, DOUBLEMATRIX *cl, DOUBLEMATRIX *dl, TOPO *top, LAND *land, SOIL *sl)
	{
	long r,c;
	short lu;
	double thetat,tau,tauc,omegar,omegac,beta; 	
	double dx,dy;
	double lth[10];
	
	dx=UV->U->co[1];                                    
	dy=UV->U->co[2];   
	
		lth[1]=0.0;
	    lth[2]=dx;
	    lth[3]=sqrt(dx*dx+dy*dy);
	    lth[4]=dy;
	    lth[5]=sqrt(dx*dx+dy*dy);
	    lth[6]=dx;
	    lth[7]=sqrt(dx*dx+dy*dy);
	    lth[8]=dy;
	    lth[9]=sqrt(dx*dx+dy*dy);	
		
		for(r=1;r<=Nr;r++){
			for(c=1;c<=Nc;c++){
				if(land->LC->co[r][c]!=NoV){ 						   
				   if(top->pixel_type->co[r][c]==0){ 
				    if(h[r][c]>0){
				   
                    //calculate the parameters used in Tc  g cm s					
                      thetat=0.2944*pow(sl->D50->co[1][r][c],-0.2834);	
                     // tauc=0.047*(sl->rho->co[1][r][c]-1000)*981*sl->D50->co[1][r][c]/1E7;					  
		              tauc=0.5*thetat*(sl->rho->co[1][r][c]-1000)*9.81*sl->D50->co[1][r][c]/1E5;  //g cm s2
				      tau=0.1*9.81*sl->rho->co[1][r][c]*h[r][c]*top->i_DD->co[r][c]/100;  
		              beta=(19-sl->D50->co[1][r][c]/30)/1.0E4;      //non-dimesion
					 
				      //tau->co[r][c]=0.0006*9.81*sl->rho->co[1][r][c]*h[r][c]*top->i_DD->co[r][c]/100;  //g  cm  s^2   0.006 is the ratio between tau surface and tau total
				      omegac=pow(tauc*v->co[r][c]*100,1.5)/pow(h[r][c]/10,2/3);
					  omegar=pow(tau*v->co[r][c]*100,1.5)/pow(h[r][c]/10,2/3);
				    //calculate Tc kg/m3	
				   if(top->i_DD->co[r][c]>0 && top->i_DD->co[r][c]<1 && q->co[r][c]>0 && h[r][c]>0.1){
			       if(100*top->i_DD->co[r][c]*v->co[r][c]>0.4){
                   TC->co[r][c]=cl->co[r][c]*pow((top->i_DD->co[r][c]*100*v->co[r][c]-0.4),dl->co[r][c])*sl->rho->co[1][r][c];
				   }else if(omegar-omegac>0){
				  //TC->co[r][c]=0;
				   TC->co[r][c]=beta/(q->co[r][c]*dx*dy/lth[top->DD->co[r][c]+1]/100)*pow(pow(omegar-omegac,0.7/5)-1,5);   //kg/m3
			       }else{TC->co[r][c]=0.0;}				   
                   }if(TC->co[r][c]<0.0)TC->co[r][c]=0.0;
				   
				   }else{TC->co[r][c]=0.0;}													
				  }			
				}
			}
		}
	}	
		
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void D_Cal(double **h, DOUBLEMATRIX *D, DOUBLEMATRIX *TC,DOUBLEMATRIX *q, DOUBLEMATRIX *Conc,TOPO *top, LAND *land, WATER *wat, SOIL *sl, PAR *par)
	{
	long r,c,R,C;
	double dx,dy;                     // the two dimensions of a pixel
	double b[10],lth[10];            // area perpendicular to the superficial flow divided by hsup
	short lu;
	double rf;                       //fraction of rill area
	double Ke,Ker,Kev,Ds,Dfp;
	
		dx=UV->U->co[1];                                      
		dy=UV->U->co[2]; 
		
        lth[1]=0.0;  
		lth[2]=dy;             
		lth[3]=sqrt(dx*dx+dy*dy); 
		lth[4]=dx;            
		lth[5]=sqrt(dx*dx+dy*dy); 
		lth[6]=dy;   
		lth[7]=sqrt(dx*dx+dy*dy); 
		lth[8]=dx;             
		lth[9]=sqrt(dx*dx+dy*dy);	
	
		
		b[1]=0.0;  
		b[2]=dy*par->wfrac;             
		b[3]=sqrt(dx*dx+dy*dy)*par->wfrac; 
		b[4]=dx*par->wfrac;            
		b[5]=sqrt(dx*dx+dy*dy)*par->wfrac; 
		b[6]=dy*par->wfrac;   
		b[7]=sqrt(dx*dx+dy*dy)*par->wfrac; 
		b[8]=dx*par->wfrac;             
		b[9]=sqrt(dx*dx+dy*dy)*par->wfrac; 
	
		for(r=1;r<=Nr;r++){
			for(c=1;c<=Nc;c++){
				if(land->LC->co[r][c]!=NoV){
                   if(h[r][c]>0){				
				    lu=(short)land->LC->co[r][c];
					
                   // Calculate KE  J*mm-1*m-2				
				    if(wat->Pnet->co[r][c]!=0){
                 	Ker=8.95+8.44*log(wat->Pnet->co[r][c]*3600);   //  I mm/h
                    Kev=15.8*pow(land->ty->co[lu][jHveg]/1000,0.5)-5.87;	
                    if(Kev<0.0)Kev=0.0;					
                    Ke=land->ty->co[lu][jcf]*Kev+(1-land->ty->co[lu][jcf])*Ker;		//	  KE  J*m-2	*mm-1 
				    }else {Ke=0.0;}
				   
				   //calculate Ds kg m-2 s-1 
				    rf=b[top->DD->co[r][c]+1]*lth[top->DD->co[r][c]+1]/dx/dy;
				    if(rf>1){rf=1;}				  
				    Ds=0.001*((0.1033/sl->cohad->co[r][c]*Ke*exp(-1.48*(h[r][c]*lth[top->DD->co[r][c]+1]*
				    b[top->DD->co[r][c]+1]/dx/dy))+3.58)*(1-rf)+(0.1033/sl->cohad->co[r][c]*
					Ker*exp(-1.48*h[r][c])+3.58)*rf)*wat->Pnet->co[r][c];    					   
				    if(Ds<0){Ds=0;}
					
				  	//calculate Min kg/s
					wat->Min->co[r][c]=0.0;
					if(top->DD->co[r][c-1]==1 && top->pixel_type->co[r][c-1]==0) 
				    wat->Min->co[r][c]+=q->co[r][c-1]*Conc->co[r][c-1]*dx*dy/1000;  //  kg/s   
				    if(top->DD->co[r+1][c-1]==2 && top->pixel_type->co[r+1][c-1]==0) 
				    wat->Min->co[r][c]+=q->co[r+1][c-1]*Conc->co[r+1][c-1]*dx*dy/1000;
					if(top->DD->co[r+1][c]==3 && top->pixel_type->co[r+1][c]==0) 
				    wat->Min->co[r][c]+=q->co[r+1][c]*Conc->co[r+1][c]*dx*dy/1000;
					if(top->DD->co[r+1][c+1]==4 && top->pixel_type->co[r+1][c+1]==0) 
				    wat->Min->co[r][c]+=q->co[r+1][c+1]*Conc->co[r+1][c+1]*dx*dy/1000;
					if(top->DD->co[r][c+1]==5 && top->pixel_type->co[r][c+1]==0) 
				    wat->Min->co[r][c]+=q->co[r][c+1]*Conc->co[r][c+1]*dx*dy/1000;
					if(top->DD->co[r-1][c+1]==6 && top->pixel_type->co[r-1][c+1]==0) 
				    wat->Min->co[r][c]+=q->co[r-1][c+1]*Conc->co[r-1][c+1]*dx*dy/1000;
					if(top->DD->co[r-1][c]==7 && top->pixel_type->co[r-1][c]==0) 
				    wat->Min->co[r][c]+=q->co[r-1][c]*Conc->co[r-1][c]*dx*dy/1000;
					if(top->DD->co[r-1][c-1]==8 && top->pixel_type->co[r-1][c-1]==0) 
				    wat->Min->co[r][c]+=q->co[r-1][c-1]*Conc->co[r-1][c-1]*dx*dy/1000;
					
                   //calculate Dfp  kg/s
				   if(top->pixel_type->co[r][c]==0){	
			        if(TC->co[r][c]>Conc->co[r][c]){
				    Dfp=(TC->co[r][c]-Conc->co[r][c])*sl->ys->co[r][c]*sl->vs->co[r][c]*b[top->DD->co[r][c]+1];
					}else if(TC->co[r][c]<Conc->co[r][c]){
				    Dfp=(TC->co[r][c]-Conc->co[r][c])*sl->vs->co[r][c]*b[top->DD->co[r][c]+1];	
                    }else if(TC->co[r][c]==Conc->co[r][c]){
                    Dfp=0;							
					}			
				    }
					
                    D->co[r][c]=Ds*dx*dy+wat->Min->co[r][c]+Dfp*lth[top->DD->co[r][c]+1];				 //kg/s	
				    /*if(D->co[r][c]>100){
				    printf("haha\n");
				    } */
				   }else{
				   D->co[r][c]=0;
				   }
                  }				   
				}
			}
		}
		
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void Erosion(double dt, double **h, DOUBLEMATRIX *Conc, DOUBLEMATRIX *q, DOUBLEMATRIX *D, TOPO *top, LAND *land, WATER *wat, SOIL *sl, PAR *par, CHANNEL *cnet,double *OCout){
		
	long r,c,R,C,ch;                                    
	static short r_DD[11]={0,0,-1,-1,-1,0,1,1,1,-9999,0}; // differential of number-pixel for rows and
	static short c_DD[11]={0,1,1,0,-1,-1,-1,0,1,-9999,0}; // columns, depending on Drainage Directions
	double dx,dy;                                         // the two dimensions of a pixel
	double b[10],lth[10];                                 // area perpendicular to the superficial flow divided by hsup
	short lu;
	double outero,outpoc,outdoc;
		
		dx=UV->U->co[1];                                      
		dy=UV->U->co[2]; 
		
        lth[1]=0.0;  
		lth[2]=dy;             
		lth[3]=sqrt(dx*dx+dy*dy); 
		lth[4]=dx;            
		lth[5]=sqrt(dx*dx+dy*dy); 
		lth[6]=dy;   
		lth[7]=sqrt(dx*dx+dy*dy); 
		lth[8]=dx;             
		lth[9]=sqrt(dx*dx+dy*dy);		
		
	    b[1]=0.0;  
		b[2]=dy*par->wfrac;             
		b[3]=sqrt(dx*dx+dy*dy)*par->wfrac; 
		b[4]=dx*par->wfrac;            
		b[5]=sqrt(dx*dx+dy*dy)*par->wfrac; 
		b[6]=dy*par->wfrac;   
		b[7]=sqrt(dx*dx+dy*dy)*par->wfrac; 
		b[8]=dx*par->wfrac;             
		b[9]=sqrt(dx*dx+dy*dy)*par->wfrac; 
		
		
		outero=0.0;
	    outpoc=0.0;
	    outdoc=0.0;
		
    //  Calculate concentration and erosion and deposition on each cell
	for(r=1;r<=Nr;r++){
	    for(c=1;c<=Nc;c++){
		   if(land->LC->co[r][c]!=NoV){			    
			 if(top->pixel_type->co[r][c]==0){
			    if(h[r][c] > 0.1 && q->co[r][c]>0){			 
     			  Conc->co[r][c]=(1000*D->co[r][c]*dt+Conc->co[r][c]*h[r][c]*lth[top->DD->co[r][c]+1]*b[top->DD->co[r][c]+1])/(lth[top->DD->co[r][c]+1]*b[top->DD->co[r][c]+1]*h[r][c]+q->co[r][c]*dx*dy*dt);
              	  if(Conc->co[r][c]<0){
				  Conc->co[r][c]=0;

				  } 
				  wat->dm->co[r][c]=wat->Min->co[r][c]*dt/(dx*dy/1.0E4)-q->co[r][c]*Conc->co[r][c]*dx*dy*dt/1.0E3/(dx*dy/1.0E4); 
				  if(wat->dm->co[r][c]>0){wat->Mp->co[r][c]+=wat->dm->co[r][c];}
				  if(wat->dm->co[r][c]<0){wat->Mc->co[r][c]+=wat->dm->co[r][c];}
				  wat->Mn->co[r][c]=-wat->Mc->co[r][c]-wat->Mp->co[r][c];  
                  sl->POC->co[r][c]=wat->Mn->co[r][c]*sl->pocc->co[r][c]*sl->er_r->co[1][r][c];
				  sl->DOC->co[r][c]+=dt*(q->co[r][c])*dx*dy*sl->docc->co[r][c]/1.0E6/(dx*dy/1.0E4);					
			      sl->EOC->co[r][c]=sl->POC->co[r][c]+sl->DOC->co[r][c];
			   }
			 } 
		   }
		}
	}
	
		
	    for(r=1;r<=Nr;r++){
		    for(c=1;c<=Nc;c++){
			    R=r+r_DD[top->DD->co[r][c]];
				C=c+c_DD[top->DD->co[r][c]];
				if(land->LC->co[r][c]!=NoV){
				 if(h[r][c] > 0){
				  if(top->pixel_type->co[r][c]==0){
			        if(top->pixel_type->co[R][C]==10){
			          	outero+=dt*q->co[r][c]*dx*dy*Conc->co[r][c]*1.0E-3;   //kg
						outdoc+=dt*(q->co[r][c])*dx*dy*sl->docc->co[r][c]*1.0E-6;   //kg
						outpoc=outero*sl->pocc->co[r][c]*sl->er_r->co[1][r][c];  //kg
					}                
			      }			     
		        }else{
				outero+=0;
				outdoc+=0;
				outpoc+=0;
				}
				/*if(outero!=outero){
				printf("error\n");
				}*/
			  }
			}
        }	
	

     		//*Eout+=outero/(par->total_pixel*dx*dy/1.0E4);  //kg/ha
			*OCout+=(outdoc+outpoc)/(par->total_pixel*dx*dy/1.0E4);   //kg/ha
			
		//  Sediment flow into channels

         for(ch=1;ch<=cnet->r->nh;ch++){
				r = cnet->r->co[ch];
				c = cnet->c->co[ch];
				
				if(r>0){
				
				    cnet->Inmass->co[ch]=0.0;
					if(top->DD->co[r][c-1]==1 && top->pixel_type->co[r][c-1]==0) 
				    cnet->Inmass->co[ch]+=q->co[r][c-1]*Conc->co[r][c-1]*dx*dy/1000;  //  kg/s   
				    if(top->DD->co[r+1][c-1]==2 && top->pixel_type->co[r+1][c-1]==0) 
				    cnet->Inmass->co[ch]+=q->co[r+1][c-1]*Conc->co[r+1][c-1]*dx*dy/1000;
					if(top->DD->co[r+1][c]==3 && top->pixel_type->co[r+1][c]==0) 
				    cnet->Inmass->co[ch]+=q->co[r+1][c]*Conc->co[r+1][c]*dx*dy/1000;
					if(top->DD->co[r+1][c+1]==4 && top->pixel_type->co[r+1][c+1]==0) 
				    cnet->Inmass->co[ch]+=q->co[r+1][c+1]*Conc->co[r+1][c+1]*dx*dy/1000;
					if(top->DD->co[r][c+1]==5 && top->pixel_type->co[r][c+1]==0) 
				    cnet->Inmass->co[ch]+=q->co[r][c+1]*Conc->co[r][c+1]*dx*dy/1000;
					if(top->DD->co[r-1][c+1]==6 && top->pixel_type->co[r-1][c+1]==0) 
				    cnet->Inmass->co[ch]+=q->co[r-1][c+1]*Conc->co[r-1][c+1]*dx*dy/1000;
					if(top->DD->co[r-1][c]==7 && top->pixel_type->co[r-1][c]==0) 
				    cnet->Inmass->co[ch]+=q->co[r-1][c]*Conc->co[r-1][c]*dx*dy/1000;
					if(top->DD->co[r-1][c-1]==8 && top->pixel_type->co[r-1][c-1]==0) 
				    cnet->Inmass->co[ch]+=q->co[r-1][c-1]*Conc->co[r-1][c-1]*dx*dy/1000;
								
				}
				
			}
		 


		
	
}	

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void channel_sedimport(double Dt, DOUBLEVECTOR *h, DOUBLEVECTOR *Mcin, DOUBLEVECTOR *qch,CHANNEL *cnet, PAR *par, TOPO *top, double *Eout)

{
    static short r_DD[11]={0,0,-1,-1,-1,0,1,1,1,-9999,0}; // differential of number-pixel for rows and
	static short c_DD[11]={0,1,1,0,-1,-1,-1,0,1,-9999,0}; // columns, depending on Drainage Directions
	long r, c, r_down, c_down, ch, ch_down;
	double dx,dy,dn;                                         // the two dimensions of a channel pixel
			
		dn = par->w_dx * UV->U->co[1];		//transversal length [m]
		dx=UV->U->co[1];                                      
		dy=UV->U->co[2];    
		
			for(ch=1;ch<=cnet->r->nh;ch++){
			        r = cnet->r->co[ch];
					c = cnet->c->co[ch];
					r_down = r+r_DD[top->DD->co[r][c]];
					c_down = c+c_DD[top->DD->co[r][c]];
					ch_down = cnet->ch->co[r_down][c_down];
					if(ch_down==0) ch_down=ch;
					
				if(h->co[ch] > 0){			 
     			  cnet->Conc->co[ch]=Mcin->co[ch]/(qch->co[ch]*Dt+h->co[ch]*1.0E-3*dn*cnet->length->co[ch]);
              	  if(cnet->Conc->co[ch]<0){
				  cnet->Conc->co[ch]=0;
				  } 
				  
				}
				if (ch_down == ch){
				
					 *Eout += cnet->Conc->co[ch]*qch->co[ch]*Dt;	// kg   total suspended sediments loss at the outlet

				//outdoc+=dt*(q->co[r][c])*dx*dy*sl->docc->co[r][c]*1.0E-6;   //kg
				//outpoc=outero*sl->pocc->co[r][c]*sl->er_r->co[1][r][c];  //kg
				}
			}			
}
