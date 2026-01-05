/* To use function from matlab, first compile by entering this into the Matlab command window:
   mex LNP_RecNet.c
   Then call the function like this:
   [s1,Vm_save,Vm_ave,IsynE_save,IsynI_save,time_save]= ...
    LNP_ClusterNet(Vm0,Input,PostID,out_degree,Wrr_c,param);  
 *
 * mesure within-cluster currents & out-of-cluster currents 
 */


#include "mex.h"
#include "math.h"
#include "time.h"
#include "matrix.h"



/* A fast approximation of the exp function */
static union 
{
	  double d;
	    struct {
#ifdef LITTLE_ENDIAN
	    int j,i;
#else 
		    int i,j;
#endif
	  } n;
} _eco;
#define EXP_A (1048576/0.69314718055994530942)
#define EXP_C 60801
#define EXP(y) (_eco.n.i = EXP_A*(y) + (1072693248 - EXP_C), _eco.d)

        
/* #define Phi(y)  0.3*fmax(y,0)*fmax(y,0);  */  /* transfer function  */ 

#define Phi(y)  30/(1+EXP(-((y+60)-15)/3));   /* transfer function  */ 
/*  param.phi=@(v) 30./(1+exp(-((v+60)-15)/3)); */


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

   mxArray *tmp,*Input[2],*rhs_rand,*lhs_input; /*  *phi[2],*lhs_phi;  */
   double A, *In,dt,*tau_m,tau_noise,tau_syn,*PostID,*out_degree,*Wrr,b; 
   double *Vm0,rho,T,*sigmaN, *idx_save,*s,*Vm_save,*Vm_ave,*IsynE_save,*IsynI_save,*eta_save,*v,*IsynE,*IsynI,*eta,*noise,*noise1,*noise2;
   double t, *Iapp, rate,*Size,*Vm,*PostID_bd,TMP,*time_save;
   int m1, m2,j,Ne,Ni, N, Nw, maxns, Nw2, Nt, Nsave,jj,ns,i,k,Nskip,Nt_save,count;
   const mxArray  *mxTmp;
   
   /******
 * Import variables from matlab
 *******/
   
 /* Number of exc neurons  */ 
 /*  b=10; 
   A=Phi(b);
   mexPrintf("\n Phi(%f)= %f \n",b,A); 
    mexErrMsgTxt("test"); */ 
   
mxTmp = mxGetField(prhs[5],0,"Ne");
Ne=(int)mxGetPr(mxTmp)[0];

/* Number of inh neurons */ 
mxTmp = mxGetField(prhs[5],0,"Ni");
Ni=(int)mxGetPr(mxTmp)[0];
  
 /* Total number of neurons  */ 
mxTmp = mxGetField(prhs[5],0,"N");
N=(int)mxGetPr(mxTmp)[0];

/* Initial conditions of Vm  */ 
Vm0 = mxGetPr(prhs[0]);
m1 = mxGetM(prhs[0]);
m2 = mxGetN(prhs[0]);
if(m1!=N){
    mexErrMsgTxt("Vm0 should be Nx1");
}

/* input function handle */ 
if( !mxIsClass( prhs[1] , "function_handle")) {
    mexErrMsgTxt("Second input argument is not a function handle.");
}
Input[0] = (mxArray *)(prhs[1]); 

/* index of post synaptic neurons */ 
PostID = mxGetPr(prhs[2]);
Nw = mxGetM(prhs[2]);
m2 = mxGetN(prhs[2]);
/* mexPrintf("\n PostID size: %d by %d \n",Nw, m2); */
/*  mexPrintf("\n PostID \n",(int)PostID[0],(int)PostID[1],(int)PostID[1]); */

/* outdegree for each neuron */ 
out_degree = mxGetPr(prhs[3]);
m1 = mxGetM(prhs[3]);
m2 = mxGetN(prhs[3]);
if(m1*m2!=N){
    mexErrMsgTxt("out_degree should be 1xN");
}

/* connection strengths */ 
Wrr = mxGetPr(prhs[4]);
Nw2 = mxGetM(prhs[4]);
m2 = mxGetN(prhs[4]);
if(Nw!=Nw2){
    mexErrMsgTxt("Wrr and PostID should be the same size");
}

/* transfer function handle */ 

/*  if( !mxIsClass( mxGetField(prhs[5],0,"phi") , "function_handle")) {
//    mexErrMsgTxt("param.phi is not a function handle.");  }
// phi[0] = (mxArray *) mxGetField(prhs[5],0,"phi");   */

 /* time step size */ 
mxTmp = mxGetField(prhs[5],0,"dt");
dt=mxGetPr(mxTmp)[0];


/* membrane potential time constants */ 
mxTmp = mxGetField(prhs[5],0,"tau_m");
tau_m=mxGetPr(mxTmp);

/* noise time constant */ 
/* mxTmp = mxGetField(prhs[5],0,"tau_noise");
tau_noise=mxGetPr(mxTmp)[0];
mexPrintf("\n tau_noise: %f  \n",tau_noise); */


/* Synaptic time constant */ 
mxTmp = mxGetField(prhs[5],0,"tau_syn");
tau_syn=mxGetPr(mxTmp)[0];

/* correlation of input noise */ 
/* mxTmp = mxGetField(prhs[5],0,"rho");
rho=mxGetPr(mxTmp)[0];
mexPrintf("\n rho: %f  \n",rho);  */ 

/* total simulation time */ 
mxTmp = mxGetField(prhs[5],0,"T");
T=mxGetPr(mxTmp)[0];

/* Noise standard deviation  */ 
/* mxTmp = mxGetField(prhs[5],0,"sigma_n");
sigmaN=mxGetPr(mxTmp); */

/* maximum number of spikes to store */ 
mxTmp = mxGetField(prhs[5],0,"maxns");
maxns=(int)mxGetPr(mxTmp)[0];

/* number of neurons to save */ 
mxTmp = mxGetField(prhs[5],0,"Nsave");
Nsave=(int)mxGetPr(mxTmp)[0];

/* Index of neurons to save */ 
mxTmp = mxGetField(prhs[5],0,"idx_save");
idx_save=mxGetPr(mxTmp);
m1 = mxGetN(mxTmp);
m2 = mxGetM(mxTmp);
if(m1*m2!=Nsave)
    mexErrMsgTxt("idx_save should be 1xNsave.");
if (mxTmp == NULL) {
    mexErrMsgTxt("Field 'idx_save' not found.");
}
if (!mxIsDouble(mxTmp)) {
    mexErrMsgTxt("idx_save must be of type double.");
}
if (mxGetNumberOfElements(mxTmp) != Nsave) {
    mexErrMsgTxt("idx_save must contain exactly Nsave elements.");
}
idx_save = mxGetPr(mxTmp);
/* number of neurons to save */ 
mxTmp = mxGetField(prhs[5],0,"Nskip");
Nskip=(int)mxGetPr(mxTmp)[0];

/* Numebr of time bins */
Nt=(int)(T/dt)+1; 

Nt_save=floor(Nt/Nskip); 

/* Allocate output vector */
plhs[0] = mxCreateDoubleMatrix(2, maxns, mxREAL);
s=mxGetPr(plhs[0]);

plhs[1] = mxCreateDoubleMatrix(Nsave, Nt_save, mxREAL);
Vm_save=mxGetPr(plhs[1]);  
m1 = mxGetM(plhs[1]);
m2 = mxGetN(plhs[1]);

plhs[2] = mxCreateDoubleMatrix(2, Nt_save, mxREAL);
Vm_ave=mxGetPr(plhs[2]);

plhs[3] = mxCreateDoubleMatrix(Nsave, Nt_save, mxREAL);
IsynE_save=mxGetPr(plhs[3]);  

plhs[4] = mxCreateDoubleMatrix(Nsave, Nt_save, mxREAL);
IsynI_save=mxGetPr(plhs[4]); 

/*  plhs[5] = mxCreateDoubleMatrix(Nsave, Nt_save, mxREAL);
eta_save=mxGetPr(plhs[5]);  */ 

plhs[5] = mxCreateDoubleMatrix(1, Nt_save, mxREAL);
time_save=mxGetPr(plhs[5]); 

rhs_rand = mxCreateDoubleMatrix(1, 2, mxREAL);
Size=mxGetPr(rhs_rand);  


/* Allocate membrane potential */
Vm = mxMalloc(N*sizeof(double));

/* synaptic currents   */ 
IsynE = mxMalloc(N*sizeof(double));  
IsynI = mxMalloc(N*sizeof(double));  

/* input noise */
/* eta = mxMalloc(N*sizeof(double));  */ 

PostID_bd = mxMalloc((N+1)*sizeof(double));;
PostID_bd[0]=0; 
TMP=0;
for(j=0;j<N;j++){
    TMP = TMP + out_degree[j]; 
    PostID_bd[j+1]=TMP;
}

if(Nw!=(int)PostID_bd[N]){
    mexErrMsgTxt(" size(PostID) should be the sum of out_degree");
}

/* Inititalize variables */

Size[0]= N;
Size[1]= 1;
mexCallMATLABWithTrap(1,&tmp,1,&rhs_rand,"randn");
noise=mxGetPr(tmp);

for(j=0;j<N;j++){
    Vm[j]=Vm0[j];
    IsynE[j]=0;
    IsynI[j]=0;
  /*  eta[j]=noise[j]*sigmaN[j]; */
}

/* save variables */
count=0;
if(Nskip==1){
    for(jj=0;jj<Nsave;jj++){
        time_save[count]=0;
        /* mexPrintf("\n jj %d idx_save[jj] %d, v %f \n",jj, (int)round(idx_save[jj]-1),v[(int)round(idx_save[jj]-1)]);  */
        if(idx_save[jj]<1 || idx_save[jj]>N){
            mexErrMsgTxt("Indices in Irecord must be between 1 and N");}
        Vm_save[jj]=Vm[(int)round(idx_save[jj]-1)];
        IsynE_save[jj]=IsynE[(int)round(idx_save[jj]-1)];
        IsynI_save[jj]=IsynI[(int)round(idx_save[jj]-1)];
     /*   eta_save[jj]=eta[(int)round(idx_save[jj]-1)];  */
    }
    
    Vm_ave[0]=0;
    Vm_ave[1]=0;
    for(j=0;j<Ne;j++){
        Vm_ave[0]=Vm_ave[0]+Vm[j];}
    for(j=Ne;j<N;j++){
        Vm_ave[1]=Vm_ave[1]+Vm[j];}
    Vm_ave[0]=Vm_ave[0]/Ne;
    Vm_ave[1]=Vm_ave[1]/Ni;
    count++;
}

ns=0;
  

for(i=1;i<Nt && ns<maxns;i++){  
    t=i*dt; 
    /*mexPrintf("\n t %f  \n",t); */ 
    Input[1]=mxCreateDoubleScalar(t);
    mexCallMATLABWithTrap(1,&lhs_input,2,Input,"feval");
    Iapp = mxGetPr(lhs_input);
  
    /*
    Input[1] = mxCreateDoubleScalar(t);
    mxArray *err = mexCallMATLABWithTrap(1, &lhs_input, 2, Input, "feval");
    if (err != NULL) {
    mexErrMsgTxt("Error in feval: failed to evaluate the input function.");
    }
    
    Iapp = mxGetPr(lhs_input);
    mwSize len = mxGetNumberOfElements(lhs_input);
    if (len < N) {
    mexPrintf("feval returned Iapp with only %d elements, but expected N = %d\n", (int)len, N);
    mexErrMsgTxt("Input function output size too small.");
}
*/
  /*  Size[0]= N+1;
    Size[1]= 1;
    mexCallMATLABWithTrap(1,&tmp,1,&rhs_rand,"randn");
    noise=mxGetPr(tmp);  */ 

  /*  for(j=0;j<N;j++){
        eta[j]=eta[j] - dt*eta[j]/tau_noise + sqrt(2/tau_noise)*sigmaN[j]*sqrt(dt)*(sqrt(1-rho) * noise[j] + sqrt(rho)*noise[N]); 
    }  */
    
    
    Size[0]= N;
    Size[1]= 1;
    mexCallMATLABWithTrap(1,&tmp,1,&rhs_rand,"rand");
    noise=mxGetPr(tmp);
  /*  m1 = mxGetM(tmp);
    m2 = mxGetN(tmp);
    mexPrintf("\n noise:  %d by %d \n",m1, m2); 
    mexPrintf("\n rand %f, %f, %f \n",noise[0], noise[1], noise[2]);  */
    for(j=0;j<N;j++){ 
       Vm[j]=Vm[j] + dt/tau_m[j]*( -(Vm[j]+60) + Iapp[j]  + IsynE[j] + IsynI[j] );}
    
    for(j=0;j<N;j++){ 
        IsynE[j] = IsynE[j] - dt*IsynE[j]/tau_syn; 
        IsynI[j] = IsynI[j] - dt*IsynI[j]/tau_syn; 
    }
    
    for(j=0;j<N;j++){ 
        /* check spike */
    /*   phi[1]=mxCreateDoubleScalar(Vm[j]);
         mexCallMATLABWithTrap(1,&lhs_phi,2,phi,"feval");
         rate=*mxGetPr(lhs_phi);    */
       rate= Phi(Vm[j]); 
     /* mexPrintf("\n Vm(j) %f, rate %f \n",Vm[j],rate);  */
       if(noise[j]<rate*dt*0.001 && ns<maxns){
     /*      mexPrintf("\n j %d, Vm %f, rate %f \n",j,Vm[j],rate);  */
           s[0+2*ns]=t;  /* spike time  */
           s[1+2*ns]=j+1;     /* neuron index 1 */
           ns++;
           /* update current */ 
           if(j<Ne){
               for(k=(int)round(PostID_bd[j]);k<(int)round(PostID_bd[j+1]);k++){
                   if(PostID[k]<1 || PostID[k]>N){
                       mexPrintf("\n j %d, k %d, postID %d \n",j,k,(int)round(PostID[k])-1);
                        mexErrMsgTxt("Indices in PostID must be between 1 and N");}
                  
                   IsynE[(int)round(PostID[k])-1]+=Wrr[k]; }}
           /*     mexPrintf("\n j %d, k %d, postID %d, Wrr(k) %f \n",j,k,(int)round(PostID[k]),Wrr[k]);  */
           else{
               for(k=(int)round(PostID_bd[j]);k<(int)round(PostID_bd[j+1]);k++){
                   if(PostID[k]<1 || PostID[k]>N){
                       mexPrintf("\n j %d, k %d, postID %d \n",j,k,(int)round(PostID[k])-1);
                        mexErrMsgTxt("Indices in PostID must be between 1 and N");}
                  
                   IsynI[(int)round(PostID[k])-1]+=Wrr[k];}}  } 
       
    }
    /* save variables */
    if(i%Nskip==0){
     /*   mexPrintf("\n i %d, count %d t %f \n",i,count,t);  */
        time_save[count]=t; 
    for(jj=0;jj<Nsave;jj++){
    /* mexPrintf("\n jj %d idx_save[jj] %d, v %f \n",jj, (int)round(idx_save[jj]-1),v[(int)round(idx_save[jj]-1)]); */
        if(idx_save[jj]<1 || idx_save[jj]>N){
            mexErrMsgTxt("Indices in Irecord must be between 1 and N");}
        Vm_save[jj+count*Nsave]=Vm[(int)round(idx_save[jj]-1)];
        IsynE_save[jj+count*Nsave]=IsynE[(int)round(idx_save[jj]-1)];
        IsynI_save[jj+count*Nsave]=IsynI[(int)round(idx_save[jj]-1)];
     /*   eta_save[jj+count*Nsave]=eta[(int)round(idx_save[jj]-1)];  */
    }
    Vm_ave[0+count*2]=0;
    Vm_ave[1+count*2]=0;
    for(j=0;j<Ne;j++){
        Vm_ave[0+count*2] += Vm[j];}
    for(j=Ne;j<N;j++){
        Vm_ave[1+count*2] += Vm[j];}
    Vm_ave[0+count*2]=Vm_ave[0+count*2]/Ne;
    Vm_ave[1+count*2]=Vm_ave[1+count*2]/Ni;
    count++;  
 /*  mexPrintf("\n Vm_ave(1,i) %f, Vm_ave(2,i) %f, %d\n",Vm_ave[0+i*2],Vm_ave[1+i*2],i);  */
    } 
    
}


 
if(ns>=maxns)
   mexWarnMsgTxt("Maximum number of spikes reached, simulation terminated.");
         

/* Clean UP  */ 
mxFree(Vm);
/* mxFree(eta); */ 
mxFree(IsynE);
mxFree(IsynI);
mxFree(PostID_bd);


mxDestroyArray(tmp);
mxDestroyArray(rhs_rand);
mxDestroyArray(lhs_input);
mxDestroyArray(Input[1]);

/*
//mxDestroyArray(phi[1]);
//mxDestroyArray(lhs_phi); */

}












