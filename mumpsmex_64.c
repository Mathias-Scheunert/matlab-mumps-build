#include "mex.h"

#define MUMPS_ARITH_d 2
#define MUMPS_ARITH_z 8


#if MUMPS_ARITH == MUMPS_ARITH_z

# include "zmumps_c.h"
# define dmumps_c       zmumps_c
# define dmumps_par     zmumps_par
# define DMUMPS_STRUC_C ZMUMPS_STRUC_C
# define DMUMPS_alloc   ZMUMPS_alloc     
# define DMUMPS_free    ZMUMPS_free
# define double2        mumps_double_complex
# define mxREAL2        mxCOMPLEX

#elif MUMPS_ARITH == MUMPS_ARITH_d

# include "dmumps_c.h"
# define double2 double
# define mxREAL2 mxREAL
# define EXTRACT_CMPLX_FROM_C_TO_MATLAB EXTRACT_FROM_C_TO_MATLAB
# define EXTRACT_CMPLX_FROM_MATLAB_TOPTR EXTRACT_FROM_MATLAB_TOPTR

#else

# error "Only d and z arithmetics are supported"

#endif

#define SYM        (prhs[0])
#define JOB        (prhs[1])
#define ICNTL_IN   (prhs[2])
#define CNTL_IN    (prhs[3])
#define PERM_IN    (prhs[4])
#define COLSCA_IN  (prhs[5])
#define ROWSCA_IN  (prhs[6])
#define RHS_IN     (prhs[7])
#define VAR_SCHUR  (prhs[8])
#define INST       (prhs[9])
#define REDRHS_IN  (prhs[10])
#define KEEP_IN    (prhs[11])
#define DKEEP_IN   (prhs[12])
#define A_IN       (prhs[13])

#define INFO_OUT   (plhs[0])
#define RINFO_OUT  (plhs[1])
#define RHS_OUT    (plhs[2])
#define INST_OUT   (plhs[3])
#define SCHUR_OUT  (plhs[4])
#define REDRHS_OUT (plhs[5])
#define PIVNUL_LIST (plhs[6])
#define PERM_OUT   (plhs[7])
#define UNS_PERM   (plhs[8])
#define ICNTL_OUT  (plhs[9])
#define CNTL_OUT   (plhs[10])
#define COLSCA_OUT (plhs[11])
#define ROWSCA_OUT (plhs[12])
#define KEEP_OUT   (plhs[13])
#define DKEEP_OUT  (plhs[14])

#define MYMALLOC(ptr,l,type)                      \
  if(!(ptr = (type *) malloc(l*sizeof(type)))){   \
    mexErrMsgTxt ("Malloc failed in mumpsmex.c"); \
  }                                               

#define MYFREE(ptr) \
    if(ptr){        \
        free(ptr);  \
        ptr = 0;    \
    }

#define EXTRACT_FROM_MATLAB_TOPTR(mxcomponent,mumpspointer,type,length)         \
  ptr_matlab = mxGetPr(mxcomponent);                                            \
  if(ptr_matlab[0] != -9999){                                                   \
    MYFREE(mumpspointer);                                                       \
    MYMALLOC(mumpspointer,length,type);                                         \
    for(i=0;i<length;i++){                                                      \
      mumpspointer[i] = ptr_matlab[i];                                          \
    }                                                                           \
  }


/* For scaling arrays, if they were previously allocated by MUMPS, touch nothing   */
/* This is not quite correct (user may want to modify MUMPS scaling and use given  */
/* scaling, or provide a new scaling vector on input after a previous call where   */
/* it was computed by MUMPS). One way to solve this might be to separate COLSCA_IN */
/* and COLSCA_OUT in the C interface (and possibly Fortran) too, but breaking      */
/* backward compatibility.                                                         */
#define EXTRACT_SCALING_FROM_MATLAB_TOPTR(mxcomponent,mumpspointer,is_a_pointer_from_mumps,length)   \
  ptr_matlab = mxGetPr(mxcomponent);                                            \
  if( ptr_matlab[0] != -9999 && ! (is_a_pointer_from_mumps) ) {                 \
    MYFREE(mumpspointer);                                                       \
    MYMALLOC(mumpspointer,length,double);                                       \
    for(i=0;i<length;i++){                                                      \
      mumpspointer[i] = ptr_matlab[i];                                          \
    }                                                                           \
  }

#define EXTRACT_FROM_MATLAB_TOARR(mxcomponent,mumpsarray,type,length)           \
  ptr_matlab = mxGetPr(mxcomponent);                                            \
  if(ptr_matlab[0] != -9999){                                                   \
    for(i=0;i<length;i++){                                                      \
      if(ptr_matlab[i] != -9998){                                               \
        mumpsarray[i] = ptr_matlab[i];                                          \
      }                                                                         \
    }                                                                           \
  }

#define EXTRACT_FROM_MATLAB_TOVAL(mxcomponent,mumpsvalue)                       \
  ptr_matlab = mxGetPr(mxcomponent);                                            \
  if(ptr_matlab[0] != -9999){                                                   \
      mumpsvalue = ptr_matlab[0];                                               \
  }

#define EXTRACT_FROM_C_TO_MATLAB(mxcomponent,mumpspointer,length)               \
  if(mumpspointer == 0){                                                        \
    mxcomponent = mxCreateDoubleMatrix (1, 1, mxREAL);                          \
    ptr_matlab = mxGetPr (mxcomponent);                                         \
    ptr_matlab[0] = -9999;                                                      \
  }else{                                                                        \
    mxcomponent = mxCreateDoubleMatrix (1,length,mxREAL);                       \
    ptr_matlab = mxGetPr (mxcomponent);                                         \
    for(i=0;i<length;i++){                                                      \
      ptr_matlab[i]=(double)(mumpspointer)[i];                                  \
    }                                                                           \
  }

#if MUMPS_ARITH == MUMPS_ARITH_z

#define EXTRACT_CMPLX_FROM_MATLAB_TOPTR(mxcomponent,mumpspointer,type,length)   \
  ptr_matlab = mxGetPr(mxcomponent);                                            \
  if(ptr_matlab[0] != -9999){                                                   \
    MYFREE(mumpspointer);                                                       \
    MYMALLOC(mumpspointer,length,double2);                                      \
    for(i=0;i<length;i++){                                                      \
      (mumpspointer[i]).r = ptr_matlab[i];                                      \
    }                                                                           \
    ptr_matlab = mxGetPi(mxcomponent);                                          \
    if(ptr_matlab){                                                             \
      for(i=0;i<length;i++){                                                    \
        (mumpspointer[i]).i = ptr_matlab[i];                                    \
      }                                                                         \
    }else{                                                                      \
      for(i=0;i<length;i++){                                                    \
        (mumpspointer[i]).i = 0.0;                                              \
      }                                                                         \
    }                                                                           \
  }


#define EXTRACT_CMPLX_FROM_C_TO_MATLAB(mxcomponent,mumpspointer,length)         \
  if(mumpspointer == 0){                                                        \
    mxcomponent = mxCreateDoubleMatrix (1, 1, mxCOMPLEX);                       \
    ptr_matlab = mxGetPr (mxcomponent);                                         \
    ptr_matlab[0] = -9999;                                                      \
    ptr_matlab = mxGetPi (mxcomponent);                                         \
    ptr_matlab[0] = -9999;                                                      \
  }else{                                                                        \
    mxcomponent = mxCreateDoubleMatrix (1,length,mxCOMPLEX);                    \
    ptr_matlab = mxGetPr (mxcomponent);                                         \
    for(i=0;i<length;i++){                                                      \
      ptr_matlab[i] = (mumpspointer[i]).r;                                      \
    }                                                                           \
    ptr_matlab = mxGetPi (mxcomponent);                                         \
    for(i=0;i<length;i++){                                                      \
      ptr_matlab[i] = (mumpspointer[i]).i;                                      \
    }                                                                           \
  }

#endif

void DMUMPS_free(DMUMPS_STRUC_C **dmumps_par){
  if(*dmumps_par){
  MYFREE( (*dmumps_par)->irn );
  MYFREE( (*dmumps_par)->jcn  );
  MYFREE( (*dmumps_par)->a );
  MYFREE( (*dmumps_par)->irn_loc );
  MYFREE( (*dmumps_par)->jcn_loc );
  MYFREE( (*dmumps_par)->a_loc );
  MYFREE( (*dmumps_par)->eltptr );
  MYFREE( (*dmumps_par)->eltvar );
  MYFREE( (*dmumps_par)->a_elt );
  MYFREE( (*dmumps_par)->perm_in );
  /* colsca/rowsca might have been allocated by
   * MUMPS but in that case the corresponding pointer
   * is already equal to 0 so line below will do nothing */
  MYFREE( (*dmumps_par)->colsca );
  MYFREE( (*dmumps_par)->rowsca  );
  MYFREE( (*dmumps_par)->pivnul_list );
  MYFREE( (*dmumps_par)->listvar_schur );
  MYFREE( (*dmumps_par)->sym_perm );
  MYFREE( (*dmumps_par)->uns_perm );
  MYFREE( (*dmumps_par)->irhs_ptr);
  MYFREE( (*dmumps_par)->irhs_sparse);
  MYFREE( (*dmumps_par)->rhs_sparse);
  MYFREE( (*dmumps_par)->rhs);
  MYFREE( (*dmumps_par)->redrhs);
  MYFREE(*dmumps_par);
  }
}

void DMUMPS_alloc(DMUMPS_STRUC_C **dmumps_par){

  MYMALLOC((*dmumps_par),1,DMUMPS_STRUC_C);
  (*dmumps_par)->irn  = NULL;
  (*dmumps_par)->jcn  = NULL;
  (*dmumps_par)->a  = NULL;
  (*dmumps_par)->irn_loc  = NULL;
  (*dmumps_par)->jcn_loc  = NULL;
  (*dmumps_par)->a_loc  = NULL;
  (*dmumps_par)->eltptr  = NULL;
  (*dmumps_par)->eltvar  = NULL;
  (*dmumps_par)->a_elt  = NULL;
  (*dmumps_par)->perm_in  = NULL;
  (*dmumps_par)->colsca  = NULL;
  (*dmumps_par)->rowsca  = NULL;
  (*dmumps_par)->rhs  = NULL;
  (*dmumps_par)->redrhs  = NULL;
  (*dmumps_par)->rhs_sparse = NULL;
  (*dmumps_par)->irhs_sparse = NULL;
  (*dmumps_par)->irhs_ptr = NULL;
  (*dmumps_par)->pivnul_list  = NULL;
  (*dmumps_par)->listvar_schur  = NULL;
  (*dmumps_par)->schur  = NULL;
  (*dmumps_par)->sym_perm  = NULL;
  (*dmumps_par)->uns_perm  = NULL;
}

void mexFunction(int nlhs, mxArray *plhs[ ],
                 int nrhs, const mxArray *prhs[ ]) { 
  
  mwSize i, j, pos; // Schleifenzähler und Positionen als mwSize (64-Bit)
  double *ptr_matlab;
#if MUMPS_ARITH == MUMPS_ARITH_z
  double *ptri_matlab;
#endif
  mwSize tmp_m, tmp_n;

  /* C pointer for input parameters */
  intptr_t inst_address; // Sicherer Typ für Pointer-Adressen auf 64-Bit-Systemen
  mwSize n, m, ne, netrue;
  MUMPS_INT job; // MUMPS-Steuervariablen müssen MUMPS_INT sein
  mwIndex *irn_in, *jcn_in;
  
  /* variable for multiple and sparse rhs */
  mwSize posrhs;
  mwSize nbrhs, ldrhs, nz_rhs;
  mwIndex *irhs_ptr, *irhs_sparse;
  double *rhs_sparse;
#if MUMPS_ARITH == MUMPS_ARITH_z
  double *im_rhs_sparse;
#endif

  DMUMPS_STRUC_C *dmumps_par;
  int dosolve = 0;
  int donullspace = 0;
  int doanalysis = 0;
  int dofactorize = 0;
  
  EXTRACT_FROM_MATLAB_TOVAL(JOB, job);

  doanalysis = (job == 1 || job == 4 || job == 6);
  dofactorize = (job == 2 || job == 4 || job == 5 || job == 6);
  dosolve = (job == 3 || job == 5 || job == 6);

  if(job == -1){
    DMUMPS_alloc(&dmumps_par);
    EXTRACT_FROM_MATLAB_TOVAL(SYM, dmumps_par->sym);
    dmumps_par->job = -1;
    dmumps_par->par = 1;
    dmumps_c(dmumps_par);
    dmumps_par->nz = -1;
    dmumps_par->nz_alloc = -1;
  }else{
    EXTRACT_FROM_MATLAB_TOVAL(INST, inst_address);
    // Direktes, sicheres Casting der gespeicherten Adresse auf die MUMPS-Struktur
    dmumps_par = (DMUMPS_STRUC_C *) inst_address;

    if(job == -2){
      dmumps_par->job = -2;
      dmumps_c(dmumps_par);
      /* If colsca/rowsca were freed by MUMPS,
         dmumps_par->colsca/rowsca are now null.
         Application of MYFREE in call below thus ok */
      DMUMPS_free(&dmumps_par);
    }else{

      /* check of input arguments */
      n = mxGetN(A_IN);
      m = mxGetM(A_IN);

      if (!mxIsSparse(A_IN) || n != m )
          mexErrMsgTxt("Input matrix must be a sparse square matrix");
      
      jcn_in = mxGetJc(A_IN);
      ne = jcn_in[n];
      irn_in = mxGetIr(A_IN);
      dmumps_par->n = (MUMPS_INT)n;
      
      if(dmumps_par->sym != 0)
          netrue = (n+ne)/2;
      else
          netrue = ne;
      
      if(dmumps_par->nz_alloc < (MUMPS_INT)netrue || dmumps_par->nz_alloc >= (MUMPS_INT)(2*netrue)){  
        MYFREE(dmumps_par->jcn);
        MYFREE(dmumps_par->irn);
        MYFREE(dmumps_par->a);
        MYMALLOC((dmumps_par->jcn), netrue, MUMPS_INT); // Nutzt MUMPS_INT (64-Bit)
        MYMALLOC((dmumps_par->irn), netrue, MUMPS_INT); // Nutzt MUMPS_INT (64-Bit)
        MYMALLOC((dmumps_par->a), netrue, double2);
        dmumps_par->nz_alloc = (MUMPS_INT)netrue;
      }

      if(dmumps_par->sym == 0){
        /* if analysis already performed then we only need to read numerical values */
        if(doanalysis){ 
          for(i=0; i<n; i++){
            for(j=jcn_in[i]; j<jcn_in[i+1]; j++){
              (dmumps_par->jcn)[j] = (MUMPS_INT)(i+1);
              (dmumps_par->irn)[j] = (MUMPS_INT)(irn_in[j]+1);
            }
          }
        }
        dmumps_par->nz = (MUMPS_INT)ne;

#if MUMPS_ARITH == MUMPS_ARITH_z
        ptr_matlab = mxGetPr(A_IN);
        for(i=0; i<(mwSize)dmumps_par->nz; i++){                                                   
          ((dmumps_par->a)[i]).r = ptr_matlab[i];
        }
        ptr_matlab = mxGetPi(A_IN);
        if(ptr_matlab){
          for(i=0; i<(mwSize)dmumps_par->nz; i++){                                                   
            ((dmumps_par->a)[i]).i = ptr_matlab[i];
          }
        }else{
          for(i=0; i<(mwSize)dmumps_par->nz; i++){                                                   
             ((dmumps_par->a)[i]).i = 0.0;
             }
        }
#else
        ptr_matlab = mxGetPr(A_IN);
        for(i=0; i<(mwSize)dmumps_par->nz; i++){                                                   
          (dmumps_par->a)[i] = ptr_matlab[i];
        }
#endif
      }else{
        /* in the symmetric case we do not need to check doanalysis */
        pos = 0;
        ptr_matlab = mxGetPr(A_IN);
#if MUMPS_ARITH == MUMPS_ARITH_z
        ptri_matlab = mxGetPi(A_IN);
#endif
        for(i=0; i<n; i++){
          for(j=jcn_in[i]; j<jcn_in[i+1]; j++){
            if(irn_in[j] >= i){
              if(pos >= netrue)
                mexErrMsgTxt("Input matrix must be symmetric");
              (dmumps_par->jcn)[pos] = (MUMPS_INT)(i+1);
              (dmumps_par->irn)[pos] = (MUMPS_INT)(irn_in[j]+1);
#if MUMPS_ARITH == MUMPS_ARITH_z
              ((dmumps_par->a)[pos]).r = ptr_matlab[j];
              if(ptri_matlab){
                ((dmumps_par->a)[pos]).i = ptri_matlab[j];
              }else{
                ((dmumps_par->a)[pos]).i = 0.0;
              }
#else
              (dmumps_par->a)[pos] = ptr_matlab[j];
#endif
              pos++;
             }
          }
        }
        dmumps_par->nz = (MUMPS_INT)pos;
      }
    
      EXTRACT_FROM_MATLAB_TOVAL(JOB, dmumps_par->job);
      EXTRACT_FROM_MATLAB_TOARR(ICNTL_IN, dmumps_par->icntl, MUMPS_INT, 60);
      EXTRACT_FROM_MATLAB_TOARR(CNTL_IN, dmumps_par->cntl, double, 15);
      EXTRACT_FROM_MATLAB_TOPTR(PERM_IN, (dmumps_par->perm_in), MUMPS_INT, n);

      /* colsca and rowsca are treated differently */
      EXTRACT_SCALING_FROM_MATLAB_TOPTR(COLSCA_IN, (dmumps_par->colsca), (dmumps_par->colsca_from_mumps), n);
      EXTRACT_SCALING_FROM_MATLAB_TOPTR(ROWSCA_IN, (dmumps_par->rowsca), (dmumps_par->rowsca_from_mumps), n);

      EXTRACT_FROM_MATLAB_TOARR(KEEP_IN, dmumps_par->keep, MUMPS_INT, 500);
      EXTRACT_FROM_MATLAB_TOARR(DKEEP_IN, dmumps_par->dkeep, double, 230);

      dmumps_par->size_schur = (MUMPS_INT)mxGetN(VAR_SCHUR);
      EXTRACT_FROM_MATLAB_TOPTR(VAR_SCHUR, (dmumps_par->listvar_schur), MUMPS_INT, dmumps_par->size_schur);
      if(!dmumps_par->listvar_schur) dmumps_par->size_schur = 0;

      ptr_matlab = mxGetPr (RHS_IN);
    if ( dmumps_par->icntl[25-1] == -1 && dmumps_par->infog[28-1] > 0 ) {
          dmumps_par->nrhs = dmumps_par->infog[28-1];
          donullspace = dosolve;
         }
      else if ( dmumps_par->icntl[25-1] > 0 && dmumps_par->icntl[25-1] <= dmumps_par->infog[28-1] ) {
           dmumps_par->nrhs = 1;
           donullspace = dosolve;
         }
      else {
           donullspace = 0;
         }
      if (donullspace) {
        nbrhs = (mwSize)dmumps_par->nrhs; ldrhs = n;
        dmumps_par->lrhs = (MUMPS_INT)n;
        MYMALLOC((dmumps_par->rhs), ((size_t)(dmumps_par->n) * (size_t)(dmumps_par->nrhs)), double2);
         }
      else if((!dosolve) || ptr_matlab == (double*)-9999 ) { /* rhs not already provided, or not used */
        EXTRACT_CMPLX_FROM_MATLAB_TOPTR(RHS_IN, (dmumps_par->rhs), double, 1);
      }else{
        nbrhs = mxGetN(RHS_IN);
        ldrhs = mxGetM(RHS_IN);
        dmumps_par->nrhs = (MUMPS_INT)nbrhs;
        dmumps_par->lrhs = (MUMPS_INT)ldrhs;
        if(ldrhs != n){
          mexErrMsgTxt ("Incompatible number of rows in RHS");
        }
        if (!mxIsSparse(RHS_IN)){ /* full rhs */
          dmumps_par->icntl[20-1] = 0;
          EXTRACT_CMPLX_FROM_MATLAB_TOPTR(RHS_IN, (dmumps_par->rhs), double, (mwSize)(dmumps_par->nrhs * dmumps_par->lrhs));
        }else{ /* sparse rhs */
          if (dmumps_par->icntl[30-1] == 0) {
            dmumps_par->icntl[20-1] = 1;
          }
          irhs_ptr = mxGetJc(RHS_IN);
          irhs_sparse = mxGetIr(RHS_IN);
          rhs_sparse = mxGetPr(RHS_IN);
#if MUMPS_ARITH == MUMPS_ARITH_z
          im_rhs_sparse = mxGetPi(RHS_IN);
#endif
          nz_rhs = irhs_ptr[nbrhs];
          dmumps_par->nz_rhs = (MUMPS_INT)nz_rhs;

          MYMALLOC((dmumps_par->irhs_ptr), (size_t)(dmumps_par->nrhs + 1), MUMPS_INT);
          MYMALLOC((dmumps_par->irhs_sparse), (size_t)dmumps_par->nz_rhs, MUMPS_INT);
          MYMALLOC((dmumps_par->rhs_sparse), (size_t)dmumps_par->nz_rhs, double2);
          MYMALLOC((dmumps_par->rhs), ((size_t)dmumps_par->nrhs * (size_t)dmumps_par->lrhs), double2);

          for(i=0; i<(mwSize)dmumps_par->nrhs; i++){
            for(j=irhs_ptr[i]; j<irhs_ptr[i+1]; j++){
              (dmumps_par->irhs_sparse)[j] = (MUMPS_INT)(irhs_sparse[j] + 1);
            }
            (dmumps_par->irhs_ptr)[i] = (MUMPS_INT)(irhs_ptr[i] + 1);
          }
          (dmumps_par->irhs_ptr)[dmumps_par->nrhs] = dmumps_par->nz_rhs + 1;
#if MUMPS_ARITH == MUMPS_ARITH_z
          if(im_rhs_sparse){
            for(i=0; i<(mwSize)dmumps_par->nz_rhs; i++){                                                   
              ((dmumps_par->rhs_sparse)[i]).r = rhs_sparse[i];
              ((dmumps_par->rhs_sparse)[i]).i = im_rhs_sparse[i];
            }
          }else{
            for(i=0; i<(mwSize)dmumps_par->nz_rhs; i++){                                                   
              ((dmumps_par->rhs_sparse)[i]).r = rhs_sparse[i];
              ((dmumps_par->rhs_sparse)[i]).i = 0.0;
            }
          }
#else
          for(i=0; i<(mwSize)dmumps_par->nz_rhs; i++){                                                   
            (dmumps_par->rhs_sparse)[i] = rhs_sparse[i];
          }
#endif
        }
      }

      if(dmumps_par->size_schur > 0){
        if (dofactorize) {
          MYMALLOC((dmumps_par->schur), ((size_t)(dmumps_par->size_schur) * (size_t)(dmumps_par->size_schur)), double2);
        }
        dmumps_par->icntl[19-1] = 1;
      }else{
        dmumps_par->icntl[19-1] = 0;
      }
       /* Reduced RHS */
       if ( dmumps_par->size_schur > 0 && dosolve ) {
          if ( dmumps_par->icntl[26-1] == 2 ) {
            /* REDRHS is on input */
            tmp_m = mxGetM(REDRHS_IN);
            tmp_n = mxGetN(REDRHS_IN);
            if (tmp_m != (mwSize)dmumps_par->size_schur || tmp_n != (mwSize)dmumps_par->nrhs) {
              mexErrMsgTxt ("bad dimensions for REDRHS in mumpsmex.c");
            }
            EXTRACT_CMPLX_FROM_MATLAB_TOPTR(REDRHS_IN, (dmumps_par->redrhs), double, (tmp_m * tmp_n));
            dmumps_par->lredrhs = dmumps_par->size_schur;
          }
          if ( dmumps_par->icntl[26-1] == 1 ) {
            /* REDRHS on output. Must be allocated before the call */
            MYFREE(dmumps_par->redrhs);
            if(!(dmumps_par->redrhs = (double2 *)malloc(((size_t)dmumps_par->size_schur * (size_t)dmumps_par->nrhs) * sizeof(double2)))){
              mexErrMsgTxt("malloc redrhs failed in intmumpsc.c");
            }
          }
       }
      dmumps_c(dmumps_par);
    }
  }
  if(nlhs > 0){
    EXTRACT_FROM_C_TO_MATLAB( INFO_OUT , (dmumps_par->infog), 80);
    EXTRACT_FROM_C_TO_MATLAB( RINFO_OUT, (dmumps_par->rinfog), 40);
    /* A-1 on output */
    if ( dmumps_par->icntl[30-1] != 0 && dosolve ) {
      RHS_OUT = mxCreateSparse((mwSize)dmumps_par->n, (mwSize)dmumps_par->n, (mwSize)dmumps_par->nz_rhs, mxREAL2);

      irhs_ptr = mxGetJc(RHS_OUT);
      irhs_sparse = mxGetIr(RHS_OUT);
      for(j=0; j<(mwSize)(dmumps_par->nrhs + 1); j++){
         irhs_ptr[j] = (mwIndex) ((dmumps_par->irhs_ptr)[j] - 1);
      }
      ptr_matlab = mxGetPr(RHS_OUT);
#if MUMPS_ARITH == MUMPS_ARITH_z
        ptri_matlab = mxGetPi(RHS_OUT);
#endif
      for(i=0; i<(mwSize)dmumps_par->nz_rhs; i++){
#if MUMPS_ARITH == MUMPS_ARITH_z
        /* complex arithmetic */
        ptr_matlab[i] = (dmumps_par->rhs_sparse)[i].r;
        ptri_matlab[i] = (dmumps_par->rhs_sparse)[i].i;
#else
        /* real arithmetic */
        ptr_matlab[i] = (dmumps_par->rhs_sparse)[i];
#endif
        irhs_sparse[i] = (mwIndex)((dmumps_par->irhs_sparse)[i] - 1);
      }

    }
    else if(dmumps_par->rhs && dosolve){
      nbrhs = (mwSize)dmumps_par->nrhs;
      RHS_OUT = mxCreateDoubleMatrix ((mwSize)dmumps_par->n, (mwSize)dmumps_par->nrhs, mxREAL2);
      ptr_matlab = mxGetPr (RHS_OUT);
#if MUMPS_ARITH == MUMPS_ARITH_z
      ptri_matlab = mxGetPi (RHS_OUT);
      for(j=0; j<(mwSize)dmumps_par->nrhs; j++){
        posrhs = j * n;
        for(i=0; i<(mwSize)dmumps_par->n; i++){
          ptr_matlab[posrhs+i] = (dmumps_par->rhs)[posrhs+i].r;
          ptri_matlab[posrhs+i] = (dmumps_par->rhs)[posrhs+i].i;
        }              
      }
#else
      for(j=0; j<(mwSize)dmumps_par->nrhs; j++){
        posrhs = j * (mwSize)dmumps_par->n;
        for(i=0; i<(mwSize)dmumps_par->n; i++){
          ptr_matlab[posrhs+i] = (dmumps_par->rhs)[posrhs+i];
        }              
      }
#endif
    }else{
      EXTRACT_CMPLX_FROM_C_TO_MATLAB( RHS_OUT, (dmumps_par->rhs), 1);
    }

    inst_address = (intptr_t) dmumps_par; // Direkte und sichere Speicherung als 64-Bit-Ganzzahl
    EXTRACT_FROM_C_TO_MATLAB( INST_OUT   , &inst_address, 1); 
    EXTRACT_FROM_C_TO_MATLAB( PIVNUL_LIST, dmumps_par->pivnul_list, (mwSize)dmumps_par->infog[27-1]);
    EXTRACT_FROM_C_TO_MATLAB( PERM_OUT   , dmumps_par->sym_perm, (mwSize)dmumps_par->n);
    EXTRACT_FROM_C_TO_MATLAB( UNS_PERM   , dmumps_par->uns_perm, (mwSize)dmumps_par->n);
    EXTRACT_FROM_C_TO_MATLAB( ICNTL_OUT  , dmumps_par->icntl, 60);
    EXTRACT_FROM_C_TO_MATLAB( CNTL_OUT   , dmumps_par->cntl, 15);
    EXTRACT_FROM_C_TO_MATLAB( ROWSCA_OUT , dmumps_par->rowsca, (mwSize)dmumps_par->n);
    EXTRACT_FROM_C_TO_MATLAB( COLSCA_OUT , dmumps_par->colsca, (mwSize)dmumps_par->n);
    EXTRACT_FROM_C_TO_MATLAB( KEEP_OUT   , dmumps_par->keep, 500);
    EXTRACT_FROM_C_TO_MATLAB( DKEEP_OUT  , dmumps_par->dkeep, 230);

    if(dmumps_par->size_schur > 0 && dofactorize){
      SCHUR_OUT = mxCreateDoubleMatrix((mwSize)dmumps_par->size_schur, (mwSize)dmumps_par->size_schur, mxREAL2);
      ptr_matlab = mxGetPr (SCHUR_OUT);
#if MUMPS_ARITH == MUMPS_ARITH_z
      ptri_matlab = mxGetPi (SCHUR_OUT);
      for(i=0; i<(mwSize)dmumps_par->size_schur; i++){
        pos = i * (mwSize)(dmumps_par->size_schur);
        for(j=0; j<(mwSize)dmumps_par->size_schur; j++){
          ptr_matlab[j+pos] = ((dmumps_par->schur)[j+pos]).r;
          ptri_matlab[j+pos] = ((dmumps_par->schur)[j+pos]).i;
        }
      }
#else
      for(i=0; i<(mwSize)dmumps_par->size_schur; i++){
        pos = i * (mwSize)(dmumps_par->size_schur);
        for(j=0; j<(mwSize)dmumps_par->size_schur; j++){
          ptr_matlab[j+pos] = (dmumps_par->schur)[j+pos];
        }
      }
#endif
    }else{
      SCHUR_OUT = mxCreateDoubleMatrix(1, 1, mxREAL2);
      ptr_matlab = mxGetPr (SCHUR_OUT);
      ptr_matlab[0] = -9999; 
#if MUMPS_ARITH == MUMPS_ARITH_z
      ptr_matlab = mxGetPi (SCHUR_OUT);
      ptr_matlab[0] = -9999;
#endif 
    }
    /* REDRHS on output */
    if ( dmumps_par->icntl[26-1]==1 && dmumps_par->size_schur > 0 && dosolve ) {
      REDRHS_OUT = mxCreateDoubleMatrix((mwSize)dmumps_par->size_schur, (mwSize)dmumps_par->nrhs, mxREAL2);
      ptr_matlab = mxGetPr(REDRHS_OUT);
#if MUMPS_ARITH == MUMPS_ARITH_z
      ptri_matlab = mxGetPi (REDRHS_OUT);
#endif
      for(i=0; i<(mwSize)(dmumps_par->nrhs * dmumps_par->size_schur); i++){
#if MUMPS_ARITH == MUMPS_ARITH_z
        ptr_matlab[i] = ((dmumps_par->redrhs)[i]).r;
        ptri_matlab[i] = ((dmumps_par->redrhs)[i]).i;
#else
        ptr_matlab[i] = ((dmumps_par->redrhs)[i]);
#endif
      }
    }else{
      REDRHS_OUT = mxCreateDoubleMatrix(1, 1, mxREAL2);
      ptr_matlab = mxGetPr (REDRHS_OUT);
      ptr_matlab[0] = -9999; 
#if MUMPS_ARITH == MUMPS_ARITH_z
      ptr_matlab = mxGetPi (REDRHS_OUT);
      ptr_matlab[0] = -9999;
#endif 
    }

    MYFREE(dmumps_par->redrhs);
    MYFREE(dmumps_par->schur);
    MYFREE(dmumps_par->irhs_ptr);
    MYFREE(dmumps_par->irhs_sparse);
    MYFREE(dmumps_par->rhs_sparse);
    MYFREE(dmumps_par->rhs);
  }
}
