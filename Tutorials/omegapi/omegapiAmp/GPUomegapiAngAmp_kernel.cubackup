//July 26th 2019, Based on DOI: 10.1016/0550-3213(84)90382-1


#include <stdio.h>

#include "GPUManager/GPUCustomTypes.h"
#include "GPUManager/CUDA-Complex.cuh"

#include "breakupMomentum.cuh"
#include "barrierFactor.cuh"

///////////////////////////////////////////////////
#include "cuda.h"

#include "GPUUtils/lorentzBoost.cuh"
#include "GPUUtils/threeVector.cuh"
#include "GPUUtils/wignerD.cuh"
#include "GPUUtils/clebsch.cuh"

///////////////////////////////////////////
#define ADD4(a,b) { a[0]+b[0], a[1]+b[1], a[2]+b[2], a[3]+b[3] }

#define MASS(v)   (G_SQRT(v[0]*v[0]-v[1]*v[1]-v[2]*v[2]-v[3]*v[3]))

#define Nterm(J)  (G_SQRT((2*J+1)/(4*M_PI)))


// Macro to ease definition of loops
#define LOOP(INDEX,START,END,INC) for (int INDEX=START;INDEX<=END;INDEX+=INC)
////////////////////////////////////////////////////////////////////////////////////////
 __device__ WCUComplex CZero = { 0, 0 };
 __device__ WCUComplex COne  = { 1, 0 };
 __device__ WCUComplex ic  = { 0, 1 };

 __device__ GDouble DegToRad = PI/180.0;
 __device__ GDouble pararray[22];

//Create array of lmLM:
 __device__ int lmLM[25][4] = {{0,0,0,0}, {0,0,2,0}, {0,0,2,1}, {0,0,2,2}, {2,0,0,0}, {2,0,2,0}, {2,0,2,1}, {2,0,2,2}, {2,1,2,0}, {2,1,2,1}, {2,1,2,2}, {2,2,2,0}, {2,2,2,1}, {2,2,2,2}, {2,1,1,1}, {0,0,1,0}, {0,0,1,1}, {2,1,1,0}, {2,1,1,1}, {2,1,2,1}, {2,1,2,2}, {2,2,2,1}, {2,2,2,2}, {2,0,1,0}, {2,0,1,1}};

__global__ void print_kernel() {
    printf("block %d, thread %d\n", blockIdx.x, threadIdx.x);
}

static __device__ void //note: 4-vector input presumed
rotateZ( GDouble* v, GDouble phi ){
  GDouble sinphi = G_SIN(phi);
  GDouble cosphi = G_COS(phi);
  GDouble tx;
  tx   = v[1] * cosphi - v[2] * sinphi;
  v[2] = v[2] * cosphi + v[1] * sinphi;
  v[1] = tx;
}

static __device__ void //note: 4-vector input presumed
rotateY ( GDouble* v, GDouble theta) {
  GDouble sinphi = G_SIN(theta);
  GDouble cosphi = G_COS(theta);
  GDouble tz;
  tz = v[3] * cosphi - v[1] * sinphi;
  v[1] = v[1] * cosphi + v[3] * sinphi;
  v[3] = tz;
}

static __device__ GDouble  //note: 3-vector input presumed
theta( GDouble* pv ){
  GDouble r= G_SQRT(pv[0]*pv[0] + pv[1]*pv[1]);
  return  G_ATAN2( r , pv[2] );
}


//static __device__ GDouble  //note: 3-vector input presumed
//phi( GDouble* pv ){
//  return  G_ATAN2( pv[1] , pv[0] );
//}


static __device__ void
MoveToRF(GDouble *parent, GDouble *daughter)
{
  GDouble *par3vec=parent+1;
  rotateZ( daughter , -phi(par3vec) );
  rotateY( daughter , -theta(par3vec) );

  GDouble beta[]={0,0, -G_SQRT(dot(par3vec,par3vec))/parent[0]};
  //** (x)  Might this be bootToRest???
  // beta is defined to boost to parent's rest frame
  // I just adapted GPUUtil boost fcn with vector beta input
  boost( daughter , beta );

}

//////////////////////////////////////////////////////////////////////////////////////
static __device__ int delta(int first, int second)
{
  int kroneckerdelta = 0;

  if (first == second) kroneckerdelta = 1;
    
  //cout << "kroneker delta =" << kroneckerdelta << endl;
  return kroneckerdelta;
}

static __device__ GDouble calpha(int alpha)
{
  GDouble normalization = 0.0;
  
  int l = lmLM[alpha][0];
  int m = lmLM[alpha][1];
  int L = lmLM[alpha][2];
  int M = lmLM[alpha][3];
  
  GDouble normalization_num = 16.0 * PI * PI;
  GDouble normalization_denum = ((2.0 * l)+1.0) * ((2.0 * L)+1.0) * (2.0-delta(m,0)) * (2.0-delta(M,0));
  normalization = normalization_denum/normalization_num;
  
  return normalization;
}

static __device__ GDouble hmoment(int alpha, GDouble* vector)
{
  GDouble loccostheta = G_COS(vector[0]);
  GDouble locphi = vector[1];
  GDouble loccosthetaH = G_COS(vector[2]);
  GDouble locphiH = vector[3];

  int l = lmLM[alpha][0];
  int m = lmLM[alpha][1];
  int L = lmLM[alpha][2];
  int M = lmLM[alpha][3];

   WCUComplex moment = {0.0, 0.0};
  if (alpha < 15)
  {
    moment = 0.5 * wignerD(L, M, m, loccostheta, locphi) * wignerD(l, m, 0, loccosthetaH, locphiH) + pow(-1.0,L+M) * wignerD(L, -M, m, loccostheta, locphi) * wignerD(l, m, 0, loccosthetaH, locphiH);
  }
  else {moment = 0.5 * wignerD(L, M, m, loccostheta, locphi) * wignerD(l, m, 0, loccosthetaH, locphiH) - pow(-1.0,L+M) * wignerD(L, -M, m, loccostheta, locphi) * wignerD(l, m, 0, loccosthetaH, locphiH);}

//    GDouble real_moment = SqrtP(moment * Conjugate(moment));
  return moment.m_dRe;
}


//Define breakup momentum (x here is the measured mass of the omegapi)
static __device__ GDouble q1(GDouble x) {
  GDouble pi0mass = 0.1349766; //mass of the pi0 in GeV
  GDouble omegamass = 0.78265; //mass of the omega in GeV
  if (x < pi0mass + omegamass)
    return 0.0;
  return 0.5 * G_SQRT((G_POW(x, 2.0) - (pi0mass + omegamass)*(pi0mass + omegamass))*(G_POW(x, 2.0) - (omegamass - pi0mass)*(omegamass - pi0mass))) / x;
}

//Define barrier factor ratio
static __device__ GDouble barrierratio(GDouble x, GDouble resonancemass, int l) {
  return barrierFactor(q1(x), l) / barrierFactor(q1(resonancemass), l);
}

static __device__ GDouble Gamma_alpha(GDouble x, GDouble resonancewidth, GDouble resonancemass, int l) {
  return resonancewidth * q1(x) / q1(resonancemass) * G_POW(barrierratio(x, resonancemass, l), 2);
}

static __device__ WCUComplex D_alpha(GDouble x, GDouble resonancemass, GDouble resonancewidth, int l) {
  WCUComplex denom = {G_POW(x, 2) - G_POW(resonancemass, 2), -1.0*resonancemass * Gamma_alpha(x, resonancewidth, resonancemass, l) };
  return resonancewidth * resonancemass / denom;
}

//Define J (spin) for a given ialpha
static __device__ int J_spin(int ialpha) { //ialpha{0,1,2} -> JP{1+,1-,0-}
  if (ialpha == 2)
    return 0;
  return 1;
}

//Define a parity function
static __device__ int eta_parity(int ialpha) {
  if (ialpha == 0)
    return 1;
  return -1;
}

//Define sum over l with Clebsch- coefficients and barrier factors
static __device__ GDouble lsum(GDouble x, GDouble resonancemass, int ialpha, int lambda, GDouble DoverS) {
  GDouble c_alpha; //partial wave amplitudes
  GDouble lsum = 0;
  for (int l = 0; l < 3; l++) { 
    if (l == 1 && (ialpha == 1 || ialpha == 2))
      c_alpha = 1;
    else if (l == 2 && ialpha == 0)
      c_alpha = DoverS;
    else if (l == 0 && ialpha == 0)
      c_alpha = 1;
    else 
      c_alpha = 0; 
    lsum += G_SQRT((2.0*l + 1.0)/(2.0*J_spin(ialpha) + 1.0)) * clebsch(l, 0, 1, lambda, J_spin(ialpha), lambda) * c_alpha * barrierratio(x, resonancemass, l);
  }
  return lsum;
}

static __device__ WCUComplex F_alpha_lambda(GDouble x, GDouble resonancemass, GDouble resonancewidth, GDouble G_alpha, int ialpha, int lambda, GDouble DoverS, int l) {
  return D_alpha(x, resonancemass, resonancewidth, l) * G_alpha * lsum(x, resonancemass, ialpha, lambda, DoverS);
}

//Define sum over omega spin states
static __device__ WCUComplex f_Llm(GDouble x, GDouble resonancemass, GDouble resonancewidth, GDouble G_alpha, int ialpha, GDouble resonancemass2, GDouble resonancewidth2, GDouble G_beta, int ibeta, GDouble DoverS, int l, int m, int L) {
  WCUComplex fsum = CZero;
  for (int lambda = -1; lambda < 2; lambda++) {
    for (int lambdaprime = -1; lambdaprime < 2; lambdaprime++) {
      if (ibeta == 2 && lambdaprime != 0)
	continue;
      if (ialpha == 2 && lambda != 0)
	continue;
      fsum += F_alpha_lambda(x, resonancemass, resonancewidth, G_alpha, ialpha, lambda, DoverS, l) * Conjugate(F_alpha_lambda(x, resonancemass2, resonancewidth2, G_beta, ibeta, lambdaprime, DoverS, l)) * clebsch(J_spin(ibeta), lambdaprime, L, m, J_spin(ialpha), lambda) * clebsch(1, lambdaprime, l, m, 1, lambda);
    }
  }
  return fsum;
}


//Define complex "amplitudes"
static __device__ WCUComplex f_0(GDouble phi0, GDouble theta) { 
  return {G_SQRT(0.5) * G_COS(phi0) * G_COS(theta), G_SQRT(0.5) * G_SIN(phi0) * G_COS(theta)};
}
static __device__ WCUComplex f_plus(GDouble phiplus, GDouble theta, GDouble psi) {
  return {G_SQRT(0.5) * G_COS(phiplus) * G_SIN(theta) * G_COS(psi), G_SQRT(0.5) * G_SIN(phiplus) * G_SIN(theta) * G_COS(psi)};
}
static __device__ WCUComplex f_minus(GDouble phiminus, GDouble theta, GDouble psi) {
  return {G_SQRT(0.5) * G_COS(phiminus) * G_SIN(theta) * G_SIN(psi), G_SQRT(0.5) * G_SIN(phiminus) * G_SIN(theta) * G_SIN(psi)};
}

//Define a helicity function
static __device__ int Lambda_H(int iH) { // {0,1,2}->{0,+1,-1}
  if (iH == 0)
    return 0;
  else if (iH == 1) 
    return 1;
  else
    return -1;
}


//Define production density matrix
static __device__ WCUComplex rho(GDouble phi0, GDouble theta, GDouble phiplus, GDouble phiminus, GDouble psi, GDouble phi02, GDouble theta2, GDouble phiplus2, GDouble phiminus2, GDouble psi2, int iH, int iHPrime, int ialpha, int ibeta) {
  WCUComplex f_alpha = CZero;
  WCUComplex f_beta = CZero; 
  WCUComplex f_alpha_neg = CZero;
  WCUComplex f_beta_neg = CZero; 
  if (iH == 0) {
    f_alpha = f_0(phi0, theta);
    f_alpha_neg = f_0(phi0, theta);
  }
  else if (iH == 1) {
    f_alpha = f_plus(phiplus, theta, psi);
    f_alpha_neg = f_minus(phiminus, theta, psi);
  }
  else { 
    f_alpha = f_minus(phiminus, theta, psi);
    f_alpha_neg = f_plus(phiplus, theta, psi);
  }
  if (iHPrime == 0) {
    f_beta = f_0(phi02, theta2);
    f_beta_neg = f_0(phi02, theta2);
  }
  else if (iHPrime == 1) {
    f_beta = f_plus(phiplus2, theta2, psi2);
    f_beta_neg = f_minus(phiminus2, theta2, psi2);
  }
  else { 
    f_beta = f_minus(phiminus2, theta2, psi2);
    f_beta_neg = f_plus(phiplus2, theta2, psi2);
  }

  return f_alpha * Conjugate(f_beta) + eta_parity(ialpha) * eta_parity(ibeta) * G_POW(-1.0, J_spin(ialpha) - J_spin(ibeta)) * G_POW(-1.0, Lambda_H(iH) - Lambda_H(iHPrime)) * f_alpha_neg * Conjugate(f_beta_neg);

}


//Define sum over helicity states
static __device__ WCUComplex HelicitySum(GDouble phi0, GDouble theta, GDouble phiplus, GDouble phiminus, GDouble psi, GDouble phi02, GDouble theta2, GDouble phiplus2, GDouble phiminus2, GDouble psi2, int ialpha, int ibeta, int L, int M) {
  WCUComplex sumH = CZero;
  for (int iH = 0; iH < 3; iH++) {
    for (int iHPrime = 0; iHPrime < 3; iHPrime++) {
      if (ialpha == 2 && iH > 0)
	continue;
      if (ibeta == 2 && iHPrime > 0)
	continue;
      sumH += rho(phi0, theta, phiplus, phiminus, psi, phi02, theta2, phiplus2, phiminus2, psi2, iH, iHPrime, ialpha, ibeta) * clebsch(J_spin(ibeta), Lambda_H(iHPrime), L, M, J_spin(ialpha), Lambda_H(iH));
    }
  }
  
      //cout << "SumH = " << sumH << endl;
  return sumH;
}

//Define t_star (sum of production density matrix)
static __device__ WCUComplex t_star_LM(GDouble phi0, GDouble theta, GDouble phiplus, GDouble phiminus, GDouble psi, GDouble phi02, GDouble theta2, GDouble phiplus2, GDouble phiminus2, GDouble psi2, int ialpha, int ibeta, int L, int M) {

  WCUComplex t_star_LM_par = G_SQRT((2.0*J_spin(ibeta) + 1.0)/(2.0*J_spin(ialpha) + 1.0)) * HelicitySum(phi0, theta, phiplus, phiminus, psi, phi02, theta2, phiplus2, phiminus2, psi2, ialpha, ibeta, L, M);
  
      //cout << "T*(LM) = " << t_star_LM_par << endl;

  return t_star_LM_par;
}

static __device__ GDouble SingleIntensity(GDouble x, GDouble *pararray, int l, int m, int L, int M, int ialpha, int ibeta) {
//    GDouble single_intensity = 0.;
    WCUComplex single_intensity = CZero;
    GDouble phiplus_alpha = 0;
    GDouble phiminus_alpha = 0;
    GDouble psi_alpha = 0;
    GDouble phiplus_beta = 0;
    GDouble phiminus_beta = 0;
    GDouble psi_beta = 0;
    if (ialpha < 2) {
      phiplus_alpha = pararray[5*ialpha + 12];
      phiminus_alpha = pararray[5*ialpha + 13];
      psi_alpha = pararray[5*ialpha + 14];
    }
    if (ibeta < 2) {
      phiplus_beta = pararray[5*ibeta + 12];
      phiminus_beta = pararray[5*ibeta + 13];
      psi_beta = pararray[5*ibeta + 14];
    }

    single_intensity = t_star_LM(pararray[5*ialpha + 10], pararray[5*ialpha + 11], phiplus_alpha, phiminus_alpha, psi_alpha, pararray[5*ibeta + 10], pararray[5*ibeta + 11], phiplus_beta, phiminus_beta, psi_beta, ialpha, ibeta, L, M) * f_Llm(x, pararray[3*ialpha + 0], pararray[3*ialpha + 1], pararray[3*ialpha + 2], ialpha, pararray[3*ibeta + 0], pararray[3*ibeta + 1], pararray[3*ibeta + 2], ibeta, pararray[9], l, m, L) * clebsch(1, 0, l, 0, 1, 0);
    //cout << "single intensity = " << single_intensity << endl;

  return single_intensity.m_dRe;
}

static __device__ GDouble Intensity(GDouble x, GDouble *pararray, int i){
  
    GDouble IntensitySum = 0;
  for (int ialpha = 0; ialpha < 3; ialpha++) { //GDouble sum over combination of JP states
    for (int ibeta = 0; ibeta < 3; ibeta++) {
      if ( (ialpha != ibeta) & (eta_parity(ialpha) == eta_parity(ibeta)) ) //Only want interference moments with opposite parity
	continue;
      IntensitySum += SingleIntensity(x, pararray, lmLM[i][0], lmLM[i][1], lmLM[i][2], lmLM[i][3], ialpha, ibeta);
    }
  }
    //cout << "IntensitySum = " << IntensitySum << endl;

  return IntensitySum;
}

///////////////////////////////////////////////////////////////////////////////////
static __device__ GDouble wdist(int k, GDouble x, GDouble Phi, GDouble* vector)
  {
    GDouble dist = 0.0;

          for (int alpha = 0; alpha < 25; alpha++)//quantum states
	  {
	    if (k == 0){dist += Intensity(x, pararray, alpha) * hmoment(alpha, vector) * calpha(alpha);}

	    if (k == 1){dist += Intensity(x, pararray, alpha) * hmoment(alpha, vector) * calpha(alpha) * G_SQRT(2.0) * G_COS(2.0 * Phi);}

	    if (k == 2){dist +=  Intensity(x, pararray, alpha) * hmoment(alpha, vector) * calpha(alpha) * G_SQRT(2.0) * G_SIN(2.0 * Phi);}

	  }//alpha loop
  // cout << "wdist(" << k << ") = " << dist << endl;

    return dist;
  }

static __device__ WCUComplex sqrtIntensity(GDouble polfrac, GDouble b1mass, GDouble prod_angle, GDouble* vector){

	GDouble intensity = 0.0;
	GDouble real_sqrt_intensity = 0.0;
	GDouble ima_sqrt_intensity = 0.0;

	intensity = wdist(0, b1mass, prod_angle, vector) - polfrac * wdist(1, b1mass, prod_angle, vector) * G_COS(2.0 * prod_angle) -  polfrac * wdist(2, b1mass, prod_angle, vector) * G_SIN(2.0 * prod_angle);
	
       if (intensity >= 0.0) real_sqrt_intensity = G_SQRT(intensity);
       if (intensity < 0.0) ima_sqrt_intensity = G_SQRT(G_FABS(intensity));
	
	WCUComplex sqrt_intensity = {real_sqrt_intensity, ima_sqrt_intensity};
	//cout << "intensity = " << intensity << ", sqrt_intensity = " << sqrt_intensity << endl;
       return sqrt_intensity;
	//return COne;
}


/////////////////////////////////////////////////////////////////////////////////////////
__global__ void
GPUomegapiAngAmp_kernel( GPU_AMP_PROTO, GDouble m_1p, GDouble m_w_1p, GDouble m_n_1p, GDouble m_1m, GDouble m_w_1m, GDouble m_n_1m, GDouble m_0m, GDouble m_w_0m, GDouble m_n_0m, GDouble m_ds_ratio, GDouble m_phi0_1p,
                            GDouble m_theta_1p, GDouble m_phip_1p, GDouble m_phim_1p, GDouble m_psi_1p,
			    GDouble m_phi0_1m, GDouble m_theta_1m, GDouble m_phip_1m, GDouble m_phim_1m,
                            GDouble m_psi_1m, GDouble m_phi0_0m, GDouble m_theta_0m, bool useCutoff, GDouble polAngle, GDouble polFraction )
{

	int iEvent = GPU_THIS_EVENT;
        
    pararray[0] = m_1p;    pararray[1] = m_w_1p;    pararray[2] = m_n_1p;    pararray[3] = m_1m;
    pararray[4] = m_w_1m;    pararray[5] = m_n_1m;    pararray[6] = m_0m;    pararray[7] = m_w_0m;
    pararray[8] = m_n_0m;    pararray[9] = m_ds_ratio;    pararray[10] = m_phi0_1p;    pararray[11] = m_theta_1p;
    pararray[12] = m_phip_1p;    pararray[13] = m_phim_1p;    pararray[14] = m_psi_1p;    pararray[15] = m_phi0_1m;
    pararray[16] = m_theta_1m;    pararray[17] = m_phip_1m;    pararray[18] = m_phim_1m;    pararray[19] = m_psi_1m;
    pararray[20] = m_phi0_0m;    pararray[21] = m_theta_0m;

//////////////////////////////////////////// Get and create 4-vectors /////////////////////////////////////////
  GDouble  beam     [4] = GPU_P4(0);
  GDouble  target   [4] = {0.938, 0.0, 0.0, 0.0};
  GDouble  recoil   [4] = GPU_P4(1);
  GDouble  b1s_pi   [4] = GPU_P4(2);
  GDouble  omegas_pi[4] = GPU_P4(3);
  GDouble  rhos_pip [4] = GPU_P4(4);
  GDouble  rhos_pim [4] = GPU_P4(5);
/*
  GDouble  beam     [4] = {0.0, 0.0, 0.0, 3.0};
  GDouble  target   [4] = {0.938, 0.0, 0.0, 0.0};
  GDouble  recoil   [4] = GPU_P4(0);
  GDouble  b1s_pi   [4] = GPU_P4(1);
  GDouble  omegas_pi[4] = GPU_P4(2);
  GDouble  rhos_pip [4] = GPU_P4(3);
  GDouble  rhos_pim [4] = GPU_P4(4);
*/
  //  Make four-vector sums
  GDouble  rho   [4] = ADD4(rhos_pip, rhos_pim );
  GDouble  omega [4] = ADD4(rho,     omegas_pi);
  GDouble  b1    [4] = ADD4(omega,   b1s_pi);
  GDouble  gammap[4] = ADD4(beam,   target);

//////////////////////////////////////////// Calculate omega pi helicity angles /////////////////////////////////////////

  //omega angles in b1 decay the daughter = omega, parent = b1, InverseOfX = beam, rf = gammap
  MoveToRF(gammap,beam);
  MoveToRF(gammap,b1);
  MoveToRF(gammap,omega);    MoveToRF(b1,omega);

  GDouble z[3] = { b1[1], b1[2], b1[3] };
  makeUnit( z );
  GDouble y[3] = { beam[1], beam[2], beam[3] };
  cross( y, z );
  makeUnit( y );
  GDouble x[3] = { y[0], y[1], y[2] };
  cross( x, z );

  GDouble *omega_3vec=omega+1;
  GDouble ang_omega[]={dot(omega_3vec, x),
                    dot(omega_3vec, y),
                    dot(omega_3vec, z)};
  GDouble omega_theta = theta(ang_omega);
  GDouble omega_phi      = phi(ang_omega);
  
  GDouble eps[3] = {G_COS(polAngle*DegToRad), G_SIN(polAngle*DegToRad), 0.0}; // beam polarization vector
  GDouble eps_cross_y[3] = { eps[0], eps[1], eps[2] };
  cross( eps_cross_y, y );
  GDouble inverseofx[3] = { -x[0], -x[1], -x[2] };

  GDouble prod_angle = G_ATAN2(dot(y, eps), dot(inverseofx, eps_cross_y) );

////////////////////////////////////////////
  
//normal to the piplus+piminus plane angles in the b1 decay the daughter = piplus, parent = omega, InverseOfX = b1, rf = gammap, seconddaughter = piminus
  MoveToRF(gammap,rhos_pip); MoveToRF(b1,rhos_pip); MoveToRF(omega,rhos_pip);
  MoveToRF(gammap,rhos_pim); MoveToRF(b1,rhos_pim); MoveToRF(omega,rhos_pim);
  
  GDouble zh[3] = { omega[1], omega[2], omega[3] };
  makeUnit( zh );
  GDouble yh[3] = { b1[1], b1[2], b1[3] };
  cross( yh, zh );
  makeUnit( yh );
  GDouble xh[3] = { yh[0], yh[1], yh[2] };
  cross( xh, zh );

  GDouble normal_3vec[3] = { rhos_pip[1], rhos_pip[2], rhos_pip[3] };
  cross( normal_3vec, rhos_pim );
  makeUnit( normal_3vec );

  GDouble ang_normal[]={dot(normal_3vec, xh),
                    dot(normal_3vec, yh),
                    dot(normal_3vec, zh)};
  GDouble normal_theta = theta(ang_normal);
  GDouble normal_phi   = phi(ang_normal);
/////////////////////////////////////////////////////////////////   

  GDouble b1mass = MASS(b1);

  //GDouble q = breakupMomentum( MASS(b1), MASS(omega), MASS(b1s_pi) );

  //GDouble alpha = phi( &(recoil[1]) );

  //GDouble cosAlpha=G_COS(alpha), sinAlpha=G_SIN(alpha);

  GDouble angles[4] = {omega_theta, omega_phi, normal_theta, normal_phi};
  ////////////////////////////////////////////////////////////////////////////////////        
    //print_kernel<<<10, 10>>>();
    //cudaDeviceSynchronize();

  pcDevAmp[iEvent] = sqrtIntensity(polFraction, b1mass, prod_angle, angles);
  
//  pcDevAmp[iEvent] = COne;
}


void
GPUomegapiAngAmp_exec( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO, 
                     GDouble m_1p, GDouble m_w_1p, GDouble m_n_1p, GDouble m_1m, GDouble m_w_1m, GDouble m_n_1m,
                            GDouble m_0m, GDouble m_w_0m, GDouble m_n_0m, GDouble m_ds_ratio, GDouble m_phi0_1p,
                            GDouble m_theta_1p, GDouble m_phip_1p, GDouble m_phim_1p, GDouble m_psi_1p,
			    GDouble m_phi0_1m, GDouble m_theta_1m, GDouble m_phip_1m, GDouble m_phim_1m,
                            GDouble m_psi_1m, GDouble m_phi0_0m, GDouble m_theta_0m, bool useCutoff, GDouble polAngle, GDouble polFraction )
{  

  GPUomegapiAngAmp_kernel<<< dimGrid, dimBlock >>>
    ( GPU_AMP_ARGS, m_1p, m_w_1p, m_n_1p, m_1m, m_w_1m, m_n_1m,
                            m_0m, m_w_0m, m_n_0m, m_ds_ratio, m_phi0_1p,
                            m_theta_1p, m_phip_1p, m_phim_1p, m_psi_1p,
			    m_phi0_1m, m_theta_1m, m_phip_1m, m_phim_1m,
                            m_psi_1m, m_phi0_0m, m_theta_0m, useCutoff, polAngle, polFraction );
}
