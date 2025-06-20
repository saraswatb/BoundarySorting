#include "header.hpp"
#include "models/TwoFluidWetting.hpp"
#include "error_msg.hpp"
#include "random.hpp"
#include "lb.hpp"
#include "tools.hpp"

using namespace std;
namespace opt = boost::program_options;

// from main.cpp:
extern unsigned nthreads, nsubsteps;
extern double time_step;

TwoFluidWetting::TwoFluidWetting(unsigned LX, unsigned LY, unsigned BC)
  : Model(LX, LY, BC, BC==0 ? GridType::Periodic : GridType::Layer)
{}
TwoFluidWetting::TwoFluidWetting(unsigned LX, unsigned LY, unsigned BC, GridType Type)
  : Model(LX, LY, BC, Type)
{}

void TwoFluidWetting::Initialize()
{
  angle = angle_deg*M_PI/180.;
  angle2 = angle_deg2*M_PI/180.;
  theta_anchor = theta_anchor_deg*M_PI/180;
  eta = 1./3.*(tau1- .5)*rho1;

  Q_fluct = (Q_kBT!=0);
  thermoAnchor = (L_thermo!=0);  
  wallAnchor = (W_anchor!=0);

  ff.SetSize(LX, LY, Type);  
  ff_tmp.SetSize(LX, LY, Type);  
  gg1.SetSize(LX, LY, Type);  
  gg1_tmp.SetSize(LX, LY, Type);
  gg2.SetSize(LX, LY, Type);  
  gg2_tmp.SetSize(LX, LY, Type);

  QQxx1.SetSize(LX, LY, Type);
  QQyx1.SetSize(LX, LY, Type);
  QQxx2.SetSize(LX, LY, Type);
  QQyx2.SetSize(LX, LY, Type);

  n1.SetSize(LX, LY, Type);
  n2.SetSize(LX, LY, Type);
  nC.SetSize(LX, LY, Type);

  phi.SetSize(LX, LY, Type);
  phi_old.SetSize(LX, LY, Type);
  dxPPhi.SetSize(LX, LY, Type);
  dyPPhi.SetSize(LX, LY, Type);

  ux1.SetSize(LX, LY, Type);
  uy1.SetSize(LX, LY, Type);
  ux2.SetSize(LX, LY, Type);
  uy2.SetSize(LX, LY, Type);  
  uxC.SetSize(LX, LY, Type);
  uyC.SetSize(LX, LY, Type);
  
  HHxx1.SetSize(LX, LY, Type);
  HHyx1.SetSize(LX, LY, Type);
  HHxx2.SetSize(LX, LY, Type);
  HHyx2.SetSize(LX, LY, Type);
  angDeg.SetSize(LX, LY, Type);

  dxQQxx1.SetSize(LX, LY, Type);
  dyQQxx1.SetSize(LX, LY, Type);
  dxQQyx1.SetSize(LX, LY, Type);
  dyQQyx1.SetSize(LX, LY, Type);
  sigmaXX1.SetSize(LX, LY, Type);
  sigmaYY1.SetSize(LX, LY, Type);
  sigmaYX1.SetSize(LX, LY, Type);
  sigmaXY1.SetSize(LX, LY, Type);
  dxQQxx2.SetSize(LX, LY, Type);
  dyQQxx2.SetSize(LX, LY, Type);
  dxQQyx2.SetSize(LX, LY, Type);
  dyQQyx2.SetSize(LX, LY, Type);
  sigmaXX2.SetSize(LX, LY, Type);
  sigmaYY2.SetSize(LX, LY, Type);
  sigmaYX2.SetSize(LX, LY, Type);
  sigmaXY2.SetSize(LX, LY, Type);

  subForceX.SetSize(LX, LY, Type);
  subForceY.SetSize(LX, LY, Type);
  BulkFE_LDG_12.SetSize(LX, LY, Type);
  FFx1.SetSize(LX, LY, Type);
  FFy1.SetSize(LX, LY, Type);
  FFx2.SetSize(LX, LY, Type);
  FFy2.SetSize(LX, LY, Type);
  
  if (thermoAnchor){
    thermoAnchor_delMu.SetSize(LX, LY, Type);
  }
  
  ux_old1.SetSize(LX, LY, Type);
  uy_old1.SetSize(LX, LY, Type);
  advX1.SetSize(LX, LY, Type);
  advY1.SetSize(LX, LY, Type);  
  ux_old2.SetSize(LX, LY, Type);
  uy_old2.SetSize(LX, LY, Type);
  advX2.SetSize(LX, LY, Type);
  advY2.SetSize(LX, LY, Type);  

  if(nsubsteps>1)
    throw error_msg("time stepping not implemented for this model"
                    ", please set nsubsteps=1.");

  
  if (dropletMode){
    bigPhi.SetSize(LX, LY, Type); 
    bigPhi_tmp.SetSize(LX, LY, Type); 
    bigMu.SetSize(LX, LY, Type); 
  }
   

  set_seed(init_seed);
}

void TwoFluidWetting::ConfigureAtNode(unsigned k)
{
  double nematicOrder1 = 1;
  double nematicOrder2 = 1;
  double theta1, theta2;
  double conc1 = 1, conc2 = 1;
  const unsigned x = GetXPosition(k);
  const unsigned y = GetYPosition(k);
  double dubLX = LX, dubLY = LY, xtemp=x-(dubLX/2.0-0.5), ytemp=y-(dubLY/2.0-0.5);  
  
  //double ux1_tilde, uy1_tilde;
  if (dropletMode){
    bigPhi[k] = 1;
  }        

  //different mould configurations
  if     (init_config == "UnifUnif")
  {
    conc1 = 1;
    conc2 = 1;
  }
  else if(init_config == "Circle" || init_config == "LGPS_Circle" || init_config == "Aster" || init_config == "Vortex"){ 
    double rad_test = sqrt(xtemp*xtemp+ytemp*ytemp);
    if (rad_test  > radius ){
      nematicOrder1 = 1.0;
      nematicOrder2 = 1.0;
      bigPhi[k] = 0;      
    }
    else {
      bigPhi[k] = 1;
      nematicOrder1 = 1.0;
      nematicOrder2 = 1.0;
    }
  }
  else if(init_config == "Rings"){ 
    double rad_test = sqrt(xtemp*xtemp+ytemp*ytemp);
    if (rad_test  > radius ){
      nematicOrder1 = 1.0;
      nematicOrder2 = 1.0;
      bigPhi[k] = 0;      
    }
    else {
      bigPhi[k] = 1;
      nematicOrder1 = 1.0;
      nematicOrder2 = 1.0;
      if (rad_test > 0.7*radius){ //the rings
        conc1 = 0.5;
        conc2 = 1.5;
      }
      else{
        conc1 = 1.5;
        conc2 = 0.5;
      }
    }
  }
  else if (init_config == "Channel"){
    if (fabs(ytemp)<radius){
      bigPhi[k] = 1;
    }
    else{
      bigPhi[k] = 0;
    }
  }
  else if(init_config == "Stripes"){
      bigPhi[k] = 1;
  }
  else if (init_config == "DogicTrack"){ //rectangle with end caps
    nematicOrder1 = 1.0;
    nematicOrder2 = 1.0;
    bigPhi[k] = 0;
    if (fabs(xtemp)<=level && fabs(ytemp)<=radius){ //horizontal track
      bigPhi[k] = 1;
    } 
    else if (xtemp>0 && ((xtemp-level)*(xtemp-level) + ytemp*ytemp <= radius*radius) ){ //right end cap
      bigPhi[k] = 1;
    }
    else if (xtemp<0 && ((xtemp+level)*(xtemp+level) + ytemp*ytemp <= radius*radius) ){ //right end cap
      bigPhi[k] = 1;
    }
  }
  else if (init_config == "FuzzySquare"){ //rectangle with end caps
    nematicOrder1 = 1.0;
    nematicOrder2 = 1.0;
    bigPhi[k] = 0;
    if (fabs(xtemp)<=radius && fabs(ytemp)<=radius){ //horizontal track
      bigPhi[k] = 1;
    } 
  }
  else
    throw error_msg("error: initial configuration '", init_config, "' unknown.");

  theta1   = bigPhi[k]*noise*M_PI*(random_real() - .5) + angle; 
  theta2   = angle+ bigPhi[k]*noise*M_PI*(random_real() - .5);
  
  //defining theta fields for different nematic configuations
  if (init_config == "Circle" || init_config == "Anticircle"){    
    theta1 = angle + atan2(ytemp, xtemp) + M_PI/2  + pow(bigPhi[k], 2)*noise*M_PI*(random_real() - .5);
    theta2 = angle + atan2(ytemp, xtemp) + M_PI/2  + pow(bigPhi[k], 2)*noise*M_PI*(random_real() - .5);
  }
  else if(init_config == "Ellipse"){
    double xtemptemp = xtemp/(1+level)/(1+level);
    theta1 = atan2(-xtemptemp, ytemp) + noise*M_PI*(random_real() - .5);
    theta2 = atan2(-xtemptemp, ytemp) + noise*M_PI*(random_real() - .5);
  }
  else if(init_config == "EllipCirc"){ 
    //this configuration is made for the droplet subclass of this code
    if (!dropletMode) return;
    else{
      if (xtemp<0){ //circle side
        theta1 = atan2(ytemp, xtemp) + M_PI/2  + noise*M_PI*(random_real() - .5);
        theta2 = atan2(ytemp, xtemp) + M_PI/2  + noise*M_PI*(random_real() - .5);
      }
      if (xtemp>0){ //ellipse side
        double xtemptemp = xtemp/(1+level)/(1+level);
        theta1 = atan2(-xtemptemp, ytemp) + noise*M_PI*(random_real() - .5);
        theta2 = atan2(-xtemptemp, ytemp) + noise*M_PI*(random_real() - .5);
      }  
      
            
    }

  }
  else if (init_config == "Vortex"){    
    theta1 = atan2(ytemp, xtemp) + M_PI/2  + noise*M_PI*(random_real() - .5);
    theta2 = atan2(ytemp, xtemp) + M_PI/2  + noise*M_PI*(random_real() - .5);
  }
  else if (init_config == "Aster"){
    theta1 = atan2(ytemp, xtemp)+ noise*M_PI*(random_real() - .5);
    theta2 = theta1 +   noise*M_PI*(random_real() - .5);
  }
  /*
  if (fix_Q1){ //impose a fixed Q for debugging
    theta1 = sin(2*M_PI*x*level/(LX*1.0)); //sin waves
  }
  */

  QQxx1[k] = (cos(2*theta1))*nematicOrder1*bigPhi[k];
  QQyx1[k] = (sin(2*theta1))*nematicOrder1*bigPhi[k];
  QQxx2[k] = (cos(2*theta2))*nematicOrder2*bigPhi[k];
  QQyx2[k] = (sin(2*theta2))*nematicOrder2*bigPhi[k];

  ux1[k] = uy1[k] = 0;
  n1[k]  = rho1*(conc1+0.05*M_PI*(random_real() - .5)); 
  ux2[k] = uy2[k] = 0;
  n2[k]  = rho2*(conc2+0.05*M_PI*(random_real() - .5));
  if (fix_v1){//impose a fixed flow for debugging
    //uy1[k] = level; //constant flow
    uy1[k] = level*sin(2*M_PI*x/radius); //shear flow
    ux1[k] = 0;
  }

  nC[k] = n1[k]+n2[k];
  uxC[k] = (n1[k]*ux1[k]+n2[k]*ux2[k])/nC[k];
  uyC[k] = (n1[k]*uy1[k]+n2[k]*uy2[k])/nC[k];
  if (fix_vc){ //impose a fixed shear flow for debugging is BC!=0
    uyC[k] = level*sin(2*M_PI*x/radius);
    uxC[k] = 0;
  }

  phi[k] = n1[k]/(n1[k]+n2[k]);
  phi_old[k] = phi[k];  
  
  ff[k]  = GetEquilibriumDistribution(uxC[k], uyC[k], nC[k]); 
  gg1[k] = GetEquilibriumDistribution(ux1[k], uy1[k], nC[k]);
  gg2[k] = GetEquilibriumDistribution(ux1[k], uy1[k], n1[k]); 
  
  ftot1  = accumulate(begin(ff[k]), end(ff[k]), ftot1);
  ftot2  = accumulate(begin(gg2[k]), end(gg2[k]), ftot2);
  ntot1  += n1[k];
  ntot2  += n2[k];
}
void TwoFluidWetting::Configure()
{
  //print flags
  cout << "\n Flags: ";
  if (fix_v1) cout<<" v1 fixed | ";
  if (fix_vc) cout<<" vc fixed | ";
  if (fix_Q1) cout<<" Q1 fixed | ";
  if (fric_track_test) cout<<" Friction tracks ON | ";
  if (backflow_on) cout<<" backflow ON | ";
  if (wallAnchor) cout<<" Wall anchoring ON | "; 
  if (thermoAnchor) cout<<" TD anchoring ON | "; 
  if (fixedBigPhi) cout<<" Fixed BigPhi ON | "; 
  cout << "\n";

  setAngDeg();

  //initialize
  dt_FD = 1.0/nsteps_FD;
  for(unsigned k=0; k<DomainSize; ++k){
    ConfigureAtNode(k);    
  }

  //and make sure the boundary conditions work too
  BoundaryConditionsLB();
  BoundaryConditionsFields();
  BoundaryConditionsFields2();  
  
  //and relax bigPhi over some timesteps
  for (int i=0; i<nPP_init; i++){
    UpdateBigPhi();
  }

  //and do preinitialization
  cout << "Preinitialization started. ... ";
  preinit_flag = true;
  for (int i = 0; i< n_preinit; i++){
    Step();
  }
  preinit_flag = false;
  cout << " ... Preinitialization done. \n ";

}

LBNode TwoFluidWetting::CorrectLBAtNode(LBNode fb, unsigned k){
  //this correction makes sure the second moment of f has the anisotropic pressure  
  
  const double dxp = dxPPhi[k];
  const double dyp = dyPPhi[k];
  const double bigP = bigPhi[k];

  const double dxp2 = dxp*dxp;
  const double dyp2 = dyp*dyp;
  const double dxpdyp = dxp*dyp;

  fb[1] = fb[1] + 1./4.*KK*abs(1-bigP)*(dxp2-dyp2);
  fb[2] = fb[2] + 1./4.*KK*abs(1-bigP)*(dxp2-dyp2);
  fb[3] = fb[3] - 1./4.*KK*abs(1-bigP)*(dxp2-dyp2);
  fb[4] = fb[4] - 1./4.*KK*abs(1-bigP)*(dxp2-dyp2);
  fb[5] = fb[5] + 1./4.*KK*abs(1-bigP)*(dxpdyp);
  fb[6] = fb[6] + 1./4.*KK*abs(1-bigP)*(dxpdyp);
  fb[7] = fb[7] - 1./4.*KK*abs(1-bigP)*(dxpdyp);
  fb[8] = fb[8] - 1./4.*KK*abs(1-bigP)*(dxpdyp);
  // from mass conservation
  //fb[0] = nn-fe[1]-fe[2]-fe[3]-fe[4]-fe[5]-fe[6]-fe[7]-fe[8];

  return fb;
}

//big phi code
void TwoFluidWetting::UpdateBigPhiFieldAtNode(unsigned k){
  const auto& d = get_neighbours(k);  
  const double del2mu = laplacian(bigMu, d, sD);
  bigPhi_tmp[k] = bigPhi[k] + gammaPP*del2mu;
}
void TwoFluidWetting::UpdateBigPhiQntyAtNode(unsigned k){
  const auto& d = get_neighbours(k);  
  const double p = bigPhi[k];
  const double del2p = laplacian(bigPhi, d, sD);
  const double mu = APP*p*(1-p)*(1-2*p) - KPP*del2p;
  bigMu[k] = mu;
}
void TwoFluidWetting::UpdateBigPhi(){
  #pragma omp parallel for num_threads(nthreads) if(nthreads)
  for(unsigned k=0; k<DomainSize; k++) {    
    UpdateBigPhiQntyAtNode(k);
  }
  
  #pragma omp parallel for num_threads(nthreads) if(nthreads)
  for(unsigned k=0; k<DomainSize; k++) {    
    UpdateBigPhiFieldAtNode(k);
  }
  
  swap(bigPhi.get_data(), bigPhi_tmp.get_data());
}

void TwoFluidWetting::UpdateNemQuantitiesAtNode(unsigned k, bool first){
  const auto& d = get_neighbours(k);
  const double p = phi[k];
  const double bigP = bigPhi[k];
  double xtemp=GetXPosition(k)-(LX/2.0-0.5);
  double ytemp=GetYPosition(k)-(LY/2.0-0.5);
  
  const double NFE_mod = modNFE;//(bigPhi[k]>=0.10?1:0); //nematic free energy modifier   
  double zetaMod1 = (bigP>0.9?1:0); //p //bigP; //
  double zetaMod2 = (bigP>0.9?1:0); //(1-p)  
  double Snem_mod = abs(bigP);
  double CC1_eff = (bigP>BPcrit?CC1:0.0);
  double LL1_eff = (bigP>BPcrit?LL1:0.01);
  
  const double Qxx1 = QQxx1[k];
  const double Qyx1 = QQyx1[k];

  const double del2Qxx1  = laplacian(QQxx1,  d, sD);
  const double dxQxx1    = derivX   (QQxx1,  d, sB);
  const double dyQxx1    = derivY   (QQxx1,  d, sB);
  const double del2Qyx1  = laplacian(QQyx1,  d, sD);
  const double dxQyx1    = derivX   (QQyx1,  d, sB);
  const double dyQyx1    = derivY   (QQyx1,  d, sB);

  const double Qxx2 = QQxx2[k];
  const double Qyx2 = QQyx2[k];

  const double del2Qxx2  = laplacian(QQxx2,  d, sD);
  const double dxQxx2    = derivX   (QQxx2,  d, sB);
  const double dyQxx2    = derivY   (QQxx2,  d, sB);
  const double del2Qyx2  = laplacian(QQyx2,  d, sD);
  const double dxQyx2    = derivX   (QQyx2,  d, sB);
  const double dyQyx2    = derivY   (QQyx2,  d, sB);
  
  const double term1 = Snem1*Snem1*Snem_mod*Snem_mod*p*p - Qxx1*Qxx1 - Qyx1*Qyx1; //*p //*bigP*bigP
  double Hxx1 = NFE_mod*(CC1_eff*term1*Qxx1 + LL1_eff*del2Qxx1);
  double Hyx1 = NFE_mod*(CC1_eff*term1*Qyx1 + LL1_eff*del2Qyx1);
  
  //anchoring at wall terms
  double FE_anchoring = 0; bool wallEdge = false;
  if (wallAnchor){ //check if we are on a wall
    double ytemp = (double)(GetYPosition(k));
    double xtemp = (double)(GetYPosition(k));
    wallEdge = (BC>=1)? (ytemp==0 || ytemp == LY-1 || 
                    ((BC>=3)? (xtemp==0 || xtemp == LX-1) : false)) :  false;
    if (wallEdge){ //only preferentially anchor nematic 1
      double S_anchor_local = p>0 ? S_anchor*sqrt(p) : 0;
      const double Qxx_anchor = S_anchor_local*cos(2*theta_anchor);
      const double Qyx_anchor = S_anchor_local*sin(2*theta_anchor);

      const double diffXX = (Qxx1-Qxx_anchor);
      const double diffYX = (Qyx1-Qyx_anchor);
      Hxx1 = Hxx1 - W_anchor*diffXX;
      Hyx1 = Hyx1 - W_anchor*diffYX;

      if (p>0.05){
        FE_anchoring = 0.5*W_anchor*(diffXX*diffXX + diffYX*diffYX );//- p/(2*p+2)*(diffXX*Qxx_anchor + diffYX*Qyx_anchor));
      }
    } 
  }
  if (thermoAnchor){ //calculate and update thermodynamic anchoring terms for Hxx and Hxy    
    const double dxp = derivX(bigPhi, d, sB);
    const double dyp = derivY(bigPhi, d, sB);
    Hxx1 = Hxx1 + L_thermo*bigP*bigP*(dxp*dxp-dyp*dyp)/2;
    Hyx1 = Hyx1 + L_thermo*bigP*bigP*dxp*dyp;    
  }
  if (LL12 != 0 || K12_NR != 0){ //Q1 coupling to Q2
    //Hxx1 = Hxx1 - (LL12+K12_NR)*del2Qxx2;
    //Hyx1 = Hyx1 - (LL12+K12_NR)*del2Qyx2;
    Hxx1 = Hxx1 + (LL12+K12_NR)*Qxx2;
    Hyx1 = Hyx1 + (LL12+K12_NR)*Qyx2;

  }

  HHxx1[k]    =  Hxx1;
  HHyx1[k]    =  Hyx1;

  dxQQxx1[k]  =  dxQxx1;
  dxQQyx1[k]  =  dxQyx1;
  dyQQxx1[k]  =  dyQxx1;
  dyQQyx1[k]  =  dyQyx1;

  const double term2 = Snem2*Snem2*Snem_mod*Snem_mod - Qxx2*Qxx2 - Qyx2*Qyx2;  //*bigP*bigP
  double Hxx2 = NFE_mod*(CC2*term2*Qxx2 + LL2*del2Qxx2);
  double Hyx2 = NFE_mod*(CC2*term2*Qyx2 + LL2*del2Qyx2);
  
  if (LL12 != 0 || K21_NR != 0){ //Q2 coupling to Q1
    Hxx2 = Hxx2 + (LL12+K21_NR)*Qxx1;
    Hyx2 = Hyx2 + (LL12+K21_NR)*Qyx1;
  }

  HHxx2[k]    =  Hxx2;
  HHyx2[k]    =  Hyx2;

  dxQQxx2[k]  =  dxQxx2;
  dxQQyx2[k]  =  dxQyx2;
  dyQQxx2[k]  =  dyQxx2;
  dyQQyx2[k]  =  dyQyx2;

  if (!first)
    return;  

  //these quantities are calculated BEFORE the first nematic update step  
  const double sigmaB1 = (backflow_on? .5*CC1_eff*term1*term1 : 0) - (wallEdge ? FE_anchoring: 0); 
  const double sigmaF1 =  - (preinit_flag ? 0 : zeta1*zetaMod1*Qxx1) 
                          + (backflow_on? 2*xi1*( (Qxx1*Qxx1-1)*Hxx1 + Qxx1*Qyx1*Hyx1 ) : 0)                          
                          + (backflow_on? LL1_eff*(dyQxx1*dyQxx1+dyQyx1*dyQyx1-dxQxx1*dxQxx1-dxQyx1*dxQyx1) : 0)  ;                          
  const double sigmaS1 =  - (preinit_flag ? 0 : zeta1*zetaMod1*Qyx1)
                          + (backflow_on? 2*xi1*(Qyx1*Qxx1*Hxx1 + (Qyx1*Qyx1-1)*Hyx1) : 0)
                          - (backflow_on? 2*LL1_eff*(dxQxx1*dyQxx1+dxQyx1*dyQyx1) : 0) ;                                           
  const double sigmaA1 = backflow_on? 2*(Qxx1*Hyx1 - Qyx1*Hxx1) : 0; 

  sigmaXX1[k] =  sigmaF1 + sigmaB1;
  sigmaYY1[k] = -sigmaF1 + sigmaB1;
  sigmaXY1[k] =  sigmaS1 + sigmaA1;
  sigmaYX1[k] =  sigmaS1 - sigmaA1;
    
  if (thermoAnchor){
    //calculate and update thermodynamic anchoring terms for the stresses
    const double dxp = derivX(bigPhi, d, sB);
    const double dyp = derivY(bigPhi, d, sB);
    sigmaXX1[k] = sigmaXX1[k] + L_thermo*( dxp*dxp*Qxx1 + dxp*dyp*Qyx1);
    sigmaYY1[k] = sigmaYY1[k] + L_thermo*(-dyp*dyp*Qxx1 + dxp*dyp*Qyx1);
    sigmaXY1[k] = sigmaXY1[k] + L_thermo*(-dxp*dyp*Qxx1 + dxp*dxp*Qyx1);
    sigmaYX1[k] = sigmaYX1[k] + L_thermo*( dxp*dyp*Qxx1 + dyp*dyp*Qyx1);
  }  

  const double sigmaB2 = backflow_on? .5*CC2*term2*term2 : 0
                          - LL12*(Qxx1*Qxx2 + Qyx1*Qyx2);                          
  const double sigmaF2 =  - (preinit_flag ? 0 : zeta2*zetaMod2*Qxx2) 
                          + (backflow_on? 2*xi2*( (Qxx2*Qxx2-1)*Hxx2 + Qxx2*Qyx2*Hyx2 ) : 0)                          
                          + (backflow_on? LL2*(dyQxx2*dyQxx2+dyQyx2*dyQyx2-dxQxx2*dxQxx2-dxQyx2*dxQyx2) : 0)  ;                                    
  const double sigmaS2 =  - (preinit_flag? 0 : zeta2*zetaMod2*Qyx2)
                          + (backflow_on? 2*xi2*(Qyx2*Qxx2*Hxx2 + (Qyx2*Qyx2-1)*Hyx2) : 0)
                          - (backflow_on? 2*LL2*(dxQxx2*dyQxx2+dxQyx2*dyQyx2) : 0) ;                                    
  const double sigmaA2 = backflow_on? 2*(Qxx2*Hyx2 - Qyx2*Hxx2) : 0;  
  
  sigmaXX2[k] =  sigmaF2 + sigmaB2;
  sigmaYY2[k] = -sigmaF2 + sigmaB2;
  sigmaXY2[k] =  sigmaS2 + sigmaA2;
  sigmaYX2[k] =  sigmaS2 - sigmaA2;

  const double dFdphi = phi_FreeE? LL1*(dxQxx1*dxQxx1+dyQxx1*dyQxx1+dxQyx1*dxQyx1+dyQyx1*dyQyx1) + CC1*term1*term1 
                            - LL2*(dxQxx2*dxQxx2+dyQxx2*dyQxx2+dxQyx2*dxQyx2+dyQyx2*dyQyx2) - CC2*term2*term2 : 0;
  BulkFE_LDG_12[k] = (backflow_on? NFE_mod*dFdphi : 0) ;   
}
void TwoFluidWetting::UpdateNemQuantities(bool first){  
  // for reduction (+:sum)
  #pragma omp parallel for num_threads(nthreads) if(nthreads)
  for(unsigned k=0; k<DomainSize; k++)
  {
    UpdateNemQuantitiesAtNode(k, first);
  }
}

void TwoFluidWetting::UpdateFluidQuantitiesAtNode_p1(unsigned k)
{  
  //const auto& d = get_neighbours(k);
  const auto& f = ff[k];
  const auto& g2 = gg2[k];

  const double nnC = f[0] + f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7] + f[8];
  double nn1 = g2[0] + g2[1] + g2[2] + g2[3] + g2[4] + g2[5] + g2[6] + g2[7] + g2[8];
  
  const double nTOL = 0.08; // phi -> [nTOL, 1-nTOL]  
  if (nn1/nnC <= nTOL) nn1 = nTOL*nnC;     
  else if (nn1/nnC >= 1-nTOL) nn1 = (1-nTOL)*nnC;

  const double p  = nn1/nnC;
  const double nn2 = nnC-nn1;
  
  n1[k]       =  nn1;
  n2[k]       =  nn2;
  phi[k]      =  p;   
  nC[k]       =  nnC;  

  //times p factor dropped
  sigmaXX1[k] =  sigmaXX1[k]*p;
  sigmaYY1[k] =  sigmaYY1[k]*p;
  sigmaXY1[k] =  sigmaXY1[k]*p;
  sigmaYX1[k] =  sigmaYX1[k]*p;

  sigmaXX2[k] =  sigmaXX2[k]*(1-p);
  sigmaYY2[k] =  sigmaYY2[k]*(1-p);
  sigmaXY2[k] =  sigmaXY2[k]*(1-p);
  sigmaYX2[k] =  sigmaYX2[k]*(1-p);

}
void TwoFluidWetting::EvaluatePhiQuantities(unsigned k)
{
  const auto& d = get_neighbours(k); 
  const double dxPhi = derivX(phi, d, sB);
  const double dyPhi = derivY(phi, d, sB);

  dxPPhi[k] = dxPhi;
  dyPPhi[k] = dyPhi;
}
void TwoFluidWetting::UpdateFluidQuantitiesAtNode_p2(unsigned k)
{  
  const auto& d = get_neighbours(k); 
  const auto& f = ff[k];
  const auto& g1 = gg1[k];

  const double nnC = nC[k];
  const double p = phi[k];  
  const double bigP = bigPhi[k];

  const double FE_mod       = (bigP<BPcrit)? modFE:1;//0.2/landA:1; //1+abs(1-bigP)*20; //LG free energy modifier
  const double rel_fric_mod = (bigP<BPcrit)? modRF:1; //1+abs(1-bigP)*10; //relative drag modifier
  double fr_mod       = (bigP<BPcrit)? modFric:friction1; //1+abs(1-bigP)*10; //substrate friction modifier
  const double landA_mod = (init_config == "LGPS_Circle" || init_config == "Rings")? ((2*bigP-1)>0?1:-1):1;
  const bool isInterface = (bigP < BPcrit) && (bigP > 1-BPcrit);

  const double dxSxx1 = derivX(sigmaXX1, d, sB);
  const double dySxy1 = derivY(sigmaXY1, d, sB);
  const double dxSyx1 = derivX(sigmaYX1, d, sB);
  const double dySyy1 = derivY(sigmaYY1, d, sB);
  
  const double Fx1 = dxSxx1 + dySxy1;
  const double Fy1 = dxSyx1 + dySyy1;

  const double dxSxx2 = derivX(sigmaXX2, d, sB);
  const double dySxy2 = derivY(sigmaXY2, d, sB);
  const double dxSyx2 = derivX(sigmaYX2, d, sB);
  const double dySyy2 = derivY(sigmaYY2, d, sB);

  const double Fx2 = dxSxx2 + dySxy2;
  const double Fy2 = dxSyx2 + dySyy2;

  FFx1[k]     = Fx1;
  FFy1[k]     = Fy1;  
  FFx2[k]     = Fx2;
  FFy2[k]     = Fy2; 
  
  const double AA_LBF_C = isGuo? 0.5 : p*tau1 + (1-p)*tau2; 
  double vxC = fix_vc? uxC[k] : (f[1] - f[2] + f[5] - f[6] - f[7] + f[8] + AA_LBF_C*(Fx1+Fx2))/(nnC + AA_LBF_C*(p*fr_mod + (1-p)*fr_mod));
  double vyC = fix_vc? uyC[k] : (f[3] - f[4] + f[5] - f[6] + f[7] - f[8] + AA_LBF_C*(Fy1+Fy2))/(nnC + AA_LBF_C*(p*fr_mod + (1-p)*fr_mod));
  
  uxC[k]      =  vxC;
  uyC[k]      =  vyC;

  const double dxPhi = dxPPhi[k];
  const double dyPhi = dyPPhi[k];
  const double nn1 = n1[k];
  const double nn2 = n2[k];

  const double d2dxPhi = (KK!=0) ? laplacian(dxPPhi, d, sD) : 0;
  const double d2dyPhi = (KK!=0) ? laplacian(dyPPhi, d, sD) : 0;
 
  double phi_eqb = (init_config == "UnevenConcs") ? level : 0.5;
  double derLan_FreeEn = (2*landA*landA_mod+ 12*landB*(p-phi_eqb)*(p-phi_eqb));    //*(bigP<BPcrit?abs(landA)/landA:1)

  const double F_LG_x  = FE_mod*(nnC*derLan_FreeEn*dxPhi- KK*abs(1-bigP)*d2dxPhi)*(isInterface?0:1);    
  const double F_LG_y  = FE_mod*(nnC*derLan_FreeEn*dyPhi- KK*abs(1-bigP)*d2dyPhi)*(isInterface?0:1); 
  
  const double F_LdG_x = derivX(BulkFE_LDG_12, d, sB)*(isInterface?0:1);
  const double F_LdG_y = derivY(BulkFE_LDG_12, d, sB)*(isInterface?0:1);  

  const double sfX_stub =  rel_fric*rel_fric_mod*vxC  - (1-p)*(F_LG_x + F_LdG_x); 
  const double sfY_stub =  rel_fric*rel_fric_mod*vyC  - (1-p)*(F_LG_y + F_LdG_y);

  const double AA_LBF = isGuo? 0.5 : tau1;
  double vx1 = fix_v1? ux1[k] : ( g1[1]/phi_old[d[opp(1)]] - g1[2]/phi_old[d[opp(2)]] 
          + g1[5]/phi_old[d[opp(5)]] - g1[6]/phi_old[d[opp(6)]] - g1[7]/phi_old[d[opp(7)]] 
          + g1[8]/phi_old[d[opp(8)]]  +AA_LBF*(sfX_stub + FFx1[k]/p))/(nnC + AA_LBF*rel_fric*rel_fric_mod + AA_LBF*fr_mod);
  double vy1 = fix_v1? uy1[k] : ( g1[3]/phi_old[d[opp(3)]] - g1[4]/phi_old[d[opp(4)]] 
          + g1[5]/phi_old[d[opp(5)]] - g1[6]/phi_old[d[opp(6)]] + g1[7]/phi_old[d[opp(7)]] 
          - g1[8]/phi_old[d[opp(8)]]  +AA_LBF*(sfY_stub + FFy1[k]/p))/(nnC + AA_LBF*rel_fric*rel_fric_mod + AA_LBF*fr_mod);
  
  const double vx2 = (nnC*vxC-nn1*vx1)/nn2;
  const double vy2 = (nnC*vyC-nn1*vy1)/nn2;
  
  ux1[k]      =  vx1;
  uy1[k]      =  vy1;
  ux2[k]      =  vx2;
  uy2[k]      =  vy2;
  
  const double sfX = sfX_stub - rel_fric*vx1;
  const double sfY = sfY_stub - rel_fric*vy1;
  
  subForceX[k] = sfX;
  subForceY[k] = sfY;

  /* //TD Anchoring
  if (thermoAnchor){
    double d2xp  = derivXX(bigPhi, d, sB);
    double d2yp  = derivYY(bigPhi, d, sB);
    double dxdyp = derivXY(bigPhi, d, sB);

    thermoAnchor_delMu[k] = -L_thermo*(dxQQxx1[k]*dxPhi + dxQQyx1[k]*dyPhi + dyQQyx1[k]*dxPhi -dyQQxx1[k]*dyPhi 
                                + 2*dxQQyx1[k]*dxdyp + QQxx1[k]*(d2xp-d2yp)); 
  }
  */
}
void TwoFluidWetting::UpdateFluidQuantities()
{
  // for reduction (+:sum)
  #pragma omp parallel for num_threads(nthreads) if(nthreads)
  for(unsigned k=0; k<DomainSize; ++k)
  {
    UpdateFluidQuantitiesAtNode_p1(k);
  }

  #pragma omp parallel for num_threads(nthreads) if(nthreads)
  for(unsigned k=0; k<DomainSize; ++k)
  {
    EvaluatePhiQuantities(k);
  }

  #pragma omp parallel for num_threads(nthreads) if(nthreads)
  for(unsigned k=0; k<DomainSize; ++k)
  {
    UpdateFluidQuantitiesAtNode_p2(k);
  }

}

void TwoFluidWetting::UpdateNemFieldsAtNode(unsigned k, bool first){
  const auto& d = get_neighbours(k);  
  
  const double bigP = bigPhi[k];
  //if (bigPhi[k]<0.8) return;
  const double NFE_mod = (bigP<BPcrit)? modNFE:1; //1+abs(1-bigPhi[k])*10; //see also, L616

  const double vx1 = ux1[k];
  const double vy1 = uy1[k];
  const double vx2 = ux2[k];
  const double vy2 = uy2[k]; 

  const double dxux1  = derivX(ux1, d, sB);
  const double dyux1  = derivY(ux1, d, sB);
  const double dxuy1  = derivX(uy1, d, sB);
  const double dyuy1  = derivY(uy1, d, sB);

  const double dxux2  = derivX(ux2, d, sB);
  const double dyux2  = derivY(ux2, d, sB);
  const double dxuy2  = derivX(uy2, d, sB);
  const double dyuy2  = derivY(uy2, d, sB); 
  
  const double Qxx1   = QQxx1[k];
  const double Qyx1   = QQyx1[k];
  const double Hxx1   = HHxx1[k];
  const double Hyx1   = HHyx1[k];

  const double dxQxx1 = dxQQxx1[k];
  const double dyQxx1 = dyQQxx1[k];
  const double dxQyx1 = dxQQyx1[k];
  const double dyQyx1 = dyQQyx1[k];

  const double Qxx2   = QQxx2[k];
  const double Qyx2   = QQyx2[k];
  const double Hxx2   = HHxx2[k];
  const double Hyx2   = HHyx2[k];

  const double dxQxx2 = dxQQxx2[k];
  const double dyQxx2 = dyQQxx2[k];
  const double dxQyx2 = dxQQyx2[k];
  const double dyQyx2 = dyQQyx2[k]; 
  
  const double expansion1 = dxux1 + dyuy1;
  const double shear1     = .5*(dxuy1 + dyux1);
  const double vorticity1 = debugMode ? 0 : .5*(dxuy1 - dyux1);
  const double traceQL1   = Qxx1*(dxux1 - dyuy1) + 2*Qyx1*shear1;
  const double Dxx1       = NFE_mod*GammaQ1*Hxx1 - 2*vorticity1*Qyx1- 1*( vx1*dxQxx1 + vy1*dyQxx1 + Qxx1*expansion1 ) 
                            + xi1*((Qxx1+1)*(2*dxux1-traceQL1) +2*Qyx1*shear1 -expansion1); 
  const double Dyx1       = NFE_mod*GammaQ1*Hyx1 + 2*vorticity1*Qxx1- 1*( vx1*dxQyx1 + vy1*dyQyx1 + expansion1*Qyx1 )
                             + xi1*( Qyx1*(expansion1-traceQL1) + 2*shear1);

  const double expansion2 = dxux2 + dyuy2;
  const double shear2     = .5*(dxuy2 + dyux2);
  const double vorticity2 = debugMode ? 0 :  .5*(dxuy2 - dyux2);
  const double traceQL2   = Qxx2*(dxux2 - dyuy2) + 2*Qyx2*shear2;
  const double Dxx2       = NFE_mod*GammaQ2*Hxx2 - 2*vorticity2*Qyx2- 1*( vx2*dxQxx2 + vy2*dyQxx2 + expansion2*Qxx2 ) 
                            + xi2*((Qxx2+1)*(2*dxux2-traceQL2) +2*Qyx2*shear2 -expansion2); 
  const double Dyx2       = NFE_mod*GammaQ2*Hyx2 + 2*vorticity2*Qxx2- 1*( vx2*dxQyx2 + vy2*dyQyx2 + expansion2*Qyx2 )
                            + xi2*( Qyx2*(expansion2-traceQL2) + 2*shear2);

  if(first)
  {
    double Qxx1_noise = 0., Qyx1_noise = 0.;    

    if(Q_fluct)
    {
      // the noise hits S too hard, but we only want changes in theta
      /*
      double Q_stren = 2*Q_kBT*dt_FD*M_PI/180.0;
      Qxx1_noise = -Qxx1*(1-cos(Q_stren)) - Qyx1*(sin(Q_stren));
      Qyx1_noise = Qxx1*cos(Q_stren) - Qyx1*(1-sin(Q_stren));
      */
      double S1 = sqrt(Qxx1*Qxx1 + Qyx1*Qyx1);
      double theta = atan2(Qyx1, Qxx1);
      double theta_new = theta + dt_FD*(Q_kBT*M_PI/180.0)*(2*random_real()-1);
      Qxx1_noise = -Qxx1 + S1*cos(2*theta_new);
      Qyx1_noise = -Qyx1 + S1*sin(2*theta_new);
      
    }

    //QNxx1[k] = QNxx1[k] + 0.5*dt_FD*Dxx1;
    //QNyx1[k] = QNyx1[k] + 0.5*dt_FD*Dyx1;    
    //QNxx2[k] = QNxx2[k] + 0.5*dt_FD*Dxx2;
    //QNyx2[k] = QNyx2[k] + 0.5*dt_FD*Dyx2;
    if (!fix_Q1){
      QQxx1[k] = QQxx1[k] + dt_FD*Dxx1 + Qxx1_noise;
      QQyx1[k] = QQyx1[k] + dt_FD*Dyx1 + Qyx1_noise;    
    }    
    QQxx2[k] = QQxx2[k] + dt_FD*Dxx2;
    QQyx2[k] = QQyx2[k] + dt_FD*Dyx2;
  }
  else
  {
    //QQxx1[k] = QQxx1[k] + 0.5*dt_FD*Dxx1;
    //QQyx1[k] = QQyx1[k] + 0.5*dt_FD*Dyx1;    
    //QQxx2[k] = QQxx2[k] + 0.5*dt_FD*Dxx2;
    //QQyx2[k] = QQyx2[k] + 0.5*dt_FD*Dyx2;
  }

}
void TwoFluidWetting::UpdateNemFields(bool first){
  #pragma omp parallel for num_threads(nthreads) if(nthreads)
  for(unsigned k=0; k<DomainSize; ++k){    
    //if (dropletMode && bigPhi[k] <= BPcrit) continue; //k outside droplet
    UpdateNemFieldsAtNode(k, first);      
  }
}

void TwoFluidWetting::UpdateFluidFieldsAtNode(unsigned k, bool first)
{
  const auto& d = get_neighbours(k);
  const double bigP = bigPhi[k];

  double fr_mod = (bigP<BPcrit)? modFric:friction1; //1+abs(1-bigPhi[k])*10;  
  const double extra_visc_mod = (bigP<BPcrit)? extraVisc:0;//abs(1-bigPhi[k])*5;
  
  double vx1 = ux1[k]; 
  double vy1 = uy1[k];
  double vxC = uxC[k];
  double vyC = uyC[k];
  double vx2 = ux2[k];
  double vy2 = uy2[k];
  const double nn1 = n1[k];
  const double nnC = nC[k];
  const double p   = phi[k];
  const double p_toggle = phi_FreeE ? p : 0.5;  

  const double del2ux1 = laplacian(ux1, d, sD);
  const double del2uy1 = laplacian(uy1, d, sD);
  const double del2ux2 = laplacian(ux2, d, sD);
  const double del2uy2 = laplacian(uy2, d, sD);
  
  const double extra_viscX1 = extra_visc_mod*eta*del2ux1; 
  const double extra_viscY1 = extra_visc_mod*eta*del2uy1; 
  const double extra_viscX2 = extra_visc_mod*eta*del2ux2; 
  const double extra_viscY2 = extra_visc_mod*eta*del2uy2; 
  
  double FxC = FFx1[k] + FFx2[k] + extra_viscX1 + extra_viscX2
                      - p_toggle*fr_mod*vx1 - (1-p_toggle)*fr_mod*vx2; 
  double FyC = FFy1[k] + FFy2[k] + extra_viscY1 + extra_viscY2
                      - p_toggle*fr_mod*vy1 - (1-p_toggle)*fr_mod*vy2;                               

  double Gx = p_toggle*(subForceX[k] + (FFx1[k]+extra_viscX1)/p_toggle) - fr_mod*vx1;
  double Gy = p_toggle*(subForceY[k] + (FFy1[k]+extra_viscY1)/p_toggle) - fr_mod*vy1;
  
  //
  if (bigP < 0.75 && bigP > 0.25){ //boundary
    double angDegLoc = angDeg[k];
    double FtC = -sin(angDegLoc)*FxC + cos(angDegLoc)*FyC;
    double GtC = -sin(angDegLoc)*Gx + cos(angDegLoc)*Gy;
    FxC = -FtC*sin(angDegLoc); 
    FyC = FtC*cos(angDegLoc);
    Gx = -GtC*sin(angDegLoc); 
    Gy = GtC*cos(angDegLoc);

    double vtC = -vxC*sin(angDegLoc) + vyC*cos(angDegLoc);
    double vt1 = -vx1*sin(angDegLoc) + vy1*cos(angDegLoc);
    vxC = -vtC*sin(angDegLoc); vyC = vtC*cos(angDegLoc);
    vx1 = -vt1*sin(angDegLoc); vy1 = vt1*cos(angDegLoc);
  }
    //

  if (thermoAnchor){
    Gx = Gx + p_toggle*(1-p_toggle)*derivX(thermoAnchor_delMu, d, sB);
    Gy = Gy + p_toggle*(1-p_toggle)*derivY(thermoAnchor_delMu, d, sB);
  }

  double tau_mod = 1;//(bigPhi[k]>=0.5?1:100);
  double tau = (p_toggle*tau1+(1-p_toggle)*tau2)*tau_mod;
  
  const double dxn1    = derivX   (n1,  d, sB);
  const double dyn1    = derivY   (n1,  d, sB);
  const double vx1_tilde = vx1 + 1./6.* dxn1/nn1; 
  const double vy1_tilde = vy1 + 1./6.* dyn1/nn1; 
  
  // calculate the equilibrium distributions fe
  const auto fe  = CorrectLBAtNode(GetEquilibriumDistribution(vxC, vyC, nnC), k);
  const auto ge1 = CorrectLBAtNode(GetEquilibriumDistribution(vx1, vy1, nn1), k); 
  const auto ge2 = CorrectLBAtNode(GetEquilibriumDistribution(vx1_tilde, vy1_tilde, nn1), k);  

  
  if(first)
  {
    for(unsigned v=0; v<lbq; ++v)
    {
      double Si_C = isGuo? w[v]*(1-0.5/tau)*(FxC*xdir(v) + FyC*ydir(v))/cs2 : 0;      
      ff[k][v] = ff[k][v] +    (fe[v]-ff[k][v])/tau +    Si_C;

      double Si_1 = isGuo? w[v]*(1-0.5/tau1)*(Gx*xdir(v) + Gy*ydir(v))/cs2: 0;      
      //second order force correction : does not seem to work
      //Si_1 = isGuo? Si_1 + w[v]*(1-0.5/tau1)*(FO2_wtXX[v]*Gx*vx1 + FO2_wtXY[v]*Gx*vy1 + FO2_wtXY[v]*Gy*vx1 + FO2_wtYY[v]*Gy*vy1 ): 0;
      
      gg1[k][v] = gg1[k][v] +    (ge1[v]-gg1[k][v])/(tau1*tau_mod) +    Si_1;  
      gg2[k][v] = gg2[k][v] +    (ge2[v]-gg2[k][v])/(tau1*tau_mod);
    }
    
  }
  else
  { 
    //npc>1 not implemented   
    /* for(unsigned v=0; v<lbq; ++v){
      double Si_C = isGuo? w[v]*(1-0.5/tau)*(FxC*xdir(v) + FyC*ydir(v))/cs2 : 0;
      double Si_1 = isGuo? w[v]*(1-0.5/tau1)*(Gx*xdir(v) + Gy*ydir(v))/cs2: 0;
      //Si_1 = isGuo? Si_1 + FO2_wtXX[v]*Gx*vx1 + FO2_wtXY[v]*Gx*vy1 + FO2_wtXY[v]*Gy*vx1 + FO2_wtYY[v]*Gy*vy1 : 0;

      ff[k][v] = fn[k][v]    + .5*(fe[v]-ff[k][v])/tau    + .5*Si_C;                  
      gg1[k][v] = gn1[k][v]  + .5*(ge1[v]-gg1[k][v])/tau1 + .5*Si_1; ;
      
    }
    */
    
  }
}
void TwoFluidWetting::UpdateFluidFields(bool first)
{
  #pragma omp parallel for num_threads(nthreads) if(nthreads)
  for(unsigned k=0; k<DomainSize; ++k){
    //if (bigPhi[k] == 0) continue;
    UpdateFluidFieldsAtNode(k, first);      
  }

  /*
  //we have the phi for this timestep - save it in phi_old before streaming happens
  #pragma omp parallel for num_threads(nthreads) if(nthreads)
  for(unsigned k=0; k<DomainSize; ++k){
    phi_old[k]= phi[k];      
  }
  */
}

void TwoFluidWetting::Move()
{
  #pragma omp parallel for num_threads(nthreads) if(nthreads)
  for(unsigned k=0; k<TotalSize; ++k)
  {
    for(unsigned v=0; v<lbq; ++v)
    {
      // advect particles
      ff_tmp[next(k, v)][v] = ff[k][v];
      //fn_tmp[next(k, v)][v] = fn[k][v];
      
      gg1_tmp[next(k, v)][v] = gg1[k][v];
      //gn1_tmp[next(k, v)][v] = gn1[k][v]; 

      gg2_tmp[next(k, v)][v] = gg2[k][v];
      //gn2_tmp[next(k, v)][v] = gn2[k][v];      
    }
  }
  // swap temp variables
  swap(ff.get_data(), ff_tmp.get_data());
  //swap(fn.get_data(), fn_tmp.get_data());
  swap(gg1.get_data(), gg1_tmp.get_data());
  //swap(gn1.get_data(), gn1_tmp.get_data());

  swap(gg2.get_data(), gg2_tmp.get_data());
  //swap(gn2.get_data(), gn2_tmp.get_data());
}

void TwoFluidWetting::BoundaryConditionsLB()
{  
  switch(BC)
  {
    // pbc without bdry layer (nothing to do)
    case 0:
      break;
    // free-slip channel
    case 1:
      ff.ApplyFreeSlipChannel();
      //fn.ApplyFreeSlipChannel();
      gg1.ApplyFreeSlipChannel();
      //gn1.ApplyFreeSlipChannel();
      gg2.ApplyFreeSlipChannel();
      //gn2.ApplyFreeSlipChannel();
      break;
    // no-slip channel
    case 2:
      ff.ApplyNoSlipChannel();
      //fn.ApplyNoSlipChannel();
      gg1.ApplyNoSlipChannel();
      //gn1.ApplyNoSlipChannel();
      gg2.ApplyNoSlipChannel();
      //gn2.ApplyNoSlipChannel();
      break;
    // free-slip box
    case 3:
      ff.ApplyFreeSlip();
      //fn.ApplyFreeSlip();
      gg1.ApplyFreeSlip();
      //gn1.ApplyFreeSlip();
      gg2.ApplyFreeSlip();
      //gn2.ApplyFreeSlip();
      break;
    // no slip box
    case 4:
      ff.ApplyNoSlip();
      //fn.ApplyNoSlip();
      gg1.ApplyNoSlip();
      //gn1.ApplyNoSlip();
      gg2.ApplyNoSlip();
      //gn2.ApplyNoSlip();
      break;
    // pbc with boundary layer
    default:
      ff.ApplyPBC();
      //fn.ApplyPBC();
      gg1.ApplyPBC();
      //gn1.ApplyPBC();
      gg2.ApplyPBC();
      //gn2.ApplyPBC();
  }
}
void TwoFluidWetting::BoundaryConditionsFields()
{  
  const double anchorSign1 = (abs(theta_anchor) < 0.01) ? 1 : (abs(M_PI/2-theta_anchor)<0.01 ? -1 : 0);
  const double anchorSign2 = 0;//anchorSign1;
  switch(BC)
  {
    // pbc without bdry layer (nothing to do)
    case 0:
      break;
    // channel
    case 1:
    case 2:
      QQxx1.ApplyNeumannChannel();
      QQyx1.ApplyNeumannChannel();
      QQxx2.ApplyNeumannChannel();
      QQyx2.ApplyNeumannChannel();
      break;
    // box
    case 3:
    case 4:
      //QQxx1.ApplyNeumann();
      //QQyx1.ApplyNeumann();
      //QQxx2.ApplyNeumann();
      //QQyx2.ApplyNeumann();
      QQxx1.AnchorAtBoxWall(anchorSign1*wallAnchor);
      QQyx1.SetBoundaryValue(0);
      QQxx2.AnchorAtBoxWall(anchorSign2*wallAnchor);
      QQyx2.SetBoundaryValue(0);
      break;
    // pbc with bdry layer
    default:
      QQxx1.ApplyPBC();
      QQyx1.ApplyPBC();
      QQxx2.ApplyPBC();
      QQyx2.ApplyPBC();
  }
}
void TwoFluidWetting::BoundaryConditionsFields2()
{
  switch(BC)
  {
    // pbc without bdry layer (nothing to do)
    case 0:
      break;
    // channel
    case 1:
    case 2:
      uy1     .ApplyNeumannChannel();
      ux1     .ApplyNeumannChannel();
      phi_old .ApplyNeumannValueChannel(h_surf/(2*KK));
      phi     .ApplyNeumannValueChannel(h_surf/(2*KK));
      uy2     .ApplyNeumannChannel();
      ux2     .ApplyNeumannChannel();
      n1      .ApplyNeumannChannel();
      n2      .ApplyNeumannChannel();
      nC      .ApplyNeumannChannel();
      uxC     .ApplyNeumannChannel();
      uyC     .ApplyNeumannChannel();

      break;
    // box
    case 3:
    case 4:
      nC      .ApplyNeumann();
      n1      .ApplyNeumann();
      n2      .ApplyNeumann();
      uxC     .ApplyNeumann();
      uyC     .ApplyNeumann();
      uy1     .ApplyNeumann();
      ux1     .ApplyNeumann();  
      phi_old .ApplyNeumannValue(h_surf/(2*KK));
      phi     .ApplyNeumannValue(h_surf/(2*KK));
      uy2     .ApplyNeumann();
      ux2     .ApplyNeumann();   
      
      /*
      sigmaXX1.ApplyNeumann();
      sigmaYY1.ApplyNeumann();
      sigmaYX1.ApplyNeumann();
      sigmaXY1.ApplyNeumann();    

      sigmaXX2.ApplyNeumann();
      sigmaYY2.ApplyNeumann();
      sigmaYX2.ApplyNeumann();
      sigmaXY2.ApplyNeumann();
      */

      break;
    // pbc with bdry layer
    default:
      nC      .ApplyPBC();
      n1      .ApplyPBC();
      n2      .ApplyPBC();
      ux1     .ApplyPBC();
      uy1     .ApplyPBC();
      ux2     .ApplyPBC();
      uy2     .ApplyPBC();
      uxC     .ApplyPBC();
      uyC     .ApplyPBC();
      phi     .ApplyPBC();
      phi_old .ApplyPBC();   
      
      /*
      sigmaXX2.ApplyPBC();
      sigmaYY2.ApplyPBC();
      sigmaYX2.ApplyPBC();
      sigmaXY2.ApplyPBC();

      sigmaXX1.ApplyPBC();
      sigmaYY1.ApplyPBC();
      sigmaYX1.ApplyPBC();
      sigmaXY1.ApplyPBC();
      */
  }
}
void TwoFluidWetting::BoundaryConditionsPhi(){  
  for (unsigned nbcphi=1; nbcphi<=2; nbcphi++){
    #pragma omp parallel for num_threads(nthreads) if(nthreads)
    for(unsigned k=0; k<DomainSize; k++) {    
      BoundaryConditionsPhiAtNode(k);
    }
  }
}

void TwoFluidWetting::BoundaryConditionsPhiAtNode(unsigned k){
  
  return;
  const double bigP = bigPhi[k];
  if (bigP > 0.9 || bigP < 0.1) return; //this function is only for the interface
  /*
  const auto& d = get_neighbours(k);
  const double dxBigP = derivX(bigPhi, d, sB);
  const double dyBigP = derivY(bigPhi, d, sB);
  const double modBigP = sqrt(dxBigP*dxBigP + dyBigP*dyBigP);
  */
  const unsigned x = GetXPosition(k); const unsigned y = GetYPosition(k);
  double x_prev = x - cos(angDeg[k]); double y_prev = y - sin(angDeg[k]);
  
  const double phi_interp = phi.Interpolate(x_prev, y_prev);
  phi[k] = phi_interp;

  //
  const double QQxx_interp = QQxx1.Interpolate(x_prev, y_prev);
  QQxx1[k] = QQxx_interp;

  const double QQyx_interp = QQyx1.Interpolate(x_prev, y_prev);
  QQyx1[k] = QQyx_interp;
  //

  /*
  const double nnC = nC[k]; const double nn1 = n1[k];
  double vxC = (ff[k][1] - ff[k][2] + ff[k][5] - ff[k][6] - ff[k][7] + ff[k][8])/nnC;
  double vyC = (ff[k][3] - ff[k][4] + ff[k][5] - ff[k][6] + ff[k][7] - ff[k][8])/nnC;
  double vx1 = (gg2[k][1] - gg2[k][2] + gg2[k][5] - gg2[k][6] - gg2[k][7] + gg2[k][8])/nnC;
  double vy1 = (gg2[k][3] - gg2[k][4] + gg2[k][5] - gg2[k][6] + gg2[k][7] - gg2[k][8])/nn1;
  */
  
  /*
  const auto fe  = CorrectLBAtNode(GetEquilibriumDistribution(vxC, vyC, nnC), k);
  const auto ge1 = CorrectLBAtNode(GetEquilibriumDistribution(vx1, vy1, nn1), k); 
  for(unsigned v=0; v<lbq; ++v)
  {
    ff[k][v] = fe[v];
    gg1[k][v] = ge1[v];
    gg2[k][v] = ge1[v];
  }
  */
  /*
  ff.ApplyFreeSlipAtNode(static_cast<unsigned>(angDeg[k]), k);
  gg1.ApplyFreeSlipAtNode(static_cast<unsigned>(angDeg[k]), k);
  gg2.ApplyFreeSlipAtNode(static_cast<unsigned>(angDeg[k]), k);
  */
  
}
void TwoFluidWetting::setAngDeg(){
  const double xCentre = (LX-1.0)/2.0, yCentre= (LY-1.0)/2.0;  

  #pragma omp parallel for num_threads(nthreads) if(nthreads)
  for(unsigned k=0; k<DomainSize; k++) {    
    unsigned x = GetXPosition(k); unsigned y = GetYPosition(k);
    double angDeg_temp = atan2(y-yCentre, x-xCentre);
    // if (angDeg_temp > - M_PI/8 && angDeg_temp <= M_PI/8) angDeg[k] = 1;
    // else if (angDeg_temp >  M_PI/8 && angDeg_temp <= 3*M_PI/8) angDeg[k] = 5;
    // else if (angDeg_temp >  3*M_PI/8 && angDeg_temp <= 5*M_PI/8) angDeg[k] = 3;
    // else if (angDeg_temp >  5*M_PI/8 && angDeg_temp <= 7*M_PI/8) angDeg[k] = 7;
    // else if (angDeg_temp >  -3*M_PI/8 && angDeg_temp <= -M_PI/8) angDeg[k] = 8;
    // else if (angDeg_temp >  -5*M_PI/8 && angDeg_temp <= -3*M_PI/8) angDeg[k] = 4;
    // else if (angDeg_temp >  -7*M_PI/8 && angDeg_temp <= -5*M_PI/8) angDeg[k] = 6;
    // else angDeg[k] = 2;
    angDeg[k] = angDeg_temp;
  }
}

void TwoFluidWetting::Step()
{
  if (my_timer < act_start_timer){
    preinit_flag = true; //turns off activity
  }
  else if (act_stop_timer > 0 && my_timer > act_stop_timer){
    preinit_flag = true;
  }
  else{
    preinit_flag = false;
  }

  BoundaryConditionsFields();
  BoundaryConditionsPhi();  
  UpdateNemQuantities(true);
  UpdateFluidQuantities();
  
  BoundaryConditionsFields2();
  BoundaryConditionsPhi();
 
  for (unsigned lcv = 1; lcv<=nsteps_FD; lcv++){
    this->UpdateNemFields(true);
    UpdateNemQuantities(true);
  }
  this->UpdateFluidFields(true);

  BoundaryConditionsLB();  
  Move();
  BoundaryConditionsPhi();
  
  swap(ux_old1.get_data(), ux1.get_data());
  swap(uy_old1.get_data(), uy1.get_data());
  swap(ux_old2.get_data(), ux2.get_data());
  swap(uy_old2.get_data(), uy2.get_data());
  swap(phi_old.get_data(), phi.get_data());


  // corrector steps
  for(unsigned n=1; n<=npc; ++n)
  { //not implemented yet - do not use
    cout << "Corrector step not implemented yet !!!! \n";
    BoundaryConditionsFields();
    UpdateNemQuantities(true);
    UpdateFluidQuantities();
    BoundaryConditionsFields2();
    for (unsigned lcv = 1; lcv<=nsteps_FD; lcv++){
      this->UpdateNemFields(false);
      UpdateNemQuantities(false);
    }    
    this->UpdateFluidFields(false);

    swap(ux_old1.get_data(), ux1.get_data());
    swap(uy_old1.get_data(), uy1.get_data());
    swap(ux_old2.get_data(), ux2.get_data());
    swap(uy_old2.get_data(), uy2.get_data());
    swap(phi_old.get_data(), phi.get_data());
  }

  my_timer +=1; //just for the Malevanets original simulation
  if (massRealloc && my_timer%nRealloc == (nRealloc-1)){ //just before printing
    //cout << "Reallocating masses ... \n";

    //calculate material to reallocate
    nInnerNodes = 0; extraMaterialIn = 0;
    #pragma omp parallel for num_threads(nthreads) if(nthreads)
    for(unsigned k=0; k<DomainSize; ++k){
      if (bigPhi[k] >= BPcrit){
        nInnerNodes+=1; //inner node
        extraMaterialIn+=(phi[k]-0.50); //outer node
      } 
    }

    //
    double phi_tilde = -extraMaterialIn/nInnerNodes;
    //rearrange uniformly
    #pragma omp parallel for num_threads(nthreads) if(nthreads)
    for(unsigned k=0; k<DomainSize; ++k){ 
      if (bigPhi[k] >= BPcrit){ //inner node
        phi[k] = (phi[k] + phi_tilde)/(1+phi_tilde); //add the extra material in phi
        gg2[k][0] += rho1*phi_tilde; //and in the n1 tracking
        gg1[k][0] += rho1*phi_tilde; //and in the n1 tracking
      } 
      else{ //outer node
        // phi[k] = 0.50;
        // auto& g2 = gg2[k];
        // g2[0] = rho1 - (g2[1] + g2[2] + g2[3] + g2[4] + g2[5] + g2[6] + g2[7] + g2[8]);
      }      
    }
    //

    //rearrange to edge nodes only
    
    /*double N_edge = (radius-1)*2*M_PI;
    double phi_edge = extraMaterialIn/N_edge;
    #pragma omp parallel for num_threads(nthreads) if(nthreads)
    for(unsigned k=0; k<DomainSize; ++k){
      if (bigPhi[k] > 0.90 && bigPhi[k]<0.95) { //edge node
        phi[k] = (phi[k] + phi_edge)/(1+phi_edge); //add the extra material in phi
        gg2[k][0] += rho1*phi_edge; //and in the n1 tracking
        gg1[k][0] += rho1*phi_edge; //and in the n1 tracking        
      }
    }
      */
  }
}

void TwoFluidWetting::RuntimeChecks()
{
  // check that the sum of f is constant
  {
    double fcheck1 = 0;
    double fcheck2 = 0;
    double fcheck3 = 0;
    for(unsigned k=0; k<DomainSize; ++k){
        fcheck1 = accumulate(begin(ff[k]), end(ff[k]), fcheck1);
        fcheck2 = accumulate(begin(gg2[k]), end(gg2[k]), fcheck2);
        fcheck3 = accumulate(begin(gg1[k]), end(gg1[k]), fcheck3);
    }
    cout << "ff::  fcheck1: " << fcheck1 << "/" << ftot1 << '\n';
    cout << "gg2:: fcheck2: " << fcheck2 << "/" << ftot2 << '\n';
    cout << "gg1:: fcheck3: " << fcheck3 << "/" << ftot2 << '\n';
    if(abs(ftot1-fcheck1)>0.01*ftot1){
      cout << "f is not conserved !!! ";
      for(unsigned k=0; k<DomainSize; ++k){
        ff[k][0] = ff[k][0] + (ftot1-fcheck1)/DomainSize;
      }
    } 
    if(abs(ftot2-fcheck2)>0.01*ftot2){
      cout << "g2 is not conserved !!! ";
      for(unsigned k=0; k<DomainSize; ++k){
        gg2[k][0] = gg2[k][0] + (ftot2-fcheck2)/DomainSize;
      }
    }
    if(abs(ftot2-fcheck3)>0.01*ftot2){
      cout << "g1 is not conserved !!! ";
      for(unsigned k=0; k<DomainSize; ++k){
        gg1[k][0] = gg1[k][0] + (ftot2-fcheck3)/DomainSize;
      }
    }
      //throw error_msg("f is not conserved !!! ");
      
  }

  // check that n is positive
  {    
    for(unsigned k=0; k<DomainSize; ++k){
        if (n1[k]<0 || n2[k]<0)
          throw error_msg("n is negative somewhere");
    } 
  } 

  /* // check that n is conserved  - turned OFF
  {
    double n1_check = 0;
    double n2_check = 0;
    for(unsigned k=0; k<DomainSize; ++k){
        n1_check += n1[k];
        n2_check += n2[k];
    }
    cout << "n1 check: " << n1_check << "/" << ntot1 << '\n';
    cout << "n2 check: " << n2_check << "/" << ntot2 << '\n';
    if(abs((ntot1-n1_check)/ntot1)>0.01 || abs((ntot2-n2_check)/ntot2)>0.01) //1% particle loss tolerance
      throw error_msg("n is not conserved ");
  } */
  /* //check that phi is conserved - turned OFF
  {
    double pcheck = 0;
    for(unsigned k=0; k<DomainSize; ++k)
        pcheck += phi[k];
    cout << "pcheck: " << pcheck << "/" << ptot << '\n';
    if(abs(ptot-pcheck)>1)
      throw error_msg("phi is not conserved (", ptot, "/", pcheck, ")");
  }
  */
}

option_list TwoFluidWetting::GetOptions()
{
  // model specific options
  opt::options_description model_options("Model options");
  model_options.add_options()
    ("GammaQ1", opt::value<double>(&GammaQ1),
     "Q-tensor mobility")
    ("xi1", opt::value<double>(&xi1),
     "tumbling/aligning parameter")
    ("GammaQ2", opt::value<double>(&GammaQ2),
     "Q-tensor mobility")
    ("xi2", opt::value<double>(&xi2),
     "tumbling/aligning parameter")
    ("tau1", opt::value<double>(&tau1),
     "viscosity 1")
    ("tau2", opt::value<double>(&tau2),
     "viscosity 2")
    ("rho1", opt::value<double>(&rho1),
     "fluid density 1")
    ("rho2", opt::value<double>(&rho2),
     "fluid density 2")
    ("rel_fric", opt::value<double>(&rel_fric),
     "relative friction between two velocities")
    ("friction1", opt::value<double>(&friction1),
     "friction from confinement")
    ("friction2", opt::value<double>(&friction2),
     "friction from confinement")
    ("CC1", opt::value<double>(&CC1),
     "coupling constant")
    ("LL1", opt::value<double>(&LL1),
     "elastic constant")    
    ("CC2", opt::value<double>(&CC2),
     "coupling constant")
    ("LL2", opt::value<double>(&LL2),
     "elastic constant")
    ("LL12", opt::value<double>(&LL12),
     "elastic coupling constant")
    ("KK", opt::value<double>(&KK),
     "binary gradient constant")
    ("APP", opt::value<double>(&APP),
     "A value for Big Phi")
    ("KPP", opt::value<double>(&KPP),
     "K value for Big Phi")
    ("gammaPP", opt::value<double>(&gammaPP),
     "gamma value for Big Phi")
    ("nPP_init", opt::value<int>(&nPP_init),
     "number of timesteps for initializing bigPhi")
    ("zeta1", opt::value<double>(&zeta1),
     "activity parameter, fluid 1")   
    ("zeta2", opt::value<double>(&zeta2),
     "activity parameter, fluid 2")
    ("zetaIso1", opt::value<double>(&zeta_ISO1),
     "isotropic activity parameter for fluid 1")   
    ("act_start_timer", opt::value<int>(&act_start_timer),
     "Activity start time")
    ("act_stop_timer", opt::value<int>(&act_stop_timer),
     "Activity stop time")
    ("landA", opt::value<double>(&landA),
     "Landau free energy A")
    ("landB", opt::value<double>(&landB),
     "Landau free energy B")
    ("nsteps_FD", opt::value<unsigned>(&nsteps_FD),
     "number of steps for the finite-difference method")
    ("n_preinit", opt::value<int>(&n_preinit),
     "number of preinitialization steps") 
    ("W_anchor", opt::value<double>(&W_anchor),
     "anchoring strength")
    ("S_anchor", opt::value<double>(&S_anchor),
     "anchhoring order parameter")
    ("theta_anchor", opt::value<double>(&theta_anchor_deg),
     "anchoring angle parameter")
    ("K12_NR", opt::value<double>(&K12_NR),
     "non-reciprocal elastic coupling constant")
    ("K21_NR", opt::value<double>(&K21_NR),
     "non-reciprocal elastic coupling constant")  
    ("isGuo", opt::value<bool>(&isGuo),
     "Is the forcing scheme Guo? If not, use Shan-Chen")
    ("backflow_on", opt::value<bool>(&backflow_on),
     "Do we have backflow?")
    ("phi_FreeE", opt::value<bool>(&phi_FreeE),
     "Does free energy depend on phi") 
    ("fix_Q1", opt::value<bool>(&fix_Q1),
     "Fix Q1, for debugging")    
    ("L_thermo", opt::value<double>(&L_thermo),
     "Strength of thermodynamic anchoring")
    ("Q_kBT", opt::value<double>(&Q_kBT),
     "magnitude of Q fluctuations (if present)")
    ("Snem1", opt::value<double>(&Snem1),
     "magnitude of nematic order in Fluid 1")
    ("Snem2", opt::value<double>(&Snem2),
     "magnitude of nematic order in Fluid 2")
    ("growthRate", opt::value<double>(&growthRate),
     "growthRate")
    ("h_surf", opt::value<double>(&h_surf),
     "Strength of wetting term")
    ("extraVisc", opt::value<double>(&extraVisc),
     "Extra viscosity of Phi=0")
    ("modFE", opt::value<double>(&modFE),
     "Stronger free energy when Phi=0")
    ("modFric", opt::value<double>(&modFric),
     "Stronger free energy when Phi=0")
    ("modRF", opt::value<double>(&modRF),
     "Stronger free energy when Phi=0")
    ("massRealloc", opt::value<bool>(&massRealloc),
     "Do we reallocate masses?")
    ("nRealloc", opt::value<int>(&nRealloc),
     "how often do we reallocate?")     
     ;

  // init config options
  opt::options_description config_options("Initial configuration options");
  config_options.add_options()
    ("config", opt::value<string>(&init_config),
     "initial configuration")
    ("level", opt::value<double>(&level),
     "starting thickness of the nematic region")
    ("init_seed", opt::value<int>(&init_seed),
     "starting seed")
    ("conc", opt::value<double>(&conc),
     "starting phi concentration of nematic region")
    ("radius", opt::value<double>(&radius),
     "radius of the initial circle")
    ("angle", opt::value<double>(&angle_deg),
     "initial angle to x direction (in degrees)")
    ("angle2", opt::value<double>(&angle_deg2),
     "initial angle to x direction (in degrees)")
    ("noise", opt::value<double>(&noise),
     "size of initial variations")
    ("GD_config", opt::value<string>(&GD_config),
     "growth death configuration") 
    ;

  return { model_options, config_options };
}

