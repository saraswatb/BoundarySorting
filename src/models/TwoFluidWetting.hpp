#ifndef MODELS_TwoFluidWetting_HPP_
#define MODELS_TwoFluidWetting_HPP_

#include "models.hpp"

class TwoFluidWetting : public Model
{
protected:
  int my_timer = 0; 
  int act_start_timer = 0;
  int act_stop_timer = -1;

  LBField ff, ff_tmp;
  LBField gg2, gg2_tmp; 
  LBField gg1, gg1_tmp;
  
  ScalarField QQxx1, QQyx1;
  ScalarField QQxx2, QQyx2;  
  ScalarField ux1, uy1;
  ScalarField ux2, uy2;  
  ScalarField n1, n2, phi, phi_old;
  ScalarField dxPPhi, dyPPhi;
  ScalarField nC, uxC, uyC, nC_tmp;  

  //master phase field, which identifies IN and OUT regions
  const bool dropletMode = true; 
  ScalarField bigPhi, bigPhi_tmp, bigMu;
  double APP=0.4, KPP=0.1, gammaPP=0.005; int nPP_init = 10000;
  bool fixedBigPhi = true; 
  //Phi < BPcrit is outside. Modifiers and parameters for both regions
  double BPcrit = 0.9, extraVisc = 10, modNFE = 1, modFE = 100; 
  double modFric = 0.1, modRF = 20;
    
  //derivatives and derived fields
  ScalarField HHxx1, HHyx1; 
  ScalarField HHxx2, HHyx2; 
  ScalarField dxQQxx1, dyQQxx1, dxQQyx1, dyQQyx1;
  ScalarField dxQQxx2, dyQQxx2, dxQQyx2, dyQQyx2;
  ScalarField sigmaXX1, sigmaYY1, sigmaYX1, sigmaXY1;
  ScalarField sigmaXX2, sigmaYY2, sigmaYX2, sigmaXY2;
  ScalarField angDeg;
  
  ScalarField subForceX, subForceY;
  ScalarField FFx1, FFy1, FFx2, FFy2;  

  //Individual fluid and nematic parameters
  double rho1 = 20., rho2 = 20., eta;
  double GammaQ1, xi1, tau1 = 1.0, friction1, LL1, CC1, zeta1;
  double GammaQ2, xi2, tau2 = 1.0, friction2, LL2, CC2, zeta2; 
  double Snem1 = 1, Snem2 = 1;
  double zeta_ISO1 = 0;

  //Fluid coupling parameters
  double landA = 0.004, landB = 0.0;
  double rel_fric = 4, KK =5;
  double LL12 = 0.0, K12_NR = 0, K21_NR = 0; 
  ScalarField BulkFE_LDG_12; 
    
  //configuration options
  double level, conc, angle_deg, angle, noise, radius, angle_deg2, angle2;
  std::string init_config;
  unsigned npc = 0, nsteps_FD = 1; double dt_FD = 1;
  int n_preinit = 0; bool preinit_flag = false; int init_seed = 2441139;
  
  //LB checks, balances, and fudges
  double ftot1 = 0, ftot2 = 0, ntot1 = 0, ntot2 = 0;
  double n_TOL = 1e-8; //phi -> 0, 1 fix
  bool isGuo = false; //if false, do by Shan-Chen instead
  //for mass reallocation code
  bool massRealloc = true;
  double extraMaterialIn; int nInnerNodes, nRealloc=100;      
  
  //effect of backflow
  bool backflow_on = false; bool phi_FreeE = false;

  //thermodynamic anchoring, wall anchoring and wetting
  bool wallAnchor = false; double W_anchor, S_anchor = 1, theta_anchor, theta_anchor_deg; 
  bool thermoAnchor = false; double L_thermo;
  ScalarField thermoAnchor_delMu;
  double h_surf = 0.0;

  // source terms, fluctuations, and special cases.
  double growthRate; std::string GD_config;
  double Q_kBT = 0; bool Q_fluct = false; 
  const bool viscChange = false, malvnt_test = false, fric_track_test = false; 
  // double u_kBT = 0; bool u_fluct = false; //velocity fluctuations not implemented yet
  
  //Fields for debugging/analysis
  ScalarField ux_old1, uy_old1, ux_old2, uy_old2; 
  ScalarField advX1, advY1, advX2, advY2;
  bool debugMode = false;
  bool fix_Q1 = false; const bool fix_vc = false, fix_v1 = false;
  
  /** Update fields using predictor-corrector method */
  virtual void UpdateNemFields(bool);
  virtual void UpdateNemQuantities(bool);
  virtual void UpdateFluidFields(bool);
  virtual void UpdateFluidQuantities();

  void UpdateNemFieldsAtNode(unsigned, bool);
  void UpdateNemQuantitiesAtNode(unsigned, bool);
  void UpdateFluidFieldsAtNode(unsigned, bool);
  void UpdateFluidQuantitiesAtNode_p1(unsigned);
  void UpdateFluidQuantitiesAtNode_p2(unsigned);
  
  void EvaluatePhiQuantities(unsigned);
  void UpdateBigPhi();
  void UpdateBigPhiFieldAtNode(unsigned);
  void UpdateBigPhiQntyAtNode(unsigned);
  LBNode CorrectLBAtNode(LBNode, unsigned);
 
  virtual void BoundaryConditionsLB();
  virtual void BoundaryConditionsFields();
  virtual void BoundaryConditionsFields2();
  
  //Emergency patch: Apply correct boundary conditions for phi and S
  virtual void BoundaryConditionsPhi();
  void BoundaryConditionsPhiAtNode(unsigned);
  void setAngDeg();
  
  //virtual void ReallocateMaterial(unsigned);

  /** Move the LB particles */
  void Move();

public:
  TwoFluidWetting() = default;
  TwoFluidWetting(unsigned, unsigned, unsigned);
  TwoFluidWetting(unsigned, unsigned, unsigned, GridType);

  /** Configure a single node
   *
   * This allows to change the way the arrays are configured in derived
   * classes, see for example TwoFluidWettingFreeBoundary.
   * */
  virtual void ConfigureAtNode(unsigned);

  // functions from base class Model
  virtual void Initialize();
  virtual void Step();
  virtual void Configure();
  virtual void RuntimeChecks();
  virtual option_list GetOptions();

  /** Serialization of parameters (do not change) */
  template<class Archive>
  void serialize_params(Archive& ar)
  {
    ar 
       & auto_name(rho1)
       & auto_name(rho2)
       & auto_name(rel_fric) 
       & auto_name(friction1)  
       & auto_name(friction2)   
      
       & auto_name(landA)
       & auto_name(landB)
       & auto_name(KK) 
       
       & auto_name(zeta1)
       & auto_name(zeta_ISO1)       
       & auto_name(zeta2)
       & auto_name(LL1)
       & auto_name(LL2)  
       & auto_name(LL12)
              
       & auto_name(GammaQ1)       
       & auto_name(xi1)             
       & auto_name(CC1)
       & auto_name(GammaQ2)
       & auto_name(xi2)
       & auto_name(CC2)

       & auto_name(nsteps_FD)
       & auto_name(init_config)
       & auto_name(W_anchor)
       & auto_name(theta_anchor)
       & auto_name(L_thermo)

       & auto_name(Snem1)
       & auto_name(Snem2)
       & auto_name(Q_kBT)
       & auto_name(extraVisc)
       & auto_name(modNFE)
       & auto_name(modFE)
       & auto_name(modFric)
       & auto_name(modRF)
       //& auto_name(growthRate)
       //& auto_name(GD_config)
       //& auto_name(h_surf)

       & auto_name(level)
       & auto_name(conc)
       & auto_name(angle)
       & auto_name(radius)
       & auto_name(noise)
       & auto_name(act_start_timer)
       & auto_name(act_stop_timer)
       & auto_name(init_seed)
       ;
  }

  /** Serialization of the current frame (time snapshot) */
  template<class Archive>
  void serialize_frame(Archive& ar)
  {
    ar & auto_name(ff)
       & auto_name(gg1)
       & auto_name(gg2)
       & auto_name(QQxx1)
       & auto_name(QQyx1)
       & auto_name(QQxx2)
       & auto_name(QQyx2)       
       & auto_name(ux1)
       & auto_name(uy1)
       & auto_name(bigPhi)
       ;
  }
};

/*  //removed parameters from "serialize_frame"
      & auto_name(FFx1)
      & auto_name(FFy1)
      & auto_name(FFx2)
      & auto_name(FFy2)
      & auto_name(subForceX)
      & auto_name(subForceY)
      
      & auto_name(n1)
      & auto_name(ux1)
      & auto_name(uy1)
*/

#endif//MODELS_TwoFluidWetting_HPP_
