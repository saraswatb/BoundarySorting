// ======================================================================
// Model declaration (gets compiled in models.cpp)
// model headers and declare_model must be consistent!


// model headers
#include "models/TwoFluidFull.hpp"
#include "models/TwoFluidWetting.hpp"


void DeclareModels()
{  

  declare_model<TwoFluidFull>(
     "TwoFluidFull",
     "A two fluid model with different velocity fields as described in Malevanets et al (1999, 2000)."
     "Modifications: Added active forces on both components separately."     
     ); 

  declare_model<TwoFluidWetting>(
     "TwoFluidWetting",
     "A two fluid model specifically for the confined case. Uses a master phase field to model confinement region."     
     ); 
         
}
