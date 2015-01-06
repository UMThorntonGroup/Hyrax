/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  15 December 2011
*
*************************************************************************/

#ifndef AUXSUPERSATURATION_H
#define AUXSUPERSATURATION_H

#include "AuxKernel.h"

//forward declaration
class AuxSupersaturation;

/**
 * AuxSupersaturation handles the supersaturation calculation (tyically C - C1) in the concurrent
 * nucleation and growth algorithm first proposed by J.P. Simmons (2000).
 */

template<>
InputParameters validParams<AuxSupersaturation>();

class AuxSupersaturation : public AuxKernel
{
public:
  AuxSupersaturation(const std::string & name, InputParameters parameters);

protected:

  /**
   * computeValue() 
   * @return returns the supersaturation (C-C1), element average value
   */

  virtual Real computeValue();
  
private:

  VariableValue & _coupled_conc;  ///< PDE variable (concentration)
  MaterialProperty<Real> & _c1;   ///< parameter from landau polynomial (first energy well position in c-space)
  Real _supersaturation;	  ///< C-C1
};

#endif //AUXSUPERSATURATION_H
