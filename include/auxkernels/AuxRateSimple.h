/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  18 April 2013
*
*************************************************************************/

#ifndef AUXRATESIMPLE_H
#define AUXRATESIMPLE_H

#include "AuxKernel.h"

class AuxRateSimple;

template<>
InputParameters validParams<AuxRateSimple>();

/**
 *  AuxRateSimple handles the nucleation rate (j_star) calculation in the concurrent
 *  nucleation and growth algorithm first proposed by J.P. Simmons (2000).
 */

class AuxRateSimple : public AuxKernel
{
public:
  AuxRateSimple(const InputParameters & parameters);

protected:
  /**
   * computeValue()
   * @return returns the nucleation rate (element average value), j_star.
   * Handles adaptivity and different mesh dimensions.
   * j_star = Kn1 * exp(-1*Kn2 / supersaturation)
   */

  virtual Real computeValue();

private:
  const unsigned int _mesh_dimension;

  VariableValue & _coupled_energy; ///< area/volume free energy change of transformation
  Real _Kn1;                       ///< First nucleation rate value
  Real _Kn2;                       ///< Second nucleation rate value

  unsigned int _n_OP_vars;
  std::vector<VariableValue *> _coupled_OP_vars;
};

#endif //AUXRATESIMPLE_H
