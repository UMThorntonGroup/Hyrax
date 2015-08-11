/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  20 March 2015
*
*************************************************************************/

#ifndef AUXAMRNUCLEATIONPROBABILITY_H
#define AUXAMRNUCLEATIONPROBABILITY_H

#include "AuxKernel.h"

class AuxAMRNucleationProbability;

template<>
InputParameters validParams<AuxAMRNucleationProbability>();

/**
 *  AuxAMRNucleationProbability handles the nucleation probability (P_nm) calculation in the concurrent
 *  nucleation and growth algorithm first proposed by J.P. Simmons (2000).
 *  Returns the nucleation probability for an adaptive mesh in 2 or 3 dimensions.
 */

class AuxAMRNucleationProbability : public AuxKernel
{
public:
  AuxAMRNucleationProbability(const InputParameters & params);

protected:
  /**
   * computeValue()
   * @return returns the nucleation probability for each element, p_nm = 1 - exp(-1*j_star*dt*vol)
   */

  virtual Real computeValue();

private:
  VariableValue & _coupled_nuc_rate;
  VariableValue & _coupled_OP;
  Real _OP_threshold;
  Real _2D_mesh_height;
  const unsigned int _mesh_dimension;

};

#endif //AUXAMRNUCLEATIONPROBABILITY_H
