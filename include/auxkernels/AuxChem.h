/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  4 September 2013
*
*************************************************************************/

#ifndef AUXCHEM_H
#define AUXCHEM_H

#include "AuxKernel.h"

//forward declaration
class AuxChem;

template<>
InputParameters validParams<AuxChem>();

class AuxChem : public AuxKernel
{
public:
  AuxChem(const InputParameters & parameters);

protected:
  virtual Real computeValue();

  virtual Real computeEnergy(Real & conserved, Real & nonconserved, bool matrix);
  virtual Real computeDifferential(Real & coupled_conserved,Real & nonconserved);

  virtual Real computeFchem(Real & conserved, Real & nonconserved);

  virtual Real computeDfchemDcons(Real & coupled_conserved, Real & coupled_nonconserved);

  virtual Real computeDfchemDnoncons(Real & coupled_conserved, Real & coupled_nonconserved);

  VariableValue & _coupled_cons;
  VariableValue & _coupled_noncons;

  Real _precip_conserved;
  Real _precip_nonconserved;

  const MaterialProperty<Real> & _a1;
  const MaterialProperty<Real> & _a2;
  const MaterialProperty<Real> & _a3;
  const MaterialProperty<Real> & _a4;
  const MaterialProperty<Real> & _c1;
  const MaterialProperty<Real> & _c2;

  // unsigned int _n_variants;
  unsigned int _noncons_var_num;

private:
};

#endif //AUXCHEM_H
