/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  10 April 2014
*
*************************************************************************/

#ifndef SORET_H
#define SORET_H

#include "Diffusion.h"

//forward declaration
class Soret;

template<>
InputParameters validParams<Soret>();

class Soret : public Diffusion
{
public:

  Soret(const InputParameters & parameters);

protected:

  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  // virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  const MaterialProperty<Real> & _L1Q;

private:
  VariableValue & _T;
  VariableGradient & _grad_T;
};

#endif //SORET_H
