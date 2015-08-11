/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  1 November 2013
*
*************************************************************************/

#ifndef HEAT_H
#define HEAT_H

#include "Diffusion.h"

//forward declaration
class Heat;

template<>
InputParameters validParams<Heat>();

class Heat : public Diffusion
{
public:

  Heat(const InputParameters & parameters);

protected:

  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  const MaterialProperty<Real> & _diffusivity;
//  const MaterialProperty<Real> & _dDiffusivity_dT;

private:
  unsigned int _w_var;
  unsigned int _c_var;
  //unsigned int _n_var;

  VariableValue & _w;
  VariableValue & _c;
  // VariableValue & _n;

  unsigned int _n_OP_vars;

  std::vector<unsigned int> _n_var;
  std::vector<VariableValue *> _coupled_OP_vars;
};

#endif //HEAT_H
