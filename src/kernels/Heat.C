/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  1 November 2013
*
*************************************************************************/

#include "Heat.h"

template<>
InputParameters validParams<Heat>()
{
  InputParameters params = validParams<Diffusion>();

  params.addRequiredCoupledVar("w", "chemical potential");
  params.addRequiredCoupledVar("c", "concentration");
  //params.addRequiredCoupledVar("n", "structural order parameter");

  params.addRequiredParam<int>("n_OP_vars", "# of coupled OP variables");
  params.addRequiredCoupledVar("OP_var_names", "Array of coupled OP variable names");

  return params;
}


Heat::Heat(const InputParameters & parameters) :
    Diffusion(parameters),
    _diffusivity(getMaterialProperty<Real>("thermal_diffusivity")),
    _w_var(coupled("w")),
    _c_var(coupled("c")),
    //_n_var(coupled("n")),
    _w(coupledValue("w")),
    _c(coupledValue("c")),
    // _n(coupledValue("n"))
    _n_OP_vars(getParam<int>("n_OP_vars"))
    //_OP_number(getParam<int>("OP_number"))
    // _dDiffusivity_dT(getMaterialProperty<Real>("dThermal_diffusivity_dT"))
{
    // Create a vector of the coupled OP variables and set = 0 the one that the kernel
  // is operating on
  if(_n_OP_vars != coupledComponents("OP_var_names"))
    mooseError("Please match the number of orientation variants to coupled OPs (ACCoupledCalphad).");

  _n_var.resize(_n_OP_vars);
  _coupled_OP_vars.resize(_n_OP_vars);

  for (unsigned int i=0; i<_n_OP_vars; i++)
  {
    _n_var[i] = coupled("OP_var_names", i);
    _coupled_OP_vars[i] = &coupledValue("OP_var_names", i);
  }
}


Real
Heat::computeQpResidual()
{
  return _grad_test[_i][_qp]*_diffusivity[_qp]*_grad_u[_qp];
}

Real
Heat::computeQpJacobian()
{
 return _grad_phi[_j][_qp]*_diffusivity[_qp]*_grad_test[_i][_qp];
}

Real
Heat::computeQpOffDiagJacobian(unsigned int jvar)
{
  for (unsigned int i=0; i< _n_OP_vars; i++)
  {
    if (jvar == _n_var[i])
      return 0;
  }

   if (jvar == _w_var)
    return 0;

  else if (jvar == _c_var)
    return 0;

   //else if (jvar == _n_var)
   //  return 0;

  else
    mooseError("Screwed up Heat::computeQpOffDiagonalJacobian.");
}
