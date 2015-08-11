/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  30 August 2013
*
*************************************************************************/

#include "ACFerroelectric.h"

template<>
InputParameters validParams<ACFerroelectric>()
{
  InputParameters params = validParams<ACBulk>();

  params.addRequiredCoupledVar("OP_var_names", "Array of coupled OP variable names");
  params.addRequiredParam<int>("n_OP_vars", "number of coupled OP variables");
  params.addRequiredParam<int>("OP_number", "# of OP for this kernel, starting from 1");

  return params;
}

ACFerroelectric::ACFerroelectric(const InputParameters & parameters) :
    ACBulk(parameters),
    _a1(getMaterialProperty<Real>("a1")),
    _a11(getMaterialProperty<Real>("a11")),
    _a12(getMaterialProperty<Real>("a12")),
    _n_OP_vars(getParam<int>("n_OP_vars")),
    _OP_number(getParam<int>("OP_number"))
{
  //make sure you didn't cock up the input file
  if(_n_OP_vars != coupledComponents("OP_var_names"))
    mooseError("Please match the # of polarizations to coupled OPs (ACFerroelectric)");


  // get the coupled polarizatons and this kernel's polarization variable set up correctly
  _coupled_OP_vars.resize(_n_OP_vars);

  for(unsigned int i=0; i< _n_OP_vars; i++)
  {
    if(i == _OP_number-1)
    {
      _coupled_OP_vars[i] = NULL;
    }
    else
    {
      _coupled_OP_vars[i] = &coupledValue("OP_var_names", i);
    }
  }
}

Real
ACFerroelectric::computeDFDOP(PFFunctionType type)
{
  Real square_sum(0.0);

  for(unsigned int i=0; i< _n_OP_vars; i++)
  {
    if(i != _OP_number-1)
      square_sum += ((*_coupled_OP_vars[i])[_qp])*((*_coupled_OP_vars[i])[_qp]);
  }

  switch (type)
  {
  case Residual:
    return _a1[_qp]*2.0*_u[_qp] + _a11[_qp]*4.0*_u[_qp]*_u[_qp]*_u[_qp]
      + _a12[_qp]*2.0*_u[_qp]*(square_sum);

  case Jacobian:
    return _phi[_j][_qp]*(2.0*_a1[_qp] + 12.0*_a11[_qp]*_u[_qp]*_u[_qp]
                          + 2.0*_a12[_qp]*square_sum);
  }

  mooseError("Invalid type passed in");
}
