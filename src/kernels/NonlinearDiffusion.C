/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  14 October 2014
*
*************************************************************************/

#include "NonlinearDiffusion.h"

template<>
InputParameters validParams<NonlinearDiffusion>()
{
  InputParameters params = validParams<Diffusion>();

  params.addRequiredParam<Real>("k", "constant coefficient for the diffusion equation");

  return params;
}

NonlinearDiffusion::NonlinearDiffusion(const InputParameters & parameters)
    :Diffusion(parameters),
     _k(getParam<Real>("k"))
{
}

NonlinearDiffusion::~NonlinearDiffusion()
{
}

Real
NonlinearDiffusion::computeQpResidual()
{
  return _grad_test[_i][_qp]*( _k*_grad_u[_qp]/_u[_qp] );
}

Real
NonlinearDiffusion::computeQpJacobian()
{
  //this might be right.
  return _grad_test[_i][_qp]*( _k*_grad_phi[_j][_qp]/_u[_qp]
                               - _k*_grad_u[_qp]/( _phi[_j][_qp]*_u[_qp]*_u[_qp] ) );
}
