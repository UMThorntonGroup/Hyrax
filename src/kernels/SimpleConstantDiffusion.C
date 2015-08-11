/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  14 October 2014
*
*************************************************************************/

#include "SimpleConstantDiffusion.h"

template<>
InputParameters validParams<SimpleConstantDiffusion>()
{
  InputParameters params = validParams<Diffusion>();

  params.addRequiredParam<Real>("k", "constant coefficient for the diffusion equation");

  return params;
}

SimpleConstantDiffusion::SimpleConstantDiffusion(const InputParameters & parameters)
    :Diffusion(parameters),
     _k(getParam<Real>("k"))
{
}

SimpleConstantDiffusion::~SimpleConstantDiffusion()
{
}

Real
SimpleConstantDiffusion::computeQpResidual()
{
  return _k*Diffusion::computeQpResidual();
}

Real
SimpleConstantDiffusion::computeQpJacobian()
{
  return _k*Diffusion::computeQpJacobian();
}
