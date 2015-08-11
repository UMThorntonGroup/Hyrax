/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*************************************************************************/

#include "AuxNucleationRate.h"

#include <cmath>

template<>
InputParameters validParams<AuxNucleationRate>()
{
  InputParameters params = validParams<AuxKernel>();

  params.addRequiredCoupledVar("coupled_aux_var","coupled auxiliary variable for free energy change");
  // I'd have line in input file, coupled_aux_var = supersaturation, or chem_elastic, or whatever

  params.addRequiredParam<Real>("gamma", "Surface energy");
  params.addParam<Real>("Kb", 1.3806503e-23, "Boltzmann's constant, make sure units all match");
  params.addParam<Real>("temperature", 473, "Temperature");
  params.addRequiredParam<Real>("scale_factor","factor to scale energy/dimensions by");

  params.addRequiredParam<Real>("Z", "Non-equilibrium Zeldovitch factor");
  params.addRequiredParam<Real>("Beta_star", "1/characteristic nucleation time");
  params.addRequiredParam<Real>("linear_density", "linear atomic density of matrix");

  return params;
}

AuxNucleationRate::AuxNucleationRate(const InputParameters & parameters)
  : AuxKernel(parameters),
    _mesh_dimension(_mesh.dimension()),
    _coupled_energy(coupledValue("coupled_aux_var")),
    _Z(getParam<Real>("Z")),
   _beta_star(getParam<Real>("Beta_star")),
    _linear_density(getParam<Real>("linear_density")),
    _gamma(getParam<Real>("gamma")),
    _Kb(getParam<Real>("Kb")),
    _temperature(getParam<Real>("temperature")),
    _scale_factor(getParam<Real>("scale_factor"))
{
}

Real
AuxNucleationRate::computeValue()
{
  Real kn1;
  Real kn2;
  Real alpha;

  // handling the dimension of the problem here and making sure we get the correct
  // (areal or volume) density
  if (_mesh_dimension == 2)
  {
    kn1 = _Z*_beta_star*pow(_linear_density, 2);
    alpha = libMesh::pi;
  }
  else if (_mesh_dimension == 3)
  {
    kn1 = _Z*_beta_star*pow(_linear_density, 3);
    alpha = (16*libMesh::pi)/3;
  }
  else
    mooseError("honky, your problem dimesion must be 2 or 3 (AuxNucleationRate");

  // correct the density to the actual element volume to get # of atoms
  kn1 *= _current_elem_volume;

  kn2 = _scale_factor*alpha*std::pow(_gamma, (int)_mesh_dimension)/(_Kb*_temperature);

  return kn1*std::exp(-1.0*kn2/std::pow(_coupled_energy[_qp], (int)_mesh_dimension-1));
}

