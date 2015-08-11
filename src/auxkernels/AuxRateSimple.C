/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  18 April 2013
*
*************************************************************************/

#include "AuxRateSimple.h"

#include <cmath>
#include <ostream>

template<>
InputParameters validParams<AuxRateSimple>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredParam<Real>("Kn1", "First nucleation rate coefficient");
  params.addRequiredParam<Real>("Kn2", "Second nucleation rate coefficient");

  params.addRequiredCoupledVar("coupled_aux_var","coupled auxiliary variable for free energy change");
  // I'd have line in input file, coupled_aux_var = supersaturation, or chem_elastic, or whatever

  params.addRequiredParam<int>("n_OP_vars", "# of coupled OP variables");
  params.addRequiredCoupledVar("OP_var_names", "Array of coupled OP variable names");

  return params;
}

AuxRateSimple::AuxRateSimple(const InputParameters & parameters)
  : AuxKernel(parameters),
    _mesh_dimension(_mesh.dimension()),
    _coupled_energy(coupledValue("coupled_aux_var")),
    _Kn1(getParam<Real>("Kn1")),
    _Kn2(getParam<Real>("Kn2")),
    _n_OP_vars(getParam<int>("n_OP_vars"))
{
  if(_n_OP_vars != coupledComponents("OP_var_names"))
    mooseError("Please match the number of orientation variants to coupled OPs (AuxRateSimple).");

  _coupled_OP_vars.resize(_n_OP_vars);

  for(unsigned int i=0; i< _n_OP_vars; i++)
    _coupled_OP_vars[i] = &coupledValue("OP_var_names", i);
}

Real
AuxRateSimple::computeValue()
{
// NOTE: AuxRateSimple WORKS ONLY FOR NON-ADAPTIVE MESHES!
  if ( _mesh.changed() && _t_step > 1)
    mooseError("AuxRateSimple cannot be used with adaptive meshing.");

  // _console<<"exponential term"<< exp(-1.0*kn2/std::pow(_coupled_energy[_qp], (int)_dim-1)) <<std::endl;
  return _Kn1*std::exp(-1.0*_Kn2/std::pow(_coupled_energy[_qp], (int)_mesh_dimension-1));
}

