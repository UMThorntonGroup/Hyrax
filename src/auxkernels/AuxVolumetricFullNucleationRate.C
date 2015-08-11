/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
* 13 January 2015
*
*************************************************************************/

#include "AuxVolumetricFullNucleationRate.h"

#include <cmath>
#include <iostream>

//NOTE: THIS ASSUMES ELEMENTS ARE SQUARE OR CUBIC!

template<>
InputParameters validParams<AuxVolumetricFullNucleationRate>()
{
  InputParameters params = validParams<AuxFullNucleationRate>();
  params.addRequiredParam<Real>("rate_volume", "volume in m^3 for which the nucleation rate is calculated");
  /*
  params.addRequiredCoupledVar("coupled_bulk_energy_change","coupled auxiliary variable for free energy change");
  // I'd have line in input file, coupled_aux_var = supersaturation, or chem_elastic, or whatever

  params.addRequiredParam<Real>("gamma", "Surface energy");
  params.addParam<Real>("Kb", 1.3806503e-23, "Boltzmann's constant, make sure units all match");
  //params.addParam<Real>("temperature", 473, "Temperature");
  params.addParam<Real>("scale_factor", 1, "factor to scale nucleation rate by");

  //params.addRequiredParam<Real>("Z", "Non-equilibrium Zeldovitch factor");
  //params.addRequiredParam<Real>("Beta_star", "1/characteristic nucleation time");
  params.addRequiredParam<Real>("linear_density", "linear atomic density of matrix");

  params.addRequiredCoupledVar("T", "temperature");
  params.addRequiredCoupledVar("X", "atomic fraction of solute");

  params.addRequiredParam<Real>("jump_distance", "atomic distance");
  params.addRequiredParam<int>("n_OP_variables", "# of coupled OP variables, >=1");
  params.addRequiredCoupledVar("OP_variable_names", "Array of coupled OP variable names");
  params.addParam<Real>("OP_threshold", 0.0001, "threshold below which not to calculate nucleation rate");
  params.addParam<Real>("length_scale_factor", "characteristic length the simulation is scaled by");
*/
  return params;
}

AuxVolumetricFullNucleationRate::AuxVolumetricFullNucleationRate(const InputParameters & parameters)
    : AuxFullNucleationRate(parameters),
      _rate_volume(getParam<Real>("rate_volume"))
      /*
      _mesh_dimension(_mesh.dimension()),
      _coupled_energy(coupledValue("coupled_bulk_energy_change")),
      _Z(),
      _beta_star(),
      _linear_density(getParam<Real>("linear_density")),
      _gamma(getParam<Real>("gamma")),
      _Kb(getParam<Real>("Kb")),
      _scale_factor(getParam<Real>("scale_factor")),
      _n_OP_variables(getParam<int>("n_OP_variables")),
      _T(coupledValue("T")),
      _X(coupledValue("X")),
      _D(getMaterialProperty<Real>("D_alpha")),
      _jump_distance(getParam<Real>("jump_distance")),
      _Omega(getMaterialProperty<Real>("molar_volume")),
      _OP_threshold(getParam<Real>("OP_threshold")),
      _length_scale(getParam<Real>("length_scale_factor"))
      */
{
  // Create a vector of the coupled OP variables and gradients
  if(_n_OP_variables != coupledComponents("OP_variable_names"))
    mooseError("Please match the # of orientation variants to coupled OPs (CHCoupledCalphad)");

  _OP.resize(_n_OP_variables);

  for(unsigned int i=0; i< _n_OP_variables; i++)
    _OP[i] = &coupledValue("OP_variable_names", i);
}

/*
Real
AuxVolumetricFullNucleationRate::computeValue()
{

  computeCriticalRadius();
  computeCriticalEnergy();
  computeZeldovichFactor();
  computeCriticalFrequency();
  computeNumAtoms();

  _console<<"scale_factor = "<<_scale_factor<<std::endl;
  _console<<"Z = "<<_Z<<std::endl;
  _console<<"N = "<<_N<<std::endl;
  _console<<"beta_star = "<<_beta_star<<std::endl;
  _console<<"G_star = "<<_G_star<<std::endl;
  _console<<"Kb = "<<_Kb<<std::endl;
  _console<<"T = "<<_T[_qp]<<std::endl;
  _console<<"exp(-Gstar/kT) = "<<std::exp( (-1*_G_star)/ (_Kb*_T[_qp]) )<<std::endl;

  Real rate = _scale_factor*( _Z*_N*_beta_star*std::exp( (-1*_G_star)/ (_Kb*_T[_qp]) ));

  if (rate < 0)
    return 0;
  else
    return rate;

    }*/

/*void
AuxFullNucleationRate::computeCriticalRadius()
{
  //this is calculated as if in 3D
  _r_star = 2*_gamma/_coupled_energy[_qp];

  // _console<<"r* = "<<_r_star<<std::endl;
  }*/

/*void
AuxFullNucleationRate::computeCriticalEnergy()
{
  //this is calculated as if in 3D
  Real alpha = 16*libMesh::pi/3;

  _G_star = alpha*std::pow(_gamma, 3)/std::pow(_coupled_energy[_qp], 2);

  // _console<<"G* in Joules = "<<_G_star<<std::endl;
  }*/


/*void
AuxFullNucleationRate::computeZeldovichFactor()
{
  Real critical_volume = 4*libMesh::pi*_r_star*_r_star*_r_star/3;

  //here I'm assuming that the molar volume is constant across phases.. meh
  Real Nc = critical_volume*6.02214E23/_Omega[_qp];

  _Z =std::sqrt( _G_star/( 3*libMesh::pi*Nc*Nc*_Kb*_T[_qp] ));

  //_console<<"Z = "<<_Z<<std::endl;
  }*/

/*void
AuxVolumetricFullNucleationRate::computeCriticalFrequency()
{
  Real Zc = 4*libMesh::pi*_r_star*_r_star*_linear_density*_linear_density;

  _console<<"Zc = "<<Zc<<std::endl;
  _console<<"X = "<<_X[_qp]<<std::endl;
  _console<<"D = "<<_D[_qp]<<std::endl;
  _console<<"jump_dist^2 = "<< ( std::pow(_jump_distance,2))<<std::endl;

  _beta_star = Zc*_X[_qp]*_D[_qp]/ ( std::pow(_jump_distance,2)) ;

  // _console<<"beta* = "<<_beta_star<<std::endl;
  }*/

void
AuxVolumetricFullNucleationRate::computeNumAtoms()
{
 //correct the density to the actual element volume to get # of atoms

 Real atomic_volume = _Omega[_qp]/6.02214E23;

 //this still assumes all the nucleation rate computation is 3D
 //this also assumes the elements are square.  MAJOR ASSUMPTION.

/*
 Real vol;

 //this volume is non-dimensionalized
 if (_mesh_dimension == 2)
   vol = std::sqrt(_current_elem_volume);
 else if (_mesh_dimension == 3 )
   vol = std::pow(_current_elem_volume, (1./3.));
 else
   mooseError("please perform this simulation in 2D or 3D (AuxFullNucleationRate)");

//this volume is dimensionalized  by the simulation length scale
 vol = vol*vol*vol*_length_scale*_length_scale*_length_scale;

 //_console<<"atomic volume = "<<atomic_volume<<std::endl;
 //_console<<"volume = "<<vol<<std::endl;

 _N = vol/atomic_volume;
*/
 _N = _rate_volume/atomic_volume;

}
