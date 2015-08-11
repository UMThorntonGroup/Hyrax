/*************************************************************************
*************************************************************************/

#include "FreeEnergy.h"

template<>
InputParameters validParams<FreeEnergy>()
{
  InputParameters params = validParams<Material>();

  params.addRequiredCoupledVar("coupled_conc", "Coupled concentration variable name");
  params.addRequiredCoupledVar("coupled_n", "coupled order parameter variable name");
  params.addParam<Real>("length_scale_factor", 1, "length scale factor in simulation");
  params.addParam<Real>("energy_scale_factor", 1, "energy scale factor in simulation");

  return params;
}

FreeEnergy::FreeEnergy(const InputParameters & parameters)
    : Material(parameters),
      _length_scale(getParam<Real>("length_scale_factor")),
      _energy_scale(getParam<Real>("energy_scale_factor")),
      _kappa_c(getMaterialProperty<Real>("kappa_c")),
      _kappa_n(getMaterialProperty<Real>("kappa_n")),
      _W(getMaterialProperty<Real>("well_height")),
      _molar_vol(getMaterialProperty<Real>("molar_volume")),
      _alpha_energy(getMaterialProperty<Real>("G_AB1CD1")),
      _delta_energy(getMaterialProperty<Real>("G_AB1CD2")),
      _c(coupledValue("coupled_conc")),
      _n(coupledValue("coupled_n")),
      _grad_c(coupledGradient("coupled_conc")),
      _grad_n(coupledGradient("coupled_n")),
      _free_energy_density(declareProperty<Real>("free_energy_density")),
      _constant_phases_energy_density(declareProperty<Real>("constant_phase_energy_density"))
{
}

void
FreeEnergy::computeQpProperties()
{

  Real solute = _c[_qp];
  if (solute < 0)
    solute = 0;

  Real OP = _n[_qp];
  if (OP < 0) OP = 0;
  if (OP > 1) OP = 1;

  Real Heaviside = 3*OP*OP - 2*OP*OP*OP;
  Real g = OP*OP*(OP-1)*(OP-1);

  //the energies are in J/mol
  Real fchem = (1 - Heaviside)*_alpha_energy[_qp] + Heaviside*_delta_energy[_qp] + _W[_qp]*g;
  Real fchem_scale_factor = (_length_scale*_length_scale*_length_scale/_molar_vol[_qp]);

  Real grad_scale_factor = _energy_scale*_length_scale*_length_scale*_length_scale;
  Real grad_c_term = 0.5*_kappa_c[_qp]*_grad_c[_qp].size_sq();
  Real grad_n_term = 0.5*_kappa_n[_qp]*_grad_n[_qp].size_sq();

  //convert fchem to J/m^3 and then nondimensionalize it
   _free_energy_density[_qp] = fchem*fchem_scale_factor
     + grad_c_term*grad_scale_factor
     + grad_n_term*grad_scale_factor;

   //  _console<<"fchem = "<<fchem*fchem_scale_factor<<std::endl;
   // _console<<"gradc = "<<grad_c_term*grad_scale_factor<<std::endl;
   //_console<<"gradn = "<<grad_n_term*grad_scale_factor<<std::endl;

  Real alpha_energy_density = 0;
  if(OP < 0.1)
    alpha_energy_density = _alpha_energy[_qp]*fchem_scale_factor;

  Real delta_energy_density = 0;
  if(OP > 0.9)
    delta_energy_density = _delta_energy[_qp]*fchem_scale_factor;

  _constant_phases_energy_density[_qp] = alpha_energy_density + delta_energy_density;
}

