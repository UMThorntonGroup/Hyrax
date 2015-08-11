/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  7 November 2013
*
*************************************************************************/

#include "CalphadAB1CD2MaterialOld.h"

template<>
InputParameters validParams<CalphadAB1CD2MaterialOld>()
{
  InputParameters params = validParams<CalphadEnergyMaterial>();

  params.addRequiredParam<std::vector<Real> >("pure_EP1_phase1_coeffs", "coeffs of pure endpoint at low composition in the first phase");

  params.addParam<Real>("low_cutoff", 0.001, "linearization cutoff, low end");
  params.addParam<Real>("high_cutoff", 0.655, "linearization cutoff, high end");
  params.addParam<Real>("precip_conc", 0.6, "concentration of precipitate");

  return params;
}

CalphadAB1CD2MaterialOld::CalphadAB1CD2MaterialOld(const InputParameters & parameters) :
    CalphadEnergyMaterial(parameters),
    _energy(),
    _low_cutoff(getParam<Real>("low_cutoff")),
    _high_cutoff(getParam<Real>("high_cutoff")),
    _precip_conc(getParam<Real>("precip_conc")),
    _pure_EP1_phase1_coeffs(getParam<std::vector<Real> >("pure_EP1_phase1_coeffs")),
    _G_AB1CD2(declareProperty<Real>("G_AB1CD2")),
    _dG_dc(declareProperty<Real>("dGAB1CD2_dc")),
    _d2G_dc2(declareProperty<Real>("d2GAB1CD2_dc2")),
    // _d3G_dc3(declareProperty<Real>("d3GAB1CD2_dc3"))
    _G_AB1CD2_precip(declareProperty<Real>("G_AB1CD2_precip"))
{
  _energy.parameterize(_R, _pure_endpoint_low_coeffs, _pure_endpoint_high_coeffs, _mixture_coeffs,
                       _L0_coeffs, _L1_coeffs, _pure_EP1_phase1_coeffs);
}

void
CalphadAB1CD2MaterialOld::computeQpProperties()
{
  Real c;

  if (_c[_qp]  < _low_cutoff)
    c = _low_cutoff;
  else if (_c[_qp] > _high_cutoff)
    c = _high_cutoff;
  else c = _c[_qp];

  /* _G_AB1CD2[_qp] = computeGMix(c);
  _dG_dc[_qp] = computeDGMixDc(c);
  _d2G_dc2[_qp] = computeD2GMixDc2();
  _d3G_dc3[_qp] = computeD3GMixDc3(); */


  _G_AB1CD2[_qp] = _energy.computeGMix(c, _T[_qp]);
  _dG_dc[_qp] = _energy.computeDGMixDc(c, _T[_qp]);
  _d2G_dc2[_qp] = _energy.computeD2GMixDc2(c, _T[_qp]);

  _G_AB1CD2_precip[_qp] = _energy.computeGMix(_precip_conc, _T[_qp]);
}
/*
Real
CalphadAB1CD2Material::calculateReference(Real c)
{
  Real first_energy = calculateFirstLatticeGminusHser();
  Real second_energy = calculateSecondLatticeGminusHser();

  return  0.5*( (2 - 3*c)*first_energy + c*second_energy );
}


Real
CalphadAB1CD2Material::calculateIdeal(Real c)
{
   return _R*_T[_qp]*( c*std::log(c/(2-2*c)) + (2-3*c)*std::log((2-3*c)/(2-2*c)) );
}

Real
CalphadAB1CD2Material::calculateExcess(Real c)
{
  Real L0 = _L0_coeffs[0] + _L0_coeffs[1]*_T[_qp];
  Real L1 = _L1_coeffs[0] + _L1_coeffs[1]*_T[_qp];

  Real polynomial = (c-1)*L0 + (1-2*c)*L1;

  Real prefactor = 3*c*c - 2*c;
  prefactor /= 4*(c-1)*(c-1);

  return prefactor*polynomial;
}

Real
CalphadAB1CD2Material::calculateSecondLatticeGminusHser()
{
  Real first_term = _mixture_coeffs[0]
    + _mixture_coeffs[1]*_T[_qp]
    + _mixture_coeffs[2]*_T[_qp]*std::log(_T[_qp]);

  Real second_term = _pure_EP1_phase1_coeffs[0]
    + _pure_EP1_phase1_coeffs[1]*_T[_qp]
    + _pure_EP1_phase1_coeffs[2]*_T[_qp]*std::log(_T[_qp])
    + _pure_EP1_phase1_coeffs[3]*_T[_qp]*_T[_qp]
    + _pure_EP1_phase1_coeffs[4]/_T[_qp];

  Real third_term = _pure_endpoint_high_coeffs[0]
    + _pure_endpoint_high_coeffs[1]*_T[_qp]
    + _pure_endpoint_high_coeffs[2]*_T[_qp]*std::log(_T[_qp])
    + _pure_endpoint_high_coeffs[3]*_T[_qp]*_T[_qp]
    + _pure_endpoint_high_coeffs[4]/_T[_qp];

 return first_term + second_term + third_term;
}

Real
CalphadAB1CD2Material::computeGMix(Real c)
{
  Real c1;

  if( _c[_qp] < 0.001)
  {
    c1 = 0.001;

    return CalphadEnergyMaterial::computeGMix(c1) + computeDGMixDc(c1)*(_c[_qp] - c1);
  }
  else if (_c[_qp] > 0.655)
  {
    c1 = 0.655;

    return CalphadEnergyMaterial::computeGMix(c1) + computeDGMixDc(c1)*(_c[_qp] - c1);
  }
  else
    return CalphadEnergyMaterial::computeGMix(_c[_qp]);
}

Real
CalphadAB1CD2Material::computeDGMixDc(Real c)
{
  Real ref;
  Real ideal;
  ref = -1.5*calculateFirstLatticeGminusHser() + 0.5*calculateSecondLatticeGminusHser();

  ideal = _R*_T[_qp]*( std::log(4*c/(1-c)) - 3*std::log((2-3*c)/(1-c)) );

  Real xs;

  Real L0 = _L0_coeffs[0] + _L0_coeffs[1]*_T[_qp];
  Real L1 = _L1_coeffs[0] + _L1_coeffs[1]*_T[_qp];

  Real polynomial = (3*c*c*c -9*c*c + 8*c -2)*L0 + (-6*c*c*c + 18*c*c -12*c +2)*L1;

  Real prefactor = 1/(4*(c-1)*(c-1)*(c-1) );

  xs = prefactor*polynomial;

  return ref + ideal + xs;
}

Real
CalphadAB1CD2Material::computeD2GMixDc2()
{
  Real ref;
  Real ideal;
  Real xs;

  //make this piecewise in concentration space
  if( _c[_qp] < 0.001)
    return 0;

  else if (_c[_qp] > 0.655)
    return 0;

  else
  {
    ref = 0;

    ideal = _R*_T[_qp]*( 2/( _c[_qp]*(3*_c[_qp]*_c[_qp] - 5*_c[_qp] + 2)) );

    Real L0 = _L0_coeffs[0] + _L0_coeffs[1]*_T[_qp];
    Real L1 = _L1_coeffs[0] + _L1_coeffs[1]*_T[_qp];

    Real polynomial = (_c[_qp] - 1)*L0 + (3 - 6*_c[_qp])*L1;

    Real prefactor = 0.5/(4*std::pow(_c[_qp]-1, 4));

    xs = prefactor*polynomial;

    return ref + ideal + xs;
  }

}

Real
CalphadAB1CD2Material::computeD3GMixDc3()
{
  Real ref;
  Real ideal;
  Real xs;

  //make this piecewise in concentration space
  if( _c[_qp] < 0.001)
    return 0;

  else if (_c[_qp] > 0.655)
    return 0;

  else
  {
    ref = 0;

    ideal = -1*_R*_T[_qp]*( (18*_c[_qp]*_c[_qp] - 20*_c[_qp] + 4)
                            / std::pow( _c[_qp]*(3*_c[_qp]*_c[_qp] - 5*_c[_qp] + 2), 2) );

    Real L0 = _L0_coeffs[0] + _L0_coeffs[1]*_T[_qp];
    Real L1 = _L1_coeffs[0] + _L1_coeffs[1]*_T[_qp];

    Real polynomial = (1 - _c[_qp])*L0 + (6*_c[_qp] - 2)*L1;

    Real prefactor = 1.5/(std::pow(_c[_qp]-1, 5));

    xs = prefactor*polynomial;

    return ref + ideal + xs;
  }

}

*/
