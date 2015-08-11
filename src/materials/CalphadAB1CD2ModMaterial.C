/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  21 April 2015
*
*************************************************************************/

#include "CalphadAB1CD2ModMaterial.h"

template<>
InputParameters validParams<CalphadAB1CD2ModMaterial>()
{
  InputParameters params = validParams<CalphadAB1CD2Material>();

//  params.addRequiredParam<std::vector<Real> >("pure_EP1_phase1_coeffs", "coeffs of pure endpoint at low composition in the first phase");

//  params.addParam<Real>("low_cutoff", 0.001, "linearization cutoff, low end");
//  params.addParam<Real>("high_cutoff", 0.655, "linearization cutoff, high end");
//  params.addParam<Real>("precip_conc", 0.6, "concentration of precipitate");
  params.addRequiredParam<Real>("mod_cutoff", "composition for mod function cutoff");
  params.addRequiredParam<std::vector<Real> >("A1", "vector of coefficients, Aa1*T + Aa2 for double-well mod function component");
  params.addRequiredParam<std::vector<Real> >("A2", "vector of coefficients, Aa1*T + Aa2 for linear mod function component");
  params.addParam<Real>("A1_to_0", 0, "Temperature at which A1 goes to zero");
  params.addParam<Real>("A2_to_0", 0, "Temperature at which A2 goes to zero");

  return params;
}

CalphadAB1CD2ModMaterial::CalphadAB1CD2ModMaterial(const InputParameters & parameters) :
    CalphadAB1CD2Material(parameters),
    //_energy(),
    //_low_cutoff(getParam<Real>("low_cutoff")),
    //_high_cutoff(getParam<Real>("high_cutoff")),
    //_precip_conc(getParam<Real>("precip_conc")),
    //_pure_EP1_phase1_coeffs(getParam<std::vector<Real> >("pure_EP1_phase1_coeffs")),
    //_G_AB1CD2(declareProperty<Real>("G_AB1CD2")),
    //_dG_dc(declareProperty<Real>("dGAB1CD2_dc")),
    //_d2G_dc2(declareProperty<Real>("d2GAB1CD2_dc2")),
    //// _d3G_dc3(declareProperty<Real>("d3GAB1CD2_dc3"))
    //_G_AB1CD2_precip(declareProperty<Real>("G_AB1CD2_precip")),
    //_d2G_dc2_precip(declareProperty<Real>("d2GAB1CD2_dc2_precip"))
    _mod_cutoff(getParam<Real>("mod_cutoff")),
    _A1(getParam<std::vector<Real> >("A1")),
    _A2(getParam<std::vector<Real> >("A2")),
    _A1_to_0(getParam<Real>("A1_to_0")),
    _A2_to_0(getParam<Real>("A2_to_0"))
{
  //error check
  if (_A1.size() != 2)
    mooseError("A1 vector needs to be size of 2 (CalphadAB1CD2ModMaterial");

  if (_A2.size() != 2)
    mooseError("A2 vector needs to be size of 2 (CalphadAB1CD2ModMaterial");

  _energy.parameterize(_R, _low_cutoff, _high_cutoff, _pure_endpoint_low_coeffs, _pure_endpoint_high_coeffs, _mixture_coeffs,
                       _L0_coeffs, _L1_coeffs, _pure_EP1_phase1_coeffs);
}

void
CalphadAB1CD2ModMaterial::computeQpProperties()
{
  //Real c;

  /*
  if (_c[_qp]  < _low_cutoff)
    c = _low_cutoff;
  else if (_c[_qp] > _high_cutoff)
    c = _high_cutoff;
  else c = _c[_qp];
  */
  computeModFunction();

  _G_AB1CD2[_qp] = _energy.computeGMix(_c[_qp], _T[_qp]) + _mod;
  _dG_dc[_qp] = _energy.computeDGMixDc(_c[_qp], _T[_qp]) + _dmoddc;
  _d2G_dc2[_qp] = _energy.computeD2GMixDc2(_c[_qp], _T[_qp]) + _d2moddc2;

  //this ASSUMES that the mod function concentration cutoff is below the precip
  //concentration used here.  DON'T Screw it up!
  _G_AB1CD2_precip[_qp] = _energy.computeGMix(_precip_conc, _T[_qp]);
  _d2G_dc2_precip[_qp] = _energy.computeD2GMixDc2(_precip_conc, _T[_qp]);
}

void
CalphadAB1CD2ModMaterial::computeModFunction()
{
  Real A1(0);
  Real A2(0);

  if (_c[_qp] >= 0.57)
  {
    _mod = 0;
    _dmoddc = 0;
    _d2moddc2 = 0;
  }
  else
  {
    if (_T[_qp] < _A1_to_0)
      A1 = _A1[0]*_T[_qp] + _A1[1];
    else
      A1 = 0;

    if (_T[_qp] < _A2_to_0)
      A2 = _A2[0]*_T[_qp] + _A2[1];
    else
      A2 = 0;

    Real g(0);
    Real dgdc(0);
    Real d2gdc2(0);

    g = std::pow( ((1/_mod_cutoff)*_c[_qp]), 2.0)
      - 2*std::pow( ((1/_mod_cutoff)*_c[_qp]), 3.0)
      + std::pow( ((1/_mod_cutoff)*_c[_qp]), 4.0);

    dgdc = (2/_mod_cutoff)*(_c[_qp]/_mod_cutoff)
      - (6/_mod_cutoff)*(std::pow( ((1/_mod_cutoff)*_c[_qp]), 2.0))
      + (4/_mod_cutoff)*(std::pow( ((1/_mod_cutoff)*_c[_qp]), 3.0));

    d2gdc2 = (2/_mod_cutoff)*(1/_mod_cutoff)
      - (12/_mod_cutoff)*(1/_mod_cutoff)*((1/_mod_cutoff)*_c[_qp])
      + (12/_mod_cutoff)*(1/_mod_cutoff)*(std::pow( ((1/_mod_cutoff)*_c[_qp]), 2.0));

    Real l(0);
    Real dldc(0);
    Real d2ldc2(0);

    l = _mod_cutoff - _c[_qp];
    dldc = -1;
    d2ldc2 = 0;

    _mod = A1*g + A2*l;
    _dmoddc = A1*dgdc + A2*dldc;
    _d2moddc2 = A1*d2gdc2 + A2*d2ldc2;
  }

}
