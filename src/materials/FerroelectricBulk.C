/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  30 August 2013
*
*************************************************************************/

#include "FerroelectricBulk.h"

template<>
InputParameters validParams<FerroelectricBulk>()
{
  InputParameters params = validParams<Material>();

  params.addRequiredParam<Real>("a1", "1st Landau coefficient");
  params.addRequiredParam<Real>("a11", "2nd Landau coefficient");
  params.addRequiredParam<Real>("a12", "3rd Landau coefficient");

  params.addRequiredParam<Real>("L", "Allen-Cahn mobility coefficient");
  params.addRequiredParam<Real>("g", "Gradient energy coefficient");

  return params;
}

FerroelectricBulk::FerroelectricBulk(const InputParameters & parameters) :
    Material(parameters),
    _a1_i(getParam<Real>("a1")),
    _a11_i(getParam<Real>("a11")),
    _a12_i(getParam<Real>("a12")),
    _L_i(getParam<Real>("L")),
    _g_i(getParam<Real>("g")),

    // declare all materials properties
    _a1(declareProperty<Real>("a1")),
    _a11(declareProperty<Real>("a11")),
    _a12(declareProperty<Real>("a12")),
    _L(declareProperty<Real>("L")),
    _g(declareProperty<Real>("g"))
{
}

void
FerroelectricBulk::computeQpProperties()
{
  _L[_qp] = _L_i;
  _g[_qp] = _g_i;

  _a1[_qp] = _a1_i;
  _a11[_qp] = _a11_i;
  _a12[_qp] = _a12_i;
}

