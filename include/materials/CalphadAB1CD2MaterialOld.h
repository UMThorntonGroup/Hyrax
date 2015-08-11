/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  7 November 2013
*
*************************************************************************/

#ifndef CALPHADAB1CD2MATERIALOLD_H
#define CALPHADAB1CD2MATERIALOLD_H

#include "CalphadEnergyMaterial.h"
#include "CalphadAB1CD2Old.h"

class CalphadAB1CD2MaterialOld;

template<>
InputParameters validParams<CalphadAB1CD2MaterialOld>();

class CalphadAB1CD2MaterialOld : public CalphadEnergyMaterial
{
public:
  CalphadAB1CD2MaterialOld(const InputParameters & parameters);

protected:
  virtual void computeQpProperties();

  /* virtual Real calculateReference(Real c);
  virtual Real calculateIdeal(Real c);
  virtual Real calculateExcess(Real c);

  virtual Real calculateSecondLatticeGminusHser();

  virtual Real computeGMix(Real c);
  virtual Real computeDGMixDc(Real c);
  virtual Real computeD2GMixDc2();
  virtual Real computeD3GMixDc3();*/

private:

  CalphadAB1CD2Old _energy;

  Real _low_cutoff;
  Real _high_cutoff;
  Real _precip_conc;

  std::vector<double> _pure_EP1_phase1_coeffs;
  MaterialProperty<Real> & _G_AB1CD2;
  MaterialProperty<Real> & _dG_dc;
  MaterialProperty<Real> & _d2G_dc2;
  //MaterialProperty<Real> & _d3G_dc3;
  MaterialProperty<Real> & _G_AB1CD2_precip;
};

#endif //CALPHADAB1CD2MATERIALOLD_H
