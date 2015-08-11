/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  21 April 2015
*
*************************************************************************/

#ifndef CALPHADAB1CD2MODMATERIAL_H
#define CALPHADAB1CD2MODMATERIAL_H

#include "CalphadAB1CD2Material.h"
#include "CalphadAB1CD2.h"

class CalphadAB1CD2ModMaterial;

template<>
InputParameters validParams<CalphadAB1CD2ModMaterial>();

class CalphadAB1CD2ModMaterial : public CalphadAB1CD2Material
{
public:
  CalphadAB1CD2ModMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpProperties();
  virtual void computeModFunction();

  /* virtual Real calculateReference(Real c);
  virtual Real calculateIdeal(Real c);
  virtual Real calculateExcess(Real c);

  virtual Real calculateSecondLatticeGminusHser();

  virtual Real computeGMix(Real c);
  virtual Real computeDGMixDc(Real c);
  virtual Real computeD2GMixDc2();
  virtual Real computeD3GMixDc3();*/

private:

  /*
  CalphadAB1CD2 _energy;

  Real _low_cutoff;
  Real _high_cutoff;
  Real _precip_conc;

  std::vector<double> _pure_EP1_phase1_coeffs;
  MaterialProperty<Real> & _G_AB1CD2;
  MaterialProperty<Real> & _dG_dc;
  MaterialProperty<Real> & _d2G_dc2;
  //MaterialProperty<Real> & _d3G_dc3;
  MaterialProperty<Real> & _G_AB1CD2_precip;
  MaterialProperty<Real> & _d2G_dc2_precip;
  */

  Real _mod_cutoff;
  std::vector<Real> _A1;
  std::vector<Real> _A2;
  Real _A1_to_0;
  Real _A2_to_0;

  Real _mod;
  Real _dmoddc;
  Real _d2moddc2;

};

#endif //CALPHADAB1CD2MODMATERIAL_H
