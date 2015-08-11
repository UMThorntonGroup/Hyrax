/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  30 August 2013
*
*************************************************************************/

#ifndef FERROELECTRICBULK_H
#define FERROELECTRICBULK_H

#include "Material.h"

//forward declarations
class FerroelectricBulk;

template<>
InputParameters validParams<FerroelectricBulk>();

class FerroelectricBulk : public Material
{
public:
  FerroelectricBulk(const InputParameters & parameters);

protected:
  virtual void computeQpProperties();

private:
  Real _a1_i;
  Real _a11_i;
  Real _a12_i;

  Real _L_i;
  Real _g_i;

//Allen-Cahn equation
  MaterialProperty<Real> & _a1;
  MaterialProperty<Real> & _a11;
  MaterialProperty<Real> & _a12;

  MaterialProperty<Real> & _L;
  MaterialProperty<Real> & _g;
};


#endif //FERROELECTRICBULK_H
