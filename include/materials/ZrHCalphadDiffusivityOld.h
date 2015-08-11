/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  13 February 2014
*
*************************************************************************/

#ifndef ZRHCALPHADDIFFUSIVITYOLD_H
#define ZRHCALPHADDIFFUSIVITYOLD_H

#include "ZrHCalphad.h"

//forward declaration
class ZrHCalphadDiffusivityOld;

template<>
InputParameters validParams<ZrHCalphad>();

class ZrHCalphadDiffusivityOld : public ZrHCalphad
{
public:
  ZrHCalphadDiffusivityOld(const InputParameters & parameters);

protected:
  virtual void computeQpProperties();
  Real computeHeaviside();


private:
  //Diffusion coefficient information
  Real _H_Zr_D0;
  Real _H_ZrH2_D0;
  Real _H_Zr_Q0;
  Real _H_ZrH2_Q0;
  Real _R;
  Real _k;
  Real _mobility_CH_scaling;

  const MaterialProperty<Real> & _d2Galpha_dc2;
  const MaterialProperty<Real> & _d2Gdelta_dc2;

  MaterialProperty<Real> & _D_alpha;
  MaterialProperty<Real> & _D_delta;

  unsigned int _n_OP_variables;
  std::vector<VariableValue *> _OP;

  VariableValue & _c;

  MaterialProperty<Real> & _L1Q;
  Real _Q_transport;
  const MaterialProperty<Real> & _d2Galpha_dcdT;

};

#endif //ZRHCALPHADDIFFUSIVITYOLD_H
