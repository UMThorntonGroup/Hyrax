/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  14 October 2014
*
*************************************************************************/

#ifndef SIMPLECONSTANTDIFFUSION_H
#define SIMPLECONSTANTDIFFUSION_H

#include "Diffusion.h"

//forward declarations
class SimpleConstantDiffusion;

template<>
InputParameters validParams<SimpleConstantDiffusion>();

class SimpleConstantDiffusion : public Diffusion
{
public:
  SimpleConstantDiffusion(const InputParameters & parameters);
  virtual ~SimpleConstantDiffusion();

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();

private:
  Real _k;

};

#endif //SIMPLECONSTANTDIFFUSION_H
