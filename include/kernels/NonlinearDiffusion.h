/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  14 October 2014
*
*************************************************************************/

#ifndef NONLINEARDIFFUSION_H
#define NONLINEARDIFFUSION_H

#include "Diffusion.h"

//forward declarations
class NonlinearDiffusion;

template<>
InputParameters validParams<NonlinearDiffusion>();

class NonlinearDiffusion : public Diffusion
{
public:
  NonlinearDiffusion(const InputParameters & parameters);
  virtual ~NonlinearDiffusion();

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();

private:
  Real _k;

};

#endif //NONLINEARDIFFUSION_H
