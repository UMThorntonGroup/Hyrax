/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  17 December 2013
*
*************************************************************************/

#ifndef PHASEFRACTIONDT_H
#define PHASEFRACTIONDT_H

#include "TimeStepper.h"
#include "PostprocessorInterface.h"

class PhaseFractionDT;

template<>
InputParameters validParams<PhaseFractionDT>();

class PhaseFractionDT :
  public TimeStepper,
  public PostprocessorInterface
{
public:
  PhaseFractionDT(const InputParameters & parameters);

protected:
  virtual Real computeInitialDT();
  virtual Real computeDT();

  Real _main_dt;

  Real _initial_dt;
  bool _has_initial_dt;
  unsigned int _n_initial_steps;

  Real _growth_factor;

  const PostprocessorValue & _volfrac;
  const PostprocessorValue & _volfrac_old;

private:
};

#endif //PHASEFRACTIONDT_H
