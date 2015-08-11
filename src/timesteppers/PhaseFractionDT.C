/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  17 December 2013
*
*************************************************************************/

#include "PhaseFractionDT.h"

template<>
InputParameters validParams<PhaseFractionDT>()
{
    InputParameters params = validParams<TimeStepper>();

  params.addParam<Real>("initial_dt", -1, "initial tiny timestep size to start simulation");
  params.addParam<bool>("has_initial_dt",false, "boolean for if you have tiny initial timestep");
  params.addParam<unsigned int>("n_initial_steps", 1, "how many initial tiny timesteps to take");

  params.addParam<Real>("growth_factor", 0.1, "value to increase dt by (0 to 1)");

  params.addRequiredParam<PostprocessorName>("postprocessor", "Name of the postprocessor timestepper uses to compute the dt");

  params.addRequiredParam<Real>("dt", "main initial timestep size to use with this timestepper");

  return params;
}

PhaseFractionDT::PhaseFractionDT(const InputParameters & parameters) :
    TimeStepper(parameters),
    PostprocessorInterface(parameters),
    _main_dt(getParam<Real>("dt")),
    _initial_dt(getParam<Real>("initial_dt")),
    _has_initial_dt(getParam<bool>("has_initial_dt")),
    _n_initial_steps(getParam<unsigned int>("n_initial_steps")),
    _growth_factor(getParam<Real>("growth_factor")),
    _volfrac(getPostprocessorValue("postprocessor")),
    _volfrac_old(getPostprocessorValueOld("postprocessor"))
{
  if(_growth_factor < 0 || _growth_factor > 1)
    mooseError("please set growth factor between 0 and 1 (PhaseFractionDT)");
}

Real
PhaseFractionDT::computeInitialDT()
{
  if(_has_initial_dt)
    return _initial_dt;
  else
    return _main_dt;
}

Real
PhaseFractionDT::computeDT()
{
  if (_t_step <= _n_initial_steps && _has_initial_dt)
    return _initial_dt;

  else if (_t_step == _n_initial_steps + 1 && _has_initial_dt)
    return _main_dt;

  else
  {
    Real volfrac_change = (_volfrac - _volfrac_old)/_volfrac_old;

    if (volfrac_change < 0)
      return _dt;

    else
      return _dt + _dt*(1-volfrac_change)*_growth_factor;

  }
}
