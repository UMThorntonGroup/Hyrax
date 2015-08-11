/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*************************************************************************/

#include "AuxNucleationProbability.h"

/**
 *  AuxNucleationProbability handles the nucleation probability (P_nm) calculation in the concurrent
 *  nucleation and growth algorithm first proposed by J.P. Simmons (2000).
 *  Returns the nucleation probability over the domain.
 */

template<>
InputParameters validParams<AuxNucleationProbability>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredCoupledVar("coupled_aux_var","coupled auxiliary variable: nucleation rate");
  // I'd have line in input file, coupled_aux_var = j_star
  params.addRequiredCoupledVar("coupled_variables", "coupled order parameter variables");
  params.addRequiredParam<int>("n_OP_vars", "# of coupled OP variables");
  params.addParam<Real>("OP_threshold", 0.01, "Threshold value for determining existence within 2nd phase");
  //params.addParam<Real>("P_threshold",1,"Threshold value above which to ignore probability");

  return params;
}

AuxNucleationProbability::AuxNucleationProbability(const InputParameters & parameters)
  : AuxKernel(parameters),
    _coupled_nuc_rate(coupledValue("coupled_aux_var")),
    _n_OP_vars(getParam<int>("n_OP_vars")),
    _OP_threshold(getParam<Real>("OP_threshold"))
    //_P_threshold(getParam<Real>("P_threshold"))
{
  if(_n_OP_vars != coupledComponents("coupled_variables"))
    mooseError("Please match the number of orientation variants to coupled OPs (AuxNucleationProbability).");

  _coupled_OP.resize(_n_OP_vars);

  for(unsigned int i=0; i<_n_OP_vars; i++)
    _coupled_OP[i] = &coupledValue("coupled_variables", i);
}

Real
AuxNucleationProbability::computeValue()
{
  for(unsigned int i=0; i<_n_OP_vars; i++)
  {
    //Real OP = (*_coupled_OP[i])[_qp];
    if ( (*_coupled_OP[i])[_qp] > _OP_threshold)
     return -1;
  }

//  if (1.0 - exp(-1.0*_coupled_nuc_rate[_qp]*_dt) > _P_threshold)
//    return 0;

    return 1.0 - exp(-1.0*_coupled_nuc_rate[_qp]*_dt);
}

