/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  24 July 2013
*
*************************************************************************/

#ifndef NUCLEATIONAUXACTION_H
#define NUCLEATIONAUXACTION_H

#include "InputParameters.h"
#include "Action.h"

//forward declaration
class NucleationAuxAction;

template<>
InputParameters validParams<NucleationAuxAction>();

class NucleationAuxAction : public Action
{
public:

  NucleationAuxAction(InputParameters params);

  virtual void act();

protected:

private:
  //for Action
  unsigned int _num_OPs;
  std::string _OP_name_base;
  std::string _bulk_energy_name_base;
  std::string _nucleation_rate_name_base;
  std::string _nucleation_probability_name_base;
  std::string _deltaGstar_name_base;

  //for AuxChem or AuxChemElastic
  // bool _use_auxchem;
  std::string _bulk_energy_name;

  std::string _coupled_conserved_var;
  Real _precip_conserved;
  Real _precip_nonconserved;

  //for NucleationRate
  Real _Z;
  Real _beta_star;
  Real _linear_density;
  Real _gamma;                 //also for DeltaGStar
  Real _Kb;
  Real _temperature;
  Real _scale_factor;

  //for NucleationProbability
  Real _OP_threshold;

};

#endif //NUCLEATIONAUXACTION_H
