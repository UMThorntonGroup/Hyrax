/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  24 July 2013
*
*************************************************************************/

#ifndef NUCLEATIONPOSTPROCESSORACTION_H
#define NUCLEATIONPOSTPROCESSORACTION_H

#include "InputParameters.h"
#include "Action.h"

//forward declaration
class NucleationPostprocessorAction;

template<>
InputParameters validParams<NucleationPostprocessorAction>();

class NucleationPostprocessorAction : public Action
{
public:

  NucleationPostprocessorAction(InputParameters params);

  virtual void act();

protected:

private:

  int _num_OPs;
  std::string _OP_name_base;
  FileName _particle_volume_name_base;
  FileName _Avrami_name_base;

  PostprocessorName _mesh_volume;
  Real _equil_fraction;

  const Real _threshold;

  std::string _nucleation_userobject;
};

#endif //NUCLEATIONPOSTPROCESSORACTION_H
