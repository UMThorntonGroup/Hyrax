/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  30 August 2013
*
*************************************************************************/

#ifndef ACFERROELECTRIC_H
#define ACFERROELECTRIC_H

#include "ACBulk.h"

//forward declarations
class ACFerroelectric;

template<>
InputParameters validParams<ACFerroelectric>();

class ACFerroelectric : public ACBulk
{
public:
  ACFerroelectric(const InputParameters & parameters);

protected:
  virtual Real computeDFDOP(PFFunctionType type);

  const MaterialProperty<Real> & _a1;
  const MaterialProperty<Real> & _a11;
  const MaterialProperty<Real> & _a12;

  unsigned int _n_OP_vars;
  unsigned int _OP_number;

  std::vector<VariableValue *> _coupled_OP_vars;

private:

};


#endif //ACFERROELECTRIC_H
