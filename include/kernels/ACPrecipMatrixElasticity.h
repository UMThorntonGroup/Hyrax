/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  14 February 2014
*
*************************************************************************/

#ifndef ACPRECIPMATRIXELASTICITY_H
#define ACPRECIPMATRIXELASTICITY_H

#include "ACBulk.h"
#include "ElasticityTensorR4.h"
#include "RankTwoTensor.h"

class ACPrecipMatrixElasticity;

template<>
InputParameters validParams<ACPrecipMatrixElasticity>();

class ACPrecipMatrixElasticity : public ACBulk
{
public:

  ACPrecipMatrixElasticity(const std::string & name, InputParameters parameters);

protected:

  /**
   * computeDFDOP()
   * @return returns the partial(elastic free energy/order parameter)
   */
  virtual Real computeDFDOP(PFFunctionType type);

  // system elasticity tensor, varies in space
  MaterialProperty<ElasticityTensorR4> & _elasticity_tensor;
  MaterialProperty<std::vector<ElasticityTensorR4> > & _dn_elasticity_tensor;

  MaterialProperty<RankTwoTensor> & _elastic_strain;
  MaterialProperty<std::vector<RankTwoTensor> > & _dn_misfit_strain;

  // orientation variant number for this kernel (1 to n)
  unsigned int _OP_number;

  Real _scaling_factor;

private:

};

#endif //ACPRECIPMATRIXELASTICITY_H
