/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  24 October 2012
*
*************************************************************************/

#ifndef ACBULKSIMMONS_H
#define ACBULKSIMMONS_H

#include "ACBulkCoupled.h"

//Forward Declarations
class ACBulkSimmons;

template<>
InputParameters validParams<ACBulkSimmons>();

/**
 * ACBulkSimmons replicates dF/deta for the free energy used in Simmons (2001) concurrent nucleation
 * and growth paper. It uses the PFMobilityLandau materials class.
 */

class ACBulkSimmons : public ACBulkCoupled
{
public:

  ACBulkSimmons(const InputParameters & parameters);

protected:

  /**
   * computeDFDOP()
   * @return returns the partial(bulk free energy/order parameter)
   */
  virtual Real computeDFDOP(PFFunctionType type);

private:

};

#endif //ACBULKSIMMONS_H

