/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  25 July 2013
*
*************************************************************************/

#include "NucleationPostprocessorAction.h"
#include "Factory.h"
#include "FEProblem.h"

template<>
InputParameters validParams<NucleationPostprocessorAction>()
{
  InputParameters params = validParams<Action>();

  params.addRequiredParam<std::string>("OP_name_base", "name base of the OP variable");
  params.addRequiredParam<int>("num_OPs", "# of order parameter variables to process on");
  params.addRequiredParam<FileName>("particle_volume_name_base", "name base of particle volume file");
  params.addParam<FileName>("Avrami_name_base", "name base of Avrami info file");
  params.addRequiredParam<PostprocessorName>("mesh_volume", "postprocessor from which to get mesh volume");
  params.addParam<Real>("equil_fraction", -1.0, "Equilibrium volume fraction of 2nd phase for Avrami analysis");
  params.addRequiredParam<Real>("threshold", "The threshold value for which a new bubble may be started");

  params.addRequiredParam<std::string>("nucleation_userobject", "The name of the userobject for nucleation event location");
  return params;
}

NucleationPostprocessorAction::NucleationPostprocessorAction(InputParameters params) :
    Action(params),
    _num_OPs(getParam<int>("num_OPs")),
    _OP_name_base(getParam<std::string>("OP_name_base")),
    _particle_volume_name_base(getParam<FileName>("particle_volume_name_base")),
    _mesh_volume(getParam<PostprocessorName>("mesh_volume")),
    _equil_fraction(getParam<Real>("equil_fraction")),
    _threshold(getParam<Real>("threshold")),
    _nucleation_userobject(getParam<std::string>("nucleation_userobject"))
{
}

void
NucleationPostprocessorAction::act()
{
  for (unsigned int i = 1; i <= _num_OPs; ++i)
  {
    //create names
    std::string variable_name = _OP_name_base;
    std::string particle_volume_name = _particle_volume_name_base;

    std::stringstream out;
    out << i;

    variable_name.append(out.str());
    particle_volume_name.append(variable_name);
    particle_volume_name.append(".csv");

    // get and set input parameters for NodalVolumeFraction
    InputParameters action_params = _factory.getValidParams("NodalVolumeFraction");

    action_params.set<std::vector<VariableName> >("variable")
      = std::vector<VariableName> (1, variable_name);
    action_params.set<FileName>("bubble_volume_file") = particle_volume_name;
    action_params.set<Real>("threshold") = _threshold;
    action_params.set<PostprocessorName>("mesh_volume") = _mesh_volume;

    //check to see if Avrami file ouput is required
    if(_pars.isParamValid("Avrami_name_base"))
    {
      std::string Avrami_name(getParam<FileName>("Avrami_name_base"));
      Avrami_name.append(variable_name);
      Avrami_name.append(".csv");

      action_params.set<FileName>("Avrami_file") = Avrami_name;
      action_params.set<Real>("equil_fraction") = _equil_fraction;
    }

    //make postprocessor name
    std::string postprocessor_name = "NodalVolumeFraction_";
    postprocessor_name.append(variable_name);

    //action_params.print();
    _problem->addPostprocessor("NodalVolumeFraction", postprocessor_name, action_params);

    //Get and set input parameters for NucleiInformation
    action_params = _factory.getValidParams("NucleiInformation");

    action_params.print();

    action_params.set<UserObjectName>("nucleation_userobject") = _nucleation_userobject;


    action_params.set<int>("OP_number") = i;

    postprocessor_name = "NucleiInformation_";
    postprocessor_name.append(variable_name);

    action_params.print();

    _problem->addPostprocessor("NucleiInformation", postprocessor_name, action_params);
  }
}
