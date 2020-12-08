#include "BiasPlugin/BiasPlugin.hh"
#include "BiasEnergy.hh"

using namespace std::string_literals;

namespace prf
{
  inline namespace plugins
  {
    BiasPlugin::BiasPlugin() {}
    BiasPlugin::~BiasPlugin() {}

    int BiasPlugin::init()
    {
      return 0;
    }

    std::string BiasPlugin::name()
    {
      return "BiasPlugin";
    }

    bool BiasPlugin::has_update(std::string nm) const
    {
      return 0;
    }

    bool BiasPlugin::has_energy_term(std::string nm) const
    {
      return nm == "BiasEnergy";
    } 
    
    prf::Update * BiasPlugin::get_new_update(std::string nm)
    {
      return nullptr;
    }

    prf::Energy * BiasPlugin::get_new_energy_term(std::string nm)
    {
      prf::Energy * en=nullptr;
      if (nm=="BiasEnergy")
	{
	  en = new prf::BiasEnergy();
	  en->Name("BiasEnergy");
	}
      
      return en;
    }
    
    std::string BiasPlugin::brief()
    {
      return "A Go-like energy term with Cbeta-Cbeta attractions";
    }
    
    std::string BiasPlugin::help()
    {
      return R"(Strength in file go_parmeters and list of native contacts and distances in file contactmap.sidechain)";
    }
    
  }
}
