#include <PluginSupport/PluginBase.hh>

namespace prf
{
  inline namespace plugins
  {
    class BiasPlugin : public PluginBase {
    public:
      BiasPlugin();
      ~BiasPlugin();
      int init() override;
      std::string name() override;
      std::string brief() override;
      std::string help() override;
      bool has_update(std::string nm) const override;
      bool has_energy_term(std::string nm) const override;
      prf::Update* get_new_update(std::string nm) override;
      prf::Energy* get_new_energy_term(std::string nm) override;
    };
  }
}
