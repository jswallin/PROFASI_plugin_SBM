module;
#if PRF_USE_IMPORT_STD
import std;
#else
#include <any>
#include <optional>
#include <string>
#endif
module BiasPlugin2;
import prf;
namespace BiasPlugin2 {

auto provides(const BiasPlugin2&) -> std::string
{
    return R"xmlout(<BiasPlugin2>
<brief>ProFASi plugin SBM for ProFASi 2.4</brief>
<help_text>This plugin returns a structure based energy term where
two sets of contacts are used to provide competing interactions
between different parts of the chain.


</help_text>
<resources>
<resource>
<request>SBM</request>
<return_type>BiasEnergy</return_type>
<wrapper>prf::FFProxy</wrapper>
<api>1</api>
<brief>Structure based term for 2 reference structures</brief>
</resource>
</resources>
<profasi_deps>profasi_combined
</profasi_deps>
<generated_classes>BiasEnergy	</generated_classes>
</BiasPlugin2>
)xmlout";
}

auto get(BiasPlugin2&, const std::string& resourcename) -> std::optional<prf::Any>
{
    std::optional<prf::Any> ans;
    if (resourcename == "BiasEnergy") {
        ans = prf::Any { prf::FFProxy { BiasEnergy {} /* Replace with a properly initialized object */ } };
    }
    return ans;
}
}
