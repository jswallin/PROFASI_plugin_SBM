
module;
#if PRF_USE_IMPORT_STD
import std;
#else
#include <any>
#include <optional>
#include <string>
#endif
export module BiasPlugin2;
import prf;
export import :BiasEnergy;
export namespace BiasPlugin2{
// The auto-generated class is left empty as a tag class
// You can fill it with any code you need. ProFASi does not depend on
// this class having any properties other than existing. Feel free to
// adapt it for the intended functionality.
class BiasPlugin2 {
};

// The 4 following functions must be present to satisfy requirements of a valid plugin
// under ProFASi plugin api version 1. The name() and api() functions are probably
// just fine as automatically generated. You can leave them alone.
inline auto name(const BiasPlugin2&) -> std::string { return "BiasPlugin2"; }
inline auto api(const BiasPlugin2&) -> int { return 1; }

// provides() must return an XML formatted string detailing what the plugin does. See source file!
auto provides(const BiasPlugin2&) -> std::string;

// get() responds to resource requests. Returns answers wrapped in optional<prf::Any>. See source file!
auto get(BiasPlugin2&, const std::string & resourcename) -> std::optional<prf::Any>;

}
