module;
#if PRF_USE_IMPORT_STD
#include <Utils/prf_fmt.hh>
import std;
import std.compat;
#else
// Standard library includes, if you want to test with a non-modular standard library.
#endif
module BiasPlugin2;
import prf;
using std::size_t;
namespace BiasPlugin2 {
// We need a parseCommand function to integrate with the
// ProFASi command parsing system. We just need to parse
// a single command. The loop over different commands in
// the settings file or the command line comes from other
// ProFASi routines.
auto BiasEnergy::parseCommand(prf::InstructionString gcmd) -> int
{
    if (gcmd.head() != name()) {
        return 0;
    }
    prf::Logger blog { 8 };
    gcmd = gcmd.tail();
    auto str = gcmd.str();
    auto opt = prf::utils_02::split_by(";")(str).to<std::vector<std::string>>();
    for (auto elem : opt) {
        elem = prf::utils::trim_str(elem);
        auto [k, v] = prf::utils_02::split(elem).to_pair();
        auto keyword = prf::utils::trim_str(k);
        auto value = prf::utils::trim_str(v);
        if ((keyword == "file" || keyword == "filename" || keyword == "add_contacts_from")
            && (not value.empty())) {
            blog << "Added new input contacts file " << value << "\n";
            infiles.push_back(std::string(value));
        } else if (keyword == "lambda_SC" && (not value.empty())) {
            lambda_SC = prf::svtod(value);
        } else if (keyword == "ksi_SC" && (not value.empty())) {
            ksi_SC = prf::svtod(value);
        }
    }
    return 1;
}

void BiasEnergy::init(const prf::Population& pop, [[maybe_unused]] prf::StateProperties& prp)
{
    if (initd)
        return;

    auto&& pc = pop.composition();
    contacts_data.clear();
    prf::Logger blog { 8 };

    size_t ifile {};
    blog << FMT("{}: {} input files for Go interaction parameters.\n", name(), infiles.size());
    for (auto&& inpfile : infiles) {
        blog << name() << ": Reading input file " << inpfile << "\n";
        try {
            if (not(prf::fs::exists(inpfile))) {
                throw prf::Exception { FMT("Input file '{}' must exist!\n", inpfile), __func__ };
            }
            auto ctxs = prf::FileAsLines::read(inpfile);
            // ctxs can now be used like a vector<string>, containing the lines of the file.
            for (auto&& line : ctxs) {
                std::istringstream ssin { line };
                std::string alabl1, alabl2;
                ssin >> alabl1 >> alabl2;
                prf::AtomURL u1, u2;
                u1.str(alabl1); // If either string is not a proper AtomURL, this str()
                u2.str(alabl2); // function will throw an exception, and we will quit
                                // the init function.
                auto opt_a1 = pc->atom(u1); // We ask the population composition to
                auto opt_a2 = pc->atom(u2); // give us the atoms for those URLs
                if (opt_a1.has_value() and opt_a2.has_value()) {
                    // opt_a1 has a value if the first URL designated an existing atom in the population
                    auto aid1 = opt_a1->globalIndex();
                    auto aid2 = opt_a2->globalIndex();
                    if (aid1 > aid2)
                        std::swap(aid1, aid2);
                    double tmp {};
                    std::optional<double> atdist, gaussian_width, gaussian_weight;
                    if (ssin >> tmp) {
                        atdist = tmp;
                        if (ssin >> tmp) {
                            gaussian_width = tmp;
                            if (ssin >> tmp) {
                                gaussian_weight = tmp;
                            }
                        }
                    }
                    if (atdist) {
                        auto loc = sr::find_if(contacts_data,
                            [aid1, aid2](auto&& ct) {
                                return ct.atom1 == aid1 && ct.atom2 == aid2;
                            });
                        SingleContactInteractionParameters pars { *atdist, gaussian_width.value_or(ksi_SC),
                            gaussian_weight.value_or(lambda_SC), ifile };
                        if (loc == contacts_data.end()) {
                            auto&& newelem = contacts_data.emplace_back(SingleContact(aid1, aid2));
                            newelem.par_sets.push_back(pars);
                        } else {
                            loc->par_sets.push_back(pars);
                        }
                    }
                }
            }
            ++ifile;
        } catch (prf::Exception& err) {
            err(FMT("Syntax error while processing contacts input file '{}'", infiles[ifile]), __func__);
            throw;
        }
    }
    blog << "Contacts list has " << contacts_data.size() << " entries.\n";
    for (size_t ic = 0UL; ic < contacts_data.size(); ++ic) {
        auto&& C = contacts_data[ic];
        blog << FMT("{}\t{}\t{}\t", ic, pc->url(C.atom1).str(), pc->url(C.atom2).str());
        for (auto&& par : C.par_sets)
            blog << infiles[par.label] << ", ";
        blog << "\n";
    }

    myobsindex = prp.next_obs(name());
    contribprof = prp.next_profile(FMT("{}_contributions", name()));
    sourceprof = prp.next_profile(FMT("{}_sources", name()));
    prp[contribprof].resize(contacts_data.size());
    prp[sourceprof].resize(contacts_data.size());
}
auto contact_energy(double rnow, double rref, double wdt, double wgt) -> double
{
    auto d = rnow - rref;
    return -wgt * std::exp(-d * d / (2. * wdt * wdt));
}

auto BiasEnergy::operator()(const prf::Population& pop, prf::StateProperties& prp) const -> double
{
    double vval {};
    auto&& ps = pop.state();
    auto&& xyz = ps.xyz();

    auto&& contribs = prp[contribprof];
    auto&& sources = prp[sourceprof];

    for (size_t ic = 0UL; ic < contacts_data.size(); ++ic) {
        auto&& C = contacts_data[ic];
        double r = prf::AC::dist(C.atom1, C.atom2, xyz);
        auto [en, lb] = sr::min(
            C.par_sets | sv::transform([r](auto&& gaussian) -> std::pair<double, size_t> {
                return {
                    contact_energy(r, gaussian.mean, gaussian.width, gaussian.weight),
                    gaussian.label
                };
            }),
            [](auto g1, auto g2) { return g1.first < g2.first; });
        vval += en;
        contribs[ic] = en;
        sources[ic] = lb;
    }

    prp[myobsindex] = vval;
    return vval;
}
}
