module;
#if PRF_USE_IMPORT_STD
import std;
import std.compat;
#else
// Standard library includes, if you want to test with a non-modular standard library.
#endif
export module BiasPlugin2:BiasEnergy;
import prf;
import prf.io;
using std::size_t;
namespace BiasPlugin2 {
// Parameters width, weight and nominal separation (mean of the Gaussian)
// for a single pair of entities in contact. This is a helper class, not
// meant for outside use. Therefore, it is not exported.
struct SingleContactInteractionParameters {
    double mean{0.}, width { 1.0 }, weight { 1 };
    size_t label {};
};
// Single Contact representation
struct SingleContact {
    prf::URL::AtomIndex atom1 {}, atom2 {};
    std::vector<SingleContactInteractionParameters> par_sets {};
};
export class BiasEnergy {
public:
    using population_type = prf::Population;
    using properties_type = prf::StateProperties;
    using transition_type = prf::ConformationChange;
    using result_type = double;

    BiasEnergy() = default;
    ~BiasEnergy() = default;

    auto name() const noexcept -> std::string { return nm; }
    void name(std::string_view gnm) { nm = gnm; }

    // An index used inside ProFASi to identify energy terms and other observables.
    // The value is queried from a StateProperties object in the function init, and remembered
    auto index() const noexcept -> prf::StateProperties::ObsIndexType { return myobsindex; }

    // The following function is called by ProFASi during initialisation. The Population
    // will be ready with all atoms, degrees of freedom etc. The StateProperties object
    // can be used to make room for auxilliary variables for the term. Anything that remains
    // constant for a population composition should be stored in the class and initialised
    // properly in init(). Anything calculated that depends on atom coordinates or degrees
    // of freedom should be placed in a private storage obtained from the StateProperties
    // object. The idea is that in a parallel (multi-threaded) simulation, many copies of
    // the population state, i.e., xyz and torsion coordinates, could be simultaneously
    // in play. Storing coordinate dependent changable quantities in the energy terms
    // creates a big problem of managing state dependent properties and keeping them in
    // sync with their properties. This is managed for you in ProFASi.
    void init(const prf::Population& pop, prf::StateProperties& prp);

    // This is the function where you should calculate the energy term. The input signature
    // needs to be exactly this, since this is how ProFASi calls it. The Population argument
    // has the coordinates you need to calculate. The result should be stored in the
    // StateProperties and returned.
    auto operator()(const prf::Population& pop, prf::StateProperties& prp) const -> double;

    // For faster MC updates, we use delta calculations in ProFASi. Only changes to
    // energy values are calculated during MC moves. When a move is proposed, ProFASi
    // pre-calculates a "ConformationChange" object, which summarizes many important
    // aspects of the move. The conformation change object gives access to the population
    // state before and after the Monte Carlo move, of course. But in addition, it
    // can be queried for things like which atoms were moved, which degrees of freedom
    // were changed etc. Since two states are involved, there are two StateProperties
    // arguments, with obvious meanings. Uncomment this (and the next two functions)
    // if you decide to implement the fast update method in the source file.
    // auto operator()(const ConformationChange& T,
    //    StateProperties& props1, StateProperties& props2) const -> double;

    // If you implement the operator() function for a ConformationChange, you will
    // probably have to store partial calculations in state properties 1 and 2. Whether
    // those calculations are to be kept or discarded depends on whether the move
    // was accepted or rejected. The following two functions are for proper book keeping.
    // void apply_move(const ConformationChange& T,
    //    StateProperties& sp1, StateProperties& sp2) const;
    // void undo_move(const ConformationChange& T,
    //    StateProperties& sp1, StateProperties& sp2) const;

    // Do you want the new term to be available for interactive evaluation inside the
    // interactive shell in StAn? If yes, uncomment and implement the following function.
    // void exec(const prf::InstructionString& s);
    // Do you want the new term to be able to auto-complete commands directed to it
    // inside an interactive shell? If yes, you should implement the following
    // function.
    // void helpComplete(std::span<std::string> parts, std::string_view cur);
    auto parseCommand(prf::InstructionString s) -> int;

private:
    std::string nm { "BiasEnergy" };
    prf::StateProperties::ObsIndexType myobsindex;
    prf::StateProperties::ProfileIndexType contribprof, sourceprof;
    std::vector<SingleContact> contacts_data;
    std::vector<std::string> infiles;
    double lambda_SC { 1. }, ksi_SC { 1. };
    bool initd { false };
};
}
template <>
inline constexpr auto prf::typestring<BiasPlugin2::BiasEnergy>
    = "BiasPlugin2::BiasEnergy";
