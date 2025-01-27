#include <rdkit/GraphMol/GraphMol.h>
#include <rdkit/GraphMol/RDKitBase.h>
#include "./reaction_engine.hpp"
#include <string>
#include <vector>
#include <stdexcept> 


ReactionEngine::ReactionEngine(int rxn_size, 
                                int num_reactants_per_rxn,
                                int first_reactant_count, 
                                bool use_names=true, 
                                bool use_sample_ids=true) {
    _rxn_size = rxn_size;
    _num_reactants_per_rxn = num_reactants_per_rxn;
    _first_reactant_count = first_reactant_count;
    _use_names = use_names;
    _use_sample_ids = use_sample_ids;

}

ReactionEngine::~ReactionEngine()

void ReactionEngine::set_reactions(const std::vector<std::string> &smarts) {
    _rxn.clear();
    int size = smarts.size();
    for (int i = 0; i < size; ++i) {
        _rxn.push_back(RDKit::ReactionFromSmarts(smarts[i]));
    }
}

void ReactionEngine::set_names(const std::vector<std::string> &names) {
    if (_use_names) {
        _names = names;
    }
    else {
        throw std::runtime_error("Use names is false.");
    }
}
void ReactionEngine::set_mols(const std::vector<std::string> &smiles) {
    _mols.clear();
    int size = smiles.size();
    for (int i = 0; i < size; ++i) {
        _mols.push_back(RDKit::ROMol(smiles[i]));
    }
}

void run_multi_double_reaction(std::vector<std::string> &products) {
}

