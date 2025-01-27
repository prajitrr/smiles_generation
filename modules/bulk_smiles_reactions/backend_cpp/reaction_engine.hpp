#ifndef REACTION_ENGINE_HPP
#define REACTION_ENGINE_HPP

#include <rdkit/GraphMol/ChemReactions.h>
#include <rdkit/GraphMol/RDKitBase.h>
#include <vector>
#include <string>

class ReactionEngine {
public:
    ReactionEngine(int rxn_size, 
                   int num_reactants,
                   int first_reactant_count, 
                   bool use_names=true, 
                   bool use_sample_ids=true);
    ~ReactionEngine();    
    void set_reactions(const std::vector<std::string> &smarts);
    void set_names(const std::vector<std::string> &names);
    void set_mols(const std::vector<std::string> &smiles);
    void set_sample_ids(const std::vector<std::string> &sample_ids);
    void set_first_reactant_count(int count);
    void run_multi_double_reaction(std::vector<std::string> &products);
    std::string run_reaction(const std::string &smiles);

private:
    std::vector<RDKit::ChemicalReaction*> _rxn;
    std::vector<RDKit::ROMol*> _mols;
    std::vector<std::string> _names;
    std::vector<std::string> _sample_ids;
    int _first_reactant_count;
    int _rxn_size;
    int _num_reactants_per_rxn;
    bool _use_names;
    bool _use_sample_ids;
};

#endif // REACTION_ENGINE_HPP