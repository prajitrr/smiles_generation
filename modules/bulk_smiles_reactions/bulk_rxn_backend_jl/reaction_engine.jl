
import RDKitMinimalLib



function get_reaction_smiles(reaction_smiles::String)
    # Parse the reaction SMILES string
    thing = RDKitMinimalLib.get_mol(reaction_smiles)
    old_thing = RDKitMinimalLib.get_smiles(thing)
    return old_thing
end

function multi_double_reaction(smarts::Vector{AbstractString}, 
                               smiles::Vector{AbstractString},
                               num_first_reactant::Int,
                               multi_react_first::Bool,
                               multi_react_second::Bool,
                               names::Vector{AbstractString} = [],
                               sample_ids::Vector{AbstractString} = []
                               )
    
    use_names::Bool = length(names) > 0
    use_ids::Bool = length(sample_ids) > 0

    reactions::Vector
     
    # Parse the reaction SMILES string
    thin = RDKitMinimalLib.get_mol(reaction_smiles)
    old_thing = RDKitMinimalLib.get_smiles(thing)
    return old_thing
end

