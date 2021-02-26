package org.biojava.nbio.structure.chem;

import java.util.HashMap;
import java.util.Map;

/**
 * A representation of the Chemical Component Dictionary.
 *
 * @author Andreas Prlic
 * @since 1.7
 * @see <a href="http://mmcif.rcsb.org/dictionaries/">link into mmCIF dictionary</a>
 *
 */
public class ChemicalComponentDictionary {
    private final Map<String, ChemComp> dictionary;
    private final Map<String,String> replaces;
    private final Map<String,String> isReplacedBy;

    public ChemicalComponentDictionary(){
        dictionary = new HashMap<>();
        replaces  = new HashMap<>();
        isReplacedBy = new HashMap<>();
    }

    public boolean isReplaced(ChemComp c) {
        return isReplaced(c.getId());
    }

    public boolean isReplaced(String id) {
        return isReplacedBy.containsKey(id);
    }

    public boolean isReplacer(ChemComp c) {
        return isReplacer(c.getId());
    }

    public boolean isReplacer(String id) {
        return replaces.containsKey(id);
    }

    /**
     * If ChemComp is replaced by another one, get the newer version otherwise return the same ChemComp again.
     * @param c
     * @return get the component that replaced ChemComp.
     */
    public ChemComp getReplacer(ChemComp c) {
        return getReplacer(c.getId());
    }

    public ChemComp getReplacer(String id) {
        if (isReplaced(id)) {
            return dictionary.get(isReplacedBy.get(id));
        }
        return dictionary.get(id);
    }

    /**
     * If ChemComp is replacing another one, get the old version otherwise return the same ChemComp again.
     * @param  c the ChemComp for which older versions should be looked up.
     */
    public ChemComp getReplaced(ChemComp c) {
        return getReplaced(c.getId());
    }

    public ChemComp getReplaced(String id){
        if (isReplacer(id)) {
            return dictionary.get(replaces.get(id));
        }
        return dictionary.get(id);
    }

    /**
     * Get the parent of a component. If component has no parent, return null
     * @param c
     * @return get the parent component or null if ChemComp has no parent.
     */
    public ChemComp getParent(ChemComp c) {
        if (c.hasParent()) {
            return dictionary.get(c.getMonNstdParentCompId());
        }
        return null;
    }

    /**
     * Add a new component to the dictionary
     * @param comp
     */
    public void addChemComp(ChemComp comp) {
        dictionary.put(comp.getId(),comp);
        String rep = comp.getPdbxReplaces();
        if (rep != null && !rep.equals("?")) {
            replaces.put(comp.getId(),rep);
        }

        String isrep = comp.getPdbxReplacedBy();
        if (isrep != null && !isrep.equals("?")) {
            isReplacedBy.put(comp.getId(), isrep);
        }
    }

    /**
     * Returns the number of ChemComps in this dictionary
     * @return nr. of ChemComps
     */
    public int size() {
        return dictionary.size();
    }

    public ChemComp getChemComp(String id) {
        return dictionary.get(id);
    }
}
