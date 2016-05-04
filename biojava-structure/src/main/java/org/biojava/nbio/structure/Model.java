package org.biojava.nbio.structure;

import java.util.ArrayList;
import java.util.List;

/** An internal utility class for StructureImpl to make it easier to manage poly and nonpoly chains.
 * Not to exposed to users through API.
 *
 *
 * Created by andreas on 5/3/16.
 */
public class Model {

    List<Chain> polyChains;
    List<Chain> nonPolyChains;

    public Model(){
        polyChains = new ArrayList<>(0);
        nonPolyChains = new ArrayList<>(0);
    }

    public List<Chain> getPolyChains() {
        return polyChains;
    }

    public List<Chain> getNonPolyChains() {
        return nonPolyChains;
    }

    public List<Chain> getChains(){
        ArrayList<Chain> chains = new ArrayList<>();

        chains.addAll(polyChains);
        chains.addAll(nonPolyChains);

        chains.trimToSize();

        return chains;
    }

    public void setChains(List<Chain> modelChains) {

        polyChains.clear();
        nonPolyChains.clear();

        for (Chain c : modelChains){
            addChain(c);
        }
    }

    public void addChain(Chain c) {
        EntityInfo info = c.getEntityInfo();
        if ( info == null || info.getType() == null) {

            polyChains.add(c);
        } else if ( info.getType().equals(EntityType.POLYMER)){
            polyChains.add(c);
        } else {
            nonPolyChains.add(c);
        }
    }

    public int size() {
        return polyChains.size() + nonPolyChains.size();
    }
}
