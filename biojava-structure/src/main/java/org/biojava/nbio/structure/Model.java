package org.biojava.nbio.structure;

import java.util.ArrayList;
import java.util.List;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/** 
 * An internal utility class for StructureImpl to make it easier to manage poly and nonpoly chains.
 * Not to exposed to users through API.
 *
 * Created by andreas on 5/3/16.
 * @author Andreas Prlic
 * @since 5.0
 */
public class Model {
	
	private static final Logger logger = LoggerFactory.getLogger(Model.class);

    private List<Chain> polyChains;
    private List<Chain> nonPolyChains;
    private List<Chain> waterChains;

    public Model(){
        polyChains = new ArrayList<>();
        nonPolyChains = new ArrayList<>();
        waterChains = new ArrayList<>();
    }

    public List<Chain> getPolyChains() {
        return polyChains;
    }

    public List<Chain> getNonPolyChains() {
        return nonPolyChains;
    }
    
    public List<Chain> getWaterChains() {
    	return waterChains;
    }

    /**
     * Get all chains: polymeric, non-polymeric and water
     * @return
     */
    public List<Chain> getChains(){
        ArrayList<Chain> chains = new ArrayList<>();

        chains.addAll(polyChains);
        chains.addAll(nonPolyChains);
        chains.addAll(waterChains);

        chains.trimToSize();

        return chains;
    }

    public void setChains(List<Chain> modelChains) {

        polyChains.clear();
        nonPolyChains.clear();
        waterChains.clear();

        for (Chain c : modelChains){
            addChain(c);
        }
    }

    public void addChain(Chain c) {
        EntityInfo info = c.getEntityInfo();
        
        if ( info == null || info.getType() == null) {
        	logger.warn("No entity info could be found while adding chain with asym id {} (author id {}). Will consider it a polymer chain.", c.getId(), c.getName());
            polyChains.add(c);
            
        } else if ( info.getType() == EntityType.POLYMER) {
            polyChains.add(c);
            
        } else if (info.getType() == EntityType.NONPOLYMER) {
            nonPolyChains.add(c);
            
        } else if (info.getType() == EntityType.WATER) {
        	waterChains.add(c);
        	
        } else if (info.getType() == EntityType.MACROLIDE) {
        	logger.warn("Chain with asym id {} (author id {}) has entity type 'macrolide', considering it non-polymeric", c.getId(), c.getName());
        	nonPolyChains.add(c);
        	
        } else {
        	logger.warn("Chain with asym id {} (author id {}) has unsupported entity type '{}'. Will not add it to the Structure.", c.getId(), c.getName(), info.getType().toString());
        	// ignore it
        	
        }
    }

    /**
     * Returns the total number of chains in this model: polymeric, non-polymeric and water
     * @return
     */
    public int size() {
        return polyChains.size() + nonPolyChains.size() + waterChains.size();
    }
}
