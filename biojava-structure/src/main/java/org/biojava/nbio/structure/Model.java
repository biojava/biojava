/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 */
package org.biojava.nbio.structure;

import java.io.Serializable;
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
public class Model implements Serializable {
	
	private static final long serialVersionUID = 5320613424668781882L;

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
        	logger.info("No entity info could be found while adding chain with asym id {} (author id {}). Will consider it a polymer chain.", c.getId(), c.getName());
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
    
    @Override
    public String toString() {
    	return "["+polyChains.size()+" poly chains, "+nonPolyChains.size()+" non-poly chains, "+waterChains.size()+" water chains]";
    }
}
