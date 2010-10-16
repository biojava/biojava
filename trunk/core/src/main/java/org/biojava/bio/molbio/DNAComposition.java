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

package org.biojava.bio.molbio;

import org.biojava.bio.BioError;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.IllegalSymbolException;

/**
 * Computes composition statistics about a DNA <code>SymbolList</code>.
 * 
 * @author Mark Schreiber
 * @since 1.6
 */
public class DNAComposition extends Composition{

    
    public DNAComposition(){
        super();
        try {
            setSymbolList(DNATools.createDNA("acgt"));
        } catch (IllegalSymbolException ex) {
            System.err.println("Severe Error. Cannot create DNA SymbolList with 'acgt");
        }
    }
    
    /**
     * Get the relative compositon of 'A'.
     * @return between 0.0 and 1.0
     */
    public double getA(){
        try {
            return getDistribution().getWeight(DNATools.a());
        } catch (IllegalSymbolException ex) {
            throw new BioError("Severe Error with DNA configuration.",ex);
        }
    }
    
    /**
     * Get the relative compositon of 'C'.
     * @return between 0.0 and 1.0
     */
    public double getC(){
        try {
            return getDistribution().getWeight(DNATools.c());
        } catch (IllegalSymbolException ex) {
            throw new BioError("Severe Error with DNA configuration.",ex);
        }
    }
    
    /**
     * Get the relative compositon of 'G'.
     * @return between 0.0 and 1.0
     */
    public double getG(){
        try {
            return getDistribution().getWeight(DNATools.g());
        } catch (IllegalSymbolException ex) {
            throw new BioError("Severe Error with DNA configuration.",ex);
        }
    }
    
    /**
     * Get the relative compositon of 'T'.
     * @return between 0.0 and 1.0
     */
    public double getT(){
        try {
            return getDistribution().getWeight(DNATools.t());
        } catch (IllegalSymbolException ex) {
            throw new BioError("Severe Error with DNA configuration.",ex);
        }
    }
}
