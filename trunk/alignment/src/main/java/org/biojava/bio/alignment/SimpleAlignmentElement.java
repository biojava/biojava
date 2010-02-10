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

package org.biojava.bio.alignment;

import org.biojava.bio.BioError;
import org.biojava.bio.BioException;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.SymbolList;


/** 
 * <p>SimpleSimpleAlignment is a simple implementation of
 * AlignmentElement.</p>
 *
 * @author David Waring
 */

public class SimpleAlignmentElement implements AlignmentElement{
    
    protected Object label;
    protected Location loc;
    protected SymbolList seq;
    
    public SimpleAlignmentElement(Object label,SymbolList seq,Location loc)throws BioException{
        if ((loc.getMax() - loc.getMin() + 1) != seq.length()){
            throw new BioException("Sequence length and location length do not match");
        }
        this.label = label;
        this.seq = seq;
        this.loc = loc;
    }
    
    public Object getLabel(){
        return label;
    }
    public Location getLoc(){
        return loc;
    }
    public SymbolList getSymbolList(){
        return seq;
    }
    public void setLoc(Location nLoc)throws BioError{
         if ((nLoc.getMax() - nLoc.getMin() + 1) != seq.length()){
            throw new BioError("Sequence length and location length do not match");
        }
        this.loc = nLoc;
    }

}
