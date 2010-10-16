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
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.SymbolList;

/** 
 * <p>AlignmentElement is a class which represents a SymbolList and its
 * location within an Alignment This is for use in
 * UnequalLengthAlignments and ARAlignments.</p>
 *
 * @author David Waring
 */

public interface AlignmentElement{
    
    public String getLabel();
    public Location getLoc();
    public SymbolList getSymbolList();
    public void setLoc(Location nLoc)throws BioError;

}
