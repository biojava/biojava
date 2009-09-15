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

import java.util.List;

import org.biojava.bio.symbol.Location;

/**
 * <p>UnequalLengthAlignment has the following behavior.  Two or more
 * SymbolLists may align in such a way that their ends do not
 * overlap.</p>
 *
 * <pre>
 *      example
 *         1         aaaaaatttcttt
 *         2               tttgtttggggggc
 * </pre>
 *
 * <p>
 *     length returns ??                                        <br>
 *     symbolAt(1,1) returns 20                                 <br>
 *     symbolAt(2,1) returns null -- NOT an exception           <br>
 *     symbolAt(2,99) throws NoSuchElementException             <br>
 *     leftMost returns 1                                       <br>
 *     rightMost returns 2                                      <br>
 *     locInAlignment (1) returns (1,13)                        <br>
 *     locInAlignment (2) returns (7,20)                        <br>
 *     alignmentRange() returns (7,13)                          <br>
 * </p>
 *
 * @author David Waring
 */

public interface UnequalLengthAlignment extends Alignment{
        
       
        /**
        * The location of an individual SymbolList relative to overall Alignment
        */        
        public Location locInAlignment(Object label);
                
        /**
        * Returns a list labels, of all seqs that cover that column
        */        
        public List<Object> labelsAt(int column);

        /**
        * Returns list of all the labels that intersect that range
        */      
        public List<Object> labelsInRange(Location loc);   
       


}
