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

/*
 *    SymbolPropertyTableIterator.java
 */
package org.biojava.bio.proteomics.aaindex;

import java.util.NoSuchElementException;

import org.biojava.bio.BioException;
import org.biojava.bio.symbol.SymbolPropertyTable;

/**
 * Iterator over {@link org.biojava.bio.symbol.SymbolPropertyTable} objects.
 * @author <a href="mailto:Martin.Szugat@GMX.net">Martin Szugat</a>
 * @version $Revision$
 */
public interface SymbolPropertyTableIterator {
    /**
     * Checks if there is a further 
     * {@link org.biojava.bio.symbol.SymbolPropertyTable} object.
     * @return <code>true</code> if a call to the {@link #nextTable()} method
     * is valid, <code>false</code> otherwise.
     */
    boolean hasNext();
    
    /**
     * Returns the next {@link org.biojava.bio.symbol.SymbolPropertyTable} 
     * object.
     * @return a symbol property table
     * @throws NoSuchElementException if there is no further symbol property 
     * table.
     * @throws BioException if the next symbol property table could not be
     * retrieved.
     */
    SymbolPropertyTable nextTable() 
    throws NoSuchElementException, BioException;
}
