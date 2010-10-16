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
 *    SymbolPropertyTableDB.java
 */
package org.biojava.bio.proteomics.aaindex;

import java.util.Set;

import org.biojava.bio.seq.db.IllegalIDException;
import org.biojava.bio.symbol.SymbolPropertyTable;

/**
 * Database of {@link org.biojava.bio.symbol.SymbolPropertyTable} objects.
 * @author <a href="mailto:Martin.Szugat@GMX.net">Martin Szugat</a>
 * @version $Revision$
 */
public interface SymbolPropertyTableDB {

    /**
     * Returns an iterator over 
     * {@link org.biojava.bio.symbol.SymbolPropertyTable} objects.
     * @return a new iterator
     */
     SymbolPropertyTableIterator tableIterator();
    
    /**
     * Returns the number of symbol property tables in the database. 
     * @return the number of tables
     */
     int numTables();
    
    /**
     * Returns the set of unique table names.
     * @return a set containing strings
     */
     Set names();
    
    /**
     * Returns the table with the specified name.
     * @param name the 
     * {@linkplain org.biojava.bio.symbol.SymbolPropertyTable#getName() name 
     * of the table}
     * @return the specified table
     * @throws IllegalIDException if no symbol property table with the specified name could be found.
     * @throws NullPointerException if <code>name</code> is 
     * <code>null</code>.
     */
     SymbolPropertyTable table(String name)
    throws IllegalIDException, NullPointerException;
}
