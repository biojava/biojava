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
 *    SimpleSymbolPropertyTableDB.java
 */
package org.biojava.bio.proteomics.aaindex;

import java.util.Hashtable;
import java.util.Iterator;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.Set;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.db.IllegalIDException;
import org.biojava.bio.symbol.SymbolPropertyTable;

/**
 * A simple implementation of a symbol property table database.
 * @author <a href="mailto:Martin.Szugat@GMX.net">Martin Szugat</a>
 * @version $Revision$
 */
public class SimpleSymbolPropertyTableDB implements SymbolPropertyTableDB {
    
    /* PRIVATE CLASSES */
    
    /**
     * Iterator over symbol property tables.
     * @author <a href="mailto:Martin.Szugat@GMX.net">Martin Szugat</a>
     * @version $Revision$
     */
    private static final class SimpleSymbolPropertyTableIterator implements
            SymbolPropertyTableIterator {
        
        /* PRIVATE FIELDS */

        /**
         * The internal iterator.
         */
        private Iterator iterator = null;
        
        /* PUBLIC CONSTRUCTORS */
        
        /**
         * Initializes the iterator.
         * @param iterator the internal iterator
         * @throws NullPointerException if <code>iterator</code> is
         * <code>null</code>.
         */
        public SimpleSymbolPropertyTableIterator(
                Iterator iterator) 
        throws NullPointerException {
            super();
            if (iterator == null) {
                throw new NullPointerException("iterator is null.");
            }
            this.iterator = iterator;
        }
        
        /* INTERFACE SymbolPropertyTableIterator */

        /**
         * {@inheritDoc}
         */
        public boolean hasNext() {
            return iterator.hasNext();
        }

        /**
         * {@inheritDoc}
         */
        public SymbolPropertyTable nextTable() throws NoSuchElementException, 
        BioException {
            return (SymbolPropertyTable) iterator.next();
        }

    }

    /* PRIVATE FIELDS */
    
    /**
     * Internal map of symbol property tables. 
     */
    private Map map = null;
    
    /* PUBLIC CONSTRUCTORS */

    /**
     * Initializes the database. 
     */
    public SimpleSymbolPropertyTableDB() {
        super();
        map = new Hashtable();
    }
    
    /**
     * Initializes the database by copying all symbol property tables from 
     * a given iterator into the database.
     * @param tableIterator an iterator over symbol property tables.
     * @throws BioException if the symbol property tables could not be
     * iterated.
     */
    public SimpleSymbolPropertyTableDB(
            SymbolPropertyTableIterator tableIterator) 
    throws BioException {
        this();
        if (tableIterator == null) {
            throw new NullPointerException("tableIterator is null.");
        }
        while (tableIterator.hasNext()) {
            addTable(tableIterator.nextTable());
        }
    }
    
    /* PUBLIC METHODS */
    
    /**
     * Adds a symbol property table to the database. Overrides an existing 
     * table entry with the same name.
     * @param table the symbol property table to add.
     * @throws NullPointerException if <code>table</code> is <code>null</code>.
     */
    public void addTable( SymbolPropertyTable table) 
    throws NullPointerException {
        if (table == null) {
            throw new NullPointerException("table is null.");
        }
        map.put(table.getName(), table);
    }
    
    /* INTERFACE SymbolPropertyTableDB */

    /**
     * {@inheritDoc}
     */
    public SymbolPropertyTableIterator tableIterator() {
        return new SimpleSymbolPropertyTableIterator(map.values().iterator());
    }

    /**
     * {@inheritDoc}
     */
    public int numTables() {
        return map.size();
    }

    /**
     * {@inheritDoc}
     */
    public SymbolPropertyTable table(String name) throws 
    IllegalIDException, NullPointerException {
        if (name == null) {
            throw new NullPointerException("name is null.");
        }
        if (!map.containsKey(name)) {
            throw new IllegalIDException("No table found with name "
                    + name + ".");
        }
        return (SymbolPropertyTable) map.get(name);
    }

    /**
     * {@inheritDoc}
     */
    public Set names() {
        return map.keySet();
    }
}
