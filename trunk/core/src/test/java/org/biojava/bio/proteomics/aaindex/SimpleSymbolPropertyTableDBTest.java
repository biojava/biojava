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
 *    SimpleSymbolPropertyTableDBTest.java
 */
package org.biojava.bio.proteomics.aaindex;

import java.util.ArrayList;
import java.util.Set;

import junit.framework.TestCase;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.ProteinTools;
import org.biojava.bio.seq.db.IllegalIDException;
import org.biojava.bio.symbol.SimpleSymbolPropertyTable;
import org.biojava.bio.symbol.SymbolPropertyTable;

/**
 * Test class for {@link org.biojava.bio.proteomics.aaindex.SimpleSymbolPropertyTableDB}.
 * @author <a href="mailto:Martin.Szugat@GMX.net">Martin Szugat</a>
 * @version $Revision$
 */
public class SimpleSymbolPropertyTableDBTest extends TestCase {

    /**
     * Main entry point.
     * @param args command line arguments
     */
    public static void main(String[] args) {
        junit.textui.TestRunner.run(SimpleSymbolPropertyTableDBTest.class);
    }

    /**
     * Test constructor to create an empty database.
     */
    public void testEmptyConstructor() {
        SymbolPropertyTableDB db = new SimpleSymbolPropertyTableDB();
        assertEquals("Database is not empty.", 0, db.numTables());
        assertFalse("Database does not return an empty iterator.", 
                db.tableIterator().hasNext());
    }


    /**
     * Test constructor to create a non-empty database.
     * @throws BioException if iterator fails.
     */
    public void testNonEmptyConstructor() throws BioException {
        SimpleSymbolPropertyTableDB dbInit = new SimpleSymbolPropertyTableDB();
        SymbolPropertyTable t1 = new SimpleSymbolPropertyTable(
                ProteinTools.getAlphabet(), "protein");
        SymbolPropertyTable t2 = new SimpleSymbolPropertyTable(
                DNATools.getDNA(), "dna");
        dbInit.addTable(t1);
        dbInit.addTable(t2);
        SymbolPropertyTableDB db = new SimpleSymbolPropertyTableDB(
                dbInit.tableIterator());
        assertEquals("Database has wrong number of tables.", 2, db.numTables());
        SymbolPropertyTableIterator iterator = db.tableIterator();
        ArrayList tables = new ArrayList(2);
        while (iterator.hasNext()) {
            tables.add(iterator.nextTable());
        }
        assertEquals("Iterator returned wrong number of tables.", 
                2, tables.size());
        assertEquals("Iterator returned wrong table.", t1, tables.get(0));
        assertEquals("Iterator returned wrong table.", t2, tables.get(1));
    }
    
    /**
     * Test for {@link SimpleSymbolPropertyTableDB#addTable(SymbolPropertyTable)}.
     * @throws IllegalIDException if {@link SymbolPropertyTableDB#table(String)}
     * fails.
     */
    public void testAddTable() throws IllegalIDException {
        SimpleSymbolPropertyTableDB db = new SimpleSymbolPropertyTableDB();
        SymbolPropertyTable t1 = new SimpleSymbolPropertyTable(
                ProteinTools.getAlphabet(), "protein");
        SymbolPropertyTable t2 = new SimpleSymbolPropertyTable(
                DNATools.getDNA(), "dna");
        db.addTable(t1);
        assertEquals("Database has wrong number of tables.", 1, db.numTables());
        db.addTable(t2);
        assertEquals("Database has wrong number of tables.", 2, db.numTables());
        assertEquals("Database returned wrong table.", t1, db.table("protein"));
        assertEquals("Database returned wrong table.", t2, db.table("dna"));
        Set names = db.names();
        assertTrue("Table is missing.", names.contains("protein"));
        assertTrue("Table is missing.", names.contains("dna"));
        assertEquals("Database has wrong number of tables.", 2, names.size());
        try {
            db.addTable(null);
            fail("addTable must throw NullPointerException.");
        } catch (NullPointerException e) {
            e.printStackTrace();
        }
    }
    
    /**
     * Test for {@link SimpleSymbolPropertyTableDB#table(String)}.
     */
    public void testTable() {
        SimpleSymbolPropertyTableDB db = new SimpleSymbolPropertyTableDB();
        try {
            db.table(null);
            fail("table must throw NullPointerException.");
        } catch (NullPointerException e) {
            e.printStackTrace();
        } catch (IllegalIDException e) {
            fail("table throwed IllegalIDException instead of NullPointerException.");
        }
        try {
            db.table("test");
            fail("table must throw IllegalIDException.");
        } catch (IllegalIDException e) {
            e.printStackTrace();
        }
    }
}
