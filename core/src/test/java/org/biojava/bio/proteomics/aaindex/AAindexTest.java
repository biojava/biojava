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
 *    AAindexTest.java
 */
package org.biojava.bio.proteomics.aaindex;

import junit.framework.TestCase;

import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.ProteinTools;
import org.biojava.bio.symbol.IllegalSymbolException;

/**
 * Test class for {@link org.biojava.bio.proteomics.aaindex.AAindex}.
 * @author <a href="mailto:Martin.Szugat@GMX.net">Martin Szugat</a>
 * @version $Revision$
 */
public class AAindexTest extends TestCase {

    /**
     * Main entry point.
     * @param args command line arguments
     */
    public static void main(String[] args) {
        junit.textui.TestRunner.run(AAindexTest.class);
    }
    
    /**
     * AAindex object.
     */
    private AAindex aaindex = null;

    /**
     * {@inheritDoc}
     */
    protected void setUp() throws Exception {
        super.setUp();
        aaindex = new AAindex("test");
    }

    /**
     * Test for {@link AAindex#accessionNumber()}.
     */
    public void testAccessionNumber() {
        assertEquals("test", aaindex.accessionNumber());
    }

    /**
     * Test for {@link AAindex#getArticleAuthors()}.
     */
    public void testGetArticleAuthors() {
        assertNull(aaindex.getArticleAuthors());
    }

    /**
     * Test for {@link AAindex#setArticleAuthors(String)}.
     */
    public void testSetArticleAuthors() {
        aaindex.setArticleAuthors("authors");
        assertEquals("authors", aaindex.getArticleAuthors());
    }

    /**
     * Test for {@link AAindex#getComment()}.
     */
    public void testGetComment() {
        assertNull(aaindex.getComment());
    }

    /**
     * Test for {@link AAindex#setComment(String)}.
     */
    public void testSetComment() {
        aaindex.setComment("comment");
        assertEquals("comment", aaindex.getComment());
    }

    /**
     * Test for {@link AAindex#getArticleTitle()}.
     */
    public void testGetArticleTitle() {
        assertNull(aaindex.getArticleTitle());
    }

    /**
     * Test for {@link AAindex#setArticleTitle(String)}.
     */
    public void testSetArticleTitle() {
        aaindex.setArticleTitle("title");
        assertEquals("title", aaindex.getArticleTitle());
    }

    /**
     * Test for {@link AAindex#getDescription()}.
     */
    public void testGetDescription() {
        assertNull(aaindex.getDescription());
    }

    /**
     * Test for {@link AAindex#setDescription(String)}.
     */
    public void testSetDescription() {
        aaindex.setDescription("description");
        assertEquals("description", aaindex.getDescription());
    }

    /**
     * Test for {@link AAindex#getJournalReference()}
     */
    public void testGetJournalReference() {
        assertNull(aaindex.getJournalReference());
    }

    /**
     * Test for {@link AAindex#setJournalReference(String)}.
     */
    public void testSetJournalReference() {
        aaindex.setJournalReference("journal");
        assertEquals("journal", aaindex.getJournalReference());
    }

    /**
     * Test for {@link AAindex#getLITDBEntryNumbers()}. 
     */
    public void testGetLITDBEntryNumbers() {
        assertNull(aaindex.getLITDBEntryNumbers());
    }

    /**
     * Test for {@link AAindex#setLITDBEntryNumbers(String[])}.
     */
    public void testSetLITDBEntryNumbers() {
        aaindex.setLITDBEntryNumbers(new String[] {"lit"});
        assertEquals("lit", aaindex.getLITDBEntryNumbers()[0]);
    }

    /**
     * Test for {@link AAindex#similarEntries()}.
     */
    public void testSimilarEntries() {
        assertEquals(0, aaindex.similarEntries().size());
        aaindex.similarEntries().put("test", new Double(0.0));
        assertEquals(new Double(0.0), aaindex.similarEntries().get("test"));
    }

    /**
     * Test for {@link org.biojava.bio.symbol.SymbolPropertyTable#getName()}.
     */
    public void testGetName() {
        assertEquals("test", aaindex.getName());
    }

    /**
     * Test for {@link org.biojava.bio.symbol.SymbolPropertyTable#getAlphabet()}.
     */
    public void testGetAlphabet() {
        assertEquals(ProteinTools.getAlphabet(), aaindex.getAlphabet());
    }

    /**
     * Test for {@link org.biojava.bio.symbol.SymbolPropertyTable#getDoubleValue(org.biojava.bio.symbol.Symbol)}.
     * @throws IllegalSymbolException if a symbol is illegal.
     */
    public void testDoubleValue() throws IllegalSymbolException {
        aaindex.setDoubleProperty(ProteinTools.gln(), "1.0");
        assertEquals(1.0, aaindex.getDoubleValue(ProteinTools.gln()), 0.0);
        try {
            aaindex.getDoubleValue(DNATools.a());
            fail("IllegalSymbolException not thrown.");
        } catch (IllegalSymbolException e) {
            e.printStackTrace();
        }
    }
}
