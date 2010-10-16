/*
 * DummyCrossReferenceResolverTest.java
 * JUnit based test
 *
 * Created on November 11, 2005, 3:54 PM
 */

package org.biojavax;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

import org.biojava.bio.seq.DNATools;
import org.biojavax.bio.seq.InfinitelyAmbiguousSymbolList;

/**
 *
 * @author Mark Schreiber
 */
public class DummyCrossReferenceResolverTest extends TestCase {
    DummyCrossReferenceResolver res;
    CrossRef xref;
    
    public DummyCrossReferenceResolverTest(String testName) {
        super(testName);
        xref = new SimpleCrossRef("db", "AC123456", 2);
    }

    protected void setUp() throws Exception {
        res = new DummyCrossReferenceResolver();
    }

    protected void tearDown() throws Exception {
        res = null;
    }

    public static Test suite() {
        TestSuite suite = new TestSuite(DummyCrossReferenceResolverTest.class);
        
        return suite;
    }

    /**
     * Test of getRemoteSymbolList method, of class org.biojavax.DummyCrossReferenceResolver.
     */
    public void testGetRemoteSymbolList() {
        System.out.println("testGetRemoteSymbolList");
        assertTrue(res.getRemoteSymbolList(xref, DNATools.getDNA()) 
            instanceof InfinitelyAmbiguousSymbolList);
    }

    /**
     * Test of getRemoteBioEntry method, of class org.biojavax.DummyCrossReferenceResolver.
     */
    public void testGetRemoteBioEntry() {
        System.out.println("testGetRemoteBioEntry");
        
        assertNull(res.getRemoteBioEntry(xref));
    }
    
}
