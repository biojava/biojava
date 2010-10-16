/*
 * SimpleRichLocationTest.java
 *
 * Created on April 1, 2007, 8:40 AM
 */

package org.biojavax.bio.seq;

import junit.framework.TestCase;
import org.biojavax.CrossRef;
import org.biojavax.SimpleCrossRef;
import org.biojavax.bio.seq.RichLocation.Strand;

/**
 *
 * @author George Waldon - equality
 */
public class SimpleRichLocationTest extends TestCase {
    
    SimpleRichLocation loc1;
    SimpleRichLocation loc2;
    SimpleRichLocation loc3;
    
    protected void setUp() throws Exception {
 
    }

    protected void tearDown() throws Exception {
    }
    
    /**Test of equals method.
     * Locations are equal if their term, min, max, strand, and crossref are
     * the same, and if their rank is the same too.
     */
    public void testEquals() {
        System.out.println("testEquals");
        CrossRef cr1 = new SimpleCrossRef("GenBank","A12345",3);
        CrossRef cr2 = new SimpleCrossRef("Another","A12345",3);
        Position p1 = new SimplePosition(23);
        Position p2 = new SimplePosition(24);
        
        
        loc1 = new SimpleRichLocation(p1, p2, 0,Strand.POSITIVE_STRAND,cr1); 
        loc2 = new SimpleRichLocation(p1, p2, 0,Strand.POSITIVE_STRAND,cr1);
        assertTrue(loc1.equals(loc2));
        
        Position p3 = new SimplePosition(22);
        loc3 = new SimpleRichLocation(p3, p2, 0,Strand.POSITIVE_STRAND,cr1);
        assertFalse(loc1.equals(loc3));
        
        p3 = new SimplePosition(25);
        loc3 = new SimpleRichLocation(p1, p3, 0,Strand.POSITIVE_STRAND,cr1);
        assertFalse(loc1.equals(loc3));
        
        loc3 = new SimpleRichLocation(p1, p2, 1,Strand.POSITIVE_STRAND,cr1);
        assertFalse(loc1.equals(loc3));
        
        loc3 = new SimpleRichLocation(p1, p2, 0,Strand.NEGATIVE_STRAND,cr1);
        assertFalse(loc1.equals(loc3));
        
        loc3 = new SimpleRichLocation(p1, p2, 0,Strand.POSITIVE_STRAND,cr2);
        assertFalse(loc1.equals(loc3));
    }
    
}
