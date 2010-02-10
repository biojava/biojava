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

package org.biojava.bio.seq;

import java.util.ArrayList;

import junit.framework.TestCase;

import org.biojava.bio.symbol.CircularLocation;
import org.biojava.bio.symbol.LocationTools;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;

/**
 * <code>CircularSequenceTest</code> tests are to ensure that the <code>CircularView</code>
 * class can be instantiated and that the coordinate system works
 *
 * @author Mark Schreiber
 */
public class CircularSequenceTest extends TestCase
{
    protected Sequence dna;
    protected CircularView circ;
    protected Symbol a = DNATools.a();
    protected Symbol c = DNATools.c();
    protected Symbol g = DNATools.g();
    protected Symbol t = DNATools.t();

    public CircularSequenceTest(String name)
    {
        super(name);
    }

    protected void setUp() throws Exception
    {

        dna = DNATools.createDNASequence("atcgctcaga","");
        circ = new CircularView(dna);
    }

    public void testLocations(){
      CircularLocation cloc = LocationTools.makeCircularLocation(1,5,dna.length());
      assertEquals("atcgc", cloc.symbols(dna).seqString());

      cloc = LocationTools.makeCircularLocation(1,10,dna.length());
      assertEquals("atcgctcaga", cloc.symbols(dna).seqString());

      cloc = LocationTools.makeCircularLocation(9,1,dna.length());
      assertTrue(cloc.get5PrimeEnd() == 9);
      assertEquals("gaa", cloc.symbols(dna).seqString());

      // a more tricky case
      ArrayList l = new ArrayList();
      //this is the most 5'
      l.add(LocationTools.makeCircularLocation(5,6, dna.length()));
      //this is not the most 5'
      l.add(LocationTools.makeCircularLocation(9,1,dna.length()));
      l.add(LocationTools.makeCircularLocation(3,3, dna.length()));

      cloc = (CircularLocation)LocationTools.union(l);
      assertTrue(cloc.get5PrimeEnd() == 5);
      assertEquals("ctgaac", cloc.symbols(dna).seqString());
    }


    public void testSubStr()
    {
        assertEquals("atc", circ.subStr(1,3));
        assertEquals("aat", circ.subStr(0,2));
        assertEquals("agaatc", circ.subStr(-2,3));
        assertEquals("gaatc", circ.subStr(9,3));
        assertEquals("gaatc", circ.subStr(9,13));
    }

    public void testSymbolAt()
    {
        assertEquals(a, circ.symbolAt(0));
        assertEquals(c, circ.symbolAt(13));
        assertEquals(g, circ.symbolAt(-1));
        assertEquals(t, circ.symbolAt(-4));
    }

    public void testSubList() throws Exception
    {
        SymbolList s = DNATools.createDNA("agaatc");
        //System.out.println(circ.subList(-2,3).seqString());
        assertTrue(s.equals(circ.subList(-2,3)));
        assertTrue(symsEq(circ.subList(-2,2), circ.subList(8,12)));
    }

    private boolean symsEq(SymbolList a, SymbolList b){
        if(a.length() != b.length()) return false;
        for(int i = 1; i <= a.length(); i++){
          if(a.symbolAt(i) != b.symbolAt(i)) return false;
        }
        return true;
    }

}
