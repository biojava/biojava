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

package org.biojava.bio.symbol;

import java.util.Iterator;

import junit.framework.TestCase;

import org.biojava.bio.seq.DNATools;

/**
 * Tests for SimpleSymbolList
 *
 * @author Thomas Down
 * @since 1.2
 */

public class SimpleSymbolListTest extends TestCase {
    protected SymbolList sl1;
    protected SymbolList sl2;

    public SimpleSymbolListTest(String name) {
        super(name);
    }

    protected void setUp()
        throws Exception
    {
        sl1 = DNATools.createDNA("gataca");
        sl2 = DNATools.createDNA("atgga");
    }

    public void testCopyConstructor()
        throws Exception
    {
        SimpleSymbolList tsl = new SimpleSymbolList(sl1);
        assertTrue(compareSymbolList(tsl, sl1));
    }

   public void testEquals() throws Exception{
     SymbolList sl = DNATools.createDNA("gataca");
     assertTrue(sl.equals(sl1));
     assertTrue(sl1.equals(sl));
     assertTrue(! sl.equals(sl2));
   }

    public void testAddSymbol()
        throws Exception
    {
        SimpleSymbolList tsl = new SimpleSymbolList(sl1);
        tsl.addSymbol(DNATools.a());
        assertTrue(compareSymbolList(tsl, DNATools.createDNA("gatacaa")));
    }

    public void testInsertAtStart()
        throws Exception
    {
        SimpleSymbolList tsl = new SimpleSymbolList(sl1);
        tsl.edit(new Edit(1, 0, sl2));
        assertTrue(compareSymbolList(tsl, DNATools.createDNA("atggagataca")));
    }

    public void testInsertInMiddle()
        throws Exception
    {
        SimpleSymbolList tsl = new SimpleSymbolList(sl1);
        tsl.edit(new Edit(4, 0, sl2));
        assertTrue(compareSymbolList(tsl, DNATools.createDNA("gatatggaaca")));
    }

    public void testInsertAtEnd() throws Exception {
      SimpleSymbolList ts1 = new SimpleSymbolList(sl1);
      ts1.edit(new Edit(7,0,sl2));
      assertTrue(compareSymbolList(ts1, DNATools.createDNA("gatacaatgga")));
      ts1.edit(new Edit(ts1.length()+1, 0, sl2));
      assertTrue(compareSymbolList(ts1, DNATools.createDNA("gatacaatggaatgga")));
    }

    public void testDeletion()
        throws Exception
    {
        SimpleSymbolList tsl = new SimpleSymbolList(sl1);
        tsl.edit(new Edit(4, 2, SymbolList.EMPTY_LIST));
        assertTrue(compareSymbolList(tsl, DNATools.createDNA("gata")));
    }

    public void testReplacement()
        throws Exception
    {
        SimpleSymbolList tsl = new SimpleSymbolList(sl1);
        tsl.edit(new Edit(4, 2, sl2));
        assertTrue(compareSymbolList(tsl, DNATools.createDNA("gatatggaa")));
    }

    public void testEditSuperlist()
        throws Exception
    {
        SimpleSymbolList tsl = new SimpleSymbolList(sl1);
        SymbolList sub_tsl = tsl.subList(2, 6);
        tsl.edit(new Edit(4, 2, sl2));
        assertTrue(compareSymbolList(tsl, DNATools.createDNA("gatatggaa")));
        assertTrue(compareSymbolList(sub_tsl, DNATools.createDNA("ataca")));
    }

    public void testEditSublist()
        throws Exception
    {
        SimpleSymbolList tsl = new SimpleSymbolList(sl1);
        SymbolList sub_tsl = tsl.subList(2, 6);
        sub_tsl.edit(new Edit(3, 2, sl2));
        assertTrue(compareSymbolList(sub_tsl, DNATools.createDNA("atatggaa")));
        assertTrue(compareSymbolList(tsl, DNATools.createDNA("gataca")));
    }

    private boolean compareSymbolList(SymbolList sl1, SymbolList sl2) {
        if (sl1.length() != sl2.length()) {
            return false;
        }

        Iterator si1 = sl1.iterator();
        Iterator si2 = sl2.iterator();
        while (si1.hasNext()) {
            if (! (si1.next() == si2.next())) {
                return false;
            }
        }

        return true;
    }
}

