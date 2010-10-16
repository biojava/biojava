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
 */


package org.biojava.bio.alignment;

import java.util.ArrayList;
import java.util.List;

import junit.framework.TestCase;

import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.GappedSequence;
import org.biojava.bio.seq.ProteinTools;
import org.biojava.bio.symbol.BasisSymbol;
import org.biojava.bio.symbol.LocationTools;

/* @author Lachlan Coin */
public class FlexibleAlignmentTest extends TestCase {

  public FlexibleAlignmentTest(String name){
    super(name);
  }

    //public static void main(String[] args) throws Exception{
    //	FlexibleAlignmentTest aat = new FlexibleAlignmentTest();
    //	aat.setUp();
    //	aat.testDNA();
    //	aat.testProtein();
    //}


    final static String[] alignment = new String[] {"A-C",
                                       "AGC",
                                       "A-A"};


    final static String[] names = new String[] {"MOUSE", "HUMAN","SCHPO"};

    final static FlexibleAlignment alignDNA = parse(names, alignment, true);
    final static FlexibleAlignment alignProt = parse(names, alignment, false);

    protected void setUp() throws Exception{

    }


  public void testDNA()
      throws Exception
  {
      List syms = ((BasisSymbol) alignDNA.symbolAt(2)).getSymbols();
      assertEquals(syms.get(0), DNATools.getDNA().getGapSymbol());
      assertEquals(syms.get(1), DNATools.g());
  }


  public void testProtein()
      throws Exception
  {
      List syms = ((BasisSymbol) alignProt.symbolAt(2)).getSymbols();
      assertEquals(syms.get(0), ProteinTools.getAlphabet().getGapSymbol());
  }



    private static FlexibleAlignment parse(String[] names, String[] alignment, boolean dna)
    {
        try
        {
            List sequences = new ArrayList();
            for(int i=0; i<alignment.length; i++)
            {
                GappedSequence seq;
                if(dna)
                    seq = DNATools.createGappedDNASequence(alignment[i],names[i]);
                else
                    seq = ProteinTools.createGappedProteinSequence(alignment[i],names[i]);
                AlignmentElement ae = new SimpleAlignmentElement(names[i], seq, LocationTools.makeLocation(1, alignment[i].length()));
                //         System.out.println(seq.seqString());
                sequences.add(ae);
            }
            FlexibleAlignment al = new FlexibleAlignment(sequences);
            return al;
        }
        catch (Throwable t)
        {
            t.printStackTrace();
            //System.exit(0);
            return null;
        }
    }



}
