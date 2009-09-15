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

import junit.framework.TestCase;

import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.impl.SimpleGappedSequence;
import org.biojava.bio.seq.impl.SimpleSequence;
import org.biojava.bio.seq.io.SymbolTokenization;

/**
 * Test for gapped symbol lists.
 *
 * @author Matthew Pocock
 * @since 1.3
 */

public class GappedSymbolListTest extends TestCase {
  private SymbolList symList;
  private SymbolList symList1;
  private SymbolList symList2;
  private SymbolList symList3;
  private SymbolList symList4;

  public GappedSymbolListTest(String name) {
    super(name);
  }

  protected void setUp()
  throws Exception {
    FiniteAlphabet dna = (FiniteAlphabet) AlphabetManager.alphabetForName("DNA");
    SymbolTokenization tok = dna.getTokenization("token");
    symList  = new SimpleSymbolList(tok,"gtgtggga");
    symList1 = new SimpleSymbolList(tok,"gtg-tggga");
    symList2 = new SimpleSymbolList(tok,"gtg--tggga");
    symList3 = new SimpleSymbolList(tok,"gtg---tggga");
    symList4 = new SimpleSymbolList(tok,"gtg----tggga");
  }

  public void testIndividualInsertBlockRemove()
  throws Exception {
    GappedSymbolList gsl = new SimpleGappedSymbolList(symList);
    assertTrue(SymbolUtils.compareSymbolLists(gsl, symList));

    gsl.addGapInSource(4);
    assertTrue(SymbolUtils.compareSymbolLists(gsl, symList1));

    gsl.addGapInSource(4);
    assertTrue(SymbolUtils.compareSymbolLists(gsl, symList2));

    gsl.addGapInSource(4);
    assertTrue(SymbolUtils.compareSymbolLists(gsl, symList3));

    gsl.addGapInSource(4);
    assertTrue(SymbolUtils.compareSymbolLists(gsl, symList4));

    gsl.removeGaps(4,4);
    assertTrue(SymbolUtils.compareSymbolLists(gsl, symList));
  }

  public void testBlockedInsertIndividualRemove()
  throws Exception {
    GappedSymbolList gsl = new SimpleGappedSymbolList(symList);
    assertTrue("Gapped same as ungapped:\n" + gsl.seqString() + "\n" + symList.seqString(), SymbolUtils.compareSymbolLists(gsl, symList));

    gsl.addGapsInSource(4,4);
    assertTrue("Four gaps:\n" + gsl.seqString() + "\n" + symList4.seqString(), SymbolUtils.compareSymbolLists(gsl, symList4));

    gsl.removeGap(4);
    assertTrue("Three gaps:\n" + gsl.seqString() + "\n" + symList3.seqString(), SymbolUtils.compareSymbolLists(gsl, symList3));

    gsl.removeGap(4);
    assertTrue("Two gaps:\n" + gsl.seqString() + "\n" + symList2.seqString(), SymbolUtils.compareSymbolLists(gsl, symList2));

    gsl.removeGap(4);
    assertTrue("One gap:\n" + gsl.seqString() + "\n" + symList1.seqString(), SymbolUtils.compareSymbolLists(gsl, symList1));

    gsl.removeGap(4);
    assertTrue("All gaps removed:\n" + gsl.seqString() + "\n" + symList.seqString(), SymbolUtils.compareSymbolLists(gsl, symList));
  }

  public void testBlockedInsertBlockRemove()
  throws Exception {
    GappedSymbolList gsl = new SimpleGappedSymbolList(symList);
    assertTrue(SymbolUtils.compareSymbolLists(gsl, symList));

    gsl.addGapsInSource(4,4);
    assertTrue(SymbolUtils.compareSymbolLists(gsl, symList4));

    gsl.removeGaps(4,4);
    assertTrue(SymbolUtils.compareSymbolLists(gsl, symList));
  }

  public void testBeginningGap()
  throws Exception {
    GappedSymbolList gsl = new SimpleGappedSymbolList(symList);
    gsl.addGapInSource(1);
    assertTrue(
      "Begining gap is the empty gap: " +
      Alphabet.EMPTY_ALPHABET.getGapSymbol() + " vs " +
      gsl.symbolAt(1),
      Alphabet.EMPTY_ALPHABET.getGapSymbol() == gsl.symbolAt(1)
    );
  }

  public void testBeginningGaps()
  throws Exception {
    GappedSymbolList gsl = new SimpleGappedSymbolList(symList);
    gsl.addGapsInSource(1, 4);
    for(int i = 1; i <=4; i++) {
      assertTrue(
        "Begining gap " + i + " is the empty gap: " +
        Alphabet.EMPTY_ALPHABET.getGapSymbol() + " vs " +
        gsl.symbolAt(i),
        Alphabet.EMPTY_ALPHABET.getGapSymbol() == gsl.symbolAt(i)
      );
    }
  }

  public void testTrailingGap()
  throws Exception {
    GappedSymbolList gsl = new SimpleGappedSymbolList(symList);
    gsl.addGapInSource(symList.length() + 1);
    assertTrue(
      "Trailing gap is the empty gap: " +
      Alphabet.EMPTY_ALPHABET.getGapSymbol() + " vs " +
      gsl.symbolAt(gsl.length()),
      Alphabet.EMPTY_ALPHABET.getGapSymbol() == gsl.symbolAt(gsl.length())
    );
  }

  public void testTrailingGaps()
  throws Exception {
    GappedSymbolList gsl = new SimpleGappedSymbolList(symList);
    gsl.addGapsInSource(symList.length() + 1, 4);
    for(int i = 1; i <=4; i++) {
      assertTrue(
        "Trailing gap " + i + " is the empty gap: " +
        Alphabet.EMPTY_ALPHABET.getGapSymbol() + " vs " +
        gsl.symbolAt(gsl.length() + 1 - i),
        Alphabet.EMPTY_ALPHABET.getGapSymbol() == gsl.symbolAt(gsl.length() + 1 - i)
      );
    }
  }

  public void testInternalGap()
  throws Exception {
    GappedSymbolList gsl = new SimpleGappedSymbolList(symList);
    gsl.addGapInSource(4);
    assertTrue(
      "Internal gap is the apropreate gap: " +
      gsl.getAlphabet().getGapSymbol() + " vs " +
      gsl.symbolAt(4),
      gsl.getAlphabet().getGapSymbol() == gsl.symbolAt(4)
    );
  }

  public void testInternalGaps()
  throws Exception {
    GappedSymbolList gsl = new SimpleGappedSymbolList(symList);
    gsl.addGapsInSource(4, 4);
    for(int i = 4; i < 8; i++) {
      assertTrue(
        "Internal gap " + i + " is the apropreate gap: " +
        gsl.getAlphabet().getGapSymbol() + " vs " +
        gsl.symbolAt(i),
        gsl.getAlphabet().getGapSymbol() == gsl.symbolAt(i)
      );
    }
  }

  public void testDavidBug()
  throws Exception {
    SymbolTokenization dnaToke1 = DNATools.getDNA().getTokenization("token");
    SymbolList symbolList = new SimpleSymbolList(dnaToke1, "ACTGGACCTAAGG");
    Sequence sequence = new SimpleSequence(symbolList, "test", "test", null);
    SimpleGappedSequence gappedSequence = new SimpleGappedSequence(sequence);

    gappedSequence.addGapsInView(4, 4);
    gappedSequence.removeGap(7);
    gappedSequence.removeGaps(4, 3);
    gappedSequence.addGapsInView(7, 2);
    gappedSequence.addGapsInView(9, 3);
    gappedSequence.addGapsInView(12, 2);
    gappedSequence.addGapsInView(14, 3);
    gappedSequence.addGapsInView(17, 2);
    gappedSequence.removeGap(18);
    gappedSequence.removeGaps(11, 6);

    String seqStr = gappedSequence.seqString();

    seqStr=seqStr==null?null:seqStr;//trick
  }
}
