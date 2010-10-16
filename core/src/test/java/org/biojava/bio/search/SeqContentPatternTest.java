package org.biojava.bio.search;

import junit.framework.TestCase;

import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.SymbolList;

public class SeqContentPatternTest
extends TestCase {
  public void testTooShort()
  throws IllegalSymbolException, IllegalAlphabetException {
    SymbolList zero = DNATools.createDNA("");
    SymbolList one = DNATools.createDNA("a");
    SymbolList five = DNATools.createDNA("aaaaa");
    SymbolList sixMatch = DNATools.createDNA("aaaaaa");
    SymbolList sixNoMatch = DNATools.createDNA("tttttt");
    SymbolList sevenMatch1 = DNATools.createDNA("aaaaaat");
    SymbolList sevenMatch2 = DNATools.createDNA("taaaaaa");
    SymbolList sevenMatchAll = DNATools.createDNA("aaaaaaa");
    SymbolList sevenNoMatch = DNATools.createDNA("ttttttt");


    SeqContentPattern scp = new SeqContentPattern(DNATools.getDNA());
    scp.setLength(6);
    scp.setMinCounts(DNATools.a(), 6);

    BioMatcher scm;

    scm = scp.matcher(zero);
    assertFalse("No hits of length 6 in seq length 0", scm.find());

    scm = scp.matcher(one);
    assertFalse("No hits of length 6 in seq length 1", scm.find());

    scm = scp.matcher(five);
    assertFalse("No hits of length 6 in seq length 5", scm.find());

    scm = scp.matcher(sixMatch);
    assertTrue("Hit found in matching 6", scm.find());
    assertFalse("Just the one hit", scm.find());

    scm = scp.matcher(sixNoMatch);
    assertFalse("No hit in non-matching 6", scm.find());

    scm = scp.matcher(sevenMatch1);
    assertTrue("Hit in 7-1", scm.find());
    assertEquals("Hit in 7-1 at 1", 1, scm.start());
    assertFalse("No more hits in 7-1", scm.find());

    scm = scp.matcher(sevenMatch2);
    assertTrue("Hit in 7-2", scm.find());
    assertEquals("Hit in 7-2 at 2", 2, scm.start());
    assertFalse("No more hits in 7-2", scm.find());

    scm = scp.matcher(sevenMatchAll);
    assertTrue("Hit in 7-all", scm.find());
    assertEquals("Hit in 7-all at 1", 1, scm.start());
    assertTrue("Hit in 7-all", scm.find());
    assertEquals("Hit in 7-all at 2", 2, scm.start());

    scm = scp.matcher(sevenNoMatch);
    assertFalse("No hit in seven no hits", scm.find());
  }

  // we should have tests for various boundary conditions
  // on the counts
  // also - check that obviously impossible count ranges never match
}
