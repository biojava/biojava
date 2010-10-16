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
package org.biojavax.ga.impl;

import junit.framework.TestCase;

import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.SymbolList;
import org.biojavax.ga.Organism;



/**
 * @author Mark Schreiber
 */
public class SimpleOrganismTest extends TestCase {
  private SimpleOrganism o;

  public SimpleOrganismTest(String s) {
    super(s);
  }

  protected void setUp() throws Exception{
    super.setUp();
    o = new SimpleOrganism();
    o.setName("org");

    SymbolList[] csomes = new SymbolList[1];
    csomes[0] = DNATools.createDNA("aaaaagggggtttttccccc");
    o.setChromosomes(csomes);
  }

  protected void tearDown() throws Exception{
    o = null;
    super.tearDown();
  }

  public void testIsHaploid() {
    assertTrue(o.isHaploid());
  }

  public void testReplicate() {
    Organism o2 = o.replicate("org2");

    assertEquals( o.getName(), "org");
    assertEquals( o2.getName(), "org2");

    assertEquals( o.getChromosomes()[0], o2.getChromosomes()[0]);
  }
}
