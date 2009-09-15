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

import java.util.Iterator;
import java.util.Set;

import junit.framework.TestCase;

import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.utils.ChangeVetoException;
import org.biojavax.ga.Organism;
import org.biojavax.ga.exception.IllegalOrganismException;



/**
 * @author Mark Schreiber
 */
public class SimplePopulationTest extends TestCase {
  private SimplePopulation pop;
  private Organism o;

  public SimplePopulationTest(String s) {
    super(s);
  }

  protected void setUp() throws Exception{
    super.setUp();
    pop = new SimplePopulation("pop");

    o = new SimpleOrganism();
    o.setName("org");

    SymbolList[] csomes = new SymbolList[1];
    csomes[0] = DNATools.createDNA("aaaaagggggtttttccccc");
    o.setChromosomes(csomes);
    pop.addOrganism(o);
  }

  protected void tearDown() throws Exception{
    pop= null;
    o = null;
    super.tearDown();
  }

  public void testGetOrganismByName() {
    String name1=  "org";
    Organism organismRet = pop.getOrganismByName(name1);
    assertEquals( organismRet, o);
  }

  public void testGetOrganisms() {
    Set setRet = pop.getOrganisms();
    assertTrue( setRet.contains(o));
    assertTrue( setRet.size() == 1);

    try {
      pop.addOrganism(o.replicate("org2"));
    }
    catch (Exception ex) {
      fail(ex.getMessage());
    }

    // System.out.println(pop.getOrganisms().size());

    assertTrue( pop.getOrganisms().size() == 2);
  }
  public void testOrganisms() {

    Iterator it = pop.organisms();
    assertTrue( it.hasNext());
    assertTrue( it.next().equals(o));
  }

  public void addOrganism(){
    try {
      pop.addOrganism(o.replicate("org2"));
      assertTrue( pop.size() == 2);
    }
    catch (IllegalOrganismException ex) {
      fail(ex.getMessage());
    }catch (ChangeVetoException ex) {
      fail(ex.getMessage());
    }
  }

  public void addAndRemoveOrganisms(){
    Organism[] orgs = new SimpleOrganism[10];

    for (int i = 0; i < orgs.length; i++) {
      orgs[i] = o.replicate(o.getName()+i);
    }

    try {
      pop.addOrganisms(orgs);
      assertTrue( pop.size() == 11);

      pop.removeOrganisms(orgs);
      assertTrue( pop.size() == 1);

      pop.removeOrganism(o);
      assertTrue( pop.size() == 0);
    }
    catch (Exception ex) {
      fail(ex.getMessage());
    }
  }
}
