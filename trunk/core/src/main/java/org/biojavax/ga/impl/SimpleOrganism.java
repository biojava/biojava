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
import org.biojava.bio.symbol.SimpleSymbolList;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.utils.ChangeVetoException;
import org.biojavax.ga.Organism;

/**
 * A Simple Haploid Organism implementation
 *
 * @author Mark Schreiber
 * @version 1.0
 * @since 1.5
 */

public class SimpleOrganism extends AbstractOrganism {

	public SimpleOrganism(){
    super();
  }
  public SimpleOrganism(String name){
    super();
    try {
      setName(name);
    }
    catch (ChangeVetoException ex) {
      //can't happen till after the organism has been made
    }
  }

  public SimpleOrganism(Organism org, String name){
    super(org, name);
  }

  protected void setChromImpl(SymbolList[] chromosomes) {
    this.chromosomes = chromosomes;
  }

  /**
   * Simple Organisms are Halpoid
   * @return true
   */
  public boolean isHaploid() {
    return true;
  }

  public Organism replicate(String name){
    SimpleOrganism o = new SimpleOrganism(name);
    SymbolList[] symls = new SymbolList[this.getChromosomes().length];
    for (int i = 0; i < symls.length; i++) {
      symls[i] = new SimpleSymbolList(this.getChromosomes()[i]);
    }

    o.setChromImpl(symls);
    return o;
  }

}