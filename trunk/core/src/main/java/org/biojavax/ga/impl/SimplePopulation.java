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

import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

import org.biojavax.ga.Organism;
import org.biojavax.ga.exception.IllegalOrganismException;

/**
 * <p>Simple concrete implementation of the <code>Population</code> interface</p>
 * <p>Internally the SimplePopulation store Organisms in a HashMap</p>
 * 
 * @author Mark Schreiber
 * @version 1.0
 * @since 1.5
 */

public final class SimplePopulation extends AbstractPopulation{
  private Map orgs;

  public SimplePopulation(String name) {
    super(name);
    this.orgs = new HashMap();
  }

  public SimplePopulation(){
    this("");
  }

  protected void addOrganismImpl(Organism orgToAdd) throws IllegalOrganismException{
    if(orgs.containsKey(orgToAdd.getName()))
       throw new IllegalOrganismException("All organisms in a population must have a unique name");
    orgs.put(orgToAdd.getName(),orgToAdd);
  }

  protected void removeOrganismImpl(Organism orgToRemove){
    orgs.remove(orgToRemove.getName());
  }

  protected void removeAllOrganismsImpl(){
    orgs = new HashMap();
  }

  public Organism getOrganismByName(String name){
    return (Organism)orgs.get(name);
  }

  public int size() {
    return orgs.size();
  }
  public Iterator organisms() {
    return orgs.values().iterator();
  }
  public Set getOrganisms() {
    return new HashSet(orgs.values());
  }

}
