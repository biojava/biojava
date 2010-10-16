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

import org.biojava.utils.AbstractChangeable;
import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeSupport;
import org.biojava.utils.ChangeVetoException;
import org.biojavax.ga.GeneticAlgorithm;
import org.biojavax.ga.Organism;
import org.biojavax.ga.Population;
import org.biojavax.ga.exception.IllegalOrganismException;

/**
 * Most Population implementations will want to inherit from here.
 * This class doesn't define how Organims are stored or accessed so inheriting classes
 * can define that themselves.
 *
 * @author Mark Schreiber
 * @version 1.0
 * @since 1.5
 */

public abstract class AbstractPopulation extends AbstractChangeable implements Population {
  String name;

  public AbstractPopulation() {
    this("");
  }

  public AbstractPopulation(String name){
    this.name = name;
  }

  public String getName() {
    return name;
  }
  public final void setName(String name) throws ChangeVetoException {
    if(!hasListeners()){
      this.name = name;
   }else{
     ChangeEvent ce = new ChangeEvent(this,
                                      Population.NAME,
                                      name,
                                      this.name
                                      );
     ChangeSupport changeSupport = super.getChangeSupport(GeneticAlgorithm.POPULATION);
     synchronized(changeSupport){
       changeSupport.firePreChangeEvent(ce);
       this.name = name;
       changeSupport.firePostChangeEvent(ce);
     }
    }
  }

  public final void addOrganism(Organism org) throws ChangeVetoException, IllegalOrganismException {
    if(!hasListeners()){
      addOrganismImpl(org);
   }else{
     ChangeEvent ce = new ChangeEvent(this,
                                      Population.ORGANISMS,
                                      org,
                                      getOrganisms()
                                      );
     ChangeSupport changeSupport = super.getChangeSupport(GeneticAlgorithm.POPULATION);
     synchronized(changeSupport){
       changeSupport.firePreChangeEvent(ce);
       addOrganismImpl(org);
       changeSupport.firePostChangeEvent(ce);
     }
    }
  }

  protected abstract void addOrganismImpl(Organism org) throws IllegalOrganismException;

  public final void addOrganisms(Organism[] orgs) throws ChangeVetoException, IllegalOrganismException {
    for (int i = 0; i < orgs.length; i++) {
      addOrganism(orgs[i]);
    }
  }

  public final void addOrganisms(Set orgs)throws ChangeVetoException, IllegalOrganismException {
    for (Iterator i = orgs.iterator(); i.hasNext(); ) {
      Object o = i.next();

      if (o instanceof Organism) {
        addOrganism((Organism)o);
      }
      else {
        throw new IllegalOrganismException(
            "Attempting to add a non Organism object to a population, object is: "+
            o.getClass().getName()
            );
      }
    }
  }
  public final void addOrganisms(Population orgs) throws ChangeVetoException, IllegalOrganismException {
    for (Iterator i = orgs.organisms(); i.hasNext(); ) {
      Object o = i.next();

      if (o instanceof Organism) {
        addOrganism((Organism)o);
      }
      else {
        throw new IllegalOrganismException(
              "Attempting to add a non Organism object to a population, object is: "+
              o.getClass().getName()
            );
      }

    }
  }





  public final void removeOrganisms(Organism[] orgs)throws ChangeVetoException{
    for (int i = 0; i < orgs.length; i++) {
      removeOrganismImpl(orgs[i]);
    }
  }

  public final void removeOrganisms(Set orgs) throws ChangeVetoException{
    for (Iterator i = orgs.iterator(); i.hasNext(); ) {
      removeOrganismImpl((Organism)i.next());
    }
  }

  public final void removeAllOrganisms() throws ChangeVetoException{
    if(!hasListeners()){
      removeAllOrganismsImpl();
   }else{
     ChangeEvent ce = new ChangeEvent(this,
                                      Population.ORGANISMS,
                                      null,
                                      getOrganisms()
                                      );
     ChangeSupport changeSupport = super.getChangeSupport(GeneticAlgorithm.POPULATION);
     synchronized(changeSupport){
       changeSupport.firePreChangeEvent(ce);
       removeAllOrganismsImpl();
       changeSupport.firePostChangeEvent(ce);
     }
    }
  }


  public final void removeOrganism(Organism org) throws ChangeVetoException {
    if(!hasListeners()){
      removeOrganismImpl(org);
   }else{
     ChangeEvent ce = new ChangeEvent(this,
                                      Population.ORGANISMS,
                                      org,
                                      getOrganisms()
                                      );
     ChangeSupport changeSupport = super.getChangeSupport(GeneticAlgorithm.POPULATION);
     synchronized(changeSupport){
       changeSupport.firePreChangeEvent(ce);
       removeOrganismImpl(org);
       changeSupport.firePostChangeEvent(ce);
     }
    }
  }

  protected abstract void removeOrganismImpl(Organism org);
  protected abstract void removeAllOrganismsImpl();
}