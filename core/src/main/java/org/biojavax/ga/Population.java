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


package org.biojavax.ga;

import java.util.Iterator;
import java.util.Set;

import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.Changeable;
import org.biojavax.ga.exception.IllegalOrganismException;

/**
 * A collection of GA organisms
 *
 * @author Mark Schreiber
 * @version 1.0
 * @since 1.5
 */

public interface Population extends Changeable{

  /**
   * @return the name of the population
   */
  public String getName();

  /**
   * Sets the name of the population
   * @param name set the name to this.
   * @throws ChangeVetoException if the name may not be changed
   */
  public void setName(String name) throws ChangeVetoException;

  /**
   * Adds an Organism to the Population
   * @param org the organism
   * @throws ChangeVetoException
   * @throws IllegalOrganismException if for some reason the organism is invalid
   */
  public void addOrganism(Organism org) throws ChangeVetoException, IllegalOrganismException;

  /**
   * Adds several organisms to the population
   * @param orgs the organisms to add
   * @throws ChangeVetoException
   * @throws IllegalOrganismException if for some reason the organism is invalid
   */
  public void addOrganisms(Organism[] orgs)throws ChangeVetoException, IllegalOrganismException;

  /**
   * Adds several organisms to the population
   * @param orgs the organisms to add
   * @throws ChangeVetoException
   * @throws IllegalOrganismException if for some reason the organism is invalid
   */
  public void addOrganisms(Set orgs)throws ChangeVetoException, IllegalOrganismException;

  /**
   * Adds the residents of one population to this one
   * @param orgs the population to add
   * @throws ChangeVetoException
   * @throws IllegalOrganismException if for some reason the organism is invalid
   */
  public void addOrganisms(Population orgs)throws ChangeVetoException, IllegalOrganismException;

  /**
   * Kills off the organism
   * @param org the organism to kill
   * @throws ChangeVetoException
   */
  public void removeOrganism(Organism org) throws ChangeVetoException;

  /**
   * Removes all the <code>Organisms</code> in <code>orgs</code>
   * @param orgs the <code>Organisms</code> to remove.
   * @throws ChangeVetoException if the change is vetoed
   */
  public void removeOrganisms(Organism[] orgs) throws ChangeVetoException;

  /**
   * Removes all the <code>Organisms</code> in <code>orgs</code>
   * @param orgs the <code>Organisms</code> to remove.
   * @throws ChangeVetoException if the change is vetoed
   */
  public void removeOrganisms(Set orgs) throws ChangeVetoException;

  /**
   * Removes all the <code>Organisms</code> in this <code>Population</code>
   * @throws ChangeVetoException if the change is vetoed
   */
  public void removeAllOrganisms() throws ChangeVetoException;

  /**
   * Gets the specified organism
   * @param name the name of the organism to retreive
   * @return the organism named or null if that organism doesn't exist.
   */
  public Organism getOrganismByName(String name);

  /**
   * Gets the Set of Organisms
   * @return a Set
   */
  public Set getOrganisms();

  /**
   *
   * @return an iterator over the set of Organisms.
   */
  public Iterator organisms();

  /**
   * Gets the Size of the population
   * @return the size
   */
  public int size();


  public static ChangeType ORGANISMS =
      new ChangeType("Organisms changed",Population.class,"ORGANISMS");

  public static ChangeType NAME =
      new ChangeType("Name changed",Population.class,"NAME");

}