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
package org.biojava.bio.taxa;

import org.biojava.utils.ChangeVetoException;

/**
 * Factory for handling a particular implementation of a Taxon.
 *
 * @author Matthew Pocock
 * @deprecated replaced by classes in {@link org.biojavax.bio.taxa org.biojavax.bio.taxa}
 */

public interface TaxonFactory {
  /**
   * <p>Name for this TaxonFactory.</p>
   *
   * @return the name of this TaxonFactory
   */
  public String getName();
  
  /**
   * <p>Import a Taxon and all its children into the implementation
   * provided by this factory.</p>
   *
   * <p>The return value of this method should be .equals() and
   * .hasCode() compatable with the taxon parameter. It may not be
   * implemented by the same underlying implementation.</p>
   *
   * @param source the Taxon to copy
   * @return a new Taxon
   */
  public Taxon importTaxon(Taxon source);
  
  /**
   * <p>Retrieve the root upon which all rooted Taxon that this
   * factory knows about are rooted.</p>
   *
   * @return the 'root' Taxon
   */
  public Taxon getRoot();
  
  /**
   * <p>Retrieve a Taxon that matches some ID.</p>
   *
   * <p>This method is here out of desperation. It's nasty and should
   * be replaced by some propper querying API. Without having
   * different methods for every TaxonFactory I don't know what to
   * do. All ideas appreciated.</p>
   *
   * @param id  the Object identifying a Taxon
   * @return the Taxon matching the ID, or null if none match
   */
  public Taxon search(Object id);
  
  /**
   * <p>Create a new orphan Taxon with a given scientific and common
   * name.</p>
   *
   * @param scientificName  the scientificName to give the Taxon
   * @param commonName  the common name to give the Taxon
   * @return a new Taxon with no parent and no children
   */
  public Taxon createTaxon(String scientificName, String commonName);
  
  
  /**
   * <p>Add a taxon as a child to a parent.</p>
   *
   * <p>The TaxonFactory may chose to add the child directly, or make
   * a new object which is .equals() compatable with child. The actual
   * Taxon instance inserted into the child set is returned by the add
   * method.</p>
   *
   * @param parent the parent Taxon to add the child to
   * @param child  the Taxon to add as a child
   * @return the Taxon object actualy present as the child
   * @throws ChangeVetoException if for any reason the child can't be added
   * @throws CircularReferenceException if child is this Taxon or any
   * of its parents
   */
  public Taxon addChild(Taxon parent, Taxon child)
  throws ChangeVetoException,
  CircularReferenceException;
  
  /**
   * <p>Remove a Taxon as a child to this one.</p>
   *
   * <p>This Taxon should attempt to remove a child that is .equals()
   * compatable with child. If it is sucessful, it should return the
   * Taxon instance that was removed. If not, it should return
   * null.</p>
   *
   * @param parent the parent Taxon to remove the child from
   * @param child  the Taxon to remove as a child
   * @return the actual Taxon removed, or null if none were removed
   * @throws ChangeVetoException if for any reason the child can't be
   * removed
   */
  public Taxon removeChild(Taxon parent, Taxon child)
  throws ChangeVetoException;
}
