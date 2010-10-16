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
 * Encapsulate the mapping between Taxon and stringified
 * representations of taxa.
 *
 * @author Matthew Pocock
 * @deprecated replaced by classes in {@link org.biojavax.bio.taxa org.biojavax.bio.taxa}
 */
public interface TaxonParser {
  /**
   * Convert a stringified Taxon into a Taxon instance.
   *
   * @param taxonFactory  the TaxonFactory used to instantiate taxa instances
   * @param taxonString  the String to parse
   * @return a Taxon instance created by the TaxonFactory from the taxonString
   */
  public Taxon parse(TaxonFactory taxonFactory, String taxonString)
  throws ChangeVetoException, CircularReferenceException;
  
  /**
   * Convert a Taxon into a stringified representation.
   *
   * @param taxon the Taxon to serialize
   * @return the stringified version of Taxon
   */
  public String serialize(Taxon taxon);
}
