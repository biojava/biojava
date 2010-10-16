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

import java.util.Collections;
import java.util.Set;

/**
 * A no-frills implementatation of Taxon.
 *
 * <p>A TaxonFactory implementation will probably wish to sub-class
 * this and add package-private accessors for the parent and children
 * fields as well as a pacakge-private constructor.</p>
 *
 * @author Matthew Pocock
 * @deprecated replaced by classes in {@link org.biojavax.bio.taxa org.biojavax.bio.taxa}
 */
public class SimpleTaxon extends AbstractTaxon {
  protected Taxon parent;
  protected Set children;
  
  protected SimpleTaxon() { super(); }
  
  /**
   * Create a new instance with no parent, no children and given
   * scientific and common names.
   */
  protected SimpleTaxon(String scientificName, String commonName) {
    super(scientificName, commonName);
  }
  
  public Taxon getParent() {
    return parent;
  }
  
  void setParent(Taxon parent) {
    this.parent = parent;
  }
  
  public Set getChildren() {
    if(children != null) {
      return children;
    } else {
      return Collections.EMPTY_SET;
    }
  }
}
