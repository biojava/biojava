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

import java.lang.ref.WeakReference;
import java.util.Collections;
import java.util.Set;

/**
 * <p>An implementation of Taxon that keeps only weak references to
 * children, but full references to parents.</p>
 *
 * <p>This may be suitable for deriving memory-savy implementations
 * of TaxonFactory.</p>
 *
 * <p>To manipulate the children set, use the getChildrenRaw and
 * setChildrenRaw methods. These 'box' the actual weak reference, but
 * recognize null to mean that there are no children currently
 * known. A code-fragment may wish to do something like this:</p>
 *
 * <pre><code>
 * Set children = weakTaxon.getChildrenRaw();
 * if(children == null) {
 *   children = new HashSet();
 *   weakTaxon.setChildrenRaw(children);
 * }
 * // do stuff to update child set e.g. add children 
 * </code></pre>
 * </p>
 *
 * @author Matthew Pocock
 * @deprecated replaced by classes in {@link org.biojavax.bio.taxa org.biojavax.bio.taxa}
 */
public class WeakTaxon extends AbstractTaxon {
  protected Taxon parent;
  private WeakReference /*Set*/ children;
  
  public WeakTaxon() {
    super();
  }
  
  public WeakTaxon(String scientificName, String commonName) {
    super(scientificName, commonName);
  }
  
  public Taxon getParent() {
    return parent;
  }
  
  void setParent(Taxon parent) {
    this.parent = parent;
  }
  
  public Set getChildren() {
    Set c = getChildrenRaw();
    if(c != null) {
      return c;
    } else {
      return Collections.EMPTY_SET;
    }
  }
  
  public Set getChildrenRaw() {
    if(children != null) {
      Set c = (Set) children.get();
      if(c != null) {
        return c;
      }
    }
    
    return null;
  }
  
  public void setChildrenRaw(Set children) {
    this.children = new WeakReference(children);
  }
}
