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
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.WeakHashMap;

import org.biojava.utils.SmallSet;

/**
 * <p>An implementation of TaxonFactory that builds a weak in-memory
 * Taxon tree.</p>
 *
 * <p>This implementation holds only weak references to the Taxon
 * instances it knows about. This means that WeakTaxonFactory may not
 * be appropriate for situations where you wish to browse the taxon
 * tree. It does, however, mean that massive taxa can be represented,
 * by effectively reflecting the currently useful rooted sub-tree in
 * memory.</p>
 *
 * @author Matthew Pocock
 * @deprecated replaced by classes in {@link org.biojavax.bio.taxa org.biojavax.bio.taxa}
 */
public class WeakTaxonFactory implements TaxonFactory {
  /**
   * The TaxonFactory that the biojava system should use for storing
   * the taxonomy used by swissprot and embl as in-memory objects.
   */
  public static final WeakTaxonFactory GLOBAL
    = new WeakTaxonFactory("GLOBAL");
  
  private final Taxon root;
  private final String name;
  private final Map taxonBySciName = new WeakHashMap();
  
  public WeakTaxonFactory(String name) {
    this.name = name;
    this.root = createTaxon("ROOT", "");
  }
  
  public Taxon getRoot() {
    return root;
  }
  
  public String getName() {
    return name;
  }
  
  public Taxon importTaxon(Taxon taxon) {
    WeakTaxon can = canonicalize(taxon);
    if(can == null) {
      can = new WeakTaxon(taxon.getScientificName(), taxon.getCommonName());
      
      for(Iterator i = taxon.getChildren().iterator(); i.hasNext(); ) {
        Taxon child = (Taxon) i.next();
        addChild(can, child);
      }
      
      return can;
    } else {
      return can;
    }
  }
  
  public Taxon createTaxon(String scientificName, String commonName) {
    Taxon taxon = new WeakTaxon(scientificName, commonName);
    taxonBySciName.put(scientificName, new WeakReference(taxon));
    return taxon;
  }
  
  public Taxon addChild(Taxon parent, Taxon child) {
    WeakTaxon sparent = (WeakTaxon) importTaxon(parent);
    WeakTaxon schild = (WeakTaxon) importTaxon(child);
    
    Set children = sparent.getChildrenRaw();
    if(children == null) {
      children = new SmallSet();
      sparent.setChildrenRaw(children);
    }
    
    children.add(schild);
    schild.setParent(sparent);
    
    return schild;
  }
  
  public Taxon removeChild(Taxon parent, Taxon child) {
    WeakTaxon sparent = canonicalize(parent);
    WeakTaxon schild = canonicalize(child);
    
    if(sparent == null) {
      throw new IllegalArgumentException("Don't know about parent taxon");
    }
    
    Set children = sparent.getChildrenRaw();
    if(
      (schild != null) &&
      (children != null) &&
      (children.remove(schild))
    ) {
      return schild;
    } else {
      return null;
    }
  }
  
  public Taxon search(Object id) {
    WeakReference wr = (WeakReference) taxonBySciName.get(id);
    if(wr != null) {
      return (Taxon) wr.get();
    } else {
      return null;
    }
  }
  
  private WeakTaxon canonicalize(Taxon taxon) {
    return (WeakTaxon) search(taxon.getScientificName());
  }
}
