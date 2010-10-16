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

import java.util.Set;

import org.biojava.bio.Annotatable;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;

/**
 * <p>A taxon within a classification.</p>
 *
 * <p>Taxa may be 'leaf' nodes specifying species, or 'internal' nodes
 * specifying kingdoms and the like.</p>
 *
 * @author Matthew Pocock
 * @deprecated replaced by classes in {@link org.biojavax.bio.taxa org.biojavax.bio.taxa}
 */
public interface Taxon extends Annotatable {
  /**
   * Change type to indicate that the common name of this Taxon is
   * changing.
   */
  public final static ChangeType CHANGE_COMMON_NAME = new ChangeType(
    "Common name change",
    Taxon.class,
    "CHANGE_COMMON_NAME"
  );
  
  /**
   * Change type to indicate that the scientific name of this Taxon is
   * changing.
   */
  public final static ChangeType CHANGE_SCIENTIFIC_NAME = new ChangeType(
    "Scientific name change",
    Taxon.class,
    "CHANGE_SCIENTIFIC_NAME"
  );

  /**
   * <p>The common name of the Taxon.</p>
   *
   * <p>This is the normal name used in common speech, such as
   * 'human'.</p>
   *
   * @return a String representing this Taxon's common name
   */
  public String getCommonName();
  
  /**
   * <p>Set the new common name of this Taxon.</p>
   *
   * @param commonName  the new common name
   * @throws ChangeVetoException if the name can't be changed at this time
   */
  public void setCommonName(String commonName)
  throws ChangeVetoException;
  
  /**
   * <p>The scientific name of this taxon.</p>
   *
   * <p>This will be the portion of the scientific classification
   * pertaining to just this node within the classifictaion. It will
   * be something like 'homo sapiens' or 'archaeal group 2', rather
   * than the full classification list.</p>
   */
  public String getScientificName();

  /**
   * Change the scientific name of this species.
   *
   * @param scientificName the new scientific name
   * @throws ChangeVetoException if the scientific name can't be
   * changed at this time
   */
  public void setScientificName(String scientificName)
  throws ChangeVetoException;
  
  /**
   * <p>The parent of this Taxon.</p>
   *
   * <p>Taxa live within a tree data-structure, so every taxon has a
   * single parent except for the root type. This has the null
   * parent.</p>
   *
   * @return the parent Taxon, or null if this is the root type.
   */
  public Taxon getParent();
  
  /**
   *  <p>The children of this Taxon.</p>
   *
   * <p>Taxa live within a tree data-structure, so every taxon has
   * zero or more children. In the case of zero children, the empty
   * set is returned.</p>
   *
   * <p>? read-only ? dynamicaly updated with taxon object ? copy of
   * data ?</p>
   *
   * @return the Set (possibly empty) of all child Taxa
   */
  public Set getChildren();
  
  /**
   * <p>Two taxa are equal if they have equivalent children, common
   * and scientific names.</p>
   *
   * <p>Two different implementations of Taxon should be able to
   * appropriately trans-class equality. The parent of a Taxon is not
   * considered in testing equality as this potentially leads to
   * combinatorial problems checking whole taxa hierachies against one
   * another.</p>
   *
   * @param o  the object to check
   * @return true if o is a Taxon instance and has the same properties
   * as this
   */
  public boolean equals(Object o);
  
  /**
   * The hash-code of a Taxon is equal to the hash-code of it's
   * scientific name.
   */
  public int hashCode();
}
