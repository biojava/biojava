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

package org.biojava.bio;

import java.util.Collection;

/**
 * Hidden implementation detail. Damn those interfaces not having private
 * inner classes & static fields.
 *
 * @author Matthew Pocock
 * @author Thomas Down
 * @since 1.4
 */
class NoneCollectionConstraint implements CollectionConstraint {
  public boolean accept(Object value) {
    return false;
  }
  
  public boolean subConstraintOf(CollectionConstraint subConstraint) {
    return subConstraint instanceof NoneCollectionConstraint;
  }
  
  public boolean validateAddValue(Collection oldcoll, Object newval)
  {
      return false;
  }
  
  public boolean validateRemoveValue(Collection oldcoll, Object victime)
  {
      return false;
  }
  
  public String toString() {
    return "NONE";
  }
}
