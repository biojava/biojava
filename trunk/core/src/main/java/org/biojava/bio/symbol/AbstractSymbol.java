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


package org.biojava.bio.symbol;

import org.biojava.bio.Annotatable;
import org.biojava.bio.Annotation;
import org.biojava.utils.AbstractChangeable;
import org.biojava.utils.ChangeForwarder;
import org.biojava.utils.ChangeSupport;
import org.biojava.utils.ChangeType;

/**
 * The base-class for Symbol implementations.
 *
 * @author Matthew Pocock
 * @since 1.1
 */
public abstract class AbstractSymbol
  extends
    AbstractChangeable
  implements
    Symbol
{
  protected transient ChangeForwarder annotationForwarder = null;

  protected ChangeSupport getChangeSupport(ChangeType changeType) {
    ChangeSupport changeSupport = super.getChangeSupport(changeType);

    if(
      (Annotatable.ANNOTATION.isMatchingType(changeType)) &&
      (annotationForwarder == null)
    ) {
      annotationForwarder =
              new ChangeForwarder.Retyper(this, changeSupport, Annotatable.ANNOTATION);
      getAnnotation().addChangeListener(annotationForwarder, Annotation.PROPERTY);
    }

    return changeSupport;
  }

  public String toString() {
    return getClass().getName() + ": " + getName();
  }
}
