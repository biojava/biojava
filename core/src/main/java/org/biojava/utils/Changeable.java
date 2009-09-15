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
 */

package org.biojava.utils;

/**
 * This is a flag interface that defines the common add/remove listener methods
 * for classes and interfaces that wish to indicate that they are sources of
 * ChangeEvents.
 *
 * @author Matthew Pocock
 */
public interface Changeable {
  /**
   * Add a listener that will be informed of all changes.
   *
   * @param cl  the ChangeListener to add
   * @deprecated  use addChangeListener(cl, ChangeType.UNKNOWN)
   */
  public void addChangeListener(ChangeListener cl);
  
  /**
   * Add a listener that will be informed of changes of a given type.
   *
   * @param cl  the ChangeListener
   * @param ct  the ChangeType it is to be informed of
   */
  public void addChangeListener(ChangeListener cl, ChangeType ct);
  
  /**
   * Remove a listener that was interested in all types of changes.
   *
   * @param cl  a ChangeListener to remove
   * @deprecated  use removeChangeListener(cl, ChangeType.UNKNOWN)
   */
  public void removeChangeListener(ChangeListener cl);

  /**
   * Remove a listener that was interested in a specific types of changes.
   *
   * @param cl  a ChangeListener to remove
   * @param ct  the ChangeType that it was interested in
   */
  public void removeChangeListener(ChangeListener cl, ChangeType ct);
  
  /**
   * <p>
   * A particular ChangeType can never be raised by this Changeable.
   * </p>
   *
   * <p>
   * If this returns true, then it is guaranteed that change events of this type
   * (and all child types) can never under any circumstances be fired by this
   * Changeable instance. If it returns false, that does not mean that this type
   * of event will or even can be raised, but that it is worth registering
   * listeners incase.
   * </p>
   *
   * @param ct  the ChangeType to check
   * @return true if ChangeEvents of this type are guaranteed to never be fired
   */
  public boolean isUnchanging(ChangeType ct);
}
