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

package org.biojava.bio.program.tagvalue;


/**
 * <p>
 * Interface for TagValueListeners that wrap other TagValueListeners
 * </p>
 *
 * <p>
 * Implementations will tend to intercept the tags or values as they stream
 * through and modify them in some manner before forwarding them to the delegate
 * listener. Using classes derived from SimpleTagValueWrapper, it is possible to build
 * up complex chains of handlers that process and collate information as it
 * streams through.
 * </p>
 *
 * @author Matthew Pocock
 * @author David Huen (conversion to interface)
 * @since 1.2
 */
public interface TagValueWrapper
  extends
    TagValueListener
{
  /**
   * get listener to which all calls will be delegated
   */  
  public TagValueListener getDelegate();

  /**
   * set listener to which all calls will be delegated
   */
  public void setDelegate(TagValueListener delegate);
}

