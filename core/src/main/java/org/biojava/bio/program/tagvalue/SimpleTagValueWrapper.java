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

import org.biojava.utils.ParserException;

/**
 * <p>
 * Helper class to wrap one TagValueListener inside another one.
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
 * @author David Huen (change of TagValueWrapper to interface)
 * @since 1.2
 */
public abstract class SimpleTagValueWrapper
  implements
    TagValueWrapper
{
  private TagValueListener delegate;
  
  /**
   * Build a SimpleTagValueWrapper that will forward everything to a delegate.
   *
   * @param delegate the SimpleTagValueWrapper to forward events to
   */
  public SimpleTagValueWrapper(TagValueListener delegate) {
    this.delegate = delegate;
  }

  public SimpleTagValueWrapper() {
    delegate = null;
  }
  
  public TagValueListener getDelegate() {
    return delegate;
  }

  public void setDelegate(TagValueListener delegate) {
    this.delegate = delegate;
  }
  
  public void startRecord()
  throws ParserException {
    delegate.startRecord();
  }
  
  public void endRecord()
  throws ParserException {
    delegate.endRecord();
  }
  
  public void startTag(Object tag)
  throws ParserException {
    delegate.startTag(tag);
  }
  
  public void endTag()
  throws ParserException {
    delegate.endTag();
  }
  
  public void value(TagValueContext ctxt, Object value)
  throws ParserException {
    delegate.value(ctxt, value);
  }
}

