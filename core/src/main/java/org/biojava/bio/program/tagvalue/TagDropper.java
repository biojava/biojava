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

import java.util.Set;

import org.biojava.utils.ParserException;
import org.biojava.utils.SmallSet;

/**
 * Silently drop all tags except those specified, and pass the rest onto a
 * delegate.
 *
 * @author Matthew Pocock
 * @since 1.2
 */
public class TagDropper
extends SimpleTagValueWrapper {
  private final Set tags;
  private boolean retain;
  private boolean propagate;

  {
    tags = new SmallSet();
    retain = true;
  }

  public TagDropper() {
    super();
  }
  
  /**
   * Create a new TagDropper that will pass on all retained tags and values to
   * tvl. Initially, no tags will be retained.
   *
   * @param tvl  the TagValueListener to inform of all surviving events
   */
  public TagDropper(TagValueListener tvl) {
    super(tvl);
  }
  
  /**
   * Add a tag to retain.
   *
   * @param tag  a tag that will be forwarded to the delegate
   */
  public void addTag(Object tag) {
    tags.add(tag);
  }
  
  /**
   * Remove a tag so that it will not be retained.
   *
   * @param tag  a tag that will not be forwarded to the delegate
   */
  public void removeTag(Object tag) {
    tags.remove(tag);
  }

  /**
   * Get the complete set of tags that are currently recognized.
   *
   * @return the Set of known tags
   */
  public Set getTags() {
    return tags;
  }
  
  /**
   * Set wether known tags are to be retained or dropped.
   *
   * <p>
   * If retain is true, then all known tags will be passed on and all other tags
   * will be discarded. If retain is true, then all known tags will be dropped
   * and all others passed on.
   * </p>
   *
   * @param retain  true if the tags are to be retained
   */
  public void setRetain(boolean retain) {
    this.retain = retain;
  }
  
  /**
   * Find out if known tags are retained or dropped.
   *
   * @return true if values are retained
   */
  public boolean doRetain() {
    return retain;
  }
  
  public void startTag(Object tag)
  throws ParserException {
    propagate = retain == tags.contains(tag);
    
    if(propagate) {
      super.startTag(tag);
    }
  }
  
  public void endTag()
  throws ParserException {
    if(propagate) {
      super.endTag();
    }
  }
  
  public void value(TagValueContext ctxt, Object value)
  throws ParserException {
    if(propagate) {
      super.value(ctxt, value);
    }
  }
}

