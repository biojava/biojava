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
 * Rename tags using a TagMapper.
 * </p>
 *
 * <p> 
 * This will rename tags as they stream into this listener using a TagMapper.
 * Once renamed, the events will be forwarded onto a delegate TagValueListener
 * for further processing.
 * </p>
 *
 * @author Matthew Pocock
 * @since 1.2
 */
public class TagRenamer extends SimpleTagValueWrapper {
  private PropertyChanger mapper;
  
  /**
   * Build a new TagRenamer with a delegate and mapper.
   *
   * @param delegate TagValueListener to pass mapped events onto
   * @param mapper TagMapper used to rename tags
   */
  public TagRenamer(TagValueListener delegate, PropertyChanger mapper) {
    super(delegate);
    this.mapper = mapper;
  }
  
  /**
   * Retrieve the mapper used to rename tags
   *
   * @return the current mapper
   */
  public PropertyChanger getMapper() {
    return mapper;
  }
  
  public void startTag(Object tag)
  throws ParserException {
    Object newTag = mapper.getNewTag(tag);
    if(newTag != null) {
      tag = newTag;
    }
    
    super.startTag(tag);
  }
}

