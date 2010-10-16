package org.biojava.bio.program.tagvalue;

import org.biojava.utils.ParserException;

/**
 * Partician multiple values for a tag into their own tag groups.
 *
 * <p>With tag-value files, it is not uncommon for information about logical
 * blocks of data to be encoded in the values as well as tags of the document.
 * For example, in swissprot entries, the comment block may be punctuated by
 * lines that clearly seperate one logical comment from another, but there may
 * be no change in the pattern of tags to indicate this. Some fields, such as
 * alternate names in Enzyme, are a ist of values but some of the values are
 * longer than a single line. Each value is terminated with a period ".", but
 * again there is no way from the tags to know the logical grouping.</p>
 *
 * <p>This class provides callbacks to allow event streams to be re-written
 * so that they contain this information. A single CC tag with multiple values
 * can be re-written as multiple CC tags with values for each logical comment.
 * This is done by presenting each value to an instance of Agregator.Observer
 * that indicates if the current value signals the end of a logcal block.
 * </p>
 *
 * @since 1.4
 * @author Matthew Pocock
 */
public class MultiTagger extends SimpleTagValueWrapper {
  private final BoundaryFinder observer;

  // state
  //
  boolean inTag;
  boolean seenValues;
  Object tag;

  public MultiTagger(TagValueListener listener, BoundaryFinder observer) {
    super(listener);
    this.observer = observer;
  }

  public BoundaryFinder getBoundaryFinder() {
    return observer;
  }
  
  public void startTag(Object tag)
  throws ParserException {
    this.tag = tag;
    inTag = false;
    seenValues = false;
  }

  public void value(TagValueContext ctxt, Object value)
  throws ParserException {
    seenValues = true;

    if(observer.isBoundaryStart(value)) {
      if(inTag) {
        super.endTag();
      }
      super.startTag(tag);
      if(!observer.dropBoundaryValues()) {
        super.value(ctxt, value);
      }
      inTag = true;
    } else if(observer.isBoundaryEnd(value)) {
      if(!inTag) {
        super.startTag(tag);
      }
      if(!observer.dropBoundaryValues()) {
        super.value(ctxt, value);
      }
      super.endTag();
      inTag = false;
    } else {
      if(!inTag) {
        super.startTag(tag);
        inTag = true;
      }
      super.value(ctxt, value);
    }
  }

  public void endTag()
  throws ParserException {
    if(inTag) {
      super.endTag();
    } else if(!seenValues) {
      // bounary condition where there are no values associated with a tag
      super.startTag(tag);
      super.endTag();
    }
    inTag = false;
    seenValues = false;
  }
}
