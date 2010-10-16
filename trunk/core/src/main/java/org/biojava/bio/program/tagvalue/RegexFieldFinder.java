package org.biojava.bio.program.tagvalue;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.biojava.utils.ParserException;

public class RegexFieldFinder
extends SimpleTagValueWrapper {
  private final Pattern pattern;
  private final String[] tags;
  private final boolean inLine;

  /**
   * Creates a new RegexFiledFinder.
   *
   * @param delegate  the TagValueListener to forward events to
   * @param pattern a Pattern to match to values
   * @param tags an array of Strings giving tag names for each group in the
   *        pattern
   * @param inLine if false, an entire sub-document will be generated for
   *        the parent tag
   */
  public RegexFieldFinder(
    TagValueListener delegate,
    Pattern pattern,
    String[] tags,
    boolean inLine
  ) {
    super(delegate);
    this.pattern = pattern;
    this.tags = tags;
    this.inLine = inLine;
  }

  public void startTag(Object tag)
  throws ParserException {
    if(!inLine) {
      super.startTag(tag);
      super.startRecord();
    }
  }

  public void endTag()
  throws ParserException {
    if(!inLine) {
      super.endRecord();
      super.endTag();
    }
  }

  public void value(TagValueContext ctxt, Object val)
  throws ParserException {
    try {
      Matcher m = pattern.matcher(val.toString());
      m.find();
      
      for(int i = 0; i < tags.length; i++) {
        super.startTag(tags[i]);
        super.value(ctxt, m.group(i + 1));
        super.endTag();
      }
    } catch (IllegalStateException ise) {
      throw new ParserException("Problem matching " + pattern.pattern() + " to " + val);
    }
  }
}   
