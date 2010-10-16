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

import java.util.Map;
import java.util.Set;

import org.biojava.utils.ParserException;
import org.biojava.utils.SmallMap;

/**
 * <p>
 * Pushes a new parser and listener, or delegate to a listener depending on the
 * tag.
 * </p>
 *
 * <p>
 * setParserListener() is used to associate a tag with a TagValueParser and
 * TagValueListener. When this tag is encountered, the pair will be pushed onto
 * the parser processing stack and will gain control of the stream until that
 * tag has ended. setListener() is used to associate a listener with a tag that
 * will be used to handle those values without pushing a sub-context.
 * The delegator is constructed with a default TagValueListener that will be
 * informed of all events for which there are no explicit delegate pairs
 * registered.
 * </p>
 *
 * @author Matthew Pocock
 * @since 1.2
 */
public class TagDelegator
  extends
    SimpleTagValueWrapper
{
  private Map parsers;
  private Map listeners;
  private TagValueParser delegateParser;
  
  private TagValueParser parser;
  private TagValueListener listener;

  {
    parsers = new SmallMap();
    listeners = new SmallMap();
  }

  public TagDelegator() {
    super();
  }
  
  public TagDelegator(TagValueListener delegate) {
    super(delegate);
    parsers = new SmallMap();
    listeners = new SmallMap();
  }

  public void setDelegateParser(TagValueParser delegateParser) {
    this.delegateParser = delegateParser;
  }

  public TagValueParser getDelegateParser() {
    return delegateParser;
  }
  
  public void startTag(Object tag)
  throws ParserException {
    parser = (TagValueParser) parsers.get(tag);
    listener = (TagValueListener) listeners.get(tag);

    if(parser == null && listener != null) {
      listener.startTag(tag);
    } else {
      super.startTag(tag);
    }
  }

  public void endTag()
  throws ParserException {
    if(parser == null && listener != null) {
      listener.endTag();
    } else {
      super.endTag();
    }
  }

  public void value(TagValueContext tvc, Object value)
  throws ParserException {
    if(parser != null) {
      tvc.pushParser(parser, listener);
    } else if(listener != null) {
      listener.value(tvc, value);
    } else if(delegateParser != null) {
      tvc.pushParser(delegateParser, getDelegate());
    } else {
      super.value(tvc, value);
    }
  }
  
  public void setParserListener(
    Object tag,
    TagValueParser parser,
    TagValueListener listener
  ) {
    parsers.put(tag, parser);
    listeners.put(tag, listener);
  }
  
  public void setListener(
    Object tag,
    TagValueListener listener
  ) {
    listeners.put(tag, listener);
  }

  public TagValueParser getParser(Object tag) {
    return (TagValueParser) parsers.get(tag);
  }

  public TagValueListener getListener(Object tag) {
    return (TagValueListener) listeners.get(tag);
  }

  public Set getTags() {
    return listeners.keySet();
  }
}
