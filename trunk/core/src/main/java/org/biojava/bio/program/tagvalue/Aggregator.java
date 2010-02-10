package org.biojava.bio.program.tagvalue;

import org.biojava.utils.ParserException;

/**
 * Joins multipel values into single values.
 *
 * <p>
 * Some properties have values spread across multiple lines. For example,
 * the properties on EMBL features can be spread across multiple lines.</p>
 *
 * <p>
 * This class provides callbacks to allow event streams to be re-written
 * so that they contain this information.
 * </p>
 *
 * @since 1.4
 * @author Matthew Pocock
 */
public class Aggregator extends SimpleTagValueWrapper {
  private BoundaryFinder observer;
  private String joiner;

  // state
  //
  private StringBuffer value;
  private boolean inValue;

  public Aggregator(TagValueListener listener, BoundaryFinder observer, String joiner) {
    super(listener);
    this.observer = observer;
    this.joiner = joiner;
    this.value = new StringBuffer();
  }

  public BoundaryFinder getBoundaryFinder() {
    return observer;
  }
  
  public void setBoundaryFinder(BoundaryFinder finder) {
    this.observer = finder;
  }
  
  public String getJoiner() {
    return joiner;
  }
  
  public void setJoiner(String joiner) {
    this.joiner = joiner;
  }
  
  public void startTag(Object tag)
  throws ParserException {
    super.startTag(tag);
    inValue = false;
    value.setLength(0);
  }
  
  public void value(TagValueContext ctxt, Object value)
  throws ParserException {
    boolean isStart = observer.isBoundaryStart(value);
    boolean isEnd = observer.isBoundaryEnd(value);
    boolean dbv = observer.dropBoundaryValues();
    
    if(isStart && isEnd) {
      if(inValue) {
        super.value(ctxt, this.value.toString());
        this.value.setLength(0);
        inValue = false;
      }
      
      if(!dbv) {
        super.value(ctxt, value);
      }
    }
    
    if(isStart && !isEnd) {
      if(inValue) {
        super.value(ctxt, this.value.toString());
        this.value.setLength(0);
      }
      this.value.append(value);
      inValue = true;
    }
    
    if(!isStart && isEnd) {
      if(!dbv) {
        this.value.append(joiner);
        this.value.append(value);
      }
      
      super.value(ctxt, this.value.toString());
      inValue = false;
      this.value.setLength(0);
    }
    
    if(!isStart && !isEnd) {
      this.value.append(joiner);
      this.value.append(value);
      inValue = true;
    }
  }

  public void endTag()
  throws ParserException {
    if(inValue) {
      super.value(null, this.value.toString());
      inValue = false;
    }
    super.endTag();
  }
}
