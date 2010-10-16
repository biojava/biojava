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

import java.util.List;
import java.util.Map;

import org.biojava.utils.ParserException;
import org.biojava.utils.SmallMap;

/**
 * <p>
 * A mapping between keys and actions to turn old values into new values.
 * </p>
 *
 * @author Matthew Pocock
 * @since 1.3
 */
public class ChangeTable {
  public static final Changer STRING_TO_INT = new Changer() {
    public Object change(Object val) {
      return new Integer((String) val);
    }
  };

  private final Map changers;
  private final Map splitters;

  public ChangeTable() {
    this.changers = new SmallMap();
    this.splitters = new SmallMap();
  }
  

  /**
   * Set the Changer to be used for all values of a particular tag.
   *
   * @param tag the tag Object which will have all values changed
   * @param changer the Changer used to change the values
   */
  public void setChanger(Object tag, Changer changer) {
    changers.put(tag, changer);
  }
  
  /**
   * Set the Splitter to be used for all values of a particular tag.
   *
   * @param tag the tag Object which will have all values split
   * @param splitter the Splitter used to split the values
   */
  public void setSplitter(Object tag, Splitter splitter) {
    splitters.put(tag, splitter);
  }
  
  /**
   * Get the Changer currently registered to handle a tag.
   *
   * @param tag  the tag Object for which values would be changed
   * @return the associated Changer or null
   */
  public Changer getChanger(Object tag) {
    return (Changer) changers.get(tag);
  }
  
  /**
   * Get the Splitter currently registered to handle a tag.
   *
   * @param tag  the tag Object for which values would be split
   * @return the associated Splitter or null
   */
  public Splitter getSplitter(Object tag) {
    return (Splitter) splitters.get(tag);
  }
  
  public Object change(Object tag, Object value)
  throws ParserException {
    Changer c = (Changer) changers.get(tag);
    if(c != null) {
      return c.change(value);
    }
    
    Splitter s = (Splitter) splitters.get(tag);
    if(s != null) {
      return s.split(value);
    }
    
    return value;
  }
  
  /**
   * Callback used to produce a new value from an old one.
   *
   * @author Matthew Pocock
   * @since 1.3
   */
  public static interface Changer {
    /**
     * <p>
     * Produce a modified value from an old value.
     * </p>
     *
     * <p>
     * It is strongly recommended that this method is re-entrant and does not
     * modify the state of the Changer in a way that would affect future return
     * -values.
     * </p>
     *
     * @param value  the old value Object
     * @return  the new value Object
     * @throws ParserException if value could not be changed
     */
    public Object change(Object value)
    throws ParserException;
  }
  
  /**
   * Callback used to produce a list of values from a single old one.
   *
   * @author Matthew Pocock
   * @since 1.3
   */
  public static interface Splitter {
    /**
     * <p>
     * Produce a list of values from an old value.
     * </p>
     *
     * <p>
     * It is strongly recommended that this method is re-entrant and does not
     * modify the state of the Splitter in a way that would affect future return
     * -values.
     * </p>
     *
     * @param value  the old value Object
     * @return  a List of value Objects produced by splitting the old value
     *          Object
     * @throws ParserException if the value could not be split
     */
    public List split(Object value)
    throws ParserException;
  }
  
  /**
   * An implementation of Changer that applies a list of Changer instances to
   * the value in turn.
   *
   * @author Matthew Pocock
   * @since 1.3
   */
  public static class ChainedChanger
  implements Changer {
    private Changer[] changers;
    
    public ChainedChanger(Changer[] changers) {
      this.changers = new Changer[changers.length];
      
      System.arraycopy(changers, 0, this.changers, 0, changers.length);
    }
    
    public Object change(Object value)
    throws ParserException {
      for(int i = 0; i < changers.length; i++) {
        value = changers[i].change(value);
      }
      
      return value;
    }
  }
}
