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
package org.biojava.bio.seq;
 
import java.io.NotSerializableException;
import java.io.ObjectStreamException;
import java.io.Serializable;

import org.biojava.bio.symbol.SymbolList;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.StaticMemberPlaceHolder;

/**
 * Adds the concept of 'strand' to features.
 * <p>
 * Strandedness only applies to some types of sequence, such as DNA. Any
 * implementation should blow chunks to avoid being added to a sequence for
 * which strand is a foreign concept.
 *
 * Strand is intrinsicly part of all {@link org.biojavax.bio.seq.RichFeature RichFeatures}
 * We strongly recommend using this interface.
 *
 * @author Matthew Pocock
 *@see org.biojavax.bio.seq.RichFeature
 */
public interface StrandedFeature extends Feature {

  /**
   * The strand of this feature is being altered.
   */
  public static final ChangeType STRAND =
      new ChangeType("Strand has altered", StrandedFeature.class, "STRAND");

  /**
   * Retrieve the strand that this feature lies upon.
   * <p>
   * This will be one of StrandedFeature.POSITIVE or NEGATIVE.
   *
   * @return one of the Strand constants
   */
  Strand getStrand();

  /**
   * Set the strand that this feature lies upon.
   * <p>
   * This will be one of StrandedFeature.POSITIVE or NEGATIVE.
   *
   * @param strand a <code>Strand</code>.
   *
   * @exception ChangeVetoException if the strand may not be
   * changed.
   */
  void setStrand(Strand strand) throws ChangeVetoException;

  /**
   * Return a list of symbols that are contained in this feature.
   * <p>
   * The symbols may not be contiguous in the original sequence, but they
   * will be concatenated together in the resulting SymbolList.
   * <p>
   * The order of the Symbols within the resulting symbol list will be 
   * according to the concept of ordering within the location object.
   * <p>
   * If the feature is on the negative strand then the SymbolList will be
   * reverse-complemented as appropriate.
   *
   * @return  a SymbolList containing each symbol of the parent sequence contained
   *          within this feature in the order they appear in the parent
   */
  SymbolList getSymbols();
  
  /**
   * Flag to indicate that a feature is on the positive strand.
   */
  static final Strand POSITIVE = new Strand("POSITIVE", +1, '+');

  /**
   * Flag to indicate that a feature is on the negative strand.
   */
  static final Strand NEGATIVE = new Strand("NEGATIVE", -1, '-');
  
  /**
   * Flag to indicate that a feature has an unknown strand.
   */
  static final Strand UNKNOWN = new Strand("UNKNOWN", 0, '.');
  
    /**
     * Template class for parameterizing the creation of a new
     * <code>StrandedFeature</code>.
     *
     * @author Matthew Pocock
     */

  public static class Template extends Feature.Template {
    public Strand strand;
  }
  
  /**
   * Class to represent the 'strandedness' of a feature.
   * <p>
   * Strandedness may be re-used in other situations, but basically what it means
   * is whether the feature has directionality, and if it does, does it travel
   * from its location min to max, or max to min.
   *
   * @author Matthew Pocock
   */
  public static final class Strand implements Serializable {
    private final String text;
    private final int value;
    private final char token;
    
    // Should be private. workaround for known javac 1.2 bug
    // http://developer.java.sun.com/developer/bugParade/bugs/4262105.html
    Strand(String text, int value, char token) {
      this.text = text;
      this.value = value;
      this.token = token;
    }
    public String toString() {
      return text;
    }
    
   /**
    * Returns the integer label for strandedness. That is, "+1", "-1",
    * or "0" for positive, negative, and unknown strands respectively.
    */
    public int getValue() {
      return value;
    }
    
   /**
    * Returns the token for strandedness. That is, "+","-","." for
    * positive, negative and unknown strands respectively.
    */
    public char getToken() {
      return token;
    }

    /**
     * Return a strand that represents flipping this onto the opposite strand.
     */
    public Strand flip() {
      if(this == POSITIVE) {
        return NEGATIVE;
      } else if(this == NEGATIVE) {
        return POSITIVE;
      } else {
        return this;
      }
    }

    private Object writeReplace() throws ObjectStreamException {
      try {
        return new StaticMemberPlaceHolder(StrandedFeature.class.getField(text));
      } catch (NoSuchFieldException nsfe) {
        throw new NotSerializableException(nsfe.getMessage());
      }
    }
  }
}
