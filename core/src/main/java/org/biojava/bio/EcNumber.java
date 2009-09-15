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
package org.biojava.bio;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * An ec (enzyme classification) number.
 *
 *
 * Implementations of this interface should be imutable. This makes them much
 * more usefull as keys in maps.
 * 
 * it is a good idea to validate that the data being passed in is a sane ec
 * number.
 *
 * @author Matthew Pocock
 * @since 1.4
 */
public interface EcNumber {
  /**
   * A Pattern that can be used to parse EC strings into the indiidual numbers.
   */
  public static final Pattern EC_PATTERN =
    Pattern.compile("((\\d*)|(-))\\.((\\d*)|(-))\\.((\\d*)|(-))\\.((\\d*)|(-))");

  /**
   * Constant that represents EC number components that are not defined. This
   * is often represented as a '-' in EC strings.
   */
  public static final int UNDEFINED = -1;

  /**
   * Constant that represents EC number components that are as yet unclassified.
   * This is often represented as 99 in EC strings.
   */
  public static final int UNCLASSIFIED = 99;

  /**
   * Get the class number associated with the particular level of the ec number.
   *
   * <p>The index can be between 0 and 3 inclusive. 0 correxpons to the top
   * level class, 1 to the sub-class and so on. A return value of UNDEFINED
   * indicates that this field is not populated.</p>
   *
   * @param level  the level in the ec classification to return the number for
   * @return the value at that level
   */
  public int getClassNumber(int level);

  /**
   * A simple implementation of EcNumber.
   *
   * @author Matthew Pocock
   * @since 1.4
   */
  public static class Impl
  implements EcNumber {
    private int[] classes;

    /**
     * Make a new EcNumber.Impl with the data provided.
     *
     * @param mainClass     the main class number
     * @param subClass      the sub class number
     * @param subSubClass   the sub-sub class number
     * @param group         the group number
     */
    public Impl(int mainClass, int subClass, int subSubClass, int group) {
      this.classes = new int[] { mainClass, subClass, subSubClass, group };
    }

    private Impl(String ecString) {
      Matcher matcher = EC_PATTERN.matcher(ecString);
      if(!matcher.matches()) {
        throw new IllegalArgumentException(
          "Can't parse ec string: " + ecString );
      }

      classes = new int[] {
        process(matcher.group(1)),
        process(matcher.group(4)),
        process(matcher.group(7)),
        process(matcher.group(10))
      };
    }

    private int process(String s) {
      if(s.length() > 0) {
        if(s.equals("-")) {
          return UNDEFINED;
        } else {
          return Integer.parseInt(s);
        }
      } else {
        return UNDEFINED;
      }
    }

    public int getClassNumber(int level) {
      return classes[level];
    }

    public String toString() {
      StringBuffer sBuf = new StringBuffer();
      sBuf.append(process(getClassNumber(0)));
      for(int i = 1; i < 4; i++) {
        sBuf.append(".");
        sBuf.append(process(getClassNumber(i)));
      }
      return sBuf.toString();
    }

    private String process(int val) {
      if(val == UNDEFINED) {
        return "";
      } else {
        return Integer.toString(val);
      }
    }

    public boolean equals(Object obj) {
      if(obj instanceof EcNumber) {
        EcNumber that = (EcNumber) obj;

        for(int i = 0; i < 4; i++) {
          if(this.getClassNumber(i) != that.getClassNumber(i)) {
            return false;
          }
        }

        return true;
      }

      return false;
    }

    public int hashCode() {
      return
        getClassNumber(0) * 1000000 +
        getClassNumber(1) * 10000 +
        getClassNumber(2) * 100 +
        getClassNumber(3);
    }

    /**
     * Process a string into an EcNumber.
     *
     * <p>
     * This method uses the {@link EcNumber#EC_PATTERN} regular expression.
     * </p>
     *
     * @param ecString  String to parse
     * @return a new EcNumber
     * @throws IllegalArgumentException  if ecString could not be parsed
     */
    public static Impl valueOf(String ecString) {
      return new Impl(ecString);
    }
  }
}
