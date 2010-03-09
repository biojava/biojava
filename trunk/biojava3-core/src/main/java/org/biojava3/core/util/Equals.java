package org.biojava3.core.util;

/**
 * A set of helper methods which return true if the two parameters are
 * equal to each other.
 *
 * @author ayates
 */
public class Equals {

  public static boolean equal(int one, int two) {
    return (one == two);
  }

  public static boolean equal(long one, long two) {
    return (one == two);
  }

  public static boolean equal(boolean one, boolean two) {
    return (one == two);
  }

  public static boolean equal(Object one, Object two) {
    boolean equal = false;
    //Both are null means they're equal
    if (one == null && two == null) {
      equal = true;
    }
    else if (one == null || two == null) {
      equal = false;
    }
    //True only if they are the same object
    else if (one == two) {
      equal = true;
    }
    //Otherwise ask the Object on their equals method.
    else {
      equal = one.equals(two);
    }
    return equal;
  }

  /**
   * This method should be called before beginning any equals methods. In order
   * to return true the method:
   *
   * <ol>
   * <li>The two given objects are the same instance using ==. This also means
   * if both Objects are null then this method will return true (well
   * technically they are equal)</li>
   * <li>Tests that neither object is null</li>
   * <li>The the two classes from the objects are equal using ==</li>
   * </ol>
   *
   * The boilerplate using this method then becomes:
   *
   * <pre>
   * boolean equals = false;
   * if (EqualsHelper.classEqual(this, obj)) {
   *   TargetClass casted = (TargetClass) obj;
   *   equals = (EqualsHelper.equal(this.getId(), casted.getId()) &amp;&amp; EqualsHelper
   *       .equal(this.getName(), casted.getName()));
   * }
   * return equals;
   * </pre>
   *
   * @param one
   *          The first object to test
   * @param two
   *          The second object to test
   * @return A boolean indicating if the logic agrees that these two objects are
   *         equal at the class level
   */
  public static boolean classEqual(Object one, Object two) {

    boolean equal = false;

    if (one == two) {
      equal = true;
    }
    else if (one == null || two == null) {
      equal = false;
    }
    else {
      equal = (one.getClass() == two.getClass());
    }

    return equal;
  }
}
