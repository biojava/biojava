package org.biojava3.core.util;

import java.lang.reflect.Array;

/**
 * Contains helper methods for generating a HashCode without having to resort to
 * the commons lang hashcode builders.
 *
 * Example where the property name is a String and the property age is an int
 *
 * <pre>
 * public int hashCode() {
 *   int result = Hashcoder.SEED;
 *   result = Hashcoder.hash(result, this.getName());
 *   result = Hashcoder.hash(result, this.getAge());
 *   return result;
 * }
 * </pre>
 *
 * @author ayates
 */
public class Hashcoder {

  /**
   * An initial value for a <code>hashCode</code>, to which we add
   * contributions from fields. Using a non-zero value decreases collisions of
   * <code>hashCode</code> values.
   */
  public static final int SEED = 9;

  /**
   * The prime number used to multiply any calculated hashcode seed by
   *
   * i.e. result = PRIME*result + c
   *
   * Where result is the result of the previous calculation (at first this will
   * be seed) and c is the calculated int to add to result
   */
  public static final int PRIME = 79;

  public static int hash(int seed, boolean b) {
    return (PRIME * seed) + (b ? 1 : 0);
  }

  public static int hash(int seed, char c) {
    return (PRIME * seed) + (int) c;
  }

  /**
   * Used for ints, bytes and shorts
   */
  public static int hash(int seed, int i) {
    return (PRIME * seed) + i;
  }

  /**
   * long support done by shifting by 32 (using unsigned shift)
   */
  public static int hash(int seed, long l) {
    return (PRIME * seed) + (int) (l ^ (l >>> 32));
  }

  /**
   * float support done via {@link Float#floatToIntBits(float)} which allows
   * the encoding of a float as an int. We can then hash as normal.
   */
  public static int hash(int seed, float f) {
    return hash(seed, Float.floatToIntBits(f));
  }

  /**
   * double support which is done using the same techinque as float hashing
   * except we convert to a long not to an int.
   */
  public static int hash(int seed, double d) {
    return hash(seed, Double.doubleToLongBits(d));
  }

  /**
   * <code>o</code> is a possibly-null object field, and possibly an
   * array.
   *
   * If <code>o</code> is an array, then each element may be a primitive
   * or a possibly-null object.
   */
  public static int hash(int seed, Object o) {
    int result = seed;
    //If it was null then this is the result 0
    if (o == null) {
      result = hash(result, 0);
    }
    //If it wasn't an array then calculate the hashcode
    else if (!o.getClass().isArray()) {
      result = hash(result, o.hashCode());
    }
    //Otherwise loop &
    else {
      int length = Array.getLength(o);
      for (int i = 0; i < length; i++) {
        Object item = Array.get(o, i);
        result = hash(result, item);// recursive call!
      }
    }
    return result;
  }
}