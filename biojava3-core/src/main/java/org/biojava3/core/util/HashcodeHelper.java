/* Original code licensed as follows under BSD:
 *
 * Copyright (c) 2002-2009, Hirondelle Systems
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *    * Redistributions of source code must retain the above copyright
 *      notice, this list of conditions and the following disclaimer.
 *    * Redistributions in binary form must reproduce the above copyright
 *      notice, this list of conditions and the following disclaimer in the
 *      documentation and/or other materials provided with the distribution.
 *    * Neither the name of Hirondelle Systems nor the
 *      names of its contributors may be used to endorse or promote products
 *      derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY HIRONDELLE SYSTEMS ''AS IS'' AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL HIRONDELLE SYSTEMS BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * Code re-licensed from the BioJava Project.
 */

package org.biojava3.core.util;

import java.lang.reflect.Array;

/**
 * Contains helper methods for generating a HashCode without having to resort to
 * the commons lang hashcode builders.
 *
 * The code was taken from http://www.javapractices.com/Topic28.cjp which
 * adhers to the rules in Effective Java pg 36 - Chapter 3 Item 8
 *
 * Example where the property name is a String and the property age is an int
 *
 * <pre>
 * public int hashCode() {
 *   int result = HashcodeHelper.SEED;
 *   result = HashcodeHelper.hash(result, this.getName());
 *   result = HashcodeHelper.hash(result, this.getAge());
 *   return result;
 * }
 * </pre>
 *
 * @author ayates
 */
public class HashcodeHelper {

  /**
   * An initial value for a <code>hashCode</code>, to which is added
   * contributions from fields. Using a non-zero value decreases collisons of
   * <code>hashCode</code> values.
   */
  public static final int SEED = 17;

  /**
   * The prime number used to multiply any calculated hascode seed by
   *
   * i.e. result = DEFAULT_PRIME_NUMBER*result + c
   *
   * Where result is the result of the previous calculation (at first this will
   * be seed) and c is the calculated int to add to result
   */
  public static final int DEFAULT_PRIME_NUMBER = 37;

  /**
   * booleans
   */
  public static int hash(int aSeed, boolean aBoolean) {
    return firstTerm(aSeed) + (aBoolean ? 1 : 0);
  }

  /**
   * chars
   */
  public static int hash(int aSeed, char aChar) {
    return firstTerm(aSeed) + (int) aChar;
  }

  /**
   * Used for ints, bytes and shorts via implicit conversion
   */
  public static int hash(int aSeed, int aInt) {
    return firstTerm(aSeed) + aInt;
  }

  /**
   * longs
   */
  public static int hash(int aSeed, long aLong) {
    return firstTerm(aSeed) + (int) (aLong ^ (aLong >>> 32));
  }

  /**
   * floats
   */
  public static int hash(int aSeed, float aFloat) {
    return hash(aSeed, Float.floatToIntBits(aFloat));
  }

  /**
   * doubles
   */
  public static int hash(int aSeed, double aDouble) {
    //This is handled via the seed(int,long) method
    return hash(aSeed, Double.doubleToLongBits(aDouble));
  }

  /**
   * <code>aObject</code> is a possibly-null object field, and possibly an
   * array.
   *
   * If <code>aObject</code> is an array, then each element may be a primitive
   * or a possibly-null object.
   */
  public static int hash(int aSeed, Object aObject) {
    int result = aSeed;
    if (aObject == null) {
      result = hash(result, 0);
    }
    else if (!isArray(aObject)) {
      result = hash(result, aObject.hashCode());
    }
    else {
      int length = Array.getLength(aObject);
      for (int idx = 0; idx < length; ++idx) {
        Object item = Array.get(aObject, idx);
        // recursive call!
        result = hash(result, item);
      }
    }
    return result;
  }

  private static int firstTerm(int aSeed) {
    return DEFAULT_PRIME_NUMBER * aSeed;
  }

  private static boolean isArray(Object aObject) {
    return aObject.getClass().isArray();
  }
}