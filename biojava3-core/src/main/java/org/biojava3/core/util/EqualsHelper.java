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

/**
 * A set of helper methods which return true if the two paramters are
 * equal to each other. It should be used with Equals Methods where object
 * generation is a problem therefore stopping you from using Commons-Lang
 * EqualsBuilder.
 *
 * @author ayates
 */
public class EqualsHelper {

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
    if (one == null && two == null) {
      equal = true;
    }
    else if (one == null || two == null) {
      equal = false;
    }
    else {
      // Same ref equals should be delt with in the equals()
      // method of the calling class and not here
      equal = one.equals(two);
    }
    return equal;
  }

  /**
   * This method should be called before begininng any equals methods. In order
   * to return true the method:
   *
   * <ol>
   * <li>The two given objects are the same instance using ==. This also means
   * if both Objects are null then this method will return true (well
   * techincally they are equal)</li>
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
      // We can replace this with
      // one.getClass().isInstance(two);
      // However the only advantage is the method version does a null check
      // which
      // we cannot be in because of the first if statement
      equal = (one.getClass() == two.getClass());
    }

    return equal;
  }
}
