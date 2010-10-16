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
package org.biojava.utils;

/**
 * utility methods for implementing the equals() and hashCode() methods of Objects.
 * All credit for this class goes to Mark Davis (Java Report 5(1), 46; Java Report 5(4), 60).
 * <p>
 * All equals() methods in this class take the two fields to compare for equality as arguments and return whether they
 * are equal. Consequently, the equals() method of class AClass should be implemented (taking advantage of the equals() 
 * methods in this class) as follows:
 * <pre>
 * public boolean equals(Object o) {
 *   if (o == this) return true;
 *   // if this class AClass is a direct sub-class of Object:
 *   if (o == null) return false;
 *   if (!o.getClass().equals(this.getClass())) return false;
 *   // else if this class AClass extends another class than Object:
 *   //if (!super.equals(o)) return false;
 *   
 *   AClass that = (AClass) o;
 *   
 *   // only compare fields of this class AClass (not of super-classes):
 *   //if (!ObjectUtil.equals(this.field, that.field)) return false;
 *   
 *   // this and that are identical if we made it 'til here
 *   return true;
 * }
 * </pre>
 * <p>
 * All hashCode() methods in this class take the current hashCode value as the first argument and the field to take into
 * account as the second argument and return the newly calculated hashCode value (which now includes the influence of
 * the given field). Consequently, the hashCode() method of class AClass should be implemented (taking advantage of the
 * hashCode() methods in this class) as follows:
 * <pre>
 * public int hashCode() {
 *   // if this class AClass is a direct sub-class of Object:
 *   int hc = 0;
 *   // else if this class AClass extends another class than Object:
 *   //int hc = super.hashCode();
 *
 *   // only take into account fields of this class AClass (not of super-class):
 *   //hc = ObjectUtil.hashCode(hc, field);
 *
 *   return hc;
 * }
 * </pre>
 *
 * @author <A href="mailto:Gerald.Loeffler@vienna.at">Gerald Loeffler</A>
 */
public final class ObjectUtil {
  /**
   * the current hashCode is always first multiplied with this prime before the hashCode value for a particular field is
   * added.
   */
  public static final int PRIME = 1000003;

  public static boolean equals(boolean b1, boolean b2) {
    return b1 == b2;
  }

  public static boolean equals(int i1,  int i2) {
    return i1 == i2;
  }

  public static boolean equals(long l1,  long l2) {
    return l1 == l2;
  }

  public static boolean equals(float f1,  float f2) {
    return ((Float.isNaN(f1) && Float.isNaN(f2)) || (f1 == f2));
  }

  public static boolean equals(double d1, double d2) {
    return ((Double.isNaN(d1) && Double.isNaN(d2)) || (d1 == d2));
  }

  public static boolean equals(boolean[] a1, boolean[] a2) {
    if (a1 == null && a2 == null) return true;
    if (!(a1 != null && a2 != null)) return false;

    final int l = a1.length;
    if (l != a2.length) return false;

    for (int i = 0; i < l; i++) {
      if (!equals(a1[i], a2[i])) return false;
    }
    return true;
  }

  public static boolean equals(int[] a1, int[] a2) {
    if (a1 == null && a2 == null) return true;
    if (!(a1 != null && a2 != null)) return false;

    final int l = a1.length;
    if (l != a2.length) return false;

    for (int i = 0; i < l; i++) {
      if (!equals(a1[i], a2[i])) return false;
    }
    return true;
  }

  public static boolean equals(long[] a1, long[] a2) {
    if (a1 == null && a2 == null) return true;
    if (!(a1 != null && a2 != null)) return false;

    final int l = a1.length;
    if (l != a2.length) return false;

    for (int i = 0; i < l; i++) {
      if (!equals(a1[i], a2[i])) return false;
    }
    return true;
  }

  public static boolean equals(float[] a1, float[] a2) {
    if (a1 == null && a2 == null) return true;
    if (!(a1 != null && a2 != null)) return false;

    final int l = a1.length;
    if (l != a2.length) return false;

    for (int i = 0; i < l; i++) {
      if (!equals(a1[i], a2[i])) return false;
    }
    return true;
  }

  public static boolean equals(double[] a1, double[] a2) {
    if (a1 == null && a2 == null) return true;
    if (!(a1 != null && a2 != null)) return false;

    final int l = a1.length;
    if (l != a2.length) return false;

    for (int i = 0; i < l; i++) {
      if (!equals(a1[i], a2[i])) return false;
    }
    return true;
  }

  public static boolean equals(Object[] a1, Object[] a2) {
    if (a1 == null && a2 == null) return true;
    if (!(a1 != null && a2 != null)) return false;

    final int l = a1.length;
    if (l != a2.length) return false;

    for (int i = 0; i < l; i++) {
      if (!equals(a1[i], a2[i])) return false;
    }
    return true;
  }

  public static boolean equals(Object o1, Object o2) {
    return ((o1 == null && o2 == null) || (o1 != null && o1.equals(o2)) || (o2 != null && o2.equals(o1)));
  }

  public static int hashCode(int currentHashCodeValue, boolean b) {
    return PRIME*currentHashCodeValue + (b ? 1 : 0);
  }

  public static int hashCode(int currentHashCodeValue, int i) {
    return PRIME*currentHashCodeValue + i;
  }

  public static int hashCode(int currentHashCodeValue, long l) {
    return PRIME*(PRIME*currentHashCodeValue + ((int) (l>>>32))) + ((int) (l&0xFFFFFFFF));
  }

  public static int hashCode(int currentHashCodeValue, float f) {
    return hashCode(currentHashCodeValue, f == 0.0F ? 0 : Float.floatToIntBits(f));
  }

  public static int hashCode(int currentHashCodeValue, double d) {
    return hashCode(currentHashCodeValue, d == 0.0 ? 0L : Double.doubleToLongBits(d));
  }

  public static int hashCode(int currentHashCodeValue, boolean[] a) {
    if (a == null) return PRIME*currentHashCodeValue;
    
    final int l = a.length;
    for (int i = 0; i < l; i++) {
      currentHashCodeValue = hashCode(currentHashCodeValue, a[i]);
    }
    return currentHashCodeValue;
  }

  public static int hashCode(int currentHashCodeValue, int[] a) {
    if (a == null) return PRIME*currentHashCodeValue;
    
    final int l = a.length;
    for (int i = 0; i < l; i++) {
      currentHashCodeValue = hashCode(currentHashCodeValue, a[i]);
    }
    return currentHashCodeValue;
  }

  public static int hashCode(int currentHashCodeValue, long[] a) {
    if (a == null) return PRIME*currentHashCodeValue;
    
    final int l = a.length;
    for (int i = 0; i < l; i++) {
      currentHashCodeValue = hashCode(currentHashCodeValue, a[i]);
    }
    return currentHashCodeValue;
  }

  public static int hashCode(int currentHashCodeValue, float[] a) {
    if (a == null) return PRIME*currentHashCodeValue;
    
    final int l = a.length;
    for (int i = 0; i < l; i++) {
      currentHashCodeValue = hashCode(currentHashCodeValue, a[i]);
    }
    return currentHashCodeValue;
  }

  public static int hashCode(int currentHashCodeValue, double[] a) {
    if (a == null) return PRIME*currentHashCodeValue;
    
    final int l = a.length;
    for (int i = 0; i < l; i++) {
      currentHashCodeValue = hashCode(currentHashCodeValue, a[i]);
    }
    return currentHashCodeValue;
  }

  public static int hashCode(int currentHashCodeValue, Object[] a) {
    if (a == null) return PRIME*currentHashCodeValue;
    
    final int l = a.length;
    for (int i = 0; i < l; i++) {
      currentHashCodeValue = hashCode(currentHashCodeValue, a[i]);
    }
    return currentHashCodeValue;
  }

  public static int hashCode(int currentHashCodeValue, Object o) {
    return PRIME*currentHashCodeValue + (o == null ? 0 : o.hashCode());
  }
}
