/*
 *                 BioJava development code
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

import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

import org.biojava.utils.ChangeVetoException;

  /**
   * This class intendes to provide some FeatureHolder utilities. 
   * Currently it is mainly providing set operators.
   * 
   * @author Markus Brosch (markus[at]brosch[dot]cc)
   */
  public class FeatureHolderUtils {

    /**
     * Operator: Union of FeatureHolder1 and FeatureHolder2
     * @param fh1 FeatureHolder1
     * @param fh2 FeatureHolder2
     * @return Union of fh1 and fh2 (corresponds to logical OR)
     * @throws ChangeVetoException
     */
    public static FeatureHolder union(FeatureHolder fh1, FeatureHolder fh2)
        throws ChangeVetoException {

      int s1 = fh1.countFeatures();
      int s2 = fh2.countFeatures();

      if (s1 < s2) {
        return unionOp(fh2, fh1);
      } else {
        return unionOp(fh1, fh2);
      }
    }

    private static FeatureHolder unionOp(FeatureHolder fh1, FeatureHolder fh2)
        throws ChangeVetoException {

      SimpleFeatureHolder res = new SimpleFeatureHolder();
      for (Iterator it = fh1.features(); it.hasNext();) {
        res.addFeature((Feature) it.next());
      }
      for (Iterator it = fh2.features(); it.hasNext();) {
        Feature f = (Feature) it.next();
        if (!res.containsFeature(f)) {
          res.addFeature(f);
        }
      }
      return res;
    }

    /**
     * Operator: Intersect FeatureHolder1 with FeatureHolder2
     * @param fh1 FeatureHolder1
     * @param fh2 FeatureHolder2
     * @return Intersection of fh1 and fh2 (corresponds to logical AND)
     * @throws ChangeVetoException
     */
    public static FeatureHolder intersect(FeatureHolder fh1, FeatureHolder fh2)
        throws ChangeVetoException {

      int s1 = fh1.countFeatures();
      int s2 = fh2.countFeatures();

      if (s1 < s2) {
        return intersectOp(fh1, fh2);
      } else {
        return intersectOp(fh2, fh1);
      }
    }

    private static FeatureHolder intersectOp(FeatureHolder fh1,
        FeatureHolder fh2)
        throws ChangeVetoException {

      SimpleFeatureHolder res = new SimpleFeatureHolder();
      for (Iterator it = fh1.features(); it.hasNext();) {
        Feature f = (Feature) it.next();
        if (fh2.containsFeature(f)) {
          res.addFeature(f);
        }
      }
      return res;
    }

    /**
     * Operator: FeatureHolder 1 NOT FeatureHolder2
     * @param fh1 FeatureHolder1
     * @param fh2 FeatureHolder2
     * @return Set of fh1 without any feature of fh2 (Not)
     * @throws ChangeVetoException
     */
    public static FeatureHolder not(FeatureHolder fh1, FeatureHolder fh2)
        throws ChangeVetoException {

      SimpleFeatureHolder res = new SimpleFeatureHolder();
      for (Iterator it = fh1.features(); it.hasNext();) {
        Feature f = (Feature) it.next();
        if (!fh2.containsFeature(f)) {
          res.addFeature(f);
        }
      }
      return res;
    }

    /**
     * Returns a FeatureHolder as a Set of Features
     * @param fh FeatureHolder you want to have as a Set
     * @return Set of FeatureHoler fh
     */
    public static Set featureHolderAsSet(FeatureHolder fh) {
      return new FeatureHolderAsSet(fh);
    }

    /**
     * FeatureHolderAsSet represents a FeatureHolder as a Set.<br>
     * You can use the FeatureHolder as a normal Set.
     */
    private final static class FeatureHolderAsSet extends HashSet {

      FeatureHolderAsSet(FeatureHolder fh) {
        for (Iterator it = fh.features(); it.hasNext();) {
          add(it.next());
        }
      }
    }

  }
