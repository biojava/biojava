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

import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.LocationTools;
import org.biojava.utils.ChangeVetoException;

/**
 * <p><code>AnnotationTools</code> is a set of static utility methods for
 * manipulating <code>Annotation</code>s and <code>AnnotationType</code>s.</p>
 *
 * <p>The methods allIn() and allOut() let you compare an Annotation to an
 * AnnotationType and produce a new Annotation with only those properties
 * explicitly constrained by or not constrained by the type. This could be
 * of use when using an Annotation as a template for some object. You could use
 * allOut to make an Annotation that has all the properties that do not fit into
 * normal constructor properties, and pass that in as the Annotation bundle.</p>
 *
 * <p>intersection(AnnotationType) and union(AnnotationType) return new
 * AnnotationType instances that will accept every Annotation instance that is
 * accepted by both or either respectively. It is particularly informative to
 * compare the result of this to the AnnotationType.NONE to see if the two types
 * are mutualy disjoint.</p>
 *
 * <p>intersection(PropertyConstraint) and union(PropertyConstraint) return new
 * PropertyConstraint instances that will accept every Object that is accepted
 * by both or either one respectively.</p>
 * 
 * FilterTools uses these methods
 * when comparing filters on features by their Annotation bundles.
 * 
 * @since 1.3
 * @author Matthew Pocock
 * @author <a href="mailto:kdj@sanger.ac.uk">Keith James</a> (docs)
 * @author Thomas Down
 *
 */
public final class AnnotationTools {
    /**
     * <p>
     * Destructive down-cast an annotation to a type.
     * </p>
     *
     * <p>
     * <code>allIn</code> returns a new <code>Annotation</code>
     * containing only those values in the <code>Annotation</code>
     * argument which are of a type specified by the
     * <code>AnnotationType</code>.
     * </p>
     *
     * @param annotation an <code>Annotation</code> to scan.
     * @param annType an <code>AnnotationType</code>.
     *
     * @return an <code>Annotation</code>.
     */
    public static Annotation allIn(Annotation annotation, AnnotationType annType) {
        Annotation res;
        if (annotation instanceof SmallAnnotation) {
            res = new SmallAnnotation();
        } else {
            res = new SimpleAnnotation();
        }

        for (Iterator i = annType.getProperties().iterator(); i.hasNext();) {
            Object tag = i.next();
            try {
                res.setProperty(tag, annotation.getProperty(tag));
            } catch (ChangeVetoException cve) {
                throw new BioError("Assertion Failure: Can't alter an annotation", cve);
            }
        }

        return res;
    }

    /**
     * <code>allOut</code> returns a new <code>Annotation</code>
     * containing only those values in the <code>Annotation</code>
     * argument which are <strong>not</strong> of a type specified by
     * the <code>AnnotationType</code>.
     *
     * @param annotation an <code>Annotation</code>.
     * @param annType an <code>AnnotationType</code>.
     *
     * @return an <code>Annotation</code> value.
     */
    public static Annotation allOut(Annotation annotation, AnnotationType annType) {
        Annotation res;
        if (annotation instanceof SmallAnnotation) {
            res = new SmallAnnotation();
        } else {
            res = new SimpleAnnotation();
        }

        Set props = annType.getProperties();
        for (Iterator i = annotation.keys().iterator(); i.hasNext();) {
            Object tag = i.next();
            if (! props.contains(tag)) {
                try {
                    res.setProperty(tag, annotation.getProperty(tag));
                } catch (ChangeVetoException cve) {
                    throw new BioError("Assertion Failure: Can't alter an annotation", cve);
                }
            }
        }

        return res;
    }

    /**
     * <p>
     * Scans an Annotation with an AnnotationType and returns all Annotation
     * instances matching a Type.
     * </p>
     *
     * <p>This differs from AnnotationType.instanceOf()
     * as it will descend into properties of an Annotation if that property is
     * itself an Annotation. This allows you to scan a tree of Annotations for
     * nodes in the tree of a particular shape.
     * </p>
     *
     * @param ann  the Annotation to scan
     * @param query  the AnnotationType to match against all nodes in the tree
     * @return the set of all annotations matching the query
     */
    public static Set searchAnnotation(Annotation ann, AnnotationType query) {
      Set hits = new HashSet();
      searchAnnotation(ann, query, hits);
      return hits;
    }

    private static void searchAnnotation(Annotation ann, AnnotationType query, Set hits) {
      if(query.instanceOf(ann)) {
        hits.add(ann);
      }

      for(Iterator i = ann.keys().iterator(); i.hasNext(); ) {
        Object prop = i.next();
        Object val = ann.getProperty(prop);
        if(val instanceof Annotation) {
          searchAnnotation((Annotation) val, query, hits);
        } else if(prop instanceof Collection) {
          for(Iterator vi = ((Collection) val).iterator(); vi.hasNext(); ) {
            Object v = vi.next();
            if(v instanceof Annotation) {
              searchAnnotation((Annotation) v, query, hits);
            }
          }
        }
      }
    }

    /**
     * Calculate an AnnotationType that matches all Annotation instances matched
     * by both types. Usually you will either use this value blind or compare it to
     * AnnotationType.NONE.
     *
     * @param ann1  the first AnnotationType
     * @param ann2  the seccond AnnotationType
     * @return the intersection AnnotationType
     */
    public static AnnotationType intersection(
      AnnotationType ann1,
      AnnotationType ann2
    ) {
      if(ann1.subTypeOf(ann2)) {
        return ann2;
      } else if(ann2.subTypeOf(ann1)) {
        return ann1;
      } else {
        Set props = new HashSet();
        props.addAll(ann1.getProperties());
        props.addAll(ann2.getProperties());

        AnnotationType.Impl intersect = new AnnotationType.Impl();
        for(Iterator i = props.iterator(); i.hasNext(); ) {
          Object key = i.next();

          CollectionConstraint pc1 = ann1.getConstraint(key);
          CollectionConstraint pc2 = ann2.getConstraint(key);
          CollectionConstraint pc = intersection(pc1, pc2);
          if (pc == CollectionConstraint.NONE) {
            return AnnotationType.NONE;
          }

          intersect.setConstraint(key, pc);
        }

        intersect.setDefaultConstraint(
          intersection(ann1.getDefaultConstraint(), ann2.getDefaultConstraint())
        );

        return intersect;
      }
    }

    /**
     * Calculate the intersection of two PropertyConstraint instances. This method is realy only interesting when comparing each property in an
     * AnnotationType in turn. Usually the return value is either compared to
     * PropertyConstraint.NONE or is used blindly.
     *
     * @param pc1 the first PropertyConstraint
     * @param pc2 the seccond PropertyConstraint
     * @return the intersection PropertyConstraint
     *
     */
    public static PropertyConstraint intersection(
      PropertyConstraint pc1,
      PropertyConstraint pc2
    ) {
      if(pc1.subConstraintOf(pc2)) {
        return pc2;
      } else if(pc2.subConstraintOf(pc1)) {
        return pc1;
      } else if(
        pc1 instanceof PropertyConstraint.ByClass &&
        pc2 instanceof PropertyConstraint.ByClass
      ) {
        PropertyConstraint.ByClass pc1c = (PropertyConstraint.ByClass) pc1;
        PropertyConstraint.ByClass pc2c = (PropertyConstraint.ByClass) pc2;
        Class c1 = pc1c.getPropertyClass();
        Class c2 = pc2c.getPropertyClass();

        if(!c1.isInterface() && !c2.isInterface()) {
          return new PropertyConstraint.And(pc1c, pc2c);
        } else {
          return PropertyConstraint.NONE;
        }
      } else if(pc2 instanceof PropertyConstraint.ByClass) {
        return intersection(pc2, pc1);
      } else if(pc1 instanceof PropertyConstraint.ByClass) {
        PropertyConstraint.ByClass pc1c = (PropertyConstraint.ByClass) pc1;

        if(pc2 instanceof PropertyConstraint.Enumeration) {
          PropertyConstraint.Enumeration pc2e = (PropertyConstraint.Enumeration) pc2;
          Set values = new HashSet();
          for(Iterator i = pc2e.getValues().iterator(); i.hasNext(); ) {
            Object val = i.next();
            if(pc1c.accept(val)) {
              values.add(val);
            }
          }
          if(values.isEmpty()) {
            return PropertyConstraint.NONE;
          } else if(values.size() == 1) {
            return new PropertyConstraint.ExactValue(values.iterator().next());
          } else {
            return new PropertyConstraint.Enumeration(values);
          }
        }

        if(pc2 instanceof PropertyConstraint.ExactValue) {
          // we've already checked for containment - we know this value is of
          // the wrong class
          return PropertyConstraint.NONE;
        }
      } else if(
        (pc1 instanceof PropertyConstraint.Enumeration ||
         pc1 instanceof PropertyConstraint.ExactValue) &&
        (pc2 instanceof PropertyConstraint.Enumeration ||
         pc2 instanceof PropertyConstraint.ExactValue)
      ) {
        if (pc1 instanceof PropertyConstraint.Enumeration && pc2 instanceof PropertyConstraint.Enumeration) {
            Set intersection = new HashSet(((PropertyConstraint.Enumeration) pc1).getValues());
            intersection.retainAll(((PropertyConstraint.Enumeration) pc2).getValues());
            if (intersection.size() == 0) {
                return PropertyConstraint.NONE;
            } else if (intersection.size() == 1) {
                return new PropertyConstraint.ExactValue(intersection.iterator().next());
            } else {
                return new PropertyConstraint.Enumeration(intersection);
            }
        } else {
            // This case already handled by subset/superset logic
            return PropertyConstraint.NONE;
        }
      } else if(
        (pc1 instanceof PropertyConstraint.ByAnnotationType &&
         !(pc2 instanceof PropertyConstraint.ByAnnotationType)) ||
        (pc2 instanceof PropertyConstraint.ByAnnotationType &&
         !(pc1 instanceof PropertyConstraint.ByAnnotationType))
      ) {
        return PropertyConstraint.NONE;
      } else if(
        pc1 instanceof PropertyConstraint.ByAnnotationType &&
        pc2 instanceof PropertyConstraint.ByAnnotationType
      ) {
        PropertyConstraint.ByAnnotationType pc1a = (PropertyConstraint.ByAnnotationType) pc1;
        PropertyConstraint.ByAnnotationType pc2a = (PropertyConstraint.ByAnnotationType) pc2;

        AnnotationType intersect = intersection(
          pc1a.getAnnotationType(),
          pc2a.getAnnotationType()
        );
        if(intersect == AnnotationType.NONE) {
          return PropertyConstraint.NONE;
        } else {
          return new PropertyConstraint.ByAnnotationType(intersect);
        }
      }

      return new PropertyConstraint.And(pc1, pc2);
    }

    /**
     * Create an AnnotationType that matches all Anntotations that are accepted
     * by two others. This method is realy not very usefull in most cases. You may wish to
     * compare the result of this to AnnotationType.ANY, or use it blindly.
     *
     * @param ann1  the first AnnotationType
     * @param ann2  the seccond AnnotationType
     * @return an AnnotationType that represents their unions
     *
     */
    public static AnnotationType union(
      AnnotationType ann1,
      AnnotationType ann2
    ) {
      if(ann1.subTypeOf(ann2)) {
        return ann1;
      } else if(ann2.subTypeOf(ann1)) {
        return ann2;
      } else {
        Set props = new HashSet();
        props.addAll(ann1.getProperties());
        props.addAll(ann2.getProperties());

        AnnotationType.Impl union = new AnnotationType.Impl();
        for(Iterator i = props.iterator(); i.hasNext(); ) {
          Object key = i.next();

          CollectionConstraint pc1 = ann1.getConstraint(key);
          CollectionConstraint pc2 = ann2.getConstraint(key);
          CollectionConstraint pc = union(pc1, pc2);

          union.setConstraint(key, pc);
        }

        return union;
      }
    }

    /**
     * Create a PropertyConstraint that matches all Objects that are accepted
     * by two others. In the general case, there is no clean way to represent the union of two
     * PropertyConstraint instances. You may get back a PropertyConstraint.Or
     * instance, or perhaps PropertyConstraint.ANY. Alternatively, there may be
     * some comparrison possible. It is a thankless task introspecting this in
     * code. You have been warned.
     *
     * @param pc1 the first PropertyConstraint
     * @param pc2 the second PropertyConstraint
     * @return the union PropertyConstraint
     */
    public static PropertyConstraint union(
      PropertyConstraint pc1,
      PropertyConstraint pc2
    ) {
      if(pc1.subConstraintOf(pc2)) {
        return pc1;
      } else if(pc2.subConstraintOf(pc1)) {
        return pc2;
      } else if(
        pc1 instanceof PropertyConstraint.ByClass &&
        pc2 instanceof PropertyConstraint.ByClass
      ) {
        return new PropertyConstraint.Or(pc1, pc2);
      } else if(pc2 instanceof PropertyConstraint.ByClass) {
        return intersection(pc2, pc1);
      } else if(pc1 instanceof PropertyConstraint.ByClass) {
        PropertyConstraint.ByClass pc1c = (PropertyConstraint.ByClass) pc1;

        if(pc2 instanceof PropertyConstraint.Enumeration) {
          PropertyConstraint.Enumeration pc2e = (PropertyConstraint.Enumeration) pc2;
          Set values = new HashSet();
          for(Iterator i = pc2e.getValues().iterator(); i.hasNext(); ) {
            Object val = i.next();
            if(!pc1c.accept(val)) {
              values.add(val);
            }
          }
          if(values.isEmpty()) {
            return pc1;
          } else if(values.size() == 1) {
            return new PropertyConstraint.Or(
              pc1,
              new PropertyConstraint.ExactValue(values.iterator().next())
            );
          } else {
            return new PropertyConstraint.Or(
              pc1,
              new PropertyConstraint.Enumeration(values)
            );
          }
        }

        if(pc2 instanceof PropertyConstraint.ExactValue) {
          // we've already checked for containment - we know this value is of
          // the wrong class
          return new PropertyConstraint.Or(pc1, pc2);
        }
      } else if(
        pc1 instanceof PropertyConstraint.ByAnnotationType &&
        pc2 instanceof PropertyConstraint.ByAnnotationType
      ) {
        PropertyConstraint.ByAnnotationType pc1a = (PropertyConstraint.ByAnnotationType) pc1;
        PropertyConstraint.ByAnnotationType pc2a = (PropertyConstraint.ByAnnotationType) pc2;

        return new PropertyConstraint.ByAnnotationType(union(
          pc1a.getAnnotationType(),
          pc2a.getAnnotationType()
        ));
      }

      return new PropertyConstraint.Or(pc1, pc2);
    }

    /**
     * Return the CollectionConstraint which accept only collections accepted by
     * both of those specified.
     *
     * @param cc1 the first CollectionConstraint
     * @param cc2 the seccond CollectionConstrant
     * @return a CollectionConstraint representing the intersection of the other
     *    two
     */

    public static CollectionConstraint intersection(CollectionConstraint cc1, CollectionConstraint cc2) {
        if (cc1.subConstraintOf(cc2)) {
            return cc2;
        } else if (cc2.subConstraintOf(cc1)) {
            return cc1;
        } else if (cc1 instanceof CollectionConstraint.AllValuesIn &&
                   cc2 instanceof CollectionConstraint.AllValuesIn)
        {
            PropertyConstraint pc1 = ((CollectionConstraint.AllValuesIn) cc1).getPropertyConstraint();
            PropertyConstraint pc2 = ((CollectionConstraint.AllValuesIn) cc2).getPropertyConstraint();
            Location card1 = ((CollectionConstraint.AllValuesIn) cc1).getCardinalityConstraint();
            Location card2 = ((CollectionConstraint.AllValuesIn) cc2).getCardinalityConstraint();
            Location card = LocationTools.intersection(card1, card2);
            if (card == Location.empty) {
                return CollectionConstraint.NONE;
            }
            PropertyConstraint pc = intersection(pc1, pc2);
            if (pc == PropertyConstraint.NONE && !card.contains(0)) {
                return CollectionConstraint.NONE;
            } else {
                return new CollectionConstraint.AllValuesIn(pc, card);
            }
        } else if (cc1 instanceof CollectionConstraint.Contains &&
                   cc2 instanceof CollectionConstraint.Contains)
        {
            PropertyConstraint pc1 = ((CollectionConstraint.Contains) cc1).getPropertyConstraint();
            PropertyConstraint pc2 = ((CollectionConstraint.Contains) cc2).getPropertyConstraint();
            Location card1 = ((CollectionConstraint.Contains) cc1).getCardinalityConstraint();
            Location card2 = ((CollectionConstraint.Contains) cc2).getCardinalityConstraint();
            Location card = LocationTools.intersection(card1, card2);
            if (card == Location.empty) {
                return CollectionConstraint.NONE;
            }
            PropertyConstraint pc = intersection(pc1, pc2);
            if (pc == PropertyConstraint.NONE && !card.contains(0)) {
                return CollectionConstraint.NONE;
            } else {
                return new CollectionConstraint.Contains(pc, card);
            }
        } else if (cc1 instanceof CollectionConstraint.Contains &&
                   cc2 instanceof CollectionConstraint.AllValuesIn)
        {
            PropertyConstraint pc1 = ((CollectionConstraint.Contains) cc1).getPropertyConstraint();
            PropertyConstraint pc2 = ((CollectionConstraint.AllValuesIn) cc2).getPropertyConstraint();
            Location card1 = ((CollectionConstraint.Contains) cc1).getCardinalityConstraint();
            Location card2 = ((CollectionConstraint.AllValuesIn) cc2).getCardinalityConstraint();
            if (card1.getMin() > card2.getMax()) {
                // Requires too many values.
                return CollectionConstraint.NONE;
            }
            PropertyConstraint pc = intersection(pc1, pc2);
            if (pc == PropertyConstraint.NONE && !card1.contains(0)) {
                return CollectionConstraint.NONE;
            } else {
                return new CollectionConstraint.Contains(pc, card1);
            }
        } else if (cc1 instanceof CollectionConstraint.AllValuesIn &&
                   cc2 instanceof CollectionConstraint.Contains)
        {
            return intersection(cc2, cc1);
        } else {
            return new CollectionConstraint.And(cc1, cc2);
        }
    }

  /**
   * Calculate a CollectionConstaint that will accept all items accepted by
   * either constraint.
   *
   * @param cc1   the first CollectionConstraint
   * @param cc2   the seccond collectionConstraint
   * @return      a CollectionConstraint representing the union of the other two
   */
    public static CollectionConstraint union(CollectionConstraint cc1, CollectionConstraint cc2) {
        if (cc1.subConstraintOf(cc2)) {
            return cc1;
        } else if (cc2.subConstraintOf(cc1)) {
            return cc2;
        } else {
            return new CollectionConstraint.Or(cc1, cc2);
        }
    }
}
