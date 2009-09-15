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
import java.util.Iterator;

import org.biojava.bio.symbol.Location;

/**
 * Used by <code>AnnotationType</code> to represent the constraint on
 * the collection of values in a property-slot.
 * CollectionConstraints usually use a <code>PropertyConstraint</code>
 * to validate the individual elements.
 *
 * 
 * Use one or more of the built-in implementations to build new
 * <code>AnnotationTypes</code>.
 *
 * @since 1.3
 * @author Thomas Down
 * @author Matthew Pocock
 */
public interface CollectionConstraint {
    /**
     * <code>accept</code> returns true if the value fulfills the
     * constraint.
     *
     * @param values a <code>Collection</code> to check.
     * @return true if the values are acceptable
     *
     * powerUser Manually compare items with the CollectionConstraint. Node:
     * this will ususaly be done for you in an AnnotationType instance
     */
    public boolean accept(Object values);

    /**
     * <p><code>subConstraintOf</code> returns true if the constraint
     * is a sub-constraint.<p>
     *
     * <p>A pair of constraints super and sub are in a
     * superConstraint/subConstraint relationship if every object
     * accepted by sub is also accepted by super. To put it another
     * way, if instanceOf was used as a set-membership indicator
     * function over some set of objects, then the set produced by
     * super would be a superset of that produced by sub.</p>
     *
     * <p>It is not expected that constraints will neccesarily
     * maintain references to super/sub types. It will be more usual
     * to infer this relationship by introspecting the constraints
     * themselves. For example,
     * <code>CollectionConstraint.ByClass</code> will infer
     * subConstraintOf by looking at the possible class of all items
     * matching subConstraint.</p>
     *
     * @param subConstraint a <code>CollectionConstraint</code> to check.
     * @return a <code>boolean</code>.
     *
     * 
     * Usefull when attempting to compare two constraints to see
     * if it is necisary to retain both. You may want to check the more
     * general or the more specific constraint only.
     */
    public boolean subConstraintOf(CollectionConstraint subConstraint);


    /**
     * Return <code>true</code> iff the Collection formed by adding
     * <code>newValue</code> to <code>current</code> would be accepted
     * by this constraint.
     *
     * 
     * Implementations may <em>not</em> assume that <code>current</code>
     * is valid.
     *
     * @param current  a Collection containing the current values
     * @param newValue the new value to add
     * @return true if adding the new value will result in an acceptable
     *    property
     */
    
    public boolean validateAddValue(Collection current, Object newValue);
    
    /**
     * Return <code>true</code> iff the Collection formed by removing
     * <code>newValue</code> from <code>current</code> would be accepted
     * by this constraint.
     *
     * 
     * Implementations may <em>not</em> assume that <code>current</code>
     * is valid.  However, <code>current</code> will already have been
     * checked to ensure that it contains <code>victim</code>.
     *
     * @param current a Collection containing the current values
     * @param victim  the value to remove
     * @return true   if removing the victim will result in an acceptable
     *    property value set
     */
    
    public boolean validateRemoveValue(Collection current, Object victim);
    
    /**
     * <code>ANY</code> is a constraint which accepts a property for
     * addition under all conditions.
     *
     * Whenever a CollectionConstraint is needed and you want to allow
     * any value there
     */
    public static final CollectionConstraint ANY = new AllValuesIn(PropertyConstraint.ANY, CardinalityConstraint.ANY);
    
    /**
     * <code>EMPTY</code> is a constraint which only accepts the empty
     * set.
     * 
     * Use this to indicate that a property must be undefined
     */
     
    public static final CollectionConstraint EMPTY = new AllValuesIn(PropertyConstraint.NONE, CardinalityConstraint.ZERO);
    
    /**
     * <code>NONE</code> is a constraint which accepts no value for a property
     * under any condition.
     *
     * This value indicates an impossible condition.  It may be
     *                 returned by methods such as <code>AnnotationTools.intersection</code>
     *                 to indicate that <code>NO</code> values of a property (include undefined)
     *                 are valid.
     */
    public static final CollectionConstraint NONE = new NoneCollectionConstraint();
  
    /**
     * CollectionConstraint which validates all members of a Collection.
     * All members must be vaild according to the supplied
     * <code>PropertyConstraint</code>, and the total number of
     * members must be acceptable by the given cardinality constraint.
     *
     * @author Thomas Down
     * @author Matthew Pocock
     */
    
    public class AllValuesIn implements CollectionConstraint {
        private PropertyConstraint pc;
        private Location card;

      /**
       * Create an AllValuesIn based upon a PropertyConstraint and a
       * cardinality.
       *
       * @param pc    the PropertyConstraint to apply to each property value
       * @param card  the cardinality constraint restricting the number of
       *    values
       */
        public AllValuesIn(PropertyConstraint pc, Location card) {
            this.pc = pc;
            this.card = card;
        }

      /**
       * Get the PropertyConstraint used to validate each property value.
       *
       * @return  the PropertyConstraint used
       */
        public PropertyConstraint getPropertyConstraint() {
            return pc;
        }

      /**
       * Get the cardinality constraint used to validate the number of property
       * values.
       *
       * @return  the cardinality constraint as a Location
       */
        public Location getCardinalityConstraint() {
            return card;
        }
        
        public boolean accept(Object o) {
            if (o instanceof Collection) {
                Collection co = (Collection) o;
                if (!card.contains(co.size())) {
                    return false;
                } else {
                    for (Iterator i = co.iterator(); i.hasNext(); ) {
                        if (!pc.accept(i.next())) {
                            return false;
                        }
                    }
                    return true;
                }
            } else {
                return card.contains(1) && pc.accept(o);
            }
        }
        
        public boolean validateAddValue(Collection oldcol, Object newValue) {
            if (!pc.accept(newValue)) {
                return false;
            }
            if (!card.contains(oldcol.size() + 1)) {
                return false;
            }
            for (Iterator i = oldcol.iterator(); i.hasNext(); ) {
                if (!pc.accept(i.next())) {
                    return false;
                }
            }
            return true;
        }
        
        public boolean validateRemoveValue(Collection oldcol, Object victim) {
            if (!card.contains(oldcol.size() - 1)) {
                return false;
            }
            
            for (Iterator i = oldcol.iterator(); i.hasNext(); ) {
                Object o = i.next();
                if (!o.equals(victim) && !pc.accept(o)) {
                    return false;
                }
            }
            return true;
        }
        
        public int hashCode() {
            return pc.hashCode() + 87;
        }
        
        public boolean equals(Object o) {
            if (o instanceof AllValuesIn) {
                AllValuesIn avo = (AllValuesIn) o;
                return avo.getCardinalityConstraint().equals(getCardinalityConstraint()) &&
                       avo.getPropertyConstraint().equals(getPropertyConstraint());
            } else {
                return false;
            }
        }
        
        public boolean subConstraintOf(CollectionConstraint cc) {
            if (cc instanceof NoneCollectionConstraint) {
                return true;
            } else if (cc instanceof CollectionConstraint.AllValuesIn) {
                AllValuesIn ccavi = (AllValuesIn) cc;
                return getCardinalityConstraint().contains(ccavi.getCardinalityConstraint()) &&
                       getPropertyConstraint().subConstraintOf(ccavi.getPropertyConstraint());
            } else if (cc instanceof CollectionConstraint.Contains) {
                if (!getCardinalityConstraint().contains(Integer.MAX_VALUE)) {
                    return false;
                } else {
                    Contains ccc = (Contains) cc;
                    return getPropertyConstraint().subConstraintOf(ccc.getPropertyConstraint());
                }
            } else if (cc instanceof CollectionConstraint.And) {
                And cca = (And) cc;
                return subConstraintOf(cca.getChild1()) || subConstraintOf(cca.getChild2());
            } else if (cc instanceof CollectionConstraint.Or) {
                Or cco = (Or) cc;
                return subConstraintOf(cco.getChild1()) && subConstraintOf(cco.getChild2());
            }
            return false;
        }
        
        public String toString() {
          return "AllValuesIn: (" + pc.toString() + ", " + card.toString() + ")";
        }
    }
    
    /**
     * CollectionConstraint which validates a portion of a Collection.
     * Accepts only collections where the number of members matching
     * the <code>PropertyConstraint</code> is in the supplied cardinality.
     *
     * <p>
     * A typical application for this would be with Annotations where
     * one property can contain a number of synonyms.
     * <code>CollectionConstraint.Contains</code> could be used as
     * a query to select instances based on one of these synonyms.
     * </p>
     *
     * @author Thomas Down
     */
    
    public class Contains implements CollectionConstraint {
        private PropertyConstraint pc;
        private Location card;
        
      /**
       * Create a Contains based upon a PropertyConstraint and a
       * cardinality.
       *
       * @param pc    the PropertyConstraint to apply to each property value
       * @param card  the cardinality constraint restricting the number of
       *    values
       */
        public Contains(PropertyConstraint pc, Location card) {
            this.pc = pc;
            this.card = card;
        }
        
      /**
       * Get the PropertyConstraint used to validate each property value.
       *
       * @return  the PropertyConstraint used
       */
        public PropertyConstraint getPropertyConstraint() {
            return pc;
        }
        
      /**
       * Get the cardinality constraint used to validate the number of property
       * values.
       *
       * @return  the cardinality constraint as a Location
       */
        public Location getCardinalityConstraint() {
            return card;
        }
        
        public boolean accept(Object o) {
            if (o instanceof Collection) {
                return card.contains(countMembers((Collection) o));
            } else {
                if (pc.accept(o)) {
                    return card.contains(1);
                } else {
                    return card.contains(0);
                }
            }
        }
        
        private int countMembers(Collection co) {
            int members = 0;
            for (Iterator i = co.iterator(); i.hasNext(); ) {
                if (pc.accept(i.next())) {
                    ++members;
                }
            }
            return members;
        }
        
        public int hashCode() {
            return pc.hashCode() + 178;
        }
        
        public boolean equals(Object o) {
            if (o instanceof Contains) {
                Contains avo = (Contains) o;
                return avo.getCardinalityConstraint().equals(getCardinalityConstraint()) &&
                       avo.getPropertyConstraint().equals(getPropertyConstraint());
            } else {
                return false;
            }
        }
        
        public boolean validateAddValue(Collection oldCol, Object newValue) {
            int members = countMembers(oldCol);
            if (pc.accept(newValue)) {
                ++members;
            }
            return card.contains(members);
        }
        
        public boolean validateRemoveValue(Collection oldCol, Object newValue) {
            int members = countMembers(oldCol);
            if (pc.accept(newValue)) {
                --members;
            }
            return card.contains(members);
        }
        
        
        public boolean subConstraintOf(CollectionConstraint cc) {
            if (cc instanceof NoneCollectionConstraint) {
                return true;
            } else if (cc instanceof CollectionConstraint.AllValuesIn) {
                AllValuesIn ccavi = (AllValuesIn) cc;
                return getCardinalityConstraint().contains(ccavi.getCardinalityConstraint()) &&
                       getPropertyConstraint().subConstraintOf(ccavi.getPropertyConstraint());
            } else if (cc instanceof CollectionConstraint.Contains) {
                Contains ccavi = (Contains) cc;
                return getCardinalityConstraint().contains(ccavi.getCardinalityConstraint()) &&
                       getPropertyConstraint().subConstraintOf(ccavi.getPropertyConstraint());
            } else if (cc instanceof CollectionConstraint.And) {
                And cca = (And) cc;
                return subConstraintOf(cca.getChild1()) || subConstraintOf(cca.getChild2());
            } else if (cc instanceof CollectionConstraint.Or) {
                Or cco = (Or) cc;
                return subConstraintOf(cco.getChild1()) && subConstraintOf(cco.getChild2());
            }
            return false;
        }
        
        public String toString() {
          return "Contains: (" + pc.toString() + ", " + card.toString() + ")";
        }
    }
    
    
    /**
     * A collection constraint that accpepts collections iff they are accepted by both
     * child constraints. This effectively matches the intersection of the items
     * matched by the two constraints.
     *
     * Use this to combine multiple constraints. You can make one
     *            or both of the children And instances if you need a tighter
     *            intersection.
     * @author Matthew Pocock
     * @author Thomas Down
     */
    public class And implements CollectionConstraint {
      private CollectionConstraint c1;
      private CollectionConstraint c2;
      
      /**
       * Create a new <code>And</code> from two child constraints.
       *
       * @param c1 the first child
       * @param c2 the seccond child
       */
      public And(CollectionConstraint c1, CollectionConstraint c2) {
        this.c1 = c1;
        this.c2 = c2;
      }
      
      /**
       * Get the first child CollectionConstraint.
       *
       * @return the first child CollectionConstraint
       *
       */
      public CollectionConstraint getChild1() {
        return c1;
      }
      
      /**
       * Get the seccond child CollectionConstraint.
       *
       * @return the seccond child CollectionConstraint
       *
       */
      public CollectionConstraint getChild2() {
        return c2;
      }
      
      public boolean accept(Object object) {
        return c1.accept(object) && c2.accept(object);
      }
      
      public boolean subConstraintOf(CollectionConstraint pc) {
        return c1.subConstraintOf(pc) && c2.subConstraintOf(pc);
      }
      
      
      public boolean validateAddValue(Collection oldcoll, Object newvalue) {
          return c1.validateAddValue(oldcoll, newvalue) && c2.validateAddValue(oldcoll, newvalue);
      }
      
      public boolean validateRemoveValue(Collection oldcoll, Object victim) {
          return c1.validateAddValue(oldcoll, victim) && c2.validateAddValue(oldcoll, victim);
      }
      
      public String toString() {
        return "And(" + c1 + ", " + c2 + ")";
      }
    }
    
    /**
     * A collection constraint that accepts items iff they are accepted by either
     * child constraints. This effectively matches the union of the items
     * matched by the two constraints. Use this to combine multiple constraints. You can make one
     *            or both of the children Or instances if you need a wider
     *            union.
     *
     * @author Matthew Pocock
     * @author Thomas Down
     */
    public class Or implements CollectionConstraint {
      private CollectionConstraint c1;
      private CollectionConstraint c2;
      
      /**
       * Create a new <code>Or</code> from two child constraints.
       *
       * @param c1 the first child
       * @param c2 the seccond child
       */
      public Or(CollectionConstraint c1, CollectionConstraint c2) {
        this.c1 = c1;
        this.c2 = c2;
      }
      
      /**
       * Get the first child CollectionConstraint.
       *
       * @return the first child CollectionConstraint
       */
      public CollectionConstraint getChild1() {
        return c1;
      }
      
      /**
       * Get the seccond child CollectionConstraint.
       *
       * @return the seccond child CollectionConstraint
       * 
       */
      public CollectionConstraint getChild2() {
        return c2;
      }
      
      public boolean accept(Object object) {
        return c1.accept(object) || c2.accept(object);
      }
      
      public boolean subConstraintOf(CollectionConstraint pc) {
        return c1.subConstraintOf(pc) || c2.subConstraintOf(pc);
      }
      
      public boolean validateAddValue(Collection oldcoll, Object newvalue) {
          return c1.validateAddValue(oldcoll, newvalue) || c2.validateAddValue(oldcoll, newvalue);
      }
      
      public boolean validateRemoveValue(Collection oldcoll, Object victim) {
          return c1.validateAddValue(oldcoll, victim) || c2.validateAddValue(oldcoll, victim);
      }
      
      public String toString() {
        return "Or(" + c1 + ", " + c2 + ")";
      }
    }
}

