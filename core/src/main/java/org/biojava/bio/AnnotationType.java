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

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

import org.biojava.bio.symbol.Location;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.SmallMap;

/**
 * A set of constraints on the data contained in an <code>Annotation</code>.
 * <code>AnnotationType</code> instances can be used to validate an
 * <code>Annotation</code> to check that it has the appropriate
 * properties and that they are of the right type.
 *
 * <h2>Common Usage</h2>
 *
 * <pre>
 * // annType is going to define ID as being a single String and
 * // AC as being a list of Strings and accept any other properties at all
 * AnnotationType annType = new AnnotationType.Impl();
 * annType.setDefaultConstraints(PropertyConstraint.ANY, Cardinality.ANY)
 * annType.setConstraints("ID",
 *                        new PropertyConstraint.ByClass(String.class),
 *                        CardinalityConstraint.ONE );
 * annType.setConstraint("AC",
 *                       new PropertyConstraint.ByClass(String.class),
 *                       CardinalityConstraint.ANY );
 *
 * Annotation ann = ...;
 *
 * if(annType.accepts(ann)) {
 *   // we have proved that ann contains one ID property that is a String
 *   System.out.println("ID: " + (String) ann.getProperty("ID"));
 *   // we know that the AC property is potentialy more than one String -
 *   // let's use getProperty on AnnotationType to make sure that the
 *   // values get unpacked cleanly
 *   System.out.println("AC:");
 *   for(Iterator i = annType.getProperty(ann, "AC"); i.hasNext(); ) {
 *     System.out.println("\t" + (String) i.next());
 *   }
 * } else {
 *   throw new IllegalArgumentException(
 *     "Expecting an annotation conforming to: "
 *     annType + " but got: " + ann
 *   );
 * }
 * </pre>
 *
 * <h2>Description</h2>
 *
 * <p>AnnotationType is a  constraint-based language for describing
 * sets of Annotation bundles. It works by assuming that any given Annotation
 * may have any set of properties defined. If it matches a particular
 * AnnotationType instance, then each defined property must be of a value
 * that is acceptable to the type, and each undefined property must
 * be allowed to be absend in the type.</p>
 *
 * <p>
 * <code>AnnotationType</code> imposes a data model on <code>Annotations</code>
 * where the value of each `slot' is considered as a <code>Collection</code>.
 * Values which aren't actually <code>Collection</code> instances are treated
 * as singleton sets.  Undefined properties are treated as <code>Collection.EMPTY_SET</code>.
 * The logic to validate a complete slot in an <code>Annotation</code> is
 * encapsulated by a <code>CollectionConstraint</code> object.  In the common case,
 * slots are validated by <code>CollectionConstraint.AllValuesIn</code> which
 * ensures that all members of the collection match a specified <code>PropertyConstraint</code>,
 * and that the size of the collection matches a cardinality constraint.  This is the
 * most important type of constraint when defining datatypes.  However, when using
 * <code>AnnotationTypes</code> to specify a query over a set of Annotations, it may be more
 * useful to ask whether a slot <em>contains</em> a given value:  For example
 * </p>
 *
 * <pre>
 * // Test that gene_name contains the value "BRCA2", ignoring other values (synonyms)
 * // which might be present in that slot.
 * AnnotationType at = new AnnotationType.Impl();
 * at.setConstraint(
 *      "gene_name",
 *      new CollectionConstraint.Contains(
 *          new PropertyConstraint.ExactValue("BRCA2"),
 *          CardinalityConstraint.ONE
 *      )
 * );
 * </pre>
 *
 * <p>It is usually left up to the AnnotationType instance to work out how
 * multiple values should be packed into a single property slot in an Annotation
 * instance. Commonly, things that are allowed a cardinality of 1 will store one
 * value directly in the slot. Things that allow multiple values (and optionaly
 * things with one value) will usualy store them within a Collection in the
 * slot. This complexity is hidden from you if you use the accessor methods
 * built into AnnotationType, setProperty(), addProperty(), removeProperty() and getProperty().
 * </p>
 *

 *
 * 
 * Using AnnotationType instances that you have been provided with
 * e.g. from UnigeneTools.LIBRARY_ANNOTATION
 *
 * 
 * <code>AnnotationType</code> instances can be used as queries
 * to select from a set of Annotations based on the value of one
 * or more properties.  Commonly, this is used in conjunction
 * with <code>FeatureFilter.ByAnnotationType</code>.
 *
 * 
 * Make AnnotationType instances that describe what should and
 * should not appear in an Annotation bundle
 *
 * 
 * Constrain FeatureFilter schemas by Annotation associated with
 * the features
 *
 * 
 * Provide meta-data to the tag-value parser for automatically
 * generating object representations of flat-files
 *
 * 
 * Implementing your own AnnotationType implementations to reflect
 * frame, schema or ontology definitions. For example, dynamically
 * reflect an RDMBS schema or DAML/Oil deffinition as an
 * AnnotationType.
 * 
 * @since 1.3
 * @author Matthew Pocock
 * @author Keith James (docs)
 * @author Thomas Down
 * 
 */
public interface AnnotationType {
    /**
     * The type that accepts all annotations and is the supertype of all
     * other annotations. Only an empty annotation is an exact instance of
     * this type.
     *
     *  Use this whenever an AnnotationType is needed by an API and you
     *       don't want to constrain anything
     */
    AnnotationType ANY = new Impl(
      PropertyConstraint.ANY,
      CardinalityConstraint.ANY
    );

    /**
     * The type that accepts no annotations at all and is the subtype of all
     * other annotations.
     *
     * Use this whenever an AnnotationType is needed by an API and you
     *       want to make sure that all Annotation objects get rejected
     */
    AnnotationType NONE = new Impl(
      PropertyConstraint.NONE,
      CardinalityConstraint.NONE
    );

    /**
     * Validate an Annotation against this AnnotationType.
     *
     *  Any time you wish to see if an Annotation bundle conforms to a
     *       type
     * @param ann the Annotation to validate.
     * @return true if ann conforms to this type and false if it doesn't.
     */
    boolean instanceOf(Annotation ann);

    /**
     * <p>See if an AnnotationType is a specialisation of this type.</p>
     *
     * <p>An AnnotationType is a sub-type if it restricts each of the
     * properties of the super-type to a type that can be cast to the
     * type in the super-type. Note that this is not always a cast in
     * the pure Java sense; it may include checks on the number and
     * type of members in collections or other criteria.</p>
     *
     * @param subType an AnnotationType to check.
     * @return true if subType is a sub-type of this type.
     *
     * 
     */
    boolean subTypeOf(AnnotationType subType);

    /**
     * <p>Retrieve the constraint that will be applied to all
     * properties with a given key.</p>
     *
     * <p>For an <code>Annotation</code> to be accepted, each key in
     * getProperties() must be present in the annotation and each of the
     * values associated with those properties must match the
     * constraint.</p>
     * If you want to find out exactly what constraints will be
     * applied to a particular propery key
     *
     * @param key the property to be validated.
     * @return PropertyConstraint the constraint by which the values
     * must be accepted.
     *
     */
    CollectionConstraint getConstraint(Object key);

    /**
     * Set the constraints associated with a property.  This method constrains
     * the value of the specified property such that all members must match
     * <code>con</code>, and the number of members must match <code>card</code>.
     * It implicitly constructs a <code>CollectionConstraint.AllValuesIn</code>
     * instance.
     * <p>When you are building your own AnnotationType</p>
     * 
     * @param key  the name of the property to constrain
     * @param con  the PropertyConstraint to enforce
     * @param card the CardinalityConstraint to enforce
     *
     * 
     */
    void setConstraints(
      Object key,
      PropertyConstraint con,
      Location card
    );

    /**
     * Specifies the constraint to apply to the specified property.
     *
     * @param key    the name of the property to constrain
     * @param con    the constraint to apply to this slot.
     *
     */
    void setConstraint(
        Object key,
        CollectionConstraint con
    );

    /**
     * Set the constraints that will apply to all properties without an
     * explicitly defined set of constraints.  This method constrains
     * the value of the specified property such that all members must match
     * <code>con</code>, and the number of members must match <code>card</code>.
     * It implicitly constructs a <code>CollectionConstraint.AllValuesIn</code>
     * instance.
     *
     * @param pc  the default PropertyConstraint
     * @param cc the default CardinalityConstraint
     *
     */
    void setDefaultConstraints(PropertyConstraint pc, Location cc);

    /**
     * Specifies the default constraint to apply to properties where no
     * other constraint is specified.
     *
     * @param cc The default constraint.
     */

    void setDefaultConstraint(CollectionConstraint cc);

    /**
     * Get the CollectionConstraint that will be applied to all properties without
     * an explicit binding. This defaults to CollectionConstraint.ALL.
     * <p>If you want to find out exactly what constraint will be
     * applied to properties with no explicitly defined constraints</p>
     *
     * @return the default CollectionConstraint
     *
     */
    CollectionConstraint getDefaultConstraint();

    /**
     * Retrieve the set of properties for which constraints have been explicity specified.
     * Discover which properties have explicit constraints
     * 
     * @return the Set of properties to validate.
     * 
     */
    Set getProperties();

    /**
     * Set the property in an annotation bundle according to the type we believe
     * it should be. This will take care of any neccisary packing or unpacking
     * to Collections.
     *
     * @param ann  the Annotation to modify
     * @param property  the property key Object
     * @param value  the property value Object
     * @throws ChangeVetoException  if the value could not be accepted by this
     *         annotation type for that property key, or if the Annotation could
     *         not be modified
     *
     */
    void setProperty(Annotation ann, Object property, Object value)
        throws ChangeVetoException;

    /**
     * Add a value to the specified property slot.
     *
     * @param ann  the Annotation to modify
     * @param property  the property key Object
     * @param value  the property value Object
     * @throws ChangeVetoException  if the value could not be accepted by this
     *         annotation type for that property key, or if the Annotation could
     *         not be modified
     *
     */
    void addProperty(Annotation ann, Object property, Object value)
        throws ChangeVetoException;

    /**
     * Get the Collection of values associated with an Annotation bundle
     * according to the type we believe it to be. This will take care of any
     * neccisary packing or unpacking to Collections. Properties with no values
     * will return empty Collections.
     *
     * @param ann  the Annotation to access
     * @param property  the property key Object
     * @return a Collection of values
     * @throws ChangeVetoException  if the value could not be removed
     *
     */
    Collection getProperty(Annotation ann, Object property)
        throws ChangeVetoException;

    /**
     * Remove a value from the specified property slot.
     *
     * @param ann  the Annotation to modify
     * @param property  the property key Object
     * @param value  the property value Object
     * @throws ChangeVetoException  if the Annotation could
     *         not be modified
     *
     */
    void removeProperty(Annotation ann, Object property, Object value)
        throws ChangeVetoException;

    /**
     * Set the comment for the whole AnnotationType.
     * This is human-readable text.
     *
     * @param comment  the new comment
     */
    void setComment(String comment);
    
    /**
     * Get the comment for the whole AnnotationType.
     * This is human-readable text.
     *
     * @return the comment
     */
    String getComment();
    
    /**
     * Set the comment for a particular property.
     * This is a human-readable description of the property.
     *
     * @param property the property to comment on
     * @param comment  the comment
     */
    void setComment(Object property, String comment);
    
    /**
     * Get the comment for a particular property.
     * This is a human-readable description of the property.
     *
     * @param property the property to get a comment for
     * @return the comment
     */
    String getComment(Object property);
    
    /**
     * <p>An abstract base class useful for implementing AnnotationType
     * instances.</p>
     *
     * <p>This provides deffinitions for the logical operators (validate(),
     * subTypeOf()), the mutators (setProperty(), getProperty() and
     * deleteProperty()) and toString() that you may not want to
     * write yourself. It leaves the data-related methods up to you.</p>
     *
     * @since 1.3
     * @author Matthew Pocock
     * @author Thomas Down
     *
     */
    abstract class Abstract implements AnnotationType {
        public void setConstraints(Object key, PropertyConstraint pc, Location cc) {
            setConstraint(key, new CollectionConstraint.AllValuesIn(pc, cc));
        }

        public void setDefaultConstraints(PropertyConstraint pc, Location cc) {
            setDefaultConstraint(new CollectionConstraint.AllValuesIn(pc, cc));
        }

        public boolean instanceOf(Annotation ann) {
          Set props = new HashSet();
          props.addAll(getProperties());
          props.addAll(ann.keys());

          for (Iterator i = getProperties().iterator(); i.hasNext();) {
            Object key = i.next();
            CollectionConstraint con = getConstraint(key);
            Object value;
            if (ann.containsProperty(key)) {
                value = ann.getProperty(key);
            } else {
                value = Collections.EMPTY_SET;
            }
            if (!con.accept(value)) {
                return false;
            }
          }

          return true;
        }

        public final void setProperty(
          Annotation ann,
          Object property,
          Object value
        ) throws ChangeVetoException
        {
            CollectionConstraint cons = getConstraint(property);
            if (cons.accept(value)) {
                ann.setProperty(property, value);
            } else {
                throw new ChangeVetoException("Setting property " + property + " to " + value + " would violate constraints in " + cons);
            }
        }

        public final Collection getProperty(Annotation ann, Object property)
        throws ChangeVetoException {
          Collection vals = null;

          if(!ann.containsProperty(property)) {
            vals = Collections.EMPTY_SET;
          } else {
            Object val = ann.getProperty(property);
            if(val instanceof Collection) {
              vals = (Collection) val;
            } else {
              vals = Collections.singleton(val);
            }
          }

          return vals;
        }

        public final void addProperty(
          Annotation ann,
          Object key,
          Object value
        ) throws ChangeVetoException
        {
            CollectionConstraint cons = getConstraint(key);
            Collection oldValue;
            if (ann.containsProperty(key)) {
                Object ov = ann.getProperty(key);
                if (ov instanceof Collection) {
                    oldValue = (Collection) ov;
                } else {
                    oldValue = new ArrayList();
                    oldValue.add(ov);
                }
            } else {
                oldValue = new ArrayList();
            }
            if (cons.validateAddValue(oldValue, value)) {
                oldValue.add(value);
                ann.setProperty(key, oldValue);
            } else {
                throw new ChangeVetoException("Adding value " + value + " to " + key + " would violate constraints");
            }
        }

        public final void removeProperty(
          Annotation ann,
          Object key,
          Object value
        ) throws ChangeVetoException
        {
            CollectionConstraint cons = getConstraint(key);
            Collection oldValue;
            if (ann.containsProperty(key)) {
                Object ov = ann.getProperty(key);
                if (ov instanceof Collection) {
                    oldValue = (Collection) ov;
                } else {
                    oldValue = new ArrayList();
                    oldValue.add(ov);
                }
            } else {
                oldValue = new ArrayList();
            }
            if (!oldValue.contains(value)) {
                throw new ChangeVetoException("Property " + key + " does not have value " + value);
            } else {
                if (cons.validateRemoveValue(oldValue, value)) {
                    oldValue.remove(value);
                    if (oldValue.size() > 0) {
                        ann.setProperty(key, oldValue);
                    } else {
                        ann.removeProperty(key);
                    }
                } else {
                    throw new ChangeVetoException("Adding value " + value + " to " + key + " would violate constraints");
                }
            }
        }

        public String toString() {
          StringBuffer sb = new StringBuffer(
            "AnnotationType: " +
            getComment() +
            " {"
          );

          for(Iterator i = getProperties().iterator(); i.hasNext(); ) {
            Object key = i.next();
            CollectionConstraint cc = getConstraint(key);

            sb.append(" [" + ", " + cc + " " + key + "]");
          }
          sb.append(" [*, " + getDefaultConstraint() + "]");

          sb.append(" }");

          return sb.toString();
        }

        public boolean subTypeOf(AnnotationType subType) {
            Set props = new HashSet();
            props.addAll(getProperties());
            props.addAll(subType.getProperties());

            for (Iterator i = props.iterator(); i.hasNext();) {
                Object key = i.next();

                CollectionConstraint thisPC = getConstraint(key);
                CollectionConstraint subPC = subType.getConstraint(key);
                if (! thisPC.subConstraintOf(subPC)) {
                    return false;
                }
            }

            return true;
        }
    }

    /**
     * <p>An implementation of <code>AnnotationType</code>.</p>
     *
     * <p>To build an instance of <code>AnnotationType.Impl</code>,
     * first invoke the no-args constructor, and then use the
     * setPropertyConstraint method to build the property->constraint
     * mapping.</p>
     * A convenient class for when you need an AnnotationType
     * instance and don't want to write your own
     *
     * @since 1.3
     * @author Matthew Pocock
     * @author Thomas Down
     *
     */
    class Impl extends AnnotationType.Abstract {
        private Map cons;
        private CollectionConstraint unknown;
        private String comment;
        private Map comments;

        /**
         * Create a new Impl with no constraints.
         */
        public Impl() {
            cons = new SmallMap();
            unknown = CollectionConstraint.ANY;
            comment = "";
            comments = new SmallMap();
        }

        /**
         * Create a new Impl with a default property and cardinality constraint.
         *
         * @param defaultPC  the default PropertyConstraint
         * @param defaultCC  the default CardinalityConstraint
         */
        public Impl(PropertyConstraint defaultPC, Location defaultCC) {
            this();
            this.setDefaultConstraints(defaultPC, defaultCC);
        }

        /**
         * Create a new Impl with a default collection constraint.
         *
         * @param unknown  the default CollectionConstraint
         */
        public Impl(CollectionConstraint unknown) {
            this();
            setDefaultConstraint(unknown);
        }

        public void setDefaultConstraint(CollectionConstraint cc) {
            this.unknown = cc;
        }

        public CollectionConstraint getDefaultConstraint() {
            return unknown;
        }

        public CollectionConstraint getConstraint(Object key) {
            CollectionConstraint pc = (CollectionConstraint) cons.get(key);
            if (pc == null) {
                pc = unknown;
            }
            return pc;
        }

        public void setConstraint(
          Object key,
          CollectionConstraint cc
        ) {
            cons.put(key, cc);
        }

        public Set getProperties() {
            return cons.keySet();
        }
        
        public void setComment(String comment) {
          this.comment = comment;
        }
        
        public String getComment() {
          return comment;
        }
        
        public void setComment(Object key, String comment) {
          comments.put(key, comment);
        }
        
        public String getComment(Object key) {
          return (String) comments.get(key);
        }
    }
}
