/*
 *                  BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *   http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *   http://www.biojava.org
 */

package org.biojava3.ontology.utils;

import java.io.Serializable;
import java.lang.reflect.Field;
import java.lang.reflect.Modifier;
import java.util.HashSet;
import java.util.Iterator;
import java.util.NoSuchElementException;
import java.util.Set;


/**
 *
 * Class for all constants which are used to indicate change
 * types.  Note that all ChangeType objects must be accessible
 * via a public static field of some class or interface.  These should
 * be specified at construction time, so that the ChangeType can
 * be properly serialized.  Typically, they should be constructed
 * using code like:
 *
 * <pre>
 * class MyClassWhichCanFireChangeEvents {
 *   public final static ChangeType CHANGE_COLOR = new ChangeType(
 *              "Color change",
 *              MyClassWhichCanFireChangeEvents.class,
 *              "CHANGE_COLOR");
 *   // Rest of the class here...
 * }
 * </pre>
 *
 * <p>
 * As of BioJava 1.2, the known ChangeTypes of a system follow a simple
 * hierarchy with single inheritance.  All ChangeTypes
 * (except ChangeType.UNKNOWN) have a parent ChangeType (defaulting
 * to ChangeType.UNKNOWN).  Generally, when a listener is registered
 * for changetype <code>foo</code>, changes of type <code>bar</code>
 * should be accepted if <code>bar</code> is a sub-type of <code>foo</code>.
 * This can be checked using an expression like:
 * </p>
 *
 * <pre>
 *     bar.isMatchingType(foo);
 * </pre>
 *
 * @author     Thomas Down
 * @author     Matthew Pocock
 * @since      1.1
 */

public final class ChangeType implements Serializable {
    private final String name;
    private final Field ourField;
    private final ChangeType superType;

  /**
   * Constant ChangeType field which indicates that a change has
   * occured which can't otherwise be represented.  Please do not
   * use this when there is another, more sensible, option. This
   * is the fallback for when you realy don't know what else to
   * do.
   *
   * <p>
   * As of BioJava 1.2, this type is the root of the ChangeType
   * hierarchy.  Listening for this type is equivalent to listening
   * for all ChangeTypes.
   * </p>
   */
  public static final ChangeType UNKNOWN;

    /**
     *  Construct a new ChangeType.
     *
     * @param  name      The name of this change.
     * @param  ourField  The public static field which contains this
     *                   ChangeType.
     * @param  superType The supertype of this type.
     *
     * @since 1.2
     */

  public ChangeType(String name, Field ourField, ChangeType superType) {
    this.name = name;
    this.ourField = ourField;
    this.superType = superType;
  }

  /**
   *  Construct a new ChangeType with superType UNKNOWN.
   *
   * @param  name      The name of this change.
   * @param  ourField  The public static field which contains this
   *      ChangeType.
   */
  public ChangeType(String name, Field ourField) {
      this(name, ourField, ChangeType.UNKNOWN);
  }

  /**
   *  Construct a new ChangeType with supertype UNKNOWN.
   *
   * @param  name   The name of this change.
   * @param  clazz  The class which is going to contain this change.
   * @param  fname
   * The name of the field in <code>clazz</code> which
   * is to contain a reference to this change.
   * @throws        BioError If the field cannot be found.
   */
  public ChangeType(String name, Class clazz, String fname) {
      this(name, clazz, fname, ChangeType.UNKNOWN);
  }

  /**
   *  Construct a new ChangeType.
   *
   * @param  name   The name of this change.
   * @param  clazz  The class which is going to contain this change.
   * @param  fname
   * The name of the field in <code>clazz</code> which
   * is to contain a reference to this change.
   * @param superType the supertype of this type.
   * @throws        BioError If the field cannot be found.
   *
   * @since 1.2
   */
  public ChangeType(String name, Class clazz, String fname, ChangeType superType) {
    this.name = name;
    this.superType = superType;
    try {
      this.ourField = clazz.getField(fname);
    }
    catch (Exception ex) {
      throw new AssertionFailure("Couldn't find field " + fname + " in class " + clazz.getName(), ex);
    }
  }

    public ChangeType(String name, String className, String fieldName, ChangeType superType) {
	this.name = name;
	this.superType = superType;
	try {
	    Class clazz = Class.forName(className);
	    this.ourField = clazz.getField(fieldName);
	} catch (Exception ex) {
	    throw new AssertionFailure(
				  "Couldn't find class or field " + className +
				  "->" + fieldName,
          ex
				  );
	}
    }

    public ChangeType(String name, String className, String fieldName) {
	this(name, className, fieldName, ChangeType.UNKNOWN);
    }

  /**
   *  Return the name of this change.
   *
   * @return    The Name value
   */
  public String getName() {
    return name;
  }

    /**
     * Return a Field object where this change type is declared.
     */

    public Field getField() {
	return ourField;
    }

  /**
   *  Return a string representation of this change.
   *
   * @return    Description of the Returned Value
   */
  public String toString() {
    return "ChangeType: " + name;
  }

  /**
   *  Make a placeholder for this object in a serialized stream.
   *
   * @return    Description of the Returned Value
   */
  private Object writeReplace() {
    return new StaticMemberPlaceHolder(ourField);
  }

  static {
    UNKNOWN = new ChangeType(
      "Unknown change",
      ChangeType.class,
      "UNKNOWN",
      null
    );
  }

    /**
     * Get all ChangeType objects defined within a class.  This
     * includes ChangeTypes defined in superclasses and interfaces.
     * Only fields declared as public [final] static ChangeType are
     * returned.
     *
     * @param clazz A class to introspect
     */

    public static Set getChangeTypes(Class clazz)
    {
	Set types = new HashSet();
	Field[] fields = clazz.getFields();
	for (int i = 0; i < fields.length; ++i) {
	    Field f = fields[i];
	    if (f.getType().equals(ChangeType.class) && (f.getModifiers() & Modifier.STATIC) != 0) {
		try {
		    types.add(f.get(null));
		} catch (Exception ex) {}
	    }
	}

	return types;
    }

    /**
     * Return the immediate supertype (internal use only)
     */

    private ChangeType getSuperType() {
	return superType;
    }

    /**
     * Return an iterator which contains this type, and all supertypes.
     *
     * @since 1.2
     */

    public Iterator matchingTypes() {
	return new Iterator() {
		ChangeType cti = ChangeType.this;

		public boolean hasNext() {
		    return cti != null;
		}

		public Object next() {
		    if (cti == null) {
			throw new NoSuchElementException("No more elements");
		    }

		    ChangeType rt = cti;
		    cti = cti.getSuperType();

		    return rt;
		}

		public void remove() {
		    throw new UnsupportedOperationException("Can't remove");
		}
	    } ;
    }

    /**
     * Return <code>true</code> iff <code>ct</code> is equal to this type
     * or any of it's supertypes (including ChangeType.UNKNOWN). If this is
     * true, then ct is more general than this.
     *
     * @since 1.2
     */

    public boolean isMatchingType(ChangeType ct) {
	for (Iterator i = matchingTypes(); i.hasNext(); ) {
	    ChangeType mt = (ChangeType) i.next();
	    if (mt == ct) {
		return true;
	    }
	}

	return false;
    }
}

