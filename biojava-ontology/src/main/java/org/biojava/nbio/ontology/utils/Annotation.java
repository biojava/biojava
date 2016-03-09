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


package org.biojava.nbio.ontology.utils;

import java.util.Map;
import java.util.NoSuchElementException;
import java.util.Set;


/**
 * <p>
 * Arbitrary annotation associated with one or more objects.
 * </p>
 *
 * <p>
 * Biological information often does not fit design patterns very well, and can
 * be a jumble of facts and relationships. Annotation objects provide a standard
 * way for you to store this mess as a property of an object.
 * </p>
 *
 * <p>
 * Annotations may contain keys that have Annotations as values. In this way,
 * annotations can be shared among multiple Annotatable objects, and you can
 * represent semi-structured data.
 * </p>
 *
 * <p>
 * It is perfectly possible to wrap up almost any tree-like or flat data
 * structure as Annotation.
 * </p>
 * Other than when using the constructor, you should be able to
 * interact with nearly all Annotation implementations via this API.
 *
 * @author Matthew Pocock
 * @author Thomas Down
 * @see org.biojavax.RichAnnotation
 *
 *
 * @since 1.0
 */
public interface Annotation  {


	/**
	 * <p>
	 * Retrieve the value of a property by key.
	 * </p>
	 *
	 * <p>
	 * Unlike the Map collections, it will complain if the key does not exist. It
	 * will only return null if the key is defined and has value null.
	 * </p> Normal raw access to the property. For cleverer access, use
	 * methods in AnnotationType.
	 *
	 * @param key  the key of the property to retrieve
	 * @return  the object associated with that key
	 * @throws NoSuchElementException if there is no property with the key
	 *
	 *
	 */
	Object getProperty(Object key) throws NoSuchElementException;

	/**
	 * <p>
	 * Set the value of a property.
	 * </p>
	 *
	 * <p>
	 * This method throws an exception if either properties can not be
	 * added to this object, or that this particular property is immutable or
	 * illegal within the implementation.
	 * </p> Normal raw access to the property. For cleverer access, use
	 * methods in AnnotationType.
	 *
	 * @param key the key object
	 * @param value the new value for this key
	 * @throws IllegalArgumentException if the property <code>key</code> is not
	 *         legal
	 * @throws ChangeVetoException if this annotation object can't be changed, or
	 *         if the change was vetoed.
	 */
	void setProperty(Object key, Object value)
			throws IllegalArgumentException;

	/**
	 * Delete a property. Normal raw access to the property. For cleverer access, use
	 * methods in AnnotationType.
	 *
	 * @param key the key object
	 * @throws NoSuchElementException if the property doesn't exist
	 * @throws ChangeVetoException if the change is vetoed
	 * @since 1.3
	 *
	 */

	public void removeProperty(Object key)
			throws NoSuchElementException;

	/**
	 * Returns whether there the property is defined. Normal raw access to the property. For cleverer access, use
	 * methods in AnnotationType.
	 *
	 * @param key the key Object to search for
	 * @return true if this Annotation knows about the key, false otherwise
	 */
	boolean containsProperty(Object key);

	/**
	 * Get a set of key objects.
	 *
	 * @return  a Set of key objects
	 */
	Set keys();

	/**
	 * Return a map that contains the same key/values as this Annotation.
	 * <p>
	 * If the annotation changes, the map may not reflect this.  The Map
	 * may be unmodifiable.
	 *
	 * @return a Map
	 */
	Map asMap();

	/**
	 * <p>
	 * A really useful empty and immutable annotation object.
	 * </p>
	 *
	 * Be careful when stooring Annotation arguments to
	 *  constructors. It is possible that you have been passed EMPTY_ANNOTATION but
	 * that code later on will access this object believing it to be
	 * mutable. For example, the SeqIO factory code clones some
	 * Annotations passed in on Feature.Template instances
	 *
	 * Use this instead of null when you really don't want an object or
	 * an implementation to have annotation even though it should implement
	 * Annotatable.
	 */
	static final Annotation EMPTY_ANNOTATION = new EmptyAnnotation();
}

