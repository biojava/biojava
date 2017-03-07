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

import java.io.Serializable;
import java.util.*;


/**
 * A utility class to ease the problem of implementing an Annotation to that of
 * providing an apropreate implementation of Map. Where possible implementations
 *
 * This class is only intended as a way to implement
 * Annotation. If you are not trying to do that, then don't read on. If you
 * are reading the documentation for an Annotation implementation that extends
 * this, then don't read on. There is nothing to see here.
 *
 * If you are still reading this, then you must be trying to
 * implement Annotation. To do that, extend this class and implement
 * <code>getProperties()</code> and <code>propertiesAllocated()</code>.
 * Where possible implementations should be backed with a
 * <code>LinkedHashMap</code> or similar so properties are iterated in the order
 * they were added.
 *
 * @author Matthew Pocock
 * @author Greg Cox
 *
 * @since 1.0
 */
public abstract class AbstractAnnotation

	implements
		Annotation,
		Serializable
{
	/**
	 *
	 */
	private static final long serialVersionUID = 2753449055959952873L;

/**
	 * Implement this to return the Map delegate. Modifying this return value will
	 * modify the properties associated with this annotation.
	 *
	 * From code in the 1.2 version of AbstractAnnotation
	 * This is required for the implementation of an Annotation that
	 *            extends AbstractAnnotation. Where possible implementations
	 *            should be backed with a
	 *            <code>LinkedHashMap</code> or similar so properties are iterated in the order
	 *            they were added.
	 *
	 * @return a Map containing all properties
	 */
	protected abstract Map getProperties();

	/**
	 * A convenience method to see if we have allocated the properties
	 * Map.
	 * This is required for the implementation of an Annotation that
	 *            extends AbstractAnnotation.
	 * @return true if the properties have been allocated, false otherwise
	 *
	 */
	protected abstract boolean propertiesAllocated();


	@Override
	public Object getProperty(Object key) throws NoSuchElementException {
		if(propertiesAllocated()) {
			Map prop = getProperties();
			if(prop.containsKey(key)) {
				return prop.get(key);
			}
		}
		throw new NoSuchElementException("Property " + key + " unknown");
	}

	@Override
	public void setProperty(Object key, Object value)
	 {

			getProperties().put(key, value);

	}

	@Override
	public void removeProperty(Object key)
		throws  NoSuchElementException
	{
		if (!getProperties().containsKey(key)) {
				throw new NoSuchElementException("Can't remove key " + key.toString());
		}


			getProperties().remove(key);

	}

	@Override
	public boolean containsProperty(Object key) {
		if(propertiesAllocated()) {
			return getProperties().containsKey(key);
		} else {
			return false;
		}
	}

	@Override
	public Set keys() {
		if(propertiesAllocated()) {
			return getProperties().keySet();
		} else {
			return Collections.EMPTY_SET;
		}
	}

	@Override
	public String toString() {
		StringBuffer sb = new StringBuffer("{");
		Map prop = getProperties();
		Iterator i = prop.keySet().iterator();
		if(i.hasNext()) {
			Object key = i.next();
			sb.append(key).append("=").append(prop.get(key));
		}
		while(i.hasNext()) {
			Object key = i.next();
			sb.append(",").append(key).append("=").append(prop.get(key));
		}
		sb.append("}");
		return sb.substring(0);
	}

	@Override
	public Map asMap() {
		return Collections.unmodifiableMap(getProperties());
	}

	/**
	 * Protected no-args constructor intended for sub-classes. This class is
	 * abstract and can not be directly instantiated.
	 */
	protected AbstractAnnotation() {
	}

	/**
	 * Copy-constructor.
	 *
	 * <p>
	 * This does a shallow copy of the annotation. The result is an annotation
	 * with the same properties and values, but which is independant of the
	 * original annotation.
	 * </p>
	 *
	 * @param ann  the Annotation to copy
	 */
	protected AbstractAnnotation(Annotation ann) {
		if(ann == null) {
			throw new NullPointerException(
				"Null annotation not allowed. Use Annotation.EMPTY_ANNOTATION instead."
			);
		}
		if(ann == Annotation.EMPTY_ANNOTATION) {
			return;
		}
		Map properties = getProperties();
		for(Iterator i = ann.keys().iterator(); i.hasNext(); ) {
			Object key = i.next();
			try {
				properties.put(key, ann.getProperty(key));
			} catch (IllegalArgumentException iae) {
				throw new RuntimeException(
					"Property was there and then disappeared: " + key, iae
				);
			}
		}
	}

	/**
	 * Create a new Annotation by copying the key-value pairs from a map. The
	 * resulting Annotation is independant of the map.
	 *
	 * @param annMap  the Map to copy from.
	 */
	public AbstractAnnotation(Map annMap) {
		if(annMap == null) {
			throw new IllegalArgumentException(
				"Null annotation Map not allowed. Use an empy map instead."
			);
		}
		if(annMap.isEmpty()) {
			return;
		}

		Map properties = getProperties();
		for(Iterator i = annMap.keySet().iterator(); i.hasNext(); ) {
			Object key = i.next();
			properties.put(key, annMap.get(key));
		}
	}


	@Override
	public int hashCode() {
		return asMap().hashCode();
	}

	@Override
	public boolean equals(Object o) {
		if(o == this){
				return true;
		}
		if (! (o instanceof Annotation)) {
			return false;
		}

		return ((Annotation) o).asMap().equals(asMap());
	}
}
