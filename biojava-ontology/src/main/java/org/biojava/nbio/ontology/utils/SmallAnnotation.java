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



/**
 * Annotation that is optimized for memory usage.  Access time
 * is linear, so SmallAnnotations are not recommended when
 * the number of entries is large.  However, they are fine for
 * small numbers of keys.
 *
 * @author Thomas Down
 * @author Matthew Pocock
 * @since 1.2
 *
 *
 * A minimal-memory alternative to SimpleAnnotation
 *
 *
 * When creating a large number of small Annotation instances, it is worth
 * instantiating SmallAnnotation. Small is anything up to at least 30 properties
 * but will vary with the JavaVM and underlying platform.
 */

public class SmallAnnotation extends AbstractAnnotation {
	private Map properties;

	@Override
	protected final Map getProperties() {
		if(!propertiesAllocated()) {
			properties = new SmallMap();
		}
		return properties;
	}

	@Override
	protected final boolean propertiesAllocated() {
		return properties != null;
	}

	/**
	 * Return a new SmallAnnotation optimised for small sets of properties.
	 */
	public SmallAnnotation() {
		super();
	}

	/**
	 * Return a new SmallAnnotation that copies all values from another annoation.
	 *
	 * @param ann  the Annoation to copy all values from
	 * @throws NullPointerException if ann is null
	 */
	public SmallAnnotation(Annotation ann) {
		super(ann);
	}

	/**
	 * Return a new SmallAnnotation that copies all values from a Map.
	 *
	 * @param map  the Map to copy values from
	 */
	public SmallAnnotation(Map map) {
		super(map);
	}
}

