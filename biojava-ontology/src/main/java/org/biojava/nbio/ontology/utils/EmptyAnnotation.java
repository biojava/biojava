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

import java.io.NotSerializableException;
import java.io.ObjectStreamException;
import java.io.Serializable;
import java.util.*;

/**
 * An always-empty Annotation.
 *
 * @author Matthew Pocock
 * @author Thomas Down
 *
 * @since 1.0 as part of Annotation
 * @since 1.4 as top-level class
 * @see org.biojavax.EmptyRichAnnotation
 */
class EmptyAnnotation

implements Annotation, Serializable {
	@Override
	public Object getProperty(Object key) throws NoSuchElementException {
		throw new NoSuchElementException(
			"There are no keys in the Empty Annotation object: " +
			key
		);
	}

	@Override
	public void setProperty(Object key, Object value){
		}

	@Override
	public void removeProperty(Object key)

	{

	}

	@Override
	public boolean containsProperty(Object key) {
		return false;
	}

	@Override
	public Set keys() {
		return Collections.EMPTY_SET;
	}

	@Override
	public Map asMap() {
		//return Collections.EMPTY_MAP; 1.3
		return new HashMap();
	}

	private Object writeReplace() throws ObjectStreamException {
		try {
			return new StaticMemberPlaceHolder(Annotation.class.getField("EMPTY_ANNOTATION"));
		} catch (NoSuchFieldException nsfe) {
			throw new NotSerializableException(nsfe.getMessage());
		}
	}

	@Override
	public int hashCode() {
		return asMap().hashCode();
	}

	@Override
	public boolean equals(Object o) {
		if (! (o instanceof Annotation)) {
			return false;
		}

		return ((Annotation) o).asMap().equals(asMap());
	}
}

