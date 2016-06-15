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

import java.io.InvalidObjectException;
import java.io.ObjectStreamException;
import java.io.Serializable;
import java.lang.reflect.Field;

/**
 * @author Matthew Pocock
 */
public class StaticMemberPlaceHolder implements Serializable {
	private String className;
	private String fieldName;

	public StaticMemberPlaceHolder(Field field) {
		this.className = field.getDeclaringClass().getName();
		this.fieldName = field.getName();
	}

	protected StaticMemberPlaceHolder() {}

	public Object readResolve() throws ObjectStreamException {
		try {
			Class c = Class.forName(className);
			Field f = c.getDeclaredField(fieldName);
			return f.get(null);
		} catch (Exception e) {
			throw new InvalidObjectException(
				"Unable to retrieve static field " + fieldName +
				"for class " + className + " because:\n" +
				e.getMessage()
			);
		}
	}
}
