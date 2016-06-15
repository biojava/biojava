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

import java.lang.ref.ReferenceQueue;
import java.lang.ref.WeakReference;

/**
 * Subclass of WeakReference which includes and extra field (the key)
 * which can be used to help cleanup once this reference has been
 * enqueued.
 *
 * @author Thomas Down
 * @since 1.3
 */

public class KeyedWeakReference extends WeakReference {
	private Object key;

	public KeyedWeakReference(Object key, Object ref) {
		super(ref);
		this.key = key;
	}

	public KeyedWeakReference(Object key, Object ref, ReferenceQueue queue) {
		super(ref, queue);
		this.key = key;
	}

	public Object getKey() {
		return key;
	}
}
