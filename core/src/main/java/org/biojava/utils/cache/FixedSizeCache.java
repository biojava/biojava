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

package org.biojava.utils.cache;

import java.lang.ref.WeakReference;
import java.util.LinkedList;
import java.util.List;

/**
 * Cache which stores up to <code>limit</code> Objects.
 *
 * @since 1.1
 * @author Thomas Down
 */

public class FixedSizeCache implements Cache {
    private List objects;
    private int sizeLimit;

    {
	objects = new LinkedList();
    }

    public FixedSizeCache(int limit) {
	sizeLimit = limit;
    }

    public CacheReference makeReference(Object o) {
	CacheReference cr = new FixedSizeCacheReference(o);
	objects.add(new WeakReference(cr));
	while (objects.size() > sizeLimit) {
	    CacheReference old = (CacheReference) ((WeakReference) objects.remove(0)).get();
	    if (old != null) {
		old.clear();
	    }
	}
	
	return cr;
    }

    public int getLimit() {
	return sizeLimit;
    }

    public void setLimit(int limit) {
	this.sizeLimit = limit;
    }

    private class FixedSizeCacheReference implements CacheReference {
	private Object o;

	private FixedSizeCacheReference(Object o) {
	    this.o = o;
	}

	public Object get() {
	    return o;
	}

	public void clear() {
	    o = null;
	}
    }
}

