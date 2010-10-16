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

package org.biojava.utils;

import java.util.AbstractMap;
import java.util.AbstractSet;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.Set;

/**
 * <p>
 * Overlap one map onto another. This allows you to have a map with local values
 * and default values. The local and default values are provided by a child and
 * parent map.
 * </p>
 *
 * @author Thomas Down
 * @author Matthew Pocock
 */

public class OverlayMap extends AbstractMap {
    private Map parent;
    private Map overlay;

    /**
     * Build a new map with default key-value pairs.
     *
     * @param parent  the default fall-through Map
     * @param overlay the overriding Map
     */
    public OverlayMap(Map parent, Map overlay) {
	super();
	this.parent = parent;
	this.overlay = overlay;
    }

    /**
     * Build a new map with default key-value pairs.
     *
     * @param parent  the default fall-through Map
     */
    public OverlayMap(Map parent) {
	super();
	this.parent = parent;
	this.overlay = new HashMap();
    }

    /**
     * Return the object containing the fallback mappings.
     * This is the actual parent map, not a copy.
     *
     * @return the parent map
     */

    public Map getParentMap() {
	return parent;
    }

    /**
     * Return the object containing the overlay mappings.
     * This is the actual child map, not a copy.
     *
     * @return the child map
     */

    public Map getOverlayMap() {
	return overlay;
    }

    //
    // Basic map operations, done explicitly
    //

    public Object get(Object key) {
	Object value = overlay.get(key);
	if (value == null)
	    value = parent.get(key);
	return value;
    }

    public Set entrySet() {
	return new OEntrySet();
    }

    public Set keySet() {
	return new OKeySet();
    }

    public boolean containsKey(Object key) {
	return overlay.containsKey(key) || parent.containsKey(key);
    }

    public Object put(Object key, Object value) {
	Object old = get(key);
	overlay.put(key, value);
	return old;
    }

    private class OKeySet extends AbstractSet {
	private Set parentKeys;
	
	private OKeySet() {
	    super();
	    parentKeys = parent.keySet();
	}
     
	public Iterator iterator() {
	    return new Iterator() {
		    Iterator oi = overlay.keySet().iterator();
		    Iterator pi = parentKeys.iterator();
		    Object peek = null;
		    
		    public boolean hasNext() {
			if (peek == null)
			    peek = nextObject();
			return (peek != null);
		    }
         
		    public Object next() {
			if (peek == null) {
			    peek = nextObject();
			}
			if (peek == null) {
			    throw new NoSuchElementException();
			}
			Object o = peek;
			peek = null;
			return o;
		    }
         
		    private Object nextObject() {
			if (oi.hasNext()) {
			    return oi.next();
			}
			Object po = null;
			while (po == null && pi.hasNext()) {
			    po = pi.next();
			    if (overlay.containsKey(po)) {
				po = null;
			    }
			}
			return po;
		    }
		    
		    public void remove() {
			throw new UnsupportedOperationException();
		    }
		};
	}
     
	public int size() {
	    int i = 0;
	    Iterator keys = iterator();
	    while(keys.hasNext()) {
		keys.next();
		++i;
	    }
	    return i;
	}
     
	public boolean contains(Object o) {
	    return overlay.containsKey(o) || parentKeys.contains(o);
	}
    }

    private class OEntrySet extends AbstractSet {
	OKeySet ks;
	
	private OEntrySet() {
	    super();
	    ks = new OKeySet();
	}
    
	public Iterator iterator() {
	    return new Iterator() {
		Iterator ksi = ks.iterator();
		    
		public boolean hasNext() {
		    return ksi.hasNext();
		}
		    
		public Object next() {
		    Object k = ksi.next();
		    Object v = get(k);
		    return new OMapEntry(k, v);
		}
        
		public void remove() {
		    throw new UnsupportedOperationException();
		}
	    };
	}
    
	public int size() {
	    return ks.size();
	}
    }
    
    private static class OMapEntry implements Map.Entry {
	private Object key;
	private Object value;
    
	private OMapEntry(Object key, Object value) {
	    this.key = key;
	    this.value = value;
	}
    
	public Object getKey() {
	    return key;
	}
    
	public Object getValue() {
	    return value;
	}
    
	public Object setValue(Object v) {
	    throw new UnsupportedOperationException();
	}
    
	public boolean equals(Object o) {
	    if (! (o instanceof Map.Entry)) {
		return false;
	    }
      
	    Map.Entry mo = (Map.Entry) o;
	    return ((key == null ? mo.getKey() == null : key.equals(mo.getKey())) &&
		    (value == null ? mo.getValue() == null : value.equals(mo.getValue())));
	}
    
	public int hashCode() {
	    return (key == null ? 0 : key.hashCode()) ^ (value == null ? 0 : value.hashCode());
	}
    }
}
