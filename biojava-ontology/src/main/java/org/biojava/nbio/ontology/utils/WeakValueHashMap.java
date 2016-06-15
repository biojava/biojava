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

import java.lang.ref.PhantomReference;
import java.lang.ref.Reference;
import java.lang.ref.ReferenceQueue;
import java.util.*;

/**
 * Map implementation which keeps weak references to values.
 * Entries are removed from the map when their value is
 * no longer reachable using normal (hard) references.  This is
 * useful for maintaining canonical copies of objects without forcing
 * these objects to remain in memory forever.
 *
 * <p>
 * Note that this is distinct from the standard library class,
 * <code>WeakHashMap</code> which has weak <em>keys</em>.
 * </p>
 *
 * @author Thomas Down
 * @since 1.3
 */

public class WeakValueHashMap extends AbstractMap {
	private final Map keyToRefMap;
	private final ReferenceQueue queue;
	private final Set iteratorRefs;
	private final ReferenceQueue iteratorRefQueue;

	public WeakValueHashMap() {
	keyToRefMap = new HashMap();
	queue = new ReferenceQueue();
	iteratorRefs = new HashSet();
	iteratorRefQueue = new ReferenceQueue();
	}

	private void diddleReferenceQueue() {
	// Avoid making behind-the-scenes modifications while iterators exist.

	if (iteratorRefs.size() > 0) {
		Reference ref;
		while ((ref = iteratorRefQueue.poll()) != null) {
		iteratorRefs.remove(ref);
		}
		if (iteratorRefs.size() > 0) {
		return;
		}
	}

	KeyedWeakReference ref;
	while ((ref = (KeyedWeakReference) queue.poll()) != null) {
		keyToRefMap.remove(ref.getKey());
	}
	}

	@Override
	public Object put(Object key, Object value) {
	diddleReferenceQueue();
	Reference oldRef = (Reference) keyToRefMap.put(key, new KeyedWeakReference(key, value, queue));
	if (oldRef != null) {
		Object oldRefVal = oldRef.get();
		oldRef.clear();
		return oldRefVal;
	} else {
		return null;
	}
	}

	@Override
	public Object get(Object key) {
	Reference ref = (Reference) keyToRefMap.get(key);
	if (ref != null) {
		return ref.get();
	} else {
		return null;
	}
	}

	@Override
	public boolean containsKey(Object o) {
	diddleReferenceQueue();
	return keyToRefMap.containsKey(o);
	}

	@Override
	public Set entrySet() {
	diddleReferenceQueue();
	return new WVEntrySet();
	}

	private class WVEntrySet extends AbstractSet {
	private Set keyRefEntrySet;

	public WVEntrySet() {
		super();
		keyRefEntrySet = keyToRefMap.entrySet();
	}

	@Override
	public int size() {
		return keyRefEntrySet.size();
	}

	@Override
	public Iterator iterator() {
		Iterator i = new WVEntryIterator(keyRefEntrySet.iterator());
		iteratorRefs.add(new PhantomReference(i, iteratorRefQueue));
		return i;
	}
	}

	private class WVEntryIterator implements Iterator {
	private Object cache;
	private Iterator keyRefIterator;

	public WVEntryIterator(Iterator keyRefIterator) {
		this.keyRefIterator = keyRefIterator;
	}

	@Override
	public boolean hasNext() {
		if (cache == null) {
		primeCache();
		}
		return cache != null;
	}

	@Override
	public Object next() {
		if (cache == null) {
		primeCache();
		}
		if (cache == null) {
		throw new NoSuchElementException();
		} else {
		Object o = cache;
		cache = null;
		return o;
		}
	}

	@Override
	public void remove() {
		if (cache != null) {
		throw new IllegalStateException("next() not called");
		} else {
		keyRefIterator.remove();
		}
	}

	private void primeCache() {
		while (keyRefIterator.hasNext()) {
		Map.Entry krme = (Map.Entry) keyRefIterator.next();
		Object ref = ((Reference) krme.getValue()).get();
		if (ref != null) {
			cache = new WVMapEntry(krme.getKey(), ref);
			return;
		}
		}
	}
	}

	private static class WVMapEntry implements Map.Entry {
	private Object key;
	private Object value;

	private WVMapEntry(Object key, Object value) {
		this.key = key;
		this.value = value;
	}

	@Override
	public Object getKey() {
		return key;
	}

	@Override
	public Object getValue() {
		return value;
	}

	@Override
	public Object setValue(Object v) {
		throw new UnsupportedOperationException();
	}

	@Override
	public boolean equals(Object o) {
		if (! (o instanceof Map.Entry)) {
		return false;
		}

		Map.Entry mo = (Map.Entry) o;
		return ((key == null ? mo.getKey() == null : key.equals(mo.getKey())) &&
			(value == null ? mo.getValue() == null : value.equals(mo.getValue())));
	}

	@Override
	public int hashCode() {
		return (key == null ? 0 : key.hashCode()) ^ (value == null ? 0 : value.hashCode());
	}
	}
}
