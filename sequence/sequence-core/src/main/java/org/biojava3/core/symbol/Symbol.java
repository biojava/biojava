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
package org.biojava3.core.symbol;

import java.io.ObjectStreamException;
import java.io.Serializable;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Set;

/**
 * The symbol class encapsulates the concept of a unique entity, such as a DNA
 * base. It can contain a reference to any object at all. The class of the object
 * inside the symbol determines whether or not the symbol is cacheable. If it is
 * cacheable, then the symbol behaves as a singleton. By default symbols based around
 * Character and Boolean objects are cacheable. If it is not cacheable, then
 * the symbol is not a singleton and may be regenerated on demand, depending on how 
 * recently it has been used. Symbols created on demand are kept in a LRU cache
 * so that frequently requested symbols are reused to save time and space.
 * @author Richard Holland
 * @since 3.0
 * @param T the type of object to hold within the symbol.
 */
public final class Symbol implements Serializable, Comparable<Symbol> {
	
	private static final long serialVersionUID = 1L;

    @SuppressWarnings("unchecked")
	private final Comparable obj;
    private static final Map<Object, Symbol> cache = new HashMap<Object, Symbol>();
    private static int MAX_LRU_SIZE = 100;
    private static final LinkedHashMap<Object, Symbol> lru = new LinkedHashMap<Object, Symbol>() {					

    	private static final long serialVersionUID = 1L;

        @Override
        protected boolean removeEldestEntry(Map.Entry<Object, Symbol> eldest) {
            return this.size() > MAX_LRU_SIZE;
        }
    };
    @SuppressWarnings("unchecked")
	private static final Set<Class<? extends Comparable>> cacheableClasses = new HashSet<Class<? extends Comparable>>();

    static {
        cacheableClasses.add(Character.class);
        cacheableClasses.add(Boolean.class);
    }

    /**
     * Sets the maximum number of non-singleton symbols to keep in the cache.
     * @param size the size of the cache.
     */
    public static void setMaxLRUSize(int size) {
        MAX_LRU_SIZE = size;
    }

    /**
     * Gets the maximum number of non-singleton symbols to keep in the cache.
     * @return the size of the cache.
     */
    public static int getMaxLRUSize() {
        return MAX_LRU_SIZE;
    }

    /**
     * Indicates whether symbols wrapping the given object type are singletons.
     * @param clazz the object type the symbols wrap.
     * @param cacheable {@code true} if they should be cached as singletons,
     * {@code false} if they should be kept in the LRU cache and recreated on
     * demand instead.
     */
    @SuppressWarnings("unchecked")
	public static void setCacheable(Class<? extends Comparable> clazz, boolean cacheable) {
        boolean alreadyCacheable = isCacheable(clazz);
        if (cacheable && !alreadyCacheable) {
            cacheableClasses.add(clazz);
        } else if (!cacheable && alreadyCacheable) {
            cacheableClasses.remove(clazz);
            for (Iterator<Map.Entry<Object, Symbol>> i = cache.entrySet().iterator(); i.hasNext();) {
                Map.Entry<Object, Symbol> entry = i.next();
                if (clazz == entry.getKey().getClass()) {
                    i.remove();
                }
            }
        }
    }

    /**
     * Checkes whether symbols wrapping the given object type are singletons.
     * @param clazz the object type the symbols wrap.
     * @return {@code true} if they should be cached as singletons,
     * {@code false} if they should be kept in the LRU cache and recreated on
     * demand instead.
     */
    @SuppressWarnings("unchecked")
	public static boolean isCacheable(Class<? extends Comparable> clazz) {
        return cacheableClasses.contains(clazz);
    }

    /**
     * Obtain a symbol wrapping the given object. If the object type is
     * cacheable, the symbol returned will be a singleton. Otherwise it will
     * be recreated on demand if required.
     * @param obj the object the symbol should wrap.
     * @return the symbol.
     */
	@SuppressWarnings("unchecked")
	public static Symbol get(Comparable obj) {
        if (obj == null) {
            throw new NullPointerException("Symbols cannot be made from null objects.");
        }
        Class<? extends Comparable> clazz = obj.getClass();
        if (isCacheable(clazz)) {
            if (!cache.containsKey(obj)) {
                cache.put(obj, new Symbol(obj));
            }
            return cache.get(obj);
        } else {
            if (!lru.containsKey(obj)) {
                lru.put(obj, new Symbol(obj));
            }
            return lru.get(obj);
        }
    }

    /**
     * Construct a new symbol around the given object.
     * @param obj the object.
     */
    @SuppressWarnings("unchecked")
	private Symbol(Comparable obj) {
        this.obj = obj;
    }

    /**
     * Get the object this symbol contains.
     * @return the object.
     */
    @SuppressWarnings("unchecked")
	public Comparable getObject() {
        return this.obj;
    }

    /**
     * Check to see if another symbol has the same object type as this one.
     * @param sym the other symbol.
     * @return {@code true} if the two symbols share the same object type,
     * based on class equality.
     */
    public boolean hasSameObjectType(Symbol sym) {
        return this.obj.getClass() == sym.obj.getClass();
    }

    @Override
    public boolean equals(Object obj) {
        if (obj == null) {
            return false;
        }
        if (this.getClass() != obj.getClass()) {
            return false;
        }
        final Symbol other = (Symbol) obj;
        return this.obj.equals(other.obj);
    }

    @Override
    public int hashCode() {
        return this.obj.hashCode();
    }
    
    @SuppressWarnings("unchecked")
	public int compareTo(Symbol o) {
        return this.obj.compareTo(o.obj);
    }

    @Override
    public String toString() {
        return this.obj.toString();
    }

    /**
     * Make sure that the symbol remains a singleton on deserialization.
     * @return the singleton symbol.
     * @throws ObjectStreamException should never happen.
     */
    protected Object readResolve() throws ObjectStreamException {
        return get(this.obj);
    }
}
