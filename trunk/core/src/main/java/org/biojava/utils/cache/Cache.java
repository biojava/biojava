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

/**
 * Interface for managing caches of objects.  Caches serve
 * as factories for CacheReference objects.  These references
 * can be cleared at any time after creation, using heuristics
 * determined by the Cache implementation.
 *
 * @author Thomas Down
 * @since 1.1
 */

public interface Cache {
    /**
     * Construct a temporary reference to an object.  The reference
     * persists until it becomes dereferenced itself, it is explicitly
     * cleared by the user, or the cache determines that it is a
     * candidate for disposal.
     */

    public CacheReference makeReference(Object o);
}
