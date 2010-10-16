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
 * Interface for a reference to an object, analogous to
 * <code>java.lang.ref.Referencce</code>, but more flexible.
 * CacheReference implementations are always obtained by calling
 * <code>makeReference</code> on a <code>Cache</code> object.
 *
 * @author Thomas Down
 * @since 1.1
 */

public interface CacheReference {
    public Object get();
    public void clear();
}

