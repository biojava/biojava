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

import java.lang.ref.Reference;
import java.lang.ref.SoftReference;

/**
 * Cache which is cleared according to memory pressure.  This
 * is simply a wrapper around java.lang.ref.SoftReference, and
 * the performance will depend on the behaviour of SoftReference
 * on your platform.
 *
 * @author Thomas Down
 * @since 1.1
 */

public class SoftReferenceCache implements Cache {
    public CacheReference makeReference(Object o) {
	    return new ReferenceReference(new SoftReference(o));
    }

    private static class ReferenceReference implements CacheReference {
      private Reference ref;
      
      ReferenceReference(Reference r) {
        ref = r;
      }
      
      public Object get() {
        if (ref != null)
          return ref.get();
        return null;
      }
      
      public void clear() {
        ref = null;
      }
    }
}
