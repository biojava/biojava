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

package org.biojava.bio.program.tagvalue;

import java.util.Map;

import org.biojava.utils.SmallMap;

/**
 * <p>
 * <code>TagMapper</code> maps arbitrary object keys to new keys.
 * </p>
 *
 * <p>
 * If there is no explicit mapping from old to new keys, then the old key will
 * be used.
 * </p>
 *
 * @since 1.2
 * @author Matthew Pocock
 * @author Keith James (docs).
 */
public class TagMapper implements PropertyChanger {
    private Map tags;

    /**
     * Creates a new, empty <code>TagMapper</code>.
     */
    public TagMapper() {
        this.tags = new SmallMap();
    }

    /**
     * <code>setNewTag</code>.
     *
     * @param oldTag an <code>Object</code> tag to be substituted.
     * @param newTag an <code>Object</code> tag to substitue for the
     * old value.
     */
    public void setNewTag(Object oldTag, Object newTag) {
        tags.put(oldTag, newTag);
    }

    public Object getNewTag(Object oldTag) {
      Object newTag = tags.get(oldTag);
      if(newTag == null) {
        newTag = oldTag;
      }
      return newTag;
    }
}

