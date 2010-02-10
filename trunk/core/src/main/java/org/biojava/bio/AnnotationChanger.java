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

package org.biojava.bio;

import java.util.AbstractMap;
import java.util.AbstractSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

import org.biojava.bio.program.tagvalue.ChangeTable;
import org.biojava.utils.ParserException;

/**
 * <code>AnnotationChanger</code> remaps the values of an
 * <code>Annotation</code> to new values specified by a
 * <code>ValueChanger</code>. This will modify the values
 * associated with properties, but not the property names.
 * 
 * For writing light-weigth adaptors to project one type of
 * Annotation to another using a ChangeTable.
 *
 * @since 1.3
 * @author Matthew Pocock
 *
 */
public class AnnotationChanger extends AbstractAnnotation
{
    private final Annotation wrapped;
    private final ChangeTable changer;
    private final Map properties;
  
    /**
     * Creates a new <code>AnnotationChanger</code> using the
     * specified <code>ValueChanger</code> to remap its values.
     *
     * @param wrapped an <code>Annotation</code>.
     * @param changer a <code>ValueChanger</code>.
     */
    public AnnotationChanger(Annotation wrapped, ChangeTable changer) {
        this.wrapped = wrapped;
        this.changer = changer;
        this.properties = new MappedHash();
    }

    /**
     * <code>getWrapped</code> returns the <code>Annotation</code>
     * being remapped.
     *
     * @return an <code>Annotation</code>.
     */
    public Annotation getWrapped() {
        return wrapped;
    }

    /**
     * <code>getMapper</code> returns the <code>ValueChanger</code> being
     * used to remap the <code>Annotation</code>.
     *
     * @return a <code>ValueChanger</code>.
     */
    public ChangeTable getChanger() {
        return changer;
    }

    /**
     * <code>getProperties</code> returns the mapped contents of the
     * underlying <code>Annotation</code> as a <code>Map</code>.
     *
     * @return a <code>Map</code>.
     */
    public Map getProperties() {
        return properties;
    }

    /**
     * <code>propertiesAllocated</code> is a convenience method to see
     * if we have allocated the properties <code>Map</code>.
     *
     * @return a <code>boolean</code> true if the properties have been
     * allocated, false otherwise.
     */
    public boolean propertiesAllocated() {
        return true;
    }

    private class MappedHash extends AbstractMap {
        public int size() {
            return wrapped.asMap().size();
        }

        public Set entrySet() {
            return new WrappedSet(wrapped.asMap().entrySet());
        }
    }

    private class WrappedSet extends AbstractSet {
        private Set entrySet;

        public WrappedSet(Set entrySet) {
            this.entrySet = entrySet;
        }

        public int size() {
            return entrySet.size();
        }

        public Iterator iterator() {
            return new Iterator() {
                    Iterator i = entrySet.iterator();

                    public boolean hasNext() {
                        return i.hasNext();
                    }

                    public Object next() {
                        final Map.Entry entry = (Map.Entry) i.next();

                        return new Map.Entry() {
                                public Object getKey() {
                                    return entry.getKey();
                                }

                                public Object getValue() {
                                  try {
                                    return changer.change(getKey(), entry.getValue());
                                  } catch (ParserException pe) {
                                    throw new BioError(pe);
                                  }
                                }

                                public Object setValue(Object value) {
                                    return entry.setValue(value);
                                }
                            };
                    }
        
                    public void remove() {
                        i.remove();
                    }
                };
        }
    }
}
