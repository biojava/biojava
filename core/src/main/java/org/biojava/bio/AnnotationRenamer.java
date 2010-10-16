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

import org.biojava.bio.program.tagvalue.PropertyChanger;

/**
 * <code>AnnotationRenamer</code> remaps the keys of an
 * <code>Annotation</code> to new keys specified by a
 * <code>TagMapper</code>. This will rename properties, but not alter their
 * values.
 * For writing light-weigth adaptors to project one type of
 * Annotation to another using a TagMapper.
 * @since 1.3
 * @author Matthew Pocock
 * @author <a href="mailto:kdj@sanger.ac.uk">Keith James</a> (docs)
 *
 */
public class AnnotationRenamer extends AbstractAnnotation
{
    private final Annotation wrapped;
    private final PropertyChanger mapper;
    private final Map properties;
  
    /**
     * Creates a new <code>AnnotationRenamer</code> using the
     * specified <code>TagMapper</code> to remap its keys.
     *
     * @param wrapped an <code>Annotation</code>.
     * @param mapper a <code>TagMapper</code>.
     */
    public AnnotationRenamer(Annotation wrapped, PropertyChanger mapper) {
        this.wrapped = wrapped;
        this.mapper = mapper;
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
     * <code>getMapper</code> returns the <code>TagMapper</code> being
     * used to remap the <code>Annotation</code>.
     *
     * @return a <code>TagMapper</code>.
     */
    public PropertyChanger getMapper() {
        return mapper;
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
     * <code>propertiesAllocated</code> Javadoc FIXME - this overrides
     * a protected method and I'm not sure why (KJ).
     *
     * @return a <code>boolean</code>.
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
                                    return mapper.getNewTag(entry.getKey());
                                }

                                public Object getValue() {
                                    return entry.getValue();
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
