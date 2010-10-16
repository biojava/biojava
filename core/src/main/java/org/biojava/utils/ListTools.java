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

import java.io.Serializable;
import java.util.AbstractList;
import java.util.AbstractMap;
import java.util.AbstractSet;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * @author Matthew Pocock
 * @author Thomas Down
 * @author Mark Schreiber
 */
public class ListTools implements Serializable{
  public static Iterator nonRemoveIterator(Iterator i) {
    final Iterator it = i;
    return new Iterator() {
      public boolean hasNext() {
        return it.hasNext();
      }
      public Object next() {
        return it.next();
      }
      public void remove() {
        throw new UnsupportedOperationException();
      }
    };
  }

  public static List createList(List l) {
    switch (l.size()) {
      case 0:
        return Collections.EMPTY_LIST;
      case 1:
        return Collections.nCopies(1, l.get(0));
      case 2:
        Doublet d = new Doublet(l.get(0), l.get(1));
        return d;
      case 3:
        Triplet t = new Triplet(l.get(0), l.get(1), l.get(2));
        return t;
      default:
        return new ArrayList(l);
    }
  }

  public static List createList(Object[] a) {
    switch (a.length) {
      case 0:
        return Collections.EMPTY_LIST;
      case 1:
        return Collections.nCopies(1, a[0]);
      case 2:
        Doublet d = new Doublet(a[0],a[1]);
        return d;
      case 3:
        Triplet t = new Triplet(a[0],a[1],a[2]);
        return t;
      default:
        return Arrays.asList(a);
    }
  }

  /**
   * Create a new SeriesList with the given leader, trailer and size.
   *
   * @param leader  the String that will prefix the index
   * @param trailer the String that will suffix the index
   * @param size  the length of the list
   * @throws NullPointerException if leader or trailer are null (use the empty
   *   string instead)
   * @throws IllegalArgumentException if the size is negative
   */
  public static SeriesList createSeriesList(
    String leader,
    String trailer,
    int size
  ) {
    return new SeriesList(leader, trailer, size);
  }

  public static class Doublet extends AbstractList implements Serializable {
    private Object a;
    private Object b;

    public Doublet() {}
    public Doublet(Object a, Object b) {
      this();
      set(a, b);
    }

    public void set(Object a, Object b) {
      this.a = a;
      this.b = b;
    }

    public void setA(Object a) {
      this.a = a;
    }

    public void setB(Object b) {
      this.b = b;
    }

    public Object getA() {
      return a;
    }

    public Object getB() {
      return b;
    }

    public int size() {
      return 2;
    }

    public Object get(int indx) {
      switch (indx) {
        case 0:
          return a;
        case 1:
          return b;
        default:
          throw new IndexOutOfBoundsException("indx must be 0 or 1");
      }
    }

    public Iterator getIterator() {
      return new Iterator() {
        private int indx = 0;

        public boolean hasNext() {
          return indx < 2;
        }

        public Object next() {
          return get(indx++);
        }

        public void remove()
        throws UnsupportedOperationException {
          throw new UnsupportedOperationException();
        }
      };
    }

    public int hashCode() {
      int hashcode = 1;
      hashcode = 31*hashcode + a.hashCode();
      hashcode = 31*hashcode + b.hashCode();
      return hashcode;
    }

    public boolean equals(Object o) {
      if(! (o instanceof List) ) {
        return false;
      }

      List other = (List) o;
      if(other.size() != 2) {
        return false;
      }

      return other.get(0).equals(a) && other.get(1).equals(b);
    }
  }

  public static class Triplet extends AbstractList implements Serializable {
    private Object a;
    private Object b;
    private Object c;

    public Triplet() {}
    public Triplet(Object a, Object b, Object c) {
      this();
      set(a, b, c);
    }

    public void set(Object a, Object b, Object c) {
      this.a = a;
      this.b = b;
      this.c = c;
    }

    public void setA(Object a) {
      this.a = a;
    }

    public void setB(Object b) {
      this.b = b;
    }

    public void setC(Object c) {
      this.c = c;
    }

    public Object getA() {
      return a;
    }

    public Object getB() {
      return b;
    }

    public Object getC() {
      return c;
    }

    public int size() {
      return 3;
    }

    public Object get(int indx) {
      switch (indx) {
        case 0:
          return a;
        case 1:
          return b;
        case 2:
          return c;
        default:
          throw new IndexOutOfBoundsException("indx must be 0 or 1");
      }
    }

    public Iterator getIterator() {
      return new Iterator() {
        private int indx = 0;

        public boolean hasNext() {
          return indx < 3;
        }

        public Object next() {
          return get(indx++);
        }

        public void remove()
        throws UnsupportedOperationException {
          throw new UnsupportedOperationException();
        }
      };
    }

    public int hashCode() {
      int hashcode = 1;
      hashcode = 31*hashcode + a.hashCode();
      hashcode = 31*hashcode + b.hashCode();
      hashcode = 31*hashcode + c.hashCode();
      return hashcode;
    }

    public boolean equals(Object o) {
      if(! (o instanceof List) ) {
        return false;
      }

      List other = (List) o;
      if(other.size() != 3) {
        return false;
      }

      return other.get(0).equals(a) && other.get(1).equals(b) && other.get(2).equals(c);
    }
  }

  /**
   * A list that represents a series of values.
   *
   * <p>This provides a simple list implementation that synthesises elements from
   * a leading and trailing string and the index into the list.</p>
   *
   * <p>For example, a SeriesList with leader "" and trailer ":" will contain
   * values like "0:", "1:", "2:" and so on. A SeriesList with leader "Chapter "
   * and trailer "" will have values like "Chapter 5".</p>
   *
   * @author Matthew Pocock
   * @since 1.4
   */
  public static class SeriesList
  extends AbstractList {
    private final String leader;
    private final String trailer;
    private final int size;

    private SeriesList(String leader, String trailer, int size) {
      if(leader == null) {
        throw new NullPointerException(
        "Leader was null. Use the empty string instead");
      }

      if(trailer == null) {
        throw new NullPointerException(
        "Trailer was null. Use the empty string instead");
      }

      if(size < 0) {
        throw new IllegalArgumentException(
          "Size must be zero or positive: " + size );
      }

      this.leader = leader;
      this.trailer = trailer;
      this.size = size;
    }

    public String getLeader() {
      return leader;
    }

    public String getTrailer() {
      return trailer;
    }

    public int size() {
      return size;
    }

    public Object get(int indx) {
      return leader + indx + trailer;
    }
  }

  public static List mapList(final List list,
                             final Mapper mapper)
  {
    return new AbstractList() {
      public Object get(int index)
      {
        return mapper.map(list.get(index));
      }

      public int size()
      {
        return list.size();
      }
    };
  }

  public static Set mapSet(final Set set,
                           final Mapper mapper)
  {
    return new AbstractSet() {
      public Iterator iterator()
      {
        return new Iterator() {
          Iterator i = set.iterator();
          public boolean hasNext()
          {
            return i.hasNext();
          }

          public Object next()
          {
            return mapper.map(i.next());
          }

          public void remove()
          {
            i.remove();
          }
        };
      }

      public int size()
      {
        return set.size();
      }
    };
  }

  public static Map mapMap(final Map map,
                           final Mapper keyMapper,
                           final Mapper valMapper)
  {
    return new AbstractMap() {
      public Set entrySet()
      {
        return mapSet(map.entrySet(), new Mapper() {
          public Object map(Object val)
          {
            final Map.Entry ent = (Map.Entry) val;
            return new Map.Entry() {
              public Object getKey()
              {
                return keyMapper.map(ent.getKey());
              }

              public Object getValue()
              {
                return valMapper.map(ent.getValue());
              }

              public Object setValue(Object value)
              {
                throw new UnsupportedOperationException();
              }
            };
          }
        });
      }
    };
  }

  /**
   * Maps one object to another.
   *
   * @author Matthew Pocock
   * @since 1.4
   */
  public interface Mapper {
    /**
     * Map the object.
     *
     * @param val   the object to map
     * @return      the new value
     */
    public Object map(Object val);
  }

  public static final Mapper NULL_MAPPER = new Mapper() {
    public Object map(Object val)
    {
      return val;
    }
  };
}
