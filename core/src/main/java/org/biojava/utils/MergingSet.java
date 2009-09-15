package org.biojava.utils;

import java.util.AbstractSet;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

/**
 *
 *
 * @author Matthew Pocock
 */
public class MergingSet
        extends AbstractSet
{
  private final Set sets;
  private Set modifiable;

  public static MergingSet merge(Set first, Set seccond)
  {
    MergingSet ms = new MergingSet();
    ms.add(first);
    ms.add(seccond);
    return ms;
  }

  public static MergingSet modifiableMerge()
  {
    MergingSet ms = new MergingSet();
    ms.modifiable = new HashSet();
    ms.addSet(ms.modifiable);
    return ms;
  }

  public static MergingSet modifiableMerge(Set original)
  {
    MergingSet ms = new MergingSet();
    ms.modifiable = new HashSet();
    ms.addSet(ms.modifiable);
    ms.addSet(original);
    return ms;
  }

  public MergingSet() {
    this.sets = new SmallSet();
  }

  public MergingSet(Set sets) {
    this.sets = new SmallSet(sets);
  }

  public void addSet(Set set) {
    sets.add(set);
  }

  public boolean removeSet(Set set) {
    if(set == modifiable) {
      throw new IllegalArgumentException(
              "Can't remove the set that contains modifications");
    }
    return sets.remove(set);
  }

  public Set getModifiable()
  {
    return modifiable;
  }

  public int size() {
    int size = 0;

    for(Iterator i = sets.iterator(); i.hasNext(); ) {
      Set s = (Set) i.next();
      size += s.size();
    }

    return size;
  }

  public boolean contains(Object o) {
    for (Iterator i = sets.iterator(); i.hasNext();) {
      Set s = (Set) i.next();
      if(s.contains(o)) {
        return true;
      }
    }

    return false;
  }

  public Iterator iterator() {
    return new MergingIterator(sets.iterator());
  }

  public boolean add(Object o)
  {
    if(modifiable == null) {
      throw new UnsupportedOperationException();
    }

    return modifiable.add(o);
  }

  public boolean remove(Object o)
  {
    if(modifiable == null) {
      throw new UnsupportedOperationException();
    }

    if(this.contains(o) && !modifiable.contains(o)) {
      throw new IllegalArgumentException(
              "Can't remove items not added to this merged view");
    }

    return modifiable.remove(o);
  }
}
