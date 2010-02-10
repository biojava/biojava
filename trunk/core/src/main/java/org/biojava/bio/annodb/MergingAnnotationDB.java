package org.biojava.bio.annodb;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;

import org.biojava.bio.AnnotationTools;
import org.biojava.bio.AnnotationType;

/**
 * <p>An AnnotationDB that provides a merged view of a list of underlying DBs.</p>
 *
 * @author Matthew Pocock
 * @since 1.3
 */
public class MergingAnnotationDB implements AnnotationDB {
  private final String name;
  private final List merged;

  /**
   * Create a new MergingAnnotationDB with a name and no DBs to merge.
   *
   * @param name  the name of this DB
   */
  public MergingAnnotationDB(String name) {
    this.name = name;
    this.merged = new ArrayList();
  }

  /**
   * Create a new MergingAnnotationDB with a name and a list of DBs to merge.
   *
   * @param name    the name of this DB
   * @param merged  a list of DBs to merge
   */
  public MergingAnnotationDB(String name, List merged) {
    this.name = name;
    this.merged = new ArrayList(merged);
  }

  /**
   * Add a DB to be merged in this view.
   *
   * @param toAdd the AnnotationDB to add
   */
  public void addAnnotationDB(AnnotationDB toAdd) {
    if(!merged.contains(toAdd)) {
      merged.add(toAdd);
    }
  }

  /**
   * Remove a DB from this view.
   *
   * @param toRemove  the AnnotationDB to remove
   */
  public void removeAnnotationDB(AnnotationDB toRemove) {
    merged.remove(toRemove);
  }

  /**
   * Return a list of merged DBs. This can be modified independantly of this DB.
   *
   * @return a List of merged DBs
   */
  public List getMerged() {
    return new ArrayList(merged);
  }
  
  public String getName() {
    return name;
  }
  
  public AnnotationType getSchema() {
    AnnotationType schema = AnnotationType.NONE;
    
    for(Iterator i = merged.iterator(); i.hasNext(); ) {
      AnnotationDB db = (AnnotationDB) i.next();
      schema = AnnotationTools.union(schema, db.getSchema());
    }
    
    return schema;
  }

  public Iterator iterator() {
    return new Iterator() {
      Iterator ii;
      Iterator ci;
      Object item;
      
      {
        ii = merged.iterator();
       EVERYTHING:
        while(item == null) {
          if(ii.hasNext()) {
            AnnotationDB adb = (AnnotationDB) ii.next();
            ci = adb.iterator();
            if(ci.hasNext()) {
              item = ci.next();
            }
          } else {
            break EVERYTHING;
          }
        }
      }
      
      public boolean hasNext() {
        return item != null;
      }
      
      public Object next() {
        Object it = item;
        item = _next();
        return it;
      }
      
      private Object _next() {
        while(!ci.hasNext()) {
          if(ii.hasNext()) {
            AnnotationDB adb = (AnnotationDB) ii.next();
            ci = adb.iterator();
          } else {
            return null;
          }
        }
        
        return ci.next();
      }
      
      public void remove() {
        throw new NoSuchElementException();
      }
    };
  }
  
  public int size() {
    int size = 0;
    
    for(Iterator dbi = merged.iterator(); dbi.hasNext(); ) {
      size += ((AnnotationDB) dbi.next()).size();
    }
    
    return size;
  }
  
  public AnnotationDB filter(AnnotationType at) {
    List anns = new ArrayList();
    
    for(Iterator i = merged.iterator(); i.hasNext(); ) {
      AnnotationDB adb = (AnnotationDB) i.next();
      AnnotationDB res = adb.filter(at);
      if(res.size() > 0) {
        anns.add(res);
      }
    }
    
    if(anns.isEmpty()) {
      return AnnotationDB.EMPTY;
    } else if(anns.size() == 1) {
      return (AnnotationDB) anns.get(0);
    } else {
      return new MergingAnnotationDB("", anns);
    }
  }
  
  public AnnotationDB search(AnnotationType at) {
    List anns = new ArrayList();
    
    for(Iterator i = merged.iterator(); i.hasNext(); ) {
      AnnotationDB adb = (AnnotationDB) i.next();
      AnnotationDB res = adb.search(at);
      if(res.size() > 0) {
        anns.add(res);
      }
    }
    
    if(anns.isEmpty()) {
      return AnnotationDB.EMPTY;
    } else if(anns.size() == 1) {
      return (AnnotationDB) anns.get(0);
    } else {
      return new MergingAnnotationDB("", anns);
    }
  }
}

