package org.biojava.bio.annodb;

import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

import org.biojava.bio.Annotation;
import org.biojava.bio.AnnotationTools;
import org.biojava.bio.AnnotationType;

/**
 * <p>A no-frills implementation of AnnotationDB.</p>
 *
 * @author Matthew Pocock
 * @since 1.3
 */
public class SimpleAnnotationDB implements AnnotationDB {
  private final String name;
  private final Set anns;
  private final AnnotationType schema;

  /**
   * Create a no-frills AnnotationDB instancec.
   *
   * @param name    the name of this Annotation DB
   * @param anns    a Set of Annotation instances it contains
   * @param schema  an AnnotationType schema that applies to them
   */
  public SimpleAnnotationDB(String name, Set anns, AnnotationType schema) {
    this.name = name;
    this.anns = anns;
    this.schema = schema;
  }
  
  public String getName() {
    return this.name;
  }
  
  public int size() {
    return anns.size();
  }
  
  public AnnotationType getSchema() {
    return schema;
  }

  public Iterator iterator() {
    return anns.iterator();
  }
  
  public AnnotationDB filter(AnnotationType at) {
    Set hits = new HashSet();
    
    AnnotationType intersection = AnnotationTools.intersection(schema, at);
    if(intersection != AnnotationType.NONE) {
      for(Iterator i = anns.iterator(); i.hasNext(); ) {
        Annotation ann = (Annotation) i.next();
        if(at.instanceOf(ann)) {
          hits.add(ann);
        }
      }
    }
    
    if(hits.isEmpty()) {
      return AnnotationDB.EMPTY;
    } else {
      return new SimpleAnnotationDB("", hits, intersection);
    }
  }
  
  public AnnotationDB search(AnnotationType at) {
    Set hits = new HashSet();
    
    for(Iterator i = anns.iterator(); i.hasNext(); ) {
      Annotation ann = (Annotation) i.next();
      hits.addAll(AnnotationTools.searchAnnotation(ann, at));
    }
    
    if(hits.isEmpty()) {
      return AnnotationDB.EMPTY;
    } else {
      return new SimpleAnnotationDB("", hits, at);
    }
  }
}

