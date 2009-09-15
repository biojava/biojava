package org.biojava.bio.annodb;

import java.util.Collections;
import java.util.Iterator;

import org.biojava.bio.AnnotationType;

class EmptyAnnotationDB implements AnnotationDB {
  public String getName() { return "EMPTY"; }
  public AnnotationType getSchema() { return AnnotationType.NONE; }
  public Iterator iterator() { return Collections.EMPTY_LIST.iterator(); }
  public int size() { return 0; }
  public AnnotationDB filter(AnnotationType at) { return this; }
  public AnnotationDB search(AnnotationType at) { return this; }
}
