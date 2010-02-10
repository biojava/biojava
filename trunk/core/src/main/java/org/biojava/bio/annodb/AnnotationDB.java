package org.biojava.bio.annodb;

import java.util.Iterator;

import org.biojava.bio.AnnotationType;

/**
 * <p>A database of Annotation instances.</p>
 *
 * <p>It is often a pain to provide explicit API for a particular file format,
 * but it is still necisary to present it as some collection of structured
 * objects. Annotation, together with AnnotationType are capable of representing
 * structured data and the tag-value parser API is uniquely suited to creating
 * these from structured text files. AnnotationDB is provided as a way to wrap
 * up a whole collection of Annotation instances so that they can be queried and
 * handled as a unit.</p>
 *
 * @author Matthew Pocock
 * @since 1.3
 */
public interface AnnotationDB {
  /**
   * An AnnotationDB that is always empty.
   */
  public static AnnotationDB EMPTY = new EmptyAnnotationDB();
  
  /**
   * <p>The name of this AnnotationDB.</p>
   *
   * @return the name of this AnnotationDB
   */
  public String getName();
  
  /**
   * <p>
   * Get an AnnotationType that accepts all Annotation instances in this DB.
   * </p>
   *
   * <p>
   * The schema should accept all Annotations in the DB. However, it may hit
   * other Annotations. So, AnnotationType.ALL is always a valid schema.
   * Obviously, the more retrictive it is, the more usefull it becomes for
   * introspection.
   * </p>
   *
   * @return  the schema AnnotationType
   */
  public AnnotationType getSchema();

  /**
   * Loop over each Annotation in this db.
   *
   * @return an Iterator over each item in the DB
   */
  public Iterator iterator();
  
  /**
   * The number of Annotation instances in the DB.
   *
   * @return the size of this DB
   */
  public int size();
  
  /**
   * Find all Annotation instances in this DB that are of a particular type.
   *
   * @param at  the AnnotationType to match
   * @return an AnnotationDB with all matching Annotation instances
   */
  public AnnotationDB filter(AnnotationType at);
  
  /**
   * Find all Annotation instances in this DB and any Annotations that are child
   * properties of these that match an AnnotationType.
   *
   * @param at  the AnnotationType to search with
   * @return an AnnotationDB with all matching Annotation instances,
   *         irregardless of how deep in the hieracy they are
   */
  public AnnotationDB search(AnnotationType at);
}

