package org.biojava.bio.program.unigene;

import org.biojava.bio.Annotatable;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.db.SequenceDB;

/**
 * <p>
 * A single unigene cluster.
 * </p>
 *
 * <p>
 * This represents all of the information available about a single unigene
 * cluster. This includes the sequences that are in it, the unique sequence that
 * is its representative, the cluster id and all the annotation available via
 * the data file. Much of the annotation may be accessible via the
 * getAnnotation() method.
 * </p>
 *
 * <p>
 * There is much more information in the Unigene clusters than just the id,
 * title and sequences. This is all stoored in the annotation associated with
 * the cluster. The annotation bundle conforms to the schema in 
 * <code>UnigeneTools.UNIGENE_ANNOTATION</code>.
 * </p>
 *
 * @author Matthew Pocock
 */
public interface UnigeneCluster
extends Annotatable {
  /**
   * The public unigene ID.
   *
   * @return get the cluster ID as a String
   */
  public String getID();
  
  /**
   * The cluster title.
   *
   * @return the cluster title as a String
   */
  public String getTitle();
  
  /**
   * All sequences that map to this cluster.
   *
   * @return a SequenceDB of all sequences mapping to this cluster.
   */
  public SequenceDB getAll();
  
  /**
   * The unique sequence that is used as a representative for this cluster.
   *
   * @return the cluster's unique Sequence
   */
  public Sequence getUnique();
}
