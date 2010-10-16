package org.biojava.bio.program.unigene;

import org.biojava.bio.BioException;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.Changeable;

/**
 * A database of Unigene clusters.
 *
 * @author Matthew Pocock
 */
public interface UnigeneDB
extends Changeable {
  /**
   * Fetch a cluster by its cluster id.
   *
   * @param clusterID  the cluster ID as a String
   * @return the UnigeneCluster for that ID
   * @throws BioException if there is no known cluster by that ID or if there
   *         was an error fetching it
   */
  public UnigeneCluster getCluster(String clusterID)
  throws BioException;

  /**
   * Add a cluster to a database.
   *
   * @param cluster  the UnigeneCluster to add
   * @return a (possibly new) UnigeneCluster that is equivalent to
   *         <code>cluster</code> but is served from this <code>UnigeneDB</code>
   *         instance
   */
  public UnigeneCluster addCluster(UnigeneCluster cluster)
  throws BioException, ChangeVetoException;
}
