package org.biojava.bio.program.unigene;

import java.util.Map;

import org.biojava.bio.BioException;
import org.biojava.utils.AbstractChangeable;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.JDBCConnectionPool;
import org.biojava.utils.cache.WeakValueHashMap;

/**
 * An implementation of UnigeneDB that fetches data from an SQL database.
 *
 * @author Matthew Pocock
 */
class SQLUnigeneDB
extends AbstractChangeable
implements UnigeneDB {
  private final Map clusterCache;
  
  public SQLUnigeneDB(JDBCConnectionPool connPool) {
    this.clusterCache = new WeakValueHashMap();
  }
  
  public UnigeneCluster getCluster(String clusterID)
  throws BioException {
    UnigeneCluster cluster = (UnigeneCluster) clusterCache.get(clusterID);
    
    if(cluster == null) {
      clusterCache.put(clusterID, cluster = fetchCluster(clusterID));
    }
    
    return cluster;
  }
  
  public UnigeneCluster addCluster(UnigeneCluster cluster)
  throws BioException, ChangeVetoException {

    return fetchCluster(cluster.getID());
  }
  
  public UnigeneCluster fetchCluster(String clusterID) {
    return null;
  }
}
