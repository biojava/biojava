package org.biojava.bio.program.unigene;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.AbstractSet;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

import org.biojava.bio.Annotation;
import org.biojava.bio.BioError;
import org.biojava.bio.BioException;
import org.biojava.bio.program.indexdb.BioStore;
import org.biojava.bio.program.indexdb.Record;
import org.biojava.bio.program.tagvalue.AnnotationBuilder;
import org.biojava.bio.program.tagvalue.Parser;
import org.biojava.bio.program.tagvalue.ParserListener;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.db.AbstractSequenceDB;
import org.biojava.bio.seq.db.CachingSequenceDB;
import org.biojava.bio.seq.db.IllegalIDException;
import org.biojava.bio.seq.db.SequenceDB;
import org.biojava.bio.seq.io.FastaFormat;
import org.biojava.bio.seq.io.SeqIOTools;
import org.biojava.bio.seq.io.SequenceFormat;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.ParserException;
import org.biojava.utils.Unchangeable;
import org.biojava.utils.cache.WeakValueHashMap;
import org.biojava.utils.io.RandomAccessReader;

/**
 * A UnigeneDB that uses flat-file indexing.
 *
 * @author Matthew Pocock
 */
class FlatFileUnigeneDB
extends Unchangeable
implements UnigeneDB {
  private final BioStore dataStore;
  private final BioStore allStore;

  private final Map clusterCache;
  private final Map allCache;
  private final SequenceDB uniqueDB;
  private final ParserListener dataPL;
  private final Parser dataParser;
  private final AnnotationBuilder dataBuilder;

  public FlatFileUnigeneDB(
    BioStore dataStore,
    BioStore liStore,
    BioStore uniqueStore,
    BioStore allStore
  ) throws BioException {
    this.dataStore = dataStore;
    this.allStore = allStore;

    try {
      clusterCache = new WeakValueHashMap();
      allCache = new WeakValueHashMap();

      FastaFormat fasta = new FastaFormat();
      uniqueDB = new CachingSequenceDB(
        new BioIndexSequenceDB(
          uniqueStore,
          fasta
        )
      );

      dataBuilder = new AnnotationBuilder(
        UnigeneTools.UNIGENE_ANNOTATION
      );
      dataPL = UnigeneTools.buildDataParser(dataBuilder);
      dataParser = new Parser();
    } catch (ParserException pe) {
      throw new BioException("Could not initialize unigene DB", pe);
    }
  }

  public UnigeneCluster getCluster(String clusterID)
  throws BioException {
    UnigeneCluster cluster = (UnigeneCluster) clusterCache.get(clusterID);
    if(cluster == null) {
      synchronized(dataParser) {
        cluster = (UnigeneCluster) clusterCache.get(clusterID);
        if(cluster == null) { // break race condition
          try {
            Record rec = dataStore.get(clusterID);
            RandomAccessReader rar = new RandomAccessReader(rec.getFile());
            rar.seek(rec.getOffset());
            BufferedReader reader = new BufferedReader(rar);
            dataParser.read(reader, dataPL.getParser(), dataPL.getListener());
          } catch (IOException ioe) {
            throw new BioException("Failed to load cluster: " + clusterID, ioe);
          } catch (ParserException pe) {
            throw new BioException("Failed to parse cluster: " + clusterID, pe);
          }
        }
        cluster = new AnnotationCluster(dataBuilder.getLast());
        clusterCache.put(clusterID, cluster);
      }
    }
    return cluster;
  }

  public SequenceDB getAll(String clusterID)
  throws BioException {
    SequenceDB db = (SequenceDB) allCache.get(clusterID);

    if(db == null) {
      synchronized(db) {
        db = (SequenceDB) allCache.get(clusterID);
        if(db == null) {
          allCache.put(clusterID, db = new AllDB(getCluster(clusterID), allStore));
        }
      }
    }

    return db;
  }

  public UnigeneCluster addCluster(UnigeneCluster cluster)
  throws BioException, ChangeVetoException {
    throw new ChangeVetoException("Can't alter a file-based unigene installation");
  }

  public Sequence getUnique(String clusterID)
  throws IllegalIDException, BioException {
    return uniqueDB.getSequence(clusterID);
  }

  private class AnnotationCluster
  extends Unchangeable
  implements UnigeneCluster {
    private Annotation ann;

    public AnnotationCluster(Annotation ann) {
      this.ann = ann;
    }

    public String getID() {
      return (String) ann.getProperty("ID");
    }

    public String getTitle() {
      return (String) ann.getProperty("TITLE");
    }

    public SequenceDB getAll() {
      try {
        return FlatFileUnigeneDB.this.getAll(getID());
      } catch (BioException be) {
        throw new BioError(be);
      }
    }

    public Sequence getUnique() {
      try {
        return FlatFileUnigeneDB.this.getUnique(getID());
      } catch (BioException be) {
        throw new BioError(be);
      }
    }

    public Annotation getAnnotation() {
      return ann;
    }
  }

  private static class BioIndexSequenceDB
  extends AbstractSequenceDB {
    private final BioStore store;
    private Set ids = null;

    public BioIndexSequenceDB(BioStore store, SequenceFormat format) {
      this.store = store;
    }

    public Set ids() {
      if(ids == null) {
        ids = new AbstractSet() {
          public int size() {
            return store.getRecordList().size();
          }

          public boolean contains(Object o) {
            return store.get((String) o) != null;
          }

          public Iterator iterator() {
            return store.getRecordList().iterator();
          }
        };
      }

      return ids;
    }

    public String getName() {
      return "UniqueStore";
    }

    public Sequence getSequence(String id)
    throws BioException {
      try {
        Record rec = store.get(id);
        RandomAccessReader rar = new RandomAccessReader(rec.getFile());
        rar.seek(rec.getOffset());
        BufferedReader reader = new BufferedReader(rar);
        return SeqIOTools.readFastaDNA(reader).nextSequence();
      } catch (IOException ioe) {
        throw new BioException(ioe);
      }
    }
  }

  private static class AllDB
  extends AbstractSequenceDB {
    private final Set ids;
    private final BioStore store;
    private final String name;

    public AllDB(UnigeneCluster cluster, BioStore store) {
      this.name = "All:" + cluster.getID();
      ids = new HashSet();
      this.store = store;

      Annotation ann = cluster.getAnnotation();
      Set seqs = (Set) ann.getProperty("SEQUENCES");
      for(Iterator i = seqs.iterator(); i.hasNext(); ) {
        Annotation sa = (Annotation) i.next();
        ids.add(sa.getProperty("ACC"));
      }
    }

    public Set ids() {
      return ids;
    }

    public String getName() {
      return name;
    }

    public Sequence getSequence(String id)
    throws BioException {
      try {
        Record rec = store.get(id);
        RandomAccessReader rar = new RandomAccessReader(rec.getFile());
        rar.seek(rec.getOffset());
        BufferedReader reader = new BufferedReader(rar);
        return SeqIOTools.readFastaDNA(reader).nextSequence();
      } catch (IOException ioe) {
        throw new BioException(ioe);
      }
    }
  }
}
