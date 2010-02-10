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

package org.biojava.bio.seq.distributed;

import java.util.Iterator;
import java.util.List;
import java.util.Set;

import org.biojava.bio.Annotation;
import org.biojava.bio.BioException;
import org.biojava.bio.BioRuntimeException;
import org.biojava.bio.SmallAnnotation;
import org.biojava.bio.seq.ComponentFeature;
import org.biojava.bio.seq.Feature;
import org.biojava.bio.seq.FeatureFilter;
import org.biojava.bio.seq.FeatureHolder;
import org.biojava.bio.seq.MergeFeatureHolder;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.projection.ProjectedFeatureHolder;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.Edit;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.utils.AbstractChangeable;
import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeListener;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;

/**
 * Sequence from the meta-DAS system.
 *
 * @author Thomas Down
 * @author Matthew Pocock
 * @since 1.2
 */

class DistributedSequence
        extends
        AbstractChangeable
        implements
        Sequence {
  private DistributedSequenceDB db;
  private DistDataSource seqSource;
  private Set featureSources;
  private String id;

  private SymbolList symbols;
  private MergeFeatureHolder mfh;

  private transient ChangeListener dsListener;

  DistributedSequence(String id,
                      DistributedSequenceDB db,
                      DistDataSource seqSource,
                      Set featureSources) {
//System.err.println("*** Constructing DistributedSequence: " + id);

    this.id = id;
    this.seqSource = seqSource;
    this.featureSources = featureSources;
    this.db = db;
    installListener();
  }

  private void installListener() {
    dsListener = new ChangeListener() {
      public void preChange(ChangeEvent cev)
              throws ChangeVetoException {
        if (cev.getPrevious() == seqSource) {
          throw new ChangeVetoException(cev, "Can't remove this datasource, since it is providing sequence data");
        }

        if (hasListeners()) {
          getChangeSupport(ChangeType.UNKNOWN).firePreChangeEvent(makeChainedEvent(cev));
        }
      }

      public void postChange(ChangeEvent cev) {
        if (cev.getType() == DistributedSequenceDB.DATASOURCE_SELECTION) {
          if (cev.getChange() != null && cev.getPrevious() == null) {
            DistDataSource added = (DistDataSource) cev.getChange();
            featureSources.add(added);
          } else if (cev.getChange() == null && cev.getPrevious() != null) {
            DistDataSource removed = (DistDataSource) cev.getChange();
            featureSources.remove(removed);
          }
        }

        mfh = null;  // C'mon, we can do better than that...

        if (hasListeners()) {
          getChangeSupport(ChangeType.UNKNOWN).firePostChangeEvent(makeChainedEvent(cev));
        }
      }

      private ChangeEvent makeChainedEvent(ChangeEvent cev) {
        return new ChangeEvent(DistributedSequence.this,
                               FeatureHolder.FEATURES,
                               null, null,
                               cev);
      }
    };
    db.addChangeListener(dsListener, DistributedSequenceDB.DATASOURCE);
  }

  public DistributedSequenceDB getSequenceDB() {
    return db;
  }

  public String getName() {
    return id;
  }

  public String getURN() {
    return id;
  }

  public Annotation getAnnotation() {
    return Annotation.EMPTY_ANNOTATION;
  }


  public Alphabet getAlphabet() {
    return getSymbols().getAlphabet();
  }

  public int length() {
    return getSymbols().length();
  }

  public Symbol symbolAt(int i) {
    return getSymbols().symbolAt(i);
  }

  public SymbolList subList(int start, int end) {
    return getSymbols().subList(start, end);
  }

  public List toList() {
    return getSymbols().toList();
  }

  public Iterator iterator() {
    return getSymbols().iterator();
  }

  public String seqString() {
    return getSymbols().seqString();
  }

  public String subStr(int start, int end) {
    return getSymbols().subStr(start, end);
  }

  public void edit(Edit e)
          throws ChangeVetoException {
    throw new ChangeVetoException("Can't edit sequence in EGADS -- or at least not yet...");
  }

  public Iterator features() {
    return getFeatures().features();
  }
  
  
    private FeatureHolder components = null;

    public FeatureHolder filter(FeatureFilter ff, boolean recurse) {
        if (recurse) {
            MergeFeatureHolder results = new MergeFeatureHolder();
            try {
                results.addFeatureHolder(getFeatures().filter(new FeatureFilter.And(new FeatureFilter.Not(new FeatureFilter.ByAncestor(new FeatureFilter.ByClass(ComponentFeature.class))), ff), recurse));
                if (components == null) {
                    components = getFeatures().filter(new FeatureFilter.And(FeatureFilter.top_level, new FeatureFilter.ByClass(ComponentFeature.class)));
                }
                for (Iterator i = components.features(); i.hasNext(); ) {
                    FeatureHolder fh = ((Feature) i.next()).filter(ff);
                    if (fh.countFeatures() > 0) {
                        results.addFeatureHolder(fh);
                    }
                }
            } catch (Exception ex) {
                throw new BioRuntimeException(ex);
            }
            return results;
        } else {
            return getFeatures().filter(ff, recurse);
        }
    }

    public FeatureHolder filter(FeatureFilter ff) {
        return filter(ff, true);
    }

  public void removeFeature(Feature f)
          throws ChangeVetoException {
    throw new ChangeVetoException("Can't edit sequence in EGADS -- or at least not yet...");
  }

  public Feature createFeature(Feature.Template f)
          throws ChangeVetoException {
    throw new ChangeVetoException("Can't edit sequence in EGADS -- or at least not yet...");
  }

  public boolean containsFeature(Feature f) {
    return getFeatures().containsFeature(f);
  }

  public int countFeatures() {
    return getFeatures().countFeatures();
  }

  protected SymbolList getSymbols() {
    if (symbols == null) {
      try {
        symbols = seqSource.getSequence(id);
      } catch (BioException ex) {
        throw new BioRuntimeException(ex);
      }
    }

    return symbols;
  }

  protected FeatureHolder getFeatures() {
    if (mfh == null) {
      mfh = new MergeFeatureHolder();
      for (Iterator i = featureSources.iterator(); i.hasNext();) {
        DistDataSource dds = (DistDataSource) i.next();
        try {
          Annotation ann = new SmallAnnotation();
          ann.setProperty("source", dds);
          mfh.addFeatureHolder(
              new ProjectedFeatureHolder(
                  new DistProjectionContext(
                      dds.getFeatures(id, FeatureFilter.all, false),
                      this,
                      ann
                  )
              )
          );
        } catch (Exception ex) {
          ex.printStackTrace();
        }
      }
    }
    return mfh;
  }

  public FeatureFilter getSchema() {
    return FeatureFilter.all; // FIXME!
  }
}
