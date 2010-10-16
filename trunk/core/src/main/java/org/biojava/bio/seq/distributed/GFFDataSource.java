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

import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

import org.biojava.bio.Annotation;
import org.biojava.bio.BioError;
import org.biojava.bio.BioException;
import org.biojava.bio.program.gff.GFFEntrySet;
import org.biojava.bio.program.gff.GFFRecord;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.FeatureFilter;
import org.biojava.bio.seq.FeatureHolder;
import org.biojava.bio.seq.MergeFeatureHolder;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.impl.SimpleSequence;
import org.biojava.bio.symbol.DummySymbolList;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.utils.ChangeVetoException;

/**
 * Use a GFFEntrySet as a DataSource for adding annotation to sequences.
 *
 * Instantiate this and add it to an instance of DistributeSequenceDB. All
 * of the GFF features that have sequence fields matching sequence IDs in the
 * db will be merged in.
 * 
 * @author Thomas Down
 * @author Matthew Pocock

 */
public class GFFDataSource implements DistDataSource {
    private GFFEntrySet gffe;
    private Set ids;
    private Map id2seq;
    private MergeFeatureHolder delegateFH;

    public GFFDataSource(GFFEntrySet gffe) {
	this.gffe = gffe;
        this.id2seq = new HashMap();
        delegateFH = new MergeFeatureHolder();
    }

    public boolean hasSequence(String id) throws BioException {
	return false;
    }

    public boolean hasFeatures(String id) throws BioException {
	return ids(false).contains(id);
    }

    public FeatureHolder getFeatures(FeatureFilter ff) throws BioException {
	return getDelegateFH(true).filter(ff);
    }

    public FeatureHolder getFeatures(String id, FeatureFilter ff, boolean recurse) throws BioException {
	if (! hasFeatures(id)) {
	    return FeatureHolder.EMPTY_FEATURE_HOLDER;
	}
	
        Sequence seq = populateDelegateFH(id);
        return seq.filter(ff, recurse);
    }

    private Sequence populateDelegateFH(String id) {
      Sequence seq = (Sequence) id2seq.get(id);

      if(seq == null) {
        SymbolList dummy = new DummySymbolList(DNATools.getDNA(), 1000000000);
        seq = new SimpleSequence(dummy, id, id, Annotation.EMPTY_ANNOTATION);

        try {
	  seq = gffe.getAnnotator().annotate(seq);
          delegateFH.addFeatureHolder(seq);
          id2seq.put(id, seq);
	} catch (ChangeVetoException cve) {
	  throw new BioError(cve);
	} catch (BioException be) {
          throw new BioError(be);
        }
      }

      return seq;
    }

    private FeatureHolder getDelegateFH(boolean populate)
    throws BioException {
      if(populate == true) {
        for(Iterator i = ids(true).iterator(); i.hasNext(); ) {
          populateDelegateFH((String) i.next());
        }
      }

      return delegateFH;
    }

    public Sequence getSequence(String id) throws BioException {
	throw new BioException();
    }

    public Set ids(boolean topLevel) throws BioException {
	if (ids == null) {
	    Set _ids = new HashSet();

	    for (Iterator i = gffe.lineIterator(); i.hasNext(); ) {
		Object o = i.next();
		if (o instanceof GFFRecord) {
		    GFFRecord rec = (GFFRecord) o;
		    _ids.add(rec.getSeqName());
		}
	    }

	    ids = Collections.unmodifiableSet(_ids);
	}

	return ids;
    }
}
