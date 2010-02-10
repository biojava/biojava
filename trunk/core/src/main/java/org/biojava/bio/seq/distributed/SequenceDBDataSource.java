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

import java.util.Set;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.FeatureFilter;
import org.biojava.bio.seq.FeatureHolder;
import org.biojava.bio.seq.FilterUtils;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.db.IllegalIDException;
import org.biojava.bio.seq.db.SequenceDB;

/**
 * Turn an entire SequenceDB instance into a DistDataSource.
 * This is very usefull for providing the 'reference' sequence and feature set
 * upon which you can layer any other features you have.
 * 
 * @author Thomas Down
 * 
 */
public class SequenceDBDataSource implements DistDataSource {
    private SequenceDB db;

    public SequenceDBDataSource(SequenceDB seqDB) {
	this.db = seqDB;
    }

    public boolean hasSequence(String id) throws BioException {
	return db.ids().contains(id);
    }

    public boolean hasFeatures(String id) throws BioException {
	return hasSequence(id);
    }

    public FeatureHolder getFeatures(FeatureFilter ff) throws BioException {
	throw new BioException();
    }

    public FeatureHolder getFeatures(String id, FeatureFilter ff, boolean recurse) throws BioException {
	FeatureHolder fh;
	try {
	    fh = db.getSequence(id);
	} catch (IllegalIDException ex) {
	    return FeatureHolder.EMPTY_FEATURE_HOLDER;
	}
	
	if (recurse == false && FilterUtils.areProperSubset(FeatureFilter.all, ff)) {
	    return fh;
	} else {
	    return fh.filter(ff, recurse);
	}
    }

    public Sequence getSequence(String id) throws BioException {
	return db.getSequence(id);
    }

    public Set ids(boolean topLevel) throws BioException {
	return db.ids();
    }
}
