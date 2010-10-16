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

package org.biojava.bio.seq.db.biosql;

import java.util.ArrayList;
import java.util.List;

import org.biojava.bio.BioError;
import org.biojava.bio.BioException;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.Feature;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.seq.io.ParseException;
import org.biojava.bio.seq.io.SeqIOAdapter;
import org.biojava.utils.ChangeVetoException;

/**
 * Class for instantiating BioSQL features.
 *
 * @author Thomas Down
 * @deprecated Use hibernate and org.biojavax.bio.db.*
 * @since 1.3
 */

abstract class BioSQLFeatureReceiver extends SeqIOAdapter {
    private List stack = new ArrayList();
    private Sequence seq;
    private BioSQLSequenceDB seqDB;

    BioSQLFeatureReceiver(Sequence seq) {
	this.seq = seq;
	this.seqDB = ((BioSQLSequenceI) seq).getSequenceDB();
    }

    BioSQLFeatureReceiver(BioSQLSequenceDB seqDB) {
	this.seq = null;
	this.seqDB = seqDB;
    }

    public void startFeature(Feature.Template templ)
	throws ParseException
    {
	try {
	    BioSQLFeature newFeature = _realizeFeature(templ);
	    stack.add(newFeature);
	} catch (BioException ex) {
	    throw new ParseException(ex, "Couldn't realize feature");
	} 
    }
    
    protected abstract void deliverTopLevelFeature(Feature f)
        throws ParseException, ChangeVetoException;

    public void endFeature()
	throws ParseException
    {
	if (stack.size() > 0) {
	    BioSQLFeature f = (BioSQLFeature) stack.remove(stack.size() - 1);
	    BioSQLFeature stackTop = getCurrent();
	    try {
		if (stackTop == null) {
		    deliverTopLevelFeature(seqDB.canonicalizeFeature(f, f._getInternalID()));
		} else {
		    stackTop._addFeature(f);
		}
	    } catch (ChangeVetoException ex) {
		throw new BioError(ex);
	    }
	} else {
	    throw new ParseException("start/end feature messages don't match");
	}
    }

    public void addSequenceProperty(Object key, Object value)
        throws ParseException
    {
	if ("_biosql_internal.bioentry_id".equals(key)) {
	    if (seq != null) {
		throw new ParseException("Attempting to set the sequence when it's already known!");
	    }
	    Integer bid = (Integer) value;
	    try {
		seq = seqDB.getSequence(null, bid.intValue());
	    } catch (Exception ex) {
		throw new ParseException("Non-existant sequence!");
	    }
	}
    }

    public void addFeatureProperty(Object key, Object value)
	throws ParseException
    {
	if ("_biosql_internal.feature_id".equals(key)) {
	    Integer fid = (Integer) value;
	    getCurrent()._setInternalID(fid.intValue());
	    // getCurrent()._setAnnotation(new BioSQLFeatureAnnotation(seqDB, fid.intValue()));
	} else if ("_biosql_internal.feature_id".equals(key)) {
	    Integer pid = (Integer) value;
	    getCurrent().setParentID(pid.intValue());
	} else if ("_biosql_internal.hint_childfree".equals(key)) {
	    getCurrent().hintChildFree();
	}
    }

    private BioSQLFeature getCurrent() {
	if (stack.size() > 0) {
	    return (BioSQLFeature) stack.get(stack.size() - 1);
	} else {
	    return null;
	}
    }

    private BioSQLFeature _realizeFeature(Feature.Template templ)
        throws BioException
    {
	if (templ instanceof StrandedFeature.Template && seq.getAlphabet() == DNATools.getDNA()) {
	    return new BioSQLStrandedFeature(seq, (StrandedFeature.Template) templ);
	} else {
	    return new BioSQLFeature(seq, templ);
	}
    }
}
