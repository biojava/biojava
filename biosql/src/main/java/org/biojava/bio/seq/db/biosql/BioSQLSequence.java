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

import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;

import org.biojava.bio.Annotation;
import org.biojava.bio.BioException;
import org.biojava.bio.BioRuntimeException;
import org.biojava.bio.seq.Feature;
import org.biojava.bio.seq.FeatureFilter;
import org.biojava.bio.seq.FeatureHolder;
import org.biojava.bio.seq.RealizingFeatureHolder;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.db.biosql.DBHelper.BioSequenceStyle;
import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.DummySymbolList;
import org.biojava.bio.symbol.Edit;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.SimpleSymbolList;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.utils.ChangeListener;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;

/**
 * Sequence keyed off a BioSQL biosequence.
 *
 * @author Thomas Down
 * @author Matthew Pocock
 * @author Simon Foote (modifications for schema version 1.0)
 * @author Richard Holland
 * @deprecated Use hibernate and org.biojavax.bio.db.*
 * @since 1.3
 */

class BioSQLSequence
  implements
    Sequence,
    RealizingFeatureHolder,
    BioSQLSequenceI
{
    private BioSQLSequenceDB seqDB;
    private String name;
    private int bioentry_id = -1;
    private Annotation annotation;
    private Alphabet alphabet;
    private RealizingFeatureHolder features;
    private SymbolList symbols;
    private int length;

    public BioSQLSequenceDB getSequenceDB() {
	return seqDB;
    }

    public int getBioEntryID() {
	return bioentry_id;
    }

    BioSQLSequence(BioSQLSequenceDB seqDB,
	           String name,
	           int bioentry_id,
                   String alphaName,
                   int length)
	throws BioException
    {
	this.seqDB = seqDB;
	this.name = name;
	this.bioentry_id = bioentry_id;
        this.length = length;

	try {
	    this.alphabet = AlphabetManager.alphabetForName(alphaName.toUpperCase());
	} catch (NoSuchElementException ex) {
	    throw new BioException("Can't load sequence with unknown alphabet " + alphaName, ex);
	}

	// features = new BioEntryFeatureSet(this, seqDB, bioentry_id);
    }

    public String getName() {
	return name;
    }

    public String getURN() {
	return name;
    }

    //
    // implements Annotatable
    //

    public Annotation getAnnotation() {
	if (annotation == null) {
	    annotation = new BioSQLSequenceAnnotation(seqDB, bioentry_id);
	}

	return annotation;
    }

    //
    // implements SymbolList
    //

    public Alphabet getAlphabet() {
	return alphabet;
    }

    public int length() {
        if (length >= 0) {
            return length;
        } else {
	    return getSymbols().length();
        }
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
        throws ChangeVetoException 
    {
	throw new ChangeVetoException("Can't edit sequence in BioSQL -- or at least not yet...");
    }    

    protected synchronized SymbolList getSymbols() {
	if (symbols == null) {
	    Connection conn = null;
	    try {
		conn = seqDB.getDataSource().getConnection();
                
                DBHelper dh = DBHelper.getDBHelper(conn);
                BioSequenceStyle bs = dh.getBioSequenceStyle();
		                
		PreparedStatement get_symbols = conn.prepareStatement("select seq " +
								      "from   biosequence " +
								      "where  bioentry_id = ?");
		get_symbols.setInt(1, bioentry_id);
		ResultSet rs = get_symbols.executeQuery();
		String seqString = null;
		if (rs.next()) {
                    
                    if (bs==DBHelper.BIOSEQUENCE_ORACLECLOB) {
                        OracleDBHelper odh = (OracleDBHelper)dh;
                        seqString = odh.clobToString(conn, rs, 1);
                    } else { // BIOSEQUENCE_GENERIC
                        seqString = rs.getString(1);  // FIXME should do something stream-y
                    }
                    
                    if (rs.wasNull()) {
                        seqString = null;
                    }
		}
                rs.close();
		get_symbols.close();

		conn.close();

		if (seqString != null) {
		    try {
			Alphabet alpha = getAlphabet();
			SymbolTokenization toke = alpha.getTokenization("token");
			symbols = new SimpleSymbolList(toke, seqString);
		    } catch (Exception ex) {
			throw new BioRuntimeException("Couldn't parse tokenized symbols", ex);
		    }
		} else {
                    if (! (length >= 0)) {
                        throw new BioRuntimeException("Length not available from database");
                    }
                                
		    symbols = new DummySymbolList((FiniteAlphabet) alphabet, length);
		}
	    } catch (SQLException ex) {
            if (conn!=null) try {conn.close();} catch (SQLException ex3) {}
		throw new BioRuntimeException("Unknown error getting symbols from BioSQL.  Oh dear.", ex);
	    }
	}

	return symbols;
    }

    //
    // implements FeatureHolder
    //

    private RealizingFeatureHolder getFeatures() {
	if (features == null) {
	    if (length() < 5000000) {
		features = new BioSQLAllFeatures(this, seqDB, bioentry_id);
	    } else {
		features = new BioSQLTiledFeatures(this, seqDB, bioentry_id, 1000000);
	    }
	}
	return features;
    }

    public FeatureFilter getSchema() {
        return getFeatures().getSchema();
    }
    
    public Iterator features() {
	return getFeatures().features();
    }

    public int countFeatures() {
	return getFeatures().countFeatures();
    }

    public boolean containsFeature(Feature f) {
	return getFeatures().containsFeature(f);
    }

    public FeatureHolder filter(FeatureFilter ff) {
        return getFeatures().filter(ff);
    }
    
    public FeatureHolder filter(FeatureFilter ff, boolean recurse) {
        return getFeatures().filter(ff, recurse);
    }

    public Feature createFeature(Feature.Template ft)
        throws ChangeVetoException, BioException
    {
	return getFeatures().createFeature(ft);
    }

    public void removeFeature(Feature f)
        throws ChangeVetoException, BioException
    {
	getFeatures().removeFeature(f);
    }

    public Feature realizeFeature(FeatureHolder parent, Feature.Template templ)
        throws BioException
    {
	return getFeatures().realizeFeature(parent, templ);
    }

    public void persistFeature(Feature f, int parent_id)
        throws BioException
    {
	seqDB.getFeaturesSQL().persistFeature(f, parent_id, bioentry_id);
    }

    public void addChangeListener(ChangeListener cl) {
	addChangeListener(cl, ChangeType.UNKNOWN);
    }
    
    public void addChangeListener(ChangeListener cl, ChangeType ct) {
	getSequenceDB().getEntryChangeHub().addListener(bioentry_id, cl, ct);
    }

    public void removeChangeListener(ChangeListener cl) {
	removeChangeListener(cl, ChangeType.UNKNOWN);
    }

    public void removeChangeListener(ChangeListener cl, ChangeType ct) {
	getSequenceDB().getEntryChangeHub().removeListener(bioentry_id, cl, ct);
    }

    public boolean isUnchanging(ChangeType ct) {
        return false;
    }
}
