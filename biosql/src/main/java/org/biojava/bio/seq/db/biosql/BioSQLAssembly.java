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
import org.biojava.bio.BioError;
import org.biojava.bio.BioException;
import org.biojava.bio.BioRuntimeException;
import org.biojava.bio.seq.ComponentFeature;
import org.biojava.bio.seq.Feature;
import org.biojava.bio.seq.FeatureFilter;
import org.biojava.bio.seq.FeatureHolder;
import org.biojava.bio.seq.MergeFeatureHolder;
import org.biojava.bio.seq.RealizingFeatureHolder;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SimpleFeatureHolder;
import org.biojava.bio.seq.impl.AssembledSymbolList;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.Edit;
import org.biojava.bio.symbol.RangeLocation;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.utils.AssertionFailure;
import org.biojava.utils.ChangeListener;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;

/**
 * Sequence keyed off a BioSQL biosequence.
 *
 * @author Thomas Down
 * @author Matthew Pocock
 * @author David Huen
 * @author Richard Holland
  * @deprecated Use hibernate and org.biojavax.bio.db.*
 * @since 1.3
 */

class BioSQLAssembly
  implements
    Sequence,
    RealizingFeatureHolder,
    BioSQLSequenceI
{
    private BioSQLSequenceDB      seqDB;
    private String                name;
    private int                   assembly_id;
    private int                   bioentry_id;
    private Annotation            annotation;
    private Alphabet              alphabet;
    private int                   length;

    private RealizingFeatureHolder features;
    private SimpleFeatureHolder    componentFeatures;
    private MergeFeatureHolder     allFeatures;

    public BioSQLSequenceDB getSequenceDB() {
	return seqDB;
    }

    public int getBioEntryID() {
	return bioentry_id;
    }

    BioSQLAssembly(BioSQLSequenceDB seqDB,
		   String name,
		   int bioentry_id,
		   int assembly_id,
		   String alphaName,
		   int length)
	throws BioException
    {
	this.seqDB = seqDB;
	this.name = name;
	this.bioentry_id = bioentry_id;
	this.assembly_id = assembly_id;
	this.length = length;

	try {
	    this.alphabet = AlphabetManager.alphabetForName(alphaName.toUpperCase());
	} catch (NoSuchElementException ex) {
	    throw new BioException("Can't load sequence with unknown alphabet " + alphaName, ex);
	}

	features = new BioSQLAllFeatures(this, seqDB, bioentry_id);
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
        throws ChangeVetoException 
    {
	throw new ChangeVetoException("Can't edit sequence in BioSQL -- or at least not yet...");
    }    

    protected synchronized SimpleFeatureHolder getComponentFeatures() {
	if (componentFeatures == null) {
	    componentFeatures = new SimpleFeatureHolder();

          Connection conn = null;
	    try {
		conn = seqDB.getDataSource().getConnection();
		
		PreparedStatement get_assembly = conn.prepareStatement("select assembly_fragment_id, fragment_name, assembly_start, assembly_end, fragment_start, fragment_end, strand " +
								       "  from assembly_fragment " +
								       " where assembly_id = ?");
		get_assembly.setInt(1, assembly_id);
		ResultSet rs = get_assembly.executeQuery();
		while (rs.next()) {
		    int assembly_fragment_id = rs.getInt(1);
		    String fragment_name = rs.getString(2);
		    int assembly_start = rs.getInt(3);
		    int assembly_end = rs.getInt(4);
		    int fragment_start = rs.getInt(5);
		    int fragment_end = rs.getInt(6);

		    ComponentFeature.Template temp = new ComponentFeature.Template();
		    temp.type = "component";
		    temp.source = "biosql";
		    temp.location = new RangeLocation(assembly_start, assembly_end);
		    temp.componentSequenceName = fragment_name;
		    temp.componentLocation = new RangeLocation(fragment_start, fragment_end);
		    componentFeatures.addFeature(new BioSQLComponentFeature(seqDB, this, temp, assembly_fragment_id));
		}
                rs.close();
                get_assembly.close();
		conn.close();
	    } catch (SQLException ex) {
            if (conn!=null) try {conn.close();} catch (SQLException ex3) {}
		throw new BioRuntimeException("Error fetching assembly data", ex);
	    } catch (ChangeVetoException ex) {
            if (conn!=null) try {conn.close();} catch (SQLException ex3) {}
		throw new AssertionFailure("Assertion failure: Couldn't modify private featureholder", ex);
	    }
	}

	return componentFeatures;
    }

    private AssembledSymbolList symbols;

    protected synchronized SymbolList getSymbols() {
	if (symbols == null) {
	    FeatureHolder components = getComponentFeatures();
	    symbols = new AssembledSymbolList();
	    symbols.setLength(length);
	    for (Iterator i = components.features(); i.hasNext(); ) {
		ComponentFeature cf = (ComponentFeature) i.next();
		symbols.putComponent(cf);
	    }
	}

	return symbols;
    }

    //
    // implements FeatureHolder
    //

    private RealizingFeatureHolder getFeatures() {
	return features;
    }

    private FeatureHolder getAllFeatures() {
	if (allFeatures == null) {
	    try {
		allFeatures = new MergeFeatureHolder();
		allFeatures.addFeatureHolder(getComponentFeatures());
		allFeatures.addFeatureHolder(getFeatures());
	    } catch (ChangeVetoException ex) {
		throw new BioError("Assertion failure: Couldn't modify private featureholder");
	    }
	}

	return allFeatures;
    }

    public Iterator features() {
	return getAllFeatures().features();
    }

    public int countFeatures() {
	return getAllFeatures().countFeatures();
    }

    public boolean containsFeature(Feature f) {
	return getAllFeatures().containsFeature(f);
    }

    public FeatureHolder filter(FeatureFilter ff, boolean recurse) {
        return getAllFeatures().filter(ff, recurse);
    }

    public FeatureHolder filter(FeatureFilter ff) {
        return getAllFeatures().filter(ff);
    }
    
    
    public FeatureFilter getSchema() {
        return getFeatures().getSchema();
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
