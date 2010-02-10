/* -*- c-basic-offset: 4; indent-tabs-mode: nil -*- */
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
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.Set;

import org.biojava.bio.Annotation;
import org.biojava.bio.BioError;
import org.biojava.bio.BioRuntimeException;
import org.biojava.bio.SmallAnnotation;
import org.biojava.bio.seq.io.OrganismParser;
import org.biojava.bio.taxa.Taxon;
import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeListener;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;

/**
 * Annotation keyed off a BioSQL comment table
 *
 * @author Thomas Down
 * @author Matthew Pocock
 * @author Len Trigg
 * @author Richard Holland
 * @deprecated Use hibernate and org.biojavax.bio.db.*
 * @since 1.3
 */
class BioSQLSequenceAnnotation implements Annotation {

    private BioSQLSequenceDB seqDB;
    private int bioentry_id;
    private Annotation underlyingAnnotation;

    int getBioentryID() {
	return bioentry_id;
    }

    BioSQLSequenceAnnotation(BioSQLSequenceDB seqDB,
			     int bioentry_id)
    {
	this.seqDB = seqDB;
	this.bioentry_id = bioentry_id;
    }


    private void initAnnotations() {
      Connection conn = null;
	try {
	    conn = seqDB.getDataSource().getConnection();
	    underlyingAnnotation = new SmallAnnotation();
            underlyingAnnotation.setProperty("bioentry_id", new Integer(bioentry_id));

	    //
	    // Handle all the hacky special cases first
	    //
            initTaxon(conn);

	    //
	    // General-purpose tagvalue data.
	    //
	    if (seqDB.isBioentryPropertySupported()) {
		PreparedStatement get_properties = conn.prepareStatement(
			"select term.name as qn, bioentry_qualifier_value.value " + /*, bioentry_property.property_rank as rank " + */
			"  from bioentry_qualifier_value, term " +
			" where bioentry_qualifier_value.bioentry_id = ? " +
			"   and term.term_id = bioentry_qualifier_value.term_id " /* + */
			/* " order by qn, rank" */);
		get_properties.setInt(1, bioentry_id);
		ResultSet rs = get_properties.executeQuery();
		while (rs.next()) {
		    String key = rs.getString(1).trim();   // HACK due to stupid schema change
		    String value = rs.getString(2);
		    if (underlyingAnnotation.containsProperty(key)) {
			Object current = underlyingAnnotation.getProperty(key);
			Collection coll;
			if (! (current instanceof Collection)) {
			    coll = new ArrayList();
			    coll.add(current);
			    underlyingAnnotation.setProperty(key, coll);
			} else {
			    coll = (Collection) current;
			}
			coll.add(value);
		    } else {
			underlyingAnnotation.setProperty(key, value);
		    }
		}
                rs.close();
                get_properties.close();
	    }
	    
	    conn.close();
	} catch (SQLException ex) {
          if (conn!=null) try {conn.close();} catch (SQLException ex3) {}
	    throw new BioRuntimeException("Error fetching annotations", ex);
	} catch (ChangeVetoException ex) {
          if (conn!=null) try {conn.close();} catch (SQLException ex3) {}
	    throw new BioError(ex);
	}
    }


    /**
     * Initialize the Taxon annotation. Reads from the taxon and
     * taxon_name table and stuffs the results into the annotation
     * object under key OrganismParser.PROPERTY_ORGANISM.
     *
     * @param conn the current <code>Connection</code>.
     */
    private void initTaxon(Connection conn) {
        PreparedStatement statement = null;
        ResultSet rs = null;
        try {
            int taxon_id = 0;
            statement = conn.prepareStatement("select taxon_id " 
                                              + "from bioentry " 
                                              + "where bioentry_id = ? ");
            statement.setInt(1, bioentry_id);
            rs = statement.executeQuery();
            if (rs.next()) {
                taxon_id = rs.getInt(1);
                if (rs.wasNull()) {
                    taxon_id = 0;
                }
            }
            rs.close();
            statement.close();

            if (taxon_id != 0) {
                // We are expected to get a taxon structure
                Taxon taxon = TaxonSQL.getDBTaxon(conn, taxon_id);
                underlyingAnnotation.setProperty(OrganismParser.PROPERTY_ORGANISM, taxon);
            }
	} catch (ChangeVetoException ex) {
	    throw new BioRuntimeException("Couldn't manipulate in-memory taxonomy", ex);
	} catch (SQLException ex) {
	    throw new BioRuntimeException("Error fetching taxonomy annotations", ex);
        } finally {
            TaxonSQL.attemptClose(rs);
            TaxonSQL.attemptClose(statement);
        }
    }


    public Object getProperty(Object key) {
	if (underlyingAnnotation == null) {
	    initAnnotations();
	}

	return underlyingAnnotation.getProperty(key);
    }

    public void setProperty(Object key, Object value)
        throws ChangeVetoException
    {
        BioSQLEntryAnnotationChangeHub entryAnnotationHub = seqDB.getEntryAnnotationChangeHub();
        synchronized (entryAnnotationHub) {
            ChangeEvent cev = new ChangeEvent(this, Annotation.PROPERTY, key);
            entryAnnotationHub.firePreChange(cev);
            _setProperty(key, value);
            entryAnnotationHub.firePostChange(cev);
        }
    }

    private void _setProperty(Object key, Object value) 
        throws ChangeVetoException
    {
        persistProperty(key, value);
        if (underlyingAnnotation != null) {
            underlyingAnnotation.setProperty(key, value);
        }
    }
    
    public void removeProperty(Object key)
        throws ChangeVetoException
    {
        if (underlyingAnnotation == null) {
            initAnnotations();
        }
        if (!underlyingAnnotation.containsProperty(key)) {
            throw new NoSuchElementException("Annotation doesn't contain property " + key.toString());
        }
        
        BioSQLEntryAnnotationChangeHub entryAnnotationHub = seqDB.getEntryAnnotationChangeHub();
        synchronized (entryAnnotationHub) {
            ChangeEvent cev = new ChangeEvent(this, Annotation.PROPERTY, key);
            entryAnnotationHub.firePreChange(cev);
            underlyingAnnotation.removeProperty(key);
            persistProperty(key, null);
            entryAnnotationHub.firePostChange(cev);
        }
    }

    private void persistProperty(Object key, Object value)
        throws ChangeVetoException
    {
	Connection conn = null;
	try {
	    conn = seqDB.getDataSource().getConnection();
	    conn.setAutoCommit(false);

	    seqDB.persistBioentryProperty(conn, bioentry_id, key, value, true, false);
	    
	    conn.commit();
	    conn.close();
	} catch (SQLException ex) {
	    boolean rolledback = false;
	    if (conn != null) {
		try {
		    conn.rollback();
		    rolledback = true;
		} catch (SQLException ex2) {}
            try {conn.close();} catch (SQLException ex3) {}
	    }
	    throw new BioRuntimeException("Error adding to BioSQL tables" + 
					(rolledback ? " (rolled back successfully)" : ""), ex);
	}
    }

    public boolean containsProperty(Object key) {
	if (underlyingAnnotation == null) {
	    initAnnotations();
	}

	return underlyingAnnotation.containsProperty(key);
    }

    public Set keys() {
	if (underlyingAnnotation == null) {
	    initAnnotations();
	}

	return underlyingAnnotation.keys();
    }

    public Map asMap() {
	if (underlyingAnnotation == null) {
	    initAnnotations();
	}

	return Collections.unmodifiableMap(underlyingAnnotation.asMap());
    }

    
    public void addChangeListener(ChangeListener cl) {
	addChangeListener(cl, ChangeType.UNKNOWN);
    }
    
    public void addChangeListener(ChangeListener cl, ChangeType ct) {
	seqDB.getEntryAnnotationChangeHub().addListener(new Integer(bioentry_id), cl, ct);
    }

    public void removeChangeListener(ChangeListener cl) {
	removeChangeListener(cl, ChangeType.UNKNOWN);
    }

    public void removeChangeListener(ChangeListener cl, ChangeType ct) {
	seqDB.getEntryAnnotationChangeHub().removeListener(new Integer(bioentry_id), cl, ct);
    }

    public boolean isUnchanging(ChangeType ct) {
	return false;
    }
}
