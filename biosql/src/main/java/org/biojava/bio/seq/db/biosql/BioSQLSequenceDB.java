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
import java.io.StringReader;
import java.sql.Connection;
import java.sql.Driver;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.sql.DataSource;

import org.biojava.bio.Annotation;
import org.biojava.bio.BioError;
import org.biojava.bio.BioException;
import org.biojava.bio.BioRuntimeException;
import org.biojava.bio.seq.Feature;
import org.biojava.bio.seq.FeatureFilter;
import org.biojava.bio.seq.FeatureHolder;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.SimpleFeatureHolder;
import org.biojava.bio.seq.db.IDMaker;
import org.biojava.bio.seq.db.IllegalIDException;
import org.biojava.bio.seq.db.SequenceDB;
import org.biojava.bio.seq.db.biosql.DBHelper.BioSequenceStyle;
import org.biojava.bio.seq.io.OrganismParser;
import org.biojava.bio.seq.io.ParseException;
import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.taxa.Taxon;
import org.biojava.ontology.Ontology;
import org.biojava.utils.AbstractChangeable;
import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.JDBCPooledDataSource;
import org.biojava.utils.cache.Cache;
import org.biojava.utils.cache.FixedSizeCache;
import org.biojava.utils.cache.WeakValueHashMap;

/**
 * SequenceDB keyed off a BioSQL database.  This is an almost-complete
 * implementation of the BioJava Sequence, SequenceDB, and Feature interfaces,
 * and can be used in a wide range of applications.
 *
 * Note: It now uses BioSQL schema version 1.0 (Post Singapore)
 * All previous versions are no longer supported.
 *
 * @author Thomas Down
 * @author Matthew Pocock
 * @author Simon Foote
 * @author Len Trigg
 * @author Mark Schreiber
 * @author Richard Holland
 * @deprecated Use hibernate and org.biojavax.bio.db.*
 * @since 1.3
 */
public class BioSQLSequenceDB extends AbstractChangeable implements SequenceDB {

    private DataSource dataSource;
    private int dbid = -1;
    private String name;
    private IDMaker idmaker = new IDMaker.ByName();
    private WeakValueHashMap sequencesByName = new WeakValueHashMap();
    private WeakValueHashMap sequencesByID = new WeakValueHashMap();
    private DBHelper helper;
    private FeaturesSQL featuresSQL;
    private OntologySQL ontologySQL;
//    private BioSQLChangeHub changeHub;
    private BioSQLEntryChangeHub entryChangeHub;
    private BioSQLEntryAnnotationChangeHub entryAnnotationChangeHub;
    private BioSQLFeatureChangeHub featureChangeHub;
    private BioSQLFeatureAnnotationChangeHub featureAnnotationChangeHub;
    private WeakValueHashMap featuresByID = new WeakValueHashMap();
    private Cache tileCache = new FixedSizeCache(10);

    DataSource getDataSource() {
        return dataSource;
    }

    DBHelper getDBHelper() {
        return helper;
    }

    FeaturesSQL getFeaturesSQL() {
        return featuresSQL;
    }
/*
    BioSQLChangeHub getChangeHub() {
        return changeHub;
    }
*/
    BioSQLEntryChangeHub getEntryChangeHub() {
        return entryChangeHub;
    }

    BioSQLEntryAnnotationChangeHub getEntryAnnotationChangeHub() {
        return entryAnnotationChangeHub;
    }

    BioSQLFeatureChangeHub getFeatureChangeHub() {
        return featureChangeHub;
    }

    BioSQLFeatureAnnotationChangeHub getFeatureAnnotationChangeHub() {
        return featureAnnotationChangeHub;
    }

    /**
     * Connect to a BioSQL database.
     *
     * @param dbDriver A JDBC database driver.  For example, <code>com.jdbc.mysql.Driver</code>
     * @param dbURL A JDBC database URL.  For example, <code>jdbc:postgresql://localhost/thomasd_biosql</code>
     * @param dbUser The username to use when connecting to the database (or an empty string).
     * @param dbPass The password to use when connecting to the database (or an empty string).
     * @param biodatabase The identifier of a namespace within the physical BioSQL database.
     * @param create If the requested namespace doesn't exist, and this flag is <code>true</code>,
     *               a new namespace will be created.
     *
     * @throws BioException if an error occurs communicating with the database
     */

    public BioSQLSequenceDB(String dbDriver,
                            String dbURL,
                            String dbUser,
                            String dbPass,
                            String biodatabase,
                            boolean create)
        throws BioException {

        try {
            dataSource = JDBCPooledDataSource.getDataSource(dbDriver, dbURL, dbUser, dbPass);
        } catch (Exception ex) {
            throw new BioException("Error getting datasource", ex);
        }
        this.initDb(biodatabase, create);
    }

    /**
     * Connect to a BioSQL database.
     *
     * @param dbURL A JDBC database URL.  For example, <code>jdbc:postgresql://localhost/thomasd_biosql</code>
     * @param dbUser The username to use when connecting to the database (or an empty string).
     * @param dbPass The password to use when connecting to the database (or an empty string).
     * @param biodatabase The identifier of a namespace within the physical BioSQL database.
     * @param create If the requested namespace doesn't exist, and this flag is <code>true</code>,
     *               a new namespace will be created.
     *
     * @throws BioException if an error occurs communicating with the database
     */

    public BioSQLSequenceDB(String dbURL,
                            String dbUser,
                            String dbPass,
                            String biodatabase,
                            boolean create)
        throws BioException {

        try {
            Driver drv = DriverManager.getDriver(dbURL);
            dataSource = JDBCPooledDataSource.getDataSource(drv.getClass().getName(), dbURL, dbUser, dbPass);
        } catch (Exception ex) {
            throw new BioException("Error getting datasource", ex);
        }
        this.initDb(biodatabase, create);
    }

    public BioSQLSequenceDB(DataSource ds,
                            String biodatabase,
                            boolean create)
        throws BioException {

        dataSource = ds;
        this.initDb(biodatabase, create);
    }

    void initDb(String biodatabase, boolean create)
        throws BioException {

        // Create helpers
        entryChangeHub = new BioSQLEntryChangeHub(this);
        entryAnnotationChangeHub = new BioSQLEntryAnnotationChangeHub(this, entryChangeHub);
        featureChangeHub = new BioSQLFeatureChangeHub(this, entryChangeHub);
        featureAnnotationChangeHub = new BioSQLFeatureAnnotationChangeHub(this, featureChangeHub);

	  Connection conn = null;
        try {
            conn = dataSource.getConnection();
            conn.setAutoCommit(false);

            // DBHelper needs to be initialized before checks and ontologies are created
            helper = DBHelper.getDBHelper(conn);

            // Check that BioSQL database schema is post-Singapore
            if (! isDbSchemaSupported()) {
		    try {conn.close();} catch (SQLException ex3) {}
                throw new BioException("This database appears to be an old (pre-Singapore) BioSQL."
                                       + " If you need to access it, try an older BioJava snapshot (1.3pre1 or earlier)");
            }

            if (! isBioentryPropertySupported()) {
		    try {conn.close();} catch (SQLException ex3) {}
                throw new BioException("This database appears to be an old (pre-Cape-Town) BioSQL."
                                       + " If you need to access it, try an older BioJava snapshot");
            }

            // Create adapters
            featuresSQL = new FeaturesSQL(this);
            try {
                //ontologySQL = new OntologySQL(dataSource, helper);
                ontologySQL = OntologySQL.getOntologySQL(dataSource, helper);
            } catch (SQLException ex) {
		    try {conn.close();} catch (SQLException ex3) {}
                throw new BioException("Error accessing ontologies", ex);
            }

            PreparedStatement getID = conn.prepareStatement("select biodatabase_id from biodatabase where name = ?");
            getID.setString(1, biodatabase);
            ResultSet rs = getID.executeQuery();
            if (rs.next()) {
                dbid = rs.getInt(1);
                name = biodatabase;
                rs.close();
                getID.close();
                conn.close();
                
            } else {
                rs.close();
                getID.close();
            
                if (create) {
                    PreparedStatement createdb = conn.prepareStatement(
                            "insert into biodatabase (name) values ( ? )");
                    createdb.setString(1, biodatabase);
                    createdb.executeUpdate();
                    conn.commit();
                    createdb.close();
		    conn.close();
                    dbid = getDBHelper().getInsertID(conn, "biodatabase", "biodatabase_id");
                } else {
		    conn.close();
                    throw new BioException("Biodatabase " + biodatabase + " doesn't exist");
                }
            }
        } catch (SQLException ex) {
		if (conn!=null) try {conn.close();} catch (SQLException ex3) {}
            throw new BioException("Error connecting to BioSQL database: " + ex.getMessage(), ex);
        }
    }

    public Ontology createOntology(String name, String description)
        throws Exception
    {
        return ontologySQL.createOntology(name, description);
    }

    public Ontology getOntology(String name)
        throws Exception
    {
        return ontologySQL.getOntology(name);
    }

    public Ontology addOntology(Ontology onto)
        throws Exception
    {
        return ontologySQL.addOntology(onto);
    }

    public String getName() {
        return name;
    }

    public void createDummySequence(String id,
                                    Alphabet alphabet,
                                    int length)
        throws ChangeVetoException, BioException
    {
        synchronized (this) {
          ChangeEvent cev = new ChangeEvent(this, SequenceDB.SEQUENCES, null);
          firePreChangeEvent(cev);
          _createDummySequence(id, alphabet, length);
          firePostChangeEvent(cev);
        }
    }

    private void _createDummySequence(String id,
                                      Alphabet seqAlpha,
                                      int length)
        throws ChangeVetoException, BioException
    {
        int version = 1;

        Connection conn = null;
        try {
            conn = dataSource.getConnection();
            conn.setAutoCommit(false);

            PreparedStatement create_bioentry = conn.prepareStatement(
                                                                      "insert into bioentry " +
                                                                      "(biodatabase_id, name, accession, version, division) " +
                                                                      "values (?, ?, ?, ?, ?)");
            create_bioentry.setInt(1, dbid);
            create_bioentry.setString(2, id);
            create_bioentry.setString(3, id);
            create_bioentry.setInt(4, version);
            create_bioentry.setString(5, "?");
            create_bioentry.executeUpdate();
            create_bioentry.close();

            int bioentry_id = getDBHelper().getInsertID(conn, "bioentry", "bioentry_id");

            PreparedStatement create_dummy = conn.prepareStatement("insert into biosequence " +
                                                                   "       (bioentry_id, version, alphabet, length) " +
                                                                   "values (?, ?, ?, ?)");
            create_dummy.setInt(1, bioentry_id);
            create_dummy.setInt(2, version);
            create_dummy.setString(3, seqAlpha.getName());
            create_dummy.setInt(4, length);
            create_dummy.executeUpdate();
            create_dummy.close();
            //int dummy_id = getDBHelper().getInsertID(conn, "biosequence", "biosequence_id");
            
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
            throw new BioRuntimeException("Error adding dummy sequence" +
                                          (rolledback ? " (rolled back successfully)" : ""), ex);
        }
    }

    public void addSequence(Sequence seq)
        throws ChangeVetoException, BioException
    {
        synchronized (this) {
            ChangeEvent cev = new ChangeEvent(this, SequenceDB.SEQUENCES, seq);
            firePreChangeEvent(cev);
            _addSequence(seq);
            firePostChangeEvent(cev);
        }
    }

    private void _addSequence(Sequence seq)
        throws ChangeVetoException, BioException {
        String seqName = idmaker.calcID(seq);
        int version = 1;

        Alphabet seqAlpha = seq.getAlphabet();
        SymbolTokenization seqToke;
        try {
            seqToke = seqAlpha.getTokenization("token");
        } catch (Exception ex) {
            throw new BioException("Can't store sequences in BioSQL unless they can be sensibly tokenized/detokenized", ex);
        }

        Connection conn = null;
        try {
            conn = dataSource.getConnection();
            conn.setAutoCommit(false);
            //ResultSet rs;

            //
            // we will need this annotation bundle for various things
            //

            Annotation ann = seq.getAnnotation();

            PreparedStatement create_bioentry = conn.prepareStatement(
                                                                      "insert into bioentry " +
                                                                      "(biodatabase_id, name, accession, version, division) " +
                                                                      "values (?, ?, ?, ?, ?)"
                                                                      );
            create_bioentry.setInt(1, dbid);
            create_bioentry.setString(2, seqName);
            create_bioentry.setString(3, seqName);
            create_bioentry.setInt(4, version);
            create_bioentry.setString(5, "?");
            create_bioentry.executeUpdate();
            create_bioentry.close();

            // System.err.println("Created bioentry");

            int bioentry_id = getDBHelper().getInsertID(conn, "bioentry", "bioentry_id");
                            
            BioSequenceStyle bs = getDBHelper().getBioSequenceStyle();

            // See if we are using CLOBs.
            if (bs==DBHelper.BIOSEQUENCE_ORACLECLOB) {
                PreparedStatement create_biosequence = conn.prepareStatement("insert into biosequence " +
                                                                         "(bioentry_id, version, length, seq, alphabet) " +
                                                                         "values (?, ?, ?, empty_clob(), ?)");
                create_biosequence.setInt(1, bioentry_id);
                create_biosequence.setInt(2, version);
                create_biosequence.setInt(3, seq.length());
                String seqstr = seqToke.tokenizeSymbolList(seq);

                create_biosequence.setString(4, seqAlpha.getName());
                create_biosequence.executeUpdate();
                create_biosequence.close();
                    
                // Now retrieve and update
                PreparedStatement retrieve_biosequence = conn.prepareStatement("select seq from biosequence " +
                                                                        "where bioentry_id = ? for update");
                retrieve_biosequence.setInt(1, bioentry_id);
                ResultSet rs = retrieve_biosequence.executeQuery();
                if (!rs.next()) throw new BioRuntimeException("Could not read newly inserted sequence!");
                    
                OracleDBHelper odh = (OracleDBHelper)getDBHelper();
                odh.stringToClob(conn, rs, 1, seqstr);

                rs.close();
                retrieve_biosequence.close();
                    
            } else { // BIOSEQUENCE_GENERIC
                PreparedStatement create_biosequence = conn.prepareStatement("insert into biosequence " +
                                                                         "(bioentry_id, version, length, seq, alphabet) " +
                                                                         "values (?, ?, ?, ?, ?)");
                create_biosequence.setInt(1, bioentry_id);
                create_biosequence.setInt(2, version);
                create_biosequence.setInt(3, seq.length());
                String seqstr = seqToke.tokenizeSymbolList(seq);
            
                create_biosequence.setCharacterStream(4, new StringReader(seqstr), seqstr.length());

                create_biosequence.setString(5, seqAlpha.getName());
                create_biosequence.executeUpdate();
                create_biosequence.close();
            }
            

            // System.err.println("Stored sequence");

            //
            // Store the features
            //

            FeatureHolder features = seq;
            int num = features.countFeatures();
            if (!isHierarchySupported()) {
                features = features.filter(FeatureFilter.all, true);
                if (features.countFeatures() != num) {
                    System.err.println("*** Warning: feature hierarchy was lost when adding sequence to BioSQL");
                }
            }
            getFeaturesSQL().persistFeatures(conn, bioentry_id, features);

            // System.err.println("Stored features");

            //
            // Store generic properties
            //

            for (Iterator i = ann.asMap().entrySet().iterator(); i.hasNext(); ) {
                Map.Entry me = (Map.Entry) i.next();
                Object key = me.getKey();
                Object value = me.getValue();
                persistBioentryProperty(conn, bioentry_id, key, value, false, true);
            }
                    
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
            throw new BioRuntimeException(
                                          "Error adding sequence: " + seq.getName() +
                                          (rolledback ? " (rolled back successfully)" : ""), ex);
        }
    }


    public Sequence getSequence(String id)
        throws BioException {
        return getSequence(id, -1);
    }
    
    public Sequence getSequence(int bioentry_id)
	    throws BioException {
      return getSequence(null, bioentry_id);
    }

    Sequence getSequence(String id, int bioentry_id)
        throws BioException {
        Sequence seq = null;
        if (id != null) {
            seq = (Sequence) sequencesByName.get(id);
        } else if (bioentry_id >= 0) {
            seq = (Sequence) sequencesByID.get(new Integer(bioentry_id));
        } else {
            throw new BioError("Neither a name nor an internal ID was supplied");
        }

        if (seq != null) {
            return seq;
        }

	  Connection conn = null;
        try {
            conn = dataSource.getConnection();

            if (bioentry_id < 0) {
                PreparedStatement get_bioentry = conn.prepareStatement("select bioentry.bioentry_id " +
                                                                       "from bioentry " +
                                                                       "where bioentry.accession = ? and " +
                                                                       "      bioentry.biodatabase_id = ?");
                get_bioentry.setString(1, id);
                get_bioentry.setInt(2, dbid);
                ResultSet rs = get_bioentry.executeQuery();
                if (rs.next()) {
                    bioentry_id = rs.getInt(1);
                }
                rs.close();
                get_bioentry.close();

                if (bioentry_id < 0) {
                    conn.close();
                    throw new IllegalIDException("No bioentry with accession " + id);
                }
            } else {
                PreparedStatement get_accession = conn.prepareStatement("select bioentry.accession from bioentry where bioentry.bioentry_id = ? and bioentry.biodatabase_id = ?");
                get_accession.setInt(1, bioentry_id);
                get_accession.setInt(2, dbid);
                ResultSet rs = get_accession.executeQuery();
                if (rs.next()) {
                    id = rs.getString(1);
                }
                rs.close();
                get_accession.close();

                if (id == null) {
                    conn.close();
                    throw new IllegalIDException("No bioentry with internal ID " + bioentry_id);
                }
            }

            if (seq == null) {
                PreparedStatement get_biosequence = conn.prepareStatement("select alphabet, length " +
                                                                          "from   biosequence " +
                                                                          "where  bioentry_id = ?");
                get_biosequence.setInt(1, bioentry_id);
                ResultSet rs = get_biosequence.executeQuery();
                if (rs.next()) {
                    // UC conversion required for lower-case alphabets from bioperl. 
                    // This is because BioSQL accepts both UC and LC.
                    String molecule = rs.getString(1).toUpperCase(); 
                    int length = rs.getInt(2);
                    if (rs.wasNull()) {
                        length = -1;
                    }
                    seq = new BioSQLSequence(this, id, bioentry_id, molecule, length);
                }
                rs.close();
                get_biosequence.close();
            }

            if (seq == null && isAssemblySupported()) {
                PreparedStatement get_assembly = conn.prepareStatement("select assembly_id, length, molecule " +
                                                                       "from   assembly " +
                                                                       "where  bioentry_id = ?");
                get_assembly.setInt(1, bioentry_id);
                ResultSet rs = get_assembly.executeQuery();
                if (rs.next()) {
                    int assembly_id = rs.getInt(1);
                    int length = rs.getInt(2);
                    String molecule = rs.getString(3);
                    seq = new BioSQLAssembly(this, id, bioentry_id, assembly_id, molecule, length);
                }
                rs.close();
                get_assembly.close();
            }

            conn.close();

            if (seq != null) {
                sequencesByName.put(id, seq);
                sequencesByID.put(new Integer(bioentry_id), seq);
                return seq;
            }
        } catch (SQLException ex) {
            if (conn!=null) try {conn.close();} catch (SQLException ex3) {}
            throw new BioException("Error accessing BioSQL tables", ex);
        }

        throw new BioException("BioEntry " + id + " exists with unknown sequence type");
    }

    public void removeSequence(String id)
        throws ChangeVetoException, BioException {

        synchronized (this) {
            ChangeEvent cev = new ChangeEvent(this, SequenceDB.SEQUENCES, null);
            firePreChangeEvent(cev);
            _removeSequence(id);
            firePostChangeEvent(cev);
        }
    }

    private void _removeSequence(String id)
        throws BioException, ChangeVetoException
    {

        Sequence seq = (Sequence) sequencesByName.get(id);

        if (seq != null) {
            seq = null;  // Don't want to be holding the reference ourselves!
            try {
                Thread.sleep(100L);
                System.gc();
            } catch (Exception ex) {
                ex.printStackTrace();
            }
            seq = (Sequence) sequencesByName.get(id);
            if (seq != null) {
                throw new BioException("There are still references to sequence with ID " + id + " from this database.");
            }
        }

        Connection conn = null;
        try {
            conn = dataSource.getConnection();
            conn.setAutoCommit(false);
            PreparedStatement get_sequence = conn.prepareStatement("select bioentry.bioentry_id " +
                                                                   "from bioentry where accession = ? and biodatabase_id = ?"
                                                                   );
            get_sequence.setString(1, id);
            get_sequence.setInt(2, dbid);
            ResultSet rs = get_sequence.executeQuery();
            boolean exists;
            if ((exists = rs.next())) {

                // For MySQL4 (default is to use InnoDB tables), delete is via a CASCADE
                // so only have to delete from bioentry to get all references to that
                // bioentry deleted
                DBHelper.DeleteStyle dstyle = getDBHelper().getDeleteStyle();
                int bioentry_id = rs.getInt(1);

                if (dstyle !=  DBHelper.DELETE_MYSQL4) {

                    PreparedStatement delete_reference = conn.prepareStatement("delete from bioentry_reference where bioentry_id = ?");
                    delete_reference.setInt(1, bioentry_id);
                    delete_reference.executeUpdate();
                    delete_reference.close();

                    String commentTableName = getCommentTableName();
                    if (commentTableName != null) {
                        PreparedStatement delete_comment = conn.prepareStatement("delete from " + commentTableName
                                                                                 + " where bioentry_id = ?");
                        delete_comment.setInt(1, bioentry_id);
                        delete_comment.executeUpdate();
                        delete_comment.close();
                    }

                    PreparedStatement delete_qv = conn.prepareStatement("delete from bioentry_qualifier_value where bioentry_id = ?");
                    delete_qv.setInt(1, bioentry_id);
                    delete_qv.executeUpdate();
                    delete_qv.close();

                    ArrayList generic_ids = null;  // default delete style will cache seqfeature_id's that need to be deleted

                    PreparedStatement delete_locs;
                    if (dstyle ==  DBHelper.DELETE_POSTGRESQL) {
                        delete_locs = conn.prepareStatement("delete from location" +
                                                            " where location.seqfeature_id = seqfeature.seqfeature_id and" +
                                                            " seqfeature.bioentry_id = ?");
                        delete_locs.setInt(1, bioentry_id);
                        delete_locs.executeUpdate();
                        delete_locs.close();
                    } else {
                        delete_locs = conn.prepareStatement("delete from location where seqfeature_id = ?");

                        PreparedStatement get_seqfeats = conn.prepareStatement("select seqfeature_id"
                                                                               + " from seqfeature"
                                                                               + " where bioentry_id = ?"
                                                                               );
                        get_seqfeats.setInt(1, bioentry_id);
                        ResultSet sfids = get_seqfeats.executeQuery();
                        generic_ids = new ArrayList();
                        while (sfids.next()) {
                            int sfid = sfids.getInt(1);
                            generic_ids.add(new Integer(sfid));
                            delete_locs.setInt(1, sfid);
                            delete_locs.executeUpdate();
                        }
                        sfids.close();
                        get_seqfeats.close();
                    }
                    delete_locs.close();

                    PreparedStatement delete_fqv;
                    if (dstyle ==  DBHelper.DELETE_POSTGRESQL) {
                        delete_fqv = conn.prepareStatement("delete from seqfeature_qualifier_value" +
                                                           " where seqfeature_qualifier_value.seqfeature_id = seqfeature.seqfeature_id" +
                                                           " and seqfeature.bioentry_id = ?");
                        delete_fqv.setInt(1, bioentry_id);
                        delete_fqv.executeUpdate();
                    } else {
                        delete_fqv = conn.prepareStatement("delete from seqfeature_qualifier_value"
                                                           + " where seqfeature_qualifier_value.seqfeature_id = ?");
                        for (int i = 0; i < generic_ids.size(); i++) {
                            int sfid = ((Integer) generic_ids.get(i)).intValue();
                            delete_fqv.setInt(1, sfid);
                            delete_fqv.executeUpdate();
                        }
                    }
                    delete_fqv.close();

                    PreparedStatement delete_rel;
                    if (dstyle ==  DBHelper.DELETE_POSTGRESQL) {
                        delete_rel = conn.prepareStatement("delete from seqfeature_relationship" +
                                                           " where object_seqfeature_id = seqfeature.seqfeature_id" +
                                                           " and seqfeature.bioentry_id = ?");
                        delete_rel.setInt(1, bioentry_id);
                        delete_rel.executeUpdate();
                    } else {
                        delete_rel = conn.prepareStatement("delete from seqfeature_relationship"
                                                           + " where object_seqfeature_id = ?");
                        for (int i = 0; i < generic_ids.size(); i++) {
                            int sfid = ((Integer) generic_ids.get(i)).intValue();
                            delete_rel.setInt(1, sfid);
                            delete_rel.executeUpdate();
                        }
                    }
                    delete_rel.close();

                    PreparedStatement delete_features = conn.prepareStatement("delete from seqfeature " +
                                                                              " where bioentry_id = ?");
                    delete_features.setInt(1, bioentry_id);
                    delete_features.executeUpdate();
                    delete_features.close();

                    PreparedStatement delete_biosequence =
                        conn.prepareStatement("delete from biosequence where bioentry_id = ?");

                    delete_biosequence.setInt(1, bioentry_id);
                    delete_biosequence.executeUpdate();
                    delete_biosequence.close();
                } // End of if for non-MYSQL4 deletion

                // All DB types must delete the bioentry via its id
                // MySQL4 only needs to delete this, as it will cascade delete all references to it
                PreparedStatement delete_entry =
                    conn.prepareStatement("delete from bioentry where bioentry_id = ?");
                delete_entry.setInt(1, bioentry_id);
                int status = delete_entry.executeUpdate();
                // Check that a deletion actually occurred, if not indicate so
                if (status < 1) {
                    System.out.println("Bioentry (ID " + bioentry_id + ") failed to delete!!");
                }
                delete_entry.close();
            }

            rs.close();
            get_sequence.close();

            conn.commit();
            conn.close();

            if (!exists) {
                throw new IllegalIDException("Sequence " + id + " didn't exist");
            }
        } catch (SQLException ex) {
            boolean rolledback = false;
            if (conn != null) {
                try {
                    conn.rollback();
                    rolledback = true;
                } catch (SQLException ex2) {}
                try {conn.close();} catch (SQLException ex3) {}
            }
            throw new BioException("Error removing from BioSQL tables" +
                                   (rolledback ? " (rolled back successfully)" : ""), ex);
        }
    }

    public Set ids() {
	  Connection conn = null;
        try {
            Set _ids = new HashSet();
            conn = dataSource.getConnection();

            PreparedStatement st = conn.prepareStatement("select bioentry.accession from bioentry where bioentry.biodatabase_id = ?");
            st.setInt(1, dbid);
            ResultSet rs = st.executeQuery();
            while (rs.next()) {
                _ids.add(rs.getString(1));
            }
            rs.close();
            st.close();

            conn.close();
            return Collections.unmodifiableSet(_ids);
        } catch (SQLException ex) {
            if (conn!=null) try {conn.close();} catch (SQLException ex3) {}
            throw new BioRuntimeException("Error reading from BioSQL tables", ex);
        }
    }

    //
    // Sequence support
    //


    void persistBioentryProperty(Connection conn,
                                 int bioentry_id,
                                 Object key,
                                 Object value,
                                 boolean removeFirst,
                                 boolean silent)
        throws SQLException
    {
        // Ought to check for special-case keys. (or just wait 'til the special case
        // tables get nuked :-)
        // ex. taxon, references, dbxrefs

        if (key.equals(OrganismParser.PROPERTY_ORGANISM)) {
            int taxon_id = TaxonSQL.putTaxon(conn, getDBHelper(), (Taxon) value);
            if (taxon_id != -1) {
                PreparedStatement set_taxon = conn.prepareStatement("update bioentry set taxon_id = ? "
                                                                    + " where bioentry_id = ?");
                set_taxon.setInt(1, taxon_id);
                set_taxon.setInt(2, bioentry_id);
                set_taxon.executeUpdate();
                set_taxon.close();
            }

        } else {
            String keyString = key.toString();

            if (!isBioentryPropertySupported()) {
                if (silent) {
                    return;
                } else {
                    throw new SQLException("Can't persist this property since the bioentry_qualifier_value table isn't available");
                }
            }

            if (removeFirst) {
                int id = intern_ontology_term(conn, keyString);
                PreparedStatement remove_old_value = conn.prepareStatement("delete from bioentry_qualifier_value " +
                                                                           " where bioentry_id = ? and term_id = ?");
                remove_old_value.setInt(1, bioentry_id);
                remove_old_value.setInt(2, id);
                remove_old_value.executeUpdate();
                remove_old_value.close();
            }

            if (value != null) {
                PreparedStatement insert_new;
                if (isSPASupported()) {
                    insert_new = conn.prepareStatement("insert into bioentry_qualifier_value " +
                                                       "       (bioentry_id, term_id, value, rank) " +
                                                       "values (?, intern_ontology_term( ? ), ?, ?)");
                    if (value instanceof Collection) {
                        int cnt = 0;
                        for (Iterator i = ((Collection) value).iterator(); i.hasNext(); ) {
                            insert_new.setInt(1, bioentry_id);
                            insert_new.setString(2, keyString);
                            insert_new.setInt(4, ++cnt);
                            insert_new.setString(3, i.next().toString());
                            insert_new.executeUpdate();
                        }
                    } else {
                        insert_new.setInt(1, bioentry_id);
                        insert_new.setString(2, keyString);
                        insert_new.setInt(3, 1);
                        insert_new.setString(3, value.toString());
                        insert_new.executeUpdate();
                    }
                } else {
                    insert_new = conn.prepareStatement("insert into bioentry_qualifier_value " +
                                                       "       (bioentry_id, term_id, rank, value) " +
                                                       "values (?, ?, ?, ?)");
                    int termID = intern_ontology_term(conn, keyString);
                    if (value instanceof Collection) {
                        int cnt = 0;
                        for (Iterator i = ((Collection) value).iterator(); i.hasNext(); ) {
                            insert_new.setInt(1, bioentry_id);
                            insert_new.setInt(2, termID);
                            insert_new.setInt(3, ++cnt);
                            insert_new.setString(4, i.next().toString());
                            insert_new.executeUpdate();
                        }
                    } else {
                        insert_new.setInt(1, bioentry_id);
                        insert_new.setInt(2, termID);
                        insert_new.setInt(3, 1);
                        insert_new.setString(4, value.toString());
                        insert_new.executeUpdate();
                    }
                }
                insert_new.close();
            }
        }
    }

    /**
     * Legacy method -- will go eventually
     */

    int intern_ontology_term(Connection conn, String s)
        throws SQLException
    {
        Ontology legacy = ontologySQL.getLegacyOntology();
        String ts = s.trim();  // Hack for schema change
        if (legacy.containsTerm(ts)) {
            return ontologySQL.termID(legacy.getTerm(ts));
            // Same term but different case causes error when try to add it for MySQL
            // These hacks prevent it.  ex. genbank can have ORGANISM & organism keys
            // Removed hack as if set term name to BINARY handles case correctly
//        } else if (legacy.containsTerm(ts.toLowerCase())) {
//            return ontologySQL.termID(legacy.getTerm(ts.toLowerCase()));
//        } else if (legacy.containsTerm(ts.toUpperCase())) {
//            return ontologySQL.termID(legacy.getTerm(ts.toUpperCase()));

        } else {
            try {
                return ontologySQL.termID(legacy.createTerm(ts, ""));
            } catch (Exception ex) {
                //ex.printStackTrace();
                //System.err.println("Term: " + ts + "   " + ex.getMessage());
                throw (SQLException) new SQLException(
                        "Couldn't create term '" + ts +
                        "' for '" + s + "' in legacy ontology namespace"
                ).initCause(ex);
            }
        }
    }

    String getOntologyTerm(int termId) {
        return ontologySQL.termForID(termId).getName();
    }

    private boolean hierarchyChecked = false;
    private boolean hierarchySupported = false;

    boolean isHierarchySupported() {
        if (!hierarchyChecked) {
            hierarchySupported = getDBHelper().containsTable(dataSource, "seqfeature_relationship");
            hierarchyChecked = true;
        }

        return hierarchySupported;
    }

    private boolean assemblyChecked = false;
    private boolean assemblySupported = false;

    boolean isAssemblySupported() {
        if (!assemblyChecked) {
            assemblySupported = getDBHelper().containsTable(dataSource, "assembly");
            assemblyChecked = true;
        }

        return assemblySupported;
    }



    private boolean bioentryPropertyChecked = false;
    private boolean bioentryPropertySupported = false;

    boolean isBioentryPropertySupported() {
        if (!bioentryPropertyChecked) {
            bioentryPropertySupported = getDBHelper().containsTable(dataSource, "bioentry_qualifier_value");
            bioentryPropertyChecked = true;
        }

        return bioentryPropertySupported;
    }


    private boolean dbSchemaChecked = false;
    private boolean dbSchemaSupported = false;

    boolean isDbSchemaSupported() {
        if (!dbSchemaChecked) {
            dbSchemaSupported = getDBHelper().containsTable(dataSource, "location");
            dbSchemaChecked = true;
        }

        return dbSchemaSupported;
    }

    private boolean commentTableNameChecked = false;
    private String commentTableName = null;

    // Get the name of the table used for comments.
    // "comment" isn't allowed in oracle (where comment is a reserved word)
    // Hilmar has said this table will be renamed to anncomment post-1.0 BioSQLe
    // We support both "comment" and "anncomment"
    String getCommentTableName() {
        if (!commentTableNameChecked) {
            if (getDBHelper().containsTable(dataSource, "comment")) {
                commentTableName = "comment";
            } else if (getDBHelper().containsTable(dataSource, "anncomment")) {
                commentTableName = "anncomment";
            }
            commentTableNameChecked = true;
        }
        return commentTableName;
    }



    private boolean spaChecked = false;
    private boolean spaSupported = false;

    boolean isSPASupported() {
        if (!spaChecked) {
            Connection conn = null;
            try {
                spaSupported = false;
                conn = dataSource.getConnection();
                PreparedStatement ps = null;
                try {
                    ps = conn.prepareStatement("select biosql_accelerators_level()");
                    ResultSet rs = ps.executeQuery();
                    if (rs.next()) {
                        int level = rs.getInt(1);
                        if (level >= 2) {
                            spaSupported = true;
                            // System.err.println("*** Accelerators present in the database: level " + level);
                        }
                    }
                    rs.close();
                } catch (SQLException ex) {
                }
                if (ps != null) {
                    ps.close();
                }
                conn.close();

                spaChecked = true;
            } catch (SQLException ex) {
		    if (conn!=null) try {conn.close();} catch (SQLException ex3) {}
                throw new BioRuntimeException(ex);
            }
        }

        return spaSupported;
    }

    private class SqlizedFilter {
        private List tables = new ArrayList();
        private String filter;
        private int used_ot = 0;
        private int used_sfs = 0;
        //private int used_sqv = 0;

        SqlizedFilter(FeatureFilter ff) {
            filter = sqlizeFilter(ff, false);
        }

        private String sqlizeFilter(FeatureFilter ff, boolean negate) {
            if (ff instanceof FeatureFilter.ByType) {
                String type = ((FeatureFilter.ByType) ff).getType();
                String tableName = "ot_" + (used_ot++);
                tables.add("term as " + tableName);
                return tableName + ".name " + eq(negate) + qw(type) + " and seqfeature.type_term_id = " + tableName + ".term_id";
            } else if (ff instanceof FeatureFilter.BySource) {
                String source = ((FeatureFilter.BySource) ff).getSource();
                String tableName = "sfs_" + (used_sfs++);
                tables.add("term as " + tableName);
                return tableName + ".name " + eq(negate) + qw(source) + " and seqfeature.source_term_id = " + tableName + ".term_id";
            } else if (ff instanceof FeatureFilter.ByAnnotation) {

            return "";

            /* FIXME disabled until Matthew works out what he's doing with AnnotationTypes

                FeatureFilter.ByAnnotation ffba = (FeatureFilter.ByAnnotation) ff;
                Object key = ffba.getKey();
                Object value = ffba.getValue();
                String keyString = key.toString();
                String valueString = value.toString();

                String otName = "ot_" + (used_ot++);
                tables.add("ontology_term as " + otName);

                String sqvName = "sqv_" + (used_sqv++);
                tables.add("seqfeature_qualifier_value as " + sqvName);

        */

                // FIXME this doesn't actually do negate quite right -- it doesn't
                // match if the property isn't defined.  Should do an outer join
                // to fix this.  But for now, we'll just punt :-(

        /*

                if (negate) {
                    return "";
                }

                return sqvName + ".qualifier_value" + eq(negate) + qw(valueString) + " and " +
                       sqvName + ".term_id = " + otName + ".term_id and " +
                       otName + ".term_name = " + qw(keyString) + " and " +
                       "seqfeature.seqfeature_id = " + sqvName + ".seqfeature_id";

               */

            } else if (ff instanceof FeatureFilter.And) {
                FeatureFilter.And and = (FeatureFilter.And) ff;
                FeatureFilter ff1 = and.getChild1();
                FeatureFilter ff2 = and.getChild2();
                String filter1 = sqlizeFilter(ff1, negate);
                String filter2 = sqlizeFilter(ff2, negate);
                if (filter1.length() > 0) {
                    if (filter2.length() > 0) {
                        return filter1 + " and " + filter2;
                    } else {
                        return filter1;
                    }
                } else {
                    if (filter2.length() > 0) {
                        return filter2;
                    } else {
                        return "";
                    }
                }
            } else  if (ff instanceof FeatureFilter.Not) {
                FeatureFilter child = ((FeatureFilter.Not) ff).getChild();
                return sqlizeFilter(child, !negate);
            } else {
                return "";
            }
        }

        private String eq(boolean negate) {
            if (negate) {
                return " <> ";
            } else {
                return "=";
            }
        }

        private String qw(String word) {
            return "'" + word + "'";
        }

        public String getQuery() {
            StringBuffer query = new StringBuffer();
            query.append("select bioentry.accession, seqfeature.seqfeature_id ");
            query.append("  from seqfeature, bioentry");
            for (Iterator i = tables.iterator(); i.hasNext(); ) {
                query.append(", ");
                query.append((String) i.next());
            }
            query.append(" where bioentry.bioentry_id = seqfeature.bioentry_id");
            query.append("   and bioentry.biodatabase_id = ?");
            if (filter.length() > 0) {
                query.append(" and ");
                query.append(filter);
            }
            query.append(" order by bioentry.accession");

            return query.substring(0);
        }
    }

    private class FilterByInternalID implements FeatureFilter {
        private int id;

        public FilterByInternalID(int id) {
            this.id = id;
        }

        public boolean accept(Feature f) {
            if (! (f instanceof BioSQLFeature)) {
                return false;
            }

            int intID = ((BioSQLFeature) f)._getInternalID();
            return (intID == id);
        }
    }

    public FeatureHolder filter(FeatureFilter ff) {
        Connection conn = null;
        try {
            SqlizedFilter sqf = new SqlizedFilter(ff);
            System.err.println("Doing BioSQL filter");
            System.err.println(sqf.getQuery());

            conn = dataSource.getConnection();
            PreparedStatement get_features = conn.prepareStatement(sqf.getQuery());
            get_features.setInt(1, dbid);
            ResultSet rs = get_features.executeQuery();

            String lastAcc = "";
            Sequence seq = null;
            SimpleFeatureHolder fh = new SimpleFeatureHolder();

            while (rs.next()) {
                String accession = rs.getString(1);
                int fid = rs.getInt(2);

                System.err.println(accession + "\t" + fid);

                if (seq == null || ! lastAcc.equals(accession)) {
                    seq = getSequence(accession);
                }

                FeatureHolder hereFeature = seq.filter(new FilterByInternalID(fid), true);
                Feature f = (Feature) hereFeature.features().next();
                if (ff.accept(f)) {
                    fh.addFeature(f);
                }
            }
            rs.close();
            get_features.close();
            conn.close();

            return fh;
        } catch (SQLException ex) {
            if (conn!=null) try {conn.close();} catch (SQLException ex3) {}
            throw new BioRuntimeException("Error accessing BioSQL tables", ex);
        } catch (ChangeVetoException ex) {
            if (conn!=null) try {conn.close();} catch (SQLException ex3) {}
            throw new BioError("Assert failed: couldn't modify internal FeatureHolder", ex);
        } catch (BioException ex) {
            if (conn!=null) try {conn.close();} catch (SQLException ex3) {}
            throw new BioRuntimeException("Error fetching sequence", ex);
        }
    }

    // SequenceIterator here, 'cos AbstractSequenceDB grandfathers in AbstractChangable :-(

  public SequenceIterator sequenceIterator() {
    return new SequenceIterator() {
      private Iterator pID = ids().iterator();

      public boolean hasNext() {
        return pID.hasNext();
      }

      public Sequence nextSequence() throws BioException {
        return getSequence((String) pID.next());
      }
    };
  }
    // change support stuff
    void firePreChangeEvent(ChangeEvent cev)
        throws ChangeVetoException
    {
        getChangeSupport(cev.getType()).firePreChangeEvent(cev);
    }

    void firePostChangeEvent(ChangeEvent cev)
    {
        getChangeSupport(cev.getType()).firePostChangeEvent(cev);
    }

    //
    // Feature canonicalization
    //

    BioSQLFeature canonicalizeFeature(BioSQLFeature f, int feature_id) {
        // System.err.println("Canonicalizing feature at " + f.getLocation());

        Integer key = new Integer(feature_id);
        BioSQLFeature oldFeature = (BioSQLFeature) featuresByID.get(key);
        if (oldFeature != null) {
            return oldFeature;
        } else {
            featuresByID.put(key, f);
            return f;
        }
    }

    private class SingleFeatureReceiver extends BioSQLFeatureReceiver {
        private Feature feature;

        private SingleFeatureReceiver() {
            super(BioSQLSequenceDB.this);
        }

        protected void deliverTopLevelFeature(Feature f)
            throws ParseException
        {
            if (feature == null) {
                feature = f;
            } else {
                throw new ParseException("Expecting only a single feature");
            }
        }

        public Feature getFeature() {
            return feature;
        }
    }

    BioSQLFeature getFeatureByID(int feature_id)
    {
        Integer key = new Integer(feature_id);
        BioSQLFeature f = (BioSQLFeature) featuresByID.get(key);
        if (f != null) {
            return f;
        }

        try {
            SingleFeatureReceiver receiver = new SingleFeatureReceiver();
            getFeaturesSQL().retrieveFeatures(-1, receiver, null, -1, feature_id);
            if (receiver.getFeature() == null) {
                throw new BioRuntimeException("Dangling internal_feature_id");
            } else {
                featuresByID.put(key, (BioSQLFeature) receiver.getFeature());
                return (BioSQLFeature) receiver.getFeature();
            }
        } catch (SQLException ex) {
            throw new BioRuntimeException("Database error", ex);
        } catch (BioException ex) {
            throw new BioRuntimeException(ex);
        }
    }
/*
    //
    // Dbxref canonicalization
    //
    BioSQLXRef canonicalizeXRef(BioSQLXRef r, int dbxref_id) {

        Integer key = new Integer(dbxref_id);
        BioSQLXRef oldXRef = (BioSQLXRef) XRefsByID.get(key);
        if (oldXRef != null) {
            return oldXRef;
        } else {
            featuresByID.put(key, r);
            return r;
        }
    }

    BioSQLXRef getXRefsByID(int dbxref_id)
    {
        Integer key = new Integer(dbxref_id);
        BioSQLFeature f = (BioSQLFeature) featuresByID.get(key);
        if (f != null) {
            return f;
        }

        try {
            SingleFeatureReceiver receiver = new SingleFeatureReceiver();
            getFeaturesSQL().retrieveFeatures(-1, receiver, null, -1, feature_id);
            if (receiver.getFeature() == null) {
                throw new BioRuntimeException("Dangling internal_feature_id");
            } else {
                featuresByID.put(key, (BioSQLFeature) receiver.getFeature());
                return (BioSQLFeature) receiver.getFeature();
            }
        } catch (SQLException ex) {
            throw new BioRuntimeException(ex, "Database error");
        } catch (BioException ex) {
            throw new BioRuntimeException(ex);
        }
    }
*/

    Cache getTileCache() {
        return tileCache;
    }
}
