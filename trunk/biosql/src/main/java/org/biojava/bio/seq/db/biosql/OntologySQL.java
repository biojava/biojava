/* -*- c-basic-offset: 2; indent-tabs-mode: nil -*- */
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
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.NoSuchElementException;

import javax.sql.DataSource;

import org.biojava.bio.BioError;
import org.biojava.bio.BioException;
import org.biojava.bio.BioRuntimeException;
import org.biojava.ontology.AlreadyExistsException;
import org.biojava.ontology.OntoTools;
import org.biojava.ontology.Ontology;
import org.biojava.ontology.OntologyException;
import org.biojava.ontology.OntologyTerm;
import org.biojava.ontology.RemoteTerm;
import org.biojava.ontology.Term;
import org.biojava.ontology.Triple;
import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeListener;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;

/**
 * Behind-the-scenes adaptor to the features sub-schema of BioSQL.
 *
 * @author Thomas Down
 * @author Matthew Pocock
 * @author Eric Haugen
 * @author Len Trigg
 * @author Richard Holland
 * @deprecated Use hibernate and org.biojavax.bio.db.*
 * @since 1.4
 */
class OntologySQL {

  private static HashMap ONTOLOGIES = new HashMap();

  static synchronized OntologySQL getOntologySQL(DataSource source,
                                                 DBHelper dbHelper)
          throws BioException, SQLException
  {
    OntologySQL ontologySQL = (OntologySQL) ONTOLOGIES.get(source);
    if (ontologySQL == null) {
      ontologySQL = new OntologySQL(source, dbHelper);
      ONTOLOGIES.put(source, ontologySQL);
    }
    return ontologySQL;
  }

  static synchronized void clearCache() {
    ONTOLOGIES.clear();
  }

  private DataSource source;
  private DBHelper dbHelper;
  private Map ontologiesByID;
  private Map ontologiesByName;
  private Map termsByID;
  private Map IDsByTerm;
  private Map monitors;

  private Map blessedExternalAliases;
  private Map blessedExternalTerms;

  private Ontology guano;

  {
    ontologiesByID = new HashMap();
    ontologiesByName = new HashMap();
    termsByID = new HashMap();
    IDsByTerm = new HashMap();
    monitors = new HashMap();

    blessedExternalAliases = new HashMap();
    blessedExternalTerms = new HashMap();
  }

  Ontology getLegacyOntology() {
    return guano;
  }

  public Ontology getOntology(String name) throws BioException {
    Object ont = ontologiesByName.get(name);
    if (ont == null) {
      throw new NoSuchElementException("Can't find ontology named " + name);
    } else if (ont instanceof OntologyPlaceholder) {
      OntologyPlaceholder op = (OntologyPlaceholder) ont;
      Ontology onto = loadOntology(op.name, op.description, op.id);
      OntologyMonitor om = new OntologyMonitor(onto);
      monitors.put(onto, om);
      onto.addChangeListener(om, ChangeType.UNKNOWN);
      return onto;
    } else {
      return (Ontology) ont;
    }
  }

  public Ontology createOntology(String name, String description)
          throws AlreadyExistsException, BioException
  {
    if (ontologiesByName.containsKey(name)) {
      throw new AlreadyExistsException("This BioSQL database already contains an ontology of name " + name);
    }
    Ontology ont = new Ontology.Impl(name, description);

    persistOntology(ont);

    OntologyMonitor om = new OntologyMonitor(ont);
    monitors.put(ont, om);
    ont.addChangeListener(om, ChangeType.UNKNOWN);

    return ont;
  }

  public Ontology addOntology(Ontology old)
          throws AlreadyExistsException
  {
    if (ontologiesByName.containsKey(old.getName())) {
      throw new AlreadyExistsException("This BioSQL database already contains an ontology of name " + old.getName());
    }

    Connection conn = null;
    try {
      conn = source.getConnection();
      conn.setAutoCommit(false);
      Ontology ont = new Ontology.Impl(old.getName(), old.getDescription());
      persistOntology(conn, ont);

      Map localTerms = new HashMap();
      for (Iterator i = old.getTerms().iterator(); i.hasNext(); ) {
        Term t = (Term) i.next();
        Term localTerm;
        if (t instanceof RemoteTerm) {
          localTerm = ont.importTerm(((RemoteTerm) t).getRemoteTerm(), null);
        } else {
          localTerm = ont.createTerm(t.getName(), t.getDescription(), t.getSynonyms());
          persistTerm(conn, localTerm);
        }
        localTerms.put(t, localTerm);
      }
      for (Iterator i = old.getTriples(null, null, null).iterator(); i.hasNext(); ) {
        Triple t = (Triple) i.next();
        Triple localT = ont.createTriple(
                (Term) localTerms.get(t.getSubject()),
                (Term) localTerms.get(t.getObject()),
                (Term) localTerms.get(t.getPredicate()),
                null, null
        );

        persistTriple(conn, ont, localT);
      }

      conn.commit();

      OntologyMonitor om = new OntologyMonitor(ont);
      monitors.put(ont, om);
      ont.addChangeListener(om, ChangeType.UNKNOWN);

      return ont;
    } catch (SQLException ex) {
      boolean rolledback = false;
      if (conn != null) {
        try {
          conn.rollback();
          rolledback = true;
        } catch (SQLException ex2) {}
        try {conn.close();} catch (SQLException ex3) {}
      }
      throw new BioRuntimeException("Error removing from BioSQL tables" + (rolledback ? " (rolled back successfully)" : ""),ex);
    } catch (AlreadyExistsException ex) {
      if (conn!=null) try {conn.close();} catch (SQLException ex3) {}
      throw new BioError("Unexpected ontology duplication error",ex);
    } catch (ChangeVetoException ex) {
      if (conn!=null) try {conn.close();} catch (SQLException ex3) {}
      throw new BioError("Unexpected veto altering internal Ontology object",ex);
    }
  }

  public void addCore(Connection conn)
          throws SQLException
  {
    System.err.println("*** Importing a core ontology -- hope this is okay");

    try {
      conn.setAutoCommit(false);
      Ontology old = OntoTools.getCoreOntology();
      Ontology ont = new Ontology.Impl("__core_ontology", "BioSQL core ontology (imported by BioJava)");
      persistOntology(conn, ont);

      Map localTerms = new HashMap();
      System.err.println("*** Importing terms");
      for (Iterator i = old.getTerms().iterator(); i.hasNext(); ) {
        Term t = (Term) i.next();
        processTerm(t, localTerms, ont, conn);
      }

      conn.commit();

      OntologyMonitor om = new OntologyMonitor(ont);
      monitors.put(ont, om);
      ont.addChangeListener(om, ChangeType.UNKNOWN);
      blessExternal(ont, old);
    } catch (AlreadyExistsException ex) {
      throw new BioError("Unexpected ontology duplication error",ex);
    } catch (ChangeVetoException ex) {
      throw new BioError("Unexpected veto altering internal Ontology object",ex);
    }
  }

  private void processTerm(Term t,
                           Map localTerms,
                           Ontology ont,
                           Connection conn)
          throws AlreadyExistsException, ChangeVetoException, SQLException
  {
    //System.err.println("Processing term: " + t);
    Term localTerm;

    if(localTerms.keySet().contains(t)) {
      // don't double-import
      return;
    }

    if (t instanceof Triple) {
      Triple trip = (Triple) t;

      processTerm(trip.getSubject(),   localTerms, ont, conn);
      processTerm(trip.getObject(),    localTerms, ont, conn);
      processTerm(trip.getPredicate(), localTerms, ont, conn);

      localTerm = ont.createTriple(
              (Term) localTerms.get(trip.getSubject()),
              (Term) localTerms.get(trip.getObject()),
              (Term) localTerms.get(trip.getPredicate()),
              null, null
      );

      persistTriple(conn, ont, (Triple) localTerm);
    } else if (t instanceof RemoteTerm) {
      localTerm = ont.importTerm(((RemoteTerm) t).getRemoteTerm(), null);
    } else {
      localTerm = ont.createTerm(t.getName(), t.getDescription(), t.getSynonyms());
      persistTerm(conn, localTerm);
    }

    localTerms.put(t, localTerm);

    //System.err.println("Done processing term: " + t);
  }

  private void loadTerms(Ontology ont, int id)
          throws SQLException, OntologyException, ChangeVetoException
  {
    Connection conn = source.getConnection();
    String query =
        " SELECT term.term_id, term.name, term.definition, " +
        "       term_relationship.term_relationship_id, " +
        "       term_relationship.subject_term_id, " +
        "       term_relationship.object_term_id, " +
        "       term_relationship.predicate_term_id ";
    if (dbHelper.getJoinStyle() == DBHelper.JOIN_ORACLE8) {
      query +=
        "FROM term, term_relationship_term, term_relationship " +
        "     WHERE term.term_id = term_relationship_term.term_id (+) " +
        "     AND term_relationship_term.term_relationship_id = term_relationship.term_relationship_id (+) " +
        "AND term.ontology_id = ? ";
    } else {
      query +=
        "FROM term LEFT OUTER JOIN term_relationship_term " +
        "     ON term.term_id = term_relationship_term.term_id LEFT OUTER JOIN term_relationship " +
        "     ON term_relationship_term.term_relationship_id = term_relationship.term_relationship_id " +
        "WHERE term.ontology_id = ? ";
    }
    query += "ORDER BY term.term_id";
    PreparedStatement get_terms = conn.prepareStatement(query);
    
    String query2 = "SELECT synonym FROM term_synonym WHERE term_id = ?";
    PreparedStatement get_synonyms = conn.prepareStatement(query2);

    get_terms.setInt(1, id);
    ResultSet rs = get_terms.executeQuery();
    while (rs.next()) {
      int term_id = rs.getInt(1);
      Integer tid = new Integer(term_id);

      String name = rs.getString(2);
      String description = rs.getString(3);
      if (description == null) {
        description = "";
      }
      int term_relation_id = rs.getInt(4);

      Term t = null;

      if(term_relation_id == 0) {
        t = ont.createTerm(name, description);
      } else {
        int subject_term_id = rs.getInt(5);
        int object_term_id = rs.getInt(6);
        int predicate_term_id = rs.getInt(7);

        //System.err.println("Loading triple: " + term_relation_id + " -> " + subject_term_id + ", " + object_term_id + ", " + predicate_term_id);

        t = ont.createTriple((Term) termsByID.get(new Integer(subject_term_id)),
                             (Term) termsByID.get(new Integer(object_term_id)),
                             (Term) termsByID.get(new Integer(predicate_term_id)),
                             name, description);
      }
                
      get_synonyms.setInt(1, term_id);
      ResultSet rs2 = get_synonyms.executeQuery();
      while (rs2.next()) t.addSynonym(rs2.getString(1));
      rs2.close();
      
      //System.err.println("Loaded term: " + tid + " -> " + t);

      termsByID.put(tid, t);
      IDsByTerm.put(t, tid);
    }
    rs.close();
    get_terms.close();
    get_synonyms.close();

    conn.close();
    conn = null;

    if (ont.getName().equals("__core_ontology")) {
      blessExternal(ont, OntoTools.getCoreOntology());
    }
  }

  public void blessExternal(Ontology internalOntology, Ontology externalOntology) {
    for (Iterator i = externalOntology.getTerms().iterator(); i.hasNext(); ) {
      Term extTerm = (Term) i.next();
      Term intTerm = internalOntology.getTerm(extTerm.getName());
      blessedExternalAliases.put(extTerm, intTerm);
      blessedExternalTerms.put(intTerm, extTerm);
    }
  }

  OntologySQL(DataSource source, DBHelper dbHelper) throws SQLException, BioException {
    this.source = source;
    this.dbHelper = dbHelper;

    Connection conn = source.getConnection();

    try {
      PreparedStatement get_onts = conn.prepareStatement(
              "select ontology_id, name, definition " +
              "  from ontology "
      );
      ResultSet rs = get_onts.executeQuery();
      while (rs.next()) {
        int id = rs.getInt(1);
        String name = rs.getString(2);
        String description = rs.getString(3);

        // Postpone loading of ontologies until they are first accessed by getOntology()
        OntologyPlaceholder op = new OntologyPlaceholder(name, description, id);
        ontologiesByID.put(new Integer(id), op);
        ontologiesByName.put(name, op);

      }
      rs.close();
      get_onts.close();

      if (!ontologiesByName.containsKey("__core_ontology")) {
        addCore(conn);
      }

      if (ontologiesByName.containsKey("__biojava_guano")) {
        guano = getOntology("__biojava_guano");
      } else {
        guano = createOntology("__biojava_guano", "Namespace for old, but still useful, shit imported from ontology-less BioJava data models");
      }
    } catch (AlreadyExistsException ex) {
      try {conn.close();} catch (SQLException ex3) {}
      throw new BioException("Duplicate term name in BioSQL",ex);
    }

    conn.close();
    conn = null;
  }

  private static class OntologyPlaceholder {
    final String name;
    final String description;
    final int id;
    public OntologyPlaceholder(String name, String description, int id) {
      this.name = name;
      this.description = description;
      this.id = id;
    }
  }

  private Ontology loadOntology(String name, String description, int id) throws BioException {
    Ontology ont = new Ontology.Impl(name, description);
    try {
      //long start = System.currentTimeMillis();
      loadTerms(ont, id);
      //System.err.println("Loading terms for " + name + " took " + (System.currentTimeMillis() - start) + "ms");
    } catch (SQLException ex) {
      throw new BioException("Error loading ontology terms", ex);
    } catch (OntologyException ex) {
      throw new BioException("Error loading ontology terms", ex);
    } catch (ChangeVetoException ex) {
      throw new BioError("Assertion failed: couldn't modify Ontology", ex);
    }
    ontologiesByID.put(new Integer(id), ont);
    ontologiesByName.put(name, ont);
    return ont;
  }


  private class OntologyMonitor implements ChangeListener {
    private Ontology ontology;

    OntologyMonitor(Ontology o) {
      this.ontology = o;
    }

    public void preChange(ChangeEvent cev)
            throws ChangeVetoException
    {
      ChangeType type = cev.getType();
      if (type.isMatchingType(Ontology.TERM)) {
        ChangeEvent chain = cev.getChainedEvent();
        if (chain != null) {
          // looks like it's an actual term which has been added.
          throw new ChangeVetoException("Biojava does not handle mutable terms once persisted");
        } else {
          if (cev.getChange() != null && cev.getPrevious() == null) {
            // Adding
            if (! (cev.getChange() instanceof Term)) {
              throw new ChangeVetoException("Can't understand this change");
            }
            Term addedTerm = (Term) cev.getChange();
            if (addedTerm instanceof OntologyTerm) {
              throw new ChangeVetoException("BioSQL does (currently) represent OntologyTerms but Biojava wants to lie about that");
            } else if (addedTerm instanceof RemoteTerm) {
              Term gopher = addedTerm;
              while (gopher instanceof RemoteTerm) {
                gopher = ((RemoteTerm) gopher).getRemoteTerm();
                if (!ontologiesByID.values().contains(gopher.getOntology()) && !blessedExternalAliases.containsKey(gopher)) {
                  throw new ChangeVetoException("BioSQL ontologies can't contain references to external ontologies");
                }
              }
            } else {
              // We presume it's a term with no special semantics
            }
          } else if (cev.getChange() == null && cev.getPrevious() != null) {
            throw new ChangeVetoException("FIXME: Biojava can't remove terms from biosql ontology");
          } else {
            throw new ChangeVetoException("Unknown TERM change");
          }
        }
      } else if (type.isMatchingType(Ontology.TRIPLE)) {
        // Any triple-changes that we can't handle?  I don't think so.
      } else {
        throw new ChangeVetoException("Biojava does not understand this change");
      }
    }

    public void postChange(ChangeEvent cev) {
      ChangeType type = cev.getType();
      if (type.isMatchingType(Ontology.TERM)) {
        if (cev.getChange() != null && cev.getPrevious() == null) {
          // Adding
          if (! (cev.getChange() instanceof Term)) {
            throw new BioError("Assertion failed: added object isn't a term: " + cev.getChange());
          }
          Term addedTerm = (Term) cev.getChange();
          if (addedTerm instanceof RemoteTerm) {
            // We're not actually persisting these explicitly.  Hope this isn't a problem
          }
          persistTerm(addedTerm);
        }
        // Don't support removal yet, but it should have been vetoed.
      } else if (type.isMatchingType(Ontology.TRIPLE)) {
        if (cev.getChange() != null && cev.getPrevious() == null) {
          // Adding
          if (! (cev.getChange() instanceof Triple)) {
            throw new BioError("Assertion failed: added object isn't a triple: " + cev.getChange());
          }
          persistTriple(ontology, (Triple) cev.getChange());
        }
      }
    }
  }

  private void persistTerm(Term term) {
    Connection conn = null;
    try {
      conn = source.getConnection();
      conn.setAutoCommit(false);

      persistTerm(conn, term);

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
      throw new BioRuntimeException("Error commiting to BioSQL tables" + (rolledback ? " (rolled back successfully)" : ""),ex);
    }
  }

  private void persistTerm(Connection conn, Term term)
          throws SQLException
  {
    try {
      //System.err.println("Persisting " + term + " using " + term.getName() + ", " + term.getDescription() + ", " + ontologyID(term.getOntology()));

      PreparedStatement import_term = conn.prepareStatement(
              "insert into term " +
              "       (name, definition, ontology_id) " +
              "values (?, ?, ?)"
      );
      import_term.setString(1, term.getName());
      import_term.setString(2, term.getDescription());
      import_term.setInt(3, ontologyID(term.getOntology()));
      import_term.executeUpdate();
      import_term.close();
      int id = dbHelper.getInsertID(conn, "term", "term_id");
      Integer tid = new Integer(id);
      termsByID.put(tid, term);
      IDsByTerm.put(term, tid);
      
      PreparedStatement stored_synonym = conn.prepareStatement(
              "insert into term_synonym (term_id, name) values (?,?) "
              );
      Object[] synonyms = term.getSynonyms();
      for (int i = 0; i < synonyms.length; i++) {
          stored_synonym.setInt(1,id);
          stored_synonym.setString(2,""+synonyms[i]);
          stored_synonym.executeUpdate();
      }
      stored_synonym.close();
      
      //System.err.println("Persisted term: " + tid + " -> " + term);
    } catch (SQLException se) {
      throw (SQLException) new SQLException(
              "Failed to persist term: " + term +
              " from ontology: " + term.getOntology() +
              " with error: " + se.getErrorCode() + " : " + se.getSQLState()
      ).initCause(se);
    }
  }

  private void persistTriple(Ontology ont, Triple triple) {
    Connection conn = null;
    try {
      conn = source.getConnection();
      conn.setAutoCommit(false);

      persistTriple(conn, ont, triple);

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
      throw new BioRuntimeException("Error adding to BioSQL tables" + (rolledback ? " (rolled back successfully)" : ""),ex);
    }
  }

  private void persistTriple(Connection conn, Ontology ont, Triple triple)
    throws SQLException {
   try {
      // Code borrowed from termID to check for triple's existence.
      Term t = triple;
      while (t instanceof RemoteTerm) t = ((RemoteTerm)t).getRemoteTerm();
      if(blessedExternalAliases.containsKey(t)) t = (Term)
         blessedExternalAliases.get(t);
      if (!IDsByTerm.containsKey(t)) persistTerm(conn,triple);
                                 
      // End of borrowed code.
      PreparedStatement import_trip = conn.prepareStatement(
            "insert into term_relationship " +
            " (subject_term_id, predicate_term_id,object_term_id, ontology_id) " +
            "values (?, ?, ?, ?)");

      import_trip.setInt(1, termID(triple.getSubject()));
      import_trip.setInt(2, termID(triple.getPredicate()));
      import_trip.setInt(3, termID(triple.getObject()));
      import_trip.setInt(4, ontologyID(ont));
      import_trip.executeUpdate();
      import_trip.close();
      int tripID = dbHelper.getInsertID(conn, "term_relationship", "term_relationship_id");
      PreparedStatement link_trip_to_term = conn.prepareStatement(
            "insert into term_relationship_term " +
            " (term_relationship_id, term_id) " +
            "values (?, ?)" );
      link_trip_to_term.setInt(1, tripID);
      link_trip_to_term.setInt(2, termID(triple));
      link_trip_to_term.executeUpdate();
      link_trip_to_term.close();
      //System.err.println("Persisted triple: " + triple);
    } catch (SQLException se) {
        throw (SQLException) new SQLException(
          "Failed to persist triple: " + triple+ " from ontology: " + 
           triple.getOntology() + " with error: " +se.getErrorCode() + 
	   " : " +se.getSQLState()).initCause(se);
    }
  }


  private void persistOntology(Ontology onto)
  {
    Connection conn = null;
    try {
      conn = source.getConnection();
      conn.setAutoCommit(false);

      persistOntology(conn, onto);

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
      throw new BioRuntimeException("Error adding to BioSQL tables" + (rolledback ? " (rolled back successfully)" : ""),ex);
    }
  }

  private void persistOntology(Connection conn, Ontology ont)
          throws SQLException
  {
    PreparedStatement insert_ontology = conn.prepareStatement(
            "insert into ontology " +
            "       (name, definition) " +
            "values (?, ?)"
    );
    insert_ontology.setString(1, ont.getName());
    insert_ontology.setString(2, ont.getDescription());
    insert_ontology.executeUpdate();
    insert_ontology.close();
    int id = dbHelper.getInsertID(conn, "ontology", "ontology_id");
    Integer tid = new Integer(id);
    ontologiesByID.put(tid, ont);
    ontologiesByName.put(ont.getName(), ont);

  }

  int ontologyID(Ontology ont) {
    for (Iterator i = ontologiesByID.entrySet().iterator(); i.hasNext(); ) {
      Map.Entry me = (Map.Entry) i.next();
      if (me.getValue().equals(ont)) {
        return ((Integer) me.getKey()).intValue();
      }
    }
    throw new BioError("Couldn't find ontology " + ont.getName());
  }

  int termID(Term t) {
    while (t instanceof RemoteTerm) {
      t = ((RemoteTerm) t).getRemoteTerm();
    }
    if (blessedExternalAliases.containsKey(t)) {
      t = (Term) blessedExternalAliases.get(t);
    }
    try {
      return ((Integer) IDsByTerm.get(t)).intValue();
    } catch (NullPointerException ex) {
      throw new BioError("Error looking up biosqlized ID for " + t.toString(),ex);
    }
  }

  Term termForID(int id) {
    Term t = (Term) termsByID.get(new Integer(id));
    if (t == null) {
      throw new BioError("Invalid term id " + id);
    }
    return t;
  }
}
