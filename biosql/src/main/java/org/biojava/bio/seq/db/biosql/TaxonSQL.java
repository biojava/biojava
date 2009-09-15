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

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.sql.Types;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.Set;
import java.util.Stack;

import org.biojava.bio.Annotation;
import org.biojava.bio.BioException;
import org.biojava.bio.BioRuntimeException;
import org.biojava.bio.taxa.CircularReferenceException;
import org.biojava.bio.taxa.EbiFormat;
import org.biojava.bio.taxa.Taxon;
import org.biojava.bio.taxa.TaxonFactory;
import org.biojava.bio.taxa.WeakTaxonFactory;
import org.biojava.utils.ChangeVetoException;

/**
 * Methods for retrieving, storing and manipulate Taxa stored in a BioSQL database.
 *
 * @author Len Trigg
 * @author Andreas Dr&auml;ger
 * @deprecated Use hibernate and org.biojavax.bio.db.*
 */
public class TaxonSQL {

    /** 
     * Attempts to get a Taxon object corresponding to the specified
     * NCBI taxon ID.
     *
     * @param conn the connection to the database
     * @param ncbi_taxon_id the NCBI taxon ID.
     * @return the corresponding Taxon (which may have already been
     * present in memory after an earlier retrieval), or null if the
     * Taxon could not be found in the database.
     */
    public static Taxon getTaxon(Connection conn, int ncbi_taxon_id) 
    {
      PreparedStatement ps = null;
      ResultSet rs = null;
      Taxon t = null;
      try {
        ps = conn.prepareStatement(
          "SELECT name " +
          "FROM   taxon_name, taxon " +
          "WHERE  taxon_name.taxon_id = taxon.taxon_id AND" +
          "       taxon_name.name_class LIKE 'scientific%' AND" +
          "       ncbi_taxon_id = ?");
        ps.setInt(1, ncbi_taxon_id);
        rs = ps.executeQuery();
        if (rs.next()) t = getTaxon(conn, rs.getString(1));
      } catch (SQLException exc) {
        throw new BioRuntimeException(exc);
      } finally {
        attemptClose(rs);
        attemptClose(ps);
      }
      return t; 
    }


    /** 
     * Attempts to get a Taxon object corresponding to the specified
     * taxon_id (i.e. the database's internal id for the taxon).
     *
     * @param conn the connection to the database
     * @param taxon_id the database-specific id for the Taxon.
     * @return the corresponding Taxon (which may have already been
     * present in memory after an earlier retrieval).
     */
    public static Taxon getDBTaxon(Connection conn, int taxon_id) throws SQLException, ChangeVetoException 
    {
      PreparedStatement statement = null;
      ResultSet rs = null;
      try {
        // Constants for our wee id array
        final int NCBI_ID = 1;
        final int TAXON_ID = 0;
          
        // First, get the taxon ids up to the root.
        statement = conn.prepareStatement(
          "SELECT ncbi_taxon_id, parent_taxon_id "+ 
          "FROM   taxon "+ 
          "WHERE  taxon_id = ? ");
          
        ArrayList path = new ArrayList();
        while (taxon_id != 0) {
          statement.setInt(1, taxon_id);
          rs = statement.executeQuery();
          if (rs.next()) {
            path.add(new int [] {taxon_id, rs.getInt(1)});
            taxon_id = rs.getInt(2);
            if (rs.wasNull()) 
              taxon_id = 0;
          } else 
            throw new BioRuntimeException("Error fetching taxonomy structure. No taxon with taxon_id=" + taxon_id);
          rs.close();
        }
        statement.close();
          
        // Traverse from the root down as far has has been created previously...
        TaxonFactory factory = WeakTaxonFactory.GLOBAL;
        Taxon taxon = factory.getRoot();
        int pos = path.size() - 1;
        int []ids = (int[]) path.get(pos--);
        Map names = getTaxonNames(conn, ids[TAXON_ID]);
        taxon.getAnnotation().setProperty(EbiFormat.PROPERTY_NCBI_TAXON, "" + ids[NCBI_ID]);
        taxon.getAnnotation().setProperty(EbiFormat.PROPERTY_TAXON_NAMES, names);
        for (; pos >= 0; pos--) {
          // Who's the next id down the path?
          ids = (int[]) path.get(pos);
          String nextID = "" + ids[NCBI_ID];
          // Now look among the children for the next child.
          Set children = taxon.getChildren();
          for (Iterator it = children.iterator(); it.hasNext(); ) {
            Taxon child = (Taxon) it.next();
            Annotation anno = child.getAnnotation();
            if (anno.containsProperty(EbiFormat.PROPERTY_NCBI_TAXON)) {
              String childID = (String) anno.getProperty(EbiFormat.PROPERTY_NCBI_TAXON);
              if (childID.equals(nextID)) {
                taxon = child;
                continue;
              }
            } else {
              throw new BioRuntimeException("Taxon has not been annotated with NCBI taxon ids.");
            }
          }
          // No child with desired ncbi_id has been found.
          break;
        }
          
        // Now create taxa from here on down.
        try {
          for (; pos >= 0; pos--) 
          {
            // Now look for the next child.
            ids = (int[]) path.get(pos);
            String nextID = "" + ids[NCBI_ID];
            names = getTaxonNames(conn, ids[TAXON_ID]);
            String sciName = (String) names.get("scientific name");
            if (sciName == null) 
              throw new BioRuntimeException("No scientific name for taxon_id=" + ids[TAXON_ID]);
            String commonName = (String) names.get("common name");
            taxon = factory.addChild(taxon, factory.createTaxon(sciName, commonName));
            taxon.getAnnotation().setProperty(EbiFormat.PROPERTY_NCBI_TAXON, nextID);
            taxon.getAnnotation().setProperty(EbiFormat.PROPERTY_TAXON_NAMES, names);
            taxon = getProperties(conn, taxon);
          }           
        } catch (CircularReferenceException ex) {
          throw new BioRuntimeException("Circular references in taxon table. taxon_id=" + ids[TAXON_ID]);
        }
      
        return getProperties(conn, taxon);
      } finally {
        attemptClose(rs);
        attemptClose(statement);
      }
    }
    

    /**
     * Look up all the names associated with a taxon_id.
     *
     * @param conn the current <code>Connection</code>.
     * @param taxon_id the NCBI taxon id for the taxon of interest.
     * @return a <code>Map</code> from name_class (e.g.: "scientific
     * name") to name.
     */
    private static Map getTaxonNames(Connection conn, int taxon_id) 
    {
      PreparedStatement statement = null;
      ResultSet rs = null;
      try {
        statement = conn.prepareStatement(
          "SELECT name_class, name "+
          "FROM   taxon_name "+ 
          "WHERE  taxon_id = ? ");
        statement.setInt(1, taxon_id);
        rs = statement.executeQuery();

        Map names = new HashMap();
        while (rs.next()) {  	
          String name_class = rs.getString(1);
          String name = rs.getString(2);
          if ((name_class.equals("scientific name")) ||
              (name_class.equals("common name")))
              names.put(name_class, name);
          else {
            Set s = new HashSet();
            s.add(name);
            if (names.containsKey(name_class))
              s.addAll((Collection) names.get(name_class));
              names.put(name_class, s);
            }
        }

        return names;
      } catch (SQLException ex) {
        throw new BioRuntimeException("Error fetching taxonomy annotations", ex);
      } finally {
        attemptClose(rs);
        attemptClose(statement);
      }
    }


    /**
     * Adds a <code>Taxon</code> (along with its parents) to the
     * database. If it is already present in the database, no action
     * is taken.  Returns the id by which the database refers to the
     * specified <code>Taxon</code> object.
     *
     * @param taxon a <code>Taxon</code>. The <code>Taxon</code> must
     * be annotated with the NCBI taxon id
     * (<code>key=EbiFormat.PROPERTY_ORGANISM</code>).
     * @param helper for the certain database system which is in use.
     * @return an <code>int</code> that corresponds to the
     * <code>Taxon</code> in the database.
     */
    public static int putTaxon(Connection conn, DBHelper helper, Taxon taxon) throws SQLException 
    {
      // Find the NCBI taxon id annotation 
      Annotation anno = taxon.getAnnotation();
      Object t  = anno.getProperty(EbiFormat.PROPERTY_NCBI_TAXON);
      if (t instanceof List) {
        t = (String) ((List) t).get(0);
      }
      int ncbi_taxon_id = Integer.parseInt((String) t);
      PreparedStatement selectTaxon = conn.prepareStatement(
        "select taxon_id " 
        + "from taxon " 
        + "where ncbi_taxon_id = ? "
      );
      selectTaxon.setInt(1, ncbi_taxon_id);
      ResultSet trs = selectTaxon.executeQuery();
      int taxon_id;
      if (trs.next()) {
        // entry exists - link to it
        taxon_id = trs.getInt(1);
      } else {
        // Taxon entry does not exist - create it
        Taxon parent = taxon.getParent();
        
        // Mein modifizierter Code:
        String rank, gencode, mitocode, left, right;
        
        rank = gencode = mitocode = left = right = "NULL";
        
        if (taxon.getAnnotation().containsProperty("rank"))
          rank = taxon.getAnnotation().getProperty("rank").toString();
        if (taxon.getAnnotation().containsProperty("genetic code id"))
          gencode = taxon.getAnnotation().getProperty("genetic code id").toString();
        if (taxon.getAnnotation().containsProperty("mitochondrial genetic code id"))
          mitocode = taxon.getAnnotation().getProperty("mitochondrial genetic code id").toString();
        if (taxon.getAnnotation().containsProperty("left value"))
          left = taxon.getAnnotation().getProperty("left value").toString();
        if (taxon.getAnnotation().containsProperty("right value"))
          right = taxon.getAnnotation().getProperty("right value").toString();
            
            
        PreparedStatement createTaxon = null;
        if (parent != null) {
          int parent_taxon_id = putTaxon(conn, helper, parent);
          createTaxon = conn.prepareStatement(
             "insert into taxon " 
             + "(ncbi_taxon_id, parent_taxon_id,"
             + " node_rank, genetic_code,"
             + " mito_genetic_code, left_value,"
             + " right_value) " 
             + "values (?, ?, ?, ?, ?, ?, ?)"
           );
          //System.out.println(ncbi_taxon_id+"\t"+rank+"\t"+gencode+"\t"+"\t"+mitocode+"\t"+left+"\t"+right);
          createTaxon.setInt(1, ncbi_taxon_id);
          createTaxon.setInt(2, parent_taxon_id);
          if (rank.equals("NULL")) 
            createTaxon.setNull(3, Types.VARCHAR);
          else createTaxon.setString(3, rank);
          if (gencode.equals("NULL"))
            createTaxon.setNull(4, Types.TINYINT);
          else createTaxon.setInt(4, new Integer(gencode).intValue());
          if (mitocode.equals("NULL"))
            createTaxon.setNull(5, Types.TINYINT);
          else createTaxon.setInt(5, new Integer(mitocode).intValue());
          if (left.equals("NULL"))
            createTaxon.setNull(6, Types.INTEGER);
          else createTaxon.setInt(6, new Integer(left).intValue());
          if (right.equals("NULL"))
            createTaxon.setNull(7, Types.INTEGER);
          else createTaxon.setInt(7, new Integer(right).intValue());
              
        } else {
          createTaxon = conn.prepareStatement(
            "insert into taxon " 
            + "(ncbi_taxon_id, node_rank,"
            + " genetic_code, mito_genetic_code,"
            + " left_value, right_value) " 
            + "values (?, ?, ?, ?, ?, ?)"
          );
          createTaxon.setInt(1, ncbi_taxon_id);
          if (rank.equals("NULL")) 
            createTaxon.setNull(2, Types.VARCHAR);
          else createTaxon.setString(2, rank);
          if (gencode.equals("NULL"))
            createTaxon.setNull(3, Types.TINYINT);
          else createTaxon.setInt(3, new Integer(gencode).intValue());
          if (mitocode.equals("NULL"))
            createTaxon.setNull(4, Types.TINYINT);
          else createTaxon.setInt(4, new Integer(mitocode).intValue());
          if (left.equals("NULL"))
            createTaxon.setNull(5, Types.INTEGER);
          else createTaxon.setInt(5, new Integer(left).intValue());
          if (right.equals("NULL"))
            createTaxon.setNull(6, Types.INTEGER);
          else createTaxon.setInt(6, new Integer(right).intValue());
        }
        createTaxon.executeUpdate();
        createTaxon.close();
        taxon_id = helper.getInsertID(conn, "taxon", "taxon_id");
        putTaxonNames(conn, (Map) taxon.getAnnotation().getProperty(EbiFormat.PROPERTY_TAXON_NAMES), taxon_id);
      }
      trs.close();
      selectTaxon.close();
      return taxon_id;
    }


    /** Attempts to put all the names of a taxon into the database. There's no need
     * to access this method, it's only for internal use.
     * @param conn the connection to the database
     * @param names a map with the <code>name_class</code> as key and a <code>String</code>
     *   as value or a <code>Set</code> of <code>String</code>s, respectively.
     * @param taxon_id the internal <code>taxon_id</code> within the database.
     */
    private static void putTaxonNames(Connection conn, Map names, int taxon_id) throws SQLException 
    {
      if (names != null) {     	
        Iterator it = names.keySet().iterator();
        while (it.hasNext()) {
          String nameClass = (String) it.next();
        
          // changes:
          if (names.get(nameClass) instanceof Set) {
            Set all_names = (Set) names.get(nameClass);
            Iterator i = all_names.iterator();
            while(i.hasNext()) {
              String name = (String) i.next();
              PreparedStatement createTaxon = conn.prepareStatement(
                 "insert into taxon_name " 
                 + "(taxon_id, name, name_class) " 
                 + "values (?, ?, ?)"
              );
              createTaxon.setInt(1, taxon_id);
              createTaxon.setString(2, name);
              createTaxon.setString(3, nameClass);
              createTaxon.executeUpdate();
              createTaxon.close();
            }	
          } else {
            // normal again
                	
            String name = (String) names.get(nameClass);
            PreparedStatement createTaxon = conn.prepareStatement(
              "insert into taxon_name " 
              + "(taxon_id, name, name_class) " 
              + "values (?, ?, ?)"
            );
            createTaxon.setInt(1, taxon_id);
            createTaxon.setString(2, name);
            createTaxon.setString(3, nameClass);
            createTaxon.executeUpdate();
            createTaxon.close();
          }
        }
      }
    }


    /** Attempt to close the Statement. Continue on if there is a problem during the close. */
    public static void attemptClose(Statement statement) 
    {
      if (statement != null) try {
        statement.close();
      } catch (SQLException se) {
        se.printStackTrace();
      }
    }


    /** Attempt to close the ResultSet. Continue on if there is a problem during the close. */
    public static void attemptClose(ResultSet resultset) 
    {
      if (resultset != null) try {
        resultset.close();
      } catch (SQLException se) {
        se.printStackTrace();
      }
    }
   
    
  /** 
   * Attempts to get a Taxon object corresponding to the specified
   * name.
   *
   * @param conn the connection to the database
   * @param name the species scientific name
   * @return the corresponding Taxon (which may have already been
   * present in memory after an earlier retrieval), or null if the
   * Taxon could not be found in the database.
   * @throws BioRuntimeException
   */
  public static Taxon getTaxon(Connection conn, String name) throws BioRuntimeException
  {
    PreparedStatement statement = null;
    ResultSet rs = null;
    try {
      int taxon_id = 0;
      statement = conn.prepareStatement(
        "SELECT taxon.taxon_id " +
        "FROM   taxon, taxon_name " +
        "WHERE  taxon.taxon_id = taxon_name.taxon_id " +
        "  AND  taxon_name.name LIKE ?");
      statement.setString(1, name);
      rs = statement.executeQuery();
      if (rs.next()) {
        taxon_id = rs.getInt(1);
        if (rs.wasNull()) 
          taxon_id = 0;
      }
      rs.close();
      statement.close();

      return (taxon_id != 0) ? getDBTaxon(conn, taxon_id) : null;
    } catch (ChangeVetoException ex) {
      throw new BioRuntimeException("Couldn't manipulate in-memory taxonomy", ex);
    } catch (SQLException ex) {
      throw new BioRuntimeException("Error fetching taxonomy annotations", ex);
    } finally {
      attemptClose(rs);
      attemptClose(statement);
    }
  }

  
  /** Returns all the <code>scientific name</code>s, which are currently stored in the
    * database.
    * 
    * @param conn connection to the database
    * @return array of lexicographically sorted <code>String</code>s
    * @throws BioRuntimeException
    */
  public static String[] getAllScientificNames(Connection conn) throws BioRuntimeException
  {
    String scientific_names[] = null;
    PreparedStatement statement = null;
    ResultSet rs = null;
    try {
      int taxa_counter = 0, i=0;
      statement = conn.prepareStatement(
        "SELECT COUNT(distinct name) " +
        "FROM   taxon_name " +
        "WHERE  name_class LIKE 'scientific name' ");
      rs = statement.executeQuery();
      if (rs.next())
        taxa_counter = rs.getInt(1);
      scientific_names = new String[taxa_counter];
      statement = conn.prepareStatement(
        "SELECT distinct name " +
        "FROM   taxon_name " +
        "WHERE  name_class LIKE 'scientific name' " +
        "ORDER BY name ASC");
      rs = statement.executeQuery();
      while (rs.next()) {
        scientific_names[i++] = rs.getString(1);
        if (rs.wasNull()) i--;
      }
      rs.close();
      statement.close();

      return (taxa_counter != 0) ? scientific_names : null;
    } catch (SQLException exc) {
      exc.printStackTrace();
    }
    return null;
  }

  
  /** Returns the annotated properties of a given taxon. Only for internal use,
    * since this is needed to construct a fully annotated <code>Taxon</code> object,
    * which will be returned by one of the getTaxon... methods.
    * @param conn database connection
    * @param taxon the taxon to be annotated
    * @return a fully annotated taxon.
    * @throws BioRuntimeException
    */
  private static Taxon getProperties(Connection conn, Taxon taxon) throws BioRuntimeException
  {
    PreparedStatement ps = null;
    ResultSet rs = null;
    
    try {
      String rank = null, gencode = null, mitocode = null, left = null, right = null;
      ps = conn.prepareStatement(
        "SELECT node_rank, genetic_code, mito_genetic_code, left_value, right_value " +
        "FROM   taxon " +
        "WHERE  taxon_id = ?");
      ps.setInt(1, getTaxonID(conn, getRealScientificName(taxon)));
      rs = ps.executeQuery();
      if (rs.next()) {
        rank = rs.getString(1);
        if (!rs.wasNull()) taxon.getAnnotation().setProperty("rank", rank);
        gencode = rs.getString(2);
        if (!rs.wasNull()) taxon.getAnnotation().setProperty("genetic code id", gencode);
        mitocode = rs.getString(3);
        if (!rs.wasNull()) taxon.getAnnotation().setProperty("mitochondrial genetic code id", mitocode);
        left = String.valueOf(rs.getInt(4));
        if (!rs.wasNull()) taxon.getAnnotation().setProperty("left value", left);
	    right = String.valueOf(rs.getInt(5));
        if (!rs.wasNull()) taxon.getAnnotation().setProperty("right value", right);
      }  
    } catch (SQLException exc) {
      throw new BioRuntimeException(exc);
    } catch (IllegalArgumentException exc) {
      throw new BioRuntimeException(exc);
    } catch (ChangeVetoException exc) {
      throw new BioRuntimeException(exc);
    } catch (BioException exc) {
      exc.printStackTrace();
    } finally {
      attemptClose(rs);
      attemptClose(ps);
    }
    
    return taxon;
  }
  
  
  /** This returns the true scientific name of a given taxon, if there is one. 
    * This is necessary because if a taxon does not have a parent node,
    * a <code>TaxonFactory</code> gives the scientific name 'ROOT' and the
    * real scientific name (if there is one) is only stored in the taxon's
    * <code>EbiFormat</code>-annotation. This name 'ROOT' applies only for
    * in memory taxon objects. In the database the real name is stored.
    * 
    * @param t the taxon
    * @return Name of the taxon.
    */
  public static String getRealScientificName(Taxon t)
  {
    String tName = "root";
    if (t.getScientificName() != null) tName = t.getScientificName();
    if (tName.toLowerCase().equals("root") && (t.getAnnotation() != null)) 
      if (t.getAnnotation().containsProperty(EbiFormat.PROPERTY_TAXON_NAMES)) 
        if (((Map) t.getAnnotation().getProperty(EbiFormat.PROPERTY_TAXON_NAMES)).containsKey("scientific name"))
          tName = ((Map) t.getAnnotation().getProperty(EbiFormat.PROPERTY_TAXON_NAMES)).get("scientific name").toString();
    return tName;
  }

  
  /** Returns a <code>Set</code> containing all internal <code>taxon_id</code>s of the database.
    * @param conn database connection
    * @return all <code>taxon_id</code>s or internal use of the database.
    * @throws SQLException
    */
  private static Set ids(Connection conn) throws SQLException
  {
    Set ids = new HashSet();
    Statement statement = conn.createStatement();
    ResultSet rs = statement.executeQuery("SELECT taxon_id FROM taxon");
    while (rs.next()) 
      ids.add(new Integer(rs.getInt(1)));
    attemptClose(rs);
    statement.close();
    
    return ids;
  }

  
  /** Returns a <code>Set</code> of all NCBI-Taxon-IDs which are currently stored in
    * the database. So it is easy to proove if a taxon is stored in the database or
    * perform other operations.
    * @param conn database connection
    * @return a <code>Set</code> containing all NCBI-IDs.
    * @throws BioRuntimeException
    */
  public static Set NCBIids(Connection conn) throws BioRuntimeException
  {
    Set ids = new HashSet();
    Statement statement = null;
    ResultSet rs = null;
    try {
      statement = conn.createStatement();
      rs = statement.executeQuery("SELECT ncbi_taxon_id FROM taxon");
      while (rs.next())
        ids.add(new Integer(rs.getInt(1)));
    } catch (SQLException exc) {
      throw new BioRuntimeException(exc);
    } finally {
      try {
        rs.close();
        statement.close();
      } catch (SQLException exc1) {
        throw new BioRuntimeException(exc1);
      }
    }
    return ids;
  }
   
  
  /** Deletes a taxon specified by one of it's names with all it's different
    * names, annotations and sequences from the database. This cannot be undone.
    * The removed taxon will be returned.
    * 
    * @param conn database connection
    * @param helper the helper for the certain database system to be used.
    * @param name one of the taxon's names
    * @return the taxon, which was successfully removed and which is not longer
    *   stored in the database.
    * @throws BioException
    * @throws SQLException
    */
  public static Taxon removeTaxon(Connection conn, DBHelper helper, String name) throws SQLException, BioException
  {    
    //if (name.contains("'"))  name = name.replace("\'", "\\\'");
    //if (name.contains("\"")) name = name.replace("\"", "\\\"");

    return removeTaxon(conn, helper, getTaxonID(conn, name));
  }

  
  /** This actually performs the delete operations.
    *  
    * @param conn
    * @param helper
    * @param taxon_id
    * @return taxon
    * @throws BioRuntimeException
    */
  private static Taxon removeTaxon(Connection conn, DBHelper helper, int taxon_id) throws BioRuntimeException
  {
    Taxon taxon = getTaxon(conn, taxon_id);
    /*
     * Children which don't have a parent anymore -> they become their own parent.
     */
    try {
      Stack children = getChildrenOf(conn, taxon);
      while (!children.isEmpty()) try {
        Taxon child = (Taxon) children.pop();
        setParent(conn, child, child);
      } catch (BioRuntimeException exc) {
        exc.printStackTrace();
      }
    } catch (BioException exc) {
      exc.printStackTrace();
    }
    try {
      Statement delete = conn.createStatement();
      delete.executeUpdate("DELETE FROM taxon_name WHERE taxon_id = "+taxon_id);
      delete.close();
      delete = (Statement) conn.createStatement();
      delete.executeUpdate("DELETE FROM bioentry WHERE taxon_id = "+taxon_id);
      delete.close();
      delete = (Statement) conn.createStatement();        
      delete.executeUpdate("DELETE FROM taxon WHERE taxon_id = "+taxon_id);
      delete.close();
    } catch (SQLException exc) {
      throw new BioRuntimeException("Could not delete successfully", exc);
    }
      
    return taxon;
  }
  
  
  /** Deletes the taxon given by it's NCBI-Taxon-ID from the database and returns the
    * removed taxon.
    * @param conn database connection
    * @param helper the helper for the database
    * @param ncbi_id the ncbi-id
    * @return the taxon wich is not stored in the database anymore.
    * @throws BioRuntimeException
    * @throws BioException
    * @throws SQLException
    */
  public static Taxon removeTaxon(Connection conn, int ncbi_id, DBHelper helper) throws BioRuntimeException, SQLException, BioException
  {
    return removeTaxon(conn, helper, getTaxonID(conn, getRealScientificName(getTaxon(conn, ncbi_id))));
  }
  
  
  /** This changes the scientific name of the given taxon and stores the new name persistent
    * in the database.
    * @param conn database connection
    * @param taxon the taxon to be changed
    * @param newName the new scientific name
    * @return the changed taxon with the new scientific name.
    * @throws SQLException
    */
  public static Taxon setScientificName(Connection conn, Taxon taxon, String newName) throws SQLException 
  {    
    PreparedStatement updateName = conn.prepareStatement(
      "UPDATE taxon_name SET name = ? " +
      "WHERE  name LIKE ? AND " +
      "       name_class LIKE ?");
    updateName.setString(1, newName);
    updateName.setString(2, taxon.getScientificName());
    updateName.setString(3, "scientific name");
    updateName.execute();
    attemptClose(updateName);

    return getTaxon(conn, newName);
  }
  
  
  /** With this method the common name of the given taxon can be changed or created,
    * if there was none before. The new common name will be stored persitently.
    * @param conn database connection
    * @param taxon the taxon to be updated
    * @param newName the new common name
    * @return the updated taxon.
    * @throws BioException
    * @throws SQLException
    */
  public static Taxon setCommonName(Connection conn, Taxon taxon, String newName) throws BioException, SQLException 
  {
    boolean noCommonName = false;
    String name_class = "common name";
    
    /* A taxon does not necessaryly have to have a common name. Let's have a look,
     * if we need to perform an update or insert operation.
     */
    PreparedStatement question = conn.prepareStatement(
      "SELECT COUNT(name) " +
      "FROM   taxon_name " +
      "WHERE  name_class LIKE ? AND " +
      "       taxon_id = ?");
    question.setString(1, name_class);
    question.setInt(2, getTaxonID(conn, taxon.getScientificName()));
    ResultSet rs = question.executeQuery();
    
    if (rs.next()) {
      if (rs.getInt(1) > 0) {
        attemptClose(rs);
        attemptClose(question);
        PreparedStatement updateName = conn.prepareStatement(
          "UPDATE taxon_name SET name = ? " +
          "WHERE  name LIKE ? AND " +
          "       name_class LIKE ?");
        updateName.setString(1, newName);
        updateName.setString(2, taxon.getScientificName());
        updateName.setString(3, name_class);
        updateName.execute();
        attemptClose(updateName);
      } else noCommonName = true;      
    } 
    if (noCommonName) {      
      attemptClose(rs);
      attemptClose(question);
      PreparedStatement insert = conn.prepareStatement(
        "INSERT INTO taxon_name (taxon_id, name, name_class) VALUES (?, ?, ?)");
      try {
        insert.setInt(1, getTaxonID(conn, taxon.getScientificName()));
        insert.setString(2, newName);
        insert.setString(3, name_class);
        insert.execute();
      } catch (BioException exc) {
        throw exc;
      } catch (SQLException exc) {
        throw exc;
      } finally {
        attemptClose(insert); 
      }
    }
        
    return getTaxon(conn, taxon.getScientificName());
  }

    
  /** Returns the internal <code>taxon_id</code> of the taxon specified by the 
    * scientific name. Only for internal use.
    * @param conn database connection
    * @param scientificName the scientific name
    * @return the internal taxon_id of the taxon.
    * @throws SQLException
    * @throws BioException
    */
  private static int getTaxonID(Connection conn, String scientificName) throws SQLException, BioException
  {
    int taxon_id = -1;   
    
    PreparedStatement taxID  = conn.prepareStatement(
      "SELECT taxon_id "        +
      "FROM   taxon_name "      +
      "WHERE  name LIKE ? AND " +
      "       name_class LIKE ?");
    taxID.setString(1, scientificName);
    taxID.setString(2, "scientific name");
    ResultSet rs = taxID.executeQuery();
    if (rs.next()) {
      taxon_id = rs.getInt(1);
      attemptClose(rs);
      attemptClose(taxID);
      return taxon_id;
    }
    attemptClose(rs);
    attemptClose(taxID);
    throw new BioException("The database does not contain a taxon named "+scientificName+".");
  }

  
  /** Adds a new name of the given <code>nameClass</code> to the taxon. However, there must
    * be exactly one <code>scientific name</code> and maximal one <code>common
    * name</code>. Otherwise an <code>Exception</code> will be thrown.  
    * 
    * @param conn database connection
    * @param taxon the taxon to be updated
    * @param nameClass the name_class of the new name.
    * @param newName the new name.
    * @return the persistently updated taxon.
    * @throws BioException
    * @throws SQLException
    * @throws BioRuntimeException
    */
  public static Taxon addName(Connection conn, Taxon taxon, String nameClass, String newName) throws BioException, SQLException, BioRuntimeException
  {
    if (nameClass.equals("scientific name") || nameClass.equals("common name")) 
      throw new BioException("There can only be one "+nameClass+".");  
    
    PreparedStatement insert = conn.prepareStatement(
      "INSERT INTO taxon_name (taxon_id, name, name_class) VALUES (?, ?, ?)");
    try {
      insert.setInt(1, getTaxonID(conn, getRealScientificName(taxon)));
      insert.setString(2, newName);
      insert.setString(3, nameClass);
      insert.executeUpdate();
    } catch (BioException exc) {
      throw exc;
    } catch (SQLException exc) {
      throw exc;
    } finally {
      attemptClose(insert);
    }
    
    return getTaxon(conn, taxon.getScientificName());
  }
  
  
  /** Deletes the specified name from of the taxon from the database. 
    * The <code>scientific name</code> has to be uniqe, so this cannot be removed by
    * this method.
    * 
    * @param conn the database connection
    * @param helper the helper for the used database system
    * @param taxon the taxon to be updated
    * @param nameClass the name_class of the name to be removed
    * @param oldName the old name, which is not needed anymore.
    * @return the updated taxon.
    * @throws BioException
    * @throws SQLException
    */
  public static Taxon removeName(Connection conn, DBHelper helper, Taxon taxon, String nameClass, String oldName) throws BioException, SQLException
  {
    if (nameClass.equals("scientific name")) 
      throw new BioException("You can't delete the "+nameClass);
        
    PreparedStatement delete = conn.prepareStatement(
      "DELETE FROM taxon_name " +
      "WHERE  name_class LIKE ? AND " +
      "       name LIKE ?"); 
    delete.setString(1, nameClass);
    delete.setString(2, oldName);
    delete.executeUpdate();
    attemptClose(delete);
    
    return getTaxon(conn, taxon.getScientificName());
  }
  
  
  /** Returns all children of the specified taxon.
    * @param conn database connection
    * @param scientificName name of the taxon which children should be searched.
    * @return a <code>Stac</code>, which contains the children sorted by theire 
    *   <code>scientific name</code>s from top to the bottom.
    * @throws BioException
    * @throws SQLException
    */
  public static Stack getChildrenOf(Connection conn, String scientificName) throws BioException
  {
    Stack children = new Stack();
    PreparedStatement ps = null, ps2 = null;
    ResultSet rs = null, rs2 = null;
    try {
      ps = conn.prepareStatement(
        "SELECT taxon.taxon_id " +
        "FROM   taxon, taxon_name " +
        "WHERE  taxon_name.taxon_id = taxon.taxon_id     AND" +
        "       taxon_name.name_class LIKE 'scientific%' AND" +
        "       parent_taxon_id = ? " +
        "ORDER BY taxon_name.name DESC");
      //if (scientificName.contains("'"))  scientificName = scientificName.replace("\'", "\\\'");
      //if (scientificName.contains("\"")) scientificName = scientificName.replace("\"", "\\\"");
      ps.setInt(1, getTaxonID(conn, scientificName));
      rs = ps.executeQuery();
      while (rs.next()) {
        ps2 = conn.prepareStatement(
          "SELECT name             " +
          "FROM   taxon_name       " +
          "WHERE  taxon_id = ? AND " +
          "       name_class LIKE 'scientific%'");
        ps2.setInt(1, rs.getInt(1));
        rs2 = ps2.executeQuery();
        if (rs2.next()) { 
          children.push(getTaxon(conn, rs2.getString(1)));
          if (rs2.wasNull()) children.pop();
        }
      }
    } catch (SQLException exc) {
      throw new BioException(exc);
    } finally {
      attemptClose(ps);
      attemptClose(rs);
      attemptClose(ps2);
      attemptClose(rs2);
    }

    return children;  
  }

  
  /** Returns the children as a <code>Stack</code> of this given taxon.
    * 
    * @param conn database connection
    * @param t the parent taxon
    * @return a sorted <code>Stack</code> of the children, which might be empty,
    *   but not <code>null</code>.
    * @throws BioException
    */
  public static Stack getChildrenOf(Connection conn, Taxon t) throws BioException
  {
    String tName = t.getScientificName();
    if (tName.toLowerCase().equals("root") && (t.getAnnotation() != null))
      if (t.getAnnotation().containsProperty(EbiFormat.PROPERTY_TAXON_NAMES)) 
        if (((Map) t.getAnnotation().getProperty(EbiFormat.PROPERTY_TAXON_NAMES)).containsKey("scientific name"))
           tName = ((Map) t.getAnnotation().getProperty(EbiFormat.PROPERTY_TAXON_NAMES)).get("scientific name").toString();
    
    return getChildrenOf(conn, tName);
  }

  
  /** Updates a taxon and sets it's rank to the specified <code>String</code>.
    * @param conn database connection.
    * @param tdb taxon to be updated
    * @param rank the new rank (like 'kingdom', 'genus' or what ever)
    * @throws BioRuntimeException
    */
  public static void setRank(Connection conn, Taxon tdb, String rank) throws BioRuntimeException
  {
    PreparedStatement ps = null;
    try {
      ps = conn.prepareStatement(
        "UPDATE taxon SET node_rank = ? WHERE taxon_id = ?");
      ps.setString(1, rank);
      ps.setInt(2, getTaxonID(conn, getRealScientificName(tdb)));
      ps.executeUpdate();
    } catch (SQLException exc) {
      throw new BioRuntimeException(exc);
    } catch (BioException exc) {
      throw new BioRuntimeException(exc);
    } finally {
      attemptClose(ps);
    }
  }

  
  /** Removes the rank persistently from the taxon in the database.
    * @param conn
    * @param helper
    * @param tdb
    * @throws BioRuntimeException
    */
  public static void removeRank(Connection conn, DBHelper helper, Taxon tdb) throws BioRuntimeException
  {
    PreparedStatement ps = null;
    try {
      ps = conn.prepareStatement("UPDATE taxon SET node_rank = NULL WHERE taxon_id = ?");
      ps.setInt(1, getTaxonID(conn, getRealScientificName(tdb)));
      ps.executeUpdate();
    } catch (SQLException exc) {
      throw new BioRuntimeException(exc);
    } catch (BioException exc) {
      throw new BioRuntimeException(exc);
    } finally {
      attemptClose(ps);
    }
  }

  
  /** Updates the taxon in the database and sets its genetic code id to the specified value.
    * @param conn
    * @param tdb
    * @param id
    * @throws BioRuntimeException
    */
  public static void setGeneticCodeID(Connection conn, Taxon tdb, int id) throws BioRuntimeException
  {
    PreparedStatement ps = null;
    try {
      ps = conn.prepareStatement(
        "UPDATE taxon SET genetic_code = ? WHERE taxon_id = ?");
      ps.setInt(1, id);
      ps.setInt(2, getTaxonID(conn, getRealScientificName(tdb)));
      ps.executeUpdate();
    } catch (SQLException exc) {
      throw new BioRuntimeException(exc);
    } catch (BioException exc) {
      throw new BioRuntimeException(exc);
    } finally {
      attemptClose(ps);
    }
  }

  
  /** Deletes the genetic code annotation from the taxon in the database.
    * @param conn
    * @param helper
    * @param tdb
    * @throws BioRuntimeException
    */
  public static void removeGeneticCodeID(Connection conn, DBHelper helper, Taxon tdb) throws BioRuntimeException
  {
    PreparedStatement ps = null;
    try {
      ps = conn.prepareStatement("UPDATE taxon SET genetic_code = NULL WHERE taxon_id = ?");
      ps.setInt(1, getTaxonID(conn, getRealScientificName(tdb)));
      ps.executeUpdate();
    } catch (SQLException exc) {
      throw new BioRuntimeException(exc);
    } catch (BioException exc) {
      throw new BioRuntimeException(exc);
    } finally {
      attemptClose(ps);
    }
  }

  
  /** Updates the given taxon and sets it's so called mitochondrial genetic code id to
    * the specified value.
    * @param conn
    * @param tdb
    * @param id
    * @throws BioRuntimeException
    */
  public static void setMitochondrialGeneticCodeID(Connection conn, Taxon tdb, int id) throws BioRuntimeException
  {
    PreparedStatement ps = null;
    try {
      ps = conn.prepareStatement(
        "UPDATE taxon SET mito_genetic_code = ? WHERE taxon_id = ?");
      ps.setInt(1, id);
      ps.setInt(2, getTaxonID(conn, getRealScientificName(tdb)));
      ps.executeUpdate();
    } catch (SQLException exc) {
      throw new BioRuntimeException(exc);
    } catch (BioException exc) {
      throw new BioRuntimeException(exc);
    } finally {
      attemptClose(ps);
    }
  }

  
  /** Deletes the so called mitochondrial genetic code annotation from the given taxon.
    * @param conn
    * @param helper
    * @param tdb
    * @throws BioRuntimeException
    */
  public static void removeMitochondrialGeneticCodeID(Connection conn, DBHelper helper, Taxon tdb) throws BioRuntimeException
  {
    PreparedStatement ps = null;
    try {
      ps = conn.prepareStatement("UPDATE taxon SET mito_genetic_code = NULL WHERE taxon_id = ?");
      ps.setInt(1, getTaxonID(conn, getRealScientificName(tdb)));
      ps.executeUpdate();
    } catch (SQLException exc) {
      throw new BioRuntimeException(exc);
    } catch (BioException exc) {
      throw new BioRuntimeException(exc);
    } finally {
      attemptClose(ps);
    }
  }

  
  /** Updates the taxon and sets the left value to the specified value.
    * @param conn
    * @param tdb
    * @param left
    * @throws BioRuntimeException
    */
  public static void setLeftValue(Connection conn, Taxon tdb, int left) throws BioRuntimeException
  {
    PreparedStatement ps = null;
    try {
      ps = conn.prepareStatement(
        "UPDATE taxon SET left_value = ? WHERE taxon_id = ?");
      ps.setInt(1, left);
      ps.setInt(2, getTaxonID(conn, getRealScientificName(tdb)));
      ps.executeUpdate();
    } catch (SQLException exc) {
      throw new BioRuntimeException(exc);
    } catch (BioException exc) {
      throw new BioRuntimeException(exc);
    } finally {
      attemptClose(ps);
    }
  }

  
  /** Deletes the left value from the specified taxon in the database.
    * @param conn
    * @param helper
    * @param tdb
    * @throws BioRuntimeException
    */
  public static void removeLeftValue(Connection conn, DBHelper helper, Taxon tdb) throws BioRuntimeException
  {
    PreparedStatement ps = null;
    try {
      ps = conn.prepareStatement("UPDATE taxon SEQ left_value = NULL WHERE taxon_id = ?");
      ps.setInt(1, getTaxonID(conn, getRealScientificName(tdb)));
      ps.executeUpdate();
    } catch (SQLException exc) {
      throw new BioRuntimeException(exc);
    } catch (BioException exc) {
      throw new BioRuntimeException(exc);
    } finally {
      attemptClose(ps);
    }
  }

  
  /** Updates the taxon in the database and sets the right value to the specified value.
    * @param conn
    * @param tdb
    * @param right
    * @throws BioRuntimeException
    */
  public static void setRightValue(Connection conn, Taxon tdb, int right) throws BioRuntimeException
  {
    PreparedStatement ps = null;
    try {
      ps = conn.prepareStatement(
        "UPDATE taxon SET right_value = ? WHERE taxon_id = ?");
      ps.setInt(1, right);
      ps.setInt(2, getTaxonID(conn, getRealScientificName(tdb)));
      ps.executeUpdate();
    } catch (SQLException exc) {
      throw new BioRuntimeException(exc);
    } catch (BioException exc) {
      throw new BioRuntimeException(exc);
    } finally {
      attemptClose(ps);
    }
  }

  
  /** Deletes the right value from the specified taxon in the database.
    * @param conn
    * @param helper
    * @param tdb
    * @throws BioRuntimeException
    */
  public static void removeRightValue(Connection conn, DBHelper helper, Taxon tdb) throws BioRuntimeException
  {
    PreparedStatement ps = null;
    try {
      ps = conn.prepareStatement("UPDATE taxon SEQ right_value = NULL WHERE taxon_id = ?");
      ps.setInt(1, getTaxonID(conn, getRealScientificName(tdb)));
      ps.executeUpdate();
    } catch (SQLException exc) {
      throw new BioRuntimeException(exc);
    } catch (BioException exc) {
      throw new BioRuntimeException(exc);
    } finally {
      attemptClose(ps);
    }
  }

  
  /** This updates the taxonomic tree in the database and sets the parent of the given 
    * child taxon to the parent taxon.
    * @param conn
    * @param child
    * @param parent
    * @throws BioRuntimeException
    */
  public static void setParent(Connection conn, Taxon child, Taxon parent) throws BioRuntimeException
  {
    PreparedStatement ps = null;
    try {
      ps = conn.prepareStatement(
        "UPDATE taxon SET parent_taxon_id = ? WHERE taxon_id = ?");
      ps.setInt(1, getTaxonID(conn, getRealScientificName(parent)));
      ps.setInt(2, getTaxonID(conn, getRealScientificName(child)));
      ps.executeUpdate();
    } catch (SQLException exc) {
      throw new BioRuntimeException(exc);
    } catch (BioException exc) {
      throw new BioRuntimeException(exc);
    } finally {
      attemptClose(ps);
    }
  }
  
  
  
  /** This method tries to perform a complete update according to the given 
    * <code>TaxonFactory</code>, which already contains the newes taxa and the files
    * available at the NCBI-FTP-Site.
    * @param conn database connection
    * @param helper helper for the database system to be used
    * @param factory the TaxonFactory containing all the new taxa
    * @param delnodes file containing the ncbi taxon ids which don't exist anymore
    * @param merged file containing the ncbi taxon ids which were merged.
    * @throws IOException
    */
  public static void automaticUpdate(Connection conn, DBHelper helper, TaxonFactory factory, File delnodes, File merged) throws IOException
  {
   /*
    * Delnodes.
    */
    BufferedReader br = new BufferedReader(new FileReader(delnodes));
    Set ids = NCBIids(conn);
    while (br.ready()) {
      String line = br.readLine().split("\t")[0];
      if (ids.contains(line)) {
        try {
          removeTaxon(conn, Integer.parseInt(line), helper);
          ids.remove(line);
        } catch (NumberFormatException exc1) {
          exc1.printStackTrace();
        } catch (BioRuntimeException exc1) {
          exc1.printStackTrace();
        } catch (BioException exc1) {
          exc1.printStackTrace();
        } catch (SQLException exc) {
          exc.printStackTrace();
        }
      } 
    }
               
    
    /*
     * TaxonFactory.
     */
    /*
     * Mergenodes.
     */
    br = new BufferedReader(new FileReader(merged));
    while (br.ready()) {
      String[] line = br.readLine().split("\t");
      if (ids.contains(line[0])) {
        // line[0] merged with [2].
        //System.out.println(line[0]+"\t"+line[2]);
        /*
         * Change children
         */
        // Only possible if the new parent exists in the database
        if (ids.contains(line[2])) try {
          Taxon parent = getTaxon(conn, Integer.parseInt(line[2]));
          Stack children = getChildrenOf(conn, getTaxon(conn, Integer.parseInt(line[0])));
          while (!children.isEmpty()) {
            setParent(conn, (Taxon) children.pop(), parent); 
          }
         /*
          * Sequences to the new node.
          */
          PreparedStatement ps = conn.prepareStatement(
            "UPDATE bioentry SET taxon_id = ? WHERE taxon_id = ?");
          ps.setInt(1, getTaxonID(conn, getRealScientificName(parent)));
          ps.setInt(2, getTaxonID(conn, getRealScientificName(getTaxon(conn, Integer.parseInt(line[0])))));
          ps.executeUpdate();
          ps.close();
        } catch (NumberFormatException exc1) {
          exc1.printStackTrace();
        } catch (BioRuntimeException exc1) {
          exc1.printStackTrace();
        } catch (BioException exc1) {
          exc1.printStackTrace();
        } catch (SQLException exc) {
          exc.printStackTrace();
        }
        /*
         * Delete
         */
        try {
          removeTaxon(conn, Integer.parseInt(line[0]), helper);
        } catch (NumberFormatException exc2) {
          exc2.printStackTrace();
        } catch (BioRuntimeException exc2) {
          exc2.printStackTrace();
        } catch (BioException exc2) {
          exc2.printStackTrace();
        } catch (SQLException exc) {
          exc.printStackTrace();
        }
        
      }
    }
    
    
    try {
      /*
       * TaxonFactory.
       */
      for (Iterator i=ids(conn).iterator(); i.hasNext(); ) {
        Taxon tDB = getTaxon(conn, ((Integer) i.next()).intValue());
        Taxon tFA = factory.search(getRealScientificName(tDB));
        boolean change = false;
        // This cannot work because the TaxonFactory doesn't provide this information.
        if ((tFA == null) && (tDB.getAnnotation().containsProperty(EbiFormat.PROPERTY_NCBI_TAXON))) 
          tFA = factory.search(tDB.getAnnotation().getProperty(EbiFormat.PROPERTY_NCBI_TAXON));
        if (tFA == null) continue;        
        
        // scientific name
        if (tFA.getScientificName() != null) {
          if (!tFA.getScientificName().equals(getRealScientificName(tDB))) try {
            setScientificName(conn, tDB, tFA.getScientificName());
          } catch (SQLException exc2) {
            exc2.printStackTrace();
          }
        }
        
        // common name
        if (tFA.getCommonName() != null) {
          if (tDB.getCommonName() == null) change = true;
          else if (!tDB.getCommonName().equals(tFA.getCommonName())) change = true;
          if (change) try {
            setCommonName(conn, tDB, tFA.getCommonName());
          } catch (BioException exc1) {
            exc1.printStackTrace();
          } catch (SQLException exc) {
            exc.printStackTrace();
          } finally { change = false; }
        } else if (tDB.getCommonName() != null) try {
          removeName(conn, helper, tDB, "common name", tDB.getCommonName());
        } catch (BioException exc1) {
          exc1.printStackTrace();
        } catch (SQLException exc) {
          exc.printStackTrace();
        }
        
        // annotation
        Annotation annoFA = tFA.getAnnotation();
        Annotation annoDB = tDB.getAnnotation();
        // rank
        if (annoFA.containsProperty("rank")) {
          if (!annoDB.containsProperty("rank")) change = true;
          else if (!annoDB.getProperty("rank").equals(annoFA.getProperty("rank"))) change = true;
          if (change) setRank(conn, tDB, annoFA.getProperty("rank").toString());
          change = false;
        } else if (annoDB.containsProperty("rank")) 
          removeRank(conn, helper, tDB);
        // genetic code id
        if (annoFA.containsProperty("genetic code id")) {
          if (!annoDB.containsProperty("genetic code id")) change = true;
          else if (!annoDB.getProperty("genetic code id").equals(annoFA)) change = true;
          if (change) setGeneticCodeID(conn, tDB, Integer.parseInt(annoFA.getProperty("genetic code id").toString()));
          change = false;
        } else if (annoDB.containsProperty("genetic code id"))
          removeGeneticCodeID(conn, helper, tDB);
        // mitochondrial genetic code id
        if (annoFA.containsProperty("mitochondrial genetic code id")) {
          if (!annoDB.containsProperty("mitochondrial genetic code id")) change = true;
          else if (!annoDB.getProperty("mitochondrial genetic code id").equals(annoFA.getProperty("mitochondrial genetic code id"))) change = true;
          if (change) setMitochondrialGeneticCodeID(conn, tDB, Integer.parseInt(annoFA.getProperty("mitochondrial genetic code id").toString()));
          change = false;
        } else if (annoDB.containsProperty("mitochondrial genetic code id"))
          removeMitochondrialGeneticCodeID(conn, helper, tDB);
        // left value
        if (annoFA.containsProperty("left value")) {
          if (!annoDB.containsProperty("left value")) change = true;
          else if (!annoDB.getProperty("left value").equals(annoFA.getProperty("left value"))) change = true;
          if (change) setLeftValue(conn, tDB, Integer.parseInt(annoFA.getProperty("left value").toString()));
          change = false;
        } else if (annoDB.containsProperty("left value"))
          removeLeftValue(conn, helper, tDB);
        // right value
        if (annoFA.containsProperty("right value")) {
          if (!annoDB.containsProperty("right value")) change = true;
          else if (!annoDB.getProperty("right value").equals(annoFA.getProperty("right value"))) change = true;
          if (change) setRightValue(conn, tDB, Integer.parseInt(annoFA.getProperty("right value").toString()));
          change = false;
        } else if (annoDB.containsProperty("right value"))
          removeRightValue(conn, helper, tDB);
        // other names:
        if (annoFA.containsProperty(EbiFormat.PROPERTY_TAXON_NAMES)) {
          if (annoFA.getProperty(EbiFormat.PROPERTY_TAXON_NAMES) instanceof Map) {
            Map names = (Map) annoFA.getProperty(EbiFormat.PROPERTY_TAXON_NAMES);
            Map dbNames = (Map) annoDB.getProperty(EbiFormat.PROPERTY_TAXON_NAMES);       
          
          /*
           * delete old names, if not valid anymore.
           */
          for (Iterator j=dbNames.keySet().iterator(); j.hasNext(); ) {
            String key = j.next().toString();
            if (dbNames.get(key) instanceof Set) {
              if (!names.containsKey(key)) // delete everything 
                for (Iterator k=((Set) dbNames.get(key)).iterator(); k.hasNext(); ) try {
                  String elem = k.next().toString();
                  removeName(conn, helper, tDB, key, elem);
                } catch (BioException exc2) {
                  exc2.printStackTrace();
                } catch (SQLException exc) {
                  exc.printStackTrace();
                } 
              // this key is still in the new taxon.
              else for (Iterator k=((Set) dbNames.get(key)).iterator(); k.hasNext(); ) {
                String elem = k.next().toString();
                // many entries for this key in the new taxon.
                if (names.get(key) instanceof Set) {
                  if (!((Set) names.get(key)).contains(elem)) try {
                    removeName(conn, helper, tDB, key, elem);
                  } catch (BioException exc3) {
                    exc3.printStackTrace();
                  } catch (SQLException exc) {
                    exc.printStackTrace();
                  }
                 
                // in the new taxon only one entry for this key
                } else if (!names.get(key).equals(elem)) try {
                  removeName(conn, helper, tDB, key, elem);
                } catch (BioException exc2) {
                  exc2.printStackTrace();
                } catch (SQLException exc) {
                  exc.printStackTrace();
                }
              }
            } else { // only one entry in the new animal
              if (!names.containsKey(key)) change = true; 
              // many entries in the new animal.
              else if (names.get(key) instanceof Set) {
                if (!((Set) names.get(key)).contains(dbNames.get(key))) change = true;
              // only one entry in the new animal (what ever)
              } else if (!names.get(key).equals(dbNames.get(key))) change = true;
              if (change) try { 
                removeName(conn, helper, tDB, key, dbNames.get(key).toString());
              } catch (BioException exc2) {
                exc2.printStackTrace();
              }
              change = false;
            }
          }
          // update:
          dbNames = (Map) getTaxon(conn, getRealScientificName(tDB)).getAnnotation().getProperty(EbiFormat.PROPERTY_TAXON_NAMES);
          
          /*
           * New names.
           */
          for (Iterator j=names.keySet().iterator(); j.hasNext(); ) {
            String key = j.next().toString();
            // there're multiple names for this key
            if (names.get(key) instanceof Set) {
              for (Iterator k=((Set) names.get(key)).iterator(); k.hasNext(); ) {
                String elem = k.next().toString();
                if (!dbNames.containsKey(key)) change = true;
                else if (dbNames.get(key) instanceof Set) {
                  if (!((Set) dbNames.get(key)).contains(elem)) change = true;
                } else if (!dbNames.get(key).equals(elem)) change = true; 
                if (change) try {
                  addName(conn, tDB, key, elem);
                } catch (BioRuntimeException exc3) {
                  exc3.printStackTrace();
                } catch (BioException exc3) {
                  exc3.printStackTrace();
                } catch (SQLException exc) {
                  exc.printStackTrace();
                }
                change = false;
              }
            // there's only one name for this key.
            } else if (!dbNames.containsKey(key)) try {
                addName(conn, tDB, key, names.get(key).toString());
              } catch (BioRuntimeException exc2) {
                exc2.printStackTrace();
              } catch (BioException exc2) {
                exc2.printStackTrace();
              } catch (SQLException exc) {
                exc.printStackTrace();
              }
            }
          }
        }
      }
    } catch (NumberFormatException exc) {
      exc.printStackTrace();
    } catch (BioRuntimeException exc) {
      exc.printStackTrace();
    } catch (NoSuchElementException exc) {
      exc.printStackTrace();
    } catch (SQLException exc) {
      exc.printStackTrace();
    }
  }

  
}
