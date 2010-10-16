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
import java.sql.SQLException;

import javax.sql.DataSource;

import org.biojava.bio.BioRuntimeException;

/**
 * Isolates all code that is specific to a particular RDBMS. To add
 * support for a new RDBMS, write a new <code>DBHelper</code> subclass
 * and ensure that it can be found by editing the
 * <code>getDBHelperForURL</code> method in this class.
 *
 * @author Thomas Down
 * @author Matthew Pocock
 * @author Len Trigg
 * @author Eric Haugen
 * @author Richard Holland
 * @deprecated Use hibernate and org.biojavax.bio.db.*
 */
public abstract class DBHelper {

    /**
     * Returns a DBHelper implementation suitable for a particular
     * database.
     *
     * @param conn a connection to the database.
     * @return a <code>DBHelper</code>.
     */
    public static DBHelper getDBHelper(Connection conn) {
        try {
            String dbType = conn.getMetaData().getURL();
            if (dbType.startsWith("jdbc:")) {
                dbType = dbType.substring(5);
            }
            if (!Character.isLetter(dbType.charAt(0))) {
                throw new IllegalArgumentException("URL must start with a letter: " + dbType);
            }
            
            int colon = dbType.indexOf(':');
            if (colon > 0) {
                String protocol = dbType.substring(0, colon);
                if (protocol.indexOf("mysql") >= 0) {
                    // Accept any string containing `mysql', to cope with Caucho driver
                    return new MySQLDBHelper(conn);
                } else if (protocol.equals("postgresql")) {
                    return new PostgreSQLDBHelper();
                } else if (protocol.equals("oracle")) {
                    return new OracleDBHelper(conn);
                } else if (protocol.equals("hsqldb")) {
                    return new HypersonicDBHelper();
                }
            }
        } catch (SQLException se) {
            se.printStackTrace();
        }
        return new UnknownDBHelper();
    }

    public static final class DeleteStyle {
        private final String name;

        private DeleteStyle(String name) {
            this.name = name;
        }

        public String toString() {
            return "DBHelper.DeleteStyle: " + name;
        }
    }

    public static final DeleteStyle DELETE_POSTGRESQL = new DeleteStyle("Postgresql");
    public static final DeleteStyle DELETE_MYSQL4 = new DeleteStyle("Mysql 4.0.* or later");
    public static final DeleteStyle DELETE_GENERIC = new DeleteStyle("Portable SQL");


    public static final class JoinStyle {
        private final String name;

        private JoinStyle(String name) {
            this.name = name;
        }

        public String toString() {
            return "DBHelper.JoinStyle: " + name;
        }
    }

    public static final JoinStyle JOIN_ORACLE8 = new JoinStyle("Oracle 8i or earlier");
    public static final JoinStyle JOIN_GENERIC = new JoinStyle("Portable SQL");


    public static final class BioSequenceStyle {
        private final String name;

        private BioSequenceStyle(String name) {
            this.name = name;
        }

        public String toString() {
            return "DBHelper.BioSequenceStyle: " + name;
        }
    }

    public static final BioSequenceStyle BIOSEQUENCE_GENERIC = new BioSequenceStyle("Standard schema (except for Oracle schemas using CLOBs)");
    public static final BioSequenceStyle BIOSEQUENCE_ORACLECLOB = new BioSequenceStyle("Oracle schema using CLOBS (but not Len Trigg's schema)");

    /**
     * Returns the id value created during the last insert
     * command. This is for tables that have an auto increment column.
     * 
     * @return the last id assigned, or -1 if the id could not be
     * found.
     */
    public abstract int getInsertID(Connection conn, String table, 
                                    String columnName) throws SQLException;


    /**
     * Returns the an object indicating the style of deletion that
     * this database should employ. Unless overridden, this is
     * DELETE_GENERIC.
     * 
     * @return the preferred deletion style.
     */
    public DeleteStyle getDeleteStyle() {
        return DELETE_GENERIC;
    }


    /**
     * Returns the an object indicating the style of table joining that
     * this database should employ.
     * 
     * @return the preferred joining style.
     */
    public JoinStyle getJoinStyle() {
        return JOIN_GENERIC;
    }


    /**
     * Returns the an object indicating the style of biosequence storage
     * that this database should employ. Generally, leave it at the default
     * unless you are using the Oracle schema, in which case you need
     * to override it to return BIOSEQUENCE_ORACLECLOB. This is because, in the
     * Oracle schema we need to use CLOBs (except when using Len Trigg's 
     * version which uses LONGs instead.)
     * 
     * @return the preferred joining style.
     */
    public BioSequenceStyle getBioSequenceStyle() {
        return BIOSEQUENCE_GENERIC;
    }

    /**
     * Detects whether a particular table is present in the database.
     *
     * @param ds a <code>DataSource</code> that can provide a connection to a database 
     * @param tablename the name of the table.
     * @return true if the table exists in the database.
     * @throws NullPointerException if pool is null.
     * @throws IllegalArgumentException if tablename is null or empty.
     */
    public boolean containsTable(DataSource ds, String tablename) {
        if (ds == null) {
            throw new NullPointerException("Require a datasource.");
        }
        if ((tablename == null) || (tablename.length() == 0)) {
            throw new IllegalArgumentException("Invalid table name given");
        } 
        Connection conn = null;
        try {
            boolean present;
            PreparedStatement ps = null;
            try {
                conn = ds.getConnection();
                ps = conn.prepareStatement("select * from " + tablename + " limit 1");
                ps.executeQuery();
                present = true;
            } catch (SQLException ex) {
                present = false;
            }
            if (ps != null) {
                ps.close();
            }
            if (conn != null) {
                conn.close();
            }
            return present;
        } catch (SQLException ex) {
            if (conn!=null) try {conn.close();} catch (SQLException ex3) {}
            throw new BioRuntimeException(ex);
        }
    }
}