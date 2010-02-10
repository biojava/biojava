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
import java.lang.reflect.Method;
import java.sql.Clob;
import java.sql.Connection;
import java.sql.DatabaseMetaData;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;

import javax.sql.DataSource;

import org.biojava.bio.BioRuntimeException;

/**
 * This is a <code>DBHelper</code> that provides support for Oracle
 * databases.
 *
 * @author Len Trigg
 * @author Eric Haugen
 * @author Richard Holland
 * @deprecated Use hibernate and org.biojavax.bio.db.*
 */
public class OracleDBHelper extends DBHelper {
    
    private final JoinStyle mJoinStyle;
    private final BioSequenceStyle mBseqStyle;
    
    public OracleDBHelper(Connection connection) {
        JoinStyle joinStyle = JOIN_GENERIC;
        BioSequenceStyle bseqStyle = BIOSEQUENCE_GENERIC;
        try {
            DatabaseMetaData metadata = connection.getMetaData();
            String version = metadata.getDatabaseProductVersion();
            if ((version != null) && version.startsWith("Oracle8")) {
                joinStyle = JOIN_ORACLE8;
            }
            // Describe the biosequence table
            Statement st = null;
            ResultSet rs = null;
            try {
                // For CLOB access, the Oracle 9i (or better) JDBC drivers are required on the ClassPath.
                // This simple test for the BioSequence seq storage type makes some basic assumptions:
                //    1. That if you are using Len Trigg's schema, you are logged into it directly and
                //       not using views onto it.
                //    2. That if you are logged into any schema using views, then they are views onto
                //       the standard CLOB based schema.
                st = connection.createStatement();
                rs = st.executeQuery("select data_type from user_tab_columns where table_name='BIOSEQUENCE' and column_name='SEQ'");
                String seqType = null;
                if (rs.next()) {
                    seqType = rs.getString(1);
                }
                // If it's missing or says CLOB, then use the CLOB interfaces.
                // Else, use BIOSEQUENCE_GENERIC (and assume it allows normal get/set calls)
                if (seqType==null || "CLOB".equals(seqType)) bseqStyle = BIOSEQUENCE_ORACLECLOB;
            } finally {
                if (rs != null) try { rs.close(); } catch (SQLException se) { }
                if (st != null) try { st.close(); } catch (SQLException se) { }
            }
        } catch (SQLException e) {
            System.err.println("Exception getting DatabaseMetaData:" +  e.getMessage());
            // Stick with generic style
        }
        mJoinStyle = joinStyle;
        mBseqStyle = bseqStyle;
    }
    
    
    // Inherit docs
    public JoinStyle getJoinStyle() {
        return mJoinStyle;
    }
    
    
    // Inherit docs
    public int getInsertID(Connection conn, String table, String columnName) throws SQLException {
        Statement st = null;
        ResultSet rs = null;
        try {
            st = conn.createStatement();
            // We assume that the Oracle BioSQL schema uses sequences for the autoincrement fields,
            // one sequence per table.
            rs = st.executeQuery("select " + table + "_pk_seq.CURRVAL from dual");
            int id = -1;
            if (rs.next()) {
                id = rs.getInt(1);
            }
            
            if (id < 1) {
                throw new SQLException("Couldn't get last insert id");
            }
            return id;
        } finally {
            if (rs != null) try { rs.close(); } catch (SQLException se) { }
            if (st != null) try { st.close(); } catch (SQLException se) { }
        }
    }
    
    // Inherit docs
    public boolean containsTable(DataSource ds, String tablename) {
        if (ds == null) {
            throw new NullPointerException("Require a datasource.");
        }
        if ((tablename == null) || (tablename.length() == 0)) {
            throw new IllegalArgumentException("Invalid table name given");
        }
        //System.err.println("Checking for table existence: " + tablename);
        Connection conn = null;
        try {
            boolean present;
            conn = ds.getConnection();
            PreparedStatement ps = conn.prepareStatement("select rownum from " + tablename + " where rownum < 1");
            try {
                ps.executeQuery();
                present = true;
            } catch (SQLException ex) {
                //System.err.println("Table " + tablename + " does not exist.");
                present = false;
            } finally {
                ps.close();
                if (conn != null) {
                    conn.close();
                }
            }
            return present;
        } catch (SQLException ex) {
            if (conn!=null) try {conn.close();} catch (SQLException ex3) {}
            throw new BioRuntimeException(ex);
        }
    }
    
    // Inherit docs
    public BioSequenceStyle getBioSequenceStyle() {
        return mBseqStyle;
    }
    
    /*
     * Use this to retrieve a CLOB value.
     * @param conn a connection to an Oracle database.
     * @param rs the ResultSet to retrieve the CLOB from.
     * @param column the number of the column in the ResultSet that the CLOB lives in.
     * @return String value of the CLOB.
     */
    public String clobToString(Connection conn, ResultSet rs, int column) {
        try {
            Clob seqclob = rs.getClob(column);
            StringBuffer buf = new StringBuffer();
            int bufSize = 1024;
            long start = 1L;
            long remain = seqclob.length();
            while (remain>0L) {
                if (bufSize>remain) bufSize=(int)remain;
                buf.append(seqclob.getSubString(start,bufSize));
                start+=bufSize;
                remain-=bufSize;
            }
            return buf.toString().trim();
        } catch (Exception ex) {
            throw new BioRuntimeException(ex);
        }
    }
    
    /*
     * Use this to set a CLOB value. OJDBC version 9i must be on the ClassPath.
     * @param conn a connection to an Oracle database.
     * @param rs the ResultSet to retrieve the CLOB from.
     * @param column the number of the column in the ResultSet that the CLOB lives in.
     * @param the value to set to the CLOB.
     */
    public void stringToClob(Connection conn, ResultSet rs, int column, String value) {
        try {
            // Can't use oracle.sql.CLOB directly as we'd need it at compile time otherwise.
            Class clob = Class.forName("oracle.sql.CLOB");
            Method putString = clob.getDeclaredMethod("putString",new Class[]{long.class,String.class});
            // Only get here if we have some data to write.
            if (value==null) value=""; // To stop null pointer exceptions. End result is the same.
            putString.invoke(rs.getClob(column), new Object[]{new Long(1L),value});
        } catch (Exception ex) {
            throw new BioRuntimeException(ex);
        }
    }
    
}
