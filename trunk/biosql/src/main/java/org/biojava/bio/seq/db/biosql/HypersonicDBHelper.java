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
import java.sql.Statement;

import javax.sql.DataSource;

import org.biojava.bio.BioRuntimeException;

/**
 * This is a <code>DBHelper</code> that provides support for the
 * Hypersonic RDBMS. See the <a href="http://hsqldb.sourceforge.net/">HSQLDB home page</a>
 *
 * @author Len Trigg
 * @author Richard Holland
 * @deprecated Use hibernate and org.biojavax.bio.db.*
 */
public class HypersonicDBHelper extends DBHelper {

    // Inherit docs
    public int getInsertID(Connection conn, String table, String columnName) throws SQLException {
        Statement st = null;
        ResultSet rs = null;
        try {
            st = conn.createStatement();
            rs = st.executeQuery("call identity()");
            int id = -1;
            if (rs.next()) {
                id = rs.getInt(1);
            }
            if (id < 0) {
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
            PreparedStatement ps = null;
            try {
                ps = conn.prepareStatement("select top 1 * from " + tablename);
                ps.executeQuery();
                present = true;
            } catch (SQLException ex) {
                //System.err.println("Table " + tablename + " does not exist.");
                present = false;
            } finally {
                if (ps != null) {
                    ps.close();
                }
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

}
