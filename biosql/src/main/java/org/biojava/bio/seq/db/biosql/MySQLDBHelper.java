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
import java.sql.DatabaseMetaData;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;

/**
 * This is a <code>DBHelper</code> that provides support for MySQL
 * databases.
 *
 * @author Thomas Down
 * @author Matthew Pocock
 * @author Len Trigg
 * @deprecated Use hibernate and org.biojavax.bio.db.*
 */
public class MySQLDBHelper extends DBHelper {

    private final DeleteStyle mDeleteStyle;

    public MySQLDBHelper(Connection connection) {
        DeleteStyle deleteStyle = DELETE_GENERIC;
        try {
            DatabaseMetaData metadata = connection.getMetaData();
            int major = metadata.getDatabaseMajorVersion();
						// Minor version irrelevant as returns 0 for 4.0.*
						// Actually need subminor version which not in metadata
						// Could parse from getDatabaseProductVersion() string if really needed
            int minor = metadata.getDatabaseMinorVersion();
            if ((major > 4) || ((major == 4) && (minor >= 0))) {
                deleteStyle = DELETE_MYSQL4;
            }
        } catch (SQLException e) {
            System.err.println("Exception getting DatabaseMetaData:" +  e.getMessage());
            // Stick with generic style
        }
        mDeleteStyle = deleteStyle;
    }

    // Inherit docs
    public int getInsertID(Connection conn,
			   String table,
			   String columnName)
	throws SQLException
    {
        Statement st = conn.createStatement();
	ResultSet rs = st.executeQuery("select last_insert_id()");
	int id = -1;
	if (rs.next()) {
	    id = rs.getInt(1);
	}
        rs.close();
	st.close();
	
	if (id < 1) {
	    throw new SQLException("Couldn't get last_insert_id()");
	}
	return id;
    }

    // Inherit docs
    public DeleteStyle getDeleteStyle() {
        return mDeleteStyle;
    }
}
