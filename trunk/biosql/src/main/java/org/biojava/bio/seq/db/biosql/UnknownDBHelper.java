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
import java.sql.SQLException;

/**
 * @author Thomas Down
 * @author Matthew Pocock
 * @deprecated Use hibernate and org.biojavax.bio.db.*
 */
public class UnknownDBHelper extends DBHelper {
    public int getInsertID(Connection conn,
			   String table,
			   String columnName)
	throws SQLException
    {
	throw new SQLException("No DBHelper for this database -- can't write objects");
    }

    public DeleteStyle getDeleteStyle() {
	return DELETE_GENERIC;
    }
}
