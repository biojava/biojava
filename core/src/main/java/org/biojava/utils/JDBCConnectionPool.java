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
package org.biojava.utils;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.LinkedList;

/**
 * Really simple connection pool for JDBC databases.
 *
 * <h2>Use:</h2>
 * <pre>
 * JDBCConnectionPool pool = new JDBCConnectionPool(jdbcURL, userName, passwd);
 * ...
 * Connection conn = pool.takeConnection();
 * // do stuff with conn
 * pool.putConnection(conn);
 * // don't use conn from here on
 *
 * Statement stmt = pool.takeStatement();
 * // do stuff with stmt
 * pool.putStatement(stmt);
 * // don't do anything else with stmt
 * </pre>
 *
 * <p>It is not a good idea to call <code>close()</code> on a connection you
 * get from a pool. This would prevent it from being re-used. Also, we have
 * seen some odd behavior with connections involved in transactions being
 * re-used. We have not yet identified exactly how you can safely use a
 * pooled connection for transaction-safe code.</p>  
 *
 * <p><em>Note:</em> We should probably be moving to a propper connection pool
 * API. Let's standardise on one soon.</p>
 *
 * @author Thomas Down
 * @author Matthew Pocock
 */

public class JDBCConnectionPool {
  private final String dbURL;
  private final String dbUser;
  private final String dbPass;

  private LinkedList connectionPool;

  {
    connectionPool = new LinkedList();
  }

  public JDBCConnectionPool(String url, String user, String pass)
  {
    dbURL = url;
    dbUser = user;
    dbPass = pass;
  }

  public JDBCConnectionPool(String url)
  {
    this(url, null, null);
  }

  //
  // Manage a pool of transactions with the database.
  //

  public Connection takeConnection()
          throws SQLException
  {
    Connection conn = null;

    synchronized (connectionPool) {
      if (connectionPool.size() > 0) {
        conn = (Connection) connectionPool.removeFirst();
      }
    }

    // We don't perform the isClosed in the synchronized block in case the
    // network is being slow.

    //if (conn != null) {
    //    if (!conn.isClosed()) {
    //            // hack for turfing out bad connections
    //            Statement stmt = conn.createStatement();
    //            stmt.execute("SELECT 1");
    //            return conn;
    //
    //    } else {
    //	// We simply drop conn on the floor.  It should be safely collected.
    //	return takeConnection();
    //    }
    //}

    // Statement-pool was empty -- let's create a new connection.

    if(dbUser != null) {
      conn = DriverManager.getConnection(dbURL, dbUser, dbPass);
    } else {
      conn = DriverManager.getConnection(dbURL);
    }
    Statement st = conn.createStatement();
    // st.execute("SET OPTION SQL_BIG_SELECTS=1");
    st.close();
    return conn;
  }

  public void putConnection(Connection c)
          throws SQLException
  {
    if(!c.isClosed()) {
      synchronized (connectionPool) {
        connectionPool.addLast(c);
      }
    }
  }

  public Statement takeStatement()
          throws SQLException
  {
    return takeConnection().createStatement();
  }

  public void putStatement(Statement st)
          throws SQLException
  {
    putConnection(st.getConnection());
    st.close();
  }
}
