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

package org.biojavax.bio.db.biosql;
import java.lang.reflect.Method;

import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.SymbolList;
import org.biojavax.CrossRef;
import org.biojavax.CrossReferenceResolver;
import org.biojavax.Namespace;
import org.biojavax.RichObjectFactory;
import org.biojavax.SimpleNamespace;
import org.biojavax.bio.BioEntry;
import org.biojavax.bio.seq.InfinitelyAmbiguousSymbolList;
import org.biojavax.bio.seq.RichSequence;

/**
 * A simple implementation of CrossReferenceResolver
 * @author Richard Holland
 * @author Mark Schreiber
 * @author David Scott
 * @since 1.5
 */
public class BioSQLCrossReferenceResolver implements CrossReferenceResolver {
    
    private Object session;
    private Method createQuery;
    private Method setParameter;
    private Method uniqueResult;
    
    /**
     * Requires a Hibernate session to work correctly. The session parameter
     * is a Hibernate Session object and must not be null. It is this session
     * that database objects will be retrieved from/persisted to.
     * @see <a href="http://www.hibernate.org/hib_docs/v3/api/org/hibernate/Session.html"> org.hibernate.Session</a>
     */
    public BioSQLCrossReferenceResolver(Object session) {
        try {
            // Lazy load the Session class from Hibernate.
            Class hibernateSession = session.getClass();
            Class realHibernateSession = Class.forName("org.hibernate.Session");
            // Test to see if our parameter is really a Session
            if (!realHibernateSession.isAssignableFrom(hibernateSession))
                throw new IllegalArgumentException("Parameter must be a org.hibernate.Session object");
            this.session = session;
            // Lookup the createQuery method
            this.createQuery = hibernateSession.getMethod("createQuery", new Class[]{String.class});
            // Lazy load the Query class from Hibernate.
            Class hibernateQuery = Class.forName("org.hibernate.Query");
            // Lookup the setParameter and uniqueQuery methods
            this.setParameter = hibernateQuery.getMethod("setParameter", new Class[]{int.class,Object.class});
            this.uniqueResult = hibernateQuery.getMethod("uniqueResult", new Class[]{});
        } catch (ClassNotFoundException e) {
            throw new RuntimeException(e);
        } catch (NoSuchMethodException e) {
            throw new RuntimeException(e);
        }
    }
    
    /**
     * {@inheritDoc}
     */
    public SymbolList getRemoteSymbolList(CrossRef cr, Alphabet a) {
        BioEntry be = this.getRemoteBioEntry(cr);
        if (be instanceof RichSequence) return (RichSequence)be;
        // If we get here we didn't find it, so we must create a dummy sequence instead
        if (!(a instanceof FiniteAlphabet)) throw new IllegalArgumentException("Cannot construct dummy symbol list for a non-finite alphabet");
        return new InfinitelyAmbiguousSymbolList((FiniteAlphabet)a);
    }
    
    /**
     * {@inheritDoc}
     */
    public BioEntry getRemoteBioEntry(CrossRef cr){
        Namespace ns = (Namespace)RichObjectFactory.getObject(SimpleNamespace.class, new Object[]{cr.getDbname()});
        try {
            // Build the query object
            String queryText = "from BioEntry where namespace = ? and accession = ? and version = ?";
            Object query = this.createQuery.invoke(this.session, new Object[]{queryText});
            // Set the parameters
            query = this.setParameter.invoke(query, new Object[]{new Integer(0), ns}); 
            query = this.setParameter.invoke(query, new Object[]{new Integer(1), cr.getAccession()}); 
            query = this.setParameter.invoke(query, new Object[]{new Integer(2), new Integer(cr.getVersion())}); 
            // Get the results
            Object result = this.uniqueResult.invoke(query, (Object[])null);
            // Return the found object, if found - null if not.
            return (BioEntry)result;
        } catch (Exception e) {
            // Throw the exception with our nice message
            throw new RuntimeException("Error while trying to locate remote cross reference "+cr,e);
        }
    }
}

