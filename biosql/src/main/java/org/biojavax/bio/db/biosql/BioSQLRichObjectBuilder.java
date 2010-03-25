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
import java.lang.reflect.Constructor;
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.biojavax.DocRefAuthor;
import org.biojavax.RichObjectBuilder;
import org.biojavax.SimpleCrossRef;
import org.biojavax.SimpleDocRef;
import org.biojavax.SimpleNamespace;
import org.biojavax.bio.taxa.SimpleNCBITaxon;
import org.biojavax.ontology.SimpleComparableOntology;


/**
 * Takes requests for RichObjects and sees if it can load them from a Hibernate
 * database. If it can, it returns the loaded objects. Else, it creates them
 * and persists them, then returns them. Doesn't retain a memory map so runs
 * a lot of Hibernate queries if used frequently, but on the plus side this
 * makes it memory-efficient.
 * @author Richard Holland
 * @author David Scott
 * @author Deepak Sheoran
 * @since 1.5
 */
public class BioSQLRichObjectBuilder implements RichObjectBuilder {
    
    private Object session;
    private Method createQuery;
    private Method setParameter;
    private Method uniqueResult;
    private Method persist;
    
    /**
     * Creates a new instance of SimpleRichObjectBuilder. The session parameter
     * is a Hibernate Session object and must not be null. It is this session
     * that database objects will be retrieved from/persisted to.
     * @see <a href="http://www.hibernate.org/hib_docs/v3/api/org/hibernate/Session.html"> org.hibernate.Session</a>
     */
    public BioSQLRichObjectBuilder(Object session) {
        try {
            // Lazy load the Session class from Hibernate.
            Class hibernateSession = session.getClass();
            Class realHibernateSession = Class.forName("org.hibernate.Session");
            // Test to see if our parameter is really a Session
            if (!realHibernateSession.isAssignableFrom(hibernateSession))
                throw new IllegalArgumentException("Parameter must be a org.hibernate.Session object");
            this.session = session;
            // Lookup the createQuery and persist methods
            this.createQuery = hibernateSession.getMethod("createQuery", new Class[]{String.class});
            this.persist = hibernateSession.getMethod("persist", new Class[]{String.class,Object.class});
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
     * Attempts to look up the details of the object in the database. If it
     * finds them it loads the object and returns it. Else, it persists the
     * wrapper object it made to do the search with and returns that.
     */
    public Object buildObject(Class clazz, List paramsList) {
        // convert the params list to remove nulls as we can't process those.
        List ourParamsList = new ArrayList(paramsList);
        for (Iterator i = ourParamsList.iterator(); i.hasNext(); ) 
        	if (i.next()==null) i.remove();
        // Create the Hibernate query to look it up with
        String queryText;
        String queryType;
        if (SimpleNamespace.class.isAssignableFrom(clazz)) {
            queryText = "from Namespace as ns where ns.name = ?";
            queryType = "Namespace";
        } else if (SimpleComparableOntology.class.isAssignableFrom(clazz)) {
            queryText = "from Ontology as o where o.name = ?";
            queryType = "Ontology";
        } else if (SimpleNCBITaxon.class.isAssignableFrom(clazz)) {
            queryText = "from Taxon as o where o.NCBITaxID = ?";
            queryType = "Taxon";
        } else if (SimpleCrossRef.class.isAssignableFrom(clazz)) {
            queryText = "from CrossRef as cr where cr.dbname = ? and cr.accession = ? and cr.version = ?";
            queryType = "CrossRef";
        } else if (SimpleDocRef.class.isAssignableFrom(clazz)) {
        	// First check if record exists with pubmed or medline id, if provided
        	if (ourParamsList.size() > 3) {
        		List crossRefParams = new ArrayList();
        		String crossRefQueryText = "from DocRef as dr where dr.crossref = (select id from CrossRef as cref where cref.dbname = ? and cref.accession = ? and cref.version = ?)";
        		crossRefParams.add(ourParamsList.get(ourParamsList.size() - 3)); // get dbname
        		crossRefParams.add(ourParamsList.get(ourParamsList.size() - 2)); // get accession
        		crossRefParams.add(ourParamsList.get(ourParamsList.size() - 1)); // get version
                // Build the query object
                Object query = this.createQuery.invoke(this.session, new Object[]{crossRefQueryText});
                // Set the parameters
                for (int i = 0; i < crossRefParams.size(); i++) {
                    query = this.setParameter.invoke(query, new Object[]{new Integer(i), crossRefParams.get(i)});
                }
                // Get the results
                Object result = this.uniqueResult.invoke(query, (Object[])null);
                // Return the found object, if found
                if (result!=null) return result;
                // If we get here, resort to looking up by author and title and remove all references
                // to pubmed/medline from params list.
                ourParamsList.remove(ourParamsList.size() - 3);
        		ourParamsList.remove(ourParamsList.size() - 2);
        		ourParamsList.remove(ourParamsList.size() - 1);
        	}
        	queryType = "DocRef";
        	// convert List constructor to String representation for query
        	ourParamsList.set(0, DocRefAuthor.Tools.generateAuthorString((List)ourParamsList.get(0), true));
        	if (ourParamsList.size()<3) {
        		queryText = "from DocRef as cr where cr.authors = ? and cr.location = ? and cr.title is null";
        	} else {
        		queryText = "from DocRef as cr where cr.authors = ? and cr.location = ? and cr.title = ?";
        	}        
        } else throw new IllegalArgumentException("Don't know how to handle objects of type "+clazz);
        // Run the query.
        try {
            // Build the query object
            Object query = this.createQuery.invoke(this.session, new Object[]{queryText});
            // Set the parameters
            for (int i = 0; i < ourParamsList.size(); i++) {
                query = this.setParameter.invoke(query, new Object[]{new Integer(i), ourParamsList.get(i)});
            }
            // Get the results
            Object result = this.uniqueResult.invoke(query, (Object[])null);
            // Return the found object, if found
            if (result!=null) return result;
            // Create, persist and return the new object otherwise
            else {
                // Load the class
                Class[] types = new Class[ourParamsList.size()];
                // Find its constructor with given params
                for (int i = 0; i < ourParamsList.size(); i++) {
                    if (ourParamsList.get(i) instanceof Set) types[i] = Set.class;
                    else if (ourParamsList.get(i) instanceof Map) types[i] = Map.class;
                    else if (ourParamsList.get(i) instanceof List) types[i] = List.class;
                    else types[i] = ourParamsList.get(i).getClass();
                }
                Constructor c = clazz.getConstructor(types);
                // Instantiate it with the parameters
                Object o = c.newInstance(ourParamsList.toArray());
                // Persist and return it.
                this.persist.invoke(this.session, new Object[]{queryType,o});
                return o;
            }
        } catch (Exception e) {
            // Write a useful message explaining what we were trying to do. It will
            // be in the form "class(param,param...)".
            StringBuffer paramsstuff = new StringBuffer();
            paramsstuff.append(clazz);
            paramsstuff.append("(");
            for (int i = 0; i < ourParamsList.size(); i++) {
                if (i>0) paramsstuff.append(",");
            	paramsstuff.append(ourParamsList.get(i).getClass());
            }
            paramsstuff.append(")");
            // Throw the exception with our nice message
            throw new RuntimeException("Error while trying to call new "+paramsstuff,e);
        }
    }
    
}
