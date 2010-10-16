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

package org.biojava.bio.search;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.biojava.bio.BioException;
import org.biojava.utils.TriState;
import org.biojava.utils.walker.Visitor;
import org.biojava.utils.walker.Walker;
import org.biojava.utils.walker.WalkerFactory;

public class FilteringContentHandler
    extends SearchContentAdapter
    implements BlastLikeSearchFilter.Node 
{
    // these store properties
    private Map searchProperty = new HashMap();
    private Map hitProperty = new HashMap();
    private Map subHitProperty = new HashMap();

    // these store the filter objects
    private List searchFilters = new ArrayList();
    private List hitFilters = new ArrayList();
    private List subHitFilters = new ArrayList();

    // these keep track of what this class has issued
    // to its downstream filters.
    private boolean emittedStartSearch;
    private boolean emittedStartHit;

    // these determine whether to skip the lower ranked levels
    private boolean skipHits =false;
    private boolean skipSubHits = false;

    // these keep track of whether the test for the outcome at a certain level has
    // been done.
    private boolean firstSubHit = true;
    private boolean firstHit = true;    // first hit in a search

    public Object getSearchProperty(Object key) { return searchProperty.get(key); }
    public Object getHitProperty(Object key) { return hitProperty.get(key); }
    public Object getSubHitProperty(Object key) { return subHitProperty.get(key); }

    private static final WalkerFactory walkerFactory = WalkerFactory.getInstance(BlastLikeSearchFilter.class);

    /**
     * Visitor class that parses the filter tree
     */
    public class FilterVisitor
        implements Visitor
    {
        public void bySearchProperty(BlastLikeSearchFilter.BySearchProperty sf)
        {
            // register this instance
            searchFilters.add(sf);
        }

        public void byHitProperty(BlastLikeSearchFilter.ByHitProperty sf)
        {
            hitFilters.add(sf);
        }

        public void bySubHitProperty(BlastLikeSearchFilter.BySubHitProperty sf)
        {
            subHitFilters.add(sf);
        }
    }

    private BlastLikeSearchFilter filter;
    private SearchContentHandler delegate;

    public FilteringContentHandler(BlastLikeSearchFilter filter, SearchContentHandler delegate)
        throws BioException
    {
        construct(filter);
        setSearchContentHandler(delegate);
    }

    public FilteringContentHandler(BlastLikeSearchFilter filter)
        throws BioException
    {
        construct(filter);
    }

    private void construct(BlastLikeSearchFilter filter)
        throws BioException
    {
        this.filter = filter;

        FilterVisitor visitor = new FilterVisitor();
        Walker walker = walkerFactory.getWalker(visitor);
        walker.walk(filter, visitor);
    }

    public void setSearchContentHandler(SearchContentHandler delegate)
    {
        this.delegate = delegate;
    }

    public void startHeader()
    {
        delegate.startHeader();
    }

    public void setDatabaseID(String id)
    {
        delegate.setDatabaseID(id);
    }

    public void endHeader()
    {
        delegate.endHeader();
    }

    public void startSearch()
    {
        emittedStartSearch = false;
        firstHit = true;
    }

    public void setQueryID(String queryID)
    {
        addSearchProperty(BlastLikeSearchFilter.KEY_QUERY_ID, queryID);
    }

    public void addSearchProperty(Object key, Object value)
    {
        searchProperty.put(key, value);
    }

    public void startHit()
    {
        emittedStartHit = false;
        firstSubHit = true;

        // determine if anything at the search level
        // will cause the rest of the hits/subhits
        // in this search to be skipped.
        // do only for first hit in a search
        if (firstHit) {
            for (Iterator searchFilterI = searchFilters.iterator(); searchFilterI.hasNext();) {
                BlastLikeSearchFilter sf = (BlastLikeSearchFilter) searchFilterI.next();
                sf.evaluate(this);
            }
    
            TriState filterStatus = filter.accept();
            if (filterStatus == TriState.FALSE) {
                skipHits = true;
                skipSubHits = true;
            }

            firstHit = false;
        }
    }

    public void addHitProperty(Object key, Object value)
    {
        if (skipHits) return;
        hitProperty.put(key, value);
    }

    public void startSubHit()
    {
        if (skipSubHits) return;


        // clear outcomes that depend on the properties
        // of this subhit
        for (Iterator subHitFilterI = hitFilters.iterator(); subHitFilterI.hasNext();) {
            BlastLikeSearchFilter sf = (BlastLikeSearchFilter) subHitFilterI.next();
            sf.reset();
        }

        // test filter at hit level here
        // determine if anything at the search level
        // will cause the rest of the hits/subhits
        // in this search to be skipped.
        // do only at first subhit in a search.
        if (firstSubHit) {
            for (Iterator hitFilterI = hitFilters.iterator(); hitFilterI.hasNext();) {
                BlastLikeSearchFilter sf = (BlastLikeSearchFilter) hitFilterI.next();
                sf.evaluate(this);
            }
    
            TriState filterStatus = filter.accept();

            if (filterStatus == TriState.FALSE) {
                skipSubHits = true;
            }

            firstSubHit = false;
        }
    }

    public void addSubHitProperty(Object key, Object value)
    {
        if (skipSubHits) return;
        subHitProperty.put(key, value);
    }

    public void endSubHit()
    {

        // test filter at subhit level here
        // this will be a decision on whether to emit events or not
        if (!skipSubHits) {
            for (Iterator subHitFilterI = subHitFilters.iterator(); subHitFilterI.hasNext();) {
                BlastLikeSearchFilter sf = (BlastLikeSearchFilter) subHitFilterI.next();
                sf.evaluate(this);
            }

            // handle emitting events to delegate
            if (filter.accept() == TriState.TRUE) {
                if (!emittedStartSearch) {
                    delegate.startSearch();
                    String queryId = (String) searchProperty.get(BlastLikeSearchFilter.KEY_QUERY_ID);
                    if (queryId != null) delegate.setQueryID(queryId);

                    // dump search properties to delegate omitting special keys
                    for (Iterator searchPropertyI = searchProperty.entrySet().iterator(); 
                        searchPropertyI.hasNext(); ) {
                        Map.Entry entry = (Map.Entry) searchPropertyI.next();
                        if (BlastLikeSearchFilter.KEY_QUERY_ID.equals(entry.getKey())) continue;

                        delegate.addSearchProperty(entry.getKey(), entry.getValue());
                    }

                    emittedStartSearch = true;
                }

                if (!emittedStartHit) {
                    delegate.startHit();
                    // dump hit properties to delegate.
                    for (Iterator hitPropertyI = hitProperty.entrySet().iterator();
                        hitPropertyI.hasNext(); ) {
                        Map.Entry entry = (Map.Entry) hitPropertyI.next();

                        delegate.addSearchProperty(entry.getKey(), entry.getValue());
                    }

                    emittedStartHit = true;
                }

                // dump subhit properties
                delegate.startSubHit();

                for (Iterator subHitPropertyI = subHitProperty.entrySet().iterator();
                    subHitPropertyI.hasNext(); ) {
                    Map.Entry entry = (Map.Entry) subHitPropertyI.next();

                    delegate.addSearchProperty(entry.getKey(), entry.getValue());
                }

                delegate.endSubHit();
            }
        }
    }

    public void endHit()
    {
        if (emittedStartHit) {
            delegate.endHit();
            emittedStartHit = false;
        }
    }

    public void endSearch()
    {
        if (emittedStartSearch) {
            delegate.endSearch();
            emittedStartSearch = false;
        }
    }

}

