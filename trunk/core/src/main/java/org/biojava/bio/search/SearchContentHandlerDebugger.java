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

import java.io.PrintStream;

/**
 * This class prints to a PrintStream
 * calls to the SearchContentHandler interface
 * in human readable form.  Use to debug parser/adaptor
 * classes that output to the SearchContentHandler interface.
 * @author David Huen
 */
public class SearchContentHandlerDebugger
    implements SearchContentHandler
{
    private String [] margin = {"", "\t", "\t\t", "\t\t\t"};
    private int nesting = 0;
    private boolean moreSearches = false;

    /**
     * Create an instance that dumps to System.out.
     */
    public SearchContentHandlerDebugger()
    {
    }

    /**
     * @param pStream Stream to dump output to.
     */
    public SearchContentHandlerDebugger(PrintStream pStream)
    {
    }

    public void addHitProperty(Object key, Object value)
    {
        System.out.println(margin[nesting] + key.toString() + "->" + value.toString());
    }

    public void addSubHitProperty(Object key, Object value)
    {
        System.out.println(margin[nesting] + key.toString() + "->" + value.toString());
    }

    public void addSearchProperty(Object key, Object value)
    {
        System.out.println(margin[nesting] + key.toString() + "->" + value.toString());
    }

    public void setMoreSearches(boolean value) { moreSearches = value; }
    public boolean getMoreSearches() { return moreSearches; }

    public void endHeader() { nesting--; }
    public void endHit() { nesting--; }
    public void endSearch() { nesting--; }
    public void endSubHit() { nesting--; }

    public void setDatabaseID(String databaseID)
    {
        System.out.println(margin[nesting] + "setDatabaseID: " + databaseID);
    }

    public void setQueryID(String queryID)
    {
        System.out.println(margin[nesting] + "setQueryID: " + queryID);
    }

    public void startHeader()
    {
        System.out.println(margin[nesting] + "Start header");
        nesting++;
    }

    public void startHit()
    {
        System.out.println(margin[nesting] + "Start hit");
        nesting++;
    }

    public void startSearch()
    {
        System.out.println(margin[nesting] + "Start search");
        nesting++;
    }

    public void startSubHit()
    {
        System.out.println(margin[nesting] + "Start subhit");
        nesting++;
    }
}

