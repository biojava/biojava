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

package org.biojava.bio.program.sax;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.HashSet;
import java.util.Set;
import java.util.StringTokenizer;

import org.biojava.bio.BioException;
import org.biojava.bio.search.SearchContentHandler;
import org.biojava.utils.ParserException;

/**
 * <p><code>FastaSearchParser</code> objects provide Fasta search
 * parsing functionality for the '-m 10' output format (see the Fasta
 * documentation). Data are passed to a SearchContentHandler which
 * coordinates its interpretation and creation of objects representing
 * the result.</p>
 *
 * <p>If the search output contains no hits only the header data will
 * be sent to the handler before the dataset ends. The handler is
 * responsible for dealing with this state.</p>
 *
 * <p>This class was originally used outside of the BioJava SAX
 * framework, but is now only used to provide functionality for
 * <code>FastaSearchSAXParser</code>.</p>
 *
 * @author Keith James
 * @author Greg Cox
 * @since 1.1
 */
class FastaSearchParser
{
    private static final int    NODATA = 0;
    private static final int  INHEADER = 1;
    private static final int     INHIT = 2;
    private static final int   INQUERY = 3;
    private static final int INSUBJECT = 4;
    private static final int   INALIGN = 5;

    // Valid line identifiers for result annotation
    private static HashSet resultAnnoTokens =
        (HashSet) fillSet(new String [] { "mp_name",   "mp_ver",
                                          "mp_argv",   "mp_extrap",
                                          "mp_stats",  "mp_KS",
                                          "pg_name",   "pg_ver",
                                          "pg_optcut", "pg_cgap" },
                          new HashSet());

    // Valid line identifiers for search parameters
    private static HashSet resultSearchParmTokens =
        (HashSet) fillSet(new String [] { "pg_matrix", "pg_ktup",
                                          "pg_gap-pen" },
                          new HashSet());

    // Valid line identifiers for hit annotation
    private static HashSet hitAnnoTokens =
        (HashSet) fillSet(new String [] { "fa_frame",   "fa_initn",
                                          "fa_init1",   "fa_opt",
                                          "fa_bits",    "sw_score",
                                          "sw_ident",   "sw_gident",
                                          "sw_overlap", "fa_ident",
                                          "fa_gident",  "fa_overlap",
                                          "fa_score" },
                          new HashSet());

    // Valid line identifiers for hit data
    private static HashSet hitDataTokens =
        (HashSet) fillSet(new String [] { "fa_expect", "fa_z-score" },
                          new HashSet());

    private int     searchStatus = NODATA;
    private boolean searchParsed = false;

    private SearchContentHandler handler;
    private String               line;
    private int                  lineNumber;

    StringBuffer   querySeqTokens = new StringBuffer(1024);
    StringBuffer subjectSeqTokens = new StringBuffer(1024);
    StringBuffer      matchTokens = new StringBuffer(1024);

    /**
     * The <code>parseSearch</code> method performs the core parsing
     * operations. It parses one result from the stream before
     * returning. The handler is informed whether or not there are
     * further searches in the stream.
     *
     * @param reader a <code>BufferedReader</code> to read from.
     * @param handler a <code>SearchContentHandler</code> to notify
     * of events.
     *
     * @exception IOException if the BufferedReader fails.
     * @exception BioException if the parser (via the registered
     * SearchContentHandler) fails to resolve a query sequence and
     * target database.
     * @exception ParserException if the parser fails to parse a
     * line.
     */
    public void parseSearch(BufferedReader       reader,
                            SearchContentHandler handler)
        throws IOException, BioException, ParserException
    {
        lineNumber = 0;
        this.handler = handler;

    LINE:
        while ((line = reader.readLine()) != null)
        {
            lineNumber++;
            // This token indicates the end of the formatted search
            // data. Some outputs don't have any alignment consensus
            // tokens, so we need to check here as well as INALIGN
            if (line.startsWith(">>><<<"))
            {
                // If we got here while INHEADER then the search had
                // no hits
                if (searchStatus == INHEADER)
                {
                    // Inform of end of header and search
                    handler.endHeader();
                    handler.endSearch();

                    searchParsed = true;
                    searchStatus = NODATA;
                    continue LINE;
                }

                // Pass final data to handler
                handler.addSubHitProperty("querySeqTokens",   querySeqTokens.substring(0));
                handler.addSubHitProperty("subjectSeqTokens", subjectSeqTokens.substring(0));
                handler.addSubHitProperty("matchTokens",      matchTokens.substring(0));

                handler.endSubHit();
                handler.endSearch();

                searchParsed = true;
                searchStatus = NODATA;

                // Continue getting lines to look for start of another
                // search. This allows setMoreSearches(boolean flag)
                // to be called on the handler. When true it will know
                // to expect more.
                continue LINE;
            }

        STATUS:
            switch (searchStatus)
            {
                case NODATA:
                    // This token marks the line describing the query
                    // sequence file and database searched. It is
                    // followed by header lines containing data about
                    // the search
                    if (line.startsWith(">>>"))
                    {
                        searchStatus = INHEADER;

                        handler.setQueryID(parseQueryID(line));
                        handler.setDatabaseID(parseDatabaseID(line));

                        handler.startSearch();
                        handler.startHeader();

                        // If we already saw an end of search token
                        // then this is the start of another
                        // dataset. We break from the loop and return
                        // that the stream is not empty
                        if (searchParsed)
                        {
                            searchParsed = false;
                            handler.setMoreSearches(true);

                            // We have parsed one result, so return
                            return;
                        }
                        else
                            // We break out and setMoreSearches(false)
                            // is called below
                            break STATUS;
                    }
                    else
                        // Continue getting lines to look for start of another
                        // search. This allows setMoreSearches(boolean flag)
                        // to be called on the handler. When true it will know
                        // to expect more.
                        continue LINE;

                case INHEADER:
                    // This token marks the line describing a hit. It
                    // is followed by header lines containing data
                    // about the hit
                    if (line.startsWith(">>"))
                    {
                        searchStatus = INHIT;

                        handler.endHeader();
                        handler.startHit();

                        querySeqTokens.setLength(0);
                        subjectSeqTokens.setLength(0);
                        matchTokens.setLength(0);

                        handler.addHitProperty("id",   parseID(line));
                        handler.addHitProperty("desc", parseDesc(line));
                    }
                    else
                    {
                        if (! parseHeaderLine(line, resultAnnoTokens))
                            if (! parseHeaderLine(line, resultSearchParmTokens))
                                throw new ParserException("Fasta parser failed to recognise line type",
                                                          null,
                                                          lineNumber,
                                                          line);
                    }
                    break STATUS;

                case INHIT:
                    // This token marks the line describing the query
                    // sequence.
                    if (line.startsWith(">"))
                    {
                        searchStatus = INQUERY;
                        handler.endHit();
                        handler.startSubHit();
                    }
                    else
                    {
                        if (! parseHitLine(line, hitAnnoTokens))
                            if (! parseHitLine(line, hitDataTokens))
                                throw new ParserException("Fasta parser failed to recognise line type",
                                                          null,
                                                          lineNumber,
                                                          line);
                    }
                    break STATUS;

                case INQUERY:
                    // This token marks the line describing the
                    // subject sequence.
                    if (line.startsWith(">"))
                    {
                        searchStatus = INSUBJECT;
                    }
                    else
                    {
                        parseQuerySequence(line);
                    }
                    break STATUS;

                case INSUBJECT:
                    // This token marks the start of lines containing
                    // the consensus symbols from the Fasta alignment,
                    // which we ignore
                    if (line.startsWith("; al_cons:"))
                    {
                        searchStatus = INALIGN;
                    }
                    else if (line.startsWith(">>"))
                    {
                        searchStatus = INHIT;

                        // Pass data to handler
                        handler.addSubHitProperty("querySeqTokens",   querySeqTokens.substring(0));
                        handler.addSubHitProperty("subjectSeqTokens", subjectSeqTokens.substring(0));
                        handler.addSubHitProperty("matchTokens",      matchTokens.substring(0));

                        handler.endSubHit();
                        handler.startHit();

                        querySeqTokens.setLength(0);
                        subjectSeqTokens.setLength(0);
                        matchTokens.setLength(0);

                        handler.addHitProperty("id",   parseID(line));
                        handler.addHitProperty("desc", parseDesc(line));
                    }
                    else
                    {
                        parseSubjectSequence(line);
                    }
                    break STATUS;

                case INALIGN:
                    if (line.startsWith(">>"))
                    {
                        searchStatus = INHIT;

                        // Pass data to handler
                        handler.addSubHitProperty("querySeqTokens",   querySeqTokens.substring(0));
                        handler.addSubHitProperty("subjectSeqTokens", subjectSeqTokens.substring(0));
                        handler.addSubHitProperty("matchTokens",      matchTokens.substring(0));

                        handler.endSubHit();
                        handler.startHit();

                        querySeqTokens.setLength(0);
                        subjectSeqTokens.setLength(0);
                        matchTokens.setLength(0);

                        handler.addHitProperty("id",   parseID(line));
                        handler.addHitProperty("desc", parseDesc(line));
                    }
                    else if (line.startsWith(">>><<<"))
                    {
                        searchStatus = NODATA;

                        // Pass final data to handler
                        handler.addSubHitProperty("querySeqTokens",   querySeqTokens.substring(0));
                        handler.addSubHitProperty("subjectSeqTokens", subjectSeqTokens.substring(0));
                        handler.addSubHitProperty("matchTokens",      matchTokens.substring(0));

                        handler.endSubHit();
                        handler.endSearch();

                        searchParsed = true;
                        handler.setMoreSearches(true);

                        continue LINE;
                    }
                    else
                    {
                        matchTokens.append(line);
                    }
                    break STATUS;

                default:
                    break STATUS;
            } // end switch
        } // end while

        // This is false if we reach here and return
        handler.setMoreSearches(false);
    }

    /**
     * The <code>fillSet</code> method populates a <code>Set</code>
     * with the elements of an array.
     *
     * @param tokenArray a <code>String []</code> array.
     * @param set a <code>Set</code> to fill.
     *
     * @return a <code>Set</code>.
     */
    private static Set fillSet(String [] tokenArray, Set set)
    {
        for (int i = 0; i < tokenArray.length; i++)
            set.add(tokenArray[i]);

        return set;
    }

    /**
     * The <code>parseID</code> method parses sequence IDs from
     * lines starting with '>' and '>>'.
     *
     * @param line a <code>String</code> to be parsed.
     *
     * @return a <code>String</code> containing the ID.
     *
     * @exception ParserException if an error occurs.
     */
    private String parseID(String line)
        throws ParserException
    {
        String trimmed = line.trim();
        int firstSpace = trimmed.indexOf(' ');

        // For Hit header lines (always start with >>)
        if (trimmed.startsWith(">>"))
        {
            if (trimmed.length() == 2)
                throw new ParserException("Fasta parser encountered a sequence with no Id",
                                          null,
                                          lineNumber,
                                          line);

            if (firstSpace == -1)
                return trimmed.substring(2);
            else
                return trimmed.substring(2, firstSpace);
        }
        // For SubHit header lines (always start with >)
        else
        {
            if (trimmed.length() == 1)
                throw new ParserException("Fasta parser encountered a sequence with no Id",
                                          null,
                                          lineNumber,
                                          line);

            if (firstSpace == -1)
                return trimmed.substring(1);
            else
                return trimmed.substring(1, firstSpace);
        }
    }

    /**
     * The <code>parseDesc</code> method parses the sequence
     * description from subject header lines.
     *
     * @param line a <code>String</code> to be parsed.
     *
     * @return a <code>String</code> containing the description.
     */
    private String parseDesc(String line)
    {
        String trimmed = line.trim();
        int firstSpace = trimmed.indexOf(' ');

        if (firstSpace == -1)
            return "No description";

        return trimmed.substring(firstSpace + 1);
    }

    /**
     * Creates a new <code>parseQueryID</code> instance.
     *
     * @param line a <code>String</code> to be parsed.
     *
     * @return a <code>String</code>.
     *
     * @exception ParserException if an error occurs
     */
    private String parseQueryID(String line)
        throws ParserException
    {
        int firstComma = line.indexOf(",");

        if (firstComma == -1)
            throw new ParserException("Fasta parser failed to parse a query ID",
                                      null,
                                      lineNumber,
                                      line);

        // return string between >>> and ,
        return line.substring(3, firstComma).trim();
    }

    /**
     * The <code>parseDatabaseID</code> method parses a database
     * filename from the relevant output line.
     *
     * @param line a <code>String</code> to be parsed.
     *
     * @return a <code>String</code>.
     *
     * @exception ParserException if an error occurs.
     */
    private String parseDatabaseID(String line)
        throws ParserException
    {
        StringTokenizer st = new StringTokenizer(line);
        String id = null;

        int count = st.countTokens();

        for (int i = 0; i < count - 1; i++)
        {
            id = st.nextToken();
        }

        if (id == null)
            throw new ParserException("Fasta parser failed to parse a database ID",
                                      null,
                                      lineNumber,
                                      line);
        return id;
    }

    private boolean parseHeaderLine(String line,
                                    Set    tokenSet)
        throws ParserException
    {
        String [] data = parseLine(line, tokenSet);

        if (data.length > 0)
        {
            handler.addSearchProperty(data[0], data[1]);
            return true;
        }
        else
        {
            return false;
        }
    }

    private boolean parseHitLine(String line,
                                 Set    tokenSet)
        throws ParserException
    {
        String [] data = parseLine(line, tokenSet);

        if (data.length > 0)
        {
            handler.addHitProperty(data[0], data[1]);
            return true;
        }

        return false;
    }

    private String [] parseLine(String line,
                                Set    tokenSet)
        throws ParserException
    {
        int idTokenStart = line.indexOf(";");
        int   idTokenEnd = line.indexOf(":");

        String   idToken = line.substring(idTokenStart + 1, idTokenEnd);
        idToken          = idToken.trim();

        String   idValue = line.substring(idTokenEnd + 1);
        idValue          = idValue.trim();

        if (tokenSet.contains(idToken))
            return new String [] { idToken, idValue };
        else
            return new String [0];
    }

    private void parseQuerySequence(String line)
    {
        String [] data = parseSequence(line);

        if (data.length > 0)
        {
            // We have a key/value pair
            handler.addSubHitProperty("query" + data[0], data[1]);
        }
        else
        {
            // We have a line of sequence tokens
            querySeqTokens.append(line);
        }
    }

    private void parseSubjectSequence(String line)
    {
        String [] data = parseSequence(line);

        if (data.length > 0)
        {
            // We have a key/value pair
            handler.addSubHitProperty("subject" + data[0], data[1]);
        }
        else
        {
            // We have a line of sequence tokens
            subjectSeqTokens.append(line);
        }
    }

    private String [] parseSequence(String line)
    {
        if (line.startsWith(";"))
        {
            // Check the sequence type given by the report
            if (line.equals("; sq_type: p"))
            {
                return new String [] { "_sq_type", "protein"};
            }
            else if (line.equals("; sq_type: D"))
            {
                return new String [] { "_sq_type", "dna"};
            }

            // Record the coordinates and offset of the alignment
            if (line.startsWith("; al_start:"))
            {
                return new String [] { "_al_start", parseCoord(line) };
            }
            else if (line.startsWith("; al_stop:"))
            {
                 return new String [] { "_al_stop", parseCoord(line) };
            }
            else if (line.startsWith("; al_display_start:"))
            {
                return new String [] {"_al_display_start", parseCoord(line) };
            }
            else if (line.startsWith("; sq_len:"))
            {
                return new String [] { "_sq_len", parseCoord(line) };
            }
            else if (line.startsWith("; sq_offset:"))
            {
                return new String [] { "_sq_offset", parseCoord(line) };
            }
        }

        return new String [0];
    }

    /**
     * The <code>parseCoord</code> method extracts integer coordinates
     * from Fasta output lines.
     *
     * @param line a <code>String</code> to parse.
     *
     * @return a <code>String</code> coordinate.
     */
    private String parseCoord(String line)
    {
        int sepIndex = line.lastIndexOf(":");
        return line.substring(sepIndex + 1).trim();
    }
}
