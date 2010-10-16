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

package org.biojava.bio.seq.db.biofetch;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.HttpURLConnection;
import java.net.URL;
import java.util.StringTokenizer;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.db.IllegalIDException;
import org.biojava.bio.seq.db.SequenceDBLite;
import org.biojava.bio.seq.io.SeqIOTools;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.Unchangeable;

/**
 * Simple SequenceDB implementation backed by a BioFetch (HTTP)
 * server.
 *
 * @author Thomas Down
 * @author Matthew Pocock
 * @author Greg Cox
 * @author Keith James
 * @since 1.3
 */
public class BioFetchSequenceDB
  extends
    Unchangeable
  implements
    SequenceDBLite {
    private final String location;
    private final String dbName;

    /**
     * Construct a BioFetchSequenceDB which connects to the specified
     * BioFetch server.
     *
     * @param location The base URL of the server.
     * @param dbName The database name to use.
     */
    public BioFetchSequenceDB(String location,
                                          String dbName) {
        this.location = location;
        this.dbName = dbName;
    }

    public String getName() {
        return dbName;
    }

    public void addSequence(Sequence seq)  throws ChangeVetoException {
        throw new ChangeVetoException("Failed to add sequence."
                                      + " Sequences may not be added"
                                      + " to a biofetch database");
    }

    public void removeSequence(String id) throws ChangeVetoException {
        throw new ChangeVetoException("Failed to add sequence."
                                      + " Sequences may not be removed"
                                      + " from a biofetch database");
    }

    public Sequence getSequence(String id)
        throws BioException, IllegalIDException
    {
        String format = "";

        if (dbName.equals("embl"))
            format = "embl";
        if (dbName.equals("genbank"))
            format = "genbank";
        else if (dbName.equals("swiss"))
            format = "swissprot";
        else if (dbName.equals("refseq"))
            throw new BioException("Sequence database "
                                   + dbName
                                   + " is not supported");

        StringBuffer uri = new StringBuffer(location);
        uri.append('?');
        uri.append("style=raw;");
        uri.append("format=");
        uri.append(format);
        uri.append(";db=");
        uri.append(dbName);
        uri.append(";id=");
        uri.append(id);

        try {
            HttpURLConnection huc =
                (HttpURLConnection) new URL(uri.substring(0)).openConnection();
            huc.connect();
            BufferedReader data =
                new BufferedReader(new InputStreamReader(huc.getInputStream()));
            data.mark(1000);

            String firstLine = data.readLine();
            if (firstLine.startsWith("Content-")) {
                data.readLine();
                firstLine = data.readLine();
            }

            StringTokenizer toke = new StringTokenizer(firstLine);
            String first = toke.nextToken();
            if ("ERROR".equals(first)) {
                int errorCode = Integer.parseInt(toke.nextToken());
                if (errorCode == 4) {
                    throw new IllegalIDException("No such ID "
                                                 + id
                                                 + " in database "
                                                 + getName());
                } else {
                    throw new BioException("Error fetching from BioFetch:"
                                           + firstLine);
                }
            }

            data.reset();

            SequenceIterator si = SeqIOTools.readEmbl(data);

            if (dbName.equals("embl"))
                si = SeqIOTools.readEmbl(data);
            else if (dbName.equals("genbank"))
                si = SeqIOTools.readGenbank(data);
            else if (dbName.equals("swiss"))
                si = SeqIOTools.readSwissprot(data);

            return si.nextSequence();

        } catch (IOException ex) {
            throw new BioException("Error reading data from BioFetch",ex);
        }
    }
}
