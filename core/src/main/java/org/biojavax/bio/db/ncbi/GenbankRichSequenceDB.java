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
package org.biojavax.bio.db.ncbi;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.net.MalformedURLException;
import java.net.Socket;
import java.net.URL;
import java.util.Iterator;
import java.util.Set;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.db.FetchURL;
import org.biojava.bio.seq.db.IllegalIDException;
import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojava.utils.ChangeVetoException;
import org.biojavax.Namespace;
import org.biojavax.RichObjectFactory;
import org.biojavax.bio.db.AbstractRichSequenceDB;
import org.biojavax.bio.db.HashRichSequenceDB;
import org.biojavax.bio.db.RichSequenceDB;
import org.biojavax.bio.db.RichSequenceDBLite;
import org.biojavax.bio.seq.RichSequence;
import org.biojavax.bio.seq.RichSequenceIterator;
import org.biojavax.bio.seq.io.RichSequenceBuilderFactory;


/**
 * This class contains functions accessing DNA sequences in Genbank format.
 * It adds methods to return RichSequences instead of plain Sequences.
 *
 * @author Lei Lai
 * @author Matthew Pocock
 * @author Laurent Jourdren
 * @author Shuvankar Mukherjee
 * @author Mark Schreiber
 * @author Richard Holland
 * @since 1.5
 */
public class GenbankRichSequenceDB extends AbstractRichSequenceDB implements RichSequenceDBLite {
    
    protected static final String urlBatchSequences = "http://www.ncbi.nlm.nih.gov:80/entrez/eutils/efetch.fcgi";

    private String email = "anonymous@biojava.org";
    private String tool = "biojavax";
    
    /**
     * The default constructor delegates to the parent class. The constructor refers
     * to RichObjectFactory.getDefaultNamespace() so make sure your factory is initialised
     * before calling this constructor.
     * Sets the default factory to THRESHOLD.
     */
    public GenbankRichSequenceDB() {
        super();
        this.setFactory(RichSequenceBuilderFactory.THRESHOLD); // threshold factory is efficient
        this.setNamespace(RichObjectFactory.getDefaultNamespace()); // default namespace
    }
    
    /**
     * Get the URL object for locating sequence object using eutils.
     * The default value of the return format of the sequence object is text.
     **/
    protected URL getAddress(String id) throws MalformedURLException {
        FetchURL seqURL = new FetchURL("Genbank", "text");
        String baseurl = seqURL.getbaseURL();
        String db = seqURL.getDB();
        String url = baseurl+db+"&id="+id+"&rettype=gb&tool="+getTool()+"&email="+getEmail();
        return new URL(url);
    }
    
    /**
     * Create the Http Post Request to fetch (in batch mode) a list of sequence
     * with Genbank.
     * @param url URL of the request
     * @param list List of sequence identifier
     * @return The Post request.
     */
    protected String makeBatchRequest(URL url, Set list) {
        StringBuffer params = new StringBuffer();
        params.append("db=nucleotide&rettype=gb&id=");
        
        for (Iterator i = list.iterator(); i.hasNext();) {
            String idSequence = (String)i.next();
            params.append(idSequence);
            if(i.hasNext()) params.append(",");
        }

        params.append("&email="+getEmail()+"&tool="+getTool());
        
        StringBuffer header = new StringBuffer();
        header.append("POST ");
        header.append(url.getPath());
        header.append(
                " HTTP/1.0\r\n"
                + "Connection: close\r\n"
                + "Accept: text/html, text/plain\r\n"
                + "Host: ");
        header.append(url.getHost());
        header.append(
                "\r\n"
                + "User-Agent: Biojava/GenbankSequenceDB\r\n"
                + "Content-Type: application/x-www-form-urlencoded\r\n"
                + "Content-Length: ");
        header.append(params.length());
        header.append("\r\n\r\n");
        
        StringBuffer request = new StringBuffer();
        request.append(header);
        request.append(params);
        
        return request.toString();
    }
    
    /**
     * Given the appropriate Genbank ID, return the matching RichSequence object.
     * @param id the Genbank ID to retrieve.
     * @return the matching RichSequence object, or null if not found.
     * @throws Exception if the sequence could not be retrieved for reasons other
     * than the identifier not being found.
     */
    public RichSequence getRichSequence(String id) throws BioException, IllegalIDException {
        try {
            URL queryURL = getAddress(id); //get URL based on ID
            
            SymbolTokenization rParser = DNATools.getDNA().getTokenization("token"); //get SymbolTokenization
            RichSequenceBuilderFactory seqFactory = this.getFactory();
            Namespace ns = this.getNamespace();
            
            DataInputStream in = new DataInputStream(queryURL.openStream());
            BufferedReader reader = new BufferedReader(new InputStreamReader(in));
            RichSequenceIterator seqI = RichSequence.IOTools.readGenbank(reader, rParser, seqFactory, ns);
            
            return seqI.nextRichSequence();
        } catch (MalformedURLException e) {
            throw new BioException("Failed to create Genbank URL",e);
        } catch (BioException e) {
            throw new BioException("Failed to read Genbank sequence",e);
        } catch (IOException e) {
            throw new BioException("IO failure whilst reading from Genbank",e);
        }
    }    
    
    /**     
     * Given the appropriate Genbank ID, return the matching RichSequence object. Additionally
     * define a new Namespace for the received RichSequence object.
     * @param id the Genbank ID to retrieve.
     * @param nsp the Namespace to define.
     * @return the matching RichSequence object, or null if not found.
     * @throws Exception if the sequence could not be retrieved for reasons other
     * than the identifier not being found.
     */
    public RichSequence getRichSequence(String id, Namespace nsp) throws BioException, IllegalIDException {
        try {
            URL queryURL = getAddress(id); //get URL based on ID
            
            SymbolTokenization rParser = DNATools.getDNA().getTokenization("token"); //get SymbolTokenization
            RichSequenceBuilderFactory seqFactory = this.getFactory();
            
            DataInputStream in = new DataInputStream(queryURL.openStream());
            BufferedReader reader = new BufferedReader(new InputStreamReader(in));
            RichSequenceIterator seqI = RichSequence.IOTools.readGenbank(reader, rParser, seqFactory, nsp);
            
            return seqI.nextRichSequence();
        } catch (MalformedURLException e) {
            throw new BioException("Failed to create Genbank URL",e);
        } catch (BioException e) {
            throw new BioException("Failed to read Genbank sequence",e);
        } catch (IOException e) {
            throw new BioException("IO failure whilst reading from Genbank",e);
        }
    }
    
    
    /**
     * Retrieve rich sequences from a Genbank
     *
     * @param list List of NCBI sequence number (GI), accession, accession.version,
     * fasta or seqid.
     * @return The rich database object (HashSequenceDB) with downloaded rich sequences.
     * You will need to cast the sequences you get from this database object into
     * RichSequence objects if you want to access their full features.
     */
    public RichSequenceDB getRichSequences(Set list) throws BioException, IllegalIDException {
        
        return getRichSequences(list, null);
    }
    
    /**
     * Retrieve rich sequences from a Genbank
     *
     * @param list List of NCBI sequence number (GI), accession, accession.version,
     * fasta or seqid.
     * @param database Where to store rich sequences. If database is null, use an
     * HashSequenceDB Object.
     * @return The database object with downloaded rich sequences.
     * You will need to cast the sequences you get from this database object into
     * RichSequence objects if you want to access their full features.
     */
    public RichSequenceDB getRichSequences(Set list, RichSequenceDB database) throws BioException, IllegalIDException {
        try {
            if (database == null) database = new HashRichSequenceDB();
            
            URL url = new URL(urlBatchSequences);
            int port = url.getPort();
            String hostname = url.getHost();
            
            //Open the connection and the streams
            Socket s = new Socket(hostname, port);
            
            InputStream sin = s.getInputStream();
            BufferedReader fromServer = new BufferedReader(new InputStreamReader(sin));
            OutputStream sout = s.getOutputStream();
            PrintWriter toServer = new PrintWriter(new OutputStreamWriter(sout));
            
            // Put the Post request to the server
            toServer.print(makeBatchRequest(url, list));
            toServer.flush();
            
            // Delete response headers
            boolean finEntete = false;
            for (String l = null; ((l = fromServer.readLine()) != null) && (!finEntete);) {
                if (l.equals("")) finEntete = true;
            }
            
            SymbolTokenization rParser = DNATools.getDNA().getTokenization("token"); //get SymbolTokenization
            RichSequenceBuilderFactory seqFactory = this.getFactory();
            Namespace ns = this.getNamespace();
            
            RichSequenceIterator seqI = RichSequence.IOTools.readGenbank(fromServer, rParser, seqFactory, ns);
            
            while (seqI.hasNext()) {
                try {
                    database.addSequence(seqI.nextRichSequence());
                } catch (ChangeVetoException ce) {
                    throw new BioException("Unexpectedly couldn't add to the supplied RichSequenceDB", ce);
                }
            }
            
            return database;
        } catch (MalformedURLException e) {
            throw new BioException("Failed to create Genbank URL",e);
        } catch (BioException e) {
            throw new BioException("Failed to read Genbank sequence",e);
        } catch (IOException e) {
            throw new BioException("IO failure whilst reading from Genbank",e);
        }
    }
    
    public String getName() {
        return "Genbank";
    }
    
    public Set ids() {
        throw new RuntimeException("Complete set of Genbank ids is unavailable.");
    }
    
    /**
     * Holds value of property factory.
     */
    private RichSequenceBuilderFactory factory;
    
    /**
     * Getter for property factory.
     * @return Value of property factory.
     */
    public RichSequenceBuilderFactory getFactory() {
        
        return this.factory;
    }
    
    /**
     * Setter for property factory.
     * @param factory New value of property factory.
     */
    public void setFactory(RichSequenceBuilderFactory factory) {
        
        this.factory = factory;
    }
    
    /**
     * Holds value of property namespace.
     */
    private Namespace namespace;
    
    /**
     * Getter for property namespace.
     * @return Value of property namespace.
     */
    public Namespace getNamespace() {
        
        return this.namespace;
    }
    
    /**
     * Setter for property namespace.
     * @param namespace New value of property namespace.
     */
    public void setNamespace(Namespace namespace) {
        
        this.namespace = namespace;
    }

    /** 
     * Set the tool identifier for Entrez. Defaults to 'biojavax'.
     * @param tool the new identifier.
     */
    public void setTool(String tool) {
        this.tool = tool;
    }

    /** 
     * Get the tool identifier for Entrez. Defaults to 'biojavax'.
     * @return the identifier.
     */
    public String getTool() {
        return this.tool;
    }

    /** 
     * Set the email for Entrez. Defaults to 'anonymous@biojava.org'.
     * @param email the new email.
     */
    public void setEmail(String email) {
        this.email = email;
    }

    /** 
     * Get the email for Entrez. Defaults to 'anonymous@biojava.org'.
     * @return the email.
     */
    public String getEmail() {
        return this.email;
    }
}
