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
 * Created on 21.05.2004
 * @author Andreas Prlic
 *
 *
 */

package org.biojava.bio.structure.io;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.net.Socket;

import org.biojava.bio.structure.Structure;

/** reads a PDB file from a local SRS installation using getz Actually
 * is the same as PDBFileReader, but instead of reading from a file stream, reads from a
 * buffered stream.
 * if no matching PDB code found, returns null
 *
 * @author Andreas Prlic
 *
*/

public class PDBSRSReader implements StructureIO {  
    

    private  BufferedReader getBufferedReader(String pdbId) 
	throws IOException
    {
	// getz -view PDBdata '[pdb:5pti] 
	Socket          client  = null;
	DataInputStream input   = null;
	PrintStream     output  = null;
	String          machine = ""  ;
	int             port    = 0   ;

	String message = "please set System properties PFETCH_machine and PFETCH_port !" ;

	try {
	    // shouldbe sent as argument ...
	    machine     = System.getProperty("PFETCH_host");
	    String p    = System.getProperty("PFETCH_port");
	    port        = Integer.parseInt(p);

	} catch ( NullPointerException e) {
	    System.err.println(message);
	    e.printStackTrace();
	    throw new IOException() ;
	} catch (IllegalArgumentException  e) {
	    System.err.println(message);
	    e.printStackTrace(); 
	    throw new IOException() ;
	}

	if (port  == 0 ) {
	    throw new IOException(message);
	}
	if ( (machine.equals(""))) {	   
	    throw new IOException(message); 
	}
	System.out.println("contacting: " + machine + " " + port);
	//Process proc = Runtime.getRuntime().exec(GETZSTRING+argument);
	client = new Socket(machine , port);
	client.setSoTimeout(10000) ; // 10 seconds
	System.out.println("socket o.k.");
	input  = new DataInputStream(client.getInputStream());
	BufferedReader buf = new BufferedReader (new InputStreamReader (input));
	
	System.out.println("sending: --pdb " + pdbId.toLowerCase());
	output = new PrintStream(client.getOutputStream());	  
	output.println("--pdb "+ pdbId.toLowerCase());
	output.flush();
	
	// check if return is O.K.
	buf.mark(100);
	String line = buf.readLine();

	buf.reset();
	if ( line.equals("no match")) {
	    System.out.println("first line: " + line );
	    throw new IOException("no pdb with code "+pdbId.toLowerCase() +" found");	    
	}
	
	return buf ;


    }


     /** load a structure from from SRS installation using wgetz
      * returns null if no structure found
     */
    public  Structure getStructureById(String pdbId) 
	throws IOException
    {
	
	BufferedReader buf ;
	//inStream = getInputStream();
	try {
	    buf = getBufferedReader(pdbId) ;
	}
	catch (IOException e) {
	    // no pdb code found that is suitable ;
	    return null ;
	}
	/*String line = buf.readLine ();	
	while (line != null) {
	    System.out.println (line);
	    line = buf.readLine ();

	}
	return null ;
	*/
	Structure s = null ;
	try{	    
	    //System.out.println("Starting to parse PDB file " + getTimeStamp());
	    PDBFileParser pdbpars = new PDBFileParser();
	    s = pdbpars.parsePDBFile(buf) ;
	    //System.out.println("Done parsing PDB file " + getTimeStamp());
	} catch(Exception ex){
	    ex.printStackTrace();
	}

	return s ;	
    }

}
