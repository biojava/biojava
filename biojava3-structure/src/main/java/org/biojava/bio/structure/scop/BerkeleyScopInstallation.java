/**
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
 * Created on Oct 12, 2011
 * Created by Andreas Prlic
 *
 * @since 3.0.2
 */
package org.biojava.bio.structure.scop;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.net.URL;

public class BerkeleyScopInstallation extends ScopInstallation implements
ScopDatabase {


	String defaultBerkeleyDownloadURL = "http://scop.berkeley.edu/downloads/parse/";
	String defaultBerkeleyScopVersion="1.75A";
	public static final String claFileName = "dir.cla.scop.";
	public static final String desFileName = "dir.des.scop.";
	public static final String hieFileName = "dir.hie.scop.";
	public static final String comFileName = "dir.com.scop.";
	
	
	public BerkeleyScopInstallation(){
		super();

		setScopVersion(defaultBerkeleyScopVersion);
		setScopDownloadURL(defaultBerkeleyDownloadURL);
	}

	
	protected void downloadClaFile() throws FileNotFoundException, IOException{
		String remoteFilename = claFileName + scopVersion + ".txt";
		URL url = new URL(scopDownloadURL + remoteFilename);

		String localFileName = getClaFilename();
		File localFile = new File(localFileName);

		downloadFileFromRemote(url, localFile);

	}


	protected void downloadDesFile() throws FileNotFoundException, IOException{
		String remoteFilename = desFileName + scopVersion + ".txt";
		URL url = new URL(scopDownloadURL + remoteFilename);

		String localFileName = getDesFilename();
		File localFile = new File(localFileName);

		downloadFileFromRemote(url, localFile);

	}

	protected void downloadHieFile() throws FileNotFoundException, IOException{
		String remoteFilename = hieFileName + scopVersion + ".txt";
		URL url = new URL(scopDownloadURL + remoteFilename);

		String localFileName = getHieFilename();
		File localFile = new File(localFileName);

		downloadFileFromRemote(url, localFile);

	}
	
	protected void downloadComFile() throws FileNotFoundException, IOException{
		String remoteFilename = comFileName + scopVersion + ".txt";
		URL url = new URL(scopDownloadURL + remoteFilename);

		String localFileName = getComFilename();
		File localFile = new File(localFileName);

		downloadFileFromRemote(url, localFile);
	}

}
