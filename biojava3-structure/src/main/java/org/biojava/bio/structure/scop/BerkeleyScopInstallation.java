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
}
