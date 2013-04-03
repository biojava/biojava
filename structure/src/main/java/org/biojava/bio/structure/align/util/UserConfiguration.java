/*
 *                    PDB web development code
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
 *
 * Created on Jul 8, 2009
 * Created by ap3
 *
 */

package org.biojava.bio.structure.align.util;

import java.io.IOException;
import java.io.PrintWriter;

import org.biojava.bio.structure.align.ce.StartupParameters;
import org.biojava.bio.structure.io.PDBFileReader;
import org.biojava.utils.xml.PrettyXMLWriter;
import org.biojava.utils.xml.XMLWriter;

/** A container to persist config to the file system
 * 
 * @author Andreas Prlic
 *
 */
public class UserConfiguration
{

	String pdbFilePath;

	boolean isSplit;

	private boolean autoFetch;

	public UserConfiguration(){
		isSplit = true;
		autoFetch = true;
		 // accessing temp. OS directory:         
        String property = "java.io.tmpdir";

        String tempdir = System.getProperty(property);
        
        if ( !(tempdir.endsWith(PDBFileReader.lineSplit) ) )
           tempdir = tempdir + PDBFileReader.lineSplit;
		pdbFilePath = tempdir;
	}
	
	public String getPdbFilePath()
	{
		return pdbFilePath;
	}

	public void setPdbFilePath(String pdbFilePath)
	{
		this.pdbFilePath = pdbFilePath;
	}


	public boolean isSplit() {
		return isSplit;
	}

	public void setSplit(boolean isSplit) {
		this.isSplit = isSplit;
	}

	public boolean getAutoFetch() {
		return autoFetch;
	}

	public void setAutoFetch(boolean autoFetch) {
		this.autoFetch = autoFetch;
	}

	/** convert Configuration to an XML file so it can be serialized
	 * 
	 * @param pw
	 * @return XMLWriter
	 * @throws IOException
	 */
	public XMLWriter toXML(PrintWriter pw) 
	throws IOException
	{

		XMLWriter     xw = new PrettyXMLWriter( pw);

		toXML(xw);
		return xw ;
	}


	/** convert Configuration to an XML file so it can be serialized
	 * add to an already existing xml file.
	 * 
	 * @param xw the XML writer to use
	 * @return the writer again
	 * @throws IOException
	 */

	public XMLWriter toXML(XMLWriter xw) 
	throws IOException
	{
		xw.printRaw("<?xml version='1.0' standalone='no' ?>");
		//xw.printRaw("<!DOCTYPE " + XML_CONTENT_TYPE + " SYSTEM '" + XML_DTD + "' >");
		xw.openTag("JFatCatConfig");

		xw.openTag("PDBFILEPATH");
		xw.attribute("path", pdbFilePath);
		xw.attribute("split", isSplit +"" );
		xw.attribute("autofetch", autoFetch+"");
		xw.closeTag("PDBFILEPATH");

		xw.closeTag("JFatCatConfig");
		return xw ;

	}

	public static UserConfiguration fromStartupParams(StartupParameters params) {
		UserConfiguration config = new UserConfiguration();
		config.setPdbFilePath(params.getPdbFilePath());
		config.setAutoFetch(params.isAutoFetch());
		config.setSplit(params.isPdbDirSplit());
		return config;
	}




}
