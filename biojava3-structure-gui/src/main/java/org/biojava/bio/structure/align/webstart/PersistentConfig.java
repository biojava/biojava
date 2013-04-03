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
 * Created on 20.09.2004
 * @author Andreas Prlic
 *
 */

package org.biojava.bio.structure.align.webstart;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.net.URL;

import javax.jnlp.BasicService;
import javax.jnlp.FileContents;
import javax.jnlp.PersistenceService;
import javax.jnlp.ServiceManager;
import javax.jnlp.UnavailableServiceException;
import javax.xml.parsers.ParserConfigurationException;
import javax.xml.parsers.SAXParser;
import javax.xml.parsers.SAXParserFactory;

import org.biojava.bio.structure.align.util.UserConfiguration;
import org.biojava3.core.util.XMLWriter;

import org.xml.sax.InputSource;
import org.xml.sax.XMLReader;


/** a class to store the config using the Java Web Start
 * PersistenService.
 * @author Andreas Prlic
 */
public class PersistentConfig
{

	PersistenceService ps; 
	BasicService bs      ; 

	public PersistentConfig()
	throws UnavailableServiceException
	{
		try {
			ps = (PersistenceService)ServiceManager.lookup("javax.jnlp.PersistenceService"); 
			bs = (BasicService)ServiceManager.lookup("javax.jnlp.BasicService");
		} catch (Exception e){
			System.err.println("Can't init webstart - persistent configuration. " + e.getMessage() );
		}
	}

	/** writes the configuration
	 * 
	 * @param config
	 */
	public void save(UserConfiguration config ) {
		if (ps != null && bs != null) { 
			// Persistent Service is available, running as javaws
			saveWebStart(config) ;
		} else {
			System.err.println("can not save using persistentservice!");
		}
	}

	private void saveWebStart(UserConfiguration config ){
		//System.out.println("saving webstart");


		try {

			// find all the muffins for our URL
			URL codebase = bs.getCodeBase(); 

			FileContents fc = null ;


			try {		
				// test if persistent storage already created

				fc = ps.get(codebase);

				ps.delete(codebase);

			} catch (IOException e){
			}

			// seems not, create it first
			ps.create(codebase,3000000);
			fc = ps.get(codebase);


			OutputStream os = fc.getOutputStream(true); 

			//StringWriter sw = new StringWriter();
			//StringWriter stw = new StringWriter(os)   ;
			PrintWriter pw = new PrintWriter(os,true);
			XMLWriter xw = config.toXML(pw);


			pw.flush();
			os.flush();

			xw.close();
			pw.close();
			os.close(); 


		} catch (Exception e) {
			System.err.println(e.getMessage());
		}
	}


	/** loads Config from PersistenceService
	 *  returns null if no PErsistenceService has been created ...
	 *  
	 *  @return WebStartConfiguration
	 */
	public UserConfiguration load() {
		if (ps != null && bs != null) { 
			// Persistent Service is available, running as javaws
			return loadWebStart() ;
		} else {
			System.err.println("can not load from persistentservice!");
		}
		return null ;
	}


	/** loads Config from PersistenceService
	 *  returns null if no PErsistenceService has been created ...
	 */
	private UserConfiguration loadWebStart() {
		UserConfiguration config = null;
		try { 
			URL codebase = bs.getCodeBase(); 

			FileContents fc = null ;

			try {

				fc = ps.get(codebase);
			} catch (IOException e){
				// has not been created, so nothing can be loaded ...
				e.printStackTrace();
				return null ;
			}	


			// parse the XML file ...
			InputStream stream = fc.getInputStream();
			config = parseConfigFile(stream);

		} catch (Exception e) {
			System.err.println(e.getMessage());
		}

		return config ;

	}


	private UserConfiguration parseConfigFile(InputStream inStream) {

		try {
			SAXParserFactory spfactory =
				SAXParserFactory.newInstance();

			SAXParser saxParser = null ;

			try{
				saxParser =
					spfactory.newSAXParser();
			} catch (ParserConfigurationException e) {
				e.printStackTrace();
			}

			XMLReader xmlreader = saxParser.getXMLReader();

			ConfigXMLHandler cont_handle = new ConfigXMLHandler();
			xmlreader.setContentHandler(cont_handle);
			xmlreader.setErrorHandler(new org.xml.sax.helpers.DefaultHandler());

			InputSource insource = new InputSource() ;
			insource.setByteStream(inStream);

			// the actual parsing starts now ...
			xmlreader.parse(insource);


			UserConfiguration config = cont_handle.getConfig();
			return config ;

		} catch (Exception e){
			e.printStackTrace();
			return null;
		}

	}
}
