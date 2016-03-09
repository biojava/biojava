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
 *
 *      http://www.biojava.org/

 * @author Andreas Prlic
 *
 */
package org.biojava.nbio.structure.align.webstart;


import org.biojava.nbio.structure.align.util.UserConfiguration;
import org.biojava.nbio.structure.io.LocalPDBDirectory.FetchBehavior;
import org.biojava.nbio.structure.io.LocalPDBDirectory.ObsoleteBehavior;
import org.xml.sax.Attributes;
import org.xml.sax.helpers.DefaultHandler;



/**
 * XML content handler for serialisation of RegistryConfiguration class
 */
public class ConfigXMLHandler extends DefaultHandler {

	UserConfiguration config ;

	/**
	 *
	 */
	public ConfigXMLHandler() {
		super();

		config         = new UserConfiguration();
	}

	@Override
	public void startElement (String uri, String name, String qName, Attributes atts){
		//System.out.println("new element >" + name + "< >" + qName+"<" + uri);
		if ( qName.equals("PDBFILEPATH")){

			String path = atts.getValue("path");
			// default path is system tmp...
			if ( path != null)
				config.setPdbFilePath(path);

			//Deprecated property; supported for backwards compatibility
			String autoFetch = atts.getValue("autoFetch");
			if(autoFetch == null || !autoFetch.equals("false")) {
				config.setFetchBehavior(FetchBehavior.DEFAULT);
			} else {
				config.setFetchBehavior(FetchBehavior.LOCAL_ONLY);
			}

			String fetchBehavior = atts.getValue("fetchBehavior");
			if(fetchBehavior == null) {
				config.setFetchBehavior(FetchBehavior.DEFAULT);
			} else {
				config.setFetchBehavior(FetchBehavior.valueOf(fetchBehavior));
			}
			String obsoleteBehavior = atts.getValue("obsoleteBehavior");
			if(obsoleteBehavior == null) {
				config.setObsoleteBehavior(ObsoleteBehavior.DEFAULT);
			} else {
				config.setObsoleteBehavior(ObsoleteBehavior.valueOf(obsoleteBehavior));
			}

			String fileFormat = atts.getValue("fileFormat");
			config.setFileFormat(UserConfiguration.PDB_FORMAT);
			if ( fileFormat != null) {
				if ( fileFormat.equals(UserConfiguration.MMCIF_FORMAT))
					config.setFileFormat(UserConfiguration.MMCIF_FORMAT);
			}

		}
	}







	public UserConfiguration getConfig() {
		return config ;
	}

}
