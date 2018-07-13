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
 * Created on Feb 22, 2012
 * Created by Andreas Prlic
 *
 * @since 3.0.2
 */
package org.biojava.nbio.structure.io.sifts;

import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.core.util.FileDownloadUtils;
import org.biojava.nbio.core.util.InputStreamProvider;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import java.util.List;

public class SiftsMappingProvider {

	private final static Logger logger = LoggerFactory.getLogger(SiftsMappingProvider.class);


	private static final String EBI_SIFTS_FILE_LOCATION = "http://ftp.ebi.ac.uk/pub/databases/msd/sifts/xml/%s.xml.gz";

	private static String fileLoc = EBI_SIFTS_FILE_LOCATION;

	public static void setFileLocation(String myFileLocation){
		fileLoc = myFileLocation;
	}

	/**
	 * Return the SIFTS mappings by getting the info from individual SIFTS xml files at URL {@value EBI_SIFTS_FILE_LOCATION}
	 * @param pdbId the pdb identifier
	 * @return
	 * @throws IOException if problems downloading or parsing the file
	 */
	public static List<SiftsEntity> getSiftsMapping(String pdbId) throws IOException{
		// grab files from here:

		AtomCache cache = new AtomCache();

		String path = cache.getCachePath();

		pdbId = pdbId.toLowerCase();

		String dirHash = pdbId.substring(1,3);
		File siftsDir = new File(path , "SIFTS");


		if ( ! siftsDir.exists()) {
			logger.info("Creating directory {}", siftsDir.toString());
			siftsDir.mkdir();
		}

		File hashDir = new File(siftsDir, dirHash);

		if ( ! hashDir.exists()){
			logger.info("Creating directory {}", hashDir.toString());
			hashDir.mkdir();
		}
		File dest = new File( hashDir, pdbId + ".sifts.xml.gz");

		logger.debug("testing SIFTS file " + dest.getAbsolutePath());


		if ( ! dest.exists()){
			String u = String.format(fileLoc,pdbId);
			URL url = new URL(u);
			logger.debug("Downloading SIFTS file {} to {}",url,dest);
			FileDownloadUtils.downloadFile(url, dest);
		}

		InputStreamProvider prov = new InputStreamProvider();
		InputStream is = prov.getInputStream(dest);
		SiftsXMLParser parser = new SiftsXMLParser();

		parser.parseXmlFile(is);

		//System.out.println(parser.getEntities());
		return parser.getEntities();


	}


}
