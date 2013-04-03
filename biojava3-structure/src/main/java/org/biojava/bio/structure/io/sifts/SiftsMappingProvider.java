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
package org.biojava.bio.structure.io.sifts;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import java.util.List;


import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.io.util.FileDownloadUtils;
import org.biojava3.core.util.InputStreamProvider;

public class SiftsMappingProvider {

	static String EBI_SIFTS_FILE_LOCATION = "ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/xml/%s.xml.gz";
	
	static String RCSB_SIFTS_FILE_LOCATION = "http://www.rcsb.org/pdb/files/%s.sifts.xml.gz";
		
	static String fileLoc = EBI_SIFTS_FILE_LOCATION;
		
	public static void main(String[] args){
		try {
			List<SiftsEntity> entities = getSiftsMapping("1gc1");
			
			for (SiftsEntity e : entities){
				System.out.println(e.getEntityId() + " " +e.getType());
				
				for ( SiftsSegment seg: e.getSegments()) {
					System.out.println(" Segment: " + seg.getSegId() + " " + seg.getStart() + " " + seg.getEnd()) ;
					
					for ( SiftsResidue res: seg.getResidues() ) {
						System.out.println("  " + res);
					}
				}
				
			}
		} catch (Exception e){
			e.printStackTrace();
		}
	}
	
	public static void setFileLocation(String myFileLocation){
		fileLoc = myFileLocation;
	}
	
	public static List<SiftsEntity> getSiftsMapping(String pdbId) throws IOException{
		// grab files from here:

		AtomCache cache = new AtomCache();
		
		String path = cache.getPath();

		String dirHash = pdbId.substring(1,3);
		String SIFTS_DIR = path + "/SIFTS/";
		
		File siftsDir = new File(SIFTS_DIR);
		
		if ( ! siftsDir.exists()) 
			siftsDir.mkdir();
		
		File hashDir = new File(SIFTS_DIR + dirHash);
		
		if ( ! hashDir.exists()){
			hashDir.mkdir();
		}
		File dest = new File( hashDir.getAbsolutePath() +"/"+ pdbId + ".sifts.xml.gz");
		
		if ( ! dest.exists()){
			String u = String.format(fileLoc,pdbId);
			URL url = new URL(u);
			FileDownloadUtils.downloadGzipCompressedFile(url, dest);
		}
		
		
		InputStreamProvider prov = new InputStreamProvider();
		InputStream is = prov.getInputStream(dest);
		SiftsXMLParser parser = new SiftsXMLParser();
		
		parser.parseXmlFile(is);
		
		//System.out.println(parser.getEntities());
		return parser.getEntities();
		
		
	}

	
}
