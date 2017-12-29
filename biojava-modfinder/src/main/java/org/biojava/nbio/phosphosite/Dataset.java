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
package org.biojava.nbio.phosphosite;

import org.biojava.nbio.structure.align.util.AtomCache;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;
import java.net.URL;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.List;

/**
 * Phosphosite is available under the PhosphoSitePlus® is licensed under Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License and is freely available for non-commercial purposes from
 *
 * http://www.phosphosite.org/staticDownloads.do
 *
 * Please acknowledge PhosphoSitePlus®, www.phosphosite.org" at appropriate locations.
 *
 * Please cite : “Hornbeck PV, Kornhauser JM, Tkachev S, Zhang B, Skrzypek E, Murray B, Latham V, Sullivan M (2012) PhosphoSitePlus: a comprehensive resource for investigating the structure and function of experimentally determined post-translational modifications in man and mouse. Nucleic Acids Res. 40(Database issue), D261–70.”.
 *
 *
 *
 * Created by ap3 on 31/10/2014.
 */
public class Dataset {

	private static final Logger logger = LoggerFactory.getLogger(Dataset.class);

	public static final String ACETYLATION = "https://www.phosphosite.org/downloads/Acetylation_site_dataset.gz";

	public static final String DISEASE_ASSOC = "https://www.phosphosite.org/downloads/Disease-associated_sites.gz";

	public static final String METHYLATION = "https://www.phosphosite.org/downloads/Methylation_site_dataset.gz";

	public static final String PHOSPHORYLATION = "https://www.phosphosite.org/downloads/Phosphorylation_site_dataset.gz";

	public static final String REGULATORY = "https://www.phosphosite.org/downloads/Regulatory_sites.gz";

	public static final String SUMOYLATION = "https://www.phosphosite.org/downloads/Sumoylation_site_dataset.gz";

	public static final String UBIQUITINATION = "https://www.phosphosite.org/downloads/Ubiquitination_site_dataset.gz";


	public Dataset(){


	}

	private String[] getRemoteFiles(){
		String[] files = new String[]{ACETYLATION,DISEASE_ASSOC,METHYLATION,PHOSPHORYLATION,REGULATORY,SUMOYLATION,UBIQUITINATION};


		return files;
	}

	public File[] getLocalFiles(){

		String[] rfiles = getRemoteFiles();


		File dir = getLocalDir();

		List<File> files = new ArrayList<File>();
		for ( String f : rfiles) {


			int slashIndex = f.lastIndexOf("/");

			String fileName = f.substring(slashIndex);

			File localFile = new File(dir+"/" + fileName);

			if (  localFile.exists()){
				files.add(localFile);
			}

		}

		return files.toArray(new File[files.size()]);
	}


	public File getLocalDir(){
		AtomCache cache = new AtomCache();

		String path = cache.getCachePath();

		File dir = new File(path+"/phosphosite");

		return dir;
	}

	public void download(){

		logger.warn("Downloading data from www.phosposite.org. Data is under CC-BY-NC-SA license. Please link to site and cite: ");
		logger.warn("Hornbeck PV, Kornhauser JM, Tkachev S, Zhang B, Skrzypek E, Murray B, Latham V, Sullivan M (2012) PhosphoSitePlus: a comprehensive resource for investigating the structure and function of experimentally determined post-translational modifications in man and mouse. Nucleic Acids Res. 40(Database issue), D261–70.");

		File dir = getLocalDir();

		if ( ! dir.exists()) {

			// need to download all...

			dir.mkdir();


		}

		String[] files = getRemoteFiles();

		for ( String f : files){

			try {


				int slashIndex = f.lastIndexOf("/");

				String fileName = f.substring(slashIndex);

				File localFile = new File(dir+"/" + fileName);

				if ( ! localFile.exists()){

					URL u = new URL(f);
					downloadFile(u, localFile);
				}


			} catch (Exception e){

				e.printStackTrace();
			}


		}

	}

	public void downloadFile(URL u, File localFile) throws IOException {

		logger.info("Downloading " + u);

		File tmp = File.createTempFile("tmp","phosphosite");

		InputStream is = u.openStream();

		BufferedInputStream in = new BufferedInputStream(is);

		FileOutputStream w = new FileOutputStream(tmp);

		int i= 0;
		byte[] bytesIn = new byte[300000];
		while ((i = in.read(bytesIn)) >= 0) {
			w.write(bytesIn,0,i);
		}
		in.close();
		w.close();


		// now copy  tmp file to localFile
		copyFile(tmp, localFile);

	}



	public static void copyFile(File src, File dst) throws IOException
	{

		Files.copy(src.toPath(), dst.toPath(), StandardCopyOption.REPLACE_EXISTING);

	}


	public static void main(String[] args) {

		Dataset ds = new Dataset();

		ds.download();

		try {

			for (File f : ds.getLocalFiles()) {

				logger.info(f.getAbsolutePath());

				List<Site> sites = Site.parseSites(f);

				logger.info("Got " + sites.size() + " sites");
				for (Site s : sites) {
					if (s.getUniprot().equals("P50225") || s.getUniprot().equals("P48025")) {
						logger.info(s.toString());
					}
				}

			}


		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
