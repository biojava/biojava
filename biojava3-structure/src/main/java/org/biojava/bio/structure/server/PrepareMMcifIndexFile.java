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
 * created at Oct 18, 2008
 */
package org.biojava.bio.structure.server;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

import org.biojava.bio.structure.io.MMCIFFileReader;
import org.biojava.bio.structure.io.StructureIOFile;

/** A class to index all files in a PDB installation
 *
 * @author Andreas Prlic
 * @since 1.7
 */
public class PrepareMMcifIndexFile extends PrepareIndexFile{

	public static Logger logger =  Logger.getLogger("org.biojava.bio.structure");


	public File[] getAllMMcif(File dir){
		if ( ! dir.isDirectory()){
			throw new IllegalArgumentException("path is not a directory " + dir);
		}


		String[] all = dir.list();
		List<File> pdbFiles = new ArrayList<File>();
		for (int i = 0 ; i < all.length;i++ ){
			// filenames are like 'pdb1234.ent.gz'
			String file = all[i];

			File f = new File(dir.getAbsolutePath()+ "/"+ file);

			if ( f.isDirectory()){
				File[] subPDBs = getAllMMcif(f);
				logger.info("got " + subPDBs.length + " files from subdir " + f );
				for (File file2 : subPDBs) {
					if (! pdbFiles.contains(file2)){
						pdbFiles.add(file2);
					}

				}
			} else {


				if ( (file.endsWith(".cif.gz")) || ( file.endsWith(".mmcif.gz"))){
					pdbFiles.add(new File(dir+File.separator + file));
				}
			}
		}

		return (File[]) pdbFiles.toArray(new File[pdbFiles.size()]);

	}


	/** prepare the index file for this installation
	 *
	 * @param installation
	 */
	@SuppressWarnings("deprecation")
	public void prepareIndexFileForInstallation(FlatFileInstallation installation)
	throws FileNotFoundException,IOException{

		File[] pdbfiles = getAllMMcif(installation.getFilePath());
		logger.info("found " + pdbfiles.length + " mmcif files");
		createMMcifInfoList(pdbfiles, installation.getPDBInfoFile(), installation.getChainInfoFile());


	}


	/** parses a set of PDB files and writes info into a file
	 * the file is tab separated and has the following columns:
	 * name length  resolution depositionDate modificationDate  technique title classification filename
	 *
	 * binaryDirectory: a directory in which binary files containing the atoms will be places, to provide a speedup
	 *
	 * This method needs to be run, before a DBSearch can be performed, since the files created by this method
	 * are required for the DBSearch
	 *
	 * @param pdbfiles
	 * @param outputFile
	 * @throws FileNotFoundException
	 * @throws IOException
	 */
	public  void createMMcifInfoList(File[] pdbfiles, File outputFile, File chainInfoFile)
	throws FileNotFoundException, IOException
	{

		//String outputfile = "/Users/ap3/WORK/PDB/rotated.pdb";

		FileOutputStream out= new FileOutputStream(outputFile);
		PrintStream p =  new PrintStream( out );
		PrintWriter pdbWriter = new PrintWriter(p);

		FileOutputStream cout = new FileOutputStream(chainInfoFile);
		PrintStream pc = new PrintStream(cout);
		PrintWriter chainWriter = new PrintWriter(pc);


		StructureIOFile pdbreader = new MMCIFFileReader();

		logPDBInfoFile(pdbWriter, chainWriter, pdbreader, pdbfiles);
	}


}
