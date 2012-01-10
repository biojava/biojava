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
 * created at Sep 8, 2007
 */
package org.biojava.bio.structure.server;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.List;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Compound;
import org.biojava.bio.structure.DBRef;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.PDBHeader;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureImpl;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.io.PDBFileReader;
import org.biojava.bio.structure.io.StructureIOFile;

public class PrepareIndexFile {

	@SuppressWarnings("deprecation")
	private static final SimpleDateFormat dateFormat = FlatFileInstallation.dateFormat;

	
	public PrepareIndexFile(){

	}

	/** prepare the index file for this installation
	 *
	 * @param installation
	 */
	@SuppressWarnings("deprecation")
	public void prepareIndexFileForInstallation(FlatFileInstallation installation)
	throws FileNotFoundException,IOException{

		File[] pdbfiles = getAllPDB(installation.getFilePath());

		createPDBInfoList(pdbfiles, installation.getPDBInfoFile(), installation.getChainInfoFile());


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
	public  void createPDBInfoList(File[] pdbfiles, File outputFile, File chainInfoFile)
	throws FileNotFoundException, IOException
	{

		//String outputfile = "/Users/ap3/WORK/PDB/rotated.pdb";

		FileOutputStream out= new FileOutputStream(outputFile);
		PrintStream p =  new PrintStream( out );
		PrintWriter pdbWriter = new PrintWriter(p);

		FileOutputStream cout = new FileOutputStream(chainInfoFile);
		PrintStream pc = new PrintStream(cout);
		PrintWriter chainWriter = new PrintWriter(pc);


		PDBFileReader pdbreader = new PDBFileReader();

		logPDBInfoFile(pdbWriter, chainWriter, pdbreader, pdbfiles);
	}

	protected void logPDBInfoFile(PrintWriter pdbWriter, PrintWriter chainWriter, StructureIOFile pdbreader, File[] pdbfiles )
	throws IOException{

		pdbWriter.println("//pdbId\tnrCAAtoms\ttechnique\tresolution\tdepDate\tmodDate\ttitle\tclassification\ttime\tpath");
		chainWriter.println("//pdbId\tchainId\tseqResLength\tatomLength\taminoLength\thetLength\tnucleotideLength\tmolName\tmolId\tdbrefs");
		int l = pdbfiles.length;


        long loopStart = System.currentTimeMillis();

		for ( int i = 0 ; i < l ; i++){
			//if (i != 36)
			//	continue;
			/* if ( i < 9424)
                continue;
            if ( i > 9426)
                break;
			*/
           // if ( i< 16055)
             //   continue;

			File f = pdbfiles[i];

            //String pdb = f.getName().substring(3,7);


			long startTime = System.currentTimeMillis();
            System.out.println("# " + i + " / " + l + " " + f );
            //System.out.println("getting " + f);
            Structure s = null;

            s = pdbreader.getStructure(f);

			long stopTime = System.currentTimeMillis();

			logPDBInfo(pdbWriter, s, startTime, stopTime,f);
			logChainInfo(chainWriter, s);

		}
        long loopEnd = System.currentTimeMillis();
		pdbWriter.flush();
		pdbWriter.close();
		chainWriter.flush();
		chainWriter.close();

        long time = ((loopEnd-loopStart) / (long) (1000*60));
        System.out.println("loop took: " + time + " minutes" );
	}

	/** 			get the matching compound for this chain
	 *
	 * @param compounds
	 * @param c
	 * @return
	 */
	private  Compound getCompoundForChain(List<Compound> compounds, Chain c){
		for(Compound comp : compounds){
			List<String>chainIds =  comp.getChainId();
			if ( chainIds == null)
				continue;
			if ( chainIds.contains(c.getChainID()))
				return comp;
		}

		return null;
	}

	private String getDBRefStringForChain(List<DBRef> dbrefs, Chain c){
		List<String> dbIds = new ArrayList<String>();
		for (DBRef dbref : dbrefs){
			if (dbref.getChainId().equals(c.getChainID())){
				dbIds.add(dbref.getDbIdCode());
			}
		}
		StringBuffer buf = new StringBuffer();
		int pos = 1;
		for ( String id : dbIds){
			buf.append(id);
			if ( pos < dbIds.size())
				buf.append(":");
		}

		return buf.toString();

	}
	private  void logChainInfo(PrintWriter chainWriter, Structure s ){

		String pdbCode = s.getPDBCode();
		List<Chain> chains = s.getChains(0);
		List<Compound> compounds = s.getCompounds();
		List<DBRef> dbrefs = s.getDBRefs();

		for ( Chain c : chains){

			List<Group> agr = c.getAtomGroups("amino");
			List<Group> hgr = c.getAtomGroups("hetatm");
			List<Group> ngr = c.getAtomGroups("nucleotide");

			Compound comp =  getCompoundForChain(compounds, c);
			String molName = "-";
			String molId   = "-";
			if (comp != null){
				molName = comp.getMolName();
				molId   = comp.getMolId();
			}
			String dbRefString = getDBRefStringForChain(dbrefs,c);

			String str = String.format("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t",
					pdbCode,
					c.getChainID(),
					c.getSeqResLength(),
					c.getAtomLength(),
					agr.size(),
					hgr.size(),
					ngr.size(),
					molName,
					molId,
					dbRefString
			);
			chainWriter.println(str);
		}



	}

	@SuppressWarnings("deprecation")
	private  void logPDBInfo(PrintWriter pdbWriter, Structure s, long startTime,long stopTime, File path ){
		// only used first model in nmrs ...
		if ( s.isNmr()){
			List<Chain> chains = s.getModel(0);
			Structure newsubject = new StructureImpl();
			newsubject.setPDBCode(s.getPDBCode());
			newsubject.setHeader(s.getHeader());
			newsubject.setPDBHeader(s.getPDBHeader());
			newsubject.setDBRefs(s.getDBRefs());
			newsubject.addModel(chains);
			newsubject.setNmr(true);
			s = newsubject;
		}

		//System.out.println((stopTime-startTime) / (float)1000);
		Atom[] ca = StructureTools.getAtomCAArray(s);

		PDBHeader header = s.getPDBHeader();
		String infoline = String.format("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t\t%s",
				s.getPDBCode(),
				ca.length,
				header.getTechnique(),
				header.getResolution(),
				dateFormat.format(header.getDepDate()),
				dateFormat.format(header.getModDate()),
				header.getTitle(),
				header.getClassification(),
				((stopTime-startTime) / (float) 1000),
				path.getAbsolutePath()
		);



		//System.out.println(i +"/" + pdbfiles.length + " " + infoline);
		//System.out.println(s.getHeader());

		pdbWriter.println(infoline);
	}


	/** get all PDBfiles from a directory
	 *
	 * @param dir the directory that contains all PDB files
	 * @return an array of PDB Files
	 */
	public  File[] getAllPDB(File dir)  {



		if ( ! dir.isDirectory()){
			throw new IllegalArgumentException("path is not a directory " + dir);
		}

		String[] all = dir.list();
		List<File> pdbFiles = new ArrayList<File>();
		for (int i = 0 ; i < all.length;i++ ){
			// filenames are like 'pdb1234.ent.gz'
			String file = all[i];
			if ( (file.endsWith(".pdb.gz")) || ( file.endsWith(".ent.gz"))){
				pdbFiles.add(new File(dir+File.separator + file));
			}
		}

		return (File[]) pdbFiles.toArray(new File[pdbFiles.size()]);
	}




}
