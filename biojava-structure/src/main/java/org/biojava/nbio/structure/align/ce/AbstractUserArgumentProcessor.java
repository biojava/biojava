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
 * Created on Nov 2, 2009
 * Author: Andreas Prlic
 *
 */

package org.biojava.nbio.structure.align.ce;


import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.align.MultiThreadedDBSearch;
import org.biojava.nbio.structure.align.StructureAlignment;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.util.*;
import org.biojava.nbio.structure.align.xml.AFPChainXMLConverter;
import org.biojava.nbio.structure.io.LocalPDBDirectory.FetchBehavior;
import org.biojava.nbio.structure.io.PDBFileReader;

import java.beans.Introspector;
import java.io.*;
import java.lang.reflect.InvocationTargetException;
import java.net.URL;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;


/**
 * Base class for a new structure alignment CLI.
 *
 * <p>To add a new StructureAlignment with a CLI similar to CE or FATCAT,
 * <ol>
 * <li>Implement StructureAlignment with the main algorithm
 * <li>Implement ConfigStrucAligParams. This provides the parameters for the GUI.
 * <li>Subclass StartupParameters (can be an inner class) with the same parameters
 *     described in the ConfigStrucAligParams.
 * <li>Subclass AbstractUserArgumentProcessor, with the getStartupParams() method
 *     returning a fresh instance of the custom StartupParameters
 * <li>Implement the getParameters() method to copy values from the StartupParameters
 *     to the ConfigStrucAligParams.
 * </ol>
 *
 * <p>Note that reflection is used in a number of places, so the CLI argument names
 * must match the get/set functions in both parameter beans.
 * <ul>
 * <li>AlignmentGUI automatically takes parameter names and types from the
 *     ConfigStrucAligParams
 * <li>AbstractUserArgumentProcessor also takes parameter names and help descriptions
 *     from ConfigStrucAligParams, but it saves arguments to the StartupParameter
 *     bean.
 * </ul>
 * This means that both beans should be kept in sync.
 *
 * @author Andreas
 * @author Spencer
 *
 */
public abstract class AbstractUserArgumentProcessor implements UserArgumentProcessor {

	public static String newline = System.getProperty("line.separator");

	protected StartupParameters params ;

	public static final List<String> mandatoryArgs= new ArrayList<String>();

	protected AbstractUserArgumentProcessor(){
		params = getStartupParametersInstance();
	}

	/**
	 * StartupParameters is a bean to store all the possible
	 * command line parameters.
	 *
	 * The `StartupParameter` class contains the basic set of CLI parameters
	 * common to all `StructureAligmnent`s. This method should return a subclass
	 * of StartupParameters which has been extended to store values for all
	 * additional parameters.
	 * @return A new instance of the correct StartupParameters subclass
	 */
	protected abstract StartupParameters getStartupParametersInstance();

	public abstract StructureAlignment getAlgorithm();
	public abstract Object getParameters();
	public abstract String getDbSearchLegend();

	@Override
	public void process(String[] argv){

		printAboutMe();

//		if(argv.length == 0 ) {
//			System.out.println(printHelp());
//			return;
//		}

		for (int i = 0 ; i < argv.length; i++){
			String arg   = argv[i];

			// help string
			if(arg.equalsIgnoreCase("-h") || arg.equalsIgnoreCase("-help")
					|| arg.equalsIgnoreCase("--help") )
			{
				System.out.println(printHelp());
				return;
			}
			// version
			if(arg.equalsIgnoreCase("-version") || arg.equalsIgnoreCase("--version")) {
				StructureAlignment alg = getAlgorithm();
				System.out.println(alg.getAlgorithmName() + " v." + alg.getVersion() );
				return;
			}

			String value = null;
			if ( i < argv.length -1)
				value = argv[i+1];

			// if value starts with - then the arg does not have a value.
			if (value != null && value.startsWith("-"))
				value = null;
			else
				i++;


			String[] tmp = {arg,value};

			//System.out.println(arg + " " + value);

			try {

				CliTools.configureBean(params, tmp);

			} catch (ConfigurationException e){
				System.err.println("Error: "+e.getLocalizedMessage());
				System.exit(1); return;
			}
		}

		if ( params.getPdbFilePath() != null){
			System.setProperty(UserConfiguration.PDB_DIR,params.getPdbFilePath());
		}

		if ( params.getCacheFilePath() != null){
			System.setProperty(UserConfiguration.PDB_CACHE_DIR,params.getCacheFilePath());
		}

		if ( params.isShowMenu()){
			System.err.println("showing menu...");
			try {
				GuiWrapper.showAlignmentGUI();
			} catch (Exception e){
				System.err.println(e.getMessage());
				e.printStackTrace();
			}
		}

		if ( params.getShowDBresult() != null){
			// user wants to view DB search results:


			System.err.println("showing DB results...");
			try {
				GuiWrapper.showDBResults(params);
			} catch (Exception e){
				System.err.println(e.getMessage());
				e.printStackTrace();
			}

		}

		String pdb1  = params.getPdb1();
		String file1 = params.getFile1();


		try {
			if (pdb1 != null || file1 != null){
				runPairwise();
				return;
			}

			if ( params.getAlignPairs() != null){
				runDBSearch();
				return;
			}

			if ( params.getSearchFile() != null){
				runDBSearch();
				return;
			}
		} catch (ConfigurationException e) {
			System.err.println(e.getLocalizedMessage());
			System.exit(1); return;
		}

		System.out.println(printHelp());
		System.err.println("Error: insufficient arguments.");
		System.exit(1); return;
	}





	public static void printAboutMe() {
		try {
			ResourceManager about = ResourceManager.getResourceManager("about");

			String version = about.getString("project_version");
			String build   = about.getString("build");

			System.out.println("Protein Comparison Tool " + version + " " + build);
		} catch (Exception e){
			e.printStackTrace();
		}


	}


	private void runDBSearch() throws ConfigurationException{


		String pdbFilePath = params.getPdbFilePath();

		if ( pdbFilePath == null || pdbFilePath.equals("")){

			UserConfiguration c = new UserConfiguration();
			pdbFilePath = c.getPdbFilePath();
			System.err.println("You did not specify the -pdbFilePath parameter. Defaulting to "+pdbFilePath+".");
		}

		String cacheFilePath = params.getCacheFilePath();

		if ( cacheFilePath == null || cacheFilePath.equals("")){
			cacheFilePath = pdbFilePath;

		}


		AtomCache cache = new AtomCache(pdbFilePath, pdbFilePath);

		String alignPairs = params.getAlignPairs();

		String searchFile = params.getSearchFile();

		if ( alignPairs == null || alignPairs.equals("")) {
			if ( searchFile == null || searchFile.equals("")){
				throw new ConfigurationException("Please specify -alignPairs or -searchFile !");
			}
		}

		String outputFile = params.getOutFile();

		if ( outputFile == null || outputFile.equals("")){
			throw new ConfigurationException("Please specify the mandatory argument -outFile!");
		}

		System.out.println("running DB search with parameters: " + params);

		if ( alignPairs != null && ! alignPairs.equals("")) {
			runAlignPairs(cache, alignPairs, outputFile);
		}  else {
			// must be a searchFile request...

			int useNrCPUs = params.getNrCPU();

			runDbSearch(cache,searchFile, outputFile, useNrCPUs, params);
		}
	}


	/** Do a DB search with the input file against representative PDB domains
	 *
	 * @param cache
	 * @param searchFile
	 * @param outputFile
	 * @throws ConfigurationException
	 */
	private void runDbSearch(AtomCache cache, String searchFile,
			String outputFile,int useNrCPUs, StartupParameters params) throws ConfigurationException {


		System.out.println("will use " + useNrCPUs + " CPUs.");

		PDBFileReader reader = new PDBFileReader();
		Structure structure1 = null ;
		try {
			structure1 = reader.getStructure(searchFile);
		} catch (IOException e) {
			throw new ConfigurationException("could not parse as PDB file: " + searchFile);
		}

		File searchF = new File(searchFile);
		String name1 = "CUSTOM";



		StructureAlignment algorithm =  getAlgorithm();

		MultiThreadedDBSearch dbSearch = new MultiThreadedDBSearch(name1,
				structure1,
				outputFile,
				algorithm,
				useNrCPUs,
				params.isDomainSplit());

		dbSearch.setCustomFile1(searchF.getAbsolutePath());

		dbSearch.run();


	}


	private void runAlignPairs(AtomCache cache, String alignPairs,
			String outputFile) {
		try {
			File f = new File(alignPairs);

			BufferedReader is = new BufferedReader (new InputStreamReader(new FileInputStream(f)));

			BufferedWriter out = new BufferedWriter(new FileWriter(outputFile, true));

			StructureAlignment algorithm =  getAlgorithm();

			String header = "# algorithm:" + algorithm.getAlgorithmName();
			out.write(header);
			out.write(newline);

			out.write("#Legend: " + newline );
			String legend = getDbSearchLegend();
			out.write(legend + newline );
			System.out.println(legend);
			String line = null;
			while ( (line = is.readLine()) != null){
				if ( line.startsWith("#"))
					continue;

				String[] spl = line.split(" ");

				if ( spl.length != 2) {
					System.err.println("wrongly formattted line. Expected format: 4hhb.A 4hhb.B but found " + line);
					continue;
				}

				String pdb1 = spl[0];
				String pdb2 = spl[1];


				Structure structure1 = cache.getStructure(pdb1);
				Structure structure2 = cache.getStructure(pdb2);

				Atom[] ca1;
				Atom[] ca2;


				ca1 = StructureTools.getRepresentativeAtomArray(structure1);
				ca2 = StructureTools.getRepresentativeAtomArray(structure2);

				Object jparams = getParameters();

				AFPChain afpChain;

				afpChain = algorithm.align(ca1, ca2, jparams);
				afpChain.setName1(pdb1);
				afpChain.setName2(pdb2);

				String result = getDbSearchResult(afpChain);
				out.write(result);
				System.out.print(result);

				checkWriteFile(afpChain,ca1,ca2,true);
			}

			out.close();
			is.close();
		} catch(Exception e){
			e.printStackTrace();
		}
	}


	private void runPairwise() throws ConfigurationException{

		String name1 = params.getPdb1();
		String file1 = params.getFile1();

		if ( name1 == null && file1 == null){
			throw new ConfigurationException("You did not specify the -pdb1 or -file1 parameter. Can not find query PDB id for alignment.");
		}

		if ( file1 == null) {
			if ( name1.length() < 4) {
				throw new ConfigurationException("-pdb1 does not look like a PDB ID. Please specify PDB code or PDB.chainName.");
			}
		}



		String name2 = params.getPdb2();
		String file2 = params.getFile2();
		if ( name2 == null && file2 == null ){
			throw new ConfigurationException("You did not specify the -pdb2 or -file2 parameter. Can not find target PDB id for alignment.");
		}

		if ( file2 == null ){
			if ( name2.length() < 4) {
				throw new ConfigurationException("-pdb2 does not look like a PDB ID. Please specify PDB code or PDB.chainName.");
			}
		}

		// first load two example structures

		Structure structure1 = null;
		Structure structure2 = null;

		String path = params.getPdbFilePath();

		if ( file1 == null || file2 == null) {
			if ( path == null){
				UserConfiguration c = new UserConfiguration();
				path = c.getPdbFilePath();
				System.err.println("You did not specify the -pdbFilePath parameter. Defaulting to "+path+".");
			}

			AtomCache cache = new AtomCache(path, path);
			if(params.isAutoFetch()) {
				cache.setFetchBehavior(FetchBehavior.DEFAULT);
			} else {
				cache.setFetchBehavior(FetchBehavior.LOCAL_ONLY);
			}
			structure1 = getStructure(cache, name1, file1);
			structure2 = getStructure(cache, name2, file2);
		} else {

			structure1 = getStructure(null, name1, file1);
			structure2 = getStructure(null, name2, file2);
		}



		if ( structure1 == null){
			System.err.println("structure 1 is null, can't run alignment.");
			System.exit(1); return;
		}

		if ( structure2 == null){
			System.err.println("structure 2 is null, can't run alignment.");
			System.exit(1); return;
		}

		if ( name1 == null) {
			name1 = structure1.getName();
		}
		if ( name2 == null) {
			name2 = structure2.getName();
		}

		//                   default:      new:
		// 1buz - 1ali : time: 8.3s eqr 68 rmsd 3.1 score 161 | time 6.4 eqr 58 rmsd 3.0 scre 168
		// 5pti - 1tap : time: 6.2s eqr 48 rmsd 2.67 score 164 | time 5.2 eqr 49 rmsd 2.9 score 151
		// 1cdg - 8tim
		// 1jbe - 1ord
		// 1nbw.A - 1kid
		// 1t4y - 1rp5


		try {

			Atom[] ca1;
			Atom[] ca2;

			ca1 = StructureTools.getRepresentativeAtomArray(structure1);
			ca2 = StructureTools.getRepresentativeAtomArray(structure2);

			StructureAlignment algorithm =  getAlgorithm();
			Object jparams = getParameters();

			AFPChain afpChain;

			afpChain = algorithm.align(ca1, ca2, jparams);
			afpChain.setName1(name1);
			afpChain.setName2(name2);

			if ( params.isShow3d()){

				if (! GuiWrapper.isGuiModuleInstalled()) {
					System.err.println("The biojava-structure-gui module is not installed. Please install!");
				} else {

					try {

						Object jmol = GuiWrapper.display(afpChain,ca1,ca2);

						GuiWrapper.showAlignmentImage(afpChain, ca1,ca2,jmol);

					} catch (Exception e){

						System.err.println(e.getMessage());
						e.printStackTrace();
					}
					//StructureAlignmentJmol jmol = algorithm.display(afpChain,ca1,ca2,hetatms1, nucs1, hetatms2, nucs2);

					//String result = afpChain.toFatcat(ca1, ca2);
					//String rot = afpChain.toRotMat();

					//DisplayAFP.showAlignmentImage(afpChain, ca1,ca2,jmol);
				}
			}


			checkWriteFile(afpChain,ca1, ca2, false);




			if ( params.isPrintXML()){
				String fatcatXML = AFPChainXMLConverter.toXML(afpChain,ca1,ca2);
				System.out.println(fatcatXML);
			}
			if ( params.isPrintFatCat()) {
				// default output is to XML on sysout...
				System.out.println(afpChain.toFatcat(ca1, ca2));
			}
			if ( params. isPrintCE()){
				System.out.println(afpChain.toCE(ca1, ca2));
			}

		} catch (IOException e) {
			e.printStackTrace();
			System.exit(1); return;
		} catch (ClassNotFoundException e) {
			e.printStackTrace();
			System.exit(1); return;
		} catch (NoSuchMethodException e) {
			e.printStackTrace();
			System.exit(1); return;
		} catch (InvocationTargetException e) {
			e.printStackTrace();
			System.exit(1); return;
		} catch (IllegalAccessException e) {
			e.printStackTrace();
			System.exit(1); return;
		} catch (StructureException e) {
			e.printStackTrace();
			System.exit(1); return;
		}
	}

	/** check if the result should be written to the local file system
	 *
	 * @param params2
	 * @param afpChain
	 * @param ca1
	 * @param ca2
	 * @throws IOException If an error occurs when writing the afpChain to XML
	 * @throws ClassNotFoundException If an error occurs when invoking jmol
	 * @throws NoSuchMethodException If an error occurs when invoking jmol
	 * @throws InvocationTargetException If an error occurs when invoking jmol
	 * @throws IllegalAccessException If an error occurs when invoking jmol
	 * @throws StructureException 
	 */
	private void checkWriteFile( AFPChain afpChain, Atom[] ca1, Atom[] ca2, boolean dbsearch) throws IOException, ClassNotFoundException, NoSuchMethodException, InvocationTargetException, IllegalAccessException, StructureException
			{
		String output = null;
		if ( params.isOutputPDB()){
			if (! GuiWrapper.isGuiModuleInstalled()) {
				System.err.println("The biojava-structure-gui module is not installed. Please install!");
				output = AFPChainXMLConverter.toXML(afpChain,ca1,ca2);
			} else {

				Structure tmp = AFPAlignmentDisplay.createArtificalStructure(afpChain, ca1, ca2);
				output = "TITLE  " + afpChain.getAlgorithmName() + " " + afpChain.getVersion()  + " ";
				output += afpChain.getName1() + " vs. " + afpChain.getName2();
				output += newline;
				output += tmp.toPDB();
			}
		} else  if ( params.getOutFile() != null) {
			// output by default is XML
			// write the XML to a file...
			output = AFPChainXMLConverter.toXML(afpChain,ca1,ca2);

		} else if ( params.getSaveOutputDir() != null){
			output = AFPChainXMLConverter.toXML(afpChain,ca1,ca2);
		}

		// no output requested.
		if ( output == null)
			return;

		String fileName = null;

		if ( dbsearch ){
			if ( params.getSaveOutputDir() != null) {

				// we currently don't have a naming convention for how to store results for custom files
				// they will be re-created on the fly
				if ( afpChain.getName1().startsWith("file:") || afpChain.getName2().startsWith("file:"))
					return;
				fileName = params.getSaveOutputDir();
				fileName += getAutoFileName(afpChain);

			} else {
				return;
			}

			//
			//else {
			//	fileName = getAutoFileName(afpChain);
			//}
		} else

			if ( params.getOutFile() != null) {
				fileName = params.getOutFile();
			}

		if (fileName == null) {
			System.err.println("Can't write outputfile. Either provide a filename using -outFile or set -autoOutputFile to true .");
			System.exit(1); return;
		}
		//System.out.println("writing results to " + fileName + " " + params.getSaveOutputDir());

		FileOutputStream out; // declare a file output object
		PrintStream p; // declare a print stream object

			// Create a new file output stream
			out = new FileOutputStream(fileName);

			// Connect print stream to the output stream
			p = new PrintStream( out );

			p.println (output);

			p.close();



	}


	private String getAutoFileName(AFPChain afpChain){
		String fileName =afpChain.getName1()+"_" + afpChain.getName2()+"_"+afpChain.getAlgorithmName();

		if (params.isOutputPDB() )
			fileName += ".pdb";
		else
			fileName += ".xml";
		return fileName;
	}


	private Structure getStructure(AtomCache cache, String name1, String file)
	{

		PDBFileReader reader = new PDBFileReader();
		if ( file != null ){
			try {
				// check if it is a URL:
				try {
					URL url = new URL(file);
					System.out.println(url);

					Structure s = reader.getStructure(url);

					return fixStructureName(s,file);

				} catch ( Exception e){
					System.err.println(e.getMessage());
				}
				File f= new File(file);
				System.out.println("file from local " + f.getAbsolutePath());
				Structure s= reader.getStructure(f);
				return fixStructureName(s, file);
			} catch (Exception e){
				System.err.println("general exception:" + e.getMessage());
				System.err.println("unable to load structure from " + file);
				return null;
			}
		}
		try {
			Structure s = cache.getStructure(name1);
			return s;
		} catch ( Exception e){
			System.err.println(e.getMessage());
			System.err.println("unable to load structure from dir: " + cache.getPath() + "/"+ name1);
			return null;
		}

	}

	/** apply a number of rules to fix the name of the structure if it did not get set during loading.
	 *
	 * @param s
	 * @param file
	 * @return
	 */
	private Structure fixStructureName(Structure s, String file) {

		if ( s.getName() != null && (! s.getName().equals("")))
			return s;

		s.setName(s.getPDBCode());

		if ( s.getName() == null || s.getName().equals("")){
			File f = new File(file);
			s.setName(f.getName());
		}
		return s;
	}

	public String getDbSearchResult(AFPChain afpChain){
		return afpChain.toDBSearchResult();
	}

	@Override
	public String printHelp() {
		StringBuffer buf = new StringBuffer();
		StructureAlignment alg = getAlgorithm();

		buf.append("-------------------").append(newline);
		buf.append(alg.getAlgorithmName() + " v." + alg.getVersion() + " help: " + newline);
		buf.append("-------------------").append(newline);
		buf.append(newline);

		buf.append(alg.getAlgorithmName()).append(" accepts the following parameters:").append(newline);
		buf.append(newline);

		buf.append("--- pairwise alignments ---").append(newline);
		buf.append(" two files to align can be specified by providing a path to a file, or a URL:").append(newline);
		buf.append("   -file1 the first file to align").append(newline);
		buf.append("   -file2 the second file to align").append(newline);
		buf.append(" alternatively you can specify PDB files by their PDB ids:").append(newline);
		buf.append("   -pdbFilePath  Path to the directory in your file system that contains the PDB files.").append(newline);
		buf.append("   -pdb1  PDB ID of target structure. Chain IDs are optional. In order to specify chain IDs write e.g: 5pti.A").append(newline);
		buf.append("   -pdb2  PDB ID of query structure. Chain IDs are optional. In order to specify chain IDs write e.g: 5pti.A").append(newline);
		buf.append(newline);

		buf.append("   -h / -help / --help : print this help string.").append(newline);
		buf.append("   -version: print version info").append(newline);
		buf.append("   -printXML true/false print the XML representation of the alignment on stdout.").append(newline);
		buf.append("   -printFatCat true/false print the original FATCAT output to stdout.").append(newline);
		buf.append("   -printCE true/false print the result in CE style").append(newline);
		buf.append("   -show3d print a 3D visualisation of the alignment (requires jmolapplet.jar in classpath)").append(newline);
		buf.append("   -outFile file to write the output to (default: writes XML representation).").append(newline);
		buf.append("   -outputPDB use this flag together with -outFile to dump the PDB file of the aligned structures, instead of the XML representation, instead of XML").append(newline);
		buf.append("   -autoFetch true/false if set to true PDB files will automatically get downloaded and stored in the right location. (default: false)").append(newline);
		buf.append("   -showMenu displays the menu that allows to run alignments through a user interface.").append(newline);
		buf.append(newline);

		buf.append("--- custom searches ---").append(newline);
		buf.append("   -alignPairs (mandatory) path to a file that contains a set of pairs to compair").append(newline);
		buf.append("   -outFile (mandatory) a file that will contain the summary of all the pairwise alignments").append(newline);
		buf.append(newline);

		buf.append("--- database searches ---").append(newline);
		buf.append("   -searchFile (mandatory) path to a PDB file that should be used in the search").append(newline);
		buf.append("   -outFile (mandatory) a directory that will contain the results of the DB search").append(newline);
		buf.append("   -nrCPU (optional) Number of CPUs to use for the database search. By default will use the all, but one CPU in the system.").append(newline);
		buf.append("   -pdbFilePath (mandatory) Path to the directory in your file system that contains the PDB files.").append(newline);
		buf.append("   -saveOutputDir (optional) a directory that will contain the detailed outputs of the alignments. By default will write XML files, if used together with -outputPDB, will write PDB files of the alignment.").append(newline);
		buf.append(newline);

		buf.append(" Once DB seaches are complete it is possible to view the results with:").append(newline);
		buf.append("   -showDBresult (optional) path to a DB outFile to show. Also provide the -pdbFilePath parameter to enable visualisation of results.").append(newline);
		buf.append(newline);

		ConfigStrucAligParams params = alg.getParameters();
		List<String> paramNames = params.getUserConfigParameters();
		List<String> paramHelp = params.getUserConfigHelp();

		assert(paramNames.size() == paramHelp.size());

		int size = Math.min(paramNames.size(), paramHelp.size());
		if(size > 0) {
			Iterator<String> namesIt = paramNames.iterator();
			Iterator<String> helpIt = paramHelp.iterator();

			buf.append("--- ").append(alg.getAlgorithmName()).append(" parameters: ---").append(newline);
			for(int i = 0; i< size; i++) {
				String name = namesIt.next();
				buf.append("   -").append(Introspector.decapitalize(name));
				buf.append(" ").append(helpIt.next());
				buf.append(newline);
			}
		}
		buf.append(newline);

		buf.append(" For boolean arguments: if neither the text >true< or >false< is provided it is assumed to mean >true<. Instead of >-argument false< it is also possible to write -noArgument.").append(newline);
		buf.append(newline);

		buf.append("--- How to specify what to align ---").append(newline);
		buf.append(" If only a PDB code is provided, the whole structure will be used for the alignment.").append(newline);
		buf.append(" To specify a particular chain write as: 4hhb.A (chain IDs are case sensitive, PDB ids are not)").append(newline);
		buf.append(" To specify that the 1st chain in a structure should be used write: 4hhb:0 .").append(newline);
		buf.append(" In order to align SCOP domains, provide pdb1/pdb2 as: d4hhba_ Note: if SCOP is not installed at the -pdbFilePath, will automatically download and install.").append(newline);
		buf.append(newline);

		return buf.toString();
	}

}
