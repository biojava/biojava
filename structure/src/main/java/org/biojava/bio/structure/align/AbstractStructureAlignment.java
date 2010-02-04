package org.biojava.bio.structure.align;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.ce.ConfigStrucAligParams;
import org.biojava.bio.structure.align.model.AFPChain;

public abstract class AbstractStructureAlignment implements StructureAlignment {

	public static String newline = System.getProperty("line.separator");

	abstract public  AFPChain align(Atom[] ca1, Atom[] ca2) throws StructureException ;

	abstract public AFPChain align(Atom[] ca1, Atom[] ca2, Object params) throws StructureException;

	abstract public String getAlgorithmName() ;

	abstract public ConfigStrucAligParams getParameters() ;

	abstract public String getVersion() ;

	abstract public void setParameters(ConfigStrucAligParams parameters);

	public String printHelp() {
		StringBuffer buf = new StringBuffer();

		buf.append("-------------------" + newline);
		buf.append(getAlgorithmName() + " v." + getVersion() + " help: " + newline);
		buf.append("-------------------" + newline);
		buf.append(newline);
		buf.append(getAlgorithmName() + " accepts the following parameters:" + newline);
		buf.append(newline);
		buf.append("--- pairwise alignments ---");
		buf.append(" two files to align can be specified by providing a path to a file, or a URL:" + newline);
		buf.append("-file1 the first file to align" + newline);
		buf.append("-file2 the second file to align"+newline);
		buf.append("alternatively you can specify PDB files by their PDB ids:" + newline);		
		buf.append("-pdbFilePath  Path to the directory in your file system that contains the PDB files." + newline);
		buf.append("-pdb1  PDB ID of target structure. Chain IDs are optional. In order to specify chain IDs write e.g: 5pti.A" + newline);
		buf.append("-pdb2  PDB ID of query structure. Chain IDs are optional. In order to specify chain IDs write e.g: 5pti.A" + newline);
		
		buf.append("-h / -help / --help : print this help string." + newline);
		buf.append("-printXML true/false print the XML representation of the alignment on stdout." + newline);
		buf.append("-printFatCat true/false print the original FATCAT output to stdout." + newline);
		buf.append("-printCE true/false print the result in CE style" + newline);
		buf.append("-show3d print a 3D visualisation of the alignment (requires jmolapplet.jar in classpath)" + newline);
		buf.append("-outFile file to write the output to (writes XML representation)." + newline);
		buf.append("-autoFetch true/false if set to true PDB files will automatically get downloaded and stored in the right location. (default: false)" + newline);		
		buf.append("-pdbDirSplit true/false the directory containing PDB files has all PDBs in one level or is split into multiple subdirs, like the ftp site. (default: true)" + newline);
		buf.append("-showMenu displays the menu that allows to run alignments through a user interface." + newline);
		buf.append("-maxGapSize (jCE specific): set the maximum gap size parameter G during AFP extension. default: 30. Set to 0 for unrestricted. " + newline);
		buf.append("-showAFPRanges (jCE specific): show the raw Aligned Fragment Pair positions, prior to optimization." + newline);
		buf.append(newline);
		buf.append("--- database searches ---");
		buf.append(newline);
		buf.append("-alignPairs (mandatory) path to a file that contains a set of pairs to compair" + newline);		
		buf.append("-outFile (mandatory) a file that will contain the summary of all the pairwise alignments" + newline);
		buf.append("-pdbFilePath (mandatory) Path to the directory in your file system that contains the PDB files." + newline);
		buf.append(newline);
		buf.append("  Once DB seaches are complete it is possible to view the results with:" +newline);
		buf.append(newline);
		buf.append("-viewDBresult (optional) path to a DB outFile to show. Also provide the -pdbFilePath parameter to enable visualisation of results." + newline);
		buf.append(newline);
		buf.append("For boolean arguments: if neither the text >true< or >false< is provided it is assumed to mean >true<. Instead of >-argument false< it is also possible to write -noArgument." + newline);		
		buf.append("--- How to specify what to align ---");
		buf.append(newline);
		buf.append(" If only a PDB code is provided, the whole structure will be used for the alignment." + newline);
		buf.append(" To specify a particular chain write as: 4hhb.A (chain IDs are case sensitive, PDB ids are not)"  + newline);
		buf.append(" To specify that the 1st chain in a structure should be used write: 4hhb:0 ." + newline);

		return buf.toString();


	}



}
