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
package org.biojava.nbio.structure.align.ce;


/** a simple bean that contains the parameters that can get set at startup
 *
 * @author Andreas Prlic
 *
 */
public class StartupParameters {



	String pdbFilePath;
	String cacheFilePath;
	String outFile;
	String pdb1;
	String pdb2;
	String file1;
	String file2;
	String showDBresult;
	boolean printXML;
	boolean printFatCat;
	boolean show3d;
	boolean autoFetch;
	boolean printCE;
	boolean showMenu;
	boolean printPDB;
	boolean isDomainSplit;


	// for DB searches
	String alignPairs;
	String searchFile;
	String saveOutputDir;
	int nrCPU;


	private static final String newline = System.getProperty("line.separator");

	public StartupParameters(){
		show3d = false;
		printXML = false;
		printPDB = false;
		printFatCat = false;
		autoFetch = false;
		showMenu = false;
		isDomainSplit = true;
		nrCPU = Runtime.getRuntime().availableProcessors() -1;
		if ( nrCPU < 1)
			nrCPU = 1;
	}

	/** An input file to be used for the DB search
	 *
	 * @return
	 */
	public String getSearchFile() {
		return searchFile;
	}
	public void setSearchFile(String searchFile) {
		this.searchFile = searchFile;
	}


	/** The file that contains a list of PDB pairs to be aligned
	 *
	 * @return
	 */
	public String getAlignPairs() {
		return alignPairs;
	}

	public void setAlignPairs(String alignPairs) {
		this.alignPairs = alignPairs;
	}

	public String getSaveOutputDir() {
		return saveOutputDir;
	}

	public void setSaveOutputDir(String saveOutputDir) {
		this.saveOutputDir = saveOutputDir;
	}

	public boolean isShowMenu() {
		return showMenu;
	}

	public void setShowMenu(boolean showMenu) {
		this.showMenu = showMenu;
	}

	/** Display the output string in CE style
	 *
	 * @return flag
	 */
	public boolean isPrintCE() {
		return printCE;
	}

	/** Display the output string in CE style
	 *
	 * @param printCE a flag
	 */
	public void setPrintCE(boolean printCE) {
		this.printCE = printCE;
	}


	public String getPdb1() {
		return pdb1;
	}
	/** mandatory argument to set the first PDB (and optionally chain ID) to be aligned.
	 *
	 * @param pdb1
	 */
	public void setPdb1(String pdb1) {
		this.pdb1 = pdb1;
	}
	public String getPdb2() {
		return pdb2;
	}

	/** mandatory argument to set the second PDB (and optionally chain ID) to be aligned.
	 *  @param pdb2
	 */
	public void setPdb2(String pdb2) {
		this.pdb2 = pdb2;
	}

	public boolean isPrintXML() {
		return printXML;
	}
	public void setPrintXML(boolean printXML) {
		this.printXML = printXML;
	}
	public boolean isPrintFatCat() {
		return printFatCat;
	}
	public void setPrintFatCat(boolean printFatCat) {
		this.printFatCat = printFatCat;
	}

	public String getPdbFilePath() {
		return pdbFilePath;
	}

	/** mandatory argument to set the location of PDB files.
	 *
	 * @param pdbFilePath
	 */
	public void setPdbFilePath(String pdbFilePath) {
		this.pdbFilePath = pdbFilePath;
	}

	public String getCacheFilePath() {
		return cacheFilePath;
	}

	public void setCacheFilePath(String cacheFilePath) {
		this.cacheFilePath = cacheFilePath;
	}

	public boolean isShow3d() {
		return show3d;
	}
	public void setShow3d(boolean show3d) {
		this.show3d = show3d;
	}
	public String getOutFile() {
		return outFile;
	}
	public void setOutFile(String outFile) {
		this.outFile = outFile;
	}
	public boolean isAutoFetch() {
		return autoFetch;
	}
	public void setAutoFetch(boolean autoFetch) {
		this.autoFetch = autoFetch;
	}
	public String getShowDBresult() {
		return showDBresult;
	}
	public void setShowDBresult(String showDBresult) {
		this.showDBresult = showDBresult;
	}

	public int getNrCPU() {
		return nrCPU;
	}
	public void setNrCPU(int nrCPU) {
		this.nrCPU = nrCPU;
	}

	public String getFile1()
	{
		return file1;
	}

	public void setFile1(String file1)
	{
		this.file1 = file1;
	}

	public String getFile2()
	{
		return file2;
	}

	public void setFile2(String file2)
	{
		this.file2 = file2;
	}



	/** When writing the results to a file, don;t write as XML but write aligned PDB file
	 *
	 * @return flag
	 */
	public boolean isOutputPDB() {
		return printPDB;
	}
	/** When writing the results to a file, don;t write as XML but write aligned PDB file
	 *
	 * @param printPDB flag to print aligned PDB
	 */
	public void setOutputPDB(boolean printPDB) {
		this.printPDB = printPDB;
	}





	public boolean isDomainSplit() {
		return isDomainSplit;
	}

	public void setDomainSplit(boolean isDomainSplit) {
		this.isDomainSplit = isDomainSplit;
	}

	@Override
	public String toString() {
		return "StartupParameters [pdbFilePath=" + pdbFilePath
				+ ", " + newline + " cacheFilePath=" + cacheFilePath + ", " + newline + " outFile=" + outFile
				+ ", " + newline + " pdb1=" + pdb1 + ", " + newline + " pdb2=" + pdb2 + ", " + newline + " file1=" + file1
				+ ", " + newline + " file2=" + file2 + ", " + newline + " showDBresult=" + showDBresult
				+ ", " + newline + " printXML=" + printXML + ", " + newline + " printFatCat=" + printFatCat
				+ ", " + newline + " show3d=" + show3d + ", " + newline + " autoFetch=" + autoFetch
				+ ", " + newline + " printCE=" + printCE + ", " + newline + " showMenu=" + showMenu
				+ ", " + newline + " printPDB=" + printPDB
				+ ", " + newline + " isDomainSplit="
				+ isDomainSplit + ", " + newline + " alignPairs=" + alignPairs
				+ ", " + newline + " searchFile=" + searchFile + ", " + newline + " saveOutputDir="
				+ saveOutputDir + ", " + newline + " nrCPU=" + nrCPU + "]";
	}







}
