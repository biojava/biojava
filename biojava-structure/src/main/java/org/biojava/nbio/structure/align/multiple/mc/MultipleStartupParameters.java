package org.biojava.nbio.structure.align.multiple.mc;

import java.util.List;

/** 
 * Bean that contains the parameters that can get set at startup for {@link MultipleStructureAligner}s.<p>
 * Note that for Multiple Alignments no DB searches are supported.<p>
 * More options can be added when needed. TODO
 * 
 * @author Aleix Lafita
 *
 */
public class MultipleStartupParameters {

	String pdbFilePath;
	String cacheFilePath;
	String outFile;
	List<String> pdbs;
	List<String> files;
	boolean printFatCat;
	boolean show3d;
	boolean autoFetch;
	boolean showMenu;
	boolean printPDB;
	
	public MultipleStartupParameters(){
		show3d = false;
		printPDB = false;
		printFatCat = false;
		autoFetch = false;
		showMenu = false;
	}

	@Override
	public String toString() {
		return "MultipleStartupParameters [pdbFilePath=" + pdbFilePath
				+ ", cacheFilePath=" + cacheFilePath + ", outFile=" + outFile
				+ ", pdbs=" + pdbs + ", files=" + files + ", printFatCat="
				+ printFatCat + ", show3d=" + show3d + ", autoFetch="
				+ autoFetch + ", showMenu=" + showMenu + ", printPDB="
				+ printPDB + "]";
	}

	public String getPdbFilePath() {
		return pdbFilePath;
	}

	public void setPdbFilePath(String pdbFilePath) {
		this.pdbFilePath = pdbFilePath;
	}

	public String getCacheFilePath() {
		return cacheFilePath;
	}

	public void setCacheFilePath(String cacheFilePath) {
		this.cacheFilePath = cacheFilePath;
	}

	public String getOutFile() {
		return outFile;
	}

	public void setOutFile(String outFile) {
		this.outFile = outFile;
	}

	public List<String> getPdbs() {
		return pdbs;
	}

	public void setPdbs(List<String> pdbs) {
		this.pdbs = pdbs;
	}

	public List<String> getFiles() {
		return files;
	}

	public void setFiles(List<String> files) {
		this.files = files;
	}

	public boolean isPrintFatCat() {
		return printFatCat;
	}

	public void setPrintFatCat(boolean printFatCat) {
		this.printFatCat = printFatCat;
	}

	public boolean isShow3d() {
		return show3d;
	}

	public void setShow3d(boolean show3d) {
		this.show3d = show3d;
	}

	public boolean isAutoFetch() {
		return autoFetch;
	}

	public void setAutoFetch(boolean autoFetch) {
		this.autoFetch = autoFetch;
	}

	public boolean isShowMenu() {
		return showMenu;
	}

	public void setShowMenu(boolean showMenu) {
		this.showMenu = showMenu;
	}

	public boolean isPrintPDB() {
		return printPDB;
	}

	public void setPrintPDB(boolean printPDB) {
		this.printPDB = printPDB;
	}
}
