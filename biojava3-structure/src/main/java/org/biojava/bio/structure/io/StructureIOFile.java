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
 * Created on 26.04.2004
 * @author Andreas Prlic
 *
 */
package org.biojava.bio.structure.io;

import java.io.File;
import java.io.IOException;

import org.biojava.bio.structure.Structure;

/**
 * interface StructureIOFile extends the StructureIO interface
 * and adds a few File specific methods.
 * @author Andreas Prlic
 */
public interface StructureIOFile extends StructureIO {

    /** Set path to file / connection string to db.
     * This is for installations of PDB/mmCif where all files are located in one directory.
     *
     * @param path  a String specifying the path value
     */
    public void setPath(String path) ;

    /** get the directory path to the files
     *
     * @return path
     */
    public String getPath();

    /** add a known File extension.
     * @param ext  a String ...
     */
    public void addExtension(String ext);

    /** clear all file extensions
     *
     */
    public void clearExtensions();

    /** open filename and returns
     * a Structure object.
     * @param filename  a String
     * @return a Structure object
     * @throws IOException ...
     */
    public Structure getStructure(String filename) throws IOException ;

    /** read file from File and returns
     * a Structure object.
     * @param file file containing a PDB or mmcif file
     * @return a Structure object
     * @throws IOException ...
     */
    public Structure getStructure(File file) throws IOException ;

    /** Fetch files automatically from FTP server. Default: false
     *
     * @return flag is true or false.
     * 
     */
    public boolean isAutoFetch();
       
    /** Tell the parser to fetch missing PDB files from the FTP server automatically.
	 *
	 * default is false. If true, new PDB files will be automatically stored in the Path and gzip compressed.
	 *
	 * @param autoFetch flag.
	 */ 
	public void setAutoFetch(boolean autoFetch);

	
	/** The PDB files are organized hierarchically (as on the PDB - FTP server. Directories are split based on the two middle characters of the files).
	 * 
	 * @param isSplit
	 * 
	 */
	public void setPdbDirectorySplit(boolean isSplit);
	
	/** The PDB files are organized hierarchically (as on the PDB - FTP server. Directories are split based on the two middle characters of the files).
	 * 
	 * @return flag
	 *
	 */
	public boolean isPdbDirectorySplit();
		
	
	/** Get a Structure based on its PDB id. The reader takes care of finding the correct file in the PATH configured in get/setPath.
	 * @return a Structure object
	 */
	public Structure getStructureById(String pdbId) throws IOException;

	
	/** Set the parameters that should be used for file parsing
	 * 
	 * @param params FileParsingParameters
	 */
	public void setFileParsingParameters(FileParsingParameters params);
	   
	
	/** Get the parameters that should be used for file parsing
     * 
     * @return the FileParsingParameters that are configuring the behavior of the parser
     */
    public FileParsingParameters getFileParsingParameters();
    
	
	
}
