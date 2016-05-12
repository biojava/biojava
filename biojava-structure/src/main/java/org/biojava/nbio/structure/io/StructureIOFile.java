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
package org.biojava.nbio.structure.io;

import org.biojava.nbio.structure.Structure;

import java.io.File;
import java.io.IOException;
import java.util.List;

/**
 *  StructureIOFile extends StructureProvider with methods specific to
 *  parsing files from the filesystem.
 * @author Andreas Prlic
 */
public interface StructureIOFile extends StructureProvider {

	/**
	 * Associates a file extension with this particular StructureIOFile,
	 * indicating that files of that type can be parsed. This is generally
	 * called only in the constructor of the implementing class.
	 * @param ext  a String ...
	 */
	public void addExtension(String ext);

	/**
	 * Returns a list of extensions supported by this class
	 * @return a (potentially empty) list of strings
	 */
	public List<String> getExtensions();

	/**
	 * Open filename and return a Structure object.
	 *
	 * Not to be confused with {@link #getStructureById(String)}
	 * @param filename  The path to the file. Must be the correct format for the
	 *  implementing class.
	 * @return a Structure object
	 * @throws IOException ...
	 */
	public Structure getStructure(String filename) throws IOException ;

	/**
	 * Read file from File and returns
	 * a Structure object.
	 * @param file file containing the structure. Must be the correct format for
	 *  the implementing class
	 * @return a Structure object
	 * @throws IOException ...
	 */
	public Structure getStructure(File file) throws IOException ;
}
