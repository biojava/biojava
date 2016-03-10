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


import java.util.List;

public interface ConfigStrucAligParams {


	/** get the list of parameters that the user can change through the user interface.
	 *  Parameter names are the same names as the corresponding Get/Set methods.
	 *
	 * @return list of parameters
	 */
	public List<String> getUserConfigParameters();

	/** The labels to be displayed to the user for each parameter
	 *
	 * @return list of parameter names
	 */
	public List<String> getUserConfigParameterNames();

	/** Get the data types of the parameters
	 *
	 * @return list of parameter classes
	 */
	@SuppressWarnings("rawtypes")
	public List<Class> getUserConfigTypes();


	/** The help text for each of these parameters.
	 *
	 * @return help strings
	 */
	public List<String> getUserConfigHelp();


	/** Set the parameters to the default.
	 *
	 */
	public void reset();
}
