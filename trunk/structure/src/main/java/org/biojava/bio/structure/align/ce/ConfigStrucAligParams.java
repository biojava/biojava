package org.biojava.bio.structure.align.ce;


import java.util.List;

public interface ConfigStrucAligParams {

	
	/** get the list of parameters that the user can change through the user interface.
	 *  Parameter names are the same names as the corresponding Get/Set methods.
	 * 
	 * @return
	 */
	public List<String> getUserConfigParameters();

	/** The labels to be displayed to the user for each parameter
	 * 
	 * @return
	 */
	public List<String> getUserConfigParameterNames();
	
	/** Get the data types of the parameters
	 * 
	 * @return
	 */
	@SuppressWarnings("unchecked")
	public List<Class> getUserConfigTypes();
	
	
	/** The help text for each of these parameters.
	 * 
	 * @return
	 */
	public List<String> getUserConfigHelp();
	
	
	/** Set the parameters to the default.
	 * 
	 */
	public void reset();
}
