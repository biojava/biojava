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
/*
 *                  BioJava development code
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

package org.biojava.nbio.structure.align.util;

import java.beans.BeanInfo;
import java.beans.Introspector;
import java.beans.PropertyDescriptor;
import java.io.*;
import java.lang.reflect.Array;
import java.lang.reflect.InvocationTargetException;
import java.util.*;

/**
 * Utilities for autoconfiguring javabeans based on command line arguments.
 *
 * @author Thomas Down
 */

public class CliTools {
	private CliTools() {
	}
	

	/**
	 * Configure a JavaBean based on a set of command line arguments.  
	 * For a command line construct such as "-foo 42", this method will use
	 * available <code>BeanInfo</code> (usually obtained by introspection)
	 * to find a property named "foo".  The argument will be interpreted
	 * according to the type of the "foo" property, then the appropriate
	 * mutator method (generally named setFoo) will be called to configure
	 * the property on the bean.
	 *
	 * <p>
	 * Currently supported property types are <code>int, double,
	 * boolean, String, File, Reader, Writer, InputStream, OutputStream, Enum</code>,
	 * plus arrays of all the above types.  In the case of arrays, the option
	 * may appear multiple times on the command line, otherwise recurrance of
	 * the same option is an error.
	 * </p>
	 *
	 * <p>
	 * For stream types, the parameter is interpreted as a filename unless it
	 * is equal to "-" in which case standard input or standard output are
	 * used as appropriate.  Each of the standard streams may only be used
	 * one.
	 * </p>
	 *
	 * <p>
	 * In the future, this method will probably be extended to handle multiple
	 * parameter occurances, and use Annotations to generate more useful help
	 * messages when something goes wrong.
	 * </p>
	 *
	 * @param bean 
	 * @param args
	 * @return A string array which contains any 'anonymous' arguments (may be empty)
	 * @throws ConfigurationException
	 */

	@SuppressWarnings({ "rawtypes", "unchecked" })
	public static String[] configureBean(Object bean, String[] args) 
	throws ConfigurationException
	{
		BeanInfo bi;
		try {
			bi = Introspector.getBeanInfo(bean.getClass());
		} catch (Exception ex) {
			throw new ConfigurationException("Couldn't get information for target bean " + ex.getMessage());
		}

		Map<String,PropertyDescriptor> propertiesByName = new HashMap<String, PropertyDescriptor>();
		for (PropertyDescriptor pd : bi.getPropertyDescriptors() ) {
			propertiesByName.put(pd.getName(), pd);
		}

		List<String> anonArgs = new ArrayList<String>();
		Map<PropertyDescriptor,List<String>> arrayProps = new HashMap<PropertyDescriptor, List<String>>();
		Set<PropertyDescriptor> usedProps = new HashSet<PropertyDescriptor>();

		boolean stdInUsed = false;
		boolean stdOutUsed = false;

		for (int i = 0; i < args.length; ++i) {
			String arg = args[i];
			
			//System.out.println("checking argument: " + arg);
			
			if ( arg == null)
				continue;
			
			if ((arg.length() > 0 ) && arg.charAt(0) == '-') {
				PropertyDescriptor pd = propertiesByName.get(arg.substring(1));

				boolean arrayMode = false;
				Object propVal = null;
				Class propType = null;

				if (pd == null) {
					if (arg.startsWith("-no")) {
						String altPropName = Introspector.decapitalize(arg.substring(3));
						pd = propertiesByName.get(altPropName);
						if (pd == null) {
							throw new ConfigurationException("No property named " + arg.substring(1) + " or " + altPropName);
						}
						propType = pd.getPropertyType();
						if (propType == Boolean.TYPE) {
							propVal = Boolean.FALSE;
						} else {
							throw new ConfigurationException("Negatory option " + arg + " does not refer to a boolean property");
						}
					} else {
						throw new ConfigurationException("No property named " + arg.substring(1));
					}
				} else {
					propType = pd.getPropertyType();

					if (propType.isArray()) {
						arrayMode = true;
						propType = propType.getComponentType();
					}

					if (propType == Integer.TYPE) {
						try {
							propVal = new Integer(args[++i]);
						} catch (Exception ex) {
							throw new ConfigurationException("Option " + arg + " requires an integer parameter");
						}
					} else if (propType == Double.TYPE || propType == Double.class ) {
						try {
							propVal = new Double(args[++i]);
						} catch (Exception ex) {
							throw new ConfigurationException("Option " + arg + " requires a numerical parameter");
						}
					} else if (propType == String.class) {
						propVal = args[++i];
					} else if (propType == Boolean.TYPE) {

						String val = args[++i];
						if ( val == null )
							propVal = Boolean.TRUE;
						else {
							if ( val.equalsIgnoreCase("true") || val.equalsIgnoreCase("t"))
								propVal = Boolean.TRUE;
							else if( val.equalsIgnoreCase("false") || val.equalsIgnoreCase("f"))
								propVal = Boolean.FALSE;
							else
								throw new ConfigurationException("Option "+arg+" requires a boolean parameter");
						}
					} else if (File.class.isAssignableFrom(propType)) {
						// can't distinguish if the file is for reading or writing
						// at the moment, so accept it without validation.
						propVal = new File(args[++i]);
					} else if (Reader.class.isAssignableFrom(propType)) {
						String name = args[++i];
						if ("-".equals(name)) {
							if (stdInUsed) {
								throw new ConfigurationException("Can't use standard input more than once");
							}
							propVal = new InputStreamReader(System.in);
							stdInUsed = true;
						} else {
							try {
								propVal = new FileReader(new File(name));
							} catch (Exception ex) {
								throw new ConfigurationException("Can't open " + name + " for input");
							}
						}
					} else if (InputStream.class.isAssignableFrom(propType)) {
						String name = args[++i];
						if ("-".equals(name)) {
							if (stdInUsed) {
								throw new ConfigurationException("Can't use standard input more than once");
							}
							propVal = System.in;
							stdInUsed = true;
						} else {
							try {
								propVal = new FileInputStream(new File(name));
							} catch (Exception ex) {
								throw new ConfigurationException("Can't open " + name + " for input");
							}
						}
					} else if (Writer.class.isAssignableFrom(propType)) {
						String name = args[++i];
						if ("-".equals(name)) {
							if (stdOutUsed) {
								throw new ConfigurationException("Can't use standard output more than once");
							}
							propVal = new OutputStreamWriter(System.out);
							stdOutUsed = true;
						} else {
							try {
								propVal = new FileWriter(new File(name));
							} catch (Exception ex) {
								throw new ConfigurationException("Can't open " + name + " for output");
							}
						}
					} else if (OutputStream.class.isAssignableFrom(propType)) {
						String name = args[++i];
						if ("-".equals(name)) {
							if (stdOutUsed) {
								throw new ConfigurationException("Can't use standard output more than once");
							}
							propVal = System.out;
							stdOutUsed = true;
						} else {
							try {
								propVal = new FileOutputStream(new File(name));
							} catch (Exception ex) {
								throw new ConfigurationException("Can't open " + name + " for output");
							}
						}
					} else if( propType.isEnum()) {
						String name = args[++i];
						try {
							propVal = Enum.valueOf(propType, name);
						} catch (Exception ex) {
							try {
								// Try with uppercase version, as common for enums
								propVal = Enum.valueOf(propType, name.toUpperCase());
							}  catch (Exception ex2) {
								//give up
								StringBuilder errMsg = new StringBuilder();
								errMsg.append("Option ").append(arg);
								errMsg.append(" requires a ").append(propType.getSimpleName());
								errMsg.append(" parameter. One of: ");
								for(Object val: propType.getEnumConstants() ) {
									Enum enumVal = (Enum) val;
									errMsg.append(enumVal.name());
									errMsg.append(" ");
								}
								throw new ConfigurationException(errMsg.toString());
							}
						}

					} else {
						System.err.println("Unsupported optionType for " + arg + " propType:" + propType);
						System.exit(1);
					}
				}
				
				//System.out.println("setting to: " + propVal + " " + propVal.getClass().getName());

				if (arrayMode) {
					List valList = arrayProps.get(pd);
					if (valList == null) {
						valList = new ArrayList();
						arrayProps.put(pd, valList);
					}
					valList.add(propVal);
				} else {
					if (usedProps.contains(pd)) {
						throw new ConfigurationException("Multiple values supplied for " + pd.getName());
					}
					try {
						pd.getWriteMethod().invoke(bean, new Object[] {propVal});
					} catch (InvocationTargetException ex) {
						throw new ConfigurationException("Error configuring '" + pd.getName() + "'");
					} catch (Exception ex) {
						throw new ConfigurationException("Error configuring '" + pd.getName() + "'");
					}
					usedProps.add(pd);
				}
			} else {
				anonArgs.add(arg);
			}
		}

		for (Iterator api = arrayProps.entrySet().iterator(); api.hasNext(); ) {
			Map.Entry me = (Map.Entry) api.next();
			PropertyDescriptor pd = (PropertyDescriptor) me.getKey();
			List vals = (List) me.getValue();

			Class compType = pd.getPropertyType().getComponentType();
			Object valArray;
			if (compType.isPrimitive()) {
				if (compType == Integer.TYPE) {
					valArray = CollectionTools.toIntArray(vals);
				} else if (compType == Double.TYPE) {
					valArray = CollectionTools.toDoubleArray(vals);
				} else {
					throw new ConfigurationException("Arrays of type " + compType.getName() + " are currently unsupported");
				}
			} else {
				valArray = vals.toArray((Object[]) Array.newInstance(compType, vals.size()));
			}
			try {
				pd.getWriteMethod().invoke(bean, new Object[] {valArray});
			} catch (InvocationTargetException ex) {
				throw new ConfigurationException("Error configuring '" + pd.getName() + "'");
			} catch (Exception ex) {
				throw new ConfigurationException("Error configuring '" + pd.getName() + "'");
			}
		}

		return anonArgs.toArray(new String[anonArgs.size()]);
	}

	/**
	 * Constructs a comma-separated list of values for an enum.
	 * 
	 * Example:
	 * > getEnumValues(ScoringStrategy.class)
	 * "CA_SCORING, SIDE_CHAIN_SCORING, SIDE_CHAIN_ANGLE_SCORING, CA_AND_SIDE_CHAIN_ANGLE_SCORING, or SEQUENCE_CONSERVATION"
	 * @param enumClass
	 * @return
	 */
	public static <T extends Enum<?>> String getEnumValuesAsString(Class<T> enumClass) {
		//ScoringStrategy[] vals = ScoringStrategy.values();
		T[] vals = enumClass.getEnumConstants();

		StringBuilder str = new StringBuilder();
		if(vals.length == 1) {
			str.append(vals[0].name());
		} else if(vals.length > 1 ) {
			for(int i=0;i<vals.length-1;i++) {
				str.append(vals[i].name());
				str.append(", ");
			}
			str.append("or ");
			str.append(vals[vals.length-1].name());
		}

		return str.toString();
	}

}
