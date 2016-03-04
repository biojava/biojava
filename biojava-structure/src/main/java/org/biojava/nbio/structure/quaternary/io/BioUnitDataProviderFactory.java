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
package org.biojava.nbio.structure.quaternary.io;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Factory to create BioUnitDataProvider instances.
 * 
 * Unlike many other BioJava Factory classes, this class does not store
 * singletons, but creates a new instance for every call of
 * {@link #getBioUnitDataProvider()}.
 */
public class BioUnitDataProviderFactory {

	private static final Logger logger = LoggerFactory.getLogger(BioUnitDataProviderFactory.class);
	
	public static final String mmcifProviderClassName     = MmCifBiolAssemblyProvider.class.getName();
	
	public static final String remoteProviderClassName    = RemoteBioUnitDataProvider.class.getName();
	
	public static final String pdbProviderClassName       = PDBBioUnitDataProvider.class.getName();
	
	public static Class<? extends BioUnitDataProvider> DEFAULT_PROVIDER_CLASS = MmCifBiolAssemblyProvider.class;
	public static final String DEFAULT_PROVIDER_CLASSNAME =  DEFAULT_PROVIDER_CLASS.getName();
	
	private static Class<? extends BioUnitDataProvider> providerClass = DEFAULT_PROVIDER_CLASS;
	
	private BioUnitDataProviderFactory(){
		
	}
	
	/**
	 * 
	 * @return A new instance of the current BioUnitDataProvider class
	 */
	public static BioUnitDataProvider getBioUnitDataProvider() {
		
		// use reflection to return a new instance...
		
		try {
			return providerClass.newInstance();
		} catch (IllegalAccessException e) {
			logger.error("Exception caught",e);
		} catch (InstantiationException e) {
			logger.error("Exception caught",e);			
		}
		
		return null;
		
	}

	/**
	 * Set the type of provider to be created
	 * @param klass A BioUnitDataProvider
	 */
	public static void setBioUnitDataProvider(Class<? extends BioUnitDataProvider> klass) {
		providerClass = klass;
	}
	/**
	 * Sets the data provider to the specified class name. Use {@link #setBioUnitDataProvider(Class)}
	 * for better type safety.
	 * @param className A class implementing BioUnitDataProvider
	 * @throws ClassNotFoundException If the class cannot be loaded
	 * @throws ClassCastException If the class does not extend BioUnitDataProvider
	 */
	public static void setBioUnitDataProvider(String className) throws ClassNotFoundException, ClassCastException {
		Class<?> cls = Class.forName(className);
		Class<BioUnitDataProvider> interfaceClass = BioUnitDataProvider.class;
		Class<? extends BioUnitDataProvider> castClass = cls.asSubclass(interfaceClass);

		setBioUnitDataProvider(castClass);
	}


}
