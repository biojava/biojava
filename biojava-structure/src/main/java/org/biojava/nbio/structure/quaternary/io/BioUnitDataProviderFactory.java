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

public class BioUnitDataProviderFactory {

	private static final Logger logger = LoggerFactory.getLogger(BioUnitDataProviderFactory.class);
	
	public static final String mmcifProviderClassName     = MmCifBiolAssemblyProvider.class.getName();
	
	public static final String remoteProviderClassName    = RemoteBioUnitDataProvider.class.getName();
	
	public static final String pdbProviderClassName       = PDBBioUnitDataProvider.class.getName();
	
	public static String DEFAULT_PROVIDER_CLASSNAME =  mmcifProviderClassName;
	
	private static String providerClassName = DEFAULT_PROVIDER_CLASSNAME;
	
	private BioUnitDataProviderFactory(){
		
	}
	
	public static BioUnitDataProvider getBioUnitDataProvider() {
		
		// use reflection to return a new instance...
		
		try {
			Class<?> cls = Class.forName(providerClassName); 
			//System.out.println("Using BioUnitProvider: " + providerClassName);
			return (BioUnitDataProvider) cls.newInstance();
		} catch (ClassNotFoundException e) {
			logger.error("Exception caught",e);
		} catch (IllegalAccessException e) {
			logger.error("Exception caught",e);
		} catch (InstantiationException e) {
			logger.error("Exception caught",e);			
		}
		
		return null;
		
	}

	public static void setBioUnitDataProvider(String className) {
		
		try {
			Class<?> cls = Class.forName(providerClassName);
			Class<?> interfaceClass = Class.forName(BioUnitDataProvider.class.getName());
			Class<?>[] ifs = cls.getInterfaces();
			boolean found = false;
			for ( Class<?> c : ifs){
				if ( c.equals(interfaceClass)){
					found = true;
					break;
				}
				
			}
			if ( ! found){
				logger.warn("The provided class {} does not implement the correct interface!", className);
				return ;
			}
		} catch (ClassNotFoundException e) {
			logger.error("Exception caught",e);
			
			// we do not set classes that don't exist!
			return;
		
		}
		providerClassName = className;
	}
	
	
}
