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
 * created at May 31, 2008
 */
package org.biojava.bio.structure.io.mmcif.model;

import java.lang.reflect.Method;
import java.util.List;

import org.biojava.bio.structure.Chain;

/** a generic class that implements the toString method for a bean
 * 
 * @author Andreas Prlic
 *
 */ 
public abstract class AbstractBean {

	@SuppressWarnings({ "rawtypes", "unchecked" })
	public String toString(){
		StringBuffer buf = new StringBuffer();
        buf.append(this.getClass().getName()).append(": ");
		/* disabled for the moment

    	 buf.append(" chains: " );
    	Iterator<Chain> iter = chainList.iterator();
    	while (iter.hasNext()){
    		Chain c = iter.next();
    		buf.append (c.getName() + " ");
    	}

		 */
		try {
			Class c = this.getClass();
			Method[] methods  = c.getMethods();

			for (int i = 0; i < methods.length; i++) {
				Method m = methods[i];     

				String name = m.getName();
				if ( name.substring(0,3).equals("get")) {				 

					Object o  = m.invoke(this, new Object[]{});
					if ( o instanceof String){
                        buf.append(name.substring(3, name.length())+": "+ o + " ");
					}
					else if ( o instanceof List){
                        buf.append(name.substring(3, name.length())).append(": ");

						List<Object>lst = (List<Object>)o;
						for (Object obj : lst){
							if ( obj instanceof Chain){
								continue;
							}
                            buf.append(obj).append(" ");
						}

					} 
					else {
						// ignore...
					}
				}

			}

		} catch (Exception e){
			e.printStackTrace();
		}


		//if ( organismScientific != null)
		//	buf.append(" organism scientific: " + organismScientific);


		return buf.toString();
	}
	
}
