/*
 *                  BioJava development code
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
 * Created on Oct 3, 2007
 * 
 */

package org.biojava.dasobert.das;

import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.Date;


/** a simple class that converts a java.util.Date into the DAS
 * Date and TIme convention
 * 
 * @author Andreas Prlic
 *
 */
public class DasTimeFormat {

	public static final String DASFORMAT="yyyy-MM-dd'T'HH:mm:ssZ";
	
	
	public static String toDASString(Date date){
		
		SimpleDateFormat format = new SimpleDateFormat(DASFORMAT);
		
		return format.format(date);
		
	}
	
	public static Date fromDASString(String dasTimeFormat) throws ParseException{
		
		SimpleDateFormat format = new SimpleDateFormat(DASFORMAT);
		
		Date d = new Date(0);
				
		d = format.parse(dasTimeFormat);
		
		return d;
	}
	
}
