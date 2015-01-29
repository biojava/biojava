/* 
 * @(#)NullOutputStream.java	1.0 June 2010
 * 
 * Copyright (c) 2010 Peter Troshin
 * 
 * JRONN version: 3.1     
 *  
 *        BioJava development code
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
 
package org.biojava3.ronn;

import java.io.IOException;
import java.io.OutputStream;

/**
 * The stream that void its input
 * 
 * @author Petr Troshin
 * @version 1.0
 * @since 3.0.2
 */
public final class NullOutputStream extends OutputStream {

    @Override
    public void write(final int b) throws IOException {
	// this methods does nothing.
	// This is an intention
    }

}
