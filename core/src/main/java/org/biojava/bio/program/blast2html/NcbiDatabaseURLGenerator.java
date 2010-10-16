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
package org.biojava.bio.program.blast2html;

import java.util.Properties;


/**
 * Simple URL generator for Entrez at the NCBI.
 *
 * Primary author -
 *                 Colin Hardman      (CAT)
 * Other authors  -
 *                 Tim Dilks          (CAT)
 *                 Simon Brocklehurst (CAT)
 *                 Stuart Johnston    (CAT)
 *                 Lawerence Bower    (CAT)
 *                 Derek Crockford    (CAT)
 *                 Neil Benn          (CAT)
 *
 * Copyright 2001 Cambridge Antibody Technology Group plc.
 *
 * This code released to the biojava project, May 2001
 * under the LGPL license.
 *
 * @author Cambridge Antibody Technology Group plc
 * @version 1.0
 */
public class NcbiDatabaseURLGenerator implements DatabaseURLGenerator {

    public String toURL( String poID,  Properties poProperties ) {

	int index = poID.indexOf( ";" );
	if ( index != -1 ) {
	    poID = poID.substring( 0, index );
	}

	return "http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Search&db=".concat( poProperties.getProperty( "db" )).
	    concat( "&term=").concat( poID );
    }

    public String toLink( String poID, Properties poProperties ) {
	return "[<A class=\"dbRetrieve\" HREF=\"".concat
	    ( this.toURL( poID, poProperties ) ).concat
	    ( "\">retrieve hit from ncbi</A>]");
    }
}
