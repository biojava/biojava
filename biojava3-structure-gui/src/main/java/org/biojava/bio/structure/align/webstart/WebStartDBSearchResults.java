/**
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
 * Created on Sep 26, 2011
 * Created by Andreas Prlic
 *
 * @since 3.0.2
 */
package org.biojava.bio.structure.align.webstart;

import java.io.File;
import java.net.URL;

import javax.swing.JOptionPane;

import org.biojava.bio.structure.align.gui.DBResultTable;
import org.biojava.bio.structure.align.util.UserConfiguration;

public class WebStartDBSearchResults {

	public static void main(String[] argv){

		if (argv.length  == 0 ) {

			JOptionPane.showMessageDialog(null,
			"Not enough arguments!");
			return;


		} else if ( argv.length == 2){
			String path =  argv[1];

			DBResultTable table = new DBResultTable();
			UserConfiguration config = WebStartMain.getDefaultConfig();
			try {
				URL u = new URL(path);

				//File f = new File(u.toURI());

				table.show(u,config);
			} catch (Exception e){
				JOptionPane.showMessageDialog(null,
						e.getMessage());
				return;
			}
		}


	}
}
