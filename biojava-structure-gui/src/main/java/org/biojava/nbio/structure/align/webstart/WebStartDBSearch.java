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
package org.biojava.nbio.structure.align.webstart;

import org.biojava.nbio.structure.align.FarmJob;
import org.biojava.nbio.structure.align.client.FarmJobParameters;
import org.biojava.nbio.structure.align.gui.GUIFarmJobRunnable;
import org.biojava.nbio.structure.align.util.CliTools;
import org.biojava.nbio.structure.align.util.ConfigurationException;

import javax.swing.*;
import java.util.Arrays;
import java.util.List;




/** A Web Start wrapper for a FarmJobRunnable.
 *
 */
public class WebStartDBSearch  {

	private static final String[] mandParams = new String[] {"pdbFilePath"};

	private static final List<String> mandatoryArgs= Arrays.asList(mandParams);

	public WebStartDBSearch(){
	}



	public static void main(String[] argv) {

		FarmJob job = new FarmJob();


		if (argv.length  == 0 ) {
			job.printHelp();
			JOptionPane.showMessageDialog(null,
					"Not enough arguments!");
			return;


		}

		if ( argv.length == 1){
			if (argv[0].equalsIgnoreCase("-h") || argv[0].equalsIgnoreCase("-help")|| argv[0].equalsIgnoreCase("--help")){
				job.printHelp();
				JOptionPane.showMessageDialog(null,
				"Help not supported...");
				return;
			}
		}

		FarmJobParameters params = new FarmJobParameters();


		for (int i = 0 ; i < argv.length; i++){
			String arg   = argv[i];

			String value = null;
			if ( i < argv.length -1)
				value = argv[i+1];

			// if value starts with - then the arg does not have a value.
			if (value != null && value.startsWith("-"))
				value = null;
			else
				i++;


			String[] tmp = {arg,value};

			try {

				CliTools.configureBean(params, tmp);

			} catch (ConfigurationException e){

				e.printStackTrace();

				if ( mandatoryArgs.contains(arg) ) {
					// there must not be a ConfigurationException with mandatory arguments.
					JOptionPane.showMessageDialog(null,
							e.getMessage());
					return;

				} else {
					// but there can be with optional ...
				}
			}
		}

		params.setRunBackground(true);
		GUIFarmJobRunnable runnable = new GUIFarmJobRunnable(params);

		//javax.swing.SwingUtilities.invokeLater(runnable);
		runnable.run();




	}


}
